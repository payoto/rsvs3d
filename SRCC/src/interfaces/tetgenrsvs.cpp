#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <unordered_map>
#include <cmath>

#include "tetgen.h"
#include "warning.hpp"
#include "arraystructures.hpp"
#include "mesh.hpp"
#include "snake.hpp"
#include "voxel.hpp"
#include "postprocessing.hpp"
#include "meshrefinement.hpp"
#include "meshprocessing.hpp"
#include "tetgenrsvs.hpp"
#include "rsvsjson.hpp"

using json = rsvsjson::json;

void tetgen::internal::MeshData2Tetgenio(const mesh &meshgeom, 
	tetgen::io_safe &tetin,	int facetOffset, int pointOffset, 
	int pointMarker, const std::vector<double> &pointMtrList,
	const std::vector<double> &facetConstr,	int facetConstrOffset){
	/*
	Writes meshdata into the tetgenio format for a single mesh
	*/
	vector<int> orderVert;
	int countI, countJ, countK, nFacetConstr;

	// Set point properties to the appropriate lists
	countI = meshgeom.verts.size();
	countJ = min(int(pointMtrList.size()), tetin.numberofpointmtrs);
	for (int i = 0; i < countI; ++i) {
		// set point coordinates
		for (int j = 0; j < 3; ++j) {
			tetin.pointlist[(i+pointOffset)*3+j] = meshgeom.verts(i)->coord[j];
		}
		// set point metrics
		for (int j = 0; j < countJ; ++j) {
			tetin.pointmtrlist[(i+pointOffset)*countJ+j] =pointMtrList[j];
		}
		// set point markers
		tetin.pointmarkerlist[i+pointOffset] = pointMarker;
	}

	// Set facet properties to the appropriate lists
	nFacetConstr = facetConstr.size();
	// std::cout << "Number of facetconstraint " << nFacetConstr << std::endl;
	for (int i = 0; i < nFacetConstr; ++i) {
		// set facet constraints
		tetin.facetconstraintlist[(facetConstrOffset+i)*2] = double(facetConstrOffset+i+1);
		tetin.facetconstraintlist[(facetConstrOffset+i)*2+1] = facetConstr[i];
	}

	// Set connectivity through facet and polygon
	countI = meshgeom.surfs.size();
	if(nFacetConstr != 1 && nFacetConstr!=countI && nFacetConstr != 0){
		std::cerr << "Warning: Number of constraints is not 1 or the "
			"number of surfaces, the code will work but the behaviour "
			"is suprising." << std::endl;
	}
	nFacetConstr = max(nFacetConstr,1);
	for (int i = 0; i < countI; ++i){

		tetin.allocatefacet(i+facetOffset,1);

		meshgeom.surfs(i)->OrderedVerts(&meshgeom,orderVert);
		tetin.allocatefacetpolygon(i+facetOffset,0,orderVert.size());
		countK = orderVert.size();
		for (int k = 0; k < countK; ++k){
			tetin.facetlist[i+facetOffset].polygonlist[0].vertexlist[k] 
				= meshgeom.verts.find(orderVert[k])+pointOffset;
		}
		// Set the facet marker list to match the mod of the facets added
		tetin.facetmarkerlist[i+facetOffset]=facetConstrOffset+(i%nFacetConstr)+1;
		//pointMarker;//
	}
}

void tetgen::internal::Mesh2Tetgenio(const mesh &meshgeom, const mesh &meshdomain,
	tetgen::io_safe &tetin, int numHoles){
	/*
	Converts a mesh into the safe allocation tetgenio format.
	
	Rules of conversion are: 
		- surfaces become a facet
		- edges are not logged
		- points come as a list
	*/
	

	tetin.firstnumber=0;

	tetin.numberoffacets = meshgeom.surfs.size() + meshdomain.surfs.size();
	tetin.numberofholes = numHoles;
	tetin.numberofpoints += meshgeom.verts.size() + meshdomain.verts.size();

	tetin.numberofpointmtrs = 1;
	tetin.numberoffacetconstraints = 2;


	tetin.allocate();

	tetgen::internal::MeshData2Tetgenio(meshgeom, tetin, 0, 0, 1, {0.0},{0.0},0);
	tetgen::internal::MeshData2Tetgenio(meshdomain, tetin, meshgeom.surfs.size(), 
		meshgeom.verts.size(), -1, {-1.0},{0.0},1);
}

void tetgen::internal::Mesh2TetgenioPoints(const mesh &meshgeom, const mesh &meshdomain,
	tetgen::io_safe &tetin){
	/*
	Converts a mesh into the safe allocation tetgenio format.
	
	Rules of conversion are: 
		- volu becomes a facet
		- surf becomes a polygon in the facet
		- edges are not logged
		- points come as a list
	*/
	
	tetin.firstnumber=0;

	tetin.numberoffacets = 0;
	tetin.numberofholes = 0;
	tetin.numberofpoints += meshgeom.verts.size() + meshdomain.verts.size();

	tetin.numberofpointmtrs = 0;
	tetin.numberoffacetconstraints = 0;


	tetin.allocate();
	int nPtsGeom = meshgeom.verts.size();
	for (int i = 0; i < nPtsGeom; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			tetin.pointlist[i*3+j]=meshgeom.verts(i)->coord[j];
		}
	}
	int nPtsDom = meshdomain.verts.size();
	for (int i = 0; i < nPtsDom; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			tetin.pointlist[(i+nPtsGeom)*3+j]=meshdomain.verts(i)->coord[j];
		}
	}
}

void tetgen::internal::PointCurvature2Metric(std::vector<double> &vertCurvature, 
	const tetgen::apiparam &inparam){

	double maxCurv = *max_element(vertCurvature.begin(), vertCurvature.end());
	double minCurv = 0.0;
	double deltaCurv = maxCurv-minCurv;
	auto lininterp = [&](double curv) double {
		return (inparam.surfedgelengths[1]-inparam.surfedgelengths[0])
			/(deltaCurv)*(curv-minCurv)+inparam.surfedgelengths[0];
		};

	int count = vertCurvature.size();
	for (int i = 0; i < count; ++i)
	{
		vertCurvature[i] = lininterp(vertCurvature[i]);
	}
}

void tetgen::input::RSVS2CFD(const snake &snakein, tetgen::io_safe &tetin,
	const tetgen::apiparam &tetgenParam){
	/*
	Processes a snake object into a tetgen input object.

	This function processes snake objects into the safe tetgenio object.
	This can then used to generate a tetrahedral mesh or to save it as the
	input file.
	*/

	mesh meshdomain, meshgeom;
	triangulation triRSVS;
	int nHoles;
	// int nPtsHole, kk, count;

	std::vector<int>  holeIndices;
	std::array<double, 3> upperB, lowerB;
	std::vector<int> vertPerSubDomain;
	std::vector<double> holeCoords;

	PrepareSnakeForCFD(snakein, tetgenParam.distanceTol, meshgeom, holeCoords);

	meshgeom.ReturnBoundingBox(lowerB, upperB);
	double distBound = tetgenParam.edgelengths.size()>0?
		tetgenParam.edgelengths[0]: tetgenParam.distanceTol;
	for (int i = 0; i < 3; ++i)
	{
		lowerB[i] -= distBound;
		upperB[i] += distBound;
	}
	meshdomain = BuildCutCellDomain(tetgenParam.lowerB, tetgenParam.upperB,
		lowerB, upperB, tetgenParam.edgelengths.size(), vertPerSubDomain);
	nHoles=holeCoords.size()/3;
	tetgen::internal::Mesh2Tetgenio(meshgeom, meshdomain, tetin, nHoles);

	// std::cout<< std::endl << "Number of holes " << nHoles << std::endl;

	// Assign the holes

	int count = holeCoords.size();
	for (int i = 0; i < count; ++i){
		tetin.holelist[i] = holeCoords[i];	
	}

	// Assign domain point metrics
	int startPnt = meshgeom.verts.size();
	if(tetgenParam.surfedgelengths[0]==tetgenParam.surfedgelengths[1]){
		tetin.SpecifyTetPointMetric(0,  startPnt,
			{tetgenParam.surfedgelengths[0]});
	} else {
		auto vertCurvature = CalculateVertexCurvature(meshgeom, 
			tetgenParam.curvatureSmoothing);
		DisplayVectorStatistics(vertCurvature);
		tetgen::internal::PointCurvature2Metric(vertCurvature, tetgenParam);
		tetin.SpecifyTetPointMetric(0,  startPnt,
			vertCurvature);
		DisplayVectorStatistics(vertCurvature);

	}
	for (int i = 0; i < int(tetgenParam.edgelengths.size()); ++i){
		tetin.SpecifyTetPointMetric(startPnt,  vertPerSubDomain[i],
			{tetgenParam.edgelengths[i]});
		startPnt = startPnt + vertPerSubDomain[i];
	}
	// Remove facet markers from the internal domain faces
	startPnt = meshgeom.surfs.size();
	int nSurfPerSubDomain = meshdomain.surfs.size()/vertPerSubDomain.size();
	for (int i = 0; i < int(vertPerSubDomain.size())-1; ++i){
		tetin.SpecifyTetFacetMetric(startPnt, nSurfPerSubDomain ,0);
		startPnt = startPnt + nSurfPerSubDomain;
	}
}

void tetgen::input::RSVSGRIDS(const mesh &meshdomain, const mesh &meshboundary,
	tetgen::io_safe &tetin, const tetgen::apiparam &tetgenParam){
	/*
	Processes a snake object into a tetgen input object.

	This function processes snake objects into the safe tetgenio object.
	This can then used to generate a tetrahedral mesh or to save it as the
	input file.
	*/

	triangulation triRSVS;
	int nHoles=0;
	// int nPtsHole, kk, count;

	std::vector<int>  holeIndices;
	std::vector<int> vertPerSubDomain;

	vertPerSubDomain.push_back(meshdomain.verts.size());
	

	tetgen::internal::Mesh2Tetgenio(meshboundary, meshdomain, tetin, nHoles);

	int startPnt = meshboundary.verts.size();
	// Assign "boundary" point metrics
	auto vertEdgeLength = CalculateVertexMinEdgeLength(meshboundary);
	int count = vertEdgeLength.size();
	for (int i = 0; i < count; ++i)
	{
		vertEdgeLength[i] = vertEdgeLength[i] * tetgenParam.surfedgelengths[0];
	}
	tetin.SpecifyTetPointMetric(0, startPnt, vertEdgeLength);
	// Assign domain point metrics
	vertEdgeLength = CalculateVertexMinEdgeLength(meshdomain);
	count = vertEdgeLength.size();
	for (int i = 0; i < count; ++i)
	{
		vertEdgeLength[i] = vertEdgeLength[i] * tetgenParam.surfedgelengths[0];
	}
	tetin.SpecifyTetPointMetric(startPnt, startPnt+vertEdgeLength.size(),
		vertEdgeLength);
}

void tetgen::input::RSVSGRIDS(const mesh &meshdomain,
	tetgen::io_safe &tetin, const tetgen::apiparam &tetgenParam){
	mesh meshboundary;
	meshboundary.Init(0,0,0,0);
	meshboundary.PrepareForUse();

	tetgen::input::RSVSGRIDS(meshdomain, meshboundary, 
		tetin, tetgenParam);
}

void tetgen::input::POINTGRIDS(const mesh &meshdomain, tetgen::io_safe &tetin,
	const tetgen::apiparam &tetgenParam, bool generateVoroBound){
	/*
	Processes a snake object into a tetgen input object.

	This function processes snake objects into the safe tetgenio object.
	This can then used to generate a tetrahedral mesh or to save it as the
	input file.
	*/

	mesh meshgeom;
	triangulation triRSVS;
	// int nPtsHole, kk, count;
	if(generateVoroBound){
		// Generate a set of points which fulfill the following:
		// - off the boundaries by dist tol using the "pinch" technique:
		// 	not on the corners or edges but on the faces that need
		// 	to be closed
		// 	For each face place 1 point in the middle.
		// 	For each edge place 2 points in the middle.
		// 	For each face place 3 points on the outside.
		int nVertsAdded = 0;
		std::vector<double> vertsToAdd;
		std::array<double, 3> vertModif = {0,0,0}, edgeModif = {0,0,0};
		#ifdef RSVS_DIAGNOSTIC_RESOLVED
		tecplotfile  tecout;
		tecout.OpenFile("../TESTOUT/rsvs_voro.plt", "a");
		#endif
		vertsToAdd.reserve(100);
		// add the points associated with the vertices
		for (int l = 0; l < 4; ++l)
		{
			vertModif = {0,0,0};
			if (l<3){
				vertModif[l] = tetgenParam.distanceTol;
			}
			for (int i = 0; i < 8; ++i)
			{	
				for (int j = 0; j < 3; ++j)
				{
					int k = pow(2,j);
					double delta = (tetgenParam.upperB[j]-tetgenParam.lowerB[j])*vertModif[j];
					vertsToAdd.push_back(
						(((i/k)%2) ? 
							tetgenParam.upperB[j]+delta 
							: tetgenParam.lowerB[j]-delta)
						);
				}
				nVertsAdded++;
			}
		}

		// add the points associated with the edges
		int nEdgeSteps = ceil(1.0/(tetgenParam.distanceTol*2));
		for (int ll = 1; ll < nEdgeSteps; ++ll)
		{
			for (int l = 0; l < 3; ++l)
			{ // for each direction of the edge (coord = 0.5)
				edgeModif = {0,0,0};
				edgeModif[l] = 1.0/double(nEdgeSteps)*ll;
				for (int j = 0; j < 2; ++j)
				{ // for each displacement off the edge
					vertModif = {0,0,0};
					vertModif[(l+j+1)%3] = tetgenParam.distanceTol;
					for (int i = 0; i < 4; ++i)
					{ // each of the edge that make a square
						
						for (int m = 0; m < 3; ++m)
						{
							if(m==l){ // this is the inactive dimension
								vertsToAdd.push_back(edgeModif[l]);
							} else {
								double delta = (tetgenParam.upperB[m]
									-tetgenParam.lowerB[m])*vertModif[m];
								int k = pow(2, (m+3-(l+1))%3);
								vertsToAdd.push_back(
									(((i/k)%2) ? 
								tetgenParam.upperB[m]+delta 
								: tetgenParam.lowerB[m]-delta)
							);
							}
						}
						nVertsAdded++;
					}
					
				}
			}
		}
		// add the points associated with the faces
		// add the points into a mesh

		meshgeom = Points2Mesh(vertsToAdd,3);
		#ifdef RSVS_DIAGNOSTIC_RESOLVED
		tecout.PrintMesh(meshgeom, 0,0,4);
		#endif
	} else{ 
		meshgeom.Init(0, 0, 0, 0);
	}
	meshgeom.PrepareForUse();
	
	tetgen::internal::Mesh2TetgenioPoints(meshgeom, meshdomain, tetin);
}
/**
Ouputs a tetgen io object to an SU2 mesh file.

File format: https://su2code.github.io/docs/Mesh-File/

@param      fileName  is a string with the path to the target mesh fileName
@param      tetout    is tetgenio object (safe)

@throws     invalid_argument  if `fileName` cannot be opened.

*/
void tetgen::output::SU2(const char* fileName, const tetgenio &tetout){
	
	std::ofstream meshFile;
	std::vector<int> diffBounds, numBounds;
	int DEINCR;

	meshFile.open(fileName);
	// De-increments indices if tetgen indices start from 1 to match SU2
	// which is 0 indexed
	DEINCR = tetout.firstnumber ? -1 : 0;

	// If file was not opened correctly
	if(meshFile.fail()){
		std::cerr << std::endl << "Error: " << fileName << std::endl <<
			 " in could not be opened" << std::endl; 
		RSVS3D_ERROR_ARGUMENT("Output mesh file could not be opened");
	}
	// Set appropriate precision for mesh point position
	meshFile.precision(16);
	// Dimensions
	meshFile << "NDIME= 3" << std::endl;
	// Points
	meshFile << "NPOIN= " << tetout.numberofpoints << std::endl;
	for (int i = 0; i < tetout.numberofpoints; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			meshFile << tetout.pointlist[i*3+j] << " ";
		}
		meshFile << std::endl;
	}
	// Elements (edges, triangles and tets)
	int nElms = tetout.numberoftetrahedra;
	meshFile << "NELEM= " << nElms << std::endl;

	for (int i = 0; i < tetout.numberoftetrahedra; ++i)
	{ // Tetrahedra
		meshFile << "10 " ;
		for (int j = 0; j < 4; ++j)
		{
			meshFile << tetout.tetrahedronlist[i*tetout.numberofcorners+j]+DEINCR<< " ";
		}
		meshFile << std::endl;
	}
	// Boundaries
	// 1 - Process number of different boundaries
	diffBounds.reserve(10);
	numBounds.reserve(10);
	diffBounds.push_back(0);
	numBounds.push_back(0);
	int nBounds = diffBounds.size();
	// Build an array of different boundary numbers.
	// 0 is first as it is likely the most common
	for (int i = 0; i < tetout.numberoftrifaces; ++i)
	{	
		bool flagSame=false;
		int j = 0;
		while ((j < nBounds) && !flagSame){
			flagSame = flagSame || (diffBounds[j]== tetout.trifacemarkerlist[i]);
			++j;
		}
		if(!flagSame){
			diffBounds.push_back(tetout.trifacemarkerlist[i]);
			numBounds.push_back(0);
			++nBounds;
			++j;
		}
		numBounds[j-1]++;
	}
	// DisplayVector(diffBounds);
	// DisplayVector(numBounds);
	// 2 - write boundaries (skipping 0 boundary)
	meshFile << "NMARK= " << nBounds-1 << std::endl;
	for (int k = 1; k < nBounds; ++k)
	{
		meshFile << "MARKER_TAG= " << diffBounds[k] << std::endl;
		meshFile << "MARKER_ELEMS= " << numBounds[k] << std::endl;
		for (int i = 0; i < tetout.numberoftrifaces; ++i)
		{ // Triangular faces
			if(tetout.trifacemarkerlist[i]==diffBounds[k]){
				meshFile << "5 " ;
				for (int j = 0; j < 3; ++j)
				{
					meshFile << tetout.trifacelist[i*3+j]+DEINCR<< " ";
				}
				meshFile << std::endl;
			}
		}
	}

	// Edges not needed
	// for (int i = 0; i < tetout.numberofedges; ++i)
	// { // Edges
	// 	meshFile << "3 " << tetout.edgelist[i*2]<< " "
	// 		<< tetout.edgelist[i*2+1] << " " << std::endl;
	// }
}

void tetgen::internal::CloseVoronoiMesh(mesh &meshout, tetgen::io_safe &tetout, 
	std::vector<int> &rayEdges, int DEINCR, tetgen::dombounds boundBox){
	/*
	Mesh is still missing some elements due to Voronoi diagrams not being
	naturally closed:
	1) Terminate rays with vertices
	2.a) Link ray-vertices with edges 
	2.b) Include edges in equivalent faces
	3.a) Close volumes by creating faces from edges  
	4) Triangulate the surface mesh to make sure it is flat
	Implementation notes:
		will need hashed vectors of the edges added

	Issue:
		A voronoi is not bounded -> need another algorithm there
		(idea -> points lying out of convex hull removed and the 
		tetout vertex becomes a ray)
	*/

	int nVerts, nEdges, nSurfs;
	int count,  n;
	double lStep, lStepMin=-0.00; 
	vert vertNew;
	edge edgeNew;
	surf surfNew;
	HashedVector<int, int> rays, newVerts, newEdges, newSurfs, rayFacesSearch;
	std::vector<int> rayCells; // nonclosed voronoi faces
	std::vector<int> &rayFaces = rayFacesSearch.vec; // nonclosed voronoi faces
	std::vector<int> pair;
	

	nVerts = meshout.verts.size();
	nEdges = meshout.edges.size();
	nSurfs = meshout.surfs.size();

	rays.vec = rayEdges;
	rays.GenerateHash();
	// Terminating rays:
	// Use a fixed length and assign connectivity to the ray
	count = rayEdges.size();
	vertNew.edgeind.reserve(6);
	vertNew.edgeind.assign(1,0);
	for (int i = 0; i < count; ++i){
		vertNew.index = nVerts+1;
		vertNew.edgeind[0] = rayEdges[i];
		meshout.edges[rayEdges[i]-1].vertind[1] = nVerts+1;
		lStep = tetgen::internal::ProjectRay(3, boundBox, 
			tetout.vedgelist[rayEdges[i]-1].vnormal, 
			meshout.verts[tetout.vedgelist[rayEdges[i]-1].v1+DEINCR].coord,
			lStepMin);
		// std::cout << lStep << " " << std::endl;
		for (int j = 0; j < 3; ++j){
			vertNew.coord[j]=lStep*tetout.vedgelist[rayEdges[i]-1].vnormal[j]
				+ meshout.verts[tetout.vedgelist[rayEdges[i]-1].v1+DEINCR].coord[j];
		}
		meshout.verts.push_back(vertNew);
		newVerts.vec.push_back(nVerts+1);
		nVerts++;
	}

	newVerts.GenerateHash();
	
	// For each face 
	// check if it is linked to rays, if it is
	// an edge will be established between the two ray vertices
	meshout.edges.HashArray();
	rayFaces = ConcatenateVectorField(meshout.edges, &edge::surfind, 
		meshout.edges.find_list(rays.vec));
	sort(rayFaces);
	unique(rayFaces);
	count = rayFaces.size();
	edgeNew.surfind.reserve(6);
	edgeNew.surfind.assign(1,0);
	for (int i = 0; i < count; ++i){
		// identify in each open face the 2 rays
		pair.clear();
		for (int a : rays.find_list(meshout.surfs[rayFaces[i]-1].edgeind)){
			if (a!=-1){
				pair.push_back(a);
			}
		}
		if (pair.size()!=2){RSVS3D_ERROR_LOGIC("Voronoi open face should have 2 rays");}
		// Connect the 2 rays and their edges
		edgeNew.index = nEdges+1;
		for (int j = 0; j < 2; ++j){
			edgeNew.vertind[j] = meshout.edges[rays.vec[pair[j]]-1].vertind[1];
			meshout.verts[edgeNew.vertind[j]-1].edgeind.push_back(edgeNew.index);
		}
		edgeNew.surfind[0] = rayFaces[i];

		meshout.edges.push_back(edgeNew);
		meshout.surfs[rayFaces[i]-1].edgeind.push_back(edgeNew.index);
		newEdges.vec.push_back(edgeNew.index);
		nEdges++;
	}

	newEdges.GenerateHash();
	rayFacesSearch.GenerateHash();
	// For each volume a surface is established closing 
	// the volume using all the face indices and looking 
	// for those which have been modified.

	meshout.surfs.HashArray();
	rayCells = ConcatenateVectorField(meshout.surfs, &surf::voluind, 
		meshout.surfs.find_list(rayFaces));
	sort(rayCells);
	unique(rayCells);
	count = rayCells.size();

	for (int i = 0; i < count; ++i){
		// identify in each open cell the faces that were open
		pair.clear();
		surfNew.edgeind.clear();
		for (auto a : rayFacesSearch.find_list(
				meshout.volus[rayCells[i]-1].surfind)){
			if (a!=-1){
				pair.push_back(a);
			}
		}
		if(pair.size()<3){
			RSVS3D_ERROR_LOGIC("At least 3 edges are required to close the face");
		}
		// identify the edge in each face which is new and assign correxponding
		// connectivities to surfNew;
		surfNew.index = nSurfs+1;
		for (auto ind : pair) {
			auto &edgeind = meshout.surfs[rayFacesSearch.vec[ind]-1].edgeind;
			n = edgeind.size();
			int j, k=0;
			for (j = 0; j < n; ++j){
				if(newEdges.find(edgeind[j])!=-1){
					surfNew.edgeind.push_back(edgeind[j]);
					meshout.edges[edgeind[j]-1].surfind.push_back(surfNew.index);
					// break;
					k++;
				}
			}
			// if(j>=n)
			if(k==0){
				RSVS3D_ERROR_LOGIC("None of the edges of the open face"
				" were recognised as new. This should not happen.");
			}
			if (k>1){
				std::cerr << "Unexpected behaviour with number of edges" 
					<< std::endl;
			}
		}
		surfNew.voluind[0] = rayCells[i];
		surfNew.voluind[1] = 0;

		meshout.volus[rayCells[i]-1].surfind.push_back(surfNew.index);
		meshout.surfs.push_back(surfNew);
		nSurfs++;
	}
}

tetgen::dombounds tetgen::output::GetBoundBox(tetgen::io_safe &tetout){

	tetgen::dombounds domBounds;
	domBounds[0] = {INFINITY,INFINITY,INFINITY};
	domBounds[1] = {-INFINITY,-INFINITY,-INFINITY};
	int nVerts = tetout.numberofpoints;
	int count = nVerts;
	int n=3;
	for (int i = 0; i < count; ++i){
		for (int j = 0; j < n; ++j){
			domBounds[0][j] = domBounds[0][j] <= tetout.pointlist[i*n+j] ?
				domBounds[0][j] : tetout.pointlist[i*n+j];
			domBounds[1][j] = domBounds[1][j] >= tetout.pointlist[i*n+j] ?
				domBounds[1][j] : tetout.pointlist[i*n+j];
		}
	}

	return domBounds;
}

mesh tetgen::output::VORO2MESH(tetgen::io_safe &tetout){
	/*
	Translates a tetgen output to the RSVS native mesh format 
	*/
	mesh meshout;

	int nVerts, nEdges, nSurfs, nVolus;
	int count,INCR, DEINCR,n;
	vector<int> rayInd, rayInd2;
	tetgen::dombounds boundBox=tetgen::output::GetBoundBox(tetout);

	nVerts = tetout.numberofvpoints;
	nEdges = tetout.numberofvedges;
	nSurfs = tetout.numberofvfacets;
	nVolus = tetout.numberofvcells;

	// INCR and DEINCR convert from tetgen array position to index and
	// index to array position respectively
	INCR = tetout.firstnumber ? 0 : 1;
	DEINCR = tetout.firstnumber ? -1 : 0;

	// if `tetrahedralize` has not been run with the correct flags connectivity is
	// unavailable, throw an error. (tests for flag -v)
	if(tetout.vpointlist == NULL){
		RSVS3D_ERROR_ARGUMENT("TET2MESH requires flag -v be passed to tetgen");
	}

	meshout.Init(nVerts, nEdges, nSurfs, nVolus);
	meshout.PopulateIndices();

	// These need to be cleared so that push_back works as expected
	for (int i = 0; i < nEdges; ++i){meshout.edges[i].vertind.clear();}
	rayInd.reserve(nEdges);

	// Assign coordinates
	count = nVerts;
	n=3;
	for (int i = 0; i < count; ++i){
		meshout.verts[i].index = i+1;
		for (int j = 0; j < n; ++j){
			meshout.verts[i].coord[j] = tetout.vpointlist[i*n+j];
		}
	}

	// assign edge-vertex connectivity, relies on the implied ordering and 
	// indexing of the output tetgenio object an of `PopulateIndices()`
	count = nEdges;
	for (int i = 0; i < count; ++i)	
	{ // Assign edge to vertex connectivity
		meshout.edges[i].vertind = {
			tetout.vedgelist[i].v1+INCR, tetout.vedgelist[i].v2+INCR
		};
		// Assign equivalent vertex to edge connectivity
		meshout.verts[tetout.vedgelist[i].v1+DEINCR].edgeind.push_back(i+1);

		if(tetout.vedgelist[i].v2>-1){
			// if it is not a ray add it to vertex v2
			meshout.verts[tetout.vedgelist[i].v2+DEINCR].edgeind.push_back(i+1);
		} else {
			// if it is a ray save it to the list of rays
			rayInd.push_back(i+1);
		}
	}

	count = nSurfs;
	int nRays = 0;
	for (int i = 0; i < count; ++i)	{ // loop through surfs
		// Surfaces to edges 
		n = tetout.vfacetlist[i].elist[0];
		for (int j = 0; j < n; ++j)	{ // loop through corresponding connections
			if(tetout.vfacetlist[i].elist[j+1]>-1){
				meshout.surfs[i].edgeind.push_back(
					tetout.vfacetlist[i].elist[j+1]+INCR);
				meshout.edges[tetout.vfacetlist[i].elist[j+1]+DEINCR]
					.surfind.push_back(i+1);
			} else {
				nRays++;
			}
		}
		//  surface to volumes
		meshout.surfs[i].voluind[0]=tetout.vfacetlist[i].c1+INCR;
		meshout.surfs[i].voluind[1]=tetout.vfacetlist[i].c2+INCR;
		#ifdef SAFE_ALGO
		if(tetout.vfacetlist[i].c1==tetout.vfacetlist[i].c2){
			RSVS3D_ERROR_LOGIC("Tetgen outputs returns a face connected to "
				"the same two cells");
		}
		#endif //SAFE_ALGO
		n=2;
		for (int j = 0; j < n; ++j)	{ // loop through corresponding connections
			meshout.volus[meshout.surfs[i].voluind[j]-1].surfind.push_back(i+1);
		}
	}
	
	// Close the mesh
	tetgen::internal::CloseVoronoiMesh(meshout, tetout, rayInd, DEINCR, boundBox);

	// Prepare mesh for use
	for (int i = 0; i < nVolus; ++i){
		meshout.volus[i].target=(double(rand()%1000001)/1000000);
		meshout.volus[i].fill = meshout.volus[i].target;
		meshout.volus[i].error = meshout.volus[i].target;
	}
	meshout.HashArray();
	meshout.TestConnectivityBiDir(__PRETTY_FUNCTION__);
	meshout.TightenConnectivity();
	meshout.PrepareForUse();
	// meshout.displight();
	
	auto groupedVertices = GroupCloseVertices(meshout, 1e-7);
	auto out = meshout.MergeGroupedVertices(groupedVertices,false);
	sort(out);
	unique(out);
	meshout.RemoveSingularConnectors(out);
	meshout.HashArray();
	meshout.PrepareForUse();
	meshout.TightenConnectivity();
	meshout.OrderEdges();
	meshout.SetBorders();

	#ifdef SAFE_ALGO
	if(meshout.TestConnectivityBiDir(__PRETTY_FUNCTION__)){
		RSVS3D_ERROR_LOGIC("Errors in connectivity detected.");
	}
	#endif

	FlattenBoundaryFaces(meshout);
	meshout.PrepareForUse();
	meshout.TightenConnectivity();
	meshout.OrderEdges();
	meshout.SetBorders();
	
	#ifdef SAFE_ALGO
	if(meshout.TestConnectivityBiDir(__PRETTY_FUNCTION__)){
		RSVS3D_ERROR_LOGIC("Errors in connectivity detected.");
	}
	#endif

	return meshout;
}

/**
 * @brief      Genrates Snaking and VOS RSVS meshes from the voronoi diagram
 * of a set of points
 *
 * @param[in]  vecPts    a vector of input points (3 coordinate) followed
 *  by a target volume fraction. Vecpts is a 1D vector with 4 values per point.
 * 
 * @param      vosMesh   The vos mesh
 * @param      snakMesh  The snaking mesh
 * @param      inparam   The tetgen interface parameter at input.
 *
 * @return     Returns the mapping of the original points to the snake 
 * mesh volumes.
 */
std::vector<int> tetgen::RSVSVoronoiMesh(const std::vector<double> &vecPts, 
	mesh &vosMesh, mesh &snakMesh,
	tetgen::apiparam &inparam){
	/*
	Generates a VOS and snaking grid from a set of points
	
	Args:
		vecPts: A 1 dimensional vector containing point coordinates and 
		a point metric which will be used as the fill of the vosMesh.
		vosMesh: a reference to an empty mesh object in which the VOS mesh
		will be stored.
		snakMesh: a reference to an empty mesh object in which a snaking mesh
		will be stored.

	Steps:
		1. Voronoi the set of points
		2. Tetrahedralize the voronoi diagram
		3. Use flooding to establish the relationship between voronoi cells
		and new cells
		4. Establish the containments of the original set of points matching
		the additional data to the cells.
		5. Prep meshes, as parent child
		6. Use that info to establish containments of the original points.

	Step 3 - detail:
		1. Mark each point in the original facets with -1.
		2. Flood all points in the tetrahedralisation with the same number
		not allowing a point with a -1 to be processed.
		3. Transfer these point blocks to the tetrahedrons.
		4. Remove all points which are marked with block numbers to generate the
		parent mesh.

	*/
	mesh voroMesh;
	std::vector<int> vertsVoro;
	std::vector<int> vertBlock, elmMapping;


	// Step 1 - 2 
	auto boundFaces = tetgen::voronoi::Points2VoroAndTetmesh(vecPts,voroMesh,
		snakMesh, inparam);


	int nBlocks = snakMesh.ConnectedVolumes(elmMapping, boundFaces);

	#ifdef SAFE_ALGO
	if(nBlocks!=voroMesh.volus.size()){
		std::cerr << "Error : nBlocks (" << nBlocks 
			<< ")!=nVolus Voromesh (" << voroMesh.volus.size() <<")" 
			<< std::endl;
		RSVS3D_ERROR_LOGIC("Voromesh and flooding do not match");
	}
	#endif //SAFE_ALGO
	auto elmMappingCopy =elmMapping;
	for (int i = 0; i < nBlocks; ++i)
	{
		int indRep=-1;
		int jStart = -1;
		do {
			jStart++;
			if(elmMappingCopy[jStart]==i+1){
				indRep = snakMesh.volus(jStart)->index;
			}
		} while(indRep==-1);
		for (int j = jStart; j < int(elmMapping.size()); ++j)
		{
			if(elmMappingCopy[j]==i+1){
				elmMapping[j]=indRep;
			}
		}
	}
	CoarsenMesh(snakMesh,vosMesh,elmMapping);
	snakMesh.AddParent(&vosMesh,elmMapping);
	snakMesh.PrepareForUse();
	snakMesh.OrientFaces();
	vosMesh.PrepareForUse();
	#ifdef SAFE_ALGO
	if(vosMesh.volus.size()!=voroMesh.volus.size()){
		std::cerr << "Error : vosMesh (" << vosMesh.volus.size() 
			<< ")!=nVolus Voromesh (" << voroMesh.volus.size() <<")" 
			<< std::endl;
		RSVS3D_ERROR_LOGIC("Voromesh and vosMesh do not match");
	}
	if(vosMesh.TestConnectivityBiDir(__PRETTY_FUNCTION__)){
		RSVS3D_ERROR_LOGIC("Invalid VOS mesh generated.");
	};
	#endif //SAFE_ALGO
	// Step 6 - Check original points containments
	std::vector<int> vecPtsMapping = snakMesh.VertexInVolume(vecPts, 4);
	#ifdef RSVS_DIAGNOSTIC_RESOLVED
	DisplayVector(vecPtsMapping);
	cout << endl;
	#endif //RSVS_DIAGNOSTIC_RESOLVED
	int count = vosMesh.volus.size();
	for (int i = 0; i < count; ++i)
	{
		vosMesh.volus[i].fill = 0.0;
	}
	int nVecMap = vecPtsMapping.size();
	for (int i = 0; i < nVecMap; ++i)
	{
		vosMesh.volus[
			vosMesh.volus.find(
				elmMapping[snakMesh.volus.find(vecPtsMapping[i])]
			, true) // turn off warning
			].fill=vecPts[i*4+3];
	}
	for (int i = 0; i < count; ++i)
	{
		vosMesh.volus[i].target = vosMesh.volus[i].fill;
		vosMesh.volus[i].error = 0.0;
	}
	vosMesh.PrepareForUse();
	vosMesh.OrientFaces();
	vosMesh.SetBorders();
	snakMesh.SetBorders();
	
	return vecPtsMapping;
}

/**
 * @brief      Genrates an SU2 mesh file from a snake.
 * 
 * Uses tetgen to generate a volume mesh around a snake and outputs it to
 * the SU2 format.
 * 
 * @param[in]  snakein   A snake which needs to be meshed
 * @param[in]  fileName  The file name
 * @param      inparam   tetgen interface parameter object. Used to define boundary
 *  growth rate and element sizes.
 */
void tetgen::SnakeToSU2(const snake &snakein, const std::string &fileName,
	tetgen::apiparam &inparam){
	tetgen::io_safe tetin, tetout;

	tetgen::input::RSVS2CFD(snakein, tetin, inparam);

	// Tetrahedralize the PLC. Switches are chosen to read a PLC (p),
	//   do quality mesh generation (q) with a specified quality bound
	//   (1.414), and apply a maximum volume constraint (a0.1).
	inparam.command = "pkqm"; 
	tetrahedralize(inparam.command.c_str(), &tetin, &tetout);

	tetgen::output::SU2(fileName.c_str(),tetout);
}

/**
 * @brief      Translates a tetgen output to the RSVS native mesh format.
 *
 * @param      tetout  the tetgenio object containing a mesh to be translated to
 *                     the native RSVS mesh format.
 *
 * @return     mesh object containing the translated grid.
 * 
 * @throw  invalid_argument if tetout was generated without passing the neighbour
 * flag to tetgen (`-nn`)
 */
mesh tetgen::output::TET2MESH(tetgen::io_safe &tetout){
	/*
	Translates a tetgen output to the RSVS native mesh format 
	*/
	mesh meshout;

	int nVerts, nEdges, nSurfs, nVolus;
	int count,INCR, DEINCR,n;

	nVerts = tetout.numberofpoints;
	nEdges = tetout.numberofedges;
	nSurfs = tetout.numberoftrifaces;
	nVolus = tetout.numberoftetrahedra;

	// INCR and DEINCR convert from tetgen array position to index and
	// index to array position respectively
	INCR = tetout.firstnumber ? 0 : 1;
	DEINCR = tetout.firstnumber ? -1 : 0;

	// if `tetrahedralize` has not been run with the correct flags connectivity
	// is unavailable, throw an error. (tests for flag -nn)
	if(tetout.tet2facelist == NULL){
		RSVS3D_ERROR_ARGUMENT("TET2MESH requires flag -nn be passed to tetgen");
	}

	meshout.Init(nVerts, nEdges, nSurfs, nVolus);
	meshout.PopulateIndices();

	// These need to be cleared so that push_back works as expected
	for (int i = 0; i < nSurfs; ++i){meshout.surfs[i].voluind.clear();}
	for (int i = 0; i < nEdges; ++i){meshout.edges[i].vertind.clear();}

	// Assign coordinates
	count = nVerts;
	n=3;
	for (int i = 0; i < count; ++i){
		for (int j = 0; j < n; ++j){
			meshout.verts[i].coord[j] = tetout.pointlist[i*n+j];
		}
	}

	// assign edge-vertex connectivity, relies on the implied ordering and 
	// indexing of the output tetgenio object an of `PopulateIndices()`
	count = nEdges;
	for (int i = 0; i < count; ++i)	
	{ // Assign edge to vertex connectivity
		meshout.edges[i].vertind = {
			tetout.edgelist[i*2]+INCR, tetout.edgelist[i*2+1]+INCR
		};
		for (int j = 0; j < 2; ++j)
		{ // Assign equivalent vertex to edge connectivity
			meshout.verts[tetout.edgelist[i*2+j]+DEINCR]
				.edgeind.push_back(i+1);
		}
	}

	// Volumes to surfaces
	count = nVolus;	n = 4;
	for (int i = 0; i < count; ++i)	{ // loop through volus
		for (int j = 0; j < n; ++j)	{ // loop through corresponding connections
			meshout.volus[i].surfind.push_back(tetout.tet2facelist[i*n+j]+INCR);
			meshout.surfs[tetout.tet2facelist[i*n+j]+DEINCR].voluind.push_back(i+1);
		}
	}
	for (int i = 0; i < nSurfs; ++i){
		if(meshout.surfs[i].voluind.size()==1){
			meshout.surfs[i].voluind.push_back(0);
		}
	}

	// Surfaces to edges
	count = nSurfs;	n = 3;
	for (int i = 0; i < count; ++i)	{ // loop through volus
		for (int j = 0; j < n; ++j)	{ // loop through corresponding connections
			meshout.surfs[i].edgeind.push_back(tetout.face2edgelist[i*n+j]+INCR);
			meshout.edges[tetout.face2edgelist[i*n+j]+DEINCR].surfind.push_back(i+1);
		}
	}

	// Prepare mesh for use
	for (int i = 0; i < nVolus; ++i){
		meshout.volus[i].target=(double(rand()%1000001)/1000000);
		meshout.volus[i].fill = meshout.volus[i].target;
		meshout.volus[i].error = meshout.volus[i].target;
	}
	meshout.TightenConnectivity();
	meshout.PrepareForUse();

	meshout.TestConnectivityBiDir(__PRETTY_FUNCTION__);
	// meshout.displight();

	return meshout;
}

/**
 * @brief      Generate points inside volume cells of a mesh
 *
 * @param[in]  meshin   The input mesh
 * @param[in]  nLevels  The number of layers of points in each cell.
 *
 * @return     A tetgen io object which can be passed directly to tetrahedralize
 * as the fourth input.
 */
tetgen::io_safe tetgen::voronoi::GenerateInternalPoints(const mesh &meshin, 
	int nLevels){

	tetgen::io_safe tetinPts;
	auto vecPts = VolumeInternalLayers(meshin, nLevels);
	mesh meshgeom = Points2Mesh(vecPts,3);
	mesh meshdomain;
	meshdomain.Init(0,0,0,0);
	meshdomain.PrepareForUse();
	meshgeom.PrepareForUse();
	
	tetgen::internal::Mesh2TetgenioPoints(meshgeom, meshdomain, tetinPts);

	return(tetinPts);
}
std::vector<bool> tetgen::voronoi::Points2VoroAndTetmesh(const std::vector<double> &vecPts,
	mesh &voroMesh, mesh &tetMesh, const tetgen::apiparam &inparam){
	/*
	Turns a set of points into the voronisation and tetrahedralisation.
	*/
	tetgen::io_safe tetinV, tetoutV, tetinT, tetoutT, tetinSupport;
	mesh ptsMesh = Points2Mesh(vecPts,4);
	std::string cmd;
	int INCR;
	std::vector<int> voroPts;
	std::vector<bool> tetSurfInVoro;
	std::vector<double> lowerB, upperB;

	tetgen::input::POINTGRIDS(ptsMesh, tetinV, inparam, true);
	cmd = "Qv"; 
	tetrahedralize(cmd.c_str(), &tetinV, &tetoutV);
	voroMesh = tetgen::output::VORO2MESH(tetoutV);
	// Crop mesh 
	lowerB.assign(3, 0.0);
	upperB.assign(3, 0.0);
	for (int i = 0; i < 3; ++i)
	{
		double delta = inparam.upperB[i]-inparam.lowerB[i];
		lowerB[i]=inparam.lowerB[i] - delta*inparam.distanceTol*0.999;
		upperB[i]=inparam.upperB[i] + delta*inparam.distanceTol*0.999;
	}
	// voroMesh.displight();
	voroMesh.TestConnectivityBiDir(__PRETTY_FUNCTION__);
	voroMesh.CropAtBoundary(lowerB, upperB);
	voroMesh.TightenConnectivity();
	voroMesh.OrderEdges();
	voroMesh.SetBorders();
	// voroMesh.displight();
	FlattenBoundaryFaces(voroMesh);
	voroMesh.PrepareForUse();
	voroMesh.TightenConnectivity(); 
	voroMesh.OrderEdges();
	voroMesh.SetBorders();
	// Input mesh
	
	std::vector<int> vertPerSubDomain;
	mesh meshdomain = BuildCutCellDomain({0,0,0}, {0,0,0}, // upper bounds are not used
		inparam.lowerB, inparam.upperB, 1, vertPerSubDomain);
	tetgen::input::RSVSGRIDS(voroMesh, tetinT, inparam);

	int nInternalLayers = 0;
	if(inparam.surfedgelengths.size()>1){
		nInternalLayers = round(inparam.surfedgelengths[1]);
	}
	tetinSupport=tetgen::voronoi::GenerateInternalPoints(voroMesh,nInternalLayers);

	cmd = "pqminnefO9/7"; 
	try{
		tetrahedralize(cmd.c_str(), &tetinT, &tetoutT, &tetinSupport);
	} catch (exception const& ex) {
		std::cerr <<  std::endl << "Input points: ";
		DisplayVector(vecPts);
		std::cerr << std::endl;

		tetinT.save_poly("dbg/tetgen_error");
		tetinT.save_nodes("dbg/tetgen_error");
		voroMesh.displight();
		throw ex;
	}
	tetMesh = tetgen::output::TET2MESH(tetoutT);

	INCR = tetoutT.firstnumber ? 0 : 1;

	for (int i = 0; i < tetoutT.numberoftrifaces; ++i)
	{
		if(tetoutT.trifacemarkerlist[i]!=0){
			for (int j = 0; j < 3; ++j)
			{
				voroPts.push_back(tetoutT.trifacelist[i*3+j]+INCR);
			}
		}
	}
	sort(voroPts);
	unique(voroPts);
	tetSurfInVoro.reserve(tetoutT.numberoftrifaces);

	for (int i = 0; i < tetoutT.numberoftrifaces; ++i)
	{
		tetSurfInVoro.push_back(tetoutT.trifacemarkerlist[i]!=0);
	}

	return(tetSurfInVoro);
}

std::vector<bool> tetgen::voronoi::BoundaryFacesFromPoints(const mesh &meshin,
	const std::vector<int> &boundaryPts){
	/*
	
	*/

	int nSurfs, nVerts;
	std::vector<bool> boundaryFaces, boundaryVertsLog;
	std::vector<int> tempVerts;

	nVerts = meshin.verts.size();
	nSurfs = meshin.surfs.size();

	boundaryVertsLog.assign(nVerts, false);
	boundaryFaces.assign(nSurfs, false);
	int count  = boundaryPts.size();
	for (int i = 0; i < count; ++i)
	{
		boundaryVertsLog[meshin.verts.find(boundaryPts[i])]=true;
	}
	for (int i = 0; i < nSurfs; ++i)
	{
		tempVerts.clear();
		meshin.surfs(i)->OrderedVerts(&meshin, tempVerts);
		bool isBound = true;
		for (auto j : tempVerts){
			isBound = isBound && boundaryVertsLog[meshin.verts.find(j)];
			if(!isBound){
				break;
			}
		}
		boundaryFaces[i] = isBound;
	}
	return boundaryFaces;
}

void tetgen::io_safe::allocate(){

	// Allocation of points and associated attributes
	if (this->numberofpoints>0){
		this->pointlist = new REAL[this->numberofpoints * 3];
		this->pointattributelist = new REAL[this->numberofpoints * this->numberofpointattributes];
		this->pointmtrlist = new REAL[this->numberofpoints * this->numberofpointmtrs];
		this->pointmarkerlist = new int[this->numberofpoints];
	}
	// Allocate tetrahedron
	if(this->numberoftetrahedra>0){
		this->tetrahedronlist = new int[this->numberoftetrahedra*this->numberofcorners];
		this->tetrahedronattributelist = new REAL[this->numberoftetrahedra];
		this->tetrahedronvolumelist = new REAL[this->numberoftetrahedra
			*this->numberoftetrahedronattributes];
		// this->neighborlist = new int[this->numberoftetrahedra*4]; //output only
	}
	// Allocation of facets
	if(this->numberoffacets>0){
		this->facetlist = new tetgenio::facet[this->numberoffacets];
		this->facetmarkerlist = new int[this->numberoffacets];
		for (int i = 0; i < this->numberoffacets; ++i){
			this->facetlist[i].numberofpolygons=0;
			this->facetlist[i].numberofholes=0;
		}
	}
	// Allocation of facet constraints
	if(this->numberoffacetconstraints>0){
		// [boundary marker, area constraint value]
		this->facetconstraintlist = new REAL[this->numberoffacetconstraints*2];
	}
	// Allocation of holes (a set of points)
	if(this->numberofholes>0){
		this->holelist = new REAL[this->numberofholes*3];
	}
	// Allocation of regions (a set of points with attributes)
	// X, Y, Z, region attribute at index, maximum volume at index
	if(this->numberofregions>0){
		this->regionlist = new REAL[this->numberofregions*5];
	}

	// Allocate constraint
	if(this->numberofsegmentconstraints>0){
		this->segmentconstraintlist = new REAL[this->numberofsegmentconstraints*3];
	}
	// Allocate triangles
	if(this->numberoftrifaces>0){
		this->trifacelist = new int[this->numberoftrifaces*3];
		this->trifacemarkerlist = new int[this->numberoftrifaces];
		this->o2facelist = new int[this->numberoftrifaces*3];
		this->face2tetlist = new int[this->numberoftrifaces*2];
	}

	// Allocate edges
	if(this->numberofedges>0){
		this->edgelist = new int[this->numberofedges*2];
		this->edgemarkerlist = new int[this->numberofedges];
		this->o2edgelist = new int[this->numberofedges];
		this->edge2tetlist = new int[this->numberofedges];
		this->face2edgelist = new int[this->numberofedges*3];
	}

	// Voronoi implementation
	if(this->numberofvedges || this->numberofvpoints 
		|| this->numberofvcells || this->numberofvfacets){
		std::cerr << "Warning : tetgen::io_safe::allocate() does not support "
			"Voronoi variables" << std::endl; 
	}
}

void tetgen::io_safe::allocatefacet(int fIndex){
	if (fIndex<this->numberoffacets){
		if(this->facetlist[fIndex].numberofpolygons>0){
			this->facetlist[fIndex].polygonlist 
				= new tetgenio::polygon[this->facetlist[fIndex].numberofpolygons];
		} else {
			this->facetlist[fIndex].polygonlist = (tetgenio::polygon *) NULL;
		}
		if(this->facetlist[fIndex].numberofholes>0){
			this->facetlist[fIndex].holelist 
				= new REAL[this->facetlist[fIndex].numberofholes];
		} else {
			this->facetlist[fIndex].holelist = (REAL *) NULL;
		}
	} else {
		std::cerr << "Error: Index passed to facet allocation out "
			"of range" << std::endl;
	}
}
void tetgen::io_safe::allocatefacetpolygon(int fIndex, int pIndex){
	if (fIndex<this->numberoffacets){
		if (pIndex<this->facetlist[fIndex].numberofpolygons){
			if(this->facetlist[fIndex].polygonlist[pIndex].numberofvertices>0){
				this->facetlist[fIndex].polygonlist[pIndex].vertexlist 
					= new int[this->facetlist[fIndex].polygonlist[pIndex].numberofvertices];
			} else {
				this->facetlist[fIndex].polygonlist[pIndex].vertexlist =
					(int *) NULL;
			}

		} else {
			std::cerr << "Error: Index passed to polygon allocation out "
				"of range" << std::endl;
		}
	} else {
		std::cerr << "Error: Index passed to polygon allocation  "
			"for  facet out of range" << std::endl;
	}
}
void tetgen::io_safe::allocatefacet(int fIndex, int numPoly){
	/*
	Allocates facets specifying the number of polygons in the facet.
	*/
	if (fIndex<this->numberoffacets){
		this->facetlist[fIndex].numberofpolygons=numPoly;
		this->allocatefacet(fIndex);

	} else {
		std::cerr << "Error: Index passed to facet allocation out "
			"of range" << std::endl;
	}
}
void tetgen::io_safe::allocatefacetpolygon(int fIndex, int pIndex, int numVerts){
	if (fIndex<this->numberoffacets){
		if (pIndex<this->facetlist[fIndex].numberofpolygons){
			this->facetlist[fIndex].polygonlist[pIndex].numberofvertices = numVerts;
			this->allocatefacetpolygon(fIndex,pIndex);
		} else {
			std::cerr << "Error: Index passed to polygon allocation out "
				"of range" << std::endl;
		}
	} else {
		std::cerr << "Error: Index passed to polygon allocation  "
			"for  facet out of range" << std::endl;
	}
}
void tetgen::io_safe::SpecifyTetPointMetric(int startPnt, int numPnt, 
	const std::vector<double> &mtrs){
	/*
	assigns the mtr list from `startPnt` index
	*/

	int count = int(mtrs.size())<this->numberofpointmtrs ? 
		mtrs.size() : this->numberofpointmtrs;



	if((startPnt+numPnt)>this->numberofpoints){
		RSVS3D_ERROR_ARGUMENT("Metrics out of range in tetgen_io (too high)");
	} else if (startPnt<0){
		RSVS3D_ERROR_ARGUMENT("Metrics out of range in tetgen_io (too low)");
	}

	for (int i = startPnt; i < startPnt+numPnt; ++i)	{
		for (int j = 0; j < count; ++j)
		{
			this->pointmtrlist[(i*this->numberofpointmtrs)+j]
				=mtrs[j];
		}

	}
}

void tetgen::io_safe::SpecifyIndividualTetPointMetric(int startPnt, int numPnt, 
	const std::vector<double> &mtrs){
	/*
	assigns the mtr list from `startPnt` index
	*/

	if((startPnt+numPnt)>this->numberofpoints){
		RSVS3D_ERROR_ARGUMENT("Metrics out of range in tetgen_io (too high)");
	} else if (startPnt<0){
		RSVS3D_ERROR_ARGUMENT("Metrics out of range in tetgen_io (too low)");
	} else if (numPnt*this->numberofpointmtrs != int(mtrs.size())){
		RSVS3D_ERROR_ARGUMENT("Metrics need to be the same size as numPnt"
			"the number of point metrics");
	}
	int count = mtrs.size();
	for (int i = 0; i < count; ++i)	{
		this->pointmtrlist[(startPnt*this->numberofpointmtrs)+i]
			=mtrs[i];
	}
}
void tetgen::io_safe::SpecifyTetFacetMetric(int startPnt, int numPnt, 
	int marker){
	/*
	assigns the marker from `startPnt` index
	*/

	if((startPnt+numPnt)>this->numberoffacets){
		RSVS3D_ERROR_ARGUMENT("Metrics out of range in tetgen_io (too high)");
	} else if (startPnt<0){
		RSVS3D_ERROR_ARGUMENT("Metrics out of range in tetgen_io (too low)");
	}

	for (int i = startPnt; i < startPnt+numPnt; ++i)	{
		this->facetmarkerlist[i] = marker;
	}
}


//================================
// Tetgen to json
//================================

void tetgen::to_json(json& j, const tetgen::apiparam& p){
	j = json{
		{"lowerB", p.lowerB},
		{"upperB", p.upperB},
		{"surfedgelengths", p.surfedgelengths},
		{"curvatureSmoothing", p.curvatureSmoothing},
		{"edgelengths", p.edgelengths},
		{"distanceTol", p.distanceTol},
		{"generateMeshInside", p.generateMeshInside},
		{"command", p.command}
	};
}
void tetgen::from_json(const json& j, tetgen::apiparam& p){
	j.at("lowerB").get_to(p.lowerB);
	j.at("upperB").get_to(p.upperB);
	p.surfedgelengths=j.at("surfedgelengths").get<std::array<double,2>>();
	j.at("curvatureSmoothing").get_to(p.curvatureSmoothing);
	p.edgelengths=j.at("edgelengths").get<std::vector<double>>();
	j.at("distanceTol").get_to(p.distanceTol);
	j.at("generateMeshInside").get_to(p.generateMeshInside);
	p.command=j.at("command").get<std::string>();
}


void tetgen::apiparam::ReadJsonString(const std::string &jsonStr){

	try {
		json j=json::parse(jsonStr);
		json jparam = *this;

		try{
			j = j.unflatten();
		} catch (std::exception const& ex) {
			// if j is already not flat catch the exception and move on
			// TODO check the correct exception is being thrown (one day)
		}
		rsvsjson::flatupdate(jparam, j, false, false);

		*this = jparam.get<tetgen::apiparam>();
	} catch (std::exception const& ex) {
		RSVS3D_ERROR_NOTHROW((std::string("Unhandled error while parsing json"
			" string into tetgen parameters: \n")+jsonStr).c_str());
		throw ex;
	}
}

// Test code
void tetgen::test::LoadData(mesh &snakeMesh, mesh &voluMesh, 
	snake &snakein, mesh &triMesh){
	/*
	Loads data for tetgen testing
	*/

	triangulation triRSVS;

	snakeMesh.read("../TESTOUT/testtetgen/SnakeMesh_181205T193158_sphere2.msh");
	voluMesh.read("../TESTOUT/testtetgen/VoluMesh_181205T193158_sphere2.msh");
	snakein.read("../TESTOUT/testtetgen/Snake_181205T193158_sphere2.3snk");
	// snakein.snakeconn.read("../TESTOUT/testtetgen/SnakeConn_181205T193158_sphere2.msh");
	snakein.snakemesh = &snakeMesh;
	
	snakeMesh.PrepareForUse();
	// snakeMesh.displight();
	
	voluMesh.PrepareForUse();
	// voluMesh.displight();

	snakein.PrepareForUse();
	// snakein.displight();
	snakein.AssignInternalVerts();

	snakeMesh.AddParent(&voluMesh);

	triRSVS.stattri.clear();
	triRSVS.trivert.clear();
	triRSVS.PrepareForUse();
	TriangulateMesh(snakeMesh,triRSVS);
	TriangulateSnake(snakein,triRSVS);
	triRSVS.PrepareForUse();
	triRSVS.CalcTriVertPos();
	MaintainTriangulateSnake(triRSVS);

	MeshTriangulation(triMesh,snakein.snakeconn,triRSVS.dynatri, triRSVS.trivert);
	// triMesh.displight();
	triMesh.PrepareForUse();
	triMesh.OrderEdges();
	triMesh.TestConnectivityBiDir(__PRETTY_FUNCTION__);

	// triMesh is the triangulation in mesh format
	// it is not cleaned up for tiny surfaces
}

int tetgen::test::api(){
	int errCount = 0;
	{
		tetgen::apiparam p1, p2;
		json j1, j2;

		p1.distanceTol=1.0;


		j1 = p1;
		p2 = j1;
		j2 = p2;

		if(j1!=j2){
			RSVS3D_ERROR_NOTHROW(" Assignement operator not symmetric");
			errCount++;
		} else {
			std::cout << "Assignement passed " << std::endl;
		}
	}
	{
		tetgen::apiparam p1(
			R"foo(
			{
				"/distanceTol":2.0,
				"/edgelengths/2":2.0
			}
			)foo");
		tetgen::apiparam p2;
		json j1, j2;
		

		if(p1.distanceTol!=2.0){
			RSVS3D_ERROR_NOTHROW(" JSON Constructor not working for"
				" distanceTol");
			errCount++;
		} else {
			std::cout << "Constructor passed " << std::endl;
		}
		if(p1.edgelengths.size()!=3){
			RSVS3D_ERROR_NOTHROW(" JSON Constructor operator not "
				"working for vector beyond bounds");
			std::cerr << "edgelengths.size " << p1.edgelengths.size() 
				<< std::endl;
			errCount++;
		} else {
			std::cout << "Constructor passed " << std::endl;
		}
	}

	return(errCount);
}

int tetgen::test::CFD()
{
	/*CFD meshing process*/
	tetgen::io_safe tetin, tetout;

	tetgen::apiparam inparam;
	mesh snakeMesh, voluMesh, triMesh;
	snake snakein;

	try {

		inparam.lowerB= {-15.0, -15.0,-15.0};
		inparam.upperB= {15.0, 15.0, 15.0};
		inparam.distanceTol = 1e-3;
		inparam.edgelengths = {0.06,0.3,1.0,3.0};
		tetgen::test::LoadData(snakeMesh, voluMesh, snakein, triMesh);
		tetgen::input::RSVS2CFD(snakein, tetin, inparam);


		tetin.save_nodes("../TESTOUT/testtetgen/rsvs_3cell_2body");
		tetin.save_poly("../TESTOUT/testtetgen/rsvs_3cell_2body");
		

		// Tetrahedralize the PLC. Switches are chosen to read a PLC (p),
		//   do quality mesh generation (q) with a specified quality bound
		//   (1.414), and apply a maximum volume constraint (a0.1).
		inparam.command = "Qpkqm"; 
		tetrahedralize(inparam.command.c_str(), &tetin, &tetout);

		// tetout.save_nodes("../TESTOUT/testtetgen/rsvsout_3cell_2body");
		// tetout.save_elements("../TESTOUT/testtetgen/rsvsout_3cell_2body");
		// tetout.save_faces("../TESTOUT/testtetgen/rsvsout_3cell_2body");
		tetgen::output::SU2("../TESTOUT/testtetgen/rsvs_3cell_2body.su2",tetout);
		std::cout << "Finished the tettgen process " << std::endl;
	} catch (exception const& ex) {
		cerr << "Exception: " << ex.what() <<endl; 
		return 1;
	}
	return 0;
}

int tetgen::test::call()
{
	/*_RSVSgrid*/
	tetgen::io_safe tetin, tetin2,tetin3, tetout, tetout2, tetout3;
	tetgen::apiparam inparam;
	mesh meshdomain, meshtet, meshvoro, meshtet2;
	snake snakein;
	std::array<std::array<double, 2>, 3> dimDomain;
	tecplotfile tecout;


	try {
		tecout.OpenFile("../TESTOUT/trianglemeshconv.plt");
		inparam.lowerB= {-0.0, -0.0,-0.0};
		inparam.upperB= {1.0, 1.0, 1.0};
		inparam.distanceTol = 1e-3;
		// inparam.edgelengths = {0.2};
		inparam.edgelengths = {0.25};

		BuildBlockGrid({3,1,1}, meshdomain);
		dimDomain[0] = {0,3.0};
		// dimDomain[0] = {0,3.0};
		dimDomain[1] = {0,1.0};
		dimDomain[2] = {0,1.0};
		meshdomain.Scale(dimDomain);
		meshdomain.PrepareForUse();
		tecout.PrintMesh(meshdomain);

		tetgen::input::RSVSGRIDS(meshdomain,tetin, inparam);


		tetin.save_nodes("../TESTOUT/testtetgen/rsvs_3cell_grid");
		tetin.save_poly("../TESTOUT/testtetgen/rsvs_3cell_grid");
		

		// Tetrahedralize the PLC. Switches are chosen to read a PLC (p),
		//   do quality mesh generation (q) with a specified quality bound
		//   (1.414), and apply a maximum volume constraint (a0.1).
		inparam.command = "Qpkqnnvefm"; 
		tetrahedralize(inparam.command.c_str(), &tetin, &tetout);

		// tetout.save_nodes("../TESTOUT/testtetgen/rsvsout_3cell_2body");
		// tetout.save_elements("../TESTOUT/testtetgen/rsvsout_3cell_2body");
		// tetout.save_faces("../TESTOUT/testtetgen/rsvsout_3cell_2body");
		std::cout << "Finished the tettgen process " << std::endl;
		meshtet = tetgen::output::TET2MESH(tetout);
		tecout.PrintMesh(meshtet);


		std::cout << " Meshed the tetrahedralization" << std::endl;
		tetgen::input::POINTGRIDS(meshtet,tetin3, inparam);
		inparam.command = "Qv"; 
		tetrahedralize(inparam.command.c_str(), &tetin3, &tetout3);
		meshvoro = tetgen::output::VORO2MESH(tetout3);
		// tecout.PrintMesh(meshvoro,0,0,2);
		tecout.PrintMesh(meshvoro);
		tecout.PrintMesh(meshvoro,0,0,3);
		tecout.PrintMesh(meshvoro,0,0,4);
		std::cout << " Meshed the voronization" << std::endl;

		inparam.edgelengths = {0.15};
		meshvoro.PrepareForUse();


		std::vector<int> subs;
		subs=ConcatenateVectorField(meshvoro.volus, &volu::surfind,
			0, meshvoro.volus.size());
		std::cout << "Number of surfs " << subs.size() << std::endl;
		sort(subs);
		unique(subs);
		std::cout << "Number of surfs " << subs.size() << std::endl;
		inparam.command = "Qpkqnnefm"; 
		tetgen::input::RSVSGRIDS(meshvoro,tetin2, inparam);

		tetin2.save_nodes("../TESTOUT/testtetgen/rsvs_voro_grid");
		tetin2.save_poly("../TESTOUT/testtetgen/rsvs_voro_grid");

		tetrahedralize(inparam.command.c_str(), &tetin2, &tetout2);
		std::cout << "Finished the tettgen process " << std::endl;
		meshtet2 = tetgen::output::TET2MESH(tetout2);
		tecout.PrintMesh(meshtet2);

	} catch (exception const& ex) {
		cerr << "Exception: " << ex.what() <<endl; 
		return 1;
	}
	return 0;
}
int tetgen::test::RSVSVORO()
{
	std::vector<double> vecPts;
	mesh vosMesh, snakMesh;
	tecplotfile tecout;

	int nErrors=0;
	const char* tecoutStr =  "../TESTOUT/rsvs_voro.plt";
	vector<int> numCells = {
		0,
		1,
		2,3,4,5,10,20,100,1000
	};
	vector<double> numEdge = {
		0.05, 0.1, 0.3
	};
	// std :: cin >> nPts;
	tecout.OpenFile(tecoutStr);

	for(auto i : numCells){
		for(auto j : numEdge){
			nErrors += tetgen::test::RSVSVOROFunc(i, j, tecoutStr);
		}
	}
	
	return nErrors;
}

int tetgen::test::RSVSVOROFunc_contain(int nPts, double distanceTol, const char* tecoutStr)
{
	std::vector<double> vecPts;
	mesh vosMesh, snakMesh, meshPts;
	tecplotfile tecout;
	tetgen::apiparam inparam;
	try {
		cout << "__________________________________________________" << endl;
		cout << "Start Voronoi mesh generation for " << nPts 
			<< " points"<< endl;
		cout << "__________________________________________________" << endl;
		for (int i = 0; i < nPts*4; ++i)
		{
			vecPts.push_back(double(abs(rand())%32767)/32767.0);
		}
		for (int i = 0; i < 8; ++i)
		{
			vecPts.push_back((i%2));	
			vecPts.push_back(((i/2)%2));	
			vecPts.push_back(((i/4)%2));	
			vecPts.push_back(0);	
		}
		// tecout.OpenFile(tecoutStr, "a");
		tecout.OpenFile(tecoutStr,"a");

		inparam.lowerB= {-0.0, -0.0,-0.0};
		inparam.upperB= {1.0, 1.0, 1.0};
		inparam.distanceTol = 1e-3;
		// inparam.edgelengths = {0.2};
		inparam.edgelengths = {0.1};
		inparam.distanceTol = distanceTol;

		auto vecMap = tetgen::RSVSVoronoiMesh(vecPts, vosMesh, snakMesh, inparam);

		meshPts.Init(nPts+8,0,0,0);
		meshPts.PopulateIndices();
		for (int i = 0; i < nPts+8; ++i)
		{	
			meshPts.verts[i].coord.assign(vecPts.begin()+i*4,
				vecPts.begin()+i*4+3);
		}
		meshPts.PrepareForUse();
		// tecout.PrintMesh(meshPts,0,nPts+distanceTol,4);
		int nVecMap = vecMap.size();
		for (int i = 0; i < nVecMap; ++i)
		{
			std::cout << snakMesh.ParentElementIndex(vecMap[i]) << " ";
			tecout.PrintMesh(vosMesh,nPts*3+1,i,5,
				{snakMesh.ParentElementIndex(vecMap[i])});
			tecout.PrintMesh(snakMesh,nPts*3+2,i,5,{vecMap[i]});
			tecout.PrintMesh(meshPts,nPts*3+3,i,4,{i+1});
		}
		cout << "__________________________________________________" << endl;

	} catch (exception const& ex) {
		cerr << "Exception: " << ex.what() <<endl; 
		cerr << "for " << __PRETTY_FUNCTION__ << "with nPoints = " << nPts << endl;
		return 1;
	}
	return 0;
}

int tetgen::test::RSVSVOROFunc(int nPts, double distanceTol, const char* tecoutStr)
{
	std::vector<double> vecPts;
	mesh vosMesh, snakMesh, meshPts;
	tecplotfile tecout;
	tetgen::apiparam inparam;
	try {
		cout << "__________________________________________________" << endl;
		cout << "Start Voronoi mesh generation for " << nPts 
			<< " points"<< endl;
		cout << "__________________________________________________" << endl;
		for (int i = 0; i < nPts*4; ++i)
		{
			vecPts.push_back(double(abs(rand())%32767)/32767.0);
		}
		for (int i = 0; i < 8; ++i)
		{
			vecPts.push_back((i%2));	
			vecPts.push_back(((i/2)%2));	
			vecPts.push_back(((i/4)%2));	
			vecPts.push_back(0);	
		}
		tecout.OpenFile(tecoutStr, "a");

		inparam.lowerB= {-0.0, -0.0,-0.0};
		inparam.upperB= {1.0, 1.0, 1.0};
		inparam.distanceTol = 1e-3;
		// inparam.edgelengths = {0.2};
		inparam.edgelengths = {0.1};
		inparam.distanceTol = distanceTol;

		tetgen::RSVSVoronoiMesh(vecPts, vosMesh, snakMesh, inparam);
		tecout.PrintMesh(snakMesh,1,nPts+distanceTol);
		tecout.PrintMesh(vosMesh,2,nPts+distanceTol);
		meshPts.Init(nPts+8,0,0,0);
		meshPts.PopulateIndices();
		for (int i = 0; i < nPts+8; ++i)
		{	
			meshPts.verts[i].coord.assign(vecPts.begin()+i*4,
				vecPts.begin()+i*4+3);
		}
		meshPts.PrepareForUse();
		tecout.PrintMesh(meshPts,3,nPts+distanceTol,4);
		cout << "__________________________________________________" << endl;

	} catch (exception const& ex) {
		cerr << "Exception: " << ex.what() <<endl; 
		cerr << "for " << __PRETTY_FUNCTION__ << "with nPoints = " << nPts << endl;
		return 1;
	}
	return 0;
}

int tetgen::test::RSVSVORO_Contain()
{
	std::vector<double> vecPts;
	mesh vosMesh, snakMesh;
	tecplotfile tecout;

	int nErrors=0;
	const char* tecoutStr =  "../TESTOUT/rsvs_voro_contain.plt";
	vector<int> numCells = {
		0,
		1,
		// 2,3,4,5,10,20,100,1000
	};
	vector<double> numEdge = {
		// 0.05,
		0.1,
		// 0.3
	};
	// std :: cin >> nPts;
	tecout.OpenFile(tecoutStr);

	for(auto i : numCells){
		for(auto j : numEdge){
			nErrors += tetgen::test::RSVSVOROFunc_contain(i, j,tecoutStr);
		}
	}
	
	return nErrors;
}
