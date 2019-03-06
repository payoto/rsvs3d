#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <unordered_map>
#include <cmath>

#include "tetgen_rsvs_api.hpp"
#include "tetgen.h"
#include "mesh.hpp"
#include "snake.hpp"
#include "snakevel.hpp"
#include "voxel.hpp"
#include "arraystructures.hpp"
#include "postprocessing.hpp"
#include "meshrefinement.hpp"
#include "warning.hpp"

void load_tetgen_testdata(mesh &snakeMesh, mesh &voluMesh, snake &snakein, mesh &triMesh){
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

mesh BuildDomain(const std::array<double,3> &lowerB, 
	const std::array<double,3> &upperB, double tolInner=0.0){
	/*
	Builds a parralelipipede domain stretching from lowerB to upperB

	`lowerB` and `upperB` must be vectors of size 3. `tolInner` specifies an
	offset expanding the box by a set amount.
	*/	
	mesh cube;
	int count; 

	if(lowerB.size()!=3 || upperB.size()!=3){
		std::cerr << "Error in " << __PRETTY_FUNCTION__ << " vectors must be"
			" of size 3" << std::endl;
		RSVS3D_ERROR_ARGUMENT("Input vectors must be of size 3");
	}

	std::array<int, 3> dimGrid={1, 1, 1};
	BuildBlockGrid(dimGrid, cube);

	count = cube.verts.size();
	for (int i = 0; i < count; ++i) {
		for (int j = 0; j < 3; ++j)	{
			cube.verts[i].coord[j] = (cube.verts[i].coord[j]*((upperB[j]+tolInner)-(lowerB[j]-tolInner))) + (lowerB[j]-tolInner);
		}
	}

	return cube;
}


void MeshData2Tetgenio(const mesh &meshgeom, tetgen::io_safe &tetin,
	int facetOffset, int pointOffset, int pointMarker, 
	const std::vector<double> &pointMtrList, const std::vector<double> &facetConstr,
	int facetConstrOffset){
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


void Mesh2Tetgenio(const mesh &meshgeom, const mesh &meshdomain,
	tetgen::io_safe &tetin, int numHoles){
	/*
	Converts a mesh into the safe allocation tetgenio format.
	
	Rules of conversion are: 
		- volu becomes a facet
		- surf becomes a polygon in the facet
		- edges are not logged
		- points come as a list
	*/
	

	tetin.firstnumber=0;

	tetin.numberoffacets = meshgeom.surfs.size() + meshdomain.surfs.size();
	tetin.numberofholes = numHoles;
	tetin.numberofpoints = meshgeom.verts.size() + meshdomain.verts.size();

	tetin.numberofpointmtrs = 1;
	tetin.numberoffacetconstraints = 2;


	tetin.allocate();

	MeshData2Tetgenio(meshgeom, tetin, 0, 0, 1, {0.03},{0.0},0);
	MeshData2Tetgenio(meshdomain, tetin, meshgeom.surfs.size(), 
		meshgeom.verts.size(), -1, {-1.0},{0.0},1);
}


void Mesh2TetgenioPoints(const mesh &meshgeom, const mesh &meshdomain,
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
	tetin.numberofpoints = meshgeom.verts.size() + meshdomain.verts.size();

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

HashedVector<int, int> GroupCloseSnaxels(const snake &snakein, double distTol){
	/*
	Merges snaxel which are `distTol` to the end of their carrying edge.
	*/
	HashedVector<int, int> closeVert;
	std::vector<int> sameSnaxs, rmvInds;
	int nSnax;

	nSnax = snakein.snaxs.size();
	closeVert.vec.assign(nSnax,-1);

	// This piece of code is going to be incomplete as it won't test
	// for multiple snaxels very close on the same edge
	// 
	// closeVert is a list of vertices close
	for (int i = 0; i < nSnax; ++i)	{
		if((snakein.snaxs(i)->d<distTol)){ 
			closeVert.vec[i] = snakein.snaxs(i)->fromvert;
		} else if (((1-snakein.snaxs(i)->d)<distTol)){
			closeVert.vec[i] = snakein.snaxs(i)->tovert;
		} 
	}

	closeVert.GenerateHash();
	return closeVert;
}

void TestVertClose(int vertIndIn, std::vector<bool> &isSnaxDone, 
	const mesh &meshin, double distTol,
	std::vector<int> &sameEdges){
	/*
	recursive function for finding edges which are short and should be collapsed
	*/
	int nextvertIndIn;
	if(!isSnaxDone[vertIndIn]){
			isSnaxDone[vertIndIn]=true;
			for (auto eInd : meshin.verts(vertIndIn)->edgeind){
				double d = 0.0;
				for(int i = 0; i< 3; ++i){
					d+=pow((meshin.verts.isearch(
							meshin.edges.isearch(eInd)->vertind[0]
						)->coord[i]-
						meshin.verts.isearch(
							meshin.edges.isearch(eInd)->vertind[1]
							)->coord[i]),2.0);
				}
				if(d<distTol){
					sameEdges.push_back(eInd);
					for(int i = 0; i< 2; ++i){
						nextvertIndIn = meshin.verts.find(
								meshin.edges.isearch(eInd)->vertind[i]
							);
						TestVertClose(nextvertIndIn, isSnaxDone, 
							meshin,  distTol, sameEdges);
					}
				}
				
			}
		}
}

HashedVector<int, int> GroupCloseVertices(const mesh &meshin, double distTol){
	/*
	Groups vertices which are `distTol` close to each other.
	*/
	std::vector<bool> isSnaxDone, isEdgeDone;
	HashedVector<int, int> closeVert;
	std::vector<int> sameEdges, rmvInds;
	int nSnax, nEdge, nGroup=1;

	distTol = distTol*distTol;
	nSnax = meshin.verts.size();
	nEdge = meshin.edges.size();
	closeVert.vec.assign(nSnax,-1);
	isSnaxDone.assign(nSnax,false);
	isEdgeDone.assign(nEdge,false);
	sameEdges.reserve(10);
	// This piece of code is going to be incomplete as it won't test
	// for multiple snaxels very close on the same edge
	// 
	// closeVert is a list of vertices close
	
	for (int i = 0; i < nSnax; ++i)	{
		sameEdges.clear();
		TestVertClose(i, isSnaxDone,meshin, distTol,sameEdges);

		if (sameEdges.size()>0){
			for (auto eInd : sameEdges)
			{
				closeVert.vec[meshin.verts.find(
						meshin.edges.isearch(eInd)->vertind[0]
					)]=nGroup;
				// std::cout << meshin.verts.find(
				// 		meshin.edges.isearch(eInd)->vertind[0]
				// 	) << " " << nGroup << "|";
				closeVert.vec[meshin.verts.find(
						meshin.edges.isearch(eInd)->vertind[1]
					)]=nGroup;
				// std::cout << meshin.verts.find(
				// 		meshin.edges.isearch(eInd)->vertind[1]
				// 	) << " " << nGroup << "|";
			}
			// std::cout << "|";
			++nGroup;
		}
	}
	// std::cout << std::endl;

	// std::cout << nGroup << std::endl;
	closeVert.GenerateHash();
	return closeVert;
}

int FindVertexHole(int vertInd, const mesh &meshin, const std::vector<bool> &vertIn,
	const HashedVector<int, int> &uncertainVert, std::vector<bool> &vertExplored){
	/*
	Starting from vertInd look for a vertex which is definitely a hole.

	From this vertex check if it is in? 
	Check if it is in closeVert 
	-> yes and no, we're done
	-> yes and yes, explore this vertices neighbours
	-> no, nothing to see here go to next from vertex
	*/

	std::vector<int> currVertList, nextEdgeList;
	int returnInd=-1;
	int ii, cntII, jj, cntJJ, vertPos, vertFound;

	currVertList.reserve(20);
	nextEdgeList.reserve(20);
	currVertList.push_back(vertInd);

	while(currVertList.size()>0){
		nextEdgeList.clear();
		cntII=currVertList.size();
		for(ii=0; ii<cntII; ++ii){
			vertPos = meshin.verts.find(currVertList[ii]);
			if(!vertExplored[vertPos] && vertIn[vertPos]){
				vertExplored[vertPos]=true;
				vertFound = uncertainVert.find(currVertList[ii]);
				if(vertFound==-1){
					// Found our point!
					returnInd = currVertList[ii];
					// These operations return out of the function
					currVertList.clear();
					break;
				} else {
					// Add the edges to a list
					cntJJ = meshin.verts(vertPos)->edgeind.size();
					for (jj=0; jj< cntJJ; ++jj){
						nextEdgeList.push_back(
							meshin.verts(vertPos)->edgeind[jj]);
					}
				}
			}
		}
		if (returnInd==-1){
			// Assign currVertList
			currVertList.clear();
			if (nextEdgeList.size()>0){
				currVertList = ConcatenateVectorField(meshin.edges,
					&edge::vertind, nextEdgeList);
			}
		}
	}

	return (returnInd);
}

std::vector<int> FindHolesInSnake(const snake &snakein,
	const HashedVector<int, int> &uncertainVert){
	/*
	Resolves the problem of finding points from which to initiate voids
	in a snake.

	Takes in a snake and a vector of indices to which snaxels were considered
	close. These have to be treated separately as there is some uncertainty
	about wether they are inside or outside. uncertain Vert is a list of
	indices where the algorithm is known to be invalid
	*/

	int ii, cntII, jj, cntJJ, kk, cntKK, ll;
	auto& vols = snakein.snakeconn.volus;
	std::vector<bool> vertExplored;
	std::vector<int> holeInds;
	int vertAct, holeInd=-1;
	// const vector<int>  &surfEdgeind=currVertList;
	vertExplored.assign(snakein.isMeshVertIn.size(),false);
	// Go through the volume objects vertices,
		// looking for a 'from vertex'
			// From this vertex check if it is in? 
			// Check if it is in closeVert 
			// -> yes and no, we're done
			// -> yes and yes, explore this vertices neighbours
			// -> no, nothing to see here go to next from vertex
	cntII = snakein.snakeconn.volus.size();
	holeInds.reserve(cntII);
	for (ii = 0; ii< cntII; ++ii){
		cntJJ = vols(ii)->surfind.size();
		for (jj = 0; jj< cntJJ; ++jj){
			auto surfEdgeind = snakein.snakeconn.surfs.isearch(
				vols(ii)->surfind[jj])->edgeind;
			cntKK = surfEdgeind.size();
			for (kk = 0; kk<cntKK; ++kk){
				for (ll = 0; ll<2; ++ll){
					vertAct = snakein.snaxs.isearch( 
						snakein.snakeconn.edges(surfEdgeind[kk])->vertind[ll]
						)->fromvert;
					holeInd = FindVertexHole(vertAct, *(snakein.snakemesh), 
						snakein.isMeshVertIn, uncertainVert, vertExplored);
					if(holeInd!=-1){break;}
				}
				if(holeInd!=-1){break;}
			}
			if(holeInd!=-1){break;}
		}
		if(holeInd!=-1){holeInds.push_back(holeInd);}
	}
	return (holeInds);
}

double DomInter(double x, double y1, double y2){

	/*Interpolation function*/
	return x*(y2-y1)+y1;
}

/**
 * @brief      Builds a series of domains with different edge 
 *  properties controlling the interpolation of the metric.
 *
 * @param[in]  outerLowerB       The outer lower b
 * @param[in]  outerUpperB       The outer upper b
 * @param[in]  innerLowerB       The inner lower b
 * @param[in]  innerUpperB       The inner upper b
 * @param[in]  nSteps            The steps
 * @param      vertPerSubDomain  The vertical per sub domain
 *
 * @return     The cut cell domain.
 * 
	These also serve to avoid having a badly conditioned initial triangulation
	with very small edges.
	
	`nSteps` is the number of total domains.
	0 will return an empty mesh, 1 will return a mesh of the inner bound
	2 will return inner and outer bounds, 

			{
				{DomInter(x, innerLowerB[0], outerLowerB[0]),
				DomInter(x, innerUpperB[0], outerUpperB[0])},
				{DomInter(x, innerLowerB[1], outerLowerB[1]),
				DomInter(x, innerUpperB[1], outerUpperB[1])},
				{DomInter(x, innerLowerB[2], outerLowerB[2]),
				DomInter(x, innerUpperB[2], outerUpperB[2])}
			}

		meshtemp = BuildDomain(
			{DomInter(x, innerLowerB[0], outerLowerB[0]),
				DomInter(x, innerLowerB[1], outerLowerB[1]),
				DomInter(x, innerLowerB[2], outerLowerB[2])},
			{DomInter(x, innerUpperB[0], outerUpperB[0]),
				DomInter(x, innerUpperB[1], outerUpperB[1]),
				DomInter(x, innerUpperB[2], outerUpperB[2])}
			);
 */
mesh BuildCutCellDomain(const std::array<double,3> &outerLowerB, 
	const std::array<double,3> &outerUpperB,
	const std::array<double,3> &innerLowerB, const std::array<double,3> &innerUpperB, int nSteps, 
	std::vector<int> &vertPerSubDomain){

	mesh meshdomain, meshtemp;
	double x;
	std::array<std::array<double, 2>,3> scaleDom;

	// Special case if asked for  0 subdomains
	if(nSteps==0){
		meshdomain.Init(0, 0, 0, 0);
		meshdomain.PrepareForUse();
		return meshdomain;
	}

	// Start from inner bound
	// Then step up 
	vertPerSubDomain.clear();
	meshdomain = BuildDomain(innerLowerB, innerUpperB, 0.1);
	meshdomain.PrepareForUse();
	meshtemp = meshdomain;
	meshtemp.PrepareForUse();
	vertPerSubDomain.push_back(meshtemp.verts.size());

	for (int i = 1; i < nSteps; ++i){
		x =double(i)/double(nSteps-1);
		for (int j=0; j<3;++j){	
			scaleDom[j] ={DomInter(x, innerLowerB[j], outerLowerB[j]),
				DomInter(x, innerUpperB[j], outerUpperB[j])};
		}

		meshtemp.Scale(scaleDom);
		meshtemp.PrepareForUse();
		meshdomain.MakeCompatible_inplace(meshtemp);
		meshdomain.Concatenate(meshtemp);
		meshdomain.PrepareForUse();

		vertPerSubDomain.push_back(meshtemp.verts.size());
	}
	
	return meshdomain;
}

void PrepareSnakeForCFD(const snake &snakein, const tetgen::apiparam &tetgenParam,
	mesh &meshgeom, std::vector<double> &holeCoords){
	/*
	Prepares the snake to be used for CFD, removes duplicate points and 
	triangulates it.

	Coordinates of holes are also returned.
	*/
	mesh  meshtemp;
	triangulation triRSVS;
	int nHoles;
	// int nPtsHole, kk, count;

	std::vector<int>  holeIndices;

	// tecplotfile tecout;

	// Cleanup the mesh
	// tecout.OpenFile("../TESTOUT/rsvs_tetgen_mesh.plt");
	meshtemp=snakein.snakeconn;
	// tecout.PrintMesh(meshtemp);
	// meshtemp.displight();

	auto groupedVertices = GroupCloseSnaxels(snakein, tetgenParam.distanceTol);
	meshtemp.MergeGroupedVertices(groupedVertices);
	meshtemp.RemoveSingularConnectors();

	meshtemp.PrepareForUse();
	meshtemp.TestConnectivityBiDir(__PRETTY_FUNCTION__);
	meshtemp.TightenConnectivity();
	meshtemp.OrderEdges();
	// tecout.PrintMesh(meshtemp);
	// meshtemp.displight();

	// Triangulation preparation
	triRSVS.stattri.clear();
	triRSVS.trivert.clear();
	triRSVS.PrepareForUse();
	TriangulateMesh(meshtemp,triRSVS);

	triRSVS.PrepareForUse();
	triRSVS.CalcTriVertPos();
	MeshTriangulation(meshgeom,meshtemp,triRSVS.stattri,
		triRSVS.trivert);
	meshgeom.PrepareForUse();
	// meshgeom.displight();
	meshgeom.TightenConnectivity();
	meshgeom.OrderEdges();

	// Find number of holes
	holeIndices=FindHolesInSnake(snakein, groupedVertices);
	nHoles = holeIndices.size();

	holeCoords.assign(nHoles*3,0.0);
	// Assign the holes
	for (int i = 0; i < nHoles; ++i){
		for (int j =0; j<3; ++j){
			holeCoords[i*3+j] = snakein.snakemesh->verts.isearch(
				holeIndices[i])->coord[j];
		}
	}
}

void TetgenInput_RSVS2CFD(const snake &snakein, tetgen::io_safe &tetin,
	const tetgen::apiparam &tetgenParam){
	/*
	Processes a snake object into a tetgen input object.

	This function processes snake objects into the safe tetgenio object.
	This can then used to generate a tetrahedral mesh or to save it as the
	input file.
	*/

	mesh meshdomain ,meshdomain2, meshgeom;
	triangulation triRSVS;
	int nHoles;
	// int nPtsHole, kk, count;

	std::vector<int>  holeIndices;
	std::array<double, 3> upperB, lowerB;
	std::vector<int> vertPerSubDomain;
	std::vector<double> holeCoords;

	PrepareSnakeForCFD(snakein, tetgenParam, meshgeom, holeCoords);

	meshgeom.ReturnBoundingBox(lowerB, upperB);
	meshdomain = BuildCutCellDomain(tetgenParam.lowerB, tetgenParam.upperB,
		lowerB, upperB, tetgenParam.edgelengths.size(), vertPerSubDomain);
	
	nHoles=holeCoords.size()/3;
	Mesh2Tetgenio(meshgeom, meshdomain, tetin, nHoles);

	// std::cout<< std::endl << "Number of holes " << nHoles << std::endl;

	// Assign the holes

	int count = holeCoords.size();
	for (int i = 0; i < count; ++i){
		tetin.holelist[i] = holeCoords[i];	
	}

	// Assign domain point metrics
	int startPnt = meshgeom.verts.size();
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



void TetgenInput_RSVSGRIDS(const mesh &meshdomain, const mesh &meshboundary,
	tetgen::io_safe &tetin, const tetgen::apiparam &tetgenParam){
	/*
	Processes a snake object into a tetgen input object.

	This function processes snake objects into the safe tetgenio object.
	This can then used to generate a tetrahedral mesh or to save it as the
	input file.
	*/

	triangulation triRSVS;
	int nHoles;
	// int nPtsHole, kk, count;

	std::vector<int>  holeIndices;
	std::vector<int> vertPerSubDomain;
	std::vector<double> holeCoords;

	holeCoords.clear();

	vertPerSubDomain.push_back(meshdomain.verts.size());
	

	nHoles=holeCoords.size()/3;
	Mesh2Tetgenio(meshboundary, meshdomain, tetin, nHoles);

	// std::cout<< std::endl << "Number of holes " << nHoles << std::endl;

	// Assign the holes

	int count = holeCoords.size();
	for (int i = 0; i < count; ++i){
		tetin.holelist[i] = holeCoords[i];	
	}

	// Assign domain point metrics
	int startPnt = meshboundary.verts.size();
	for (int i = 0; i < int(tetgenParam.edgelengths.size()); ++i){
		tetin.SpecifyTetPointMetric(startPnt,  vertPerSubDomain[i],
			{tetgenParam.edgelengths[i]});
		startPnt = startPnt + vertPerSubDomain[i];
	}

	// tecout.PrintMesh(meshdomain);
}

void TetgenInput_RSVSGRIDS(const mesh &meshdomain,
	tetgen::io_safe &tetin, const tetgen::apiparam &tetgenParam){
	mesh meshboundary;
	meshboundary.Init(0,0,0,0);
	meshboundary.PrepareForUse();

	TetgenInput_RSVSGRIDS(meshdomain, meshboundary, 
		tetin, tetgenParam);
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
void TetgenInput_POINTGRIDS(const mesh &meshdomain, tetgen::io_safe &tetin,
	const tetgen::apiparam &tetgenParam, bool generateVoroBound=false){
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
		std::vector<double> vertsToAdd;
		std::array<double, 3> vertModif = {0,0,0};
		vertsToAdd.reserve(100);
		// add the points associated with the vertices
		for (int l = 0; l < 3; ++l)
		{
			vertModif = {0,0,0};
			vertModif[l] = 1.0;
			for (int i = 0; i < 8; ++i)
			{	
				for (int j = 0; j < 3; ++j)
				{
					int k = pow(2,j);
					double delta = (tetgenParam.upperB[j]-tetgenParam.lowerB[j])*vertModif[l];
					vertsToAdd.push_back(
						(((i/k)%2) ? 
							tetgenParam.upperB[j]+delta 
							: tetgenParam.lowerB[j]-delta)
						);
				}
			}
		}
		// add the points associated with the edges



	} else{ 
		meshgeom.Init(0, 0, 0, 0);
	}
	meshgeom.PrepareForUse();
	
	Mesh2TetgenioPoints(meshgeom, meshdomain, tetin);
}
#pragma GCC diagnostic pop

void TetgenOutput_SU2(const char* fileName, const tetgenio &tetout){
	/*
	Ouputs a tetgen io object to an SU2 mesh file.
	
	Args:
		fileName: is a string with the path to the target mesh fileName
		tetout: is tetgenio object (safe) 

	Raises:
		- invalid_argument if `fileName` cannot be opened. 

	File format:
	https://su2code.github.io/docs/Mesh-File/
	*/
	std::ofstream meshFile;
	std::vector<int> diffBounds, numBounds;
	int DEINCR;

	meshFile.open(fileName);
	// De-increments indices if tetgen indices start from 1 to match SU2
	// which is 0 indexed
	DEINCR = tetout.firstnumber ? -1 : 0;

	// If file was not opened correctly
	if(meshFile.fail()){
		std::cerr << "Error: " << fileName << std::endl <<
			 " in " << __PRETTY_FUNCTION__ << std::endl; 
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


void CloseVoronoiMesh(mesh &meshout, tetgen::io_safe &tetout, 
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
		lStep = tetgen::voronoi::ProjectRay(3, boundBox, 
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

void FlattenBoundaryFaces(mesh &meshin){
	/*
	Flattens non flat border faces of a mesh by triangulating them

	*/
	mesh meshout;
	vector<int> subList, vertList;
	std::array<bool, 3> isFlat;
	triangulation triangleRSVS;
	int nSurf;
	double tol=1e-9;

	nSurf = meshin.surfs.size();
	subList.reserve(nSurf);

	triangleRSVS.stattri.clear();
	triangleRSVS.trivert.clear();
	triangleRSVS.PrepareForUse();
	// Find Non flat boundary faces (only checks for one coordinate being the same)
	for (int i = 0; i < nSurf; ++i)
	{
		if(meshin.surfs(i)->isBorder  && meshin.surfs(i)->edgeind.size()>3){
			meshin.surfs(i)->OrderedVerts(&meshin,vertList);
			vertList = meshin.verts.find_list(vertList);

			isFlat.fill(true);
			int count = vertList.size();
			for (int j = 1; j < count; ++j)
			{
				for (int k = 0; k < 3; ++k)
				{
					isFlat[k] = isFlat[k] && (fabs(meshin.verts(vertList[j])->coord[k]-
										meshin.verts(vertList[0])->coord[k])<tol);
				}
			}
			if(!(isFlat[0] || isFlat[1] || isFlat[2])){
				subList.push_back(i);
			}
		}
	}
	// triangulate the container
	if(subList.size()>0){
		TriangulateContainer(meshin,triangleRSVS , 1, subList); 
		triangleRSVS.meshDep = &meshin;


		triangleRSVS.PrepareForUse();
		triangleRSVS.CalcTriVertPos();
		triangleRSVS.PrepareForUse();
		// meshin.displight();
		auto newSurfs = ConcatenateVectorField(meshin.surfs,
			&surf::edgeind, subList);
		// cerr << "Num new surfs " << newSurfs.size()-subList.size() 
		// 	<< " from " << subList.size() << endl;
		MeshTriangulation(meshout,meshin,triangleRSVS.stattri,
			triangleRSVS.trivert);
		meshin=meshout;
		meshin.PrepareForUse();
		// meshin.displight();
		meshin.TightenConnectivity();
		meshin.OrderEdges();
		// meshin.displight();
	}
}

tetgen::dombounds TetgenOutput_GetBB(tetgen::io_safe &tetout){

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

mesh TetgenOutput_VORO2MESH(tetgen::io_safe &tetout){
	/*
	Translates a tetgen output to the RSVS native mesh format 
	*/
	mesh meshout;

	int nVerts, nEdges, nSurfs, nVolus;
	int count,INCR, DEINCR,n;
	vector<int> rayInd, rayInd2;
	tetgen::dombounds boundBox=TetgenOutput_GetBB(tetout);

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
	CloseVoronoiMesh(meshout, tetout, rayInd, DEINCR, boundBox);

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

mesh TetgenOutput_TET2MESH(tetgen::io_safe &tetout){
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

	// if `tetrahedralize` has not been run with the correct flags connectivity is
	// unavailable, throw an error. (tests for flag -nn)
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

namespace voronoimesh {
	mesh Points2Mesh(const std::vector<double> &vecPts){
		/*
		Takes in a set of points and returns a mesh of points ready for
		voronisation.
		*/
		mesh meshpts;
		int nProp = 3+1; // Number of properties per point
		if((vecPts.size()%nProp)!=0){
			std::cerr << "Error in: " << __PRETTY_FUNCTION__ << std::endl;
			RSVS3D_ERROR_ARGUMENT("Vector of points is an invalid size");
		}
		int nPts = vecPts.size()/nProp;
		meshpts.Init(nPts, 0, 0, 0);

		for (int i = 0; i < nPts; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				meshpts.verts[i].coord[j] = vecPts[i*nProp+j]; 
			}
		}
		meshpts.verts.PrepareForUse();
		return (meshpts);
	}

	std::vector<bool> Points2VoroAndTetmesh(const std::vector<double> &vecPts,
		mesh &voroMesh, mesh &tetMesh, const tetgen::apiparam &inparam){
		/*
		Turns a set of points into the voronisation and tetrahedralisation.
		*/
		tetgen::io_safe tetinV, tetoutV, tetinT, tetoutT;
		mesh ptsMesh = Points2Mesh(vecPts);
		std::string cmd;
		int INCR;
		std::vector<int> voroPts;
		std::vector<bool> tetSurfInVoro;
		std::vector<double> lowerB, upperB;

		TetgenInput_POINTGRIDS(ptsMesh, tetinV, inparam);
		cmd = "Qv"; 
		tetrahedralize(cmd.c_str(), &tetinV, &tetoutV);
		voroMesh = TetgenOutput_VORO2MESH(tetoutV);
		// Crop mesh 
		lowerB.assign(3, 0.0);
		upperB.assign(3, 0.0);
		for (int i = 0; i < 3; ++i)
		{
			double delta = inparam.upperB[i]-inparam.lowerB[i];
			lowerB[i]=inparam.lowerB[i] + delta/2.0;
			upperB[i]=inparam.upperB[i] + delta/2.0;
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

		TetgenInput_RSVSGRIDS(voroMesh, meshdomain, tetinT, inparam);
		tetinT.save_nodes("../TESTOUT/testtetgen/voro/rsvs");
		tetinT.save_poly("../TESTOUT/testtetgen/voro/rsvs");
		cmd = "QpqmnnefO9/7"; 
		tetrahedralize(cmd.c_str(), &tetinT, &tetoutT);
		tetMesh = TetgenOutput_TET2MESH(tetoutT);

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
	std::vector<bool> BoundaryFacesFromPoints(const mesh &meshin,
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

}


void RSVSVoronoiMesh(const std::vector<double> &vecPts,mesh &vosMesh, mesh &snakMesh){
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
	tetgen::apiparam inparam;
	std::vector<int> vertsVoro;
	std::vector<int> vertBlock, elmMapping;


	inparam.lowerB= {-0.0, -0.0,-0.0};
	inparam.upperB= {1.0, 1.0, 1.0};
	inparam.distanceTol = 1e-3;
	// inparam.edgelengths = {0.2};
	inparam.edgelengths = {0.1};

	// Step 1 - 2 
	auto boundFaces = voronoimesh::Points2VoroAndTetmesh(vecPts,vosMesh, snakMesh, inparam);
	voroMesh = vosMesh;
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored  "-Wunused-variable"
	int nBlocks = snakMesh.ConnectedVolumes(elmMapping, boundFaces);
	#pragma GCC diagnostic pop
	#ifdef SAFE_ALGO
	if(nBlocks!=voroMesh.volus.size()){
		std::cerr << "Error : nBlocks (" << nBlocks 
			<< ")!=nVolus Voromesh (" << voroMesh.volus.size() <<")" 
			<< std::endl;
		RSVS3D_ERROR_LOGIC("Voromesh and flooding do not match");
	}
	#endif //SAFE_ALGO

	CoarsenMesh(snakMesh,vosMesh,elmMapping);
	snakMesh.AddParent(&vosMesh,elmMapping);
	vosMesh=voroMesh;
	snakMesh.PrepareForUse();
	snakMesh.OrientFaces();
	vosMesh.PrepareForUse();

	// Step 6 - Check original points containments
	std::vector<int> vecPtsMapping = snakMesh.VertexInVolume(vecPts, 4);
	DisplayVector(vecPtsMapping);
	cout << endl;
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

// Test code
int test_tetgenapi(){

	return(tetcall());
}

int tetcall_CFD()
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
		load_tetgen_testdata(snakeMesh, voluMesh, snakein, triMesh);
		TetgenInput_RSVS2CFD(snakein, tetin, inparam);


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
		TetgenOutput_SU2("../TESTOUT/testtetgen/rsvs_3cell_2body.su2",tetout);
		std::cout << "Finished the tettgen process " << std::endl;
	} catch (exception const& ex) {
		cerr << "Exception: " << ex.what() <<endl; 
		return 1;
	}
	return 0;
}

int tetcall()
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

		TetgenInput_RSVSGRIDS(meshdomain,tetin, inparam);


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
		meshtet = TetgenOutput_TET2MESH(tetout);
		tecout.PrintMesh(meshtet);


		std::cout << " Meshed the tetrahedralization" << std::endl;
		TetgenInput_POINTGRIDS(meshtet,tetin3, inparam);
		inparam.command = "Qv"; 
		tetrahedralize(inparam.command.c_str(), &tetin3, &tetout3);
		meshvoro = TetgenOutput_VORO2MESH(tetout3);
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
		TetgenInput_RSVSGRIDS(meshvoro,tetin2, inparam);

		tetin2.save_nodes("../TESTOUT/testtetgen/rsvs_voro_grid");
		tetin2.save_poly("../TESTOUT/testtetgen/rsvs_voro_grid");

		tetrahedralize(inparam.command.c_str(), &tetin2, &tetout2);
		std::cout << "Finished the tettgen process " << std::endl;
		meshtet2 = TetgenOutput_TET2MESH(tetout2);
		tecout.PrintMesh(meshtet2);

	} catch (exception const& ex) {
		cerr << "Exception: " << ex.what() <<endl; 
		return 1;
	}
	return 0;
}
int tetcall_RSVSVORO()
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
	// std :: cin >> nPts;
	tecout.OpenFile(tecoutStr);

	for(auto i : numCells){
		nErrors += tetcall_RSVSVOROFunc(i, tecoutStr);
	}
	
	return nErrors;
}

int tetcall_RSVSVOROFunc(int nPts, const char* tecoutStr)
{
	std::vector<double> vecPts;
	mesh vosMesh, snakMesh, meshPts;
	tecplotfile tecout;

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

		RSVSVoronoiMesh(vecPts, vosMesh, snakMesh);
		tecout.PrintMesh(snakMesh,1,nPts);
		tecout.PrintMesh(vosMesh,2,nPts);
		meshPts.Init(nPts+8,0,0,0);
		meshPts.PopulateIndices();
		for (int i = 0; i < nPts+8; ++i)
		{	
			meshPts.verts[i].coord.assign(vecPts.begin()+i*4,
				vecPts.begin()+i*4+3);
		}
		meshPts.PrepareForUse();
		tecout.PrintMesh(meshPts,3,nPts,4);
		cout << "__________________________________________________" << endl;

	} catch (exception const& ex) {
		cerr << "Exception: " << ex.what() <<endl; 
		cerr << "for " << __PRETTY_FUNCTION__ << "with nPoints = " << nPts << endl;
		return 1;
	}
	return 0;
}
