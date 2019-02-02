#include <iostream>
#include <vector>
#include <array>

#include "tetgen_rsvs_api.hpp"
#include "tetgen.h"
#include "mesh.hpp"
#include "snake.hpp"
#include "snakevel.hpp"
#include "voxel.hpp"
#include "arraystructures.hpp"
#include "postprocessing.hpp"
// void tetrahedralize(char *switches, tetgenio *in, tetgenio *out, 
//                     tetgenio *addin, tetgenio *bgmin)


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
	snakeMesh.displight();
	
	voluMesh.PrepareForUse();
	voluMesh.displight();

	snakein.PrepareForUse();
	snakein.displight();
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
	triMesh.displight();
	triMesh.PrepareForUse();
	triMesh.OrderEdges();
	triMesh.TestConnectivityBiDir();

	// triMesh is the triangulation in mesh format
	// it is not cleaned up for tiny surfaces
}

mesh BuildDomain(const std::array<double,3> &lowerB, const std::array<double,3> &upperB){
	/*
	Builds a parralelipipede domain stretching from lowerB to upperB

	`lowerB` and `upperB` must be vectors of size 3.
	*/	
	mesh cube;
	int count; 

	if(lowerB.size()!=3 || upperB.size()!=3){
		std::cerr << "Error in " << __PRETTY_FUNCTION__ << " vectors must be"
			" of size 3" << std::endl;
		throw invalid_argument("Input vectors must be of size 3");
	}

	std::array<int, 3> dimGrid={1, 1, 1};
	BuildBlockGrid(dimGrid, cube);

	count = cube.verts.size();
	for (int i = 0; i < count; ++i) {
		for (int j = 0; j < 3; ++j)	{
			cube.verts[i].coord[j] = (cube.verts[i].coord[j]*(upperB[j]-lowerB[j])) + lowerB[j];
		}
	}

	return cube;
}


void MeshData2Tetgenio(const mesh &meshgeom, tetgen::io_safe &tetin,
	int facetOffset, int pointOffset, int pointMarker, const std::vector<double> &pointMtrList, const std::vector<double> &facetConstr, int facetConstrOffset){
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
		tetin.pointmarkerlist[i] = pointMarker;
	}

	// Set facet properties to the appropriate lists
	nFacetConstr = facetConstr.size();
	std::cout << "Number of facetconstraint " << nFacetConstr << std::endl;
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

	tetin.numberofpointmtrs = 0;
	tetin.numberoffacetconstraints = 2;


	tetin.allocate();

	// MeshData2Tetgenio(meshgeom, tetin, facetOffset, pointOffset, 
	// pointMarker, pointMtrList);
	// MeshData2Tetgenio(meshgeom, tetin, 0, 0, 1, {},{0.05*0.05},0);
	// MeshData2Tetgenio(meshdomain, tetin, meshgeom.surfs.size(), 
	// 	meshgeom.verts.size(), 2, {},{0.01},1);
	MeshData2Tetgenio(meshgeom, tetin, 0, 0, 1, {},{0.003},0);
	MeshData2Tetgenio(meshdomain, tetin, meshgeom.surfs.size(), 
		meshgeom.verts.size(), 2, {},{0.1},1);

}

HashedVector<int, int> GroupCloseSnaxels(const snake &snakein, double distTol){
	/*
	Merges snaxel which are `distTol` to the end of their carrying edge.
	*/
	std::vector<bool> isSnaxDone;
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


void TetgenInput_RSVS(const snake &snakein, tetgen::io_safe &tetin,
	const tetgen::apiparam &tetgenParam){
	/*
	Processes a snake object into a tetgen input object.

	This function processes snake objects into the safe tetgenio object.
	This can then used to generate a tetrahedral mesh or to save it as the
	input file.
	*/

	mesh meshdomain, meshgeom, meshtemp;
	triangulation triRSVS;
	int nHoles;
	// int nPtsHole, kk, count;
	std::vector<int>  holeIndices;

	tecplotfile tecout;

	// Cleanup the mesh
	tecout.OpenFile("../TESTOUT/rsvs_tetgen_mesh.plt");
	meshtemp=snakein.snakeconn;
	tecout.PrintMesh(meshtemp);
	meshtemp.displight();

	auto groupedVertices = GroupCloseSnaxels(snakein, tetgenParam.distanceTol);
	meshtemp.MergeGroupedVertices(groupedVertices);
	meshtemp.RemoveSingularConnectors();

	meshtemp.PrepareForUse();
	meshtemp.TestConnectivityBiDir();
	meshtemp.TightenConnectivity();
	meshtemp.OrderEdges();
	tecout.PrintMesh(meshtemp);
	meshtemp.displight();

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
	meshgeom.displight();
	meshgeom.TightenConnectivity();
	meshgeom.OrderEdges();

	tecout.PrintMesh(meshgeom);


	// Find number of holes
	holeIndices=FindHolesInSnake(snakein, groupedVertices);
	nHoles = holeIndices.size();
	// count = snakein.isMeshVertIn.size();
	// for (int i = 0; i < count; ++i) {
	// 	nHoles+=int(snakein.isMeshVertIn[i]);
	// }

	meshdomain = BuildDomain(tetgenParam.lowerB, tetgenParam.upperB);
	meshdomain.PrepareForUse();
	Mesh2Tetgenio(meshgeom, meshdomain, tetin, nHoles);

	std::cout<< std::endl << "Number of holes " << nHoles << std::endl;

	// Assign the holes
	
	for (int i = 0; i < nHoles; ++i){
		for (int j =0; j<3; ++j){
			tetin.holelist[i*3+j] = snakein.snakemesh->verts.isearch(
				holeIndices[i])->coord[j];
		}
	}


	tecout.PrintMesh(meshdomain);

	
}

int tetcall()
{
	tetgen::io_safe tetin, tetout;
	tetgen::apiparam inparam;
	mesh snakeMesh, voluMesh, triMesh;
	snake snakein;

	try {

		inparam.lowerB= {-2.0, -2.0,-2.0};
		inparam.upperB= {5.0, 5.0, 5.0};
		inparam.distanceTol = 1e-3;
		load_tetgen_testdata(snakeMesh, voluMesh, snakein, triMesh);
		TetgenInput_RSVS(snakein, tetin, inparam);


		tetin.save_nodes("rsvs_3cell_2body");
		tetin.save_poly("rsvs_3cell_2body");
		

		// Tetrahedralize the PLC. Switches are chosen to read a PLC (p),
		//   do quality mesh generation (q) with a specified quality bound
		//   (1.414), and apply a maximum volume constraint (a0.1).
		inparam.command = "pkqm"; 
		tetrahedralize(inparam.command.c_str(), &tetin, &tetout);

		// tetout.save_nodes("rsvsout_3cell_2body");
		// tetout.save_elements("rsvsout_3cell_2body");
		// tetout.save_faces("rsvsout_3cell_2body");
		std::cout << "Finished the tettgen process " << std::endl;
	} catch (exception const& ex) {
		cerr << "Exception: " << ex.what() <<endl; 
		return -1;
	}
	return 0;
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
		}
		if(this->facetlist[fIndex].numberofholes>0){
			this->facetlist[fIndex].holelist 
				= new REAL[this->facetlist[fIndex].numberofholes];
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

// Test code
int test_tetgenapi(){

	tetcall();

	return(0);
}
