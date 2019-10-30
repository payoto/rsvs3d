/**
 * Interface between the RSVS project and [tetgen](tetgen.org).
 *  
 *@file
 */

#ifndef TETGEN_RSVS_API_H_INCLUDED 
#define TETGEN_RSVS_API_H_INCLUDED 



//=================================
// forward declared dependencies
// 		class foo; //when you only need a pointer not the actual object
// 		and to avoid circular dependencies


//=================================
// included dependencies

#include <array>
#include <vector>
#include <string>
#include <algorithm>

#include "tetgen.h"
#include "mesh.hpp"


//-------------------------------------------------------------------------
// Tetgen namespace containing the interfacing functions
//

namespace tetgen {

	/**
	 * Type defining domain boundaries.
	 * 
	 * Simple short hand for a matrix of 2*3 doubles.
	 */	
	typedef std::array<std::array<double, 3>, 2> dombounds;

	/**
	 * @brief      Class for memory safe interface with tetgen.h
	 * 
	 * This class provides a method called `allocate` which allocates the
	 * memory for the io arrays using the new command. Command `deallocate`
	 * can be used to free the memory before destruction, or otherwise it is
	 * called uppon when object goes out of scope.
	 */
	class io_safe : public tetgenio {
		// Adds a semblance of prebuilt allocation functions
		// to the io class build into tetgen
		public:
		// REAL *pointlist;
		// REAL *pointattributelist;
		// REAL *pointmtrlist;
		// int  *pointmarkerlist;

		using tetgenio::pointlist;

		using tetgenio::pointattributelist;
		using tetgenio::pointmtrlist;
		using tetgenio::pointmarkerlist;
		using tetgenio::numberofpointmtrs;

		// int  *tetrahedronlist;
		// REAL *tetrahedronattributelist;
		// REAL *tetrahedronvolumelist;
		// int  *neighborlist;
		using tetgenio::tetrahedronlist;
		using tetgenio::tetrahedronattributelist;
		using tetgenio::tetrahedronvolumelist;
		using tetgenio::neighborlist;

		// tetgenio::facet *facetlist;
		// int *facetmarkerlist;
		using tetgenio::facetlist;
		using tetgenio::facetmarkerlist;
		// REAL *facetconstraintlist;
	    // int numberoffacetconstraints;
	    using tetgenio::facetconstraintlist;
	    using tetgenio::numberoffacetconstraints;
		// 	REAL *holelist;
		// 	REAL *regionlist;
		using tetgenio::holelist;
		using tetgenio::regionlist;

		// REAL *segmentconstraintlist;
		using tetgenio::segmentconstraintlist;

		// int *edgelist;
		// int *edgemarkerlist;
		// int *o2edgelist;
		// int *edgeadjtetlist;
		using tetgenio::edgelist;
		using tetgenio::edgemarkerlist;
		using tetgenio::o2edgelist;
		using tetgenio::edge2tetlist;
		using tetgenio::face2edgelist;
		using tetgenio::face2tetlist;

		// REAL *vpointlist;
		// tetgenio::voroedge *vedgelist;
		// tetgenio::vorofacet *vfacetlist;
		// int **vcelllist;
		using tetgenio::vpointlist;
		using tetgenio::vedgelist;
		using tetgenio::vfacetlist;
		using tetgenio::vcelllist;

		using tetgenio::initialize;

		
		void allocate();
		void allocatefacet(int fIndex);
		void allocatefacet(int fIndex, int numPoly);
		void allocatefacetpolygon(int fIndex, int pIndex);
		void allocatefacetpolygon(int fIndex, int pIndex, int numVerts);

		void SpecifyTetPointMetric(int startPnt, int numPnt, 
			const std::vector<double> &mtrs);
		void SpecifyIndividualTetPointMetric(int startPnt, int numPnt, 
			const std::vector<double> &mtrs);
		void SpecifyTetFacetMetric(int startPnt, int numPnt, 
			int marker);

		void displaystats();

		io_safe(){initialize();}
	};

	class apiparam {
	public:
		/// Lower domain bound
		std::array<double,3> lowerB;
		/// Upper domain bound 
		std::array<double,3> upperB; 
		/// Controls the surface edgelengths in CFD in the order:
		/// {point of lowest curvature, point of highest curvature}
		std::array<double, 2> surfedgelengths;
		int curvatureSmoothing;
		/// Controls the edgelengths at regular intervals 
		std::vector<double> edgelengths;
		/// Distance tolerance 
		double distanceTol; 
		// mesh the inside of the geometry? (or the outside)
		bool generateMeshInside;

		// Commmand line string to be passed to tetgen
		std::string command;

		void ReadJsonString(const std::string &jsonStr);
		apiparam(){
			this->lowerB = {0.0, 0.0, 0.0};
			this->upperB = {1.0, 1.0, 1.0};
			this->surfedgelengths = {0.02, 0.005};
			this->curvatureSmoothing = 4;
			this->edgelengths = {0.03, 1.0};
			this->distanceTol = 0.3;
			this->generateMeshInside = false;
			this->command = "";
		}
		/**
		 * @brief      Constructs the object from a json string.
		 *
		 * @param[in]  jsonStr  The json string
		 */
		apiparam(const std::string &jsonStr) : apiparam() {
			this->ReadJsonString(jsonStr);
		}
	};


	std::vector<int> RSVSVoronoiMesh(const std::vector<double> &vecPts, 
		mesh &vosMesh, mesh &snakMesh,
		tetgen::apiparam &inparam);
	
	void SnakeToSU2(const snake &snakein, const std::string &fileName,
		tetgen::apiparam &inparam);

	tetgenmesh* rsvstetrahedralize(const char *switches, tetgen::io_safe *in, 
		tetgen::io_safe *out, tetgen::io_safe *addin = NULL,
		tetgen::io_safe *bgmin = NULL);

	namespace input {
		void POINTGRIDS(const mesh &meshdomain, 
			tetgen::io_safe &tetin, const tetgen::apiparam &tetgenParam,
			bool generateVoroBound=false);
		void RSVSGRIDS(const mesh &meshdomain,
			tetgen::io_safe &tetin, const tetgen::apiparam &tetgenParam);
		void RSVSGRIDS(const mesh &meshdomain, 
			const mesh &meshboundary, tetgen::io_safe &tetin,
			const tetgen::apiparam &tetgenParam);
		void RSVS2CFD(const snake &snakein,
			tetgen::io_safe &tetin,	const tetgen::apiparam &tetgenParam);
	}
	namespace output {
		mesh VORO2MESH(tetgen::io_safe &tetout);
		void SU2(const char* fileName, const tetgenio &tetout);
		dombounds GetBoundBox(io_safe &tetout);
		mesh TET2MESH(tetgen::io_safe &tetout);
	}
	namespace internal {
		void CloseVoronoiMesh(mesh &meshout, tetgen::io_safe &tetout, 
			std::vector<int> &rayEdges, int DEINCR, 
			tetgen::dombounds boundBox);
		template<class T, class V>
		double ProjectRay(int count, const tetgen::dombounds &boundBox,
			const T &dir, const V &orig, double minDist=0.0);
		void MeshData2Tetgenio(const mesh &meshgeom, tetgen::io_safe &tetin,
			int facetOffset, int pointOffset, int pointMarker, 
			const std::vector<double> &pointMtrList, 
			const std::vector<double> &facetConstr,
			int facetConstrOffset);
		void Mesh2Tetgenio(const mesh &meshgeom, const mesh &meshdomain,
			tetgen::io_safe &tetin, int numHoles);
		void Mesh2TetgenioPoints(const mesh &meshgeom, 
			const mesh &meshdomain,	tetgen::io_safe &tetin);
		void PointCurvature2Metric(std::vector<double> &vertCurvature, 
			const tetgen::apiparam &inparam);
	}
	namespace voronoi {
		void GenerateInternalPoints(const mesh &meshin, int nLevels, 
			tetgen::io_safe &tetinPts);
		std::vector<bool> Points2VoroAndTetmesh(const std::vector<double> &vecPts,
			mesh &voroMesh, mesh &tetMesh, const tetgen::apiparam &inparam);
		std::vector<bool> BoundaryFacesFromPoints(const mesh &meshin,
			const std::vector<int> &boundaryPts);
	}

	namespace test {
		void LoadData(mesh &snakeMesh,
			mesh &voluMesh, snake &snakein, mesh &triMesh);
		int api();
		int call();
		int CFD();
		int RSVSVORO();
		int RSVSVORO_Contain();
		int RSVSVOROFunc(const std::vector<double> &vecPts={}, double distanceTol=0.26,
			const char* tecoutStr="../TESTOUT/rsvs_voro.plt");
		int RSVSVOROFunc_contain(int nPts=0, double distanceTol=0.26, 
			const char* tecoutStr="../TESTOUT/rsvs_voro_contain.plt");
	}
}


/**
 * @brief      Project voronoi diagram rays to the bounding Box
 *
 * @param[in]  count     number of coordinates
 * @param[in]  boundBox  The bounds of the domain (array<array<double,3>,2>)
 * @param[in]  dir       vector with direction (pointing in)
 * @param[in]  orig      The origin of the ray
 * @param[in]  minDist   The minimum allowable stretch distance
 *
 * @tparam     T         type of `dir`: an iterable of size `count` 
 * @tparam     V         type of `orig`: an iterable of size `count`
 *
 * @return     Distance along the ray at which the boundBox is encountered.
 */
template<class T, class V>
double tetgen::internal::ProjectRay(int count, const tetgen::dombounds &boundBox,
	const T &dir, const V &orig, double minDist){
	/*
	Calculates the distance to project a single ray to reach the bounding box

	It is templated to accept any types of containers with `count` elements, 
	accessible by operator[];
	*/

	double l=-INFINITY;

	for (int i = 0; i < count; ++i)
	{
		l = std::max(l,
			std::min(
				(boundBox[0][i]-orig[i])/dir[i] ,
			(boundBox[1][i]-orig[i])/dir[i] ));
	}
	return std::min(l, minDist);
}
int Test_RSVSvoro_init();
#endif