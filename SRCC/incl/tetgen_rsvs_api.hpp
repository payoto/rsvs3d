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
#include "tetgen.h"
//=================================
// Code
// 
// This file defines parameters to be used in other parts of the 
// RSVS snaking framework. Default values are defined in "parameters.cpp"
//
// Substructure names are all 4-5 letters
namespace tetgen {

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

		// using tetgenio::~tetgenio;

		void allocate();
		void allocatefacet(int fIndex);
		void allocatefacet(int fIndex, int numPoly);
		void allocatefacetpolygon(int fIndex, int pIndex);
		void allocatefacetpolygon(int fIndex, int pIndex, int numVerts);
		void SpecifyTetPointMetric(int startPnt, int numPnt, 
			const std::vector<double> &mtrs);
		// void deallocate();
		// ~io_safe(){
		// 	this->deallocate();
		// }
	};

	class apiparam {
	public:
		// Builds a cuboid bound
		std::array<double,3> lowerB={-3.0,-3.0,-3.0}; // Lower domain bound
		std::array<double,3> upperB={4.0,4.0,4.0}; // Upper domain bound
		// Controls the edgelengths 
		std::vector<double> edgelengths={0.03,1.0}; 
		double distanceTol = 0.3; // Distance tolerance
		// mesh the inside of the geometry? (or the outside)
		bool generateMeshInside = false;

		// Commmand line string to be run
		std::string command;
	};
}

// Test functions

int test_tetgenapi();

#endif