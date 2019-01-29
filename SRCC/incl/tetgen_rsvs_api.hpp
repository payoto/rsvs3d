#ifndef TETGEN_RSVS_API_H_INCLUDED 
#define TETGEN_RSVS_API_H_INCLUDED 



//=================================
// forward declared dependencies
// 		class foo; //when you only need a pointer not the actual object
// 		and to avoid circular dependencies


//=================================
// included dependencies

#include "tetgen.h"

//=================================
// Code
// 
// This file defines parameters to be used in other parts of the 
// RSVS snaking framework. Default values are defined in "parameters.cpp"
//
// Substructure names are all 4-5 letters

class tetgenio_safe : public tetgenio {

public:
	// REAL *pointlist;
	// REAL *pointattributelist;
	// REAL *pointmtrlist;
	// int  *pointmarkerlist;

	using tetgenio::pointlist;
	using tetgenio::pointattributelist;
	using tetgenio::pointmtrlist;
	using tetgenio::pointmarkerlist;

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
	using tetgenio::edgeadjtetlist;

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
};


// Test functions

int test_tetgenapi();

#endif