
//===============================================
// Include Guards
#ifndef SNAKSTRUCT_H_INCLUDED
#define SNAKSTRUCT_H_INCLUDED

//===============================================
// Levels of debuging Guards
#ifdef DEBUGLVL2 // All Debugging calls
#define DEBUGLVL1

#define 
#endif

#ifdef DEBUGLVL1 // Debugging of new features.
#define TEST_SNAKEENGINE
#endif

//=================================
// forward declared dependencies
// 		class foo; //when you only need a pointer not the actual object
// 		and to avoid circular dependencies


//=================================
// included dependencies
#include <iostream>
#include <vector>
#include <unordered_map>

#include "arraystructures.hpp"
#include "snakstruct.hpp"

//==================================
// Code
// NOTE: function in a class definition are IMPLICITELY INLINED 
//       ie replaced by their code at compile time
using namespace std;


// Function prototypes
void SpawnAtVertex(snake &snakein,int indVert);
void SpawnAtVertexVert(snake& newsnake, int nVert,int indVert, int subVert, 
	const vector<int> &surfInds, const vector<int> &edgeInds, const vector<int> &edgeSubs,
	unordered_multimap<int,int> &hashSurfInds);
void SpawnAtVertexEdge(snake& newsnake,int nEdge,const vector<int> &surfInds,
	const vector<int> &edgeInds,const vector<int> &voluInds,const vector<int> &surfSubs,
	unordered_multimap<int,int> &hashEdgeInds, unordered_multimap<int,int> &hashVoluInds);

void SpawnAtVertexSurf(snake& newsnake,int nSurf,const vector<int> &surfInds,
 const vector<int> &voluInds,const vector<int> &voluSubs,unordered_multimap<int,int> &hashSurfInds);
void SpawnAtVertexVolu(snake& newsnake, int nSurf);

// Test Function prototypes
int Test_SnakeStructures();
int Test_coordvec();
int Test_snax();
int Test_snaxedge();
int Test_snake();
int Test_snakeinit();
#endif //SNAKSTRUCT_H_INCLUDED