
//===============================================
// Include Guards
#ifndef SNAKSTRUCT_H_INCLUDED
#define SNAKSTRUCT_H_INCLUDED

//===============================================
// Levels of debuging Guards
#ifdef DEBUGLVL2 // All Debugging calls
#define DEBUGLVL1


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


class ConnecRemv {
public:
	int keepind;
	int typeobj;
	vector<int> rmvind;

};



// Function prototypes
void SpawnAtVertex(snake &snakein,int indVert);
void SpawnAtVertexVert(snake& newsnake, int nVert,int indVert, int subVert, 
	const vector<int> &surfInds, const vector<int> &edgeInds, const vector<int> &edgeSubs,
	unordered_multimap<int,int> &hashSurfInds);
void SpawnAtVertexEdge(snake& newsnake,int nEdge,const vector<int> &surfInds,
	const vector<int> &edgeInds,const vector<int> &voluInds,const vector<int> &surfSubs,
	unordered_multimap<int,int> &hashEdgeInds, unordered_multimap<int,int> &hashVoluInds);

void SpawnAtVertexSurf3D(snake& newsnake,int nSurf,const vector<int> &surfInds,
 const vector<int> &voluInds,const vector<int> &voluSubs,unordered_multimap<int,int> &hashSurfInds);
void SpawnAtVertexSurf2D(snake& newsnake,int nEdge,const vector<int> &voluInds);
void SpawnAtVertexVolu(snake& newsnake, int nSurf);
void MergeAllContactVertices(snake &fullsnake, vector<int> &isImpact);

void SpawnArrivedSnaxels(snake &fullsnake, const vector<int> &isImpact);
void SpawnArrivedSnaxelsDir(const snake &fullsnake,snake &partSnake,const  vector<int> &isImpact,int dir);
// Test Function prototypes
#endif //SNAKSTRUCT_H_INCLUDED