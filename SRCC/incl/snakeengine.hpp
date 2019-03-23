/**
 * Functions needed to evolve the r-surface snake.
 * 
 *@file
 */

//===============================================
// Include Guards
#ifndef SNAKEENGINE_H_INCLUDED
#define SNAKEENGINE_H_INCLUDED

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
#include "snake.hpp"
#include "postprocessing.hpp"

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

void SpawnAtVertexSurf3D(snake& newsnake,int nSurf,const vector<int> &surfInds,
	const vector<int> &voluInds,const vector<int> &voluSubs,unordered_multimap<int,int> &hashSurfInds);
void SpawnAtVertexSurf2D(snake& newsnake,int nEdge,const vector<int> &voluInds);
void SpawnAtVertexVolu(snake& newsnake, int nSurf);
void MergeAllContactVertices(snake &fullsnake, vector<int> &isImpact);

void SpawnArrivedSnaxels(snake &fullsnake, const vector<int> &isImpact);
void SpawnArrivedSnaxelsDir(snake &fullsnake,snake &partSnake,const  vector<int> &isImpact,int dir, HashedVector<int,int> &vertNoSpawn);

//
void MergeCleanSnake(snake &fullsnake, vector<int> &isImpact);
void CleanupSnakeConnec(snake &snakein);
void IdentifyMergEdgeSameSurfConnec(const snake &snakein, vector<ConnecRemv> &connecEdit);
void IndentifyEdgeSameSurf(const snake &snakein,int currSub, int &stepCheck,vector<int> &tempSub,
	vector<int> &tempSub2,vector<int> &tempSub3,HashedVector<int,int> &tempIndHash,
	HashedVector<int,int> &edge2Surf,vector<int> tempCount);
void IdentifyMergEdgeConnec(const snake &snakein, vector<ConnecRemv> &connecEdit);
void IdentifyMergeEdgeGeneral(const snake &snakein, vector<bool> &isObjDone,
	vector<ConnecRemv> &connecEdit, ConnecRemv &tempConnec,  ConnecRemv &tempConnec2,
	vector<int> &tempSub,vector<int> &tempSub2, vector<int> &tempCount, 
	HashedVector<int,int> &tempIndHash) ;
void IdentifyMergeEdgeGeneralChain(const snake &snakein, vector<bool> &isObjDone,vector<ConnecRemv> &connecEdit, ConnecRemv &tempConnec,  ConnecRemv &tempConnec2,vector<int> &tempSub,vector<int> &tempSub2, vector<int> &tempCount, HashedVector<int,int> &tempIndHash, int jjStart);

void IdentifyMergSurfConnec(const snake &snakein, vector<ConnecRemv> &connecEdit);
void IdentifyMergeSurfGeneral(const snake &snakein, vector<bool> &isObjDone,vector<ConnecRemv> &connecEdit, 
	ConnecRemv &tempConnec,vector<int> &tempSub,vector<int> &tempSub2,
	vector<int> &tempCount,HashedVector<int,int> &edge2Surf, HashedVector<int,int> &tempIndHash) ;
void IdentifyMergeSurfRecursive(const snake &snakein,vector<bool> &isObjDone, vector<int> &tempCount,const HashedVector<int,int> &edge2Surf, const HashedVector<int,int> &tempIndHash, ConnecRemv &tempConnec, const vector<int> &tempSub, const vector<int> &tempSub2, int excludeSub);
void ModifyMergVoluConnec(snake &snakein, vector<ConnecRemv> &connecEdit, 
	const vector<int> &indRmvVert);

void ModifyMergSurf2DConnec(snake &snakein, vector<ConnecRemv> &connecEdit);
void SnaxEdgeConnecDetection(snake &snakein, vector<ConnecRemv> &connecEdit);
void SnaxNoConnecDetection(const mesh &snakeconn, vector<ConnecRemv> &connecEdit);
void dispconnrmv(vector<ConnecRemv> conn);

void CheckSnakeRemovalsVert(const snake &snakein, const vector<int> &indRmvVert);
void CheckSnakeRemovalsEdge(const snake &snakein, const vector<int> &indRmvEdge);
void CheckSnakeRemovalsSurf(const snake &snakein, const vector<int> &indRmvSurf);
void CheckSnakeRemovalsVolu(const snake &snakein, const vector<int> &indRmvVolu);
// Test Function prototypes
#endif //SNAKSTRUCT_H_INCLUDED