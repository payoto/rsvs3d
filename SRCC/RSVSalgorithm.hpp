#ifndef RSVSALGORITHM_H_INCLUDED 
#define RSVSALGORITHM_H_INCLUDED 


//=================================
// forward declared dependencies
// 		class foo; //when you only need a pointer not the actual object
// 		and to avoid circular dependencies


//=================================
// included dependencies

#include "mesh.hpp"
#include "snake.hpp"

//==================================
// Code
// NOTE: function in a class definition are IMPLICITELY INLINED 
//       ie replaced by their code at compile time
	
void FindSpawnVerts(const mesh &meshin, vector<int> &vertList,
	vector<int> &voluOutList, int outerBorder=1);
void SpawnRSVS(snake &snakein, int outerBorder=1);
void RemoveSnakeInVolu(snake &snakein, vector<int> &voluInd, int outerBorder);

int Test_RSVSalgo_init();

#endif