/**
 * Functions which are part of the RSVS algorithm but not core to the snaking
 * process.
 *
 *@file
 */

#ifndef RSVSALGORITHM_H_INCLUDED
#define RSVSALGORITHM_H_INCLUDED

//=================================
// forward declared dependencies
// 		class foo; //when you only need a pointer not the actual object
// 		and to avoid circular dependencies

class mesh;
class snake;

//=================================
// included dependencies

#include <vector>

//==================================
// Code
// NOTE: function in a class definition are IMPLICITELY INLINED
//       ie replaced by their code at compile time

std::vector<int> FindSpawnVerts(const mesh &meshin, std::vector<int> &vertList, std::vector<int> &voluOutList,
                                int outerBorder = 1);
void SpawnRSVS(snake &snakein, int outerBorder = 1);
void RemoveSnakeInVolu(snake &snakein, std::vector<int> &voluInd, int outerBorder);
void RemoveSnakeInSurf(snake &snakein, std::vector<int> &voluInd, int outerBorder);
void SpawnSnakeAndMove(snake &snakein, std::vector<int> vertSpawn);

int Test_RSVSalgo_init();
int Test_RSVSalgo();
int Test_RSVSalgoflat();

#endif
