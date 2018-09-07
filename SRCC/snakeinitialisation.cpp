#include <iostream>
#include <cstdlib>

#include "snakeinitialisation.hpp"



#include "snake.hpp"

void FindSpawnVerts(const mesh &meshin, vector<int> &vertList){
	// Function which identifies spawn points
	// Spawn points are:
	//  - Any point part of a cell which touches the void
	//    Put that is not part of a surface that is on the void 
	//    itself.
	//  - Points That are on the border of two cells one with some VOS 
	//    One without.


}


// Function which handles
// - spawning 
// - growing the snake
// - Identifying snaxels to remove
// - update connectivity:
//     + find snaxsurfs in invalid snake volus(surfs)
//     + find snaxedges in these surfs
//     + find snaxels
//     * invalid snakevolus=[border cell, empty cell]