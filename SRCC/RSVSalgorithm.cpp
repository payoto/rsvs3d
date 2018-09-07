#include <iostream>
#include <cstdlib>

#include "RSVSalgorithm.hpp"
#include "snake.hpp"
#include "snakeengine.hpp"

void FindSpawnVerts(const mesh &meshin, vector<int> &vertList, int outerBorder){
	// Function which identifies spawn points
	// Spawn points are:
	//  - Any point part of a cell which touches the void
	//    Put that is not part of a surface that is on the void 
	//    itself.
	//  - Points That are on the border of two cells one with some VOS 
	//    One without.

	int ii, ni, jj, nj;
	vector<int> offBorderVert, internalSurf;


	if (outerBorder==1){ // spawn at outer borders
		meshin.GetOffBorderVert(offBorderVert,1);
		meshin.SurfOnParentBound(internalSurf,false,true);
	} else if (outerBorder==0){
		meshin.GetOffBorderVert(offBorderVert,0);
		meshin.SurfOnParentBound(internalSurf,false,false);
	} else {
		throw invalid_argument("outerBorder has an unknown value");
	}
	ni=meshin.verts.size();
	vertList.reserve(ni);
	vertList.clear();
	ni=internalSurf.size();
	for (ii=0; ii< ni; ++ii){
		nj=meshin.surfs(ii)->edgeind.size();
		for(jj=0; jj<nj; jj++){
			vertList.push_back(meshin.edges.isearch(meshin.surfs(ii)->edgeind[jj])->vertind[0]);
			vertList.push_back(meshin.edges.isearch(meshin.surfs(ii)->edgeind[jj])->vertind[1]);
		}
	}
	ni=offBorderVert.size();
	for (ii=0; ii< ni; ++ii){
		vertList.push_back(offBorderVert[ii]);
	}

	sort(vertList);
	unique(vertList);

}

void SpawnRSVS(snake &snakein){
	// Function which handles
	// - spawning 
	// - growing the snake
	// - Identifying snaxels to remove
	// - update connectivity:
	//     + find snaxsurfs in invalid snake volus(surfs)
	//     + find snaxedges in these surfs
	//     + find snaxels
	//     * invalid snakevolus=[border cell, empty cell]
	int ii,ni;
	vector<int> vertSpawn;

	FindSpawnVerts(*(snakein.snakemesh), vertSpawn,1);
	ni=vertSpawn.size();
	for(ii=0; ii< ni; ++ii){
		SpawnAtVertex(snakein, vertSpawn[ii]);
	}
}