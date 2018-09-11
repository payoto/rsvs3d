#include <iostream>
#include <cstdlib>

#include "RSVSalgorithm.hpp"
#include "snake.hpp"
#include "snakeengine.hpp"

void FindSpawnVerts(const mesh &meshin, vector<int> &vertList,
	vector<int> &voluOutList, int outerBorder){
	// Function which identifies spawn points
	// Spawn points are:
	//  - Any point part of a cell which touches the void
	//    Put that is not part of a surface that is on the void 
	//    itself.
	//  - Points That are on the border of two cells one with some VOS 
	//    One without.

	int ii, ni, jj, nj, surfSurb, vert_i;
	vector<int> offBorderVert, internalSurf, voluIndIntern;


	if (outerBorder==1){ // spawn at outer borders
		meshin.GetOffBorderVert(offBorderVert,voluOutList,1);
		meshin.SurfOnParentBound(internalSurf,voluIndIntern,false,true);
	} else if (outerBorder==0){
		meshin.GetOffBorderVert(offBorderVert,voluOutList,0);
		meshin.SurfOnParentBound(internalSurf,voluIndIntern,false,false);
	} else {
		throw invalid_argument("outerBorder has an unknown value");
	}
	ni=meshin.verts.size();
	vertList.reserve(ni);
	vertList.clear();
	ni=internalSurf.size();

	cout << "vertices to output internalSurf " << ni << endl;
	for (ii=0; ii< ni; ++ii){
		surfSurb=meshin.surfs.find(internalSurf[ii]);
		nj=meshin.surfs(surfSurb)->edgeind.size();
		for(jj=0; jj<nj; jj++){
			vert_i=meshin.edges.isearch(meshin.surfs(surfSurb)->edgeind[jj])->vertind[0];
			if(!meshin.verts.isearch(vert_i)->isBorder){
				vertList.push_back(meshin.verts.isearch(vert_i)->index);
			}
			vert_i=meshin.edges.isearch(meshin.surfs(surfSurb)->edgeind[jj])->vertind[1];
			if(!meshin.verts.isearch(vert_i)->isBorder){
				vertList.push_back(meshin.verts.isearch(vert_i)->index);
			}
		}
	}

	sort(offBorderVert);
	unique(offBorderVert);
	ni=offBorderVert.size();
	cout << "vertices to output offBorderVert " << ni << endl;
	for (ii=0; ii< ni; ++ii){
		vertList.push_back(offBorderVert[ii]);
	}

	sort(vertList);
	unique(vertList);

	ni=voluIndIntern.size();
	for (ii=0; ii< ni; ++ii){
		voluOutList.push_back(voluIndIntern[ii]);
	}

	sort(voluOutList);
	unique(voluOutList);



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
	vector<int> voluSnaxDelete;
	vector<int> isImpact;
	vector<double> dt;

	FindSpawnVerts(*(snakein.snakemesh), vertSpawn,voluSnaxDelete,1);
	ni=vertSpawn.size();
	cout << "vertices to output " << ni << endl;
	for(ii=0; ii< ni; ++ii){
		SpawnAtVertex(snakein, vertSpawn[ii]);
	}
	// Move to half distances
	ni=snakein.snaxs.size();
	for(ii=0; ii<ni; ++ii){
		snakein.snaxs[ii].v=1.0;
	}
	snakein.CalculateTimeStep(dt,0.6);
	snakein.UpdateDistance(dt);
	snakein.UpdateCoord();
	snakein.PrepareForUse();
	// detect and remove impacts
	snakein.SnaxImpactDetection(isImpact);
	MergeAllContactVertices(snakein, isImpact);
	snakein.PrepareForUse();
	CleanupSnakeConnec(snakein);

	// Remove one of the 'snakes'
}