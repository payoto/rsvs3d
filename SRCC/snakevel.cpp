
#include <iostream>
#include <cstdlib>
#include "arraystructures.hpp"
#include "snakstruct.hpp" 

void CalculateSnakeVel(snake &snakein){

	int ii=0;

	for (ii=0;ii<int(snakein.snaxs.size());++ii){
		if (snakein.snaxs(ii)->isfreeze==1){
			snakein.snaxs[ii].v=(0.5-snakein.snaxs[ii].d)*0.3;
			snakein.snaxs[ii].isfreeze=0;
		}
		snakein.snaxs[ii].v=(double(rand()%1001)/1000.0+0.5)*snakein.snaxs[ii].v;
	}
}