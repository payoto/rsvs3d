#ifdef DBG_MEMLEAK
#define _CRTDBG_MAP_ALLOC  
#include <stdlib.h>  
#include <crtdbg.h>
#endif //DBG_MEMLEAK

#include <iostream>
#include <cstdlib>
#include <string>
#include <functional>

#include "test.hpp"
#include "arraystructures.hpp"

#ifdef TEST_ALL
#include "voxel.hpp"
#include "postprocessing.hpp"
#endif //TEST_ALL
#include "snakstruct.hpp"
#include "snakeengine.hpp"

int main(){
	customtest gridTest;
	#ifdef TEST_ALL
	gridTest.Run(Test_ArrayStructures,"arraystructures");
	
	gridTest.Run(Test_BuildBlockGrid_noout,"Voxel");

	gridTest.Run(Test_tecplotfile,"post-processing");

	gridTest.Run(Test_tecplotfile,"post-processing class");

	gridTest.Run(Test_MeshOut,"Mesh output");
	gridTest.Run(TestTemplate_ArrayStruct<vert>,"TestTemplate_ArrayStruct<vert>");
	gridTest.Run(TestTemplate_ArrayStruct<edge>,"TestTemplate_ArrayStruct<edge>");
	gridTest.Run(TestTemplate_ArrayStruct<surf>,"TestTemplate_ArrayStruct<surf>");
	gridTest.Run(TestTemplate_ArrayStruct<volu>,"TestTemplate_ArrayStruct<volu>");
	gridTest.Run(TestTemplate_ArrayStruct<snax>,"TestTemplate_ArrayStruct<snax>");
	gridTest.Run(TestTemplate_ArrayStruct<snaxedge>,"TestTemplate_ArrayStruct<snaxedge>");
	gridTest.Run(TestTemplate_ArrayStruct<snaxsurf>,"TestTemplate_ArrayStruct<snaxsurf>");
	gridTest.Run(Test_snakeOrderEdges,"Snake Order error");
	gridTest.Run(Test_SnakeStructures,"Snake containers");
	#endif //TEST_ALL

	gridTest.Run(Test_snakeinit,"Snake spawning");
	gridTest.Run(Test_snakeinitflat,"Snake spawning 2D");


	gridTest.PrintSummary();



	return(0);
}


