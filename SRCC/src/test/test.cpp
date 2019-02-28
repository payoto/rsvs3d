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
#include "postprocessing.hpp"
#include "voxel.hpp"
#include "mesh.hpp"
#include "snake.hpp"
#include "snakeengine.hpp"
#include "meshrefinement.hpp"
#include "snakevel.hpp"
#include "RSVSalgorithm.hpp"
#ifdef TEST_ALL
#endif //TEST_ALL
#include "RSVSintegration.hpp"
#include "parameters.hpp"
#include "tetgen_rsvs_api.hpp"

int main(){
	customtest gridTest;

	#ifdef TEST_ALL
	// BASE templatess
	gridTest.Run(Test_ArrayStructures,"arraystructures");
	gridTest.Run(TestTemplate_ArrayStruct<vert>,"TestTemplate_ArrayStruct<vert>");
	gridTest.Run(TestTemplate_ArrayStruct<edge>,"TestTemplate_ArrayStruct<edge>");
	gridTest.Run(TestTemplate_ArrayStruct<surf>,"TestTemplate_ArrayStruct<surf>");
	gridTest.Run(TestTemplate_ArrayStruct<volu>,"TestTemplate_ArrayStruct<volu>");
	gridTest.Run(TestTemplate_ArrayStruct<snax>,"TestTemplate_ArrayStruct<snax>");
	gridTest.Run(TestTemplate_ArrayStruct<snaxedge>,"TestTemplate_ArrayStruct<snaxedge>");
	gridTest.Run(TestTemplate_ArrayStruct<snaxsurf>,"TestTemplate_ArrayStruct<snaxsurf>");
	// Building blocks
	gridTest.Run(Test_BuildBlockGrid_noout,"Voxel");
	gridTest.Run(Test_tecplotfile,"post-processing");
	gridTest.Run(Test_tecplotfile,"post-processing class");
	gridTest.Run(Test_snakeOrderEdges,"Snake Order error");
	gridTest.Run(Test_SnakeStructures,"Snake containers");
	gridTest.Run(Test_MeshOut,"Mesh output"); 
	gridTest.Run(Test_surfcentre,"test SurfCentroid"); 
	gridTest.Run(Test_MeshRefinement,"Multi-Level Meshes");
	gridTest.Run(Test_MeshOrient,"Output mesh orientation");
	gridTest.Run(Test_Crop,"test cropping of meshes"); 
	// Snakstruct 3D tests
	gridTest.Run(Test_RSVSalgo_init,"RSVS spawn");
	gridTest.Run(Test_snakeRSVS,"Snake RSVS");
	gridTest.Run(Test_snakeinit,"Snake rand velocity"); 
	gridTest.Run(Test_RSVSalgo,"Snake RSVS from spawn");
	gridTest.Run(Test_snakeRSVS_singlevol,"Snake RSVS single vol");
	gridTest.Run(Test_RSVSalgo_singlevol,"Snake RSVS single vol");
	// Snakstruct 2D tests
	gridTest.Run(Test_snakeinitflat,"Snake spawning 2D");
	// Parameter and JSON tests
	gridTest.Run(param::test::base,"parameter implementation");
	gridTest.Run(param::test::io,"parameter read write");
	gridTest.Run(param::test::ioflat,"parameter read write flat format");
	gridTest.Run(param::test::ipartialread,"parameter partial read write flat format");
	gridTest.Run(param::test::autoflat,"Algorithm for automatic determination of flat json");
	// RSVS and integration tests
	gridTest.Run(integrate::test::Prepare,"Mesh integration function");
	gridTest.Run(integrate::test::All,"Test full integration");
	gridTest.Run(Test_RSVSalgoflat,"RSVS 2D"); // Non working test - Maths not finished
	// Tetgen interface tests
	gridTest.Run(tetcall_CFD,"tegen API testing - CFD meshing"); 
	gridTest.Run(tetcall,"tegen API testing - RSVS meshing"); 
	#endif //TEST_ALL
	gridTest.Run(tetcall_RSVSVORO,"tegen API testing - Voro to RSVS"); // working test
	gridTest.Run(Test_Crop,"test cropping of meshes"); 

	
	gridTest.PrintSummary();

	return(0);
}


