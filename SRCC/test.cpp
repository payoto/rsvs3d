#include <iostream>
#include <cstdlib>
#include <string>
#include <functional>

#include "test.hpp"
#include "arraystructures.hpp"
#include "voxel.hpp"
#include "postprocessing.hpp"
#include "snakstruct.hpp"

int main(){
	customtest gridTest;

	gridTest.Run(Test_ArrayStructures,"arraystructures");
	
	gridTest.Run(Test_BuildBlockGrid_noout,"Voxel");

	gridTest.Run(Test_tecplotfile,"post-processing");

	gridTest.Run(Test_tecplotfile,"post-processing class");

	gridTest.Run(Test_MeshOut,"Mesh output");
	
	gridTest.Run(Test_SnakeStructures,"Snake containers");

	gridTest.PrintSummary();
	return(0);
}


