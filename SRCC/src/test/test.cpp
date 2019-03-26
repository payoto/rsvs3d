#include <iostream>
#include <cstdlib>
#include <string>
#include <functional>
#include <sstream>

#include "test.hpp"

#include "arraystructures.hpp"
#include "postprocessing.hpp"
#include "voxel.hpp"
#include "mesh.hpp"
#include "snake.hpp"
#include "snakeengine.hpp"
#include "meshrefinement.hpp"
#include "triangulate.hpp"
#include "RSVSalgorithm.hpp"
#include "RSVSintegration.hpp"
#include "parameters.hpp"
#include "tetgenrsvs.hpp"

using namespace rsvstest;

int rsvstest::maintest(){ 
	customtest gridTest("3-D RSVS tests - passing");

	// BASE templatess 
	gridTest.RunSilent(Test_ArrayStructures,"arraystructures");
	gridTest.RunSilent(TestTemplate_ArrayStruct<vert>,
		"TestTemplate_ArrayStruct<vert>");
	gridTest.RunSilent(TestTemplate_ArrayStruct<edge>,
		"TestTemplate_ArrayStruct<edge>");
	gridTest.RunSilent(TestTemplate_ArrayStruct<surf>,
		"TestTemplate_ArrayStruct<surf>");
	gridTest.RunSilent(TestTemplate_ArrayStruct<volu>,
		"TestTemplate_ArrayStruct<volu>");
	gridTest.RunSilent(TestTemplate_ArrayStruct<snax>,
		"TestTemplate_ArrayStruct<snax>");
	gridTest.RunSilent(TestTemplate_ArrayStruct<snaxedge>,
		"TestTemplate_ArrayStruct<snaxedge>");
	gridTest.RunSilent(TestTemplate_ArrayStruct<snaxsurf>,
		"TestTemplate_ArrayStruct<snaxsurf>");
	// Building blocks
	gridTest.RunSilent(Test_BuildBlockGrid_noout,"Voxel");
	gridTest.RunSilent(Test_tecplotfile,"post-processing");
	gridTest.RunSilent(Test_tecplotfile,"post-processing class");
	gridTest.RunSilent(Test_snakeOrderEdges,"Snake Order error");
	gridTest.RunSilent(Test_SnakeStructures,"Snake containers");
	gridTest.RunSilent(Test_MeshOut,"Mesh output"); 
	gridTest.RunSilent(Test_surfcentre,"test SurfCentroid"); 
	gridTest.RunSilent(Test_MeshRefinement,"Multi-Level Meshes");
	gridTest.RunSilent(Test_MeshOrient,"Output mesh orientation");
	gridTest.RunSilent(Test_Crop,"test cropping of meshes"); 
	// Snakstruct 3D tests
	gridTest.RunSilent(Test_RSVSalgo_init,"RSVS spawn");
	gridTest.RunSilent(Test_RSVSvoro_init,"Snake spawning voronoi"); 
	gridTest.RunSilent(Test_snakeinit,"Snake rand velocity"); 
	gridTest.RunSilent(Test_RSVSalgo,"Snake RSVS from spawn");
	gridTest.RunSilent(Test_snakeRSVS,"Snake RSVS");
	gridTest.RunSilent(Test_snakeRSVS_singlevol,"Snake RSVS single vol");
	gridTest.RunSilent(Test_RSVSalgo_singlevol,"Snake RSVS algorithm single vol");
	// Parameter and JSON tests
	gridTest.RunSilent(param::test::base,"parameter implementation");
	gridTest.RunSilent(param::test::symmetry,"Test json internal symmetry");
	gridTest.RunSilent(param::test::io,"parameter read write");
	gridTest.RunSilent(param::test::ioflat,"parameter read write flat format");
	gridTest.RunSilent(param::test::ipartialread,"parameter partial read"
		" write flat format");
	gridTest.RunSilent(param::test::autoflat,"Algorithm for automatic "
		"determination of flat json");
	gridTest.RunSilent(tetgen::test::api,"tegen test api parameter class");
	// RSVS and integration tests
	gridTest.RunSilent(integrate::test::Prepare,"Mesh integration function");
	gridTest.RunSilent(integrate::test::All,"Test full integration");

	// Tetgen interface tests
	gridTest.RunSilent(tetgen::test::CFD,"tegen API testing - CFD meshing"); 
	gridTest.RunSilent(tetgen::test::call,"tegen API testing - RSVS meshing"); 
	gridTest.RunSilent(tetgen::test::RSVSVORO,"tegen API testing - Voro to RSVS"); 
	gridTest.RunSilent(tetgen::test::RSVSVORO_Contain,"tegen API testing - Voro to RSVS");

	rsvstest::newtest();

	return(0);
}

int rsvstest::newtest(){
	customtest gridTest("3-D RSVS tests: New and breaking ");

	#ifdef TEST_KNOWN_FAILURES
	gridTest.RunSilent(Test_snakeinitflat,"Snake spawning 2D");
	gridTest.Run(Test_RSVSalgoflat,"RSVS 2D"); // Non working test - Maths not finished	
	#endif
	
	// Parameter and JSON tests
	gridTest.RunSilent(param::test::base,"parameter implementation");
	gridTest.RunSilent(param::test::symmetry,"Test json internal symmetry");
	gridTest.RunSilent(param::test::io,"parameter read write");
	gridTest.RunSilent(param::test::ioflat,"parameter read write flat format");
	gridTest.RunSilent(param::test::ipartialread,"parameter partial read"
		" write flat format");
	gridTest.RunSilent(param::test::autoflat,"Algorithm for automatic "
		"determination of flat json");
	gridTest.RunSilent(param::test::symmetry_makefillactive,
		"makefill.active symmetrical?");


	return(0);
}


int customtest::Run(function<int ()> test, const char *funcName){
	cout << "-------------------------------------------------------------"
		"---------------------------" << endl;
	cout << "-------------------------------------------------------------"
		"---------------------------" << endl;
	cout << "      Start testing " << funcName << endl;
	cout << "-------------------------------------------------------------"
		"---------------------------" << endl;
	this->prevTime=clock();
	errFlag=test();
	this->lastRunTime = clock() - this->prevTime;

	runTotal += this->lastRunTime;
	++testCount;
	cout << "-------------------------------------------------------------"
		"---------------------------" << endl;
	if (errFlag!=0){
		++errCount;
		cout << "Finished testing " << funcName << endl;
		cout << "                      - Caught Error: " << errFlag << endl;
	} else {
		cout << "Finished testing "<< funcName << endl;
		cout << "                      - No Error" << endl;
	}
	cout << "                Execution time : " << 
		double(this->lastRunTime)/double(CLOCKS_PER_SEC)*1000 << " ms" << endl; 
	cout << "-------------------------------------------------------------"
		"---------------------------" << endl;
	cout << "-------------------------------------------------------------"
		"---------------------------" << endl;
	cout << endl;
	return(errFlag);
}

/**
 * @brief      Runs a test function silently except if it returns an error.
 *
 * @param[in]  test      The test function
 * @param[in]  funcName  string descriptor for the test.
 *
 * @return     int number of errors captured.
 */
int customtest::RunSilent(function<int ()> test, const char *funcName){
	
	stringstream streamOut;
	int errFlag;

	std::cout << "START test " << this->testCount+1 
		<< ": " << funcName << endl;

	auto coutBuff=std::cout.rdbuf(streamOut.rdbuf());
	auto cerrBuff=std::cerr.rdbuf(streamOut.rdbuf());
	try{
		try{
			errFlag = this->Run(test, funcName);
		} catch(exception const& ex){
			cerr << "Exception: " << ex.what() <<endl;
			throw ex;
		}
	} catch (...){
		std::cerr << "-------------------------------------------------------------"
		"---------------------------" << endl;
		errFlag=1;
		this->unhandledError++;
		std::cerr << "Unhandled error in " << funcName << endl;
		std::cerr << "-------------------------------------------------------------"
			"---------------------------" << endl;
		std::cerr << "-------------------------------------------------------------"
			"---------------------------" << endl;
	}
	std::cout.rdbuf(coutBuff);
	std::cerr.rdbuf(cerrBuff);
	if (errFlag!=0){
		std::cerr << streamOut.str();
	} else {
		std::cout << " DONE, finished in " 
			<< double(this->lastRunTime)/double(CLOCKS_PER_SEC)*1000 
			<< " ms" << endl;
	}
	return errFlag;
}

void customtest::PrintSummary(){
	cout << endl;

	cout << "-------------------------------------------------------------"
		"---------------------------" << endl;
	cout << "-------------------------------------------------------------"
		"---------------------------" << endl;
	cout << "      Summary of Tests: (" << this->testName << ")" << endl;
	cout << "         " << testCount << " tests completed" << endl;
	cout << "         " << errCount << "  detected errors  " << endl;
	cout << "         " << unhandledError << " errors handled by the test class " << endl;
	cout << "      Total run time:" << endl;
	cout << "         " << double(this->runTotal)/double(CLOCKS_PER_SEC) 
		<< "  seconds  " << endl;
	cout << "-------------------------------------------------------------"
		"---------------------------" << endl;
	cout << "-------------------------------------------------------------"
		"---------------------------" << endl;
	#ifdef DBG_MEMLEAK
	_CrtDumpMemoryLeaks();
	#endif
}