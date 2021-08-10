#include <cstdlib>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>

#include "test.hpp"

#include "RSVSalgorithm.hpp"
#include "RSVSintegration.hpp"
#include "RSVSmath.hpp"
#include "arraystructures.hpp"
#include "matrixtools.hpp"
#include "mesh.hpp"
#include "meshrefinement.hpp"
#include "parameters.hpp"
#include "polyscopersvs.hpp"
#include "postprocessing.hpp"
#include "snake.hpp"
#include "snakeengine.hpp"
#include "tetgenrsvs.hpp"
#include "triangulate.hpp"
#include "voxel.hpp"
using namespace rsvstest;
using namespace std;

int rsvstest::maintest()
{

    customtest gridTest("all tests (expected run time: 15 minutes)");
    gridTest.Run(rsvstest::shorttest, "Short tests");
    gridTest.Run(rsvstest::longevolution, "Long tests");
    gridTest.Run(rsvstest::integrationprocesses, "Various integration and derivative tests");
    gridTest.Run(rsvstest::newtest, "New tests");

    return (gridTest.ReturnErrCount());
}

int rsvstest::shorttest()
{

    customtest gridTest("short tests (expected run time: 7 minutes)");

    gridTest.RunSilent(rsvstest::arraystructtemplates, "Test Arraystructures templates", 1);
    gridTest.RunSilent(rsvstest::meshprocesses, "Test meshes", 15);
    gridTest.RunSilent(rsvstest::snakeprocesses, "Test snake processes", 41);
    gridTest.Run(rsvstest::RSVSprocesses, "Test RSVS process", 150);
    gridTest.Run(rsvstest::tetgenprocesses, "Tetgen interface tests", 200);
    gridTest.RunSilent(rsvstest::JSONprocesses, "Parameter and JSON tests", 1);

    return (gridTest.ReturnErrCount());
}

int rsvstest::failingtest()
{
    customtest gridTest("failing tests (expected run time: 1-mn to 7h)");
    gridTest.RunSilent(rsvstest::RSVS2Dprocesses, "Test 2D-RSVS process", 6);
    gridTest.Run(tetgen::test::RSVSVORO, "tegen API testing - Voro to RSVS");
    gridTest.Run(tetgen::test::RSVSVORO_Contain, "tegen API testing - Voro to RSVS");

    return (gridTest.ReturnErrCount());
}

int rsvstest::arraystructtemplates()
{
    customtest gridTest("Test Arraystructures templates");

    // BASE templatess
    gridTest.RunSilent(Test_ArrayStructures, "arraystructures");
    gridTest.RunSilent(TestTemplate_ArrayStruct<vert>, "TestTemplate_ArrayStruct<vert>");
    gridTest.RunSilent(TestTemplate_ArrayStruct<edge>, "TestTemplate_ArrayStruct<edge>");
    gridTest.RunSilent(TestTemplate_ArrayStruct<surf>, "TestTemplate_ArrayStruct<surf>");
    gridTest.RunSilent(TestTemplate_ArrayStruct<volu>, "TestTemplate_ArrayStruct<volu>");
    gridTest.RunSilent(TestTemplate_ArrayStruct<snax>, "TestTemplate_ArrayStruct<snax>");
    gridTest.RunSilent(TestTemplate_ArrayStruct<snaxedge>, "TestTemplate_ArrayStruct<snaxedge>");
    gridTest.RunSilent(TestTemplate_ArrayStruct<snaxsurf>, "TestTemplate_ArrayStruct<snaxsurf>");
    return (gridTest.ReturnErrCount());
}

int rsvstest::meshprocesses()
{
    customtest gridTest("Test meshes");

    // Building blocks
    gridTest.RunSilent(Test_BuildBlockGrid_noout, "Voxel");
    gridTest.RunSilent(Test_tecplotfile, "post-processing class");
    gridTest.RunSilent(Test_MeshOut, "Mesh output");
    gridTest.RunSilent(Test_MeshOut2, "Functionalised Mesh output");
    gridTest.RunSilent(Test_MeshReOrder, "Mesh Reordering");
    gridTest.RunSilent(Test_surfcentre, "test SurfCentroid");
    gridTest.RunSilent(Test_MeshRefinement, "Multi-Level Meshes");
    gridTest.RunSilent(Test_MeshOrient, "Output mesh orientation");
    gridTest.RunSilent(Test_Crop, "test cropping of meshes");

    return (gridTest.ReturnErrCount());
}

int rsvstest::snakeprocesses()
{
    customtest gridTest("Test snake processes");

    // Snakstruct 3D tests
    gridTest.RunSilent(Test_snakeOrderEdges, "Snake Order error");
    gridTest.RunSilent(Test_SnakeStructures, "Snake containers");
    gridTest.RunSilent(Test_RSVSalgo_init, "RSVS spawn");
    gridTest.RunSilent(Test_RSVSvoro_init, "Snake spawning voronoi");
    gridTest.RunSilent(Test_snakeinit_random_short, "(short) Snake rand velocity");
    gridTest.RunSilent(Test_snakeinit_unit_short, "(short) Snake unit velocity (reflected)");
    gridTest.RunSilent(Test_snakeinit_unitnoreflect_short, "(short) Snake unit velocity "
                                                           "(no reflection)");
    return (gridTest.ReturnErrCount());
}

int rsvstest::longevolution()
{
    customtest gridTest("Test process with long processing");
    gridTest.RunSilent(Test_snakeinit_random, "Snake rand velocity");
    gridTest.RunSilent(Test_snakeinit_unit, "Snake unit velocity (reflected)");
    gridTest.RunSilent(Test_snakeinit_unitnoreflect, "Snake unit velocity "
                                                     "(no reflection)");
    return (gridTest.ReturnErrCount());
}
int rsvstest::RSVSprocesses()
{
    customtest gridTest("Test RSVS process");

    gridTest.RunSilent(Test_RSVSalgo, "Snake RSVS from spawn", 90);
    gridTest.RunSilent(Test_snakeRSVS, "Snake RSVS", 5);
    gridTest.RunSilent(Test_snakeRSVS_singlevol, "Snake RSVS single vol", 5);
    gridTest.RunSilent(Test_RSVSalgo_singlevol_fullmath, "Snake RSVS algorithm (full maths)", 25);
    gridTest.RunSilent(Test_RSVSalgo_singlevol_sparse, "Snake RSVS algorithm (sparse)", 12);
    return (gridTest.ReturnErrCount());
}

int rsvstest::RSVS2Dprocesses()
{
    customtest gridTest("Test 2D-RSVS process");

    gridTest.RunSilent(Test_snakeinitflat, "Snake spawning 2D");
    gridTest.Run(Test_RSVSalgoflat, "RSVS 2D"); // Non working test - Maths not finished

    return (gridTest.ReturnErrCount());
}

int rsvstest::tetgenprocesses()
{
    customtest gridTest("Tetgen interface tests");

    // Tetgen interface tests
    gridTest.RunSilent(tetgen::test::CFD, "tegen API testing - CFD meshing", 26);
    gridTest.RunSilent(tetgen::test::call, "tegen API testing - RSVS meshing", 30);
    // gridTest.RunSilent(tetgen::test::RSVSVORO,"tegen API testing - Voro to RSVS", 78);
    // gridTest.RunSilent(tetgen::test::RSVSVORO_Contain,"tegen API testing - Voro to RSVS", 42);
    return (gridTest.ReturnErrCount());
}

int rsvstest::polyscopeprocesses()
{
    customtest gridTest("Test polyscope processes");
    gridTest.RunSilent(polyscopersvs::test::init, "Initialise polyscope");
    gridTest.RunSilent(polyscopersvs::test::init, "Try to initialise polyscope twice");
    gridTest.RunSilent(polyscopersvs::test::show, "Show the polyscope window");
    gridTest.RunSilent(polyscopersvs::test::meshShow, "Print a mesh");
    return (gridTest.ReturnErrCount());
}

int rsvstest::JSONprocesses()
{
    customtest gridTest("Parameter and JSON tests");

    gridTest.RunSilent(param::test::base, "parameter implementation");
    gridTest.RunSilent(param::test::symmetry, "Test json internal symmetry");
    gridTest.RunSilent(param::test::io, "parameter read write");
    gridTest.RunSilent(param::test::ioflat, "parameter read write flat format");
    gridTest.RunSilent(param::test::ipartialread, "parameter partial read"
                                                  " write flat format");
    gridTest.RunSilent(param::test::autoflat, "Algorithm for automatic "
                                              "determination of flat json");
    gridTest.RunSilent(tetgen::test::api, "tegen test api parameter class");
    return (gridTest.ReturnErrCount());
}

int rsvstest::integrationprocesses()
{

    customtest gridTest("Various integration and derivative tests");

    // RSVS and integration tests
    gridTest.RunSilent(integrate::test::Prepare, "Mesh integration function");
    gridTest.RunSilent(integrate::test::All, "Test full integration");
    gridTest.Run(Test_SurfCentreDerivatives, "Surf centroid");
    gridTest.Run(Test_Matrix3D, "'3D' matrix maths");
    gridTest.Run(integrate::test::CompareSurfCentreDerivatives, "Test Surf centroid");

    gridTest.Run(integrate::test::CompareDerivativesSpike, "Test spike issues");
    gridTest.Run(integrate::test::CompareDerivativesSpikeNoDPos, "Test spike issues no dPos");
    gridTest.Run(integrate::test::StudyDerivatives, "Generate Derivative study");

    return (gridTest.ReturnErrCount());
}

int rsvstest::newtest()
{
    customtest gridTest("3-D RSVS tests: New and breaking ");

#ifdef TEST_KNOWN_FAILURES
    gridTest.RunSilent(Test_snakeinitflat, "Snake spawning 2D");
    gridTest.Run(Test_RSVSalgoflat, "RSVS 2D"); // Non working test - Maths not finished
    gridTest.Run(integrate::test::CompareDerivativesSpike, "Test spike issues");
    gridTest.Run(integrate::test::CompareDerivativesSpikeNoDPos, "Test spike issues no dPos");
    gridTest.Run(integrate::test::StudyDerivatives, "Generate Derivative study");
#endif
#ifdef IGNORE_THESE_TESTS
    // Parameter and JSON tests
    gridTest.Run(tetgen::test::CFD, "tegen API testing - CFD meshing");
    gridTest.Run(tetgen::test::call, "tegen API testing - RSVS meshing");
    gridTest.Run(tetgen::test::RSVSVORO, "tegen API testing - Voro to RSVS");
    gridTest.Run(tetgen::test::RSVSVORO_Contain, "tegen API testing - Voro to RSVS");
#endif
    gridTest.Run(rsvstest::polyscopeprocesses, "Test polyscope processes");
    return (gridTest.ReturnErrCount());
}

int customtest::Run(function<int()> test, const char *funcName, int expectedTime)
{
    cout << "-------------------------------------------------------------"
            "---------------------------"
         << endl;
    cout << "-------------------------------------------------------------"
            "---------------------------"
         << endl;
    cout << "      Start testing " << funcName;
    if (expectedTime >= 0)
    {
        std::cout << "  (expected execution time: " << expectedTime << " s)";
    }
    std::cout << std::endl;
    cout << "-------------------------------------------------------------"
            "---------------------------"
         << endl;
    this->prevTime = clock();
    errFlag = test();
    this->lastRunTime = clock() - this->prevTime;

    runTotal += this->lastRunTime;
    ++testCount;
    cout << "-------------------------------------------------------------"
            "---------------------------"
         << endl;
    if (errFlag != 0)
    {
        ++errCount;
        cout << "Finished testing " << funcName << endl;
        cout << "                      - Caught Error: " << errFlag << endl;
    }
    else
    {
        cout << "Finished testing " << funcName << endl;
        cout << "                      - No Error" << endl;
    }
    cout << "                Execution time : " << double(this->lastRunTime) / double(CLOCKS_PER_SEC) * 1000 << " ms"
         << endl;
    cout << "-------------------------------------------------------------"
            "---------------------------"
         << endl;
    cout << "-------------------------------------------------------------"
            "---------------------------"
         << endl;
    cout << endl;
    return (errFlag);
}

/**
 * @brief      Runs a test function silently except if it returns an error.
 *
 * @param[in]  test      The test function
 * @param[in]  funcName  string descriptor for the test.
 *
 * @return     int number of errors captured.
 */
int customtest::RunSilent(function<int()> test, const char *funcName, int expectedTime)
{

    stringstream streamOut;
    int errFlag;

    std::cout << "START test " << this->testCount + 1 << ": " << funcName;
    if (expectedTime >= 0)
    {
        std::cout << "  (expected execution time: " << expectedTime << " s)";
    }
    std::cout << std::endl;

    auto coutBuff = std::cout.rdbuf(streamOut.rdbuf());
    auto cerrBuff = std::cerr.rdbuf(streamOut.rdbuf());
    try
    {
        try
        {
            errFlag = this->Run(test, funcName, expectedTime);
        }
        catch (exception const &ex)
        {
            cerr << "Exception: " << ex.what() << endl;
            throw ex;
        }
    }
    catch (...)
    {
        std::cerr << "-------------------------------------------------------------"
                     "---------------------------"
                  << endl;
        errFlag = 1;
        this->unhandledError++;
        std::cerr << "Unhandled error in " << funcName << endl;
        std::cerr << "-------------------------------------------------------------"
                     "---------------------------"
                  << endl;
        std::cerr << "-------------------------------------------------------------"
                     "---------------------------"
                  << endl;
    }
    std::cout.rdbuf(coutBuff);
    std::cerr.rdbuf(cerrBuff);
    if (errFlag != 0)
    {
        std::cerr << streamOut.str();
    }
    else
    {
        std::cout << " DONE, finished in " << double(this->lastRunTime) / double(CLOCKS_PER_SEC) * 1000 << " ms"
                  << endl;
    }
    return errFlag;
}

void customtest::PrintSummary()
{
    cout << endl;

    cout << "-------------------------------------------------------------"
            "---------------------------"
         << endl;
    cout << "-------------------------------------------------------------"
            "---------------------------"
         << endl;
    cout << "      Summary of Tests: (" << this->testName << ")" << endl;
    cout << "         " << testCount << " tests completed" << endl;
    cout << "         " << errCount << "  detected errors  " << endl;
    cout << "         " << unhandledError << " errors handled by the test class " << endl;
    cout << "      Total run time:" << endl;
    cout << "         " << double(this->runTotal) / double(CLOCKS_PER_SEC) << "  seconds  " << endl;
    cout << "-------------------------------------------------------------"
            "---------------------------"
         << endl;
    cout << "-------------------------------------------------------------"
            "---------------------------"
         << endl;
#ifdef DBG_MEMLEAK
    _CrtDumpMemoryLeaks();
#endif
}
