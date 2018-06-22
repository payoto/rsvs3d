#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>

#include "snake.hpp"
#include "snakeengine.hpp"
#include "snakevel.hpp"
#include "postprocessing.hpp"
#include "meshrefinement.hpp" 

using namespace std;

// Implementation in snakstruct.cpp

// -------------------------------------------------------------------------------------------
// TEST CODE
// -------------------------------------------------------------------------------------------

int Test_SnakeStructures() { 
   // Test the functionality provided by arraystructures

	int errFlag,errTest;


	errFlag=0;

	cout << "--------------------------------------------" << endl;
	cout << "      testing coordvec" << endl;
	cout << "--------------------------------------------" << endl;
	errTest=Test_coordvec();
	errFlag= errFlag | (errTest!=0);

	cout << "--------------------------------------------" << endl;
	cout << "      testing snax" << endl;
	cout << "--------------------------------------------" << endl;
	errTest=Test_snax();
	errFlag= errFlag | (errTest!=0);

	cout << "--------------------------------------------" << endl;
	cout << "      testing snaxedge" << endl;
	cout << "--------------------------------------------" << endl;
	errTest=Test_snaxedge();
	errFlag= errFlag | (errTest!=0);

	cout << "--------------------------------------------" << endl;
	cout << "      testing Snake" << endl;
	cout << "--------------------------------------------" << endl;
	errTest=Test_snake();
	errFlag= errFlag | (errTest!=0);

	return(errFlag);
} 


int Test_coordvec(){

	coordvec testCoord,unitCoord;
	try {
		testCoord.assign(1.0,2.0,3.0);

		cout << "base vector: ";
		testCoord.disp();


		unitCoord=testCoord.Unit();
		cout << "unit vector: ";
		unitCoord.disp();

		cout << "unit access: ";
		cout << "coord vec [" << testCoord.Unit(0) << ","<< testCoord.Unit(1)<< ","<< testCoord.Unit(2) << "] norm 1" << endl;

		cout << "base oper(): ";
		cout << "coord vec [" << testCoord(0) << ","<< testCoord(1)<< ","<< testCoord(2) << "] " << endl;
		cout << "base oper(): ";testCoord.disp();
		cout << "base oper[]: ";
		cout << "coord vec [" << testCoord[0] << ","<< testCoord[1]<< ","<< testCoord[2] << "] " << endl;
		cout << "base oper[]: ";testCoord.disp();

		cout << "base ope()=: {compile error}";
		//testCoord(0)=0;
		testCoord.disp();
		cout << "base ope[]=: ";
		testCoord[0]=0;
		testCoord.disp();
		testCoord.GetNorm();
		cout << "base getnor: ";testCoord.disp();

	} catch (exception const& ex) { 
		cerr << "Exception: " << ex.what() <<endl; 
		return -1;
	} 
	return(0);

}

int Test_snax(){
	int i;

	snaxarray snaxStack,snaxStack2;  // stack of ints 
	snax singleSnax;

	try {
		singleSnax.disp();
		snaxStack.assign(5,singleSnax);
		snaxStack.disp();
		snaxStack.PopulateIndices();
		cout << endl;
		snaxStack[0].index=10;
		snaxStack.disp();

		snaxStack.PrepareForUse();
		i=snaxStack.find(10);
		cout << "found position " << i << "  Index " << snaxStack(i)->index <<  endl; 
		

		snaxStack2=snaxStack;
		snaxStack2.disp();
		snaxStack[0]=singleSnax;
		snaxStack.disp();

		cout << "Are the Same " << CompareDisp(snaxStack,snaxStack2) << endl;



	} catch (exception const& ex) { 
		cerr << "Exception: " << ex.what() <<endl; 
		return -1;
	} 
	return(0);

}

int Test_snaxedge(){
	snaxedgearray snaxStack,snaxStack2;  // stack of ints 
	snaxedge singleSnax;
	try {	
		singleSnax.normvector[1]=2;
		snaxStack.Init(4);
		snaxStack.PopulateIndices();
		snaxStack[1]=singleSnax;

		snaxStack.PrepareForUse();

		
		snaxStack.disp();

	} catch (exception const& ex) { 
		cerr << "Exception: " << ex.what() <<endl; 
		return -1;	
	} 
	return(0);

}


int Test_snake(){
	snake testSnake,testSnake2,testSnake3;
	mesh snakeMesh,snakeMesh2;

	bool errFlag;
	int errTest=0;
	//int nVert,nSurf,nEdge,nVolu;
	try {
		snakeMesh.Init(8,12,6,1);
		snakeMesh2.Init(8,12,6,1);

		testSnake.Init(&snakeMesh,8,12,6,1);
		testSnake2=testSnake;

		errFlag=CompareDisp(testSnake,testSnake2);
		cout << "Compare snakes after = assignement: 1=" << errFlag << endl ; 
		errTest=errTest+int(!errFlag);
		testSnake.ChangeIndices( 1, 2, 3,4);
		cout << "Succesfully changed indices (ChangeIndices)" << endl ; 
		testSnake.ChangeIndicesSnakeMesh( 5,6,7,8);
		cout << "Succesfully changed indices (ChangeIndicesSnakeMesh)" << endl ;

		testSnake.PrepareForUse();

		cout << "-----------------------testSnake1-----------------------" << endl;
		//testSnake.disp();
		testSnake.displight();

		testSnake3=testSnake.MakeCompatible(testSnake2);
		cout << "-----------------------testSnake2-----------------------" << endl;
		//testSnake2.disp();
		cout << "-----------------------testSnake3-----------------------" << endl;
		//testSnake3.disp();
		testSnake.Concatenate(testSnake3);
		testSnake.PrepareForUse();


		testSnake3.Init(&snakeMesh2,8,12,6,1);
		testSnake.MakeCompatible_inplace(testSnake3);
		
		try {
			testSnake.Concatenate(testSnake3);
			cerr << "Error : Concatenation between snakes on different meshes was allowed" << endl;
		} catch (exception const& ex){
			cout << "Succesfully generated failure" << endl;
		}
		cout << "-----------------------testSnake-----------------------" << endl;
		testSnake.displight();

	} catch (exception const& ex) { 
		cerr << "Exception: " << ex.what() <<endl; 
		return -1;
	} 
	return(errTest);

}

int Test_snakeinit(){
	snake testSnake;
	mesh snakeMesh;
	const char *fileToOpen;
	tecplotfile outSnake;
	double totT=0.0;
	vector<double> dt;
	vector<int> isImpact;
	int start_s,stop_s,ii;
	//bool errFlag;
	int errTest=0;
	

	try {
		fileToOpen="..\\TESTOUT\\TestSnake.plt";

		outSnake.OpenFile(fileToOpen);
		errTest+=snakeMesh.read("..\\TESTOUT\\mesh203010.dat");
		snakeMesh.PrepareForUse();
		testSnake.snakemesh=&snakeMesh;
		outSnake.PrintMesh(*(testSnake.snakemesh));

		start_s=clock();
		testSnake.PrepareForUse();
		
		SpawnAtVertex(testSnake,1022);
		SpawnAtVertex(testSnake,674);
		SpawnAtVertex(testSnake,675);
		SpawnAtVertex(testSnake,728);
		SpawnAtVertex(testSnake,729);
		SpawnAtVertex(testSnake,731);
		testSnake.displight();
		stop_s=clock();
		cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << "ms" << endl;

		start_s=clock();
		testSnake.PrepareForUse();
		for(ii=0;ii<20;++ii){
			cout << ii << " ";
			if(testSnake.snaxs.size()>0){
				//testSnake.snakeconn.TightenConnectivity();
				outSnake.PrintMesh(testSnake.snakeconn,1,totT);
			}

			Test_stepalgo(testSnake, dt, isImpact,outSnake);
			
			totT=totT+1;
		}

		if(testSnake.snaxs.size()>0){
			outSnake.PrintMesh(testSnake.snakeconn,1,totT);
		}
		stop_s=clock();
		cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << "ms" << endl;
		testSnake.displight();
	// the code you wish to time goes here
		
		//testSnake.disp();


	} catch (exception const& ex) { 
		cerr << "Exception: " << ex.what() <<endl; 
		return -1;
	} 
	return(errTest);

}
int Test_snakeinitflat(){
	snake testSnake;
	mesh snakeMesh;
	const char *fileToOpen;
	tecplotfile outSnake;
	int start_s,stop_s,ii;
	vector<double> dt;
	vector<int> isImpact;
	double totT=0.0;
	//bool errFlag;
	int errTest=0;

	try {
		fileToOpen="..\\TESTOUT\\TestSnake2.plt";

		outSnake.OpenFile(fileToOpen);
		errTest+=snakeMesh.read("..\\TESTOUT\\mesh230.dat");
		snakeMesh.PrepareForUse();
		snakeMesh.SetBorders();
		snakeMesh.PrepareForUse();

		testSnake.snakemesh=&snakeMesh;
		outSnake.PrintMesh(*(testSnake.snakemesh));

		start_s=clock();
		testSnake.PrepareForUse();
		SpawnAtVertex(testSnake,13);
		SpawnAtVertex(testSnake,16);
		SpawnAtVertex(testSnake,32);
		testSnake.displight();


		testSnake.PrepareForUse();
		for(ii=0;ii<200;++ii){
			cout << ii << " ";
			if(testSnake.snaxs.size()>0){
				outSnake.PrintMesh(testSnake.snakeconn,1,totT);
			}

			Test_stepalgo(testSnake, dt, isImpact,outSnake);
			
			totT=totT+1;
		}
		//testSnake.disp();
	// the code you wish to time goes here
		stop_s=clock();
		cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << "ms" << endl;

		
		//testSnake.disp();


	} catch (exception const& ex) { 
		cerr << "Exception: " << ex.what() <<endl; 
		return -1;
	} 
	return(errTest);

}

void Test_stepalgo(snake &testSnake, vector<double> dt, vector<int> isImpact, tecplotfile &outSnake){

	int start_s,stop_s,start_f;

	start_s=clock();
	start_f=clock();

	CleanupSnakeConnec(testSnake);

	#ifdef SAFE_ALGO
	if (testSnake.Check3D()){
		testSnake.snakeconn.TestConnectivityBiDir();
	}
	#endif
	stop_s=clock();
	cout << "cleanup: " << double(stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << "ms  " ;
	start_s=clock();
	CalculateSnakeVel(testSnake);
	testSnake.CalculateTimeStep(dt,0.25);
	testSnake.UpdateDistance(dt);
	testSnake.UpdateCoord();

	stop_s=clock();
	cout << "dt dist coord: " << double(stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << "ms  " ;
	start_s=clock();

	testSnake.PrepareForUse();

	stop_s=clock();
	cout << "PrepareForUse: " << double(stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << "ms  " ;
	start_s=clock();

	
	testSnake.SnaxImpactDetection(isImpact);
	MergeAllContactVertices(testSnake, isImpact);

	stop_s=clock();
	cout << "Merge: " << double(stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << "ms  " ;
	start_s=clock();
	testSnake.PrepareForUse();
	CleanupSnakeConnec(testSnake);
	#ifdef SAFE_ALGO
	if (testSnake.Check3D()){
		testSnake.snakeconn.TestConnectivityBiDir();
	}
	#endif
	SpawnArrivedSnaxels(testSnake,isImpact);


	stop_s=clock();
	cout << "Spawn: " << double(stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << "ms  " ;
	start_s=clock();

	testSnake.SnaxImpactDetection(isImpact);
	testSnake.PrepareForUse();
	#ifdef SAFE_ALGO
	if (testSnake.Check3D()){
		testSnake.snakeconn.TestConnectivityBiDir();
	}
	#endif
	MergeAllContactVertices(testSnake, isImpact);
	testSnake.PrepareForUse();
	#ifdef SAFE_ALGO
	if (testSnake.Check3D()){
		testSnake.snakeconn.TestConnectivityBiDir();
	}
	#endif
	//testSnake.ForceCloseContainers();
	//testSnake.snakeconn.TightenConnectivity();
	//testSnake.PrepareForUse();
	stop_s=clock();
	cout << "Impact merge: " << double(stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << "ms  " ;
	cout << "Total: " << double(stop_s-start_f)/double(CLOCKS_PER_SEC)*1000 << "ms  " << endl;


}

int Test_snakeOrderEdges(){

	mesh snakeMesh;
	const char *fileToOpen;
	tecplotfile outSnake;
	
	//bool errFlag;
	int errTest=0;

	try {
		fileToOpen="..\\TESTOUT\\TestOrderEdges.plt";

		outSnake.OpenFile(fileToOpen);
		errTest+=snakeMesh.read("..\\TESTOUT\\mesh234.dat");
		snakeMesh.HashArray();
		snakeMesh.SetMaxIndex();
		
		outSnake.PrintMesh(snakeMesh);

		snakeMesh.OrderEdges();
		outSnake.PrintMesh(snakeMesh);
		//testSnake.disp();


	} catch (exception const& ex) { 
		cerr << "Exception: " << ex.what() <<endl; 
		return -1;
	} 
	return(errTest);

}

int Test_MeshRefinement(){

	mesh snakeMesh, voluMesh;
	const char *fileToOpen;
	tecplotfile outSnake;
	vector<int> elmMapping,dims;
	int ii;
	
	//bool errFlag;
	int errTest=0;

	try {
		fileToOpen="..\\TESTOUT\\Test_Multimesh.plt";

		outSnake.OpenFile(fileToOpen);
		errTest+=snakeMesh.read("..\\TESTOUT\\mesh203010.dat");
		snakeMesh.PrepareForUse();
		outSnake.PrintMesh(snakeMesh);

		snakeMesh.OrderEdges();
		outSnake.PrintMesh(snakeMesh);
		//testSnake.disp();
		for (ii=0;ii<snakeMesh.volus.size();++ii){
			elmMapping.push_back(1);
		}
		dims.assign(3,0);
		dims[0]=2;dims[1]=3;dims[2]=1;
		CartesianMapping(snakeMesh,  elmMapping, dims);
		DisplayVector(elmMapping);
		CoarsenMesh(snakeMesh,voluMesh,elmMapping);
		for (ii=0;ii<voluMesh.volus.size();++ii){
			voluMesh.volus[ii].fill=(double(rand()%1001)/1000.0+0.5);
		}
		voluMesh.PrepareForUse();
		outSnake.PrintMesh(voluMesh);
	} catch (exception const& ex) { 
		cerr << "Exception: " << ex.what() <<endl; 
		return -1;
	} 
	return(errTest);

}
