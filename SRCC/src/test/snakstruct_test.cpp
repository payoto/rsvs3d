#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>

#include "snake.hpp"
#include "snakeengine.hpp"
#include "snakevel.hpp"
#include "postprocessing.hpp"
#include "meshrefinement.hpp" 
#include "RSVSmath.hpp"
#include "RSVSinterface.hpp"
#include "RSVSalgorithm.hpp"
#include "RSVSintegration.hpp"

using namespace std;

// Implementation in snakstruct.cpp


void PrintRSVSSnake(tecplotfile &outSnake, snake &testSnake, double totT, triangulation &testTriangle,
	mesh &triMesh, triangulation &triRSVS, mesh &voluMesh, int nVoluZone, int ii){
	vector<int> vertList;
	int jj;
	if(testSnake.snaxs.size()>0){
		//testSnake.snakeconn.TightenConnectivity();
		outSnake.PrintMesh(testSnake.snakeconn,1,totT);

		testSnake.snakeconn.PrepareForUse();
		testTriangle.stattri.clear();
		testTriangle.trivert.clear();
		testTriangle.PrepareForUse();
		TriangulateMesh(testSnake.snakeconn,testTriangle);
		MeshTriangulation(triMesh,testSnake.snakeconn,testTriangle.stattri, testTriangle.trivert);
		outSnake.PrintMesh(triMesh,2,totT);
		MeshTriangulation(triMesh,testSnake.snakeconn,triRSVS.dynatri, triRSVS.trivert);
		outSnake.PrintMesh(triMesh,3,totT);
		outSnake.PrintTriangulation(triRSVS,&triangulation::dynatri,4,totT);
		if (ii==0){
			outSnake.PrintTriangulation(triRSVS,&triangulation::dynatri,5,totT,3);
			outSnake.PrintTriangulation(triRSVS,&triangulation::dynatri,6,totT,3);
			outSnake.PrintTriangulation(triRSVS,&triangulation::dynatri,7,totT,3);
		}
		outSnake.PrintTriangulation(triRSVS,&triangulation::intertri,5,totT,3);
		outSnake.PrintTriangulation(triRSVS,&triangulation::trisurf,6,totT,3);
		if (int(triRSVS.acttri.size())>0){
			outSnake.PrintTriangulation(triRSVS,&triangulation::stattri,7,totT,3,triRSVS.acttri);
		}
		
		vertList.clear();
		for(jj=0;jj<int(testSnake.isMeshVertIn.size()); ++jj){
			if(testSnake.isMeshVertIn[jj]){
				vertList.push_back(testSnake.snakemesh->verts(jj)->index);
			}
		}
		if(int(testSnake.isMeshVertIn.size())==0){
			vertList.push_back(testSnake.snakemesh->verts(0)->index);
		}
		outSnake.PrintMesh(*(testSnake.snakemesh),8,totT,4,vertList);
		outSnake.PrintVolumeDat(voluMesh,nVoluZone,9,totT);
		outSnake.PrintSnake(testSnake, 10, totT);

	}
}


void PrepareMultiLvlSnake(mesh &snakeMesh, mesh &voluMesh, snake &testSnake,
	vector<int> &dims, triangulation &triRSVS){
	vector<int> elmMapping;
	SQPcalc calcVolus;
	int ii;

	snakeMesh.PrepareForUse();
	snakeMesh.OrientFaces();
	///// Generate Coarser Volume Mesh
	testSnake.snakemesh=&snakeMesh;
	//testSnake.disp();
	for (ii=0;ii<snakeMesh.volus.size();++ii){
		elmMapping.push_back(1);
	}
	CartesianMapping(snakeMesh,  elmMapping, dims);
	CoarsenMesh(snakeMesh,voluMesh,elmMapping);
	snakeMesh.AddParent(&voluMesh,elmMapping);
	
	sort(elmMapping);
	unique(elmMapping);
	DisplayVector(elmMapping);
	for (ii=0;ii<voluMesh.volus.size();++ii){
		voluMesh.volus[ii].target=(double(rand()%1001)/1000.0);
	}
	voluMesh.PrepareForUse();
	voluMesh.OrientFaces();

	calcVolus.CalculateMesh(voluMesh);
	calcVolus.ReturnConstrToMesh(voluMesh,&volu::volume);

	triRSVS.stattri.clear();
	triRSVS.trivert.clear();
	triRSVS.PrepareForUse();
	TriangulateMesh(snakeMesh,triRSVS);

	testSnake.PrepareForUse();

}


void Test_randvelstep(snake &testSnake, vector<double> dt, vector<int> &isImpact){
	CalculateSnakeVelRand(testSnake);
	testSnake.CalculateTimeStep(dt,0.25);
	testSnake.UpdateDistance(dt);
	testSnake.UpdateCoord();
	testSnake.PrepareForUse();
	Test_stepalgo(testSnake, isImpact);
}

void Test_randvelstep_mc(snake &testSnake, vector<double> dt, vector<int> &isImpact){
	CalculateSnakeVelRand(testSnake);
	testSnake.CalculateTimeStep(dt,0.25);
	testSnake.UpdateDistance(dt);
	testSnake.UpdateCoord();
	testSnake.PrepareForUse();
	Test_stepalgo_mergeclean(testSnake, isImpact);
}


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
		cout << "coord vec [" << testCoord.Unit(0) << ","<< testCoord.Unit(1)<< ","<< 
		testCoord.Unit(2) << "] norm 1" << endl;

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
	mesh snakeMesh, triMesh;
	triangulation testTriangle;
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
		
		snakeMesh.OrientFaces();
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
		for(ii=0;ii<300;++ii){
			cout << ii << " ";
			
			if(testSnake.snaxs.size()>0){
				//testSnake.snakeconn.TightenConnectivity();
				outSnake.PrintMesh(testSnake.snakeconn,1,totT);

				testSnake.snakeconn.PrepareForUse();
				testTriangle.stattri.clear();
				testTriangle.trivert.clear();
				testTriangle.PrepareForUse();
				TriangulateMesh(testSnake.snakeconn,testTriangle);
				testTriangle.CalcTriVertPos();
				MeshTriangulation(triMesh,testSnake.snakeconn,testTriangle.stattri, testTriangle.trivert);
				outSnake.PrintMesh(triMesh,2,totT);
			}
			Test_randvelstep(testSnake, dt, isImpact);
			cout << endl;
			
			totT=totT+1;
		}

		if(testSnake.snaxs.size()>0){
			outSnake.PrintMesh(testSnake.snakeconn,1,totT);
			testTriangle.stattri.clear();
			testTriangle.trivert.clear();
			testTriangle.PrepareForUse();
			TriangulateMesh(testSnake.snakeconn,testTriangle);
			testTriangle.CalcTriVertPos();

			MeshTriangulation(triMesh,testSnake.snakeconn,testTriangle.stattri, testTriangle.trivert);
			outSnake.PrintMesh(triMesh,2,totT);
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


int Test_snakeinit_MC(){ 
	snake testSnake;
	mesh snakeMesh, triMesh;
	triangulation testTriangle;
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
		
		snakeMesh.OrientFaces();
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
		for(ii=0;ii<100;++ii){
			cout << ii << " ";
			
			if(testSnake.snaxs.size()>0){
				//testSnake.snakeconn.TightenConnectivity();
				outSnake.PrintMesh(testSnake.snakeconn,1,totT);

				testSnake.snakeconn.PrepareForUse();
				testTriangle.stattri.clear();
				testTriangle.trivert.clear();
				testTriangle.PrepareForUse();
				TriangulateMesh(testSnake.snakeconn,testTriangle);
				testTriangle.CalcTriVertPos();
				MeshTriangulation(triMesh,testSnake.snakeconn,testTriangle.stattri, testTriangle.trivert);
				outSnake.PrintMesh(triMesh,2,totT);
			}
			if(ii==30){
				cout << "break here" << endl;
			}
			Test_randvelstep_mc(testSnake, dt, isImpact);
			cout << endl;
			
			totT=totT+1;
		}

		if(testSnake.snaxs.size()>0){
			outSnake.PrintMesh(testSnake.snakeconn,1,totT);
			testTriangle.stattri.clear();
			testTriangle.trivert.clear();
			testTriangle.PrepareForUse();
			TriangulateMesh(testSnake.snakeconn,testTriangle);
			testTriangle.CalcTriVertPos();

			MeshTriangulation(triMesh,testSnake.snakeconn,testTriangle.stattri, testTriangle.trivert);
			outSnake.PrintMesh(triMesh,2,totT);
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

			Test_stepalgo(testSnake, isImpact);
			cout << endl;
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

void Test_stepalgo(snake &testSnake,  vector<int> &isImpact){

	SnakeConnectivityUpdate(testSnake,  isImpact);

}

void Test_stepalgo_mergeclean(snake &testSnake, vector<int> &isImpact){

	int start_s;

	start_s=clock();

	

	start_s=TimeStamp("position: ", start_s);

	testSnake.PrepareForUse();
	MergeCleanSnake(testSnake, isImpact);
	testSnake.PrepareForUse();
	start_s=TimeStamp("MergeClean: ", start_s);

	testSnake.SnaxImpactDetection(isImpact);
	SpawnArrivedSnaxels(testSnake,isImpact);
	start_s=TimeStamp("Spawn: ", start_s);

	testSnake.PrepareForUse();
	MergeCleanSnake(testSnake, isImpact);
	testSnake.PrepareForUse();

	testSnake.OrientFaces();
	start_s=TimeStamp("Clean: ", start_s);

	

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

	mesh snakeMesh, voluMesh, triMesh;
	const char *fileToOpen;
	tecplotfile outSnake;
	vector<int> elmMapping,dims;
	triangulation testTriangle;
	SQPcalc calcObj,calcObj2;
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
		snakeMesh.OrientFaces();
		outSnake.PrintMesh(snakeMesh);
		//testSnake.disp();
		for (ii=0;ii<snakeMesh.volus.size();++ii){
			elmMapping.push_back(1);
		}
		dims.assign(3,0);
		dims[0]=2;dims[1]=3;dims[2]=1;
		CartesianMapping(snakeMesh,  elmMapping, dims);
		CoarsenMesh(snakeMesh,voluMesh,elmMapping);
		snakeMesh.AddParent(&voluMesh,elmMapping);

		sort(elmMapping);
		unique(elmMapping);
		DisplayVector(elmMapping);
		for (ii=0;ii<voluMesh.volus.size();++ii){
			voluMesh.volus[ii].fill=(double(rand()%1001)/1000.0);
			voluMesh.volus[ii].target=voluMesh.volus[ii].fill;
			voluMesh.volus[ii].error=voluMesh.volus[ii].fill;
		}
		voluMesh.PrepareForUse();
		voluMesh.TightenConnectivity();
		voluMesh.OrderEdges();
		snakeMesh.OrientFaces();
		voluMesh.OrientFaces();
		triMesh.OrientFaces();

		testTriangle.PrepareForUse();
		TriangulateMesh(snakeMesh,testTriangle);
		testTriangle.PrepareForUse();
		testTriangle.CalcTriVertPos();

		testTriangle.acttri.clear();
		testTriangle.acttri.reserve(testTriangle.stattri.size());
		for (ii=0; ii< testTriangle.stattri.size(); ++ii){
			testTriangle.acttri.push_back(testTriangle.stattri(ii)->index);
		}

		testTriangle.PrepareForUse();
		calcObj.CalculateTriangulation(testTriangle);
		calcObj.ReturnConstrToMesh(testTriangle);
		calcObj.Print2Screen();

		calcObj2.CalculateMesh(voluMesh);
		calcObj2.ReturnConstrToMesh(voluMesh,&volu::error);
		calcObj2.Print2Screen();

		MeshTriangulation(triMesh,voluMesh,testTriangle.stattri, testTriangle.trivert);

		outSnake.PrintMesh(voluMesh);
		outSnake.PrintMesh(triMesh);


	} catch (exception const& ex) { 
		cerr << "Exception: " << ex.what() <<endl; 
		return -1;
	} 
	return(errTest);

}

int Test_surfcentre(){ 
	// int ii,n;
	// vector<int> vertind;
	vector<vector<double> const *> veccoord;
	SurfCentroid tempCalc;
	vector<double> v1 = {0.0 , 0.0, 0.0 };
	vector<double> v2 = {1.0 , 1.0, 0.0 };
	vector<double> v3 = {0.0 , 1.0, 1.0 };
	vector<double> v4 = {1.0 , 0.0, 0.0 };
	vector<double> v5 = {1.0 , 0.0, 0.0 };
	// ArrayVec<double> tempCoord,jac,hes;

	// coord.assign(0,0,0);
	// n=int(surfin.edgeind.size());

	// veccoord.reserve(n);
	// ConnVertFromConnEdge(meshin, surfin.edgeind,vertind);

	// for(ii=0; ii<n; ++ii){
	// 	veccoord.push_back(&(meshin.verts.isearch(vertind[ii])->coord));
	// }
	veccoord.push_back(&v1);
	veccoord.push_back(&v2);
	veccoord.push_back(&v3);
	veccoord.push_back(&v4);
	veccoord.push_back(&v5);
	veccoord.push_back(&v5);
	veccoord.push_back(&v5);
	veccoord.push_back(&v5);
	tempCalc.assign(veccoord);
	tempCalc.Calc();

	// tempCalc.ReturnDat(tempCoord,jac,hes);
	// coord.assign(tempCoord[0][0],tempCoord[1][0],tempCoord[2][0]);

	return(0);
}

int Test_RSVSalgo_init(){
	// int nVoluZone, ii;

	snake testSnake, testSnake2;
	mesh snakeMesh,  voluMesh, voluMesh2;
	// mesh triMesh;
	triangulation testTriangle,triRSVS, triRSVS2;
	vector<int> dims;
	const char *fileToOpen;
	tecplotfile outSnake;
	// double totT=0.0;
	// vector<double> dt;
	// vector<int> isImpact;
	int start_s,stop_s;
	//bool errFlag;
	int errTest=0;
	

	dims.assign(3,0);
	dims[0]=2;dims[1]=3;dims[2]=1;
	try {
		fileToOpen="..\\TESTOUT\\TestAlgoRSVS.plt";

		outSnake.OpenFile(fileToOpen);
		errTest+=snakeMesh.read("..\\TESTOUT\\mesh203010.dat");
		
		PrepareMultiLvlSnake(snakeMesh,voluMesh,testSnake,dims,triRSVS);
		PrepareMultiLvlSnake(snakeMesh,voluMesh2,testSnake2,dims,triRSVS2);
		voluMesh.volus[0].target=0.0;
		voluMesh.volus[3].target=0.0;
		voluMesh.volus[4].target=1.0;
		voluMesh.PrepareForUse();
		outSnake.PrintMesh(*(testSnake.snakemesh));
		outSnake.PrintMesh(voluMesh);
		// nVoluZone=outSnake.ZoneNum();
		
		start_s=clock();
		SpawnRSVS(testSnake2,0);
		SpawnRSVS(testSnake,1);
		testSnake.PrepareForUse();
		testSnake2.PrepareForUse();

		stop_s=clock();
		cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << "ms" << endl;
		testSnake.displight();
		outSnake.PrintMesh(testSnake.snakeconn);
		outSnake.PrintSnakeInternalPts(testSnake);
		outSnake.PrintMesh(testSnake2.snakeconn);
		outSnake.PrintSnakeInternalPts(testSnake2);

	} catch (exception const& ex) { 
		cerr << "Exception: " << ex.what() <<endl; 
		return -1;
	} 
	return(errTest);
}

int Test_RSVSalgo(){
	// int nVoluZone, ii;

	snake testSnake;
	mesh snakeMesh,  voluMesh, triMesh;
	// mesh triMesh;
	triangulation testTriangle,triRSVS;
	vector<int> dims;
	const char *fileToOpen;
	tecplotfile outSnake, outSnake2;
	// double totT=0.0;
	// vector<double> dt;
	// vector<int> isImpact;
	int start_s,stop_s;
	//bool errFlag;
	int errTest=0;
	SQPcalc calcObj;
	

	dims.assign(3,0);
	dims[0]=2;dims[1]=3;dims[2]=1;
	try {
		fileToOpen="..\\TESTOUT\\TestAlgoRSVSstep.plt";
		outSnake.OpenFile(fileToOpen);
		fileToOpen="..\\TESTOUT\\TestAlgoRSVSstep_snake.plt";
		outSnake2.OpenFile(fileToOpen);
		errTest+=snakeMesh.read("..\\TESTOUT\\mesh203010.dat");
		
		PrepareMultiLvlSnake(snakeMesh,voluMesh,testSnake,dims,triRSVS);
		calcObj.BuildMathArrays(1,1);
		voluMesh.volus[0].target=0.0;
		voluMesh.volus[3].target=0.0;
		voluMesh.volus[4].target=1.0;
		voluMesh.PrepareForUse();
		outSnake.PrintMesh(*(testSnake.snakemesh));
		outSnake.PrintMesh(voluMesh);
		int nVoluZone;
		nVoluZone=outSnake.ZoneNum();
		
		start_s=clock();

		SpawnRSVS(testSnake,1);
		testSnake.PrepareForUse();

		triRSVS.PrepareForUse();

		TriangulateSnake(testSnake,triRSVS);
		triRSVS.PrepareForUse();
		triRSVS.CalcTriVertPos();
		int ii;
		double totT=0.0;
		vector<double> dt;
		vector<int> isImpact;
		MaintainTriangulateSnake(triRSVS);

		for(ii=0;ii<5;++ii){
			cout << ii << " ";
			PrintRSVSSnake(outSnake, testSnake, totT, testTriangle,
				triMesh, triRSVS, voluMesh, nVoluZone, ii);
			stop_s=clock();
			Test_stepalgoRSVS(testSnake,triRSVS, dt, isImpact, calcObj, outSnake2, totT);
			stop_s=TimeStamp("Total: ", stop_s);
			cout << endl;
			totT=totT+1;
		}

		PrintRSVSSnake(outSnake, testSnake, totT, testTriangle,
				triMesh, triRSVS, voluMesh, nVoluZone, ii);

		stop_s=clock();
		cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << "ms" << endl;
		testSnake.displight();
		outSnake.PrintMesh(testSnake.snakeconn);
		outSnake.PrintSnakeInternalPts(testSnake);



	} catch (exception const& ex) { 
		cerr << "Exception: " << ex.what() <<endl; 
		return -1;
	} 
	return(errTest);
}

int Test_snakeRSVS(){
	int nVoluZone;
	snake testSnake;
	mesh snakeMesh, triMesh,voluMesh;
	triangulation testTriangle,triRSVS;
	vector<int> dims;
	const char *fileToOpen;
	tecplotfile outSnake, outSnake2;
	double totT=0.0;
	vector<double> dt;
	vector<int> isImpact;
	int start_s,stop_s,ii;
	//bool errFlag;
	int errTest=0;
	SQPcalc calcObj;
	

	dims.assign(3,0);
	dims[0]=2;dims[1]=3;dims[2]=1;
	// try {
		fileToOpen="..\\TESTOUT\\TestSnakeRSVS.plt";
		outSnake.OpenFile(fileToOpen);
		fileToOpen="..\\TESTOUT\\TestSnakeRSVS_snake.plt";
		outSnake2.OpenFile(fileToOpen);
		
		errTest+=snakeMesh.read("..\\TESTOUT\\mesh203010.dat");
		PrepareMultiLvlSnake(snakeMesh,voluMesh,testSnake,dims,triRSVS);

		outSnake.PrintMesh(*(testSnake.snakemesh));
		outSnake.PrintMesh(voluMesh);
		nVoluZone=outSnake.ZoneNum();
		
		// SpawnAtVertex(testSnake,1022);
		// SpawnAtVertex(testSnake,674);
		// SpawnAtVertex(testSnake,675);
		// SpawnAtVertex(testSnake,728);
		// SpawnAtVertex(testSnake,729);
		SpawnAtVertex(testSnake,731);
		testSnake.displight();

		start_s=clock();
		testSnake.PrepareForUse();
		triRSVS.PrepareForUse();

		TriangulateSnake(testSnake,triRSVS);
		triRSVS.PrepareForUse();
		triRSVS.CalcTriVertPos();
		MaintainTriangulateSnake(triRSVS);
		
		for(ii=0;ii<100;++ii){
			cout << ii << " ";
			PrintRSVSSnake(outSnake, testSnake, totT, testTriangle,
				triMesh, triRSVS, voluMesh, nVoluZone, ii);
			stop_s=clock();
			Test_stepalgoRSVS(testSnake,triRSVS, dt, isImpact, calcObj, outSnake2, totT);
			stop_s=TimeStamp("Total: ", stop_s);
			cout << endl;
			totT=totT+1;
		}

		PrintRSVSSnake(outSnake, testSnake, totT, testTriangle,
				triMesh, triRSVS, voluMesh, nVoluZone, ii);

		stop_s=clock();
		cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << "ms" << endl;
		testSnake.displight();


	// } catch (exception const& ex) { 
	// 	cerr << "Exception: " << ex.what() <<endl; 
	// 	return -1;
	// } 
	return(errTest);

}

void repositionatquarteredge(snake &snakein){
	/*Assumes there are only 2 snakes per edge*/

	for(int ii=0; ii<int(snakein.snaxs.size()); ++ii){
		if(snakein.snaxs(ii)->orderedge>1) {
			snakein.snaxs[ii].d=0.75;
		} else {
			snakein.snaxs[ii].d=0.25;
		}
	}


	snakein.UpdateCoord();
	snakein.PrepareForUse();
}
void Test_mathRSVS_FD(snake &testSnake,triangulation &RSVStri , vector<double> &dt,
	vector<int> &isImpact, SQPcalc &calcObj, tecplotfile &outSnake2, double totT);
int Test_RSVSalgo_singlevol(){

	// int nVoluZone, ii;

	snake testSnake;
	mesh snakeMesh,  voluMesh, triMesh, triMesh2, triMesh3, triMesh4;
	// mesh triMesh;
	triangulation testTriangle,triRSVS;
	vector<int> dims;
	const char *fileToOpen;
	tecplotfile outSnake, outSnake2;
	// double totT=0.0;
	// vector<double> dt;
	// vector<int> isImpact;
	int start_s,stop_s;
	//bool errFlag;
	int errTest=0;
	SQPcalc calcObj;
	int ii;
	double totT=0.0;
	vector<double> dt;
	vector<int> isImpact;
	

	dims.assign(3,0);
	dims[0]=1;dims[1]=1;dims[2]=3;
	try {
		fileToOpen="..\\TESTOUT\\TestAlgoRSVSstep.plt";

		outSnake.OpenFile(fileToOpen);

		fileToOpen="..\\TESTOUT\\TestAlgoRSVSstep_snake.plt";
		outSnake2.OpenFile(fileToOpen);
		errTest+=snakeMesh.read("..\\TESTOUT\\mesh203010.dat");
		
		PrepareMultiLvlSnake(snakeMesh,voluMesh,testSnake,dims,triRSVS);
		for(ii=0;ii<voluMesh.volus.size();++ii){
			voluMesh.volus[ii].target=0.8;
		}
		voluMesh.volus[1].target=0.1;
		voluMesh.volus[2].target=0.3;
		// voluMesh.volus[0].target=0.001;//0.05;//0.0001;
		// voluMesh.volus[1].target=0.05;//0.05;
		// voluMesh.volus[2].target=0.001;//0.05;//0.0001;
		// voluMesh.volus[3].target=0.3;//0.05;//0.0001;
		// voluMesh.volus[4].target=0.3;//0.05;//0.0001;
		// voluMesh.volus[3].target=1;//0.05;//0.0001;
		// voluMesh.volus[4].target=1;//0.05;
		// voluMesh.volus[5].target=1;//0.05;//0.0001;
		// voluMesh.volus[3].target=1;//0.05;//0.0001;
		voluMesh.PrepareForUse();
		outSnake.PrintMesh(*(testSnake.snakemesh));
		outSnake.PrintMesh(voluMesh);
		int nVoluZone;
		nVoluZone=outSnake.ZoneNum();
		
		start_s=clock();

		SpawnRSVS(testSnake,1);
		testSnake.PrepareForUse();

		triRSVS.PrepareForUse();

		TriangulateSnake(testSnake,triRSVS);
		triRSVS.PrepareForUse();
		triRSVS.CalcTriVertPos();
		MaintainTriangulateSnake(triRSVS);

		for(ii=0;ii<300;++ii){
			cout << ii << " ";
			// if (ii%2==0){
				PrintRSVSSnake(outSnake, testSnake, totT, testTriangle,
					triMesh, triRSVS, voluMesh, nVoluZone, ii);
			// }
			stop_s=clock();
			Test_stepalgoRSVS(testSnake,triRSVS, dt, isImpact, calcObj, outSnake2, totT);
			stop_s=TimeStamp("Total: ", stop_s);
			cout << endl;
			totT=totT+1;
		}

		PrintRSVSSnake(outSnake, testSnake, totT, testTriangle,
				triMesh, triRSVS, voluMesh, nVoluZone, ii);

		stop_s=clock();
		cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << "ms" << endl;
		testSnake.displight();
		// outSnake.PrintMesh(testSnake.snakeconn);
		// outSnake.PrintSnakeInternalPts(testSnake);
		
		repositionatquarteredge(testSnake);
		triRSVS.CalcTriVertPos();
		triRSVS.PrepareForUse();
		triMesh2= TriarrayToMesh(triRSVS, triRSVS.intertri);
		FILE *fid;
		fid = fopen("..\\TESTOUT\\triintertestbig.dat","w");
		triMesh2.write(fid);//(triMesh,3,totT);
		fclose(fid);
		/*fid = fopen("..\\TESTOUT\\snakeconnouttri.dat","w");
		triMesh3= TriarrayToMesh(triRSVS, triRSVS.dynatri);
		// testSnake.snakeconn.write(fid);//(triMesh,3,totT);
		triMesh3.write(fid);//(triMesh,3,totT);
		fclose(fid);

		fid = fopen("..\\TESTOUT\\snakeconnout.dat","w");		
		testSnake.snakeconn.write(fid);//(triMesh,3,totT);
		fclose(fid);

		fid = fopen("..\\TESTOUT\\stattri.dat","w");
		triMesh4= TriarrayToMesh(triRSVS, triRSVS.stattri);
		triMesh4.write(fid);//(triMesh,3,totT);
		fclose(fid);*/

	} catch (exception const& ex) { 
		cerr << "Exception: " << ex.what() <<endl; 
		return -1;
	} 
	return(errTest);
}
 
int Test_snakeRSVS_singlevol(){
	int nVoluZone;
	snake testSnake;
	mesh snakeMesh, triMesh,voluMesh;
	triangulation testTriangle,triRSVS;
	vector<int> dims;
	const char *fileToOpen;
	tecplotfile outSnake, outSnake2;
	double totT=0.0;
	vector<double> dt;
	vector<int> isImpact;
	int start_s,stop_s,ii;
	//bool errFlag;
	int errTest=0;
	SQPcalc calcObj;
	

	dims.assign(3,0);
	dims[0]=1;dims[1]=1;dims[2]=1;
	// try {
		fileToOpen="..\\TESTOUT\\TestSnakeRSVS.plt";
		outSnake.OpenFile(fileToOpen);
		fileToOpen="..\\TESTOUT\\TestSnakeRSVS_snake.plt";
		outSnake2.OpenFile(fileToOpen);

		errTest+=snakeMesh.read("..\\TESTOUT\\mesh203010.dat");
		
		PrepareMultiLvlSnake(snakeMesh,voluMesh,testSnake,dims,triRSVS);
		voluMesh.volus[0].target=0.01;
		voluMesh.volus.PrepareForUse();
		outSnake.PrintMesh(*(testSnake.snakemesh));
		outSnake.PrintMesh(voluMesh);
		nVoluZone=outSnake.ZoneNum();
		
		// SpawnAtVertex(testSnake,1022);
		// SpawnAtVertex(testSnake,674);
		// SpawnAtVertex(testSnake,675);
		// SpawnAtVertex(testSnake,728);
		// SpawnAtVertex(testSnake,729);
		SpawnAtVertex(testSnake,731);
		testSnake.displight();

		start_s=clock();
		testSnake.PrepareForUse();
		triRSVS.PrepareForUse();

		TriangulateSnake(testSnake,triRSVS);
		triRSVS.PrepareForUse();
		triRSVS.CalcTriVertPos();
		MaintainTriangulateSnake(triRSVS);
		for(ii=0;ii<25;++ii){
			cout << ii << " ";
			PrintRSVSSnake(outSnake, testSnake, totT, testTriangle,
				triMesh, triRSVS, voluMesh, nVoluZone, ii);
			stop_s=clock();
			Test_stepalgoRSVS(testSnake,triRSVS, dt, isImpact, calcObj, outSnake2, totT);
			stop_s=TimeStamp("Total: ", stop_s);
			cout << endl;
			totT=totT+1;
		}

		PrintRSVSSnake(outSnake, testSnake, totT, testTriangle,
				triMesh, triRSVS, voluMesh, nVoluZone, ii);

		stop_s=clock();
		cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << "ms" << endl;
		testSnake.displight();


	// } catch (exception const& ex) { 
	// 	cerr << "Exception: " << ex.what() <<endl; 
	// 	return -1;
	// } 
	return(errTest);

}

int Test_MeshOrient(){

	// int nVoluZone, ii;

	snake testSnake;
	mesh snakeMesh,  voluMesh;
	// mesh triMesh;
	triangulation testTriangle,triRSVS;
	vector<int> dims;
	const char *fileToOpen;
	tecplotfile outSnake;
	// double totT=0.0;
	// vector<double> dt;
	// vector<int> isImpact;

	//bool errFlag;
	int errTest=0;
	SQPcalc calcObj;
	FILE *fid;
	

	dims.assign(3,0);
	dims[0]=1;dims[1]=3;dims[2]=1;
	
		fileToOpen="..\\TESTOUT\\TestAlgoRSVSstep.plt";

		outSnake.OpenFile(fileToOpen);
		errTest+=snakeMesh.read("..\\TESTOUT\\mesh6612.dat");
		
		PrepareMultiLvlSnake(snakeMesh,voluMesh,testSnake,dims,triRSVS);
		voluMesh.volus[0].target=1;//0.05;//0.0001;
		voluMesh.volus[1].target=1;//0.05;
		voluMesh.volus[2].target=1;//0.05;//0.0001;
		// voluMesh.volus[3].target=1;//0.05;//0.0001;
		voluMesh.PrepareForUse();
		outSnake.PrintMesh(*(testSnake.snakemesh));
		outSnake.PrintMesh(voluMesh);

		fileToOpen = "..\\TESTOUT\\volumesh6612.dat";
		fid=fopen(fileToOpen,"w");
		if(fid!=NULL){
			voluMesh.PrepareForUse();
			voluMesh.SetBorders();
			voluMesh.OrientFaces();
			voluMesh.write(fid);
		}
	return(errTest);
}

void Test_stepalgoRSVS(snake &testSnake,triangulation &RSVStri , vector<double> &dt,
	vector<int> &isImpact, SQPcalc &calcObj, tecplotfile &outSnake2, double totT){
	int start_s;

	 
	start_s=clock();

	// calcObj.Print2Screen(1);
	// calcObj.Print2Screen(2);
	// CalculateSnakeVel(testSnake);
	// Small step away from edge without crossover.
	// Need to develop that.
	// Check if impact detect crossovers
	start_s=TimeStamp(" triangulate:", start_s);
	// calcObj.limLag=10000.0;
	calcObj.CalculateTriangulation(RSVStri);
	calcObj.ReturnConstrToMesh(RSVStri);
	start_s=TimeStamp(" deriv:", start_s);
	calcObj.CheckAndCompute();
	calcObj.ReturnVelocities(RSVStri);
	start_s=TimeStamp(" solve:", start_s);
	
	calcObj.Print2Screen();
	// calcObj.Print2Screen(2);
	CalculateNoNanSnakeVel(testSnake);
	outSnake2.PrintSnake(testSnake, 1, totT);
	testSnake.CalculateTimeStep(dt,0.9);
	testSnake.UpdateDistance(dt,0.9);
	testSnake.UpdateCoord();
	testSnake.PrepareForUse();

	SnakeConnectivityUpdate(testSnake, isImpact);
	MaintainTriangulateSnake(RSVStri);
	

}

void Test_mathRSVS_FD(snake &testSnake,triangulation &RSVStri , vector<double> &dt,
	vector<int> &isImpact, SQPcalc &calcObj, tecplotfile &outSnake2, double totT){
	int start_s;

	SQPcalc calcObj2;
	double fdStep=1e-3;
	start_s=clock();

	// calcObj.Print2Screen(1);
	// calcObj.Print2Screen(2);
	// CalculateSnakeVel(testSnake);
	// Small step away from edge without crossover.
	// Need to develop that.
	// Check if impact detect crossovers
	start_s=TimeStamp(" triangulate:", start_s);
	// calcObj.limLag=10000.0;
	// for(int ii=0; ii<testSnake.snakemesh->verts.size();++ii){
	// 	for(int jj=0; jj<3;++jj){
	// 		testSnake.snakemesh->verts[ii].coord[jj] +=1;
	// 	}
	// }
	// testSnake.PrepareForUse();
	// testSnake.UpdateCoord();
	// testSnake.PrepareForUse();
	calcObj.CalculateTriangulation(RSVStri, 2);
	testSnake.snaxs[137].d += fdStep;
	testSnake.PrepareForUse();
	testSnake.UpdateCoord();
	testSnake.PrepareForUse();
	calcObj2.CalculateTriangulation(RSVStri);

	cout << endl; 
	for(int ii = 0; ii<calcObj.numConstr();++ii){
		cout << calcObj2.constr[ii]- calcObj.constr[ii]  << " " << 
			(calcObj2.constr[ii]- calcObj.constr[ii])/fdStep << endl;
	}

	calcObj.ReturnConstrToMesh(RSVStri);
	start_s=TimeStamp(" deriv:", start_s);
	calcObj.CheckAndCompute();
	calcObj.ReturnVelocities(RSVStri);
	start_s=TimeStamp(" solve:", start_s);
	outSnake2.PrintSnake(testSnake, 1, totT);
	
	calcObj.Print2Screen();
	// calcObj.Print2Screen(2);
	CalculateNoNanSnakeVel(testSnake);
	testSnake.CalculateTimeStep(dt,0.5);
	testSnake.UpdateDistance(dt,0.34);
	testSnake.UpdateCoord();
	testSnake.PrepareForUse();

	SnakeConnectivityUpdate(testSnake, isImpact);
	MaintainTriangulateSnake(RSVStri);
	

}

