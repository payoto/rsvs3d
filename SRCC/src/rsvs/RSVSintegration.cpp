#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <vector>
#include <string>


#include<tuple>
#include "RSVSclass.hpp"
#include "RSVSintegration.hpp"
#include "snake.hpp"
#include "snakeengine.hpp"
#include "parameters.hpp"
#include "voxel.hpp"
#include "meshrefinement.hpp"
#include "RSVScalc.hpp"
#include "RSVSalgorithm.hpp"
#include "postprocessing.hpp"

#include "filesystem.hpp"

int SAFE_ALGO_TestConn(snake &snakein){
	int ret=0;

	if (snakein.Check3D()){
		#ifdef SAFE_ALGO
		ret = snakein.snakeconn.TestConnectivityBiDir();
		#endif //SAFE_ALGO
	}

	return(ret);
}

void SnakeConnectivityUpdate_legacy(snake &snakein,  vector<int> &isImpact){


	int start_s;

	start_s=clock();

	

	start_s=TimeStamp("position: ", start_s);

	snakein.SnaxImpactDetection(isImpact);
	MergeAllContactVertices(snakein, isImpact);
	snakein.PrepareForUse();


	start_s=TimeStamp("Merge: ", start_s);

	CleanupSnakeConnec(snakein);
	

	start_s=TimeStamp("Clean: ", start_s);
	snakein.SnaxImpactDetection(isImpact);
	SpawnArrivedSnaxels(snakein,isImpact);


	

	start_s=TimeStamp("Spawn: ", start_s);

	snakein.SnaxImpactDetection(isImpact);
	snakein.SnaxAlmostImpactDetection(isImpact, 0.01);
	snakein.PrepareForUse();
	SAFE_ALGO_TestConn(snakein);
	MergeAllContactVertices(snakein, isImpact);
	snakein.PrepareForUse();
	SAFE_ALGO_TestConn(snakein);
	

	start_s=TimeStamp("Impact: ", start_s);

	CleanupSnakeConnec(snakein);

	SAFE_ALGO_TestConn(snakein);
	snakein.OrientFaces();

	start_s=TimeStamp("Clean: ", start_s);
}

void SnakeConnectivityUpdate_robust(snake &snakein,  vector<int> &isImpact){
	/*
	Performs the snake step except the movement of the snake.

	This one performs it in two steps:
	 1) Impact Merge Clean
	 2) Impact Spawn Impact Merge Clean

	This function might be better in snakeengine.cpp
	*/
	double impactAlmostRange = 0.2;

	int start_s, start_f;
	start_f=clock();

	//===============================
	// Impact on edge
	// ======================
	// Impact
	SAFE_ALGO_TestConn(snakein);
	snakein.SnaxImpactDetection(isImpact);
	snakein.PrepareForUse();
	start_s=TimeStamp("Impact: ", start_f);
	// ======================
	// Merge
	MergeAllContactVertices(snakein, isImpact);
	snakein.PrepareForUse();
	start_s=TimeStamp("Merge: ", start_s);
	// ======================
	// Clean
	CleanupSnakeConnec(snakein);
	snakein.PrepareForUse();
	SAFE_ALGO_TestConn(snakein);
	snakein.OrientFaces();
	start_s=TimeStamp("Clean: ", start_s);


	//===============================
	// Spawn
	// ======================
	// Impact
	SAFE_ALGO_TestConn(snakein);
	snakein.SnaxImpactDetection(isImpact);
	snakein.SnaxAlmostImpactDetection(isImpact, impactAlmostRange);
	snakein.PrepareForUse();
	start_s=TimeStamp("Impact: ", start_s);
	// ======================
	// Spawn
	SAFE_ALGO_TestConn(snakein);
	snakein.SnaxImpactDetection(isImpact);
	SpawnArrivedSnaxels(snakein,isImpact);
	start_s=TimeStamp("Spawn: ", start_s);
	// ======================
	// Impact
	SAFE_ALGO_TestConn(snakein);
	snakein.SnaxImpactDetection(isImpact);
	snakein.PrepareForUse();
	start_s=TimeStamp("Impact: ", start_s);
	// ======================
	// Merge
	MergeAllContactVertices(snakein, isImpact);
	snakein.PrepareForUse();
	start_s=TimeStamp("Merge: ", start_s);
	// ======================
	// Clean
	CleanupSnakeConnec(snakein);
	snakein.PrepareForUse();
	SAFE_ALGO_TestConn(snakein);
	snakein.OrientFaces();
	start_s=TimeStamp("Clean: ", start_s);

	TimeStamp(" - Connec Update: ", start_f);
}

void SnakeConnectivityUpdate(snake &snakein,  vector<int> &isImpact,
	double impactAlmostRange){
	/*
	Performs the snake step except the movement of the snake.
	This one performs it in a 'single' step:
	 Impact Spawn Impact Merge Clean

	This function might be better in snakeengine.cpp
	*/

	int start_s, start_f;
	start_f=clock();

	//===============================
	// Spawn
	// ======================
	// Impact
	SAFE_ALGO_TestConn(snakein);
	snakein.SnaxImpactDetection(isImpact);
	snakein.SnaxAlmostImpactDetection(isImpact, impactAlmostRange);
	snakein.PrepareForUse();
	start_s=TimeStamp("Impact: ", start_f);
	// ======================
	// Spawn
	SAFE_ALGO_TestConn(snakein);
	snakein.SnaxImpactDetection(isImpact);
	SpawnArrivedSnaxels(snakein,isImpact);
	start_s=TimeStamp("Spawn: ", start_s);
	// ======================
	// Impact
	SAFE_ALGO_TestConn(snakein);
	snakein.SnaxImpactDetection(isImpact);
	snakein.PrepareForUse();
	start_s=TimeStamp("Impact: ", start_s);
	// ======================
	// Merge
	MergeAllContactVertices(snakein, isImpact);
	snakein.PrepareForUse();
	start_s=TimeStamp("Merge: ", start_s);
	// ======================
	// Clean
	CleanupSnakeConnec(snakein);
	snakein.PrepareForUse();
	SAFE_ALGO_TestConn(snakein);
	snakein.OrientFaces();
	start_s=TimeStamp("Clean: ", start_s);

	TimeStamp(" - Connec Update: ", start_f);
}

void SnakeConnectivityUpdate_2D(snake &snakein,  vector<int> &isImpact){
	/*
	Performs the snake step except the movement of the snake.
	This one performs it in a 'single' step:
	 Impact Spawn Impact Merge Clean

	This function might be better in snakeengine.cpp
	*/
	double impactAlmostRange = 0.2;

	int start_s, start_f;
	start_f=clock();


	//===============================
	// Spawn
	// ======================
	// Impact
	SAFE_ALGO_TestConn(snakein);
	snakein.SnaxImpactDetection(isImpact);
	snakein.SnaxAlmostImpactDetection(isImpact, impactAlmostRange);
	snakein.PrepareForUse();
	start_s=TimeStamp("Impact: ", start_f);
	// ======================
	// Spawn
	SAFE_ALGO_TestConn(snakein);
	snakein.SnaxImpactDetection(isImpact);
	SpawnArrivedSnaxels(snakein,isImpact);
	start_s=TimeStamp("Spawn: ", start_s);
	// ======================
	// Impact
	SAFE_ALGO_TestConn(snakein);
	snakein.SnaxImpactDetection(isImpact);
	snakein.PrepareForUse();
	start_s=TimeStamp("Impact: ", start_s);
	// ======================
	// Merge
	MergeAllContactVertices(snakein, isImpact);
	snakein.PrepareForUse();
	start_s=TimeStamp("Merge: ", start_s);
	// ======================
	// Clean
	CleanupSnakeConnec(snakein);
	snakein.PrepareForUse();
	SAFE_ALGO_TestConn(snakein);
	snakein.OrientFaces();
	start_s=TimeStamp("Clean: ", start_s);

	TimeStamp(" - Connec Update: ", start_f);
}

int TimeStamp(const char* str,int start_s){
	int stop_s=clock();
	#ifdef TIME_EXEC
	cout << str << " " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << "ms; ";
	#endif
	return(stop_s);
}

// ====================
// integrate
// 		prepare
// ====================

void integrate::Prepare(integrate::RSVSclass &RSVSobj){
	// likely inputs (now in RSVSclass)
	// param::parameters paramconf;
	// mesh snakeMesh;
	// mesh voluMesh;
	// snake rsvsSnake;
	// triangulation rsvsTri;
	// tecplotfile outSnake;
	// Locally defined
	param::parameters origconf;

	origconf = RSVSobj.paramconf;
	RSVSobj.paramconf.PrepareForUse();

	integrate::prepare::Mesh(RSVSobj.paramconf.grid, RSVSobj.snakeMesh, RSVSobj.voluMesh);
	integrate::prepare::Snake(RSVSobj.paramconf.snak, RSVSobj.paramconf.rsvs,
		RSVSobj.snakeMesh, RSVSobj.voluMesh, RSVSobj.rsvsSnake);
	integrate::prepare::Triangulation(RSVSobj.snakeMesh, RSVSobj.rsvsSnake, RSVSobj.rsvsTri);
	integrate::prepare::Output(RSVSobj.paramconf, origconf, RSVSobj.outSnake, RSVSobj.logFile,
		RSVSobj.coutFile, RSVSobj.cerrFile);
}

void integrate::prepare::Mesh(
	const param::grid &gridconf,
	mesh &snakeMesh,
	mesh &voluMesh
	){
	/*prepares the snake and volume meshes gor the RSVS process*/
	// Local declaration
	int ii;
	auto gridSize = gridconf.voxel.gridsizebackground;
	vector<int> elmMapping, backgroundGrid;
	RSVScalc calcVolus;

	backgroundGrid.reserve(int(gridSize.size()));
	for (int i = 0; i < int(gridSize.size()); ++i)
	{
		backgroundGrid.push_back(gridSize[i]);
		gridSize[i] = gridSize[i] * gridconf.voxel.gridsizesnake[i];
	}

	// Initial build of the grid
	BuildBlockGrid(gridSize, snakeMesh);
	snakeMesh.Scale(gridconf.voxel.domain);
	snakeMesh.PrepareForUse();
	snakeMesh.OrientFaces();

	// map elements to coarse grid
	if(snakeMesh.WhatDim()==3){
		for (ii=0;ii<snakeMesh.volus.size();++ii){
			elmMapping.push_back(1);
		}
	} else if (snakeMesh.WhatDim()==2){
		for (ii=0;ii<snakeMesh.surfs.size();++ii){
			elmMapping.push_back(1);
		}
	} else { 
		throw invalid_argument("Incorrect dimension");
	}
	CartesianMapping(snakeMesh,  elmMapping, backgroundGrid);
	CoarsenMesh(snakeMesh,voluMesh,elmMapping);
	snakeMesh.AddParent(&voluMesh,elmMapping);

	voluMesh.PrepareForUse();
	voluMesh.OrientFaces();
	calcVolus.CalculateMesh(voluMesh);
	calcVolus.ReturnConstrToMesh(voluMesh,&volu::volume);
}

void integrate::prepare::Snake(
	const param::snaking &snakconf, 
	const param::rsvs &rsvsconf, 
	mesh &snakeMesh, // non const as it is passed to the snake as a pointer
	mesh &voluMesh,
	snake &rsvsSnake
	){

	// go through the rsvs conf figuring out which fill option to use.
	int nElms;
	if(voluMesh.WhatDim()==3){
		nElms = voluMesh.volus.size();
		if(rsvsconf.filefill.active){
			voluMesh.LoadTargetFill(rsvsconf.filefill.fill);
		} else if(rsvsconf.makefill.active){
			// TODO add a fill builder
		} else if(rsvsconf.cstfill.active){
			for(int i=0; i< nElms; ++i){
				voluMesh.volus[i].target=rsvsconf.cstfill.fill;
			}
		}
	} else if(voluMesh.WhatDim()==2){
		nElms = voluMesh.surfs.size();
		if(rsvsconf.filefill.active){
			voluMesh.LoadTargetFill(rsvsconf.filefill.fill);
		} else if(rsvsconf.makefill.active){
			// TODO add a fill builder
		} else if(rsvsconf.cstfill.active){
			for(int i=0; i< nElms; ++i){
				voluMesh.surfs[i].target=rsvsconf.cstfill.fill;
			}
		}
	}
	voluMesh.PrepareForUse();


	rsvsSnake.snakemesh = &snakeMesh;
	SpawnRSVS(rsvsSnake,snakconf.initboundary);
	rsvsSnake.PrepareForUse();
}

void integrate::prepare::Triangulation(
	mesh &snakeMesh,
	snake &rsvsSnake,
	triangulation &rsvsTri
	){

	rsvsTri.PrepareForUse();
	TriangulateMesh(snakeMesh,rsvsTri);
	rsvsTri.PrepareForUse();
	TriangulateSnake(rsvsSnake,rsvsTri);
	rsvsTri.PrepareForUse();
	rsvsTri.CalcTriVertPos();
	rsvsTri.PrepareForUse();
	MaintainTriangulateSnake(rsvsTri);
	rsvsTri.PrepareForUse();
}

void integrate::prepare::Output(
	const param::parameters &paramconf,
	const param::parameters &origconf,
	tecplotfile &outSnake,
	std::ofstream &logFile,
	std::ofstream &coutFile,
	std::ofstream &cerrFile
	){
	std::string outSnakeName;
	
	if (paramconf.files.ioout.outdir.size()!=0){
		#ifdef USE_CSTD_FILESYSTEM
		// use c++17 std library filesystem
		const std::filesystem::path pathout=paramconf.files.ioout.outdir;
		std::filesystem::create_directories(pathout);
		#else
		// else use boost
		boost::filesystem::create_directories(paramconf.files.ioout.outdir);
		#endif
	}
	outSnakeName = paramconf.files.ioout.outdir + "/";
	outSnakeName += "rsvs3D_" + paramconf.files.ioout.pattern + ".plt";
	outSnake.OpenFile(outSnakeName.c_str());

	outSnakeName =  paramconf.files.ioout.outdir + "/";
	outSnakeName += "config_call_" + paramconf.files.ioout.pattern + ".json";
	param::io::writeflat(outSnakeName, origconf);
	outSnakeName =  paramconf.files.ioout.outdir + "/";
	outSnakeName += "config_active_" + paramconf.files.ioout.pattern + ".json";
	param::io::writeflat(outSnakeName, paramconf);
	// Open a text log file
	outSnakeName =  paramconf.files.ioout.outdir + "/";
	outSnakeName += "convergence_" + paramconf.files.ioout.pattern + ".log";
	logFile.open(outSnakeName);
	logFile.precision(16);
	logFile << std::scientific;

	if (paramconf.files.ioout.redirectcout){
		// auto ;
		outSnakeName =  paramconf.files.ioout.outdir + "/";
		outSnakeName += "cout_" + paramconf.files.ioout.pattern + ".txt";
		coutFile.open(outSnakeName);
		std::cout.rdbuf(coutFile.rdbuf());
	}
	if (paramconf.files.ioout.redirectcerr){
		// auto ;
		outSnakeName =  paramconf.files.ioout.outdir + "/";
		outSnakeName += "cerr_" + paramconf.files.ioout.pattern + ".txt";
		cerrFile.open(outSnakeName);
		std::cerr.rdbuf(cerrFile.rdbuf());
	}

}


// ====================
// integrate
// 		execute
// ====================

void integrate::execute::All(integrate::RSVSclass &RSVSobj){


	auto coutbuff=std::cout.rdbuf();
	auto cerrbuff=std::cout.rdbuf();

	integrate::Prepare(RSVSobj);
	auto iterateInfo=integrate::execute::RSVSiterate(RSVSobj);

	integrate::execute::PostProcessing(RSVSobj,
		iterateInfo.timeT, iterateInfo.nVoluZone, iterateInfo.stepNum);
	std::cout.rdbuf(coutbuff);
	std::cerr.rdbuf(cerrbuff);
	// std::cout << std::endl <<  " cout Buffer restored" << std::endl;
	// std::cerr << " cerr Buffer restored" << std::endl;
}

integrate::iteratereturns integrate::execute::RSVSiterate(integrate::RSVSclass &RSVSobj){

	vector<double> dt;
	vector<int> isImpact;
	int stepNum, maxStep, nVoluZone;
	double totT=0;

	auto start_s=clock();
	// for n Steps

	RSVSobj.voluMesh.PrepareForUse();
	RSVSobj.outSnake.PrintMesh(RSVSobj.snakeMesh);
	RSVSobj.outSnake.PrintMesh(RSVSobj.voluMesh);
	nVoluZone=RSVSobj.outSnake.ZoneNum();

	maxStep = RSVSobj.paramconf.snak.maxsteps;
	for(stepNum=0; stepNum<maxStep; ++stepNum){
		start_s=clock(); 
		// calcObj.limLag=10000.0;
		std::cout << std::endl << "Step " << stepNum << " ";
		RSVSobj.calcObj.CalculateTriangulation(RSVSobj.rsvsTri);
		start_s=TimeStamp(" deriv:", start_s);
		RSVSobj.calcObj.CheckAndCompute(
			RSVSobj.paramconf.rsvs.solveralgorithm);
		RSVSobj.calcObj.ReturnConstrToMesh(RSVSobj.rsvsTri);
		RSVSobj.calcObj.ReturnVelocities(RSVSobj.rsvsTri);
		start_s=TimeStamp(" solve:", start_s);
		
		CalculateNoNanSnakeVel(RSVSobj.rsvsSnake);

		integrate::execute::Logging(RSVSobj,
			totT, nVoluZone, stepNum);

		RSVSobj.rsvsSnake.CalculateTimeStep(dt,
			RSVSobj.paramconf.snak.snaxtimestep);
		RSVSobj.rsvsSnake.UpdateDistance(dt,
			RSVSobj.paramconf.snak.snaxdiststep);
		totT += *max_element(dt.begin(), dt.end());
		RSVSobj.rsvsSnake.UpdateCoord();
		RSVSobj.rsvsSnake.PrepareForUse();

		SnakeConnectivityUpdate(RSVSobj.rsvsSnake, isImpact,
			RSVSobj.paramconf.snak.multiarrivaltolerance);
		start_s=clock();
		MaintainTriangulateSnake(RSVSobj.rsvsTri);
		start_s=TimeStamp(" triangulate:", start_s);
	}
	std::cout << std::endl << "RSVS iteration finished" << std::endl;
	integrate::iteratereturns retStruct(nVoluZone,stepNum, totT);
	return(retStruct);
}

void integrate::execute::Logging(integrate::RSVSclass &RSVSobj,
	double totT, int nVoluZone, int stepNum){
	// Simple function which directs to the correct output

	if (0 < RSVSobj.paramconf.files.ioout.logginglvl){
		// calcObj logging outputs different amounts of data
		// depending.
		RSVSobj.logFile << "> step" << stepNum << " :," ;
		RSVSobj.logFile << totT << endl;
		integrate::execute::logging::Log(
			RSVSobj.logFile, RSVSobj.calcObj,
			RSVSobj.paramconf.files.ioout.logginglvl
			);
	}

	if (1 < RSVSobj.paramconf.files.ioout.logginglvl){
		integrate::execute::logging::Snake(
			RSVSobj.outSnake, RSVSobj.rsvsSnake,
			RSVSobj.voluMesh, totT, nVoluZone);
	}

	if (2 < RSVSobj.paramconf.files.ioout.logginglvl){
		integrate::execute::logging::FullTecplot(
			RSVSobj.outSnake, RSVSobj.rsvsSnake,
			RSVSobj.rsvsTri, RSVSobj.voluMesh,
			totT, nVoluZone, stepNum);
	}
}

void integrate::execute::PostProcessing(integrate::RSVSclass &RSVSobj,
	double totT, int nVoluZone, int stepNum){

	if (0 < RSVSobj.paramconf.files.ioout.outputlvl){
		RSVSobj.logFile << "> final step" << stepNum << " :," ;
			RSVSobj.logFile << totT << endl;
		integrate::execute::postprocess::Log(
			RSVSobj.logFile, RSVSobj.calcObj,
			RSVSobj.paramconf.files.ioout.logginglvl
			);
	}

	if (1 < RSVSobj.paramconf.files.ioout.logginglvl){
		integrate::execute::postprocess::Snake(
			RSVSobj.outSnake, RSVSobj.rsvsSnake,
			RSVSobj.voluMesh, totT, nVoluZone,
			RSVSobj.paramconf);
	}

	if (2 < RSVSobj.paramconf.files.ioout.logginglvl){
		integrate::execute::postprocess::FullTecplot(
			RSVSobj.outSnake, RSVSobj.rsvsSnake,
			RSVSobj.rsvsTri, RSVSobj.voluMesh,
			totT, nVoluZone, stepNum);
	}
}
// ====================
// integrate
// 		execute
// 			logging
// ====================

void integrate::execute::logging::Log(
	std::ofstream &logFile,
	RSVScalc &calcObj,
	int loglvl
	){

	// Make a logging function for 
	// volume convergence and velocity convergence
	calcObj.ConvergenceLog(logFile, loglvl);
}

void integrate::execute::logging::Snake(
	tecplotfile &outSnake, snake &rsvsSnake,
	mesh &voluMesh, double totT, int nVoluZone
	){

	outSnake.PrintVolumeDat(voluMesh,nVoluZone,1,totT);
	outSnake.PrintSnake(rsvsSnake, 2, totT);
}

void integrate::execute::logging::FullTecplot(
	tecplotfile &outSnake, snake &rsvsSnake,
	triangulation &rsvsTri, mesh &voluMesh,
	double totT, int nVoluZone, int stepNum){
	vector<int> vertList;
	int jj;
	if(rsvsSnake.snaxs.size()>0){
		//rsvsSnake.snakeconn.TightenConnectivity();
		outSnake.PrintMesh(rsvsSnake.snakeconn,1,totT);

		rsvsSnake.snakeconn.PrepareForUse();
		
		outSnake.PrintTriangulation(rsvsTri,&triangulation::dynatri,2,totT);
		if (stepNum==0){
			outSnake.PrintTriangulation(rsvsTri,&triangulation::dynatri,3,totT,3);
			outSnake.PrintTriangulation(rsvsTri,&triangulation::dynatri,4,totT,3);
			outSnake.PrintTriangulation(rsvsTri,&triangulation::dynatri,5,totT,3);
		}
		outSnake.PrintTriangulation(rsvsTri,&triangulation::intertri,3,totT,3);
		outSnake.PrintTriangulation(rsvsTri,&triangulation::trisurf,4,totT,3);
		if (int(rsvsTri.acttri.size())>0){
			outSnake.PrintTriangulation(rsvsTri,&triangulation::stattri,5,totT,3,rsvsTri.acttri);
		}
		
		vertList.clear();
		for(jj=0;jj<int(rsvsSnake.isMeshVertIn.size()); ++jj){
			if(rsvsSnake.isMeshVertIn[jj]){
				vertList.push_back(rsvsSnake.snakemesh->verts(jj)->index);
			}
		}
		if(int(rsvsSnake.isMeshVertIn.size())==0){
			vertList.push_back(rsvsSnake.snakemesh->verts(0)->index);
		}
		outSnake.PrintMesh(*(rsvsSnake.snakemesh),6,totT,4,vertList);
		outSnake.PrintVolumeDat(voluMesh,nVoluZone,7,totT);
		outSnake.PrintSnake(rsvsSnake, 8, totT);

	}
}

// ====================
// integrate
// 		execute
// 			postprocess
// ====================
void integrate::execute::postprocess::Log(
	std::ofstream &logFile,
	RSVScalc &calcObj,
	int loglvl
	){

	integrate::execute::logging::Log(
			logFile, calcObj, loglvl+2);
}

void integrate::execute::postprocess::Snake(
	tecplotfile &outSnake, snake &rsvsSnake,
	mesh &voluMesh, double totT, int nVoluZone,
	param::parameters &paramconf
	){

	string fileToOpen;

	integrate::execute::logging::Snake(
		outSnake, rsvsSnake,
		voluMesh, totT, nVoluZone);

	fileToOpen=paramconf.files.ioout.outdir + "/";
	fileToOpen += "VoluMesh_" + paramconf.files.ioout.pattern + ".msh";
	voluMesh.write(fileToOpen.c_str());

	fileToOpen=paramconf.files.ioout.outdir + "/";
	fileToOpen += "SnakeMesh_" + paramconf.files.ioout.pattern + ".msh";
	rsvsSnake.snakemesh->write(fileToOpen.c_str());

	fileToOpen=paramconf.files.ioout.outdir + "/";
	fileToOpen += "SnakeConn_" + paramconf.files.ioout.pattern + ".msh";
	rsvsSnake.snakeconn.write(fileToOpen.c_str());

	fileToOpen=paramconf.files.ioout.outdir + "/";
	fileToOpen += "Snake_" + paramconf.files.ioout.pattern + ".3snk";
	rsvsSnake.write(fileToOpen.c_str());

}

void integrate::execute::postprocess::FullTecplot(
	tecplotfile &outSnake, snake &rsvsSnake,
	triangulation &rsvsTri, mesh &voluMesh,
	double totT, int nVoluZone, int stepNum){

	integrate::execute::logging::FullTecplot(
		outSnake, rsvsSnake, rsvsTri, voluMesh,
		totT, nVoluZone, stepNum);
}

// ===================
// Tests
// ===================

int integrate::test::Prepare(){
	param::parameters paramconf, origconf;
	mesh snakeMesh;
	mesh voluMesh;
	snake rsvsSnake;
	triangulation rsvsTri;
	tecplotfile outSnake;
	std::ofstream logFile;
	std::ofstream coutFile;
	std::ofstream cerrFile;

	origconf = paramconf;
	paramconf.PrepareForUse();
	try {
		integrate::prepare::Mesh(paramconf.grid, snakeMesh, voluMesh);
	} catch (exception const& ex) { 
		cerr << "integrate::prepare::Mesh(paramconf.grid, snakeMesh, voluMesh);" << endl;
		cerr << "Exception: " << ex.what() <<endl; 
		return -1;
	} 
	
	try {
		integrate::prepare::Snake(paramconf.snak, paramconf.rsvs, snakeMesh,voluMesh,  rsvsSnake);
	} catch (exception const& ex) { 
		cerr << "integrate::prepare::Snake(paramconf.snak, snakeMesh, rsvsSnake);" << endl;
		cerr << "Exception: " << ex.what() <<endl; 
		return -1;
	} 

	try {
		integrate::prepare::Triangulation(snakeMesh, rsvsSnake, rsvsTri);
	} catch (exception const& ex) { 
		cerr << "integrate::prepare::Triangulation" << endl;
		cerr << "Exception: " << ex.what() <<endl; 
		return -1;
	} 

	try {
		integrate::prepare::Output(paramconf, origconf, outSnake,logFile,
			coutFile, cerrFile);
	} catch (exception const& ex) { 
		cerr << "integrate::prepare::Output" << endl;
		cerr << "Exception: " << ex.what() <<endl; 
		return -1;
	} 


	return(0);
}

int integrate::test::All(){

	integrate::RSVSclass RSVSobj;
	auto coutbuff=std::cout.rdbuf();
	auto cerrbuff=std::cout.rdbuf();

	integrate::Prepare(RSVSobj);
	auto iterateInfo=integrate::execute::RSVSiterate(RSVSobj);

	integrate::execute::PostProcessing(RSVSobj,
		iterateInfo.timeT, iterateInfo.nVoluZone, iterateInfo.stepNum);
	std::cout.rdbuf(coutbuff);
	std::cerr.rdbuf(cerrbuff);
	std::cout << std::endl <<  " cout Buffer restored" << std::endl;
	std::cerr << " cerr Buffer restored" << std::endl;
	return 0;
}

