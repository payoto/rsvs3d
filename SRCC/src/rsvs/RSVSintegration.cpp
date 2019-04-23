#include <iostream>
#include <fstream>
#include <cmath>
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
#include "meshprocessing.hpp"
#include "RSVScalc.hpp"
#include "RSVSalgorithm.hpp"
#include "postprocessing.hpp"
#include "warning.hpp"
#include "tetgenrsvs.hpp"

#include "filesystem.hpp"
  
int SAFE_ALGO_TestConn(snake &snakein){
	int ret=0;

	if (snakein.Check3D()){
		#ifdef SAFE_ALGO
		ret = snakein.snakeconn.TestConnectivityBiDir(__PRETTY_FUNCTION__);
		#endif //SAFE_ALGO
	}

	return(ret);
}

void SnakeConnectivityUpdate_legacy(snake &snakein,  vector<int> &isImpact){


	int start_s;

	start_s=clock();

	start_s=rsvs3d::TimeStamp("position: ", start_s);

	snakein.SnaxImpactDetection(isImpact);
	MergeAllContactVertices(snakein, isImpact);
	snakein.PrepareForUse();

	start_s=rsvs3d::TimeStamp("Merge: ", start_s);

	CleanupSnakeConnec(snakein);
	
	start_s=rsvs3d::TimeStamp("Clean: ", start_s);
	snakein.SnaxImpactDetection(isImpact);
	SpawnArrivedSnaxels(snakein,isImpact);

	start_s=rsvs3d::TimeStamp("Spawn: ", start_s);

	snakein.SnaxImpactDetection(isImpact);
	snakein.SnaxAlmostImpactDetection(isImpact, 0.01);
	snakein.PrepareForUse();
	SAFE_ALGO_TestConn(snakein);
	MergeAllContactVertices(snakein, isImpact);
	snakein.PrepareForUse();
	SAFE_ALGO_TestConn(snakein);
	
	start_s=rsvs3d::TimeStamp("Impact: ", start_s);

	CleanupSnakeConnec(snakein);

	SAFE_ALGO_TestConn(snakein);
	snakein.OrientFaces();

	start_s=rsvs3d::TimeStamp("Clean: ", start_s);
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
	start_s=rsvs3d::TimeStamp("Impact: ", start_f);
	// ======================
	// Merge
	MergeAllContactVertices(snakein, isImpact);
	snakein.PrepareForUse();
	start_s=rsvs3d::TimeStamp("Merge: ", start_s);
	// ======================
	// Clean
	CleanupSnakeConnec(snakein);
	snakein.PrepareForUse();
	SAFE_ALGO_TestConn(snakein);
	snakein.OrientFaces();
	start_s=rsvs3d::TimeStamp("Clean: ", start_s);


	//===============================
	// Spawn
	// ======================
	// Impact
	SAFE_ALGO_TestConn(snakein);
	snakein.SnaxImpactDetection(isImpact);
	snakein.SnaxAlmostImpactDetection(isImpact, impactAlmostRange);
	snakein.PrepareForUse();
	start_s=rsvs3d::TimeStamp("Impact: ", start_s);
	// ======================
	// Spawn
	SAFE_ALGO_TestConn(snakein);
	snakein.SnaxImpactDetection(isImpact);
	SpawnArrivedSnaxels(snakein,isImpact);
	start_s=rsvs3d::TimeStamp("Spawn: ", start_s);
	// ======================
	// Impact
	SAFE_ALGO_TestConn(snakein);
	snakein.SnaxImpactDetection(isImpact);
	snakein.PrepareForUse();
	start_s=rsvs3d::TimeStamp("Impact: ", start_s);
	// ======================
	// Merge
	MergeAllContactVertices(snakein, isImpact);
	snakein.PrepareForUse();
	start_s=rsvs3d::TimeStamp("Merge: ", start_s);
	// ======================
	// Clean
	CleanupSnakeConnec(snakein);
	snakein.PrepareForUse();
	SAFE_ALGO_TestConn(snakein);
	snakein.OrientFaces();
	start_s=rsvs3d::TimeStamp("Clean: ", start_s);

	rsvs3d::TimeStamp(" - Connec Update: ", start_f);
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
	// double stepLength = 1e-5;
	//===============================
	// Spawn
	// ======================
	// Impact
	SAFE_ALGO_TestConn(snakein);
	snakein.SnaxImpactDetection(isImpact);
	snakein.SnaxAlmostImpactDetection(isImpact, impactAlmostRange);
	snakein.PrepareForUse();
	start_s=rsvs3d::TimeStamp("Impact: ", start_f);
	// ======================
	// Spawn
	SAFE_ALGO_TestConn(snakein);
	snakein.SnaxImpactDetection(isImpact);
	snakein.snaxs.SetMaxIndex();
	// int maxIndPreSpawn = snakein.snaxs.GetMaxIndex();
	SpawnArrivedSnaxels(snakein,isImpact);
	start_s=rsvs3d::TimeStamp("Spawn: ", start_s);
	// snakein.PrepareForUse();
	// ==============
	// step away from edge
	// snakein.TakeSpawnStep(maxIndPreSpawn, stepLength);
	// snakein.UpdateCoord();
	// snakein.PrepareForUse();
	// ======================
	// Impact
	SAFE_ALGO_TestConn(snakein);
	snakein.SnaxImpactDetection(isImpact);
	snakein.PrepareForUse();
	start_s=rsvs3d::TimeStamp("Impact: ", start_s);
	// ======================
	// Merge
	MergeAllContactVertices(snakein, isImpact);
	snakein.PrepareForUse();
	start_s=rsvs3d::TimeStamp("Merge: ", start_s);
	// ======================
	// Clean
	CleanupSnakeConnec(snakein);
	// Take Spawn step
	snakein.PrepareForUse();
	SAFE_ALGO_TestConn(snakein);
	snakein.OrientFaces();
	start_s=rsvs3d::TimeStamp("Clean: ", start_s);
	rsvs3d::TimeStamp(" - Connec Update: ", start_f);
}
void SnakeSpawnStep(snake &snakein, int maxIndPreSpawn, double stepLength){
	// double stepLength = 1e-5;

	snakein.TakeSpawnStep(maxIndPreSpawn, stepLength);
	snakein.UpdateCoord();
	snakein.PrepareForUse();

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
	start_s=rsvs3d::TimeStamp("Impact: ", start_f);
	// ======================
	// Spawn
	SAFE_ALGO_TestConn(snakein);
	snakein.SnaxImpactDetection(isImpact);
	SpawnArrivedSnaxels(snakein,isImpact);
	start_s=rsvs3d::TimeStamp("Spawn: ", start_s);
	// ======================
	// Impact
	SAFE_ALGO_TestConn(snakein);
	snakein.SnaxImpactDetection(isImpact);
	snakein.PrepareForUse();
	start_s=rsvs3d::TimeStamp("Impact: ", start_s);
	// ======================
	// Merge
	MergeAllContactVertices(snakein, isImpact);
	snakein.PrepareForUse();
	start_s=rsvs3d::TimeStamp("Merge: ", start_s);
	// ======================
	// Clean
	CleanupSnakeConnec(snakein);
	snakein.PrepareForUse();
	SAFE_ALGO_TestConn(snakein);
	snakein.OrientFaces();
	start_s=rsvs3d::TimeStamp("Clean: ", start_s);

	rsvs3d::TimeStamp(" - Connec Update: ", start_f);
}

double SnakePositionUpdate(snake &rsvsSnake, std::vector<double> &dt,
	double snaxtimestep, double snaxdiststep){

	rsvsSnake.CalculateTimeStep(dt, snaxtimestep, snaxdiststep);
	rsvsSnake.UpdateDistance(dt,snaxdiststep);
	
	rsvsSnake.UpdateCoord();
	rsvsSnake.PrepareForUse();
	return *max_element(dt.begin(), dt.end());
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

	integrate::prepare::Mesh(RSVSobj.paramconf.grid, 
		RSVSobj.paramconf.files.ioin, RSVSobj.snakeMesh, RSVSobj.voluMesh);
	integrate::prepare::Snake(RSVSobj.paramconf.snak, RSVSobj.paramconf.rsvs,
		RSVSobj.paramconf.files.ioin,
		RSVSobj.snakeMesh, RSVSobj.voluMesh, RSVSobj.rsvsSnake);
	integrate::prepare::Triangulation(RSVSobj.snakeMesh, RSVSobj.rsvsSnake, RSVSobj.rsvsTri);
	integrate::prepare::Output(RSVSobj.paramconf, origconf, RSVSobj.outSnake,
		RSVSobj.outgradientsnake, RSVSobj.logFile, RSVSobj.coutFile,
		RSVSobj.cerrFile);
}

void integrate::prepare::Mesh(
	const param::grid &gridconf,
	const param::ioin &ioinconf,
	mesh &snakeMesh,
	mesh &voluMesh
	){
	/*prepares the snake and volume meshes gor the RSVS process*/
	// Local declaration
	
	RSVScalc calcVolus, calcSnakVolu;

	if(gridconf.activegrid.compare("voxel")==0){
		integrate::prepare::grid::Voxel(gridconf, snakeMesh, voluMesh);
	} else if(gridconf.activegrid.compare("voronoi")==0) {
		integrate::prepare::grid::Voronoi(gridconf, snakeMesh, voluMesh);
	} else if(gridconf.activegrid.compare("load")==0) {
		integrate::prepare::grid::Load(ioinconf, snakeMesh, voluMesh);
	} else {
		RSVS3D_ERROR_ARGUMENT(
			(gridconf.activegrid + " not recognised. "
			"Invalid parameter activegrid passed to "
			"RSVS process.").c_str());
	}

	snakeMesh.PrepareForUse();
	snakeMesh.OrientFaces();
	calcSnakVolu.CalculateMesh(snakeMesh);
	calcSnakVolu.ReturnConstrToMesh(snakeMesh,&volu::volume);

	voluMesh.PrepareForUse();
	voluMesh.OrientFaces();
	calcVolus.CalculateMesh(voluMesh);
	calcVolus.ReturnConstrToMesh(voluMesh,&volu::volume);

	std::cout << "Meshes prepared..." << std::endl;
}

void integrate::prepare::grid::Voxel(
	const param::grid &gridconf,
	mesh &snakeMesh,
	mesh &voluMesh
	){

	int ii;
	auto gridSize = gridconf.voxel.gridsizebackground;
	vector<int> elmMapping, backgroundGrid;

	backgroundGrid.reserve(int(gridSize.size()));
	for (int i = 0; i < int(gridSize.size()); ++i)
	{
		backgroundGrid.push_back(gridSize[i]);
		gridSize[i] = gridSize[i] * gridconf.voxel.gridsizesnake[i];
	}

	// Initial build of the grid
	BuildBlockGrid(gridSize, snakeMesh);
	snakeMesh.Scale(gridconf.domain);
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
		RSVS3D_ERROR_ARGUMENT("Incorrect dimension");
	}
	CartesianMapping(snakeMesh,  elmMapping, backgroundGrid);
	CoarsenMesh(snakeMesh,voluMesh,elmMapping);
	snakeMesh.AddParent(&voluMesh,elmMapping);
}

void integrate::prepare::grid::Voronoi(
	const param::grid &gridconf,
	mesh &snakeMesh,
	mesh &voluMesh
	){
	// Vector points are already loaded
	tetgen::apiparam inparam;
	inparam.edgelengths={gridconf.voronoi.snakecoarseness,
		double(gridconf.voronoi.vorosnakelayers)
	};
	inparam.distanceTol = gridconf.voronoi.distancebox;
	for (int i = 0; i < 3; ++i)
	{
		inparam.lowerB[i] = gridconf.domain[i][0];
		inparam.upperB[i] = gridconf.domain[i][1];
	}
	tetgen::RSVSVoronoiMesh(
		gridconf.voronoi.inputpoints,
		voluMesh, snakeMesh, inparam);
}

void integrate::prepare::grid::Load(
	const param::ioin &ioinconf,
	mesh &snakeMesh,
	mesh &voluMesh)
{
	snakeMesh.read(ioinconf.snakemeshname.c_str());
	voluMesh.read(ioinconf.volumeshname.c_str());
	snakeMesh.PrepareForUse(true);
	voluMesh.PrepareForUse(true);

	// need to find mesh lineage
	auto snakMeshVolPts = CoordInVolume(snakeMesh);
	auto elmMapping = voluMesh.VertexInVolume(snakMeshVolPts,3);

	snakeMesh.AddParent(&voluMesh, elmMapping);

	snakeMesh.PrepareForUse();
	voluMesh.PrepareForUse();
}

void integrate::prepare::Snake(
	const param::snaking &snakconf, 
	const param::rsvs &rsvsconf, 
	const param::ioin &ioinconf, 
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
	if(!ioinconf.snakefile.empty()){
		rsvsSnake.read(ioinconf.snakefile.c_str());
		rsvsSnake.snakemesh = &snakeMesh;
		rsvsSnake.PrepareForUse(true);
		rsvsSnake.OrientFaces();
		rsvsSnake.AssignInternalVerts();
	} else {
		SpawnRSVS(rsvsSnake,snakconf.initboundary);
	}
	rsvsSnake.PrepareForUse(true);
	rsvsSnake.OrientFaces();
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
	tecplotfile &outgradientsnake,
	std::ofstream &logFile,
	std::ofstream &coutFile,
	std::ofstream &cerrFile
	){
	std::string outSnakeName;
	
	if (paramconf.files.ioout.outdir.size()!=0){
		filesystem::create_directories(paramconf.files.ioout.outdir);
	}

	outSnakeName = utils::OutputFileName(paramconf, constants::tecplotsnake,".plt");
	outSnake.OpenFile(outSnakeName.c_str());

	if(paramconf.files.ioout.logginglvl>3){
		outSnakeName = utils::OutputFileName(paramconf, constants::tecplotgradient,".plt");
		outgradientsnake.OpenFile(outSnakeName.c_str());
	}

	outSnakeName = utils::OutputFileName(paramconf, "config_call_",".json");
	param::io::writeflat(outSnakeName, origconf);

	outSnakeName = utils::OutputFileName(paramconf, "config_active_",".json");
	param::io::writeflat(outSnakeName, paramconf);
	// Open a text log file

	outSnakeName = utils::OutputFileName(paramconf, "convergence_",".log");
	logFile.open(outSnakeName);
	logFile.precision(16);
	logFile << std::scientific;

	if (paramconf.files.ioout.redirectcout){
		outSnakeName = utils::OutputFileName(paramconf, "cout_",".txt");
		coutFile.open(outSnakeName);
		std::cout.rdbuf(coutFile.rdbuf());
	}
	if (paramconf.files.ioout.redirectcerr){
		outSnakeName = utils::OutputFileName(paramconf, "cerr_",".txt");
		cerrFile.open(outSnakeName);
		std::cerr.rdbuf(cerrFile.rdbuf());
	}
	integrate::utils::SpecialiseTemplateFiles(paramconf);
}


// ====================
// integrate
// 		execute
// ====================

void integrate::execute::All(integrate::RSVSclass &RSVSobj){

	auto coutbuff=std::cout.rdbuf();
	auto cerrbuff=std::cerr.rdbuf();
	std::cout << "Start RSVS preparation" << std::endl;
	integrate::Prepare(RSVSobj);

	std::cout << "Preparation finished - start iteration" << std::endl;
	auto iterateInfo=integrate::execute::RSVSiterate(RSVSobj);

	std::cout << "Iteration finished - start PostProcessing" << std::endl;
	integrate::execute::PostProcessing(RSVSobj,
		iterateInfo.timeT, iterateInfo.nVoluZone, iterateInfo.stepNum);

	std::cout << "PostProcessing finished - start Exporting" << std::endl;
	integrate::execute::Exporting(RSVSobj);

	std::cout << "Exporting finished - close." << std::endl;
	std::cout.rdbuf(coutbuff);
	std::cerr.rdbuf(cerrbuff);
	// std::cout << std::endl <<  " cout Buffer restored" << std::endl;
	// std::cerr << " cerr Buffer restored" << std::endl;
}

integrate::iteratereturns integrate::execute::RSVSiterate(
	integrate::RSVSclass &RSVSobj){

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
		// RSVSobj.calcObj.limLag=10000.0;
		std::cout << std::endl << "Step " << stepNum << " ";
		RSVSobj.calcObj.CalculateTriangulation(RSVSobj.rsvsTri);
		start_s=rsvs3d::TimeStamp(" deriv:", start_s);
		RSVSobj.calcObj.CheckAndCompute(
			RSVSobj.paramconf.rsvs.solveralgorithm);
		// Second cycle
		// start_s=rsvs3d::TimeStamp(" solve:", start_s);
		// RSVSobj.calcObj.CalculateTriangulation(RSVSobj.rsvsTri);
		// start_s=rsvs3d::TimeStamp(" deriv:", start_s);
		// RSVSobj.calcObj.CheckAndCompute(
		// 	RSVSobj.paramconf.rsvs.solveralgorithm);
		// End of Second cycle
		RSVSobj.calcObj.ReturnConstrToMesh(RSVSobj.rsvsTri);
		RSVSobj.calcObj.ReturnVelocities(RSVSobj.rsvsTri);
		start_s=rsvs3d::TimeStamp(" solve:", start_s);
		
		CalculateNoNanSnakeVel(RSVSobj.rsvsSnake);

		integrate::execute::Logging(RSVSobj,
			totT, nVoluZone, stepNum);

		totT += SnakePositionUpdate(RSVSobj.rsvsSnake, dt,
			RSVSobj.paramconf.snak.snaxtimestep,
			RSVSobj.paramconf.snak.snaxdiststep);

		int maxIndPreSpawn = RSVSobj.rsvsSnake.snaxs.GetMaxIndex();
		SnakeConnectivityUpdate(RSVSobj.rsvsSnake, isImpact,
			RSVSobj.paramconf.snak.multiarrivaltolerance);
		SnakeSpawnStep(RSVSobj.rsvsSnake, maxIndPreSpawn,
			RSVSobj.paramconf.snak.spawnposition);
		start_s=clock();
		MaintainTriangulateSnake(RSVSobj.rsvsTri);
		start_s=rsvs3d::TimeStamp(" triangulate:", start_s);
	}
	std::cout << std::endl << "RSVS iteration finished" << std::endl;
	integrate::iteratereturns retStruct(nVoluZone,stepNum, totT);
	return(retStruct);
}

void integrate::execute::Logging(integrate::RSVSclass &RSVSobj,
	double totT, int nVoluZone, int stepNum){
	// Simple function which directs to the correct output
	int logLvl = RSVSobj.paramconf.files.ioout.logginglvl;
	if (0 < logLvl){
		// calcObj logging outputs different amounts of data
		// depending.
		RSVSobj.logFile << "> step" << stepNum << " :," ;
		RSVSobj.logFile << totT << endl;
		integrate::execute::logging::Log(
			RSVSobj.logFile, RSVSobj.calcObj,
			RSVSobj.paramconf.files.ioout.logginglvl
			);
	}

	if (logLvl==2 || logLvl==5){
		integrate::execute::logging::Snake(
			RSVSobj.outSnake, RSVSobj.rsvsSnake,
			RSVSobj.voluMesh, totT, nVoluZone);
	} else if (2 < logLvl && logLvl!=5){
		integrate::execute::logging::FullTecplot(
			RSVSobj.outSnake, RSVSobj.rsvsSnake,
			RSVSobj.rsvsTri, RSVSobj.voluMesh,
			totT, nVoluZone, stepNum);
	}
	if(3 < logLvl){
		integrate::execute::logging::Gradients(RSVSobj.calcObj,
			RSVSobj.rsvsTri, RSVSobj.outgradientsnake, totT);
	}
}

void integrate::execute::PostProcessing(integrate::RSVSclass &RSVSobj,
	double totT, int nVoluZone, int stepNum){

	
	int logLvl = RSVSobj.paramconf.files.ioout.outputlvl;
	if (RSVSobj.paramconf.files.ioout.logginglvl!=5){
		logLvl = max(RSVSobj.paramconf.files.ioout.outputlvl,
			RSVSobj.paramconf.files.ioout.logginglvl);
	}

	if (logLvl>1){
		integrate::execute::postprocess::Snake(
			RSVSobj.rsvsSnake,
			RSVSobj.voluMesh,
			RSVSobj.paramconf);
	}

	if (0 < logLvl){
		RSVSobj.logFile << "> final step" << stepNum << " :," ;
			RSVSobj.logFile << totT << endl;
		integrate::execute::postprocess::Log(
			RSVSobj.logFile, RSVSobj.calcObj,logLvl
			);
	}
	if (logLvl==2 || logLvl==5){
		integrate::execute::logging::Snake(
			RSVSobj.outSnake, RSVSobj.rsvsSnake,
			RSVSobj.voluMesh, totT, nVoluZone);
	} else if (2 < logLvl && logLvl!=5){
		integrate::execute::postprocess::FullTecplot(
			RSVSobj.outSnake, RSVSobj.rsvsSnake,
			RSVSobj.rsvsTri, RSVSobj.voluMesh,
			totT, nVoluZone, stepNum);
	}
	if(3 < logLvl){
		integrate::execute::postprocess::Gradients(RSVSobj.calcObj,
			RSVSobj.rsvsTri, RSVSobj.outgradientsnake, totT);
	}
}

void integrate::execute::Exporting(integrate::RSVSclass &RSVSobj){

	RSVSobj.rsvsSnake.Scale(RSVSobj.paramconf.grid.physdomain);
	auto &paramconf = RSVSobj.paramconf;
	auto outSnakeName =  utils::OutputFileName(paramconf, "rsvs3D_export_",".plt");

	tecplotfile outSnake;
	outSnake.OpenFile(outSnakeName.c_str());

	outSnake.PrintSnake(RSVSobj.rsvsSnake);
	outSnake.PrintMesh(*RSVSobj.rsvsSnake.snakemesh);
	outSnake.PrintMesh(RSVSobj.voluMesh);

	for (auto exportType : RSVSobj.paramconf.files.exportconfig){

		if(exportType.first.compare("su2")==0)
		{
			integrate::execute::exporting::SU2(
				exportType.second, RSVSobj.rsvsSnake, RSVSobj.paramconf
				);
		} else if(exportType.first.compare("")==0 
			&& exportType.second.compare("")==0) 
		{
			/*Do nothing, this is the default*/
		}  else if(exportType.first.compare("snake")==0) 
		{
			RSVSobj.rsvsSnake.Scale(RSVSobj.paramconf.grid.domain);
			auto fileToOpen= utils::OutputFileName(paramconf, "SnakeConnExport_",".msh");
			auto tecFileToOpen = utils::OutputFileName(paramconf, "SnakeSensitivity_",".plt");

			RSVSobj.rsvsSnake.snakeconn.write(fileToOpen.c_str());
			tecplotfile tecsens;
			tecsens.OpenFile(tecFileToOpen.c_str());
			RSVSobj.calcObj.CalculateTriangulation(RSVSobj.rsvsTri);
			RSVSobj.calcObj.CheckAndCompute(
				RSVSobj.paramconf.rsvs.solveralgorithm, true);
			RSVSobj.rsvsSnake.Scale(RSVSobj.paramconf.grid.physdomain);
			tecsens.PrintSnakeSensitivityTime(RSVSobj.rsvsTri, RSVSobj.calcObj);
			
		} else {
			RSVS3D_ERROR_NOTHROW((
				std::string("The export argument '") + exportType.first
				+ "' is not known. Check <input file>.json for parameter : "
				"/files/exportconfig/..."
				).c_str());
		}

	}

	RSVSobj.rsvsSnake.Scale(RSVSobj.paramconf.grid.domain);
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
	rsvsSnake.snakeconn.PrepareForUse();
	rsvsSnake.snakeconn.OrientFaces();
	outSnake.PrintVolumeDat(voluMesh,nVoluZone,1,totT);
	outSnake.PrintSnake(rsvsSnake, 2, totT,
		rsvs3d::constants::tecplot::polygon);
}

void integrate::execute::logging::FullTecplot(
	tecplotfile &outSnake, snake &rsvsSnake,
	triangulation &rsvsTri, mesh &voluMesh,
	double totT, int nVoluZone, int stepNum){
	vector<int> vertList;
	int jj;
	if(rsvsSnake.snaxs.size()>0){

		rsvsSnake.snakeconn.PrepareForUse();
		outSnake.PrintSnake(rsvsSnake, 1, totT, 
			rsvs3d::constants::tecplot::polygon);
		outSnake.PrintTriangulation(rsvsTri,&triangulation::dynatri,2,totT, 
			rsvs3d::constants::tecplot::polygon);
		outSnake.PrintTriangulation(rsvsTri,&triangulation::intertri,3,totT,
			rsvs3d::constants::tecplot::line);
		if (stepNum==0 && rsvsTri.intertri.size()==0){
			outSnake.PrintTriangulation(rsvsTri,&triangulation::dynatri,3,totT,
				rsvs3d::constants::tecplot::line, {rsvsTri.dynatri(0)->index});
		}
		outSnake.PrintTriangulation(rsvsTri,&triangulation::trisurf,4,totT,
			rsvs3d::constants::tecplot::line);
		if (stepNum==0 && rsvsTri.trisurf.size()==0){
			outSnake.PrintTriangulation(rsvsTri,&triangulation::dynatri,4,totT,
				rsvs3d::constants::tecplot::line, {rsvsTri.dynatri(0)->index});
		}
		if (int(rsvsTri.acttri.size())>0){
			outSnake.PrintTriangulation(rsvsTri,&triangulation::stattri,5,totT,
				rsvs3d::constants::tecplot::line,rsvsTri.acttri);
		} else if (stepNum==0){
			outSnake.PrintTriangulation(rsvsTri,&triangulation::dynatri,5,totT,
				rsvs3d::constants::tecplot::line, {rsvsTri.dynatri(0)->index});
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
		outSnake.PrintMesh(*(rsvsSnake.snakemesh),6,totT,
			rsvs3d::constants::tecplot::point,vertList);
		outSnake.PrintVolumeDat(voluMesh,nVoluZone,7,totT);
		outSnake.PrintSnake(rsvsSnake, 8, totT);

	}
}

void integrate::execute::logging::Gradients(const RSVScalc &calcObj,
	const triangulation &rsvsTri,	tecplotfile &outgradientsnake, double totT){

	outgradientsnake.PrintSnakeGradients(rsvsTri, calcObj, 1, totT);
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
	snake &rsvsSnake, mesh &voluMesh, param::parameters &paramconf){

	std::string fileToOpen;
	
	fileToOpen= utils::OutputFileName(paramconf, "VoluMesh_",".msh");
	voluMesh.write(fileToOpen.c_str());
	paramconf.files.ioin.volumeshname = fileToOpen;

	fileToOpen= utils::OutputFileName(paramconf, "SnakeMesh_",".msh");
	rsvsSnake.snakemesh->write(fileToOpen.c_str());
	paramconf.files.ioin.snakemeshname = fileToOpen;

	fileToOpen= utils::OutputFileName(paramconf, "SnakeConn_",".msh");
	rsvsSnake.snakeconn.write(fileToOpen.c_str());

	fileToOpen= utils::OutputFileName(paramconf, "Snake_",".3snk");
	rsvsSnake.write(fileToOpen.c_str());
	paramconf.files.ioin.snakefile = fileToOpen;

	paramconf.files.ioin.casename = "restart_" + paramconf.files.ioout.pattern;
	fileToOpen= utils::OutputFileName(paramconf, "config_restart_",".json");
	param::io::writeflat(fileToOpen, paramconf);
	
}

void integrate::execute::postprocess::FullTecplot(
	tecplotfile &outSnake, snake &rsvsSnake,
	triangulation &rsvsTri, mesh &voluMesh,
	double totT, int nVoluZone, int stepNum){

	integrate::execute::logging::FullTecplot(
		outSnake, rsvsSnake, rsvsTri, voluMesh,
		totT, nVoluZone, stepNum);
}

void integrate::execute::postprocess::Gradients(const RSVScalc &calcObj,
	const triangulation &rsvsTri,	tecplotfile &outgradientsnake, double totT){

	integrate::execute::logging::Gradients(calcObj,rsvsTri,
		outgradientsnake, totT);
}
// ====================
// integrate
// 		execute
// 			exporting
// ====================

void integrate::execute::exporting::SU2(std::string su2ConfStr,
	snake &rsvsSnake, param::parameters &paramconf){

	std::string su2Path = "", patternReplace = "<pattern>";
	std::string su2FileStr = su2ConfStr.substr(0, su2ConfStr.find(","));
	// Make correct file path relative to the output directory 
	su2Path += paramconf.files.ioout.outdir + "/";

	auto patternIt = su2FileStr.find(patternReplace);
	if(patternIt < su2FileStr.size()){
		su2Path += su2FileStr.substr(0, patternIt) 
			+ paramconf.files.ioout.pattern
			+ su2FileStr.substr(patternIt+patternReplace.size());
	} else {
		su2Path += su2FileStr;
	}

	const filesystem::path pathout=su2Path;
	filesystem::create_directories(pathout.parent_path());

	tetgen::apiparam su2MeshParam(su2ConfStr.substr(su2ConfStr.find(",")+1));
	tetgen::SnakeToSU2(rsvsSnake, su2Path, su2MeshParam);
}


// ====================
// 
//  Utilities
// 
// ====================

void integrate::utils::WriteModifiedTemplate(const std::string &fileIn, 
	const std::string &fileOut,	const std::string &oldLine,
	const std::string newLine){

	ifstream in(fileIn, ios::in);
	if(in.fail()){
		RSVS3D_ERROR_NOTHROW((fileIn + " was not found, template could not be"
			" specialised.").c_str());
		return;
	}
	ofstream out(fileOut, ios::out);
	std::string tempString;
	while (!in.eof()){
		getline(in,tempString);
		auto pos = tempString.find(oldLine);
		if(pos!=std::string::npos){
			out << tempString.substr(0, pos) << newLine << std::endl;
			break;
		} else {
			out << tempString << std::endl;
		}
	}

	out << in.rdbuf();

}

void integrate::utils::SpecialiseTemplateFiles(const param::parameters &paramconf){

	// param::tecplottemplate& 
	auto& tecconfig = paramconf.files.ioout.tecplot;
	int loglvl = paramconf.files.ioout.logginglvl;

	if (loglvl<=3){
		SpecialiseTemplateFile(tecconfig, loglvl, paramconf.files.ioout, 
			constants::tecplotsnake);
	} else if (loglvl!=5){
		SpecialiseTemplateFile(tecconfig, 3, paramconf.files.ioout, 
			constants::tecplotsnake);
	}
	if(loglvl>=4){
		SpecialiseTemplateFile(tecconfig, 5, paramconf.files.ioout, 
			constants::tecplotgradient);
	}
	if (loglvl>5){
		SpecialiseTemplateFile(tecconfig, loglvl, paramconf.files.ioout, 
			constants::tecplotsnake);
	}

}

void integrate::utils::SpecialiseTemplateFile(const param::tecplottemplate &tecconfig, 
	int logLvl,
	const param::ioout &ioout, std::string fileName){

	// Build the name of the correct file from logging lvl and 
	std::string templateLayout, outputLayout; 
	templateLayout = tecconfig.TemplateLogging(logLvl);
	outputLayout = ioout.outdir + "/";
	outputLayout += tecconfig.TemplateLogging(logLvl,
		false,std::string("_") + ioout.pattern);
	std::string lineOut = tecconfig.filenameregex +"'\""+
		integrate::utils::OutputFileName("", ioout.pattern, 
			fileName, ".plt")+"\"'";
	integrate::utils::WriteModifiedTemplate(templateLayout, outputLayout,
		tecconfig.filenameregex, lineOut);
	

}

std::string integrate::utils::OutputFileName(const param::parameters &paramconf, 
	std::string fileName, std::string extension){
	return OutputFileName(paramconf.files.ioout.outdir,
		paramconf.files.ioout.pattern, fileName, extension);
}

std::string integrate::utils::OutputFileName(const std::string rootDirectory, 
	const std::string &filePattern,	std::string fileName, 
	std::string extension){

	if (rootDirectory.compare("")==0){
		return fileName + filePattern + extension;
	}
	return rootDirectory + "/" + fileName + filePattern + extension;
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
	tecplotfile outgradientSnake;
	std::ofstream logFile;
	std::ofstream coutFile;
	std::ofstream cerrFile;

	origconf = paramconf;
	paramconf.PrepareForUse();
	try {
		integrate::prepare::Mesh(paramconf.grid, paramconf.files.ioin,
			snakeMesh, voluMesh);
		voluMesh.volus.disp();
	} catch (exception const& ex) { 
		cerr << "integrate::prepare::Mesh(paramconf.grid, snakeMesh, "
			"voluMesh);" << endl;
		cerr << "Exception: " << ex.what() <<endl; 
		return -1;
	} 
	
	try {
		integrate::prepare::Snake(paramconf.snak, paramconf.rsvs, 
			paramconf.files.ioin,
			snakeMesh,voluMesh,  rsvsSnake);
	} catch (exception const& ex) { 
		cerr << "integrate::prepare::Snake(paramconf.snak, snakeMesh,"
			" rsvsSnake);" << endl;
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
		integrate::prepare::Output(paramconf, origconf, outSnake,
			outgradientSnake,logFile,
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

