#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <ctime>
#include <cstdlib>

#include "snake.hpp"
#include "snakeengine.hpp"
#include "triangulate.hpp"
#include "postprocessing.hpp"
#include "meshrefinement.hpp" 
#include "RSVSmath.hpp"
#include "RSVScalc.hpp"
#include "RSVSalgorithm.hpp"
#include "RSVSintegration.hpp"
#include "warning.hpp"
#include "tetgenrsvs.hpp"
#include "matrixtools.hpp"
#include "main.hpp"
#include "parameters.hpp"
#include "RSVSclass.hpp"
#include "RSVSmath_automatic.hpp"

std::pair<double, Eigen::MatrixXd> test_RSVScalcvsFD(integrate::RSVSclass &RSVSobj, 
	bool SetUseSurfCentreDeriv,
	int testVert, double fdStep){


	auto & sn = RSVSobj.rsvsSnake;
	RSVSobj.calcObj.SetUseSurfCentreDeriv(SetUseSurfCentreDeriv);
	RSVSobj.calcObj.CalculateTriangulation(RSVSobj.rsvsTri);
	int dvNum = RSVSobj.calcObj.dvMap.find(testVert);
	auto analyticalDObj = RSVSobj.calcObj.dObj;
	auto analyticalDConstr = RSVSobj.calcObj.dConstr;
	auto constrPreFD = RSVSobj.calcObj.constr;
	auto objPreFd = RSVSobj.calcObj.obj;
	sn.snaxs[sn.snaxs.find(testVert)].d += fdStep;
	sn.snaxs.PrepareForUse();
	sn.UpdateCoord();

	if(SetUseSurfCentreDeriv){
		RSVSobj.rsvsTri.CalcTriVertPos();
		RSVSobj.rsvsTri.PrepareForUse();
	}

	RSVSobj.calcObj.CalculateTriangulation(RSVSobj.rsvsTri);
	auto constrPostFD = RSVSobj.calcObj.constr;
	auto objPostFd = RSVSobj.calcObj.obj;
	auto fdDObj = (objPostFd-objPreFd)/fdStep;
	auto fdDconstr = (constrPostFD-constrPreFD)/fdStep;

	double diffDObj = (analyticalDObj[dvNum]-fdDObj)/analyticalDObj[dvNum];
	auto diffDConstr = (analyticalDConstr.block(0,dvNum,RSVSobj.calcObj.numConstr(),1)
		-fdDconstr).array()
		/ analyticalDConstr.block(0,dvNum,RSVSobj.calcObj.numConstr(),1).array();
	
	sn.snaxs[sn.snaxs.find(testVert)].d -= fdStep;
	sn.snaxs.PrepareForUse();
	sn.UpdateCoord();
	if(SetUseSurfCentreDeriv){
		RSVSobj.rsvsTri.CalcTriVertPos();
		RSVSobj.rsvsTri.PrepareForUse();
	}
	std::cout << "------------------------------------------" << std::endl
		<< "Analytical gradient vs Finite difference (" << fdStep << ")"<< std::endl;
	std::cout << "SetUseSurfCentreDeriv : " << SetUseSurfCentreDeriv  << std::endl;
	std::cout << "snax->d : " << sn.snaxs.isearch(testVert)->d  << std::endl;
	std::cout << "diffDObj: (ana-fd)/ana = delta " << std::endl
		<< analyticalDObj[dvNum] << "   " << fdDObj 
		<< "  =  " << diffDObj << std::endl << std::endl;
	std::cout << "diffDConstr (ana-fd)/ana = delta: " 
		<< std::endl  << analyticalDConstr.block(0,dvNum,
			RSVSobj.calcObj.numConstr(),1).transpose() << std::endl
		<< std::endl << fdDconstr.transpose() << std::endl
		<< std::endl << diffDConstr.transpose() << std::endl;

	return {diffDObj, diffDConstr}; 
}

void RunCommandFromString(integrate::RSVSclass &RSVSobj,
	std::vector<std::string> &commands){
	std::stringstream strstream;
	param::parameters paramconf;

	auto parseOut = parse::StringParser(commands, paramconf);


	if (parseOut.execFlow<0){
		RSVS3D_ERROR("Argument parsing returned a non-runable structure");
	}
	RSVSobj.paramconf = paramconf;

	/*Running of the core iteration process*/

	auto coutbuff=std::cout.rdbuf();
	auto cerrbuff=std::cerr.rdbuf();
	std::cout << "Start RSVS preparation" << std::endl;
	integrate::Prepare(RSVSobj);

	std::cout << "Preparation finished - start iteration" << std::endl;
	std::cout.rdbuf(strstream.rdbuf());
	auto iterateInfo=integrate::execute::RSVSiterate(RSVSobj);
	std::cout.rdbuf(coutbuff);
	std::cout << "Iteration finished - start PostProcessing" << std::endl;
	integrate::execute::PostProcessing(RSVSobj,
		iterateInfo.timeT, iterateInfo.nVoluZone, iterateInfo.stepNum);

	std::cout << "PostProcessing finished - start Exporting" << std::endl;
	integrate::execute::Exporting(RSVSobj);

	std::cout << "Exporting finished - close." << std::endl;
	std::cout.rdbuf(coutbuff);
	std::cerr.rdbuf(cerrbuff);
}
void Test_CompareSnaxDeriv_FDANA(integrate::RSVSclass &RSVSobj, int testVert){

	std::vector<std::pair<double, Eigen::MatrixXd>> derivDiffNoCentre, derivDiffCentre;
	std::vector<double> fdSteps = {1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10}; 
	for (double fdStep : fdSteps){
		derivDiffNoCentre.push_back(
			test_RSVScalcvsFD(RSVSobj,false, testVert, fdStep));
		derivDiffCentre.push_back(
			test_RSVScalcvsFD(RSVSobj, true, testVert, fdStep));
	}


	std::cout << "--------------------------------------------------" << std::endl
		<< "Sumarry of results: " << std::endl;
	std::cout << std::endl;
	for (int i = 0; i < int(fdSteps.size()); ++i)
	{
		std::cout << fdSteps[i] << " - Not centroid: obj " 
			<< derivDiffNoCentre[i].first << " constr ";
		StreamStatistics(derivDiffNoCentre[i].second, std::cout);
		std::cout << fdSteps[i] << " - Yes centroid: obj " 
			<< derivDiffCentre[i].first << " constr ";
		StreamStatistics(derivDiffCentre[i].second, std::cout);
	}

}

int integrate::test::CompareSurfCentreDerivatives(){
	// Test a snake 
	integrate::RSVSclass RSVSobj;
	std::vector<std::string> commands;

	commands.push_back("Test_CompareSurfCentreDerivatives");
	commands.push_back("-l");
	commands.push_back("config/NURBStest.json");
	commands.push_back("-p");
	commands.push_back("/snak/maxsteps:50");
	RunCommandFromString(RSVSobj, commands);

	/*Perform tests, compare */
	// Choose a snaxel pass all  it through the RSVScalc object
	// double fdStep = 1e-6;
	auto & sc = RSVSobj.rsvsSnake.snakeconn;
	auto & sn = RSVSobj.rsvsSnake;
	int testVert = rsvs3d::constants::__notfound;
	int surfTest = 0;
	int nSurf = sc.surfs.size();
	for (int i = 0; i < nSurf; ++i)
	{
		if(sc.surfs(i)->edgeind.size()>4){
			surfTest = sc.surfs(i)->index;
			auto testVerts = sc.surfs.isearch(surfTest)->vertind(sc);
			for (auto isTestVert : testVerts){
				auto dTest = sn.snaxs.isearch(isTestVert)->d;
				if(dTest>0.10 && dTest< 0.90){
				// if(dTest>fdStep*10 && dTest< (1.0 - fdStep*10)){
					testVert = isTestVert;
					break;
				}
			}
			if(testVert!=rsvs3d::constants::__notfound){
				break;
			}
		}
	}
	if(testVert==rsvs3d::constants::__notfound){
		RSVS3D_ERROR("Index was supposed to be above 0");
	}
	rsvsmath_automatic_eps_surf = 0.0;
	Test_CompareSnaxDeriv_FDANA(RSVSobj, testVert);
	return 0;
}