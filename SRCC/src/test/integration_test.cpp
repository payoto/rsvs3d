#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <tuple>

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

typedef std::tuple<double, Eigen::MatrixXd, double, double> processed_derivative;

processed_derivative test_RSVScalcvsFD(
	integrate::RSVSclass &RSVSobj, 
	bool SetUseSurfCentreDeriv,
	int testVert, double fdStep,
	int constrFocus){

	std::tuple<double, Eigen::MatrixXd, double> testtuple;

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
	Eigen::VectorXd diffDConstr = (analyticalDConstr.block(0,dvNum,RSVSobj.calcObj.numConstr(),1)
		-fdDconstr).array()
		/ analyticalDConstr.block(0,dvNum,RSVSobj.calcObj.numConstr(),1).array();

	int nRows = diffDConstr.size();
	for (int i = 0; i < nRows; ++i)
	{
		if (IsAproxEqual(0.0, analyticalDConstr(i,dvNum))){
			diffDConstr[i] = 0.0;
		}
	}
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

	return {diffDObj, diffDConstr, analyticalDObj[dvNum], 
		analyticalDConstr(constrFocus, dvNum)}; 
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

void StreamProcessedDerivative(processed_derivative &in, double fdStep,
	std::string modifier, int constrIntrest){
	std::cout << fdStep << " (";
	cout.unsetf(std::ios::fixed | std::ios::scientific);
	std::cout << std::setw(3) << floor(log10(DBL_EPSILON/fdStep));
	std::cout << std::scientific;
	std::cout << ")" << " - " 
		<< std::setw(10) << modifier << " centroid - obj ("
		<< std::setw(std::cout.precision()+7) << std::get<2>(in) 
		<< "): " << std::setw(std::cout.precision()+7) 
		<< std::get<0>(in) << " - constr ("
		<< std::setw(std::cout.precision()+7) << std::get<3>(in) 
		<< "): ";
	std::cout << std::setw(std::cout.precision()+7) 
		<< std::get<1>(in)(constrIntrest,0) << " | stats : ";
	StreamStatistics(std::get<1>(in), std::cout);
}

std::stringstream Test_CompareSnaxDeriv_FDANA(integrate::RSVSclass &RSVSobj, 
	int testVert, int constrFocus=0){

	std::vector<processed_derivative> derivDiffNoCentre, derivDiffCentre;
	std::vector<double> fdSteps = {1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10}; 

	std::stringstream strstream;
	auto coutbuff=std::cout.rdbuf();
	auto cerrbuff=std::cerr.rdbuf();

	std::cout.rdbuf(strstream.rdbuf());

	for (double fdStep : fdSteps){
		derivDiffNoCentre.push_back(
			test_RSVScalcvsFD(RSVSobj,false, testVert, fdStep, constrFocus));
		derivDiffCentre.push_back(
			test_RSVScalcvsFD(RSVSobj, true, testVert, fdStep, constrFocus));
	}
	std::cout.rdbuf(coutbuff);
	std::cerr.rdbuf(cerrbuff);


	std::cout << "--------------------------------------------------" << std::endl
		<< "Summary of results: " << std::endl;
	std::cout << "rsvsmath_automatic_eps_surf : " 
		<< rsvsmath_automatic_eps_surf << std::endl;
	std::cout << std::endl;
	std::cout << std::scientific;
	RSVSobj.rsvsSnake.snaxs.isearch(testVert)->disp();
	RSVSobj.rsvsSnake.snakeconn.verts.isearch(testVert)->disp();
	for (int i = 0; i < int(fdSteps.size()); ++i)
	{

		StreamProcessedDerivative(derivDiffNoCentre[i], fdSteps[i],"without",
			constrFocus);
		StreamProcessedDerivative(derivDiffCentre[i], fdSteps[i], "with",
			constrFocus);

	}
	cout.unsetf(std::ios::fixed | std::ios::scientific);
	return strstream;
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
	double  rsvsmath_automatic_eps_surfprev = rsvsmath_automatic_eps_surf;
	rsvsmath_automatic_eps_surf = 0.0;
	Test_CompareSnaxDeriv_FDANA(RSVSobj, testVert);

	rsvsmath_automatic_eps_surf=rsvsmath_automatic_eps_surfprev;
	return 0;
}


void CompareDerivativesSpikeStepNum(std::string numstr){
	integrate::RSVSclass RSVSobj;
	std::vector<std::string> commands;

	commands.push_back("CompareDerivativesSpike"+numstr);
	commands.push_back("-l");
	commands.push_back("config/522_fbullet.json");
	commands.push_back("-p");
	commands.push_back("/snak/maxsteps:"+numstr);
	commands.push_back("-p");
	commands.push_back("/files/ioout/logginglvl:5");
	RunCommandFromString(RSVSobj, commands);

	/*Perform tests, compare */
	// Choose a snaxel pass all  it through the RSVScalc object
	// double fdStep = 1e-6;
	int testVert = RSVSobj.rsvsSnake.snaxs(543)->index;
	std::stringstream sumarries;
	auto coutbuff = std::cout.rdbuf(sumarries.rdbuf());
	double  rsvsmath_automatic_eps_surfprev = rsvsmath_automatic_eps_surf;
	// rsvsmath_automatic_eps_surf = 0.0;
	auto stream0 = Test_CompareSnaxDeriv_FDANA(RSVSobj, testVert,2);
	rsvsmath_automatic_eps_surf = 0.0;
	auto streamEps =Test_CompareSnaxDeriv_FDANA(RSVSobj, testVert,2);
	rsvsmath_automatic_eps_surf=rsvsmath_automatic_eps_surfprev;

	std::cout.rdbuf(coutbuff);

	// std::cout << stream0.str();
	// std::cout << streamEps.str();
	
	std::cout << sumarries.str();
	std::cout << std::endl << "Step number: " << numstr << "+1" << std::endl;
	std::cout  << "DeltaDv: " 
		<< RSVSobj.calcObj.deltaDV[RSVSobj.calcObj.dvMap.find(testVert)] 
		<< std::endl;

	std::cout << std::endl << "Lagrangian multipliers:" << std::endl;
	std::cout << RSVSobj.calcObj.lagMult.transpose();
	std::cout << std::endl << "-------------------------------------------" 
		<< std::endl;
}

int integrate::test::CompareDerivativesSpike(){
	CompareDerivativesSpikeStepNum("94");
	CompareDerivativesSpikeStepNum("95");
	CompareDerivativesSpikeStepNum("96");
	CompareDerivativesSpikeStepNum("97");
	CompareDerivativesSpikeStepNum("98");
	return 0;
}