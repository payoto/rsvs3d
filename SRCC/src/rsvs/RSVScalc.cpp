
#include <cstdlib>
#include <iostream>
#include <fstream>

#include "snake.hpp"
#include "snakevel.hpp"
#include "vectorarray.hpp"
#include "RSVScalc.hpp"
#include "matrixtools.hpp"

using namespace std; 

//silent functions
  

/*
Implementation of the interfaces for RSVScalc
This file's object can be compiled quickly
*/

void RSVScalc::PrepTriangulationCalc(const triangulation &triRSVS){

	int ii;
	int nDv, nConstr;
	vector<int> vecin;
	// prepare the SQP object
	nConstr=triRSVS.meshDep->CountVoluParent(); 
	if(triRSVS.snakeDep!=NULL){
		nDv=triRSVS.snakeDep->snaxs.size();
	} else {
		nDv=0;
	} 
	BuildMathArrays(nDv, nConstr);

	// TODO this needs to be supported by mapping each volume to the constraint position
	// There can be more than one constraint for each cell.
	BuildConstrMap(triRSVS);

	vecin.clear();
	vecin.reserve(nDv);
	for(ii=0; ii<nDv; ++ii){
		vecin.push_back(triRSVS.snakeDep->snaxs(ii)->index);
	}
	BuildDVMap(vecin);
}

void RSVScalc::CalculateTriangulation(const triangulation &triRSVS, int derivMethod){

	int ii,ni;


	this->returnDeriv=true;

	// prepare the SQP object
	this->PrepTriangulationCalc(triRSVS);

	// Calculate the SQP object
	
	auto calcTriFunc=&RSVScalc::CalcTriangle;
	switch(derivMethod){
		case 1:
		calcTriFunc = &RSVScalc::CalcTriangleFD;
		break;
		case 2:
		calcTriFunc = &RSVScalc::CalcTriangleDirectVolume;
		break;
	}
	ni=triRSVS.dynatri.size();
	for(ii = 0; ii< ni ; ii++){
		(this->*calcTriFunc)(*(triRSVS.dynatri(ii)), triRSVS, true, true, true);
	} 
	ni=triRSVS.intertri.size();
	for(ii = 0; ii< ni ; ii++){
		(this->*calcTriFunc)(*(triRSVS.intertri(ii)), triRSVS, false, true, true);
	} 
	// PrintMatrixFile(this->dConstr, "matrix_dConstr.txt");
	ni=triRSVS.acttri.size();
	for(ii = 0; ii< ni ; ii++){
		(this->*calcTriFunc)(*(triRSVS.stattri.isearch(triRSVS.acttri[ii])), triRSVS, 
			false, true, false);
	} 
	// Output some data to check it makes sense
}

void RSVScalc::ReturnVelocities(triangulation &triRSVS){

	int ii, ni;

	ni=triRSVS.snakeDep->snaxs.size();
	for(ii=0; ii<ni; ii++){
		triRSVS.snakeDep->snaxs[ii].v = 
			deltaDV[dvMap.find(triRSVS.snakeDep->snaxs(ii)->index)];
	}
	triRSVS.snakeDep->snaxs.PrepareForUse();
}

void RSVScalc::CalculateMesh(mesh &meshin){

	int ii,ni, jj, nj;
	int nDv, nConstr;
	//vector<int> vecin;
	triangulation triRSVS(meshin);
	// prepare the SQP object
	triRSVS.PrepareForUse();
	this->returnDeriv=false;
	nConstr=meshin.volus.size(); 
	
	nDv=0;
	
	BuildMathArrays(nDv, nConstr);

	// TODO this needs to be supported by mapping each volume to the constraint position
	// There can be more than one constraint for each cell.
	//BuildConstrMap(triRSVS);
	BuildConstrMap(meshin);

	// Calculate the SQP object
	nj=meshin.surfs.size();
	for(jj = 0; jj< nj ; jj++){
		TriangulateSurface(*(meshin.surfs(jj)),meshin,triRSVS.stattri, 
					triRSVS.trivert, 1, 1);
		triRSVS.PrepareForUse();
		triRSVS.CalcTriVertPos();
		triRSVS.PrepareForUse();
		ni=triRSVS.stattri.size();
		for(ii = 0; ii< ni ; ii++){
			CalcTriangle(*(triRSVS.stattri(ii)), triRSVS, true, true, false);
		} 
		triRSVS.stattri.clear();
		triRSVS.trivert.clear();
	}

	// Output some data to check it makes sense
}

void RSVScalc::Print2Screen(int outType)const {
	cout << "Math result: obj " << obj << " false access " << falseaccess << endl;
	cout << "Constr " ;
	for (int i = 0; i < nConstr; ++i)
	{
		cout << constr[i] << " ";
	}
	cout << "constrTarg :" <<  endl ;	
	PrintMatrix(constrTarg);
	cout << endl;
	if((nConstr<10 && nDv < 20 && outType==1) || outType==3){
		cout << "constr :" <<  endl ;	
		PrintMatrix(constr);
		cout << "dObj :" <<  endl ;	
		PrintMatrix(dObj);
		cout << "lagMult :" <<  endl ;	
		PrintMatrix(lagMult);

		cout << "dConstr :" <<  endl ;	
		PrintMatrix(dConstr);
		cout << "HConstr :" <<  endl ;	
		PrintMatrix(HConstr);
		cout << "HObj :" <<  endl ;	
		PrintMatrix(HObj);
		cout << "constrTarg :" <<  endl ;	
		PrintMatrix(constrTarg);
	}
	if (outType==2){
		cout << "constrTarg :" <<  endl ;	
		PrintMatrix(constrTarg);
		cout << endl;
		for(int i=0; i<nDv; ++i){
			cout << deltaDV[i] << " ";
		}
		cout << endl;
		for(int i=0; i<nConstr; ++i){
			cout << lagMult[i] << " ";
		}
		cout << endl;
	}
	if (outType==4){
		const char* file="matrices\\dumpmatout.txt";
		cout << "constr :" <<  endl ;	
		PrintMatrixFile(constr,file);
		cout << "dObj :" <<  endl ;	
		PrintMatrixFile(dObj,file);
		cout << "lagMult :" <<  endl ;	
		PrintMatrixFile(lagMult,file);

		cout << "dConstr :" <<  endl ;	
		PrintMatrixFile(dConstr,file);
		cout << "HConstr :" <<  endl ;	
		PrintMatrixFile(HConstr,file);
		cout << "HObj :" <<  endl ;	
		PrintMatrixFile(HObj,file);
		cout << "constrTarg :" <<  endl ;	
		PrintMatrixFile(constrTarg,file);
	}
}

void RSVScalc::ReturnConstrToMesh(triangulation &triRSVS) const {
	
	int ii, ni;
	vector<double> temp;
	ni=constr.size();
	temp.reserve(ni);

	for(ii=0; ii<ni; ii++){
		temp.push_back(constr[ii]);
	}

	triRSVS.meshDep->MapVolu2Parent(temp, this->constrList, &volu::fill);
}

void RSVScalc::ReturnConstrToMesh(mesh &meshin, double volu::*mp) const {
	
	int ii, ni;
	vector<double> temp;
	ni=constr.size();
	temp.reserve(ni);

	for(ii=0; ii<ni; ii++){
		temp.push_back(constr[ii]);
	}
	meshin.MapVolu2Self(temp, constrMap.vec, mp);
}

void RSVScalc::BuildMathArrays(int nDvIn, int nConstrIn){
	// Builds the target math arrays
	
	nDv=nDvIn;
	nConstr=nConstrIn;
	isConstrAct.clear();
	isDvAct.clear();
	isConstrAct.assign(nConstr,false);
	isDvAct.assign(nDv,false);


	constr.setZero(nConstr);
	dConstr.setZero(nConstr,nDv);
	HConstr.setZero(nDv,nDv); 

	obj=0.0;
	dObj.setZero(nDv);
	HObj.setZero(nDv,nDv); 

	if(nConstr!=lagMult.size()){
		lagMult.setZero(nConstr);
	}
	deltaDV.setZero(nDv);
	// constrTarg.setZero(nConstr);

	dvCallConstr.setZero(nDv,1);
}

void RSVScalc::BuildConstrMap(const triangulation &triangleRSVS){

	// explore parents of mesh adding 1 by 1 elemind
	// for each parent
	// for each snakemesh.volu
	// Assign to constrMap.targ = the position in parent.volu of the parentconn
	vector<double> voluVals;
	triangleRSVS.meshDep->ReturnParentMap(constrMap.vec,constrMap.targ,constrList,voluVals);

	constrMap.GenerateHash();
	constrTarg.setZero(voluVals.size());
	for (int i = 0; i < int(voluVals.size()); ++i)
	{
		constrTarg[i] = voluVals[i];
	}
}

void RSVScalc::BuildConstrMap(const mesh &meshin){
	int ni, ii;
	ni=meshin.volus.size();
	constrMap.vec.clear();
	constrMap.targ.clear();
	constrMap.vec.reserve(ni);
	constrMap.targ.reserve(ni);
	for(ii=0; ii<ni; ++ii){
		constrMap.vec.push_back(meshin.volus(ii)->index);
		constrMap.targ.push_back(ii);
	}


	constrMap.GenerateHash();
}

void RSVScalc::BuildDVMap(const vector<int> &vecin){
	dvMap.vec=vecin;
	sort(dvMap.vec);
	unique(dvMap.vec);
	dvMap.GenerateHash();
}


