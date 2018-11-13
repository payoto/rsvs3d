
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath> 
#include "snake.hpp"
#include "snakevel.hpp"
#include "vectorarray.hpp"
#include "RSVSinterface.hpp"
#include <Eigen>

using namespace std; 
using namespace Eigen; 

//silent functions
template <class T> void PrintMatrixFile(const vector<T> &mat, const char * name);


void SQPcalc::PrepTriangulationCalc(const triangulation &triRSVS){

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

void SQPcalc::CalculateTriangulation(const triangulation &triRSVS, int derivMethod){

	int ii,ni;


	this->returnDeriv=true;

	// prepare the SQP object
	this->PrepTriangulationCalc(triRSVS);

	// Calculate the SQP object
	
	auto calcTriFunc=&SQPcalc::CalcTriangle;
	switch(derivMethod){
		case 1:
		calcTriFunc = &SQPcalc::CalcTriangleFD;
		break;
		case 2:
		calcTriFunc = &SQPcalc::CalcTriangleDirectVolume;
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

void SQPcalc::ReturnVelocities(triangulation &triRSVS){

	int ii, ni;

	ni=triRSVS.snakeDep->snaxs.size();
	for(ii=0; ii<ni; ii++){
		triRSVS.snakeDep->snaxs[ii].v = 
			deltaDV[dvMap.find(triRSVS.snakeDep->snaxs(ii)->index)];
	}
	triRSVS.snakeDep->snaxs.PrepareForUse();


}

void SQPcalc::CalculateMesh(mesh &meshin){

	int ii,ni, jj, nj;
	int nDv, nConstr;
	//vector<int> vecin;
	triangulation triRSVS(meshin);
	// prepare the SQP object
	
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

void SQPcalc::Print2Screen(int outType)const {
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

void SQPcalc::ReturnConstrToMesh(triangulation &triRSVS) const {
	
	int ii, ni;
	vector<double> temp;
	ni=constr.size();
	temp.reserve(ni);

	for(ii=0; ii<ni; ii++){
		temp.push_back(constr[ii]);
	}

	triRSVS.meshDep->MapVolu2Parent(temp, this->constrList, &volu::fill);


}
void SQPcalc::ReturnConstrToMesh(mesh &meshin, double volu::*mp) const {
	
	int ii, ni;
	vector<double> temp;
	ni=constr.size();
	temp.reserve(ni);

	for(ii=0; ii<ni; ii++){
		temp.push_back(constr[ii]);
	}
	meshin.MapVolu2Self(temp, constrMap.vec, mp);
}

void SQPcalc::BuildMathArrays(int nDvIn, int nConstrIn){
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


void SQPcalc::BuildConstrMap(const triangulation &triangleRSVS){

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

void SQPcalc::BuildConstrMap(const mesh &meshin){
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

void SQPcalc::BuildDVMap(const vector<int> &vecin){
	dvMap.vec=vecin;
	sort(dvMap.vec);
	unique(dvMap.vec);
	dvMap.GenerateHash();
}

void SQPcalc::CalcTriangle(const triangle& triIn, const triangulation &triRSVS,
	bool isObj, bool isConstr, bool isDeriv){


	int ii,ni,jj,nj,kk,ll,nCellTarg;
	int isCentre,subTemp,subTemp1,subTemp2,subTemp3,nDvAct;
	SurfCentroid centreCalc;
	Volume VolumeCalc;
	Area AreaCalc;
	vector<int> dvList,subTempVec;
	HashedVector<int,int> dvListMap;
	vector<vector<double> const *> veccoord;
	MatrixXd HPos, dPos;
	MatrixXd HVal, dVal;
	MatrixXd dConstrPart,HConstrPart, HObjPart;
	MatrixXd dObjPart;
	double constrPart, objPart;
	ArrayVec<double>* HValpnt=NULL, *dValpnt=NULL;
	double * retVal;


	veccoord.reserve(3);
	isCentre=0;
	ni=3;
	for(ii=0; ii<ni; ++ii){
		if(triIn.pointtype[ii]==1){
			veccoord.push_back(&(triRSVS.meshDep->verts.isearch(triIn.pointind[ii])->coord));
		} else if (triIn.pointtype[ii]==2){
			veccoord.push_back(&(triRSVS.snakeDep->snakeconn.verts.isearch(triIn.pointind[ii])->coord));
		} else if (triIn.pointtype[ii]==3){
			veccoord.push_back((triRSVS.trivert.isearch(triIn.pointind[ii])->coord.retPtr()));
			isCentre++;
		}
	}

	// Constr and objective
	VolumeCalc.assign(veccoord[0],veccoord[1],veccoord[2]);
	AreaCalc.assign(veccoord[0],veccoord[1],veccoord[2]);

	// Active DV lists

	for(ii=0; ii<ni; ++ii){
		if (triIn.pointtype[ii]==2){
			dvList.push_back(triIn.pointind[ii]);
		} else if (triIn.pointtype[ii]==3 && false){
			subTemp=triRSVS.trivert.find(triIn.pointind[ii]);
			nj=triRSVS.trisurf(subTemp)->indvert.size();
			for(jj=0;jj<nj;++jj){
				if (triRSVS.trisurf(subTemp)->typevert[jj]==2){
					dvList.push_back(triRSVS.trisurf(subTemp)->indvert[jj]);
				}
			}
		}
	}
	dvListMap.vec=dvList;
	sort(dvListMap.vec);
	unique(dvListMap.vec);
	dvListMap.GenerateHash();
	nDvAct=dvListMap.vec.size();

	// Positional Derivatives

	// HERE -> function to calculate SurfCentroid (dc/dd)^T Hm (dc/dd)

	HPos.setZero(nDvAct,nDvAct);
	dPos.setZero(9,nDvAct);
	for(ii=0; ii<ni; ++ii){
		if (triIn.pointtype[ii]==2){
			subTemp=triRSVS.snakeDep->snaxs.find(triIn.pointind[ii]);
			subTemp1=triRSVS.meshDep->verts.find(triRSVS.snakeDep->snaxs(subTemp)->fromvert);
			subTemp2=triRSVS.meshDep->verts.find(triRSVS.snakeDep->snaxs(subTemp)->tovert);
			subTemp3=dvListMap.find(triIn.pointind[ii]);
			for(jj=0;jj<3;++jj){
				dPos(ii*3+jj,subTemp3)+=triRSVS.meshDep->verts(subTemp2)->coord[jj]
					-triRSVS.meshDep->verts(subTemp1)->coord[jj];
			}
		} else if (triIn.pointtype[ii]==3 && false){

		}
	}
	// Total

	HVal.setZero(9,9);
	dVal.setZero(1,9);
	dConstrPart.setZero(1,nDvAct);
	HConstrPart.setZero(nDvAct,nDvAct);

	if(isConstr){
		VolumeCalc.Calc();
		VolumeCalc.ReturnDatPoint(&retVal, &dValpnt, &HValpnt); 
		ArrayVec2MatrixXd(*HValpnt, HVal);
		ArrayVec2MatrixXd(*dValpnt, dVal);
		Deriv1stChainScalar(dVal, dPos,dConstrPart);
		Deriv2ndChainScalar(dVal,dPos,HVal,HPos,HConstrPart);
		// if (isDeriv){

		// cout << dConstrPart.rows() << " " << dConstrPart.cols() << " ";
		// }
		constrPart=*retVal;
		// Assign Constraint
		// and constraint derivative
		// and Hessian
		nCellTarg=triIn.connec.celltarg.size(); 
		for(ii=0; ii< nCellTarg;++ii){
			subTempVec=this->constrMap.findall(triIn.connec.celltarg[ii]);
			nj=subTempVec.size();
			for(jj=0; jj< nj; ++jj){
				if (subTempVec[jj]!=-1){
					this->constr[subTempVec[jj]] += triIn.connec.constrinfluence[ii]*constrPart;
					if(isDeriv){
						for(kk=0; kk< nDvAct; ++kk){
							this->dConstr(subTempVec[jj],this->dvMap.find(dvListMap.vec[kk])) += 
								triIn.connec.constrinfluence[ii]*dConstrPart(0,kk);
							dvCallConstr(this->dvMap.find(dvListMap.vec[kk]),0)++;
							for(ll=0; ll< nDvAct; ++ll){
								// TODO cross product with lagrangian
								this->HConstr(this->dvMap.find(dvListMap.vec[ll]),
									 this->dvMap.find(dvListMap.vec[kk])) += 
									 triIn.connec.constrinfluence[ii]*
									 HConstrPart(ll,kk)
									 *this->lagMult[subTempVec[jj]]*0;
							}
						}
					}
				} else {
					this->falseaccess++;
				}
			}
		}

	}

	HVal.setZero(9,9);
	dVal.setZero(1,9);
	dObjPart.setZero(1,nDvAct);	HObjPart.setZero(nDvAct,nDvAct);
	
	if(isObj){
		AreaCalc.Calc();
		AreaCalc.ReturnDatPoint(&retVal, &dValpnt, &HValpnt); 
		ArrayVec2MatrixXd(*HValpnt, HVal);
		ArrayVec2MatrixXd(*dValpnt, dVal);
		Deriv1stChainScalar(dVal, dPos,dObjPart);
		Deriv2ndChainScalar(dVal,dPos,HVal,HPos,HObjPart);
		objPart=*retVal;

		// Assign to main part of the object
		// assign objective function
		this->obj+=objPart; 
		// Assign objective derivative
		if(isDeriv){
			for(ii=0; ii< nDvAct; ++ii){
				this->dObj[this->dvMap.find(dvListMap.vec[ii])] += dObjPart(0,ii);
				for(jj=0; jj< nDvAct; ++jj){
					this->HObj(dvMap.find(dvListMap.vec[jj]),
						 this->dvMap.find(dvListMap.vec[ii])) += HObjPart(jj,ii);
				}
			}
		}
	}
	// Update active lists of design variables
	for(ii=0; ii< nDvAct; ++ii){
		this->isDvAct.at(dvMap.find(dvListMap.vec[ii])) = true;
	}
	nCellTarg=triIn.connec.celltarg.size(); 
	for(ii=0; ii< nCellTarg;++ii){
		subTempVec=constrMap.findall(triIn.connec.celltarg[ii]);
		nj=subTempVec.size();
		for(jj=0; jj< nj; ++jj){
			if (subTempVec[jj]!=-1){
				this->isConstrAct.at(subTempVec[jj]) = true;
			}
		}
	}

	// Assign Objective Hessian

	// Assign Constraint Hessian
}

void SQPcalc::CalcTriangleFD(const triangle& triIn, const triangulation &triRSVS,
	bool isObj, bool isConstr, bool isDeriv){
	/*Same as calctriangle but the volume derivative is calculated using a 
	Finite Difference.
	
	helps isolate errors in RSVSmath.
	*/

	int ii,ni,jj,nj,kk,ll,nCellTarg;
	int isCentre,subTemp,subTemp1,subTemp2,subTemp3,nDvAct;
	SurfCentroid centreCalc;
	Volume VolumeCalc, VolumeCalcFD;
	Area AreaCalc;
	vector<int> dvList,subTempVec;
	HashedVector<int,int> dvListMap;
	vector<vector<double> const *> veccoord;
	MatrixXd HPos, dPos;
	MatrixXd HVal, dVal;
	MatrixXd dConstrPart,HConstrPart, HObjPart;
	MatrixXd dObjPart;
	double constrPart, objPart;
	ArrayVec<double>* HValpnt=NULL, *dValpnt=NULL;
	double *retVal;


	veccoord.reserve(3);

	isCentre=0;
	ni=3;
	for(ii=0; ii<ni; ++ii){
		if(triIn.pointtype[ii]==1){
			veccoord.push_back(&(triRSVS.meshDep->verts.isearch(triIn.pointind[ii])->coord));
		} else if (triIn.pointtype[ii]==2){
			veccoord.push_back(&(triRSVS.snakeDep->snakeconn.verts.isearch(triIn.pointind[ii])->coord));
		} else if (triIn.pointtype[ii]==3){
			veccoord.push_back((triRSVS.trivert.isearch(triIn.pointind[ii])->coord.retPtr()));
			isCentre++;
		}
	}

	// Constr and objective
	VolumeCalc.assign(veccoord[0],veccoord[1],veccoord[2]);
	AreaCalc.assign(veccoord[0],veccoord[1],veccoord[2]);

	// Active DV lists

	for(ii=0; ii<ni; ++ii){
		if (triIn.pointtype[ii]==2){
			dvList.push_back(triIn.pointind[ii]);
		} else if (triIn.pointtype[ii]==3 && false){
			subTemp=triRSVS.trivert.find(triIn.pointind[ii]);
			nj=triRSVS.trisurf(subTemp)->indvert.size();
			for(jj=0;jj<nj;++jj){
				if (triRSVS.trisurf(subTemp)->typevert[jj]==2){
					dvList.push_back(triRSVS.trisurf(subTemp)->indvert[jj]);
				}
			}
		}
	}
	dvListMap.vec=dvList;
	sort(dvListMap.vec);
	unique(dvListMap.vec);
	dvListMap.GenerateHash();
	nDvAct=dvListMap.vec.size();

	// Positional Derivatives

	// HERE -> function to calculate SurfCentroid (dc/dd)^T Hm (dc/dd)

	HPos.setZero(nDvAct,nDvAct);
	dPos.setZero(9,nDvAct);
	for(ii=0; ii<ni; ++ii){
		if (triIn.pointtype[ii]==2){
			subTemp=triRSVS.snakeDep->snaxs.find(triIn.pointind[ii]);
			subTemp1=triRSVS.meshDep->verts.find(triRSVS.snakeDep->snaxs(subTemp)->fromvert);
			subTemp2=triRSVS.meshDep->verts.find(triRSVS.snakeDep->snaxs(subTemp)->tovert);
			subTemp3=dvListMap.find(triIn.pointind[ii]);
			for(jj=0;jj<3;++jj){
				dPos(ii*3+jj,subTemp3)+=triRSVS.meshDep->verts(subTemp2)->coord[jj]
					-triRSVS.meshDep->verts(subTemp1)->coord[jj];
			}
		} else if (triIn.pointtype[ii]==3 && false){

		}
	}
	// Total

	HVal.setZero(9,9);
	dVal.setZero(1,9);
	dConstrPart.setZero(1,nDvAct);
	HConstrPart.setZero(nDvAct,nDvAct);

	if(isConstr){
		VolumeCalc.CalcFD();
		VolumeCalc.ReturnDatPoint(&retVal, &dValpnt, &HValpnt); 
		ArrayVec2MatrixXd(*HValpnt, HVal);
		ArrayVec2MatrixXd(*dValpnt, dVal);
		Deriv1stChainScalar(dVal, dPos,dConstrPart);
		Deriv2ndChainScalar(dVal,dPos,HVal,HPos,HConstrPart);

		constrPart=*retVal;
		// Assign Constraint
		// and constraint derivative
		// and Hessian
		nCellTarg=triIn.connec.celltarg.size(); 
		for(ii=0; ii< nCellTarg;++ii){
			subTempVec=this->constrMap.findall(triIn.connec.celltarg[ii]);
			nj=subTempVec.size();
			for(jj=0; jj< nj; ++jj){
				if (subTempVec[jj]!=-1){
					this->constr[subTempVec[jj]] += triIn.connec.constrinfluence[ii]*constrPart;
					if(isDeriv){
						for(kk=0; kk< nDvAct; ++kk){
							this->dConstr(subTempVec[jj],this->dvMap.find(dvListMap.vec[kk])) += 
								triIn.connec.constrinfluence[ii]*dConstrPart(0,kk);
							dvCallConstr(this->dvMap.find(dvListMap.vec[kk]),0)++;
							for(ll=0; ll< nDvAct; ++ll){
								// TODO cross product with lagrangian
								this->HConstr(this->dvMap.find(dvListMap.vec[ll]),
									 this->dvMap.find(dvListMap.vec[kk])) += 
									 triIn.connec.constrinfluence[ii]*
									 HConstrPart(ll,kk)
									 *this->lagMult[subTempVec[jj]];
							}
						}
					}
				} else {
					this->falseaccess++;
				}
			}
		}

	}

	HVal.setZero(9,9);
	dVal.setZero(1,9);
	dObjPart.setZero(1,nDvAct);	HObjPart.setZero(nDvAct,nDvAct);
	
	if(isObj){
		AreaCalc.Calc();
		AreaCalc.ReturnDatPoint(&retVal, &dValpnt, &HValpnt); 
		ArrayVec2MatrixXd(*HValpnt, HVal);
		ArrayVec2MatrixXd(*dValpnt, dVal);
		Deriv1stChainScalar(dVal, dPos,dObjPart);
		Deriv2ndChainScalar(dVal,dPos,HVal,HPos,HObjPart);
		objPart=*retVal;

		// Assign to main part of the object
		// assign objective function
		this->obj+=objPart; 
		// Assign objective derivative
		// cout << endl << "dObjPart " << nDvAct << " " ;
		// PrintMatrix(dObjPart);
		// cout <<  endl << "done" <<endl;
		if(isDeriv){
			for(ii=0; ii< nDvAct; ++ii){
				this->dObj[this->dvMap.find(dvListMap.vec[ii])] += dObjPart(0,ii);
				for(jj=0; jj< nDvAct; ++jj){
					this->HObj(dvMap.find(dvListMap.vec[jj]),
						 this->dvMap.find(dvListMap.vec[ii])) += HObjPart(jj,ii);
				}
			}
		}
	}
	// Update active lists of design variables
	for(ii=0; ii< nDvAct; ++ii){
		this->isDvAct.at(dvMap.find(dvListMap.vec[ii])) = true;
	}
	nCellTarg=triIn.connec.celltarg.size(); 
	for(ii=0; ii< nCellTarg;++ii){
		subTempVec=constrMap.findall(triIn.connec.celltarg[ii]);
		nj=subTempVec.size();
		for(jj=0; jj< nj; ++jj){
			if (subTempVec[jj]!=-1){
				this->isConstrAct.at(subTempVec[jj]) = true;
			}
		}
	}

	// Assign Objective Hessian

	// Assign Constraint Hessian
}


void SQPcalc::CalcTriangleDirectVolume(const triangle& triIn, const triangulation &triRSVS,
	bool isObj, bool isConstr, bool isDeriv){
	/*Same as calctriangle but the volume derivative is calculated without intermediate
	steps
	
	helps isolate errors on dPos.
	*/

	int ii,ni,jj,nj,kk,ll,nCellTarg;
	int isCentre,subTemp,subTemp1,subTemp2,subTemp3,nDvAct;
	SurfCentroid centreCalc;
	Volume2 VolumeCalc2;
	Volume VolumeCalc;
	Area AreaCalc;
	vector<int> dvList,subTempVec, dvOrder;
	HashedVector<int,int> dvListMap;
	vector<vector<double> const *> veccoord, veccoordvol;
	vector<double> dvec;
	MatrixXd HPos, dPos;
	MatrixXd HVal, dVal;
	MatrixXd dConstrPart,HConstrPart, HObjPart;
	MatrixXd dObjPart;
	double constrPart, objPart;
	ArrayVec<double>* HValpnt=NULL, *dValpnt=NULL;
	ArrayVec<double>* HValpnt2=NULL, *dValpnt2=NULL;
	double * retVal=NULL;
	double * retVal2=NULL;


	veccoord.reserve(3);
	veccoordvol.assign(7,NULL);
	dvec.reserve(3);
	veccoordvol[0] = &dvec;
	isCentre=0;
	ni=3;
	for(ii=0; ii<ni; ++ii){
		if(triIn.pointtype[ii]==1){
			veccoord.push_back(&(triRSVS.meshDep->verts.isearch(triIn.pointind[ii])->coord));
		} else if (triIn.pointtype[ii]==2){
			veccoord.push_back(&(triRSVS.snakeDep->snakeconn.verts.isearch(triIn.pointind[ii])->coord));
		} else if (triIn.pointtype[ii]==3){
			veccoord.push_back((triRSVS.trivert.isearch(triIn.pointind[ii])->coord.retPtr()));
			isCentre++;
		}
	}

	// Volume calc assignement
	dvOrder.assign(3,0);
	ni = 3;
	for(ii=0; ii<ni; ++ii){
		if(triIn.pointtype[ii]==1){


			dvec[ii] = 0;
			veccoordvol[1+ii] = &(triRSVS.meshDep->verts.isearch(triIn.pointind[ii])->coord);
			veccoordvol[1+3+ii] = &(triRSVS.meshDep->verts.isearch(triIn.pointind[ii])->coord);

		} else if (triIn.pointtype[ii]==2){

			subTemp = triRSVS.snakeDep->snakeconn.verts.find(triIn.pointind[ii]);
			dvOrder[ii] = triIn.pointind[ii];

			dvec[ii] = triRSVS.snakeDep->snaxs(subTemp)->d;
			veccoordvol[1+ii] = &(triRSVS.meshDep->verts.isearch(
				triRSVS.snakeDep->snaxs(subTemp)->fromvert)->coord);
			veccoordvol[1+ii+3] = &(triRSVS.meshDep->verts.isearch(
				triRSVS.snakeDep->snaxs(subTemp)->tovert)->coord);

		} else if (triIn.pointtype[ii]==3){

			dvec[ii] = 0;
			veccoordvol[1+ii] = (triRSVS.trivert.isearch(triIn.pointind[ii])->coord.retPtr());
			veccoordvol[1+ii+3] = (triRSVS.trivert.isearch(triIn.pointind[ii])->coord.retPtr());
			isCentre++;
		}
	}

	// Constr and objective
	VolumeCalc2.assign(veccoordvol); /// <-----Change assignement 
	VolumeCalc.assign(veccoord[0],veccoord[1],veccoord[2]); /// <-----Change assignement 
	AreaCalc.assign(veccoord[0],veccoord[1],veccoord[2]);

	// Active DV lists
	for(ii=0; ii<ni; ++ii){
		if (triIn.pointtype[ii]==2){
			dvList.push_back(triIn.pointind[ii]);
		} else if (triIn.pointtype[ii]==3 && false){
			subTemp=triRSVS.trivert.find(triIn.pointind[ii]);
			nj=triRSVS.trisurf(subTemp)->indvert.size();
			for(jj=0;jj<nj;++jj){
				if (triRSVS.trisurf(subTemp)->typevert[jj]==2){
					dvList.push_back(triRSVS.trisurf(subTemp)->indvert[jj]);
				}
			}
		}
	}
	dvListMap.vec=dvList;
	sort(dvListMap.vec);
	unique(dvListMap.vec);
	dvListMap.GenerateHash();
	nDvAct=dvListMap.vec.size();

	// Positional Derivatives

	// HERE -> function to calculate SurfCentroid (dc/dd)^T Hm (dc/dd)

	HPos.setZero(nDvAct,nDvAct);
	dPos.setZero(9,nDvAct);

	for(ii=0; ii<ni; ++ii){
		if (triIn.pointtype[ii]==2){
			subTemp=triRSVS.snakeDep->snaxs.find(triIn.pointind[ii]);
			subTemp1=triRSVS.meshDep->verts.find(triRSVS.snakeDep->snaxs(subTemp)->fromvert);
			subTemp2=triRSVS.meshDep->verts.find(triRSVS.snakeDep->snaxs(subTemp)->tovert);
			subTemp3=dvListMap.find(triIn.pointind[ii]);
			for(jj=0;jj<3;++jj){
				dPos(ii*3+jj,subTemp3)+=triRSVS.meshDep->verts(subTemp2)->coord[jj]
					-triRSVS.meshDep->verts(subTemp1)->coord[jj];
			}
		} else if (triIn.pointtype[ii]==3 && false){

		}
	}
	// Total

	HVal.setZero(9,9);
	dVal.setZero(1,9);
	dConstrPart.setZero(1,nDvAct);
	HConstrPart.setZero(nDvAct,nDvAct);

	if(isConstr){
		VolumeCalc2.Calc();
		VolumeCalc.Calc();
		VolumeCalc.ReturnDatPoint(&retVal, &dValpnt, &HValpnt); 
		VolumeCalc2.ReturnDatPoint(&retVal2, &dValpnt2, &HValpnt2); 
		ArrayVec2MatrixXd(*HValpnt2, HConstrPart);
		ArrayVec2MatrixXd(*dValpnt2, dConstrPart);
		// Deriv1stChainScalar(dVal, dPos,dConstrPart);/// <-----Change NOT NEEDED
		// Deriv2ndChainScalar(dVal,dPos,HVal,HPos,HConstrPart);/// <-----Change NOT NEEDED
		// if (isDeriv){

		// cout << dConstrPart.rows() << " " << dConstrPart.cols() << " ";
		// }
		if(fabs(*retVal-*retVal2)>1e-10){
			// cout << *retVal << " " << *retVal2 << endl;
			this->falseaccess++;
		}
		constrPart=*retVal2;
		
		if(dvListMap.find(1044)!=-1 && isDeriv){
			cout << endl;
			DisplayVector(dvOrder);
			cout << endl;
			PrintMatrix(dConstrPart);
			cout << endl;
		}
		// if((dvListMap.find(1044)!=-1 )  && isDeriv){
		// 	std::vector<double> v={constrPart};
		// 	PrintMatrixFile(triIn.connec.celltarg,"matrices\\matrix_dConstrPart_inter2.txt");
		// 	PrintMatrixFile(triIn.connec.constrinfluence,"matrices\\matrix_dConstrPart_inter2.txt");
		// 	PrintMatrixFile(dvListMap.vec,"matrices\\matrix_dConstrPart_inter2.txt");
		// 	PrintMatrixFile(dConstrPart,"matrices\\matrix_dConstrPart_inter2.txt");
		// 	PrintMatrixFile(v,"matrices\\matrix_dConstrPart_inter2.txt");
		// 	PrintMatrixFile(dVal,"matrices\\matrix_dVal_inter2.txt");
		// 	PrintMatrixFile(dPos,"matrices\\matrix_dPos_inter2.txt");
		// 	// cout << endl;
		// 	// DisplayVector(*veccoord[0]);
		// 	// cout << endl;
		// 	// DisplayVector(*veccoord[1]);
		// 	// cout << endl;
		// 	// DisplayVector(*veccoord[2]);
		// 	// cout << endl;
		// 	// throw invalid_argument("");
		// }
		// Assign Constraint
		// and constraint derivative
		// and Hessian
		nCellTarg=triIn.connec.celltarg.size(); 
		for(ii=0; ii< nCellTarg;++ii){
			subTempVec=this->constrMap.findall(triIn.connec.celltarg[ii]);
			nj=subTempVec.size();
			for(jj=0; jj< nj; ++jj){
				if (subTempVec[jj]!=-1){
					this->constr[subTempVec[jj]] += triIn.connec.constrinfluence[ii]*constrPart;
					if(isDeriv){
						for(kk=0; kk< 3; ++kk){
							if(dvOrder[kk]!=0){

								this->dConstr(subTempVec[jj],this->dvMap.find(dvOrder[kk])) += 
									triIn.connec.constrinfluence[ii]*dConstrPart(0,kk);
								dvCallConstr(this->dvMap.find(dvOrder[kk]),0)++;
								for(ll=0; ll< 3; ++ll){
									if(dvOrder[ll]!=0){
									// TODO cross product with lagrangian
										this->HConstr(this->dvMap.find(dvOrder[ll]),
											 this->dvMap.find(dvOrder[kk])) += 
											 triIn.connec.constrinfluence[ii]*
											 HConstrPart(ll,kk)
											 *this->lagMult[subTempVec[jj]];
									}
								}
							}
						}
					}
				} else {
					this->falseaccess++;
				}
			}
		}

	}

	HVal.setZero(9,9);
	dVal.setZero(1,9);
	dObjPart.setZero(1,nDvAct);	HObjPart.setZero(nDvAct,nDvAct);
	
	if(isObj){
		AreaCalc.Calc();
		AreaCalc.ReturnDatPoint(&retVal, &dValpnt, &HValpnt); 
		ArrayVec2MatrixXd(*HValpnt, HVal);
		ArrayVec2MatrixXd(*dValpnt, dVal);
		Deriv1stChainScalar(dVal, dPos,dObjPart);
		Deriv2ndChainScalar(dVal,dPos,HVal,HPos,HObjPart);
		objPart=*retVal;

		// Assign to main part of the object
		// assign objective function
		this->obj+=objPart; 
		// Assign objective derivative
		// cout << endl << "dObjPart " << nDvAct << " " ;
		// PrintMatrix(dObjPart);
		// cout <<  endl << "done" <<endl;
		if(isDeriv){
			for(ii=0; ii< nDvAct; ++ii){
				this->dObj[this->dvMap.find(dvListMap.vec[ii])] += dObjPart(0,ii);
				for(jj=0; jj< nDvAct; ++jj){
					this->HObj(dvMap.find(dvListMap.vec[jj]),
						 this->dvMap.find(dvListMap.vec[ii])) += HObjPart(jj,ii);
				}
			}
		}
	}
	// Update active lists of design variables
	for(ii=0; ii< nDvAct; ++ii){
		this->isDvAct.at(dvMap.find(dvListMap.vec[ii])) = true;
	}
	nCellTarg=triIn.connec.celltarg.size(); 
	for(ii=0; ii< nCellTarg;++ii){
		subTempVec=constrMap.findall(triIn.connec.celltarg[ii]);
		nj=subTempVec.size();
		for(jj=0; jj< nj; ++jj){
			if (subTempVec[jj]!=-1){
				this->isConstrAct.at(subTempVec[jj]) = true;
			}
		}
	}

	// Assign Objective Hessian

	// Assign Constraint Hessian
}

bool SQPcalc::PrepareMatricesForSQP(
	MatrixXd &dConstrAct,
	MatrixXd &HConstrAct, 
	MatrixXd &HObjAct,
	RowVectorXd &dObjAct,
	VectorXd &constrAct,
	VectorXd &lagMultAct
	){

	

	int ii, jj, nDvAct, nConstrAct;

	subDvAct.reserve(nDv);
	subDvAct.clear();
	nDvAct=0;
	for(ii=0; ii<nDv; ++ii){
		nDvAct += int(isDvAct.at(ii));
		if(isDvAct.at(ii)){
			subDvAct.push_back(ii);
		}
	}

	subConstrAct.reserve(nConstr);
	subConstrAct.clear();
	nConstrAct=0;
	for(ii=0; ii<nConstr; ++ii){
		nConstrAct += int(isConstrAct.at(ii));
		if(isConstrAct.at(ii)){
			subConstrAct.push_back(ii);
		}
	}


	dConstrAct.setZero(nConstrAct,nDvAct);
	for(ii=0; ii<nConstrAct;++ii){
		for(jj=0; jj<nDvAct;++jj){
			dConstrAct(ii,jj)=dConstr(subConstrAct[ii],subDvAct[jj]);
		}
	}
	HConstrAct.setZero(nDvAct,nDvAct); 
	for(ii=0; ii<nDvAct;++ii){
		for(jj=0; jj<nDvAct;++jj){
			HConstrAct(ii,jj)=HConstr(subDvAct[ii],subDvAct[jj]);
		}
	}
	HObjAct.setZero(nDvAct,nDvAct);  
	for(ii=0; ii<nDvAct;++ii){
		for(jj=0; jj<nDvAct;++jj){
			HObjAct(ii,jj)=HObj(subDvAct[ii],subDvAct[jj]);
		}
	}
	dObjAct.setZero(nDvAct);
	for(jj=0; jj<nDvAct;++jj){
		dObjAct[jj]=dObj[subDvAct[jj]];
	}
	constrAct.setZero(nConstrAct);
	for(jj=0; jj<nConstrAct;++jj){
		constrAct[jj]=constr[subConstrAct[jj]] -constrTarg[subConstrAct[jj]];
	}
	lagMultAct.setZero(nConstrAct);

	// DisplayVector(isDvAct);
	// DisplayVector(isConstrAct);
	HLag = HObjAct;
	return(nDvAct>0);
}

void SQPcalc::CheckAndCompute(){
	int ii;
	bool computeFlag;
	MatrixXd dConstrAct,HConstrAct, HObjAct;
	RowVectorXd dObjAct;
	VectorXd constrAct, lagMultAct, deltaDVAct;

	computeFlag = PrepareMatricesForSQP(
		dConstrAct,HConstrAct, HObjAct,dObjAct,
		constrAct,lagMultAct
		);
	computeFlag=true;
	if (computeFlag){
		ComputeSQPstep(
			dConstrAct,dObjAct,
			constrAct,lagMultAct
			);
	} else {
		for (ii=0; ii<nDv; ++ii){
			deltaDV[ii]=0.0;
			
		}
	}

}

void SQPcalc::ComputeSQPstep(
	MatrixXd &dConstrAct,
	RowVectorXd &dObjAct,
	VectorXd &constrAct,
	VectorXd &lagMultAct
	){

	VectorXd  deltaDVAct;
	MatrixXd temp1, temp2;
	bool isNan, isLarge;
	int ii, ni;

	ColPivHouseholderQR<MatrixXd> HLagSystem(HLag);
	// HouseholderQR<MatrixXd> HLagSystem(HLag);
	// LLT<MatrixXd> HLagSystem(HLag);
	// PartialPivLU<MatrixXd> HLagSystem(HLag);

	
	temp1 = HLagSystem.solve(dConstrAct.transpose());
	temp2 = HLagSystem.solve(dObjAct.transpose());

	lagMultAct = (
			dConstrAct*(temp1)
		).colPivHouseholderQr().solve(
		// ).householderQr().solve(
		// ).llt().solve(
		// ).partialPivLu().solve()
			constrAct - (dConstrAct*(temp2))
		);

	isLarge = false;
	isNan = false;
	ni = lagMultAct.size();
	for (ii=0; ii<ni; ++ii){
		if(lagMultAct[ii]<-limLag){
			lagMultAct[ii]=-limLag;
			isLarge=true;
		}else if(lagMultAct[ii]>limLag){
			lagMultAct[ii]=limLag;
			isLarge=true;
		} else if(isnan(lagMultAct[ii])){
			lagMultAct[ii]=0.0;
			isNan=true;
		}
	}
	// if (isNan){
	// isLarge = false;
	// 	deltaDVAct = -dConstrAct.transpose()*lagMultAct;
	// }else 
	
	if(isLarge) {

		// PrintMatrixFile(dConstrAct, "matrix_dConstrAct.txt");
	 	deltaDVAct = -dConstrAct.bdcSvd(ComputeThinU | ComputeThinV).solve(constrAct);

	} else {

		deltaDVAct = - (HLagSystem.solve(dObjAct.transpose() 
						+ dConstrAct.transpose()*lagMultAct));

	}

	ni = subDvAct.size();
	for (ii=0; ii<ni; ++ii){
		deltaDV[subDvAct[ii]]=deltaDVAct[ii];
		// deltaDV[subDvAct[ii]]=dConstrAct.col(subDvAct[ii])[1];
		
	}
	
	ni = subConstrAct.size();
	lagMult.setZero(nConstr);
	if(isLarge){
		for (ii=0; ii<nConstr; ++ii){
			lagMult[subConstrAct[ii]]=0;
		}

	} else {
		for (ii=0; ii<ni; ++ii){
			lagMult[subConstrAct[ii]]=lagMultAct[ii];
			isNan = isNan || isnan(lagMultAct[ii]);
		}
	}
	if (false){
		Print2Screen(3);
		DisplayVector(isConstrAct);
		DisplayVector(isDvAct);
	}

}

void ArrayVec2MatrixXd(const ArrayVec<double> &arrayIn, MatrixXd &matOut){

	int nR=0,nC=0;
	arrayIn.size(nR,nC);
	matOut.setZero(nR,nC);

	for(int i=0; i< nR; ++i){
		for(int j=0; j< nC; ++j){
			matOut(i,j)=arrayIn[i][j];
		}
	}

}

void VecBy3DimArray(const MatrixXd &vec, const MatrixXd &arr3dim, MatrixXd &retArray){
	/*
	3 Dimensional array is expected as 
	 					xc 				 					yc 				 				zc 				
	  x0 x1  xn y0 y1  yn  z0 z1  zn	x0 x1  xn y0 y1  yn  z0 z1  zn	    x0 x1  xn y0 y1  yn  z0 z1  zn
 x0 [                               ]  [                               ]  [                               ] 
 x1 [                               ]  [                               ]  [                               ] 
 xn [                               ]  [                               ]  [                               ] 
 y0 [                               ]  [                               ]  [                               ] 
 y1 [                               ]  [                               ]  [                               ] 
 yn [                               ]  [                               ]  [                               ] 
 z0 [                               ]  [                               ]  [                               ] 
 z1 [                               ]  [                               ]  [                               ] 
 zn [                               ]  [                               ]  [                               ] 
	vector expected as a row vector
	[xc yc zc]

	will sum the three Hessians according to the values in the vector
	*/

	int nRow,nCol,nVec,nColFin;
	int ii,jj,kk;
	nRow=arr3dim.rows();
	nCol=arr3dim.cols();
	nVec=vec.rows();
	nColFin=nCol/nVec;

	#ifdef SAFE_ALGO
	// size checks
	if ((nVec*nColFin)!=nCol){
		throw invalid_argument("Sizes do not match in 3D array collapse");
	}
	#endif

	retArray.setZero(nRow,nColFin);
	// needs to add checks for matching sizes
	for (ii=0;ii<nRow;ii++){
		for (jj=0;jj<nColFin;jj++){
			for (kk=0;kk<nVec;kk++){
				retArray(ii,jj)+=arr3dim(ii,jj*nVec+kk)*vec(0,kk);
			}	
		}	
	}



}



void Deriv1stChainScalar(const MatrixXd &dSdc,const MatrixXd &dcdd, MatrixXd &dSdd){
	dSdd=dSdc*dcdd;
}

void Deriv2ndChainScalar(const MatrixXd &dSdc,const MatrixXd &dcdd,
	const MatrixXd &HSc,const MatrixXd &Hcd,MatrixXd &HSd){
	
	VecBy3DimArray(dSdc, Hcd, HSd);
	HSd = HSd + (dcdd.transpose()*HSc*dcdd);

}

void PrintMatrix(const MatrixXd mat){
	int ii,jj, ni, nj;

	ni=mat.rows();
	nj=mat.cols();
	for(ii=0;ii<ni;++ii){
		for(jj=0;jj<nj;++jj){
			cout << mat(ii,jj) << " ";
		}
		cout << endl;
	}
}
void PrintMatrixFile(const MatrixXd mat, const char * name){
	int ii,jj, ni, nj;
	ofstream myfile;
	ni=mat.rows();
	nj=mat.cols();
	myfile.open(name, ios::app);
	for(ii=0;ii<ni;++ii){
		for(jj=0;jj<nj;++jj){
			myfile  << mat(ii,jj) << " ";
		}
		myfile  << endl;
	}
	myfile  << endl;
	myfile.close();

}
template <class T> void PrintMatrixFile(const vector<T> &mat, const char * name){
	int ii,ni;
	ofstream myfile;
	ni=mat.size();

	myfile.open(name, ios::app);
	myfile.precision(16);
	myfile << std::scientific;
	for(ii=0;ii<ni;++ii){
		myfile  << mat[ii] << " ";
	}
	myfile  << endl;
	myfile.close();

}
void PrintMatrix(const VectorXd mat){
	int ii, ni;

	ni=mat.size();
	for(ii=0;ii<ni;++ii){
		cout << mat[ii] << " ";
		cout << endl;
	}
}
void PrintMatrix(const RowVectorXd mat){
	int ii, ni;

	ni=mat.size();
	for(ii=0;ii<ni;++ii){
		cout << mat[ii] << " ";
		
	}
	cout << endl;
}

/// INSPIRATION
/*
void SnakeSurfaceCentroid_fun(coordvec &coord,const surf &surfin, const mesh& meshin){ 
	int ii,n;
	vector<int> vertind;
	vector<vector<double> const *> veccoord;
	SurfCentroid tempCalc;
	ArrayVec<double> tempCoord,jac,hes;

	coord.assign(0,0,0);
	n=int(surfin.edgeind.size());

	veccoord.reserve(n);
	ConnVertFromConnEdge(meshin, surfin.edgeind,vertind);

	for(ii=0; ii<n; ++ii){
		veccoord.push_back(&(meshin.verts.isearch(vertind[ii])->coord));
	}

	tempCalc.assign(veccoord);
	tempCalc.Calc();

	tempCalc.ReturnDat(tempCoord,jac,hes);
	coord.assign(tempCoord[0][0],tempCoord[1][0],tempCoord[2][0]);
}

void HybridSurfaceCentroid_fun(coordvec &coord,const trianglesurf &surfin, const mesh& meshin,
	const mesh& snakeconn){ 
	int ii,n;
	vector<int> vertind;
	vector<vector<double> const *> veccoord;
	SurfCentroid tempCalc;
	ArrayVec<double> tempCoord,jac,hes;

	coord.assign(0,0,0);

	n=surfin.indvert.size();
	for(ii=0; ii<n; ++ii){
		if(surfin.typevert[ii]==1){
			veccoord.push_back(&(meshin.verts.isearch(surfin.indvert[ii])->coord));
		} else if (surfin.typevert[ii]==2){
			veccoord.push_back(&(snakeconn.verts.isearch(surfin.indvert[ii])->coord));
		}
	}

	tempCalc.assign(veccoord);
	tempCalc.Calc();

	tempCalc.ReturnDat(tempCoord,jac,hes);
	coord.assign(tempCoord[0][0],tempCoord[1][0],tempCoord[2][0]);
}
*/

