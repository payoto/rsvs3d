#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath> 
#include <Eigen>

#include "snake.hpp"
#include "snakevel.hpp"
#include "vectorarray.hpp"
#include "RSVScalc.hpp"
#include "matrixtools.hpp"

using namespace std; 
using namespace Eigen; 

/*
Implementation of the core calculation functions of
RSVScalc, this is done to reduce the size of the objects and 
keep compilation time manageable.

*/


//==========================================
// Core functions and SQPstep template
//==========================================

void ResizeLagrangianMultiplier(const RSVScalc &calcobj, 
	VectorXd &lagMultAct, 
	bool &isLarge, bool &isNan){
	/*
	Resizes the lagrangian multiplier LagMultAct based on whether any of its values
	are nan or too large.

	This uses the RSVScalc object to guide it.
	*/
	int ii, ni;

	isLarge = false;
	isNan = false;
	ni = lagMultAct.size();
	for (ii=0; ii<ni; ++ii){
		if(lagMultAct[ii]<-calcobj.limLag){
			lagMultAct[ii]=-calcobj.limLag;
			isLarge=true;
		}else if(lagMultAct[ii]>calcobj.limLag){
			lagMultAct[ii]=calcobj.limLag;
			isLarge=true;
		} else if(isnan(lagMultAct[ii])){
			lagMultAct[ii]=0.0;
			isNan=true;
		}
	}
}
/*
void SQPstep(const RSVScalc &calcobj,
	const MatrixXd &dConstrAct, const RowVectorXd &dObjAct,
	const VectorXd &constrAct, VectorXd &lagMultAct,
	VectorXd &deltaDVAct, bool &isNan, bool &isLarge){


	MatrixXd temp1, temp2;

	// ColPivHouseholderQR<MatrixXd> HLagSystem(HLag);
	HouseholderQR<MatrixXd> HLagSystem(calcobj.HLag);
	// LLT<MatrixXd> HLagSystem(HLag);
	// PartialPivLU<MatrixXd> HLagSystem(HLag);


	
	temp1 = HLagSystem.solve(dConstrAct.transpose());
	temp2 = HLagSystem.solve(dObjAct.transpose());

	lagMultAct = (
			dConstrAct*(temp1)
		// ).colPivHouseholderQr().solve(
		).householderQr().solve(
		// ).llt().solve(
		// ).partialPivLu().solve(
			constrAct - (dConstrAct*(temp2))
		);

	ResizeLagrangianMultiplier(calcobj, lagMultAct, isLarge, isNan);
	
	if(isLarge) {

		// PrintMatrixFile(dConstrAct, "matrix_dConstrAct.txt");
	 	deltaDVAct = -dConstrAct.bdcSvd(ComputeThinU | ComputeThinV).solve(constrAct);

	} else {

		deltaDVAct = - (HLagSystem.solve(dObjAct.transpose() 
						+ dConstrAct.transpose()*lagMultAct));
	}

}*/



//==========================================
// Core class functions
//==========================================

void RSVScalc::CalcTriangle(const triangle& triIn,
	const triangulation &triRSVS,
	bool isObj, bool isConstr, bool isDeriv){


	int ii,ni,jj,nj,kk,ll,nCellTarg;
	int subTemp,subTemp1,subTemp2,subTemp3,nDvAct;
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
	ni=3;
	for(ii=0; ii<ni; ++ii){
		if(triIn.pointtype[ii]==1){
			veccoord.push_back(&(triRSVS.meshDep->verts.isearch(triIn.pointind[ii])->coord));
		} else if (triIn.pointtype[ii]==2){
			veccoord.push_back(&(triRSVS.snakeDep->snakeconn.verts.isearch(triIn.pointind[ii])->coord));
		} else if (triIn.pointtype[ii]==3){
			veccoord.push_back((triRSVS.trivert.isearch(triIn.pointind[ii])->coord.retPtr()));
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

void RSVScalc::CalcTriangleFD(const triangle& triIn,
	const triangulation &triRSVS,
	bool isObj, bool isConstr, bool isDeriv){
	/*Same as calctriangle but the volume derivative is calculated using a 
	Finite Difference.
	
	helps isolate errors in RSVSmath.
	*/

	int ii,ni,jj,nj,kk,ll,nCellTarg;
	int subTemp,subTemp1,subTemp2,subTemp3,nDvAct;
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

	ni=3;
	for(ii=0; ii<ni; ++ii){
		if(triIn.pointtype[ii]==1){
			veccoord.push_back(&(triRSVS.meshDep->verts.isearch(triIn.pointind[ii])->coord));
		} else if (triIn.pointtype[ii]==2){
			veccoord.push_back(&(triRSVS.snakeDep->snakeconn.verts.isearch(triIn.pointind[ii])->coord));
		} else if (triIn.pointtype[ii]==3){
			veccoord.push_back((triRSVS.trivert.isearch(triIn.pointind[ii])->coord.retPtr()));
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

void RSVScalc::CalcTriangleDirectVolume(const triangle& triIn,
	const triangulation &triRSVS,
	bool isObj, bool isConstr, bool isDeriv){
	/*Same as calctriangle but the volume derivative is calculated without intermediate
	steps
	
	helps isolate errors on dPos.
	*/

	int ii,ni,jj,nj,kk,ll,nCellTarg;
	int subTemp,subTemp1,subTemp2,subTemp3,nDvAct;
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
	ni=3;
	for(ii=0; ii<ni; ++ii){
		if(triIn.pointtype[ii]==1){
			veccoord.push_back(&(triRSVS.meshDep->verts.isearch(triIn.pointind[ii])->coord));
		} else if (triIn.pointtype[ii]==2){
			veccoord.push_back(&(triRSVS.snakeDep->snakeconn.verts.isearch(triIn.pointind[ii])->coord));
		} else if (triIn.pointtype[ii]==3){
			veccoord.push_back((triRSVS.trivert.isearch(triIn.pointind[ii])->coord.retPtr()));
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

			dvec[ii] = 0; // dvec used as a pointer
			veccoordvol[1+ii] = (triRSVS.trivert.isearch(triIn.pointind[ii])->coord.retPtr());
			veccoordvol[1+ii+3] = (triRSVS.trivert.isearch(triIn.pointind[ii])->coord.retPtr());
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
		// 	PrintMatrixFile(triIn.connec.celltarg,"matrices/matrix_dConstrPart_inter2.txt");
		// 	PrintMatrixFile(triIn.connec.constrinfluence,"matrices/matrix_dConstrPart_inter2.txt");
		// 	PrintMatrixFile(dvListMap.vec,"matrices/matrix_dConstrPart_inter2.txt");
		// 	PrintMatrixFile(dConstrPart,"matrices/matrix_dConstrPart_inter2.txt");
		// 	PrintMatrixFile(v,"matrices/matrix_dConstrPart_inter2.txt");
		// 	PrintMatrixFile(dVal,"matrices/matrix_dVal_inter2.txt");
		// 	PrintMatrixFile(dPos,"matrices/matrix_dPos_inter2.txt");
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