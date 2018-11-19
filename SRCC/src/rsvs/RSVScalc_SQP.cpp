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
// SQP calculation functions
//==========================================

bool RSVScalc::PrepareMatricesForSQP(
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

void RSVScalc::CheckAndCompute(int calcMethod){
	int ii;
	bool computeFlag;
	MatrixXd dConstrAct,HConstrAct, HObjAct;
	RowVectorXd dObjAct;
	VectorXd constrAct, lagMultAct, deltaDVAct;

	computeFlag = PrepareMatricesForSQP(
		dConstrAct,HConstrAct, HObjAct,dObjAct,
		constrAct,lagMultAct
		);
	// computeFlag=true;
	if (computeFlag){
		ComputeSQPstep(calcMethod,
			dConstrAct,dObjAct,
			constrAct,lagMultAct
			);
	} else {
		for (ii=0; ii<nDv; ++ii){
			deltaDV[ii]=0.0;
			
		}
	}
}

void RSVScalc::ComputeSQPstep(
	int calcMethod,
	MatrixXd &dConstrAct,
	RowVectorXd &dObjAct,
	VectorXd &constrAct,
	VectorXd &lagMultAct
	){

	VectorXd  deltaDVAct;
	bool isNan, isLarge;
	int ii, ni;
	cout << calcMethod << endl;
	switch(calcMethod){
		case 1:
			SQPstep<Eigen::HouseholderQR>(*this, dConstrAct, dObjAct,
				constrAct, lagMultAct,
				deltaDVAct, isNan, isLarge);
			// cout << "case 1:" << endl;

		break;
		case 2:
			SQPstep<Eigen::ColPivHouseholderQR>(*this, dConstrAct, dObjAct,
				constrAct, lagMultAct,
				deltaDVAct, isNan, isLarge);
			// cout << "case 2:" << endl;
		break;
		case 3: 
			// Eigen::LLT is inconveniently a 2 parameter template so
			// a full type is passed
			SQPstep<Eigen::LLT<MatrixXd>>(*this, dConstrAct, dObjAct,
				constrAct, lagMultAct,
				deltaDVAct, isNan, isLarge);
			// cout << "case 3:" << endl;
		break;
		case 4:
			SQPstep<Eigen::PartialPivLU>(*this, dConstrAct, dObjAct,
				constrAct, lagMultAct,
				deltaDVAct, isNan, isLarge);
			// cout << "case 4:" << endl;
		break;
		default:
			SQPstep<Eigen::HouseholderQR>(*this, dConstrAct, dObjAct,
				constrAct, lagMultAct,
				deltaDVAct, isNan, isLarge);
			// cout << "case default:" << endl;
		break;
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