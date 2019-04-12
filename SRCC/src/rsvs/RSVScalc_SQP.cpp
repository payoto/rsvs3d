#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath> 
#include <Eigen>

#include "snake.hpp"
#include "triangulate.hpp"
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
			dConstrAct(ii,jj)=dConstr(subConstrAct[ii],subDvAct[jj])
				// /constrTarg[subConstrAct[ii]]
				;
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
		constrAct[jj]=(constr[subConstrAct[jj]] -constrTarg[subConstrAct[jj]])
			// /constrTarg[subConstrAct[jj]]
			;
	}
	lagMultAct.setZero(nConstrAct);

	// DisplayVector(isDvAct);
	// DisplayVector(isConstrAct);
	HLag = HObjAct;
	return(nDvAct>0);
}

bool RSVScalc::PrepareMatricesForSQPSensitivity(
	const MatrixXd &dConstrAct,
	const MatrixXd &HConstrAct, 
	const MatrixXd &HObjAct,
	MatrixXd &sensMult,
	MatrixXd &sensInv,
	MatrixXd &sensRes
	) const {

	int nDvAct, nConstrAct;
	MatrixXd HLagAct;

	// Find active design variables and active constraint
	nDvAct = this->subDvAct.size();
	nConstrAct  = this->subConstrAct.size();


	/// sensitivity requires : [H_d L , J_d constr; J_d constr ^T, 0]^-1
	/// [Dd^TD_v L (=0s) ; D_V constr (=eye(lagMult))]
 	
	sensMult.resize(nDvAct+nConstrAct, nDvAct+nConstrAct);
	sensInv.resize(nDvAct+nConstrAct, nConstrAct);
	sensRes.resize(nDvAct+nConstrAct, nConstrAct);

	sensInv << HConstrAct+HObjAct , dConstrAct.transpose(), 
		dConstrAct, MatrixXd::Zero(nConstrAct, nConstrAct);

	sensMult.setZero();
	for (int i = 0; i < nConstrAct; ++i)
	{
		sensMult(nDvAct+i, i) = this->lagMult[this->subConstrAct[i]];
	}
	sensRes.setZero();

	// do some size checks on the constructed matrices.

	return(nDvAct>0 && nConstrAct>0);
}

void RSVScalc::CheckAndCompute(int calcMethod, bool sensCalc){
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
	if(sensCalc){ 
		MatrixXd sensMult, sensInv, sensRes;
		bool computeSens = this->PrepareMatricesForSQPSensitivity(dConstrAct,
			HConstrAct, HObjAct, sensMult, sensInv, sensRes);
		this->sensDv.setZero(this->nDv, this->nConstr);
		if (computeSens)
		{
			this->ComputeSQPsens(calcMethod, sensMult, sensInv, sensRes);
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

	int ni = this->subDvAct.size();
	VectorXd  deltaDVAct(ni);
	bool isNan, isLarge, rerunDefault=true, attemptConstrOnly;
	int ii;

	/*
	This while loop while expensive allows to use the most stable algorithm
	*/
	deltaDVAct.setZero();
	attemptConstrOnly=false;
	if((calcMethod%10)!=calcMethod){
		calcMethod=calcMethod%10;
		attemptConstrOnly =true;
	}
	while(rerunDefault){
		switch(calcMethod){
			case 1:
				rerunDefault=SQPstep<Eigen::HouseholderQR>(*this, dConstrAct, 
					dObjAct, constrAct, lagMultAct,
					deltaDVAct, isNan, isLarge,attemptConstrOnly);
				// cout << "case 1:" << endl;

			break;
			case 2:
				rerunDefault=SQPstep<Eigen::ColPivHouseholderQR>(*this,
					dConstrAct, dObjAct, constrAct, lagMultAct,
					deltaDVAct, isNan, isLarge,attemptConstrOnly);
				// cout << "case 2:" << endl;
			break;
			case 3: 
				// Eigen::LLT is inconveniently a 2 parameter template so
				// a full type is passed
				rerunDefault=SQPstep<Eigen::LLT<MatrixXd>>(*this, dConstrAct, 
					dObjAct, constrAct, lagMultAct,
					deltaDVAct, isNan, isLarge,attemptConstrOnly);
				// cout << "case 3:" << endl;
			break;
			case 4:
				rerunDefault=SQPstep<Eigen::PartialPivLU>(*this, dConstrAct, 
					dObjAct, constrAct, lagMultAct,
					deltaDVAct, isNan, isLarge,attemptConstrOnly);
				// cout << "case 4:" << endl;
			break;
			default:
				SQPstep<Eigen::ColPivHouseholderQR>(*this, dConstrAct, dObjAct,
					constrAct, lagMultAct,
					deltaDVAct, isNan, isLarge,true);
				rerunDefault=false;
				// cout << "case default:" << endl;
			break;
		}
		if(attemptConstrOnly){
			rerunDefault=false;
		}
		if(rerunDefault){
			calcMethod=0;
		}
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

void RSVScalc::ComputeSQPsens(
	int calcMethod,
	const Eigen::MatrixXd &sensMult,
	const Eigen::MatrixXd &sensInv,
	Eigen::MatrixXd &sensRes
	){


	switch(calcMethod){
		case 1:
			SQPsens<Eigen::HouseholderQR>(sensMult, sensInv, sensRes);
		break;
		case 2:
			SQPsens<Eigen::ColPivHouseholderQR>(sensMult, sensInv, sensRes);
		break;
		case 3: 
			// Eigen::LLT is inconveniently a 2 parameter template so
			// a full type is passed
			SQPsens<Eigen::LLT<MatrixXd>>(sensMult, sensInv, sensRes);
		break;
		case 4:
			SQPsens<Eigen::PartialPivLU>(sensMult, sensInv, sensRes);
		break;
		default:
			SQPsens<Eigen::ColPivHouseholderQR>(sensMult, sensInv, sensRes);
		break;
	}

	int nConstrAct, nDvAct;
	nConstrAct = this->subConstrAct.size();
	nDvAct = this->subDvAct.size();
	for (int i = 0; i < nDvAct; ++i)
	{
		for (int j = 0; j < nConstrAct; ++j)
		{
			this->sensDv(this->subDvAct[i],this->subConstrAct[i]) = 
				sensRes(i, j);
		}
	}
}