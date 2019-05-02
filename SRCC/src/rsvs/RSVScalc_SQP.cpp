#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath> 
#include <Eigen>

// #include "snake.hpp"
// #include "triangulate.hpp"
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
	MatrixXd_sparse &dConstrAct_sparse,
	MatrixXd_sparse &HConstrAct_sparse,
	MatrixXd_sparse &HObjAct_sparse,
	RowVectorXd &dObjAct,
	VectorXd &constrAct,
	VectorXd &lagMultAct
	){

	int ii, jj, nDvAct, nConstrAct;


	this->subDvAct.reserve(nDv);
	this->subDvAct.clear();
	nDvAct=0;
	for(ii=0; ii<nDv; ++ii){
		nDvAct += int(isDvAct.at(ii));
		if(isDvAct.at(ii)){
			this->subDvAct.push_back(ii);
		}
	}

	this->subConstrAct.reserve(nConstr);
	this->subConstrAct.clear();
	nConstrAct=0;
	for(ii=0; ii<nConstr; ++ii){
		nConstrAct += int(isConstrAct.at(ii));
		if(isConstrAct.at(ii)){
			this->subConstrAct.push_back(ii);
		}
	}

	dObjAct.setZero(nDvAct);
	for(jj=0; jj<nDvAct;++jj){
		dObjAct[jj]=dObj[this->subDvAct[jj]];
	}
	constrAct.setZero(nConstrAct);
	for(jj=0; jj<nConstrAct;++jj){
		constrAct[jj]=(constr[this->subConstrAct[jj]] 
			- constrTarg[this->subConstrAct[jj]]);
	}
	lagMultAct.setZero(nConstrAct);

	this->subConstrAct.GenerateHash();
	this->subDvAct.GenerateHash();

	if(this->UseFullMath()){
		// Full maths
		this->PrepareMatricesForSQPFull(dConstrAct, HConstrAct, HObjAct);
	} else {
		this->PrepareMatricesForSQPSparse(dConstrAct_sparse, HConstrAct_sparse,
			HObjAct_sparse);
		
	}

	return(nDvAct>0);
}

void RSVScalc::PrepareMatricesForSQPFull(
	Eigen::MatrixXd &dConstrAct,
	Eigen::MatrixXd &HConstrAct, 
	Eigen::MatrixXd &HObjAct){

	int ii, jj;
	int nConstrAct = this->subConstrAct.size();
	int nDvAct = this->subDvAct.size();

	dConstrAct.setZero(nConstrAct,nDvAct);
	for(ii=0; ii<nConstrAct;++ii){
		for(jj=0; jj<nDvAct;++jj){
			dConstrAct(ii,jj)=dConstr(this->subConstrAct[ii],this->subDvAct[jj])
				// /constrTarg[this->subConstrAct[ii]]
				;
		}
	}
	HConstrAct.setZero(nDvAct,nDvAct); 
	for(ii=0; ii<nDvAct;++ii){
		for(jj=0; jj<nDvAct;++jj){
			HConstrAct(ii,jj)=HConstr(this->subDvAct[ii],this->subDvAct[jj]);
		}
	}
	HObjAct.setZero(nDvAct,nDvAct);  
	for(ii=0; ii<nDvAct;++ii){
		for(jj=0; jj<nDvAct;++jj){
			HObjAct(ii,jj)=HObj(this->subDvAct[ii],this->subDvAct[jj]);
		}
	}

	this->dLag = this->dObj + this->lagMult.transpose()*this->dConstr;
	this->HLag = HObjAct+HConstrAct;

}
void RSVScalc::PrepareMatricesForSQPSparse(
	MatrixXd_sparse &dConstrAct_sparse,
	MatrixXd_sparse &HConstrAct_sparse,
	MatrixXd_sparse &HObjAct_sparse){

	int ii, jj, count;
	int nConstrAct = this->subConstrAct.size();
	int nDvAct = this->subDvAct.size();
	// Sparse maths 
	dConstrAct_sparse.resize(nConstrAct,nDvAct);
	count = this->dConstr_sparse.nonZeros();
	dConstrAct_sparse.reserve(count);
	for (int k=0; k<count; ++k){
		auto it = dConstr_sparse[k];
		if(isConstrAct[it.row()] && isDvAct[it.col()]){
			ii = subConstrAct.find(it.row());
			jj = subDvAct.find(it.col());
			if(ii==rsvs3d::constants::__notfound 
				|| jj==rsvs3d::constants::__notfound){
				std::cerr << std::endl << "rows " << it.row()
					<< ", cols " << it.col()
					<< ", value " << it.value() << std::endl;
				RSVS3D_ERROR_LOGIC("It should not be possible to get here.");
			}

			dConstrAct_sparse.insert(ii,jj) = it.value();
		}
	}

	HConstrAct_sparse.resize(nDvAct,nDvAct); 
	count = this->HConstr_sparse.nonZeros();
	HConstrAct_sparse.reserve(count);
	for (int k=0; k<count; ++k){
		auto it = HConstr_sparse[k];
		if(isDvAct[it.row()] && isDvAct[it.col()]){
			ii = subDvAct.find(it.row());
			jj = subDvAct.find(it.col());
			HConstrAct_sparse.insert(ii,jj) = it.value();
		}
	}

	HObjAct_sparse.resize(nDvAct,nDvAct);  
	count = this->HObj_sparse.nonZeros();
	HObjAct_sparse.reserve(count);
	for (int k=0; k<count; ++k){
		auto it = HObj_sparse[k];
		if(isDvAct[it.row()] && isDvAct[it.col()]){
			ii = subDvAct.find(it.row());
			jj = subDvAct.find(it.col());
			HObjAct_sparse.insert(ii,jj) = it.value();
		}
	}
	this->dConstr_sparse.SetEqual(this->dConstr);
	this->dLag = this->dObj + this->lagMult.transpose()*this->dConstr;
	this->HLag_sparse = HObjAct_sparse+HConstrAct_sparse;

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


	// Find active design variables and active constraint
	nDvAct = this->subDvAct.size();
	nConstrAct  = this->subConstrAct.size();


	/// sensitivity requires : [H_d L , J_d constr; J_d constr ^T, 0]^-1
	/// [Dd^TD_v L (=0s) ; D_V constr (=eye(lagMult))]
 	
	sensInv.resize(nDvAct+nConstrAct, nDvAct+nConstrAct);
	sensMult.resize(nDvAct+nConstrAct, nConstrAct);
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


bool RSVScalc::PrepareMatricesForSQPSensitivity(
	const MatrixXd_sparse &dConstrAct,
	const MatrixXd_sparse &HConstrAct, 
	MatrixXd_sparse &HObjAct,
	Eigen::MatrixXd &sensMult,
	MatrixXd_sparse &sensInv,
	Eigen::MatrixXd &sensRes) const {

	int nDvAct, nConstrAct;

	// Find active design variables and active constraint
	nDvAct = this->subDvAct.size();
	nConstrAct  = this->subConstrAct.size();


	/// sensitivity requires : [H_d L , J_d constr; J_d constr ^T, 0]^-1
	/// [Dd^TD_v L (=0s) ; D_V constr (=eye(lagMult))]
 	
	sensInv.resize(nDvAct+nConstrAct, nDvAct+nConstrAct);
	sensMult.resize(nDvAct+nConstrAct, nConstrAct);
	sensRes.resize(nDvAct+nConstrAct, nConstrAct);

	// Assign sensInv
	HObjAct = HConstrAct+HObjAct;
	MatrixXd_sparse& HlagAct = HObjAct;
	sensInv.reserve(HlagAct.nonZeros() + 2 * dConstrAct.nonZeros());
	for (int k=0; k<HlagAct.outerSize(); ++k){
		for (SparseMatrix<double>::InnerIterator it(HlagAct,k); it; ++it)
		{
			sensInv.insert(it.row(),it.col()) = it.value();
		}
	}
	for (int k=0; k<dConstrAct.outerSize(); ++k){
		for (SparseMatrix<double>::InnerIterator it(dConstrAct,k); it; ++it)
		{
			sensInv.insert(it.col(), it.row()+nDvAct) = it.value();
			sensInv.insert(it.row()+nDvAct, it.col()) = it.value();
		}
	}

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
	MatrixXd_sparse dConstrAct_sparse,HConstrAct_sparse, HObjAct_sparse;
	RowVectorXd dObjAct;
	VectorXd constrAct, lagMultAct, deltaDVAct;

	computeFlag = this->PrepareMatricesForSQP(
		dConstrAct,HConstrAct, HObjAct,
		dConstrAct_sparse,HConstrAct_sparse, HObjAct_sparse,
		dObjAct,constrAct,lagMultAct
		);
	// computeFlag=true;
	if (computeFlag){
		if(this->UseFullMath()){
			this->ComputeSQPstep(calcMethod,dConstrAct,dObjAct,
				constrAct,lagMultAct);
		} else {
			this->ComputeSQPstep(calcMethod,dConstrAct_sparse,dObjAct,
				constrAct,lagMultAct);
		}
	} else {
		for (ii=0; ii<nDv; ++ii){
			deltaDV[ii]=0.0;
		}
	}
	if(sensCalc){ 
		if(this->UseFullMath()){
			MatrixXd sensMult, sensInv, sensRes;
			bool computeSens = this->PrepareMatricesForSQPSensitivity(dConstrAct,
				HConstrAct, HObjAct, sensMult, sensInv, sensRes);
			if (computeSens)
			{
				this->ComputeSQPsens(calcMethod, sensMult, sensInv, sensRes);
			}
		} else {
			MatrixXd_sparse sensInv;
			MatrixXd sensRes, sensMult;
			bool computeSens = this->PrepareMatricesForSQPSensitivity(
				dConstrAct_sparse, HConstrAct_sparse, HObjAct_sparse,
				sensMult, sensInv, sensRes);
			if (computeSens)
			{
				sensInv.makeCompressed();
				this->ComputeSQPsens(calcMethod, sensMult, sensInv, sensRes);
			}
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

	ni = this->subDvAct.size();
	for (ii=0; ii<ni; ++ii){
		this->deltaDV[subDvAct[ii]]=deltaDVAct[ii];		
	}
	
	ni = this->subConstrAct.size();
	this->lagMult.setZero(nConstr);
	if(isLarge){
		for (ii=0; ii<nConstr; ++ii){
			lagMult[subConstrAct[ii]]=0;
		}

	} else {
		for (ii=0; ii<ni; ++ii){
			this->lagMult[this->subConstrAct[ii]]=lagMultAct[ii];
			isNan = isNan || isnan(lagMultAct[ii]);
		}
	}
	if (false){
		Print2Screen(3);
		DisplayVector(isConstrAct);
		DisplayVector(isDvAct);
	}
}


void RSVScalc::ComputeSQPstep(
	int calcMethod,
	MatrixXd_sparse &dConstrAct,
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
	this->HLag_sparse.makeCompressed();
	dConstrAct.makeCompressed();
	while(rerunDefault){
		switch(calcMethod){
			default:
				SQPstep_sparse<Eigen::SparseLU<MatrixXd_sparse>>
					(*this, dConstrAct, dObjAct, constrAct, lagMultAct, 
						deltaDVAct, isNan, isLarge,true);
				rerunDefault=false;
			break;
		}
		if(attemptConstrOnly){
			rerunDefault=false;
		}
		if(rerunDefault){
			calcMethod=0;
		}
	}

	ni = this->subDvAct.size();
	for (ii=0; ii<ni; ++ii){
		this->deltaDV[subDvAct[ii]]=deltaDVAct[ii];		
	}
	
	ni = this->subConstrAct.size();
	this->lagMult.setZero(nConstr);
	if(isLarge){
		for (ii=0; ii<nConstr; ++ii){
			lagMult[subConstrAct[ii]]=0;
		}

	} else {
		for (ii=0; ii<ni; ++ii){
			this->lagMult[this->subConstrAct[ii]]=lagMultAct[ii];
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

	this->sensDv.setZero(this->nDv, this->nConstr);
	int nConstrAct, nDvAct;
	nConstrAct = this->subConstrAct.size();
	nDvAct = this->subDvAct.size();
	for (int i = 0; i < nDvAct; ++i)
	{
		for (int j = 0; j < nConstrAct; ++j)
		{
			this->sensDv(this->subDvAct[i],this->subConstrAct[j]) = 
				sensRes(i, j);
		}
	}
}

bool SQPsens_sparse2(
	const Eigen::MatrixXd &sensMult,
	const MatrixXd_sparse &sensInv,
	Eigen::MatrixXd &sensRes){

	Eigen::SparseLU<MatrixXd_sparse> HLagSystem;
	HLagSystem.compute(sensInv);
	if(HLagSystem.info()!=Eigen::Success) {
		// decomposition failed
		RSVS3D_ERROR_NOTHROW("Failed to decompose sensitivity system.");
		return(false);
	}

	sensRes = - (HLagSystem.solve(sensMult));

	return(true);
}

void RSVScalc::ComputeSQPsens(
	int calcMethod,
	Eigen::MatrixXd &sensMult,
	MatrixXd_sparse &sensInv,
	Eigen::MatrixXd &sensRes){


	switch(calcMethod){
		default:
			SQPsens_sparse<Eigen::SparseLU<MatrixXd_sparse>>
				(sensMult, sensInv, sensRes);
		break;
	}

	this->sensDv.setZero(this->nDv, this->nConstr);
	int nConstrAct, nDvAct;
	nConstrAct = this->subConstrAct.size();
	nDvAct = this->subDvAct.size();
	for (int i = 0; i < nDvAct; ++i)
	{
		for (int j = 0; j < nConstrAct; ++j)
		{
			this->sensDv(this->subDvAct[i],this->subConstrAct[j]) = 
				sensRes(i,j);
		}
	}
}