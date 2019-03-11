#ifndef RSVSINTERFACE_H_INCLUDED 
#define RSVSINTERFACE_H_INCLUDED 


//=================================
// forward declared dependencies
// 		class foo; //when you only need a pointer not the actual object
// 		and to avoid circular dependencies


//=================================
// included dependencies
#include <iostream>
#include <fstream>
#include <vector> 
#include "vectorarray.hpp" 
#include "RSVSmath.hpp"
#include "mesh.hpp"
#include "snake.hpp"
#include "triangulate.hpp"
#include <Eigen>

using namespace std; 
using namespace Eigen; 

class RSVScalc {
protected:
	int nDv=0;
	int nConstr=0;
	int falseaccess=0;
	bool returnDeriv=true;
public :
	MatrixXd dConstr,HConstr, HObj, HLag;
	RowVectorXd dObj;
	VectorXd constr, lagMult, deltaDV, constrTarg;
	double obj=0.0;
	double limLag = INFINITY;

	vector<bool> isConstrAct, isDvAct;
	vector<int> subConstrAct, subDvAct;
	HashedVector<int, int> dvMap;
	HashedMap<int,int,int> constrMap; // maps snakemesh volu onto constr
	vector<pair<int,int>> constrList; // keeps pairs with parentindex and voluindex
	MatrixXd dvCallConstr;

	void BuildMathArrays(int nDv, int nConstr);
	void BuildConstrMap(const triangulation &triangleRSVS);
	void BuildConstrMap(const mesh &meshin);
	void BuildDVMap(const vector<int> &vecin);
	bool SnakDVcond(const triangulation &triRSVS, int ii);
	void PrepTriangulationCalc(const triangulation &triRSVS);
	// Calculate derivatives wrapper
	void CalculateMesh(mesh &meshin);
	void CalculateTriangulation(const triangulation &triRSVS, int derivMethod=0);
	// Calculate derivatives
	void CalcTriangle(
		const triangle& triIn, const triangulation &triRSVS,
		bool isObj=true, bool isConstr=true, bool isDeriv=true
		);
	void CalcTriangleFD(
		const triangle& triIn, const triangulation &triRSVS,
		bool isObj=true, bool isConstr=true, bool isDeriv=true
		);
	void CalcTriangleDirectVolume(
		const triangle& triIn, const triangulation &triRSVS,
		bool isObj=true, bool isConstr=true, bool isDeriv=true
		);
	void CalcTriangleEdgeLength(
		const triangle& triIn, const triangulation &triRSVS,
		bool isObj=true, bool isConstr=true, bool isDeriv=true);

	void ReturnConstrToMesh(triangulation &triRSVS) const ;
	void ReturnConstrToMesh(mesh &meshin, double volu::*mp=&volu::volume) const ;
	void CheckAndCompute(int calcMethod=0);
	void ComputeSQPstep(int calcMethod,
		MatrixXd &dConstrAct,
		RowVectorXd &dObjAct,
		VectorXd &constrAct,
		VectorXd &lagMultAct
		);
	bool PrepareMatricesForSQP(
		MatrixXd &dConstrAct,
		MatrixXd &HConstrAct, 
		MatrixXd &HObjAct,
		RowVectorXd &dObjAct,
		VectorXd &constrAct,
		VectorXd &lagMultAct
		);
	void ReturnVelocities(triangulation &triRSVS);
	int numConstr(){return(this->nConstr);}
	// Output functions
	void Print2Screen(int outType=0)const;
	void ConvergenceLog(ofstream &out, int loglvl=3) const;
};

//==================================
// Code
// NOTE: function in a class definition are IMPLICITELY INLINED 
//       ie replaced by their code at compile time
	

void ResizeLagrangianMultiplier(const RSVScalc &calcobj, 
	VectorXd &lagMultAct, 
	bool &isLarge, bool &isNan);
template<class T>
bool SQPstep(const RSVScalc &calcobj,
	const MatrixXd &dConstrAct, const RowVectorXd &dObjAct,
	const VectorXd &constrAct, VectorXd &lagMultAct,
	VectorXd &deltaDVAct, bool &isNan, bool &isLarge, bool attemptConstrOnly=true);
template<template<typename> class T>
bool SQPstep(const RSVScalc &calcobj,
	const MatrixXd &dConstrAct, const RowVectorXd &dObjAct,
	const VectorXd &constrAct, VectorXd &lagMultAct,
	VectorXd &deltaDVAct, bool &isNan, bool &isLarge, bool attemptConstrOnly=true);


// Code needs to be included as it is a templated functions

 
template<template<typename> class T>
bool SQPstep(const RSVScalc &calcobj,
	const MatrixXd &dConstrAct, const RowVectorXd &dObjAct,
	const VectorXd &constrAct, VectorXd &lagMultAct,
	VectorXd &deltaDVAct, bool &isNan, bool &isLarge, bool attemptConstrOnly){
	/*
	This template cannot be deduced and needs the developer to
	pass the required solver template class when it is called.

	This accepts any single parameter template for instantiation

	Instantiation options:
	Eigen::HouseholderQR
	Eigen::ColPivHouseholderQR
	Eigen::LLT<MatrixXd> (*) <- needs a full type to be defined (see below)
	Eigen::PartialPivLU

	For stability info
	https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html
	https://eigen.tuxfamily.org/dox/group__DenseDecompositionBenchmark.html

	*see : template<class T> void SQPstep
	*/

	return(SQPstep<T<MatrixXd>>(calcobj, dConstrAct, dObjAct,
				constrAct, lagMultAct,
				deltaDVAct, isNan, isLarge,attemptConstrOnly));

}

template<class T>
bool SQPstep(const RSVScalc &calcobj,
	const MatrixXd &dConstrAct, const RowVectorXd &dObjAct,
	const VectorXd &constrAct, VectorXd &lagMultAct,
	VectorXd &deltaDVAct, bool &isNan, bool &isLarge, bool attemptConstrOnly){
	/*
	This template cannot be deduced and needs the developer to
	pass the required solver class when it is called.

	Instantiation options:
	Eigen::HouseholderQR<MatrixXd>
	Eigen::ColPivHouseholderQR<MatrixXd>
	Eigen::LLT<MatrixXd>
	Eigen::PartialPivLU<MatrixXd>

	For stability info
	https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html
	https://eigen.tuxfamily.org/dox/group__DenseDecompositionBenchmark.html
	*/

	MatrixXd temp1, temp2;
	T HLagSystem(calcobj.HLag);

	temp1 = HLagSystem.solve(dConstrAct.transpose());
	temp2 = HLagSystem.solve(dObjAct.transpose());
	
	T LagMultSystem(dConstrAct*(temp1));

	lagMultAct = LagMultSystem.solve(
			constrAct - (dConstrAct*(temp2))
		);

	ResizeLagrangianMultiplier(calcobj, lagMultAct, isLarge, isNan);
	isLarge=false;
	if(!attemptConstrOnly && (isLarge || isNan)){
		return(isLarge || isNan);
	}
	if(isLarge || isNan) {
		// Use a least squared solver if only using the constraint
		std::cout << "(constrmov) " ;
	 	deltaDVAct = -dConstrAct.bdcSvd(ComputeThinU | ComputeThinV).solve(constrAct);
	} else {
		deltaDVAct = - (HLagSystem.solve(dObjAct.transpose() 
						+ dConstrAct.transpose()*lagMultAct));
	}
	// cout << __PRETTY_FUNCTION__<< endl;
	return(isLarge || isNan);
}

#endif