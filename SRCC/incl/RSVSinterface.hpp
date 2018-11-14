#ifndef RSVSINTERFACE_H_INCLUDED 
#define RSVSINTERFACE_H_INCLUDED 


//=================================
// forward declared dependencies
// 		class foo; //when you only need a pointer not the actual object
// 		and to avoid circular dependencies


//=================================
// included dependencies
#include <iostream>
#include <vector> 
#include "vectorarray.hpp" 
#include "RSVSmath.hpp"
#include "mesh.hpp"
#include "snake.hpp"
#include "snakevel.hpp"
#include <Eigen>

using namespace std; 
using namespace Eigen; 

class SQPcalc {
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

	void ReturnConstrToMesh(triangulation &triRSVS) const ;
	void ReturnConstrToMesh(mesh &meshin, double volu::*mp=&volu::volume) const ;
	void CheckAndCompute();
	void ComputeSQPstep(
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
	void ConvergenceLog(ostream out);
};

//==================================
// Code
// NOTE: function in a class definition are IMPLICITELY INLINED 
//       ie replaced by their code at compile time
	

#endif