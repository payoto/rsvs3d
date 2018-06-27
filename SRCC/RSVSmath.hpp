#ifndef RSVSMATH_H_INCLUDED 
#define RSVSMATH_H_INCLUDED 


//=================================
// forward declared dependencies
// 		class foo; //when you only need a pointer not the actual object
// 		and to avoid circular dependencies


//=================================
// included dependencies
#include <vector> 
#include <cmath> 
#include "vectorarray.hpp" 
#include "snakevel.hpp" 
#include "RSVSmath_automatic.hpp"
using namespace std; 


//==================================
// Code
// NOTE: function in a class definition are IMPLICITELY INLINED 
//       ie replaced by their code at compile time
	

class TriFunc {
protected:
	vector<double>* p0=NULL;
	vector<double>* p1=NULL;
	vector<double>* p2=NULL;

	double fun;
	ArrayVec<double> jac;
	ArrayVec<double> hes;

	bool isReady;
	bool isCalc;
	int nTarg; // target length of vectors

	bool MakeValidField(vector<double>* TriFunc::*mp);

public:
	// Check validity
	bool CheckValid();
	bool MakeValid();
	void PreCalc();
	// Build a valid object
	void assign(vector<double> &in0,vector<double> &in1,vector<double> &in2);
	void assign(int pRepI,vector<double> &pRep);
	virtual void Calc() = 0; // Virtual function that calculates 

	TriFunc(){
		fun=0;
		nTarg=3;
		jac.assign(1,3*nTarg,fun);
		hes.assign(3*nTarg,3*nTarg,fun);
		isReady=false;
		isCalc=false;
	}
	TriFunc(int nTarg){
		fun=0;
		jac.assign(1,3*nTarg,fun);
		hes.assign(3*nTarg,3*nTarg,fun);
		isReady=false;
		isCalc=false;
	}

};

class Volume : public TriFunc {
	using TriFunc::TriFunc;
	using TriFunc::PreCalc;
	using TriFunc::p0;
	using TriFunc::p1;
	using TriFunc::p2;
	using TriFunc::fun;
	using TriFunc::jac;
	using TriFunc::hes;

public:
	void Calc();
};

class Area : public TriFunc {
	using TriFunc::TriFunc;
	using TriFunc::PreCalc;
	using TriFunc::p0;
	using TriFunc::p1;
	using TriFunc::p2;
	using TriFunc::fun;
	using TriFunc::jac;
	using TriFunc::hes;

public:
	void Calc();

};

class LengthEdge : public TriFunc {
	using TriFunc::TriFunc;
	using TriFunc::PreCalc;
	using TriFunc::p0;
	using TriFunc::p1;
	using TriFunc::p2;
	using TriFunc::fun;
	using TriFunc::jac;
	using TriFunc::hes;

public:
	void Calc();
};

class SurfCentroid {


};

#endif