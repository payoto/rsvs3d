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
#include "RSVSmath_automatic.hpp"
using namespace std; 


//==================================
// Code
// NOTE: function in a class definition are IMPLICITELY INLINED 
//       ie replaced by their code at compile time
	

class TriFunc {
protected:
	vector<double> const * p0=NULL;
	vector<double> const * p1=NULL;
	vector<double> const * p2=NULL;

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


class CoordFunc {
protected:
	vector<vector<double> const *> coords;

	double fun;
	ArrayVec<double> funA;
	ArrayVec<double> jac;
	ArrayVec<double> hes;

	bool isReady;
	bool isCalc;

	int nDim; // target length of vectors
	int nCoord; // target length of vectors
	int nFun;

	bool MakeValidField(vector<double> const* mp);
	void InitialiseArrays();
public:
	// Check validity
	bool CheckValid();
	bool MakeValid();
	void PreCalc();
	// Build a valid object
	void assign(vector<vector<double> const*> &coords);
	void assign(int pRepI,vector<double> &pRep);

	virtual void Calc() = 0; // Virtual function that calculates the function

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
	void ResetDim(int nDim){InitialiseArrays();}
	void ResetNCoord(int nCoord){InitialiseArrays();}
	void ResetNFun(int nFun){InitialiseArrays();}
	CoordFunc(){
		nDim=3;
		nCoord=3;
		nFun=1;
		InitialiseArrays();
	}
	CoordFunc(int nDim){
		nCoord=3;
		nFun=1;
		InitialiseArrays();
	}
	CoordFunc(int nDim,int nCoord){
		nFun=1;
		InitialiseArrays();
	}
	CoordFunc(int nDim,int nCoord,int nFun){
		InitialiseArrays();
	}
#pragma GCC diagnostic pop
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

class LengthEdge : public CoordFunc {

	using CoordFunc::PreCalc;
	using CoordFunc::coords;
	using CoordFunc::fun;
	using CoordFunc::jac;
	using CoordFunc::hes;

public:
	void Calc();
	LengthEdge() : CoordFunc(3,2){}
};

class SurfCentroid : public CoordFunc {
protected:
	//using CoordFunc::PreCalc;
	using CoordFunc::coords;
	using CoordFunc::fun;
	using CoordFunc::jac;
	using CoordFunc::hes;
	using CoordFunc::nCoord;

	vector<double> centroid;
	double edgeLength=0.0;
public:
	
	void Calc();
	void assigncentroid(const vector<double> &vecin);
	SurfCentroid() : CoordFunc(3,4,3){centroid.assign(nDim,0);};
	SurfCentroid(int nCoord) : CoordFunc(3,nCoord,3){centroid.assign(nDim,0);};


};

#endif