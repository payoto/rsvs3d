#ifndef MATRIXTOOLS_H_INCLUDED 
#define MATRIXTOOLS_H_INCLUDED 


//=================================
// forward declared dependencies
// 		class foo; //when you only need a pointer not the actual object
// 		and to avoid circular dependencies


//=================================
// included dependencies

#include <vector> 
#include <Eigen>
#include <iostream>
#include <fstream>
#include "vectorarray.hpp"

using namespace std; 
using namespace Eigen; 

template <class T> void PrintMatrixFile(const vector<T> &mat, const char * name);
void Deriv1stChainScalar(const MatrixXd &dSdc,const MatrixXd &dcdd, MatrixXd &dSdd);
void Deriv2ndChainScalar(const MatrixXd &dSdc,const MatrixXd &dcdd,const MatrixXd &HSc,const MatrixXd &Hcd,MatrixXd &HSd);
void VecBy3DimArray(const MatrixXd &vec, const MatrixXd &arr3dim, MatrixXd &retArray);
void ArrayVec2MatrixXd(const ArrayVec<double> &arrayIn, MatrixXd &matOut);
void PrintMatrix(const MatrixXd mat);
void PrintMatrixFile(const MatrixXd mat, const char * name);
void PrintMatrix(const RowVectorXd mat);
void PrintMatrix(const VectorXd mat);

//==================================
// Code
// NOTE: function in a class definition are IMPLICITELY INLINED 
//       ie replaced by their code at compile time

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

#endif