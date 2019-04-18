/**
 * Tools to support conversion, display and derivatives of Eigen matrices.
 * @file
 */

#ifndef MATRIXTOOLS_H_INCLUDED 
#define MATRIXTOOLS_H_INCLUDED 

//=================================
// forward declared dependencies
// 		class foo; //when you only need a pointer not the actual object
// 		and to avoid circular dependencies



//=================================
// included dependencies

#include <vector> 
#include <iostream>
#include <fstream>
#include <string>
#include <Eigen>
#include "vectorarray.hpp"


template <class T> void PrintMatrixFile(const std::vector<T> &mat,
	const char * name);
void Deriv1stChainScalar(const Eigen::MatrixXd &dSdc,
	const Eigen::MatrixXd &dcdd, Eigen::MatrixXd &dSdd);
void Deriv2ndChainScalar(const Eigen::MatrixXd &dSdc,
	const Eigen::MatrixXd &dcdd, 
	const Eigen::MatrixXd &HSc,
	const Eigen::MatrixXd &Hcd,
	Eigen::MatrixXd &HSd);
void VecBy3DimArray(const Eigen::MatrixXd &vec, 
	const Eigen::MatrixXd &arr3dim, Eigen::MatrixXd &retArray);
void ArrayVec2MatrixXd(const ArrayVec<double> &arrayIn,
	Eigen::MatrixXd &matOut);
void PrintMatrix(const Eigen::MatrixXd &mat);
void PrintMatrixFile(const Eigen::MatrixXd &mat,
	const char * name);
void PrintMatrixFile(const Eigen::MatrixXd &mat, ofstream& myfile);
void PrintMatrix(const Eigen::RowVectorXd &mat);
void PrintMatrix(const Eigen::VectorXd &mat);
double StreamStatistics(const Eigen::VectorXd &&vec,
	std::ofstream &out, 
	const std::string &&sep=std::string(", "));
void StreamOutVector(const Eigen::VectorXd &&vec, std::ofstream &out, 
	const std::string &&sep=std::string(", "));

int Test_Matrix3D();


//==================================
// Code
// NOTE: function in a class definition are IMPLICITELY INLINED 
//       ie replaced by their code at compile time

template <class T> void PrintMatrixFile(const std::vector<T> &mat, 
	const char * name){
	int ii,ni;
	std::ofstream myfile;
	ni=mat.size();

	myfile.open(name, std::ios::app);
	myfile.precision(16);
	myfile << std::scientific;
	for(ii=0;ii<ni;++ii){
		myfile  << mat[ii] << " ";
	}
	myfile  << std::endl;
	myfile.close();

}

#endif