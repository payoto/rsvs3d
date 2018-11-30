#include <iostream>
#include <vector>
#include <fstream>
#include <Eigen>
#include "vectorarray.hpp"
// #include "matrixtools.hpp"

using namespace std;
using namespace Eigen; 

void ArrayVec2MatrixXd(const ArrayVec<double> &arrayIn, MatrixXd &matOut){

	int nR=0,nC=0;
	arrayIn.size(nR,nC);
	matOut.setZero(nR,nC);

	for(int i=0; i< nR; ++i){
		for(int j=0; j< nC; ++j){
			matOut(i,j)=arrayIn[i][j];
		}
	}

}

void VecBy3DimArray(const MatrixXd &vec, const MatrixXd &arr3dim, MatrixXd &retArray){
	/*
	3 Dimensional array is expected as 
	 					xc 				 					yc 				 				zc 				
	  x0 x1  xn y0 y1  yn  z0 z1  zn	x0 x1  xn y0 y1  yn  z0 z1  zn	    x0 x1  xn y0 y1  yn  z0 z1  zn
 x0 [                               ]  [                               ]  [                               ] 
 x1 [                               ]  [                               ]  [                               ] 
 xn [                               ]  [                               ]  [                               ] 
 y0 [                               ]  [                               ]  [                               ] 
 y1 [                               ]  [                               ]  [                               ] 
 yn [                               ]  [                               ]  [                               ] 
 z0 [                               ]  [                               ]  [                               ] 
 z1 [                               ]  [                               ]  [                               ] 
 zn [                               ]  [                               ]  [                               ] 
	vector expected as a row vector
	[xc yc zc]

	will sum the three Hessians according to the values in the vector
	*/

	int nRow,nCol,nVec,nColFin;
	int ii,jj,kk;
	nRow=arr3dim.rows();
	nCol=arr3dim.cols();
	nVec=vec.rows();
	nColFin=nCol/nVec;

	#ifdef SAFE_ALGO
	// size checks
	if ((nVec*nColFin)!=nCol){
		throw invalid_argument("Sizes do not match in 3D array collapse");
	}
	#endif

	retArray.setZero(nRow,nColFin);
	// needs to add checks for matching sizes
	for (ii=0;ii<nRow;ii++){
		for (jj=0;jj<nColFin;jj++){
			for (kk=0;kk<nVec;kk++){
				retArray(ii,jj)+=arr3dim(ii,jj*nVec+kk)*vec(0,kk);
			}	
		}	
	}

}


void Deriv1stChainScalar(const MatrixXd &dSdc,const MatrixXd &dcdd, MatrixXd &dSdd){
	dSdd=dSdc*dcdd;
}

void Deriv2ndChainScalar(const MatrixXd &dSdc,const MatrixXd &dcdd,
	const MatrixXd &HSc,const MatrixXd &Hcd,MatrixXd &HSd){
	
	VecBy3DimArray(dSdc, Hcd, HSd);
	HSd = HSd + (dcdd.transpose()*HSc*dcdd);

}

void PrintMatrix(const MatrixXd &mat){
	int ii,jj, ni, nj;

	ni=mat.rows();
	nj=mat.cols();
	for(ii=0;ii<ni;++ii){
		for(jj=0;jj<nj;++jj){
			cout << mat(ii,jj) << " ";
		}
		cout << endl;
	}
}
void PrintMatrixFile(const MatrixXd &mat, const char * name){
	int ii,jj, ni, nj;
	ofstream myfile;
	ni=mat.rows();
	nj=mat.cols();
	myfile.open(name, ios::app);
	for(ii=0;ii<ni;++ii){
		for(jj=0;jj<nj;++jj){
			myfile  << mat(ii,jj) << " ";
		}
		myfile  << endl;
	}
	myfile  << endl;
	myfile.close();

}

void PrintMatrix(const VectorXd &mat){
	int ii, ni;

	ni=mat.size();
	for(ii=0;ii<ni;++ii){
		cout << mat[ii] << " ";
		cout << endl;
	}
}
void PrintMatrix(const RowVectorXd &mat){
	int ii, ni;

	ni=mat.size();
	for(ii=0;ii<ni;++ii){
		cout << mat[ii] << " ";
		
	}
	cout << endl;
}

void StreamStatistics(const VectorXd &&vec, ofstream &out, const string &&sep){
	/*
	Uses a rvalue reference to allow operations to be passed 
	directly in the statistics
	*/
	out << vec.norm() << sep;
	out << vec.mean() << sep;
	out << sqrt((vec.array() - vec.mean()).square().sum()/(vec.size()-1))<< sep;
	out << vec.maxCoeff() << sep;
	out << vec.minCoeff() << sep;
	out << endl;
}


void StreamOutVector(const VectorXd &&vec, ofstream &out, const string &&sep){
	/*
	Uses a rvalue reference to allow operations to be passed 
	directly in the statistics
	*/
	int i, n;
	n = vec.size();
	for(i = 0; i<n; ++i){
		out << vec[i] << sep;
	}
	out << endl;
}

