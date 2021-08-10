#include "matrixtools.hpp"

#include <Eigen>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "vectorarray.hpp"
#include "warning.hpp"
// #include "matrixtools.hpp"

using namespace std;
using namespace Eigen;

void ArrayVec2MatrixXd(const ArrayVec<double> &arrayIn, MatrixXd &matOut)
{
    int nR = 0, nC = 0;
    arrayIn.size(nR, nC);
    matOut.setZero(nR, nC);

    for (int i = 0; i < nR; ++i)
    {
        for (int j = 0; j < nC; ++j)
        {
            matOut(i, j) = arrayIn[i][j];
        }
    }
}

void VecBy3DimArray(const MatrixXd &vec, const MatrixXd &arr3dim, MatrixXd &retArray)
{
    /*
    3 Dimensional array is expected as
                                            xc
  yc 				 				zc x0 x1  xn y0
  y1  yn  z0 z1  zn	x0 x1  xn y0 y1  yn  z0 z1  zn	    x0 x1  xn y0 y1  yn
  z0 z1  zn x0 [                               ]  [ ]  [ ] x1 [ ]  [ ]  [ ] xn [
  ]  [                               ]  [                               ] y0 [ ]
  [                               ]  [                               ] y1 [ ]  [
  ]  [                               ] yn [                               ]  [ ]
  [                               ] z0 [                               ]  [ ]  [ ]
  z1 [                               ]  [                               ]  [ ] zn
  [                               ]  [                               ]  [ ] vector
  expected as a row vector [xc yc zc]

    will sum the three Hessians according to the values in the vector
    */

    int nRow, nCol, nVec, nColFin;
    int ii, jj, kk;
    nRow = arr3dim.rows();
    nCol = arr3dim.cols();
    nVec = vec.cols();
    nColFin = nCol / nVec;

#ifdef SAFE_ALGO
    // size checks
    if ((nVec * nColFin) != nCol)
    {
        std::cerr << std::endl << "nVec (" << nVec << ") ; nColFin (" << nColFin << ") ; nCol (" << nCol << ")";
        RSVS3D_ERROR_ARGUMENT("Sizes do not match in 3D array collapse");
    }
#endif

    retArray.setZero(nRow, nColFin);
    // needs to add checks for matching sizes
    for (ii = 0; ii < nRow; ii++)
    {
        for (jj = 0; jj < nColFin; jj++)
        {
            for (kk = 0; kk < nVec; kk++)
            {
                retArray(ii, jj) += arr3dim(ii, jj + kk * nColFin) * vec(0, kk);
            }
        }
    }
}

void Deriv1stChainScalar(const MatrixXd &dSdc, const MatrixXd &dcdd, MatrixXd &dSdd)
{
    dSdd = dSdc * dcdd;
}

void Deriv2ndChainScalar(const MatrixXd &dSdc, const MatrixXd &dcdd, const MatrixXd &HSc, const MatrixXd &Hcd,
                         MatrixXd &HSd)
{
    VecBy3DimArray(dSdc, Hcd, HSd);
#ifdef SAFE_ALGO
    bool errFlag = false;
    errFlag |= dcdd.rows() != HSc.cols();
    errFlag |= dcdd.cols() != HSd.cols();
    errFlag |= HSd.rows() != HSd.cols();
    if (errFlag)
    {
        std::cerr << std::endl << "dSdc: [" << dSdc.rows() << ", " << dSdc.cols() << "]" << std::endl;
        PrintMatrix(dSdc);
        std::cerr << std::endl << "dcdd: [" << dcdd.rows() << ", " << dcdd.cols() << "]" << std::endl;
        PrintMatrix(dcdd);
        std::cerr << std::endl << "HSc: [" << HSc.rows() << ", " << HSc.cols() << "]" << std::endl;
        PrintMatrix(HSc);
        std::cerr << std::endl << "HSd: [" << HSd.rows() << ", " << HSd.cols() << "]" << std::endl;
        PrintMatrix(HSd);

        std::cerr << std::endl << "dSdc: [" << dSdc.rows() << ", " << dSdc.cols() << "]";
        std::cerr << std::endl << "dcdd: [" << dcdd.rows() << ", " << dcdd.cols() << "]";
        std::cerr << std::endl << "HSc: [" << HSc.rows() << ", " << HSc.cols() << "]";
        std::cerr << std::endl << "HSd: [" << HSd.rows() << ", " << HSd.cols() << "]" << std::endl;

        RSVS3D_ERROR_LOGIC("Matrix sizes do not match");
    }
#endif

    HSd = HSd + (dcdd.transpose() * HSc * dcdd);
}

void PrintMatrix(const MatrixXd &mat)
{
    int ii, jj, ni, nj;

    ni = mat.rows();
    nj = mat.cols();
    for (ii = 0; ii < ni; ++ii)
    {
        for (jj = 0; jj < nj; ++jj)
        {
            cout << mat(ii, jj) << " ";
        }
        cout << endl;
    }
}
void PrintMatrixFile(const MatrixXd &mat, const char *name)
{
    ofstream myfile;
    myfile.open(name, ios::app);
    PrintMatrixFile(mat, myfile);
    myfile.close();
}

void PrintMatrixFile(const MatrixXd &mat, std::ostream &myfile)
{
    int ii, jj, ni, nj;

    ni = mat.rows();
    nj = mat.cols();

    for (ii = 0; ii < ni; ++ii)
    {
        for (jj = 0; jj < nj; ++jj)
        {
            myfile << mat(ii, jj) << " ";
        }
        myfile << endl;
    }
    myfile << endl;
}

void PrintMatrix(const VectorXd &mat)
{
    int ii, ni;

    ni = mat.size();
    for (ii = 0; ii < ni; ++ii)
    {
        cout << mat[ii] << " ";
        cout << endl;
    }
}
void PrintMatrix(const RowVectorXd &mat)
{
    int ii, ni;

    ni = mat.size();
    for (ii = 0; ii < ni; ++ii)
    {
        cout << mat[ii] << " ";
    }
    cout << endl;
}

double StreamStatistics(const VectorXd &&vec, ostream &out, const string &&sep)
{
    /*
    Uses a rvalue reference to allow operations to be passed
    directly in the statistics
    */
    double norm = vec.norm();
    out << norm << sep;
    if (vec.size() > 0)
    {
        out << std::setw(out.precision() + 7) << vec.mean() << sep;
        out << std::setw(out.precision() + 7) << sqrt((vec.array() - vec.mean()).square().sum() / (vec.size() - 1))
            << sep;
        out << std::setw(out.precision() + 7) << vec.maxCoeff() << sep;
        out << std::setw(out.precision() + 7) << vec.minCoeff() << sep;
    }
    out << endl;
    return (norm);
}

void StreamOutVector(const VectorXd &&vec, std::ostream &out, const string &&sep)
{
    /*
    Uses a rvalue reference to allow operations to be passed
    directly in the statistics
    */
    int i, n;
    n = vec.size();
    for (i = 0; i < n; ++i)
    {
        out << vec[i] << sep;
    }
    out << endl;
}
