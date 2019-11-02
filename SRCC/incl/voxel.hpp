/**
 * Generation of cartesian grids.
 * 
 *@file
 */

//===============================================
// Include Guards
#ifndef VOXEL_H_INCLUDED
#define VOXEL_H_INCLUDED

//===============================================
// Levels of debuging Guards
#ifdef DEBUGLVL2 // All Debugging calls
#define DEBUGLVL1
#define TEST_VOXEL_SURF
#define TEST_VOXEL_VOLU
#define TEST_EIGENEXT
#define TEST_VOXEL_EDGE
#define TEST_VOXEL_VERT
#define TEST_VOXEL
#endif

#ifdef DEBUGLVL1 // Debugging of new features.
//#define TEST_EIGENEXT
#endif

//=================================
// forward declared dependencies
// 		class foo; //when you only need a pointer not the actual object
// 		and to avoid circular dependencies


//=================================
// included dependencies

#include <iostream>
#include <numeric>      // std::partial_sum
#include <Eigen>

#include "arraystructures.hpp"
#include "postprocessing.hpp"




//==================================
// Code
// NOTE: function in a class definition are IMPLICITELY INLINED 
//       ie replaced by their code at compile time
// Templates

/**
 * @brief      template which applies cumulative sum to Eigen Matrix.
 *
 *	Cumprod is applied row-wise for d=0 and col-wise for d=1
 *
 * @param[in]  matIn  The matrix input
 * @param[in]  d      dimension to use 0-row wise, 1 col-wise
 *
 * @tparam     T      Eigen type
 *
 * @return     The cumulative sum.
 */
template <class T> T cumsum(const T &matIn, int d) {
	T mat=matIn;
	if (!d) {
		for (int i=0;i<mat.rows();++i){
			for (int j=1;j<mat.cols();++j){
				mat(i,j)+=mat(i,j-1);
			}
		}
	}else{
		for (int i=1;i<mat.rows();++i){
			for (int j=0;j<mat.cols();++j){
				mat(i,j)+=mat(i-1,j);
			}
		}
	}
	return(mat);

}

/**
 * @brief      template which applies cumulative product to Eigen Matrix.
 *
 *	Cumprod is applied row-wise for d=0 and col-wise for d=1
 *
 * @param[in]  matIn  The matrix input
 * @param[in]  d      dimension to use 0-row wise, 1 col-wise
 *
 * @tparam     T      Eigen type
 *
 * @return     The cumulative product.
 */
template <class T> T cumprod(const T &matIn, int d) {
	T mat=matIn;
	if (!d) {
		for (int i=0;i<mat.rows();++i){
			for (int j=1;j<mat.cols();++j){
				mat(i,j)*=mat(i,j-1);
			}
		}
	}else{
		for (int i=1;i<mat.rows();++i){
			for (int j=0;j<mat.cols();++j){
				mat(i,j)*=mat(i-1,j);
			}
		}
	}
	return(mat);

}

// functions
int BuildBlockGrid(std::array<int, 3> &dimGrid, mesh& blockGrid);
int BuildBlockGrid(const Eigen::RowVector3i  &dimGrid, mesh& blockGrid);
int BuildBlockVert(const Eigen::RowVector3i  &dimGrid, mesh& blockGrid, int nVert, 
	const Eigen::Matrix3i  &edgeProp, const Eigen::RowVector3i  &nEdgeDim);

int BuildBlockEdge(const Eigen::RowVector3i  &dimGrid, mesh& blockGrid, int nEdge ,
	const Eigen::RowVector3i  &nEdgeDim,  const Eigen::RowVector3i  &nSurfDim, const Eigen::Matrix3i  &edgeProp,
	const Eigen::Matrix3i  &surfProp );

int BuildBlockSurf(const Eigen::RowVector3i  &dimGrid, int nSurf ,mesh& blockGrid ,
	const Eigen::Matrix3i  &surfProp, const Eigen::Matrix3i  &edgeProp, const Eigen::RowVector3i  &nSurfDim,
	 const Eigen::RowVector3i  &nEdgeDim);

int BuildBlockVolu(const Eigen::RowVector3i  &dimGrid, int nVolu , mesh& blockGrid,
	const Eigen::RowVector3i  &nSurfDim, const Eigen::Matrix3i  &surfProp);
//test functions
int Test_BuildBlockGrid_noout();
int Test_MeshOut();

#endif // VOXEL_H_INCLUDED

