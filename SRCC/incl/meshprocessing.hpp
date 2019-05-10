/**
 * Tools for the mathematical processing of meshes. 
 *  
 *@file
 */

//===============================================
// Include Guards
#ifndef MESHPROCESSING_H_INCLUDED
#define MESHPROCESSING_H_INCLUDED

//===============================================
// Levels of debuging Guards
#ifdef DEBUGLVL2 // All Debugging calls
#define DEBUGLVL1
#define TEST_RANGE // deprecated
#define TEST_ARRAYSTRUCT

#endif

#ifdef DEBUGLVL1 // Debugging of new features.
#ifndef SAFE_ACCESS
#define SAFE_ACCESS
#endif
#endif


//=================================
// forward declared dependencies
// 		class foo; //when you only need a pointer not the actual object
// 		and to avoid circular dependencies


//=================================
// included dependencies

#include <vector>

#include "arraystructures.hpp"
#include "mesh.hpp"
#include "snake.hpp"
// 

//==================================
// Code
// NOTE: function in a class definition are IMPLICITELY INLINED 
//       ie replaced by their code at compile time
//       

void FlattenBoundaryFaces(mesh &meshin);
void TriangulateAllFaces(mesh &meshin);
std::vector<int> FindHolesInSnake(const snake &snakein,
	const HashedVector<int, int> &uncertainVert);
void PrepareSnakeForCFD(const snake &snakein, double distanceTol,
	mesh &meshgeom, std::vector<double> &holeCoords);
HashedVector<int, int> GroupCloseSnaxels(const snake &snakein,
 	double distTol);
void TestVertClose(int vertIndIn, std::vector<bool> &isSnaxDone, 
	const mesh &meshin, double distTol,
	std::vector<int> &sameEdges);
HashedVector<int, int> GroupCloseVertices(const mesh &meshin,
	double distTol);
int FindVertexHole(int vertInd, const mesh &meshin, 
	const std::vector<bool> &vertIn,
	const HashedVector<int, int> &uncertainVert,
	std::vector<bool> &vertExplored);
double DomInter(double x, double y1, double y2);
mesh BuildDomain(const std::array<double,3> &lowerB, 
	const std::array<double,3> &upperB, double tolInner=0.0);
mesh BuildCutCellDomain(const std::array<double,3> &outerLowerB, 
	const std::array<double,3> &outerUpperB,
	const std::array<double,3> &innerLowerB, 
	const std::array<double,3> &innerUpperB, int nSteps, 
	std::vector<int> &vertPerSubDomain);
double PseudoSurfaceAngle(const mesh &meshin, 
	const std::array<int, 2> &surfInds);
std::vector<double> CalculateEdgeCurvature(const mesh &meshin);
std::vector<double> CalculateVertexCurvature(const mesh &meshin,
	int smoothingSteps);
std::vector<double> CalculateVertexMinEdgeLength(const mesh &meshin);
std::vector<double> CalculateVertexMeanEdgeLength(const mesh &meshin);
std::vector<double> CalculateEdgeLengths(const mesh &meshin);
std::vector<double> CoordInVolume(const mesh &meshin);
std::vector<double> VolumeCentroids(const mesh &meshin);
std::vector<double> VolumeInternalLayers(const mesh &meshin, int nLayers);
std::vector<double> SurfaceCentroids(const mesh &meshin);
std::vector<double> SurfaceInternalLayers(const mesh &meshin, int nLayers);
double CotanLaplacianWeight(const std::vector<double> &centre,
	const std::vector<double> &p_im1, const std::vector<double> &p_i, 
	const std::vector<double> &p_ip1, coordvec &temp1, coordvec &temp2);
int VertexLaplacianVector(const mesh& meshin, const vert* vertsmooth,
	coordvec &lapVec, bool isOrdered=false);
coordvec VertexLaplacianVector(const mesh& meshin, int vertIndex);
std::array<double, 2> IntersectLineSphere(const coordvec &lineVec, 
	const coordvec &offset, double sphereRadius);
std::array<double, 2> IntersectLineSphere(const coordvec &lineVec,
	const std::vector<double> &linePoint, const coordvec &sphereCentre,
	double sphereRadius);
std::vector<double> MeshUnitNormals(const mesh& meshin);
std::vector<double> MeshLaplacians(const mesh& meshin);
// Forward declared templated functions

// Base class



#endif // MESH_H_INCLUDED