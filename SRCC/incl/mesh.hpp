/**
 * Provides all the mesh tools used for the generation of 3D grids and 
 * geometries.
 * 
 * This file provides the mesh class and it's associated sub-component. These 
 * can be used to robustly control changes in geometry.
 *  
 *@file
 */

//===============================================
// Include Guards
#ifndef MESH_H_INCLUDED
#define MESH_H_INCLUDED

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

//------------------------------------------------------------------------------
// forward declared dependencies class foo; //when you only need a pointer not
// the actual object and to avoid circular dependencies


//------------------------------------------------------------------------------
// included dependencies
#ifdef DBG_MEMLEAK
#define _CRTDBG_MAP_ALLOC  
#include <stdlib.h>  
#include <crtdbg.h>
#endif //DBG_MEMLEAK

#include <iostream>
#include <array>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <string>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <functional>
#include <cmath>
#include <cfloat>

#include "warning.hpp"
#include "arraystructures.hpp"

//------------------------------------------------------------------------------
// Code NOTE: function in a class definition are IMPLICITELY INLINED ie replaced
// by their code at compile time


class meshpart;
class mesh;
class meshdependence; 
class ConnecRemv;

class volu;
class surf;
class vert;
class edge;
class coordvec;

//typedef ArrayStruct<surf> surfarray;
typedef ModiftrackArray<surf> surfarray;
typedef ArrayStruct<volu> voluarray;
typedef ModiftrackArray<edge> edgearray;
typedef ArrayStruct<vert> vertarray;
namespace grid {
	
	typedef std::array<std::array<double, 2>,3> limits;
	/// Defines a linear transformation to the mesh where for each dimension:
	/// {new minimum, old minimum , scaling}
	typedef std::array<std::array<double, 3>,3> transformation;
	/// Defines a list of coordinates
	typedef std::vector<const std::vector<double> *> coordlist;
}

namespace rsvs3d {
	namespace constants {
		static const auto __issetlength =[](double l)->bool {return l>=-0.0;};
		static const double __unsetlength=-1.0;
		namespace ordering {
			static const int ordered=0;
			static const int truncated=1;
			static const int open=-1;
			static const int error=2;
			static const auto __isordered =[](int l)->bool {return l==ordered;};
		}
	}
}


/**
 * Handles the use and norm of a vector for which the norm and the unit value
 *  might be needed.
 *  
 *  Implements some simple mathematical operations for coordinate (3-D) vectors.
*/
class coordvec {
protected:
	std::vector<double> elems;
	double norm;
	int isuptodate;
	

public:
	double CalcNorm();
	double GetNorm();
	double GetNorm() const ;
	void PrepareForUse();
	coordvec Unit() const ;
	double Unit(const int a) const;
	double Normalize();
	void assign(double a, double b, double c);
	double& operator[](int a);
	double operator()(int a) const;
	void disp() const;
	bool isready() const {return(bool(isuptodate));};
	const std::vector<double>& usedata() const {return(elems);}
	const std::vector<double>* retPtr() const {return(&elems);}
	// Math and logical operations (element wise)
	void flipsign();
	void max(const std::vector<double> &vecin);
	void min(const std::vector<double> &vecin);
	void add(const std::vector<double> &vecin);
	void substract(const std::vector<double> &vecin);
	void substractfrom(const std::vector<double> &vecin);
	void div(const std::vector<double> &vecin);
	void div(double scalin);
	void mult(const std::vector<double> &vecin);
	void mult(double scalin);
	void swap(std::vector<double> &vecin);
	void swap(coordvec &coordin);
	std::vector<double> cross(const std::vector<double> &vecin) const ;
	double dot(const std::vector<double> &vecin) const ;
	double angle(const coordvec &coordin) const ;

	coordvec(){
		elems.reserve(3); // reserves 3 as this is the size of the array
		elems.assign(3,0);
		norm=0;
		isuptodate=0;
		#ifdef TEST_SNAKSTRUCT
		std::cout << "constructor called for coordvec" << std::endl;
		#endif
	}
	void operator=(const std::vector<double> &a){
		if(int(a.size())!=3){
			RSVS3D_ERROR_NOTHROW("Warning : Coordinate std::vector is being a "
				"std::vector other than 3 long");
		}
		elems=a;
		isuptodate=0;
	}
};


/**
 * @brief      /Abstract class to ensure mesh interfaces are correct.
 */
class meshpart : public ArrayStructpart { 
public : 
	virtual void disptree(const mesh &meshin, int n) const =0 ;

};


/**
 * @brief      Class for volume cell objects in a mesh.
 */
class volu: public meshpart {
public:
	
	double fill,target,error, volume;
	std::vector<int> surfind;

	std::vector<int> vertind(const mesh &meshin) const;
	void ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu);
	void disp() const;
	void disptree(const mesh &meshin, int n) const;
	void PrepareForUse(){};
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wunused-parameter"
	bool isready(bool isInMesh) const {return(true);}
	#pragma GCC diagnostic pop
	void read(FILE * fid);
	void write(FILE * fid) const;
	void TightenConnectivity() {sort(surfind);unique(surfind);};
	coordvec PseudoCentroid(const mesh &meshin) const;

	volu(){ // Constructor
		this->index=0;
		this->fill=0;
		this->target=1;
		this->error=1;

		#ifdef TEST_ARRAYSTRUCT
		std::cout << "volu #" << index << " Was created " << surfind.size() << std::endl;
		#endif
	}
	volu(const volu& oldVolu){ // Copy-Constructor
		this->index=oldVolu.index;
		this->fill=oldVolu.fill;
		this->target=oldVolu.target;
		this->error=oldVolu.error;
		this->surfind=oldVolu.surfind;

		#ifdef TEST_ARRAYSTRUCT
		std::cout << "copyvolu #" << index << " Was created " 
			<< surfind.size() << std::endl;
		#endif
	}
	~volu(){ // Destructor
		surfind.clear();

		#ifdef TEST_ARRAYSTRUCT
		std::cout << "volu #" << index << " Was deleted " << surfind.size() << std::endl;
		#endif

	}
	void operator=(const volu* other){
		index=other->index;
		fill=other->fill;
		target=other->target;
		error=other->error;
		surfind=other->surfind;

		#ifdef TEST_ARRAYSTRUCT
		std::cout << "OTHER: " ; other->disp();
		#endif
	}

	int Key() const {return(index);}
};



/**
 * @brief      Class for surface object in a mesh.
 */
class surf: public meshpart , public modiftrackpart {
protected:
	bool isordered;
	//bool isModif=true;
public:
	friend class mesh;
	friend surfarray;
	// friend void mesh::SwitchIndex(int typeInd, int oldInd, int newInd,
	// vector<int> scopeInd); friend void mesh::RemoveIndex(int typeInd, int
	// oldInd);

	double fill,target,error,area;
	std::vector<int> edgeind;
	std::vector<int> voluind;
	 // reserves 2 as this is the size of the array
	
	std::vector<int> vertind(const mesh &meshin) const;
	void disp() const;
	void disptree(const mesh &meshin, int n) const;
	void ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu);
	void PrepareForUse(){};	
	bool isready(bool isInMesh) const {return(isInMesh? isordered : true);}
	bool IsOrdered() const {return(isordered);}
	void read(FILE * fid);
	void write(FILE * fid) const;
	int OrderEdges(mesh *meshin);
	int SplitSurface(mesh &meshin, const std::vector<int> &fullEdgeInd);
	void OrderedVerts(const mesh *meshin, std::vector<int> &vertList) const;
	std::vector<int> OrderedVerts(const mesh *meshin) const ;
	void TightenConnectivity() {sort(voluind);unique(voluind);
		sort(edgeind);unique(edgeind);isordered=false;};
	void FlipVolus();
	bool edgeconneq(const surf &other, bool recurse=true) const;
	coordvec PseudoCentroid(const mesh &meshin) const;
	int PseudoCentroid(const mesh &meshin, coordvec &coord) const;


	surf(){ // Constructor
		index=0;
		fill=1;
		target=1;
		error=1;
		voluind.reserve(2); // reserves 2 as this is the size of the array
		voluind.assign(2,0);
		isordered=false;
	}
	~surf(){ // Destructor

		edgeind.clear();
		voluind.clear();

	}
	surf(const surf& oldSurf){ // Copy-Constructor
		index=oldSurf.index;
		fill=oldSurf.fill;
		target=oldSurf.target;
		error=oldSurf.error;
		edgeind=oldSurf.edgeind;
		voluind=oldSurf.voluind;
		isordered=oldSurf.isordered;
	}
	void operator=(const surf* other){
		index=other->index;
		fill=other->fill;
		error=other->error;
		target=other->target;
		edgeind=other->edgeind;
		voluind=other->voluind;
		isordered=other->isordered;
	}

	int Key() const {return(index);}

};


/**
 * @brief      Class for an edge object in a mesh.
 */
class edge: public meshpart , public modiftrackpart {
protected:
	double length=rsvs3d::constants::__unsetlength;
public:
	friend class mesh;
	friend edgearray;
	

	std::vector<int> vertind;
	std::vector<int> surfind;
	 // reserves 2 as this is the size of the array

	void ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu);
	void disp() const;
	void disptree(const mesh &meshin, int n) const;
	void PrepareForUse(){};
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wunused-parameter"
	bool isready(bool isInMesh) const {return(true);}
	#pragma GCC diagnostic pop
	void read(FILE * fid);
	void write(FILE * fid) const;
	void TightenConnectivity() {
		if(vertind.size()>2){
			RSVS3D_ERROR_ARGUMENT("vertind should be size 2");
		}
		sort(surfind);unique(surfind);
	};
	void GeometricProperties(const mesh *meshin, coordvec &centre,
		double &length) const;
	double Length(const mesh &meshin) const;
	double SetLength(const mesh &meshin){
		this->length = this->Length(meshin);
		return this->length;};
	double GetLength(bool warn=true) const {
		if(!rsvs3d::constants::__issetlength(this->length) && warn){
			RSVS3D_ERROR_NOTHROW("Length is accessed but not set. "
				"Run SetEdgeLengths before call.");
		}
		return this->length;
	}
	void InvalidateLength() {this->length=rsvs3d::constants::__unsetlength;}
	double LengthSquared(const mesh &meshin) const;
	bool IsLength0(const mesh &meshin, double eps=__DBL_EPSILON__) const;
	bool vertconneq(const edge &other) const;
	edge(){ // Constructor
		index=0;
		vertind.reserve(2);
		vertind.assign(2,0); // reserves 2 as this is the size of the array

	}
	edge(const edge& oldEdge){ // Copy-Constructor
		index=oldEdge.index;
		vertind=oldEdge.vertind;
		surfind=oldEdge.surfind;
	}
	~edge(){ // Destructor

		vertind.clear();
		surfind.clear();

	}
	void operator=(const edge* other){
		index=other->index;

		vertind=other->vertind;
		surfind=other->surfind;
	}

	int Key() const {return(index);}
};

/**
 * @brief      Class for a vertex in a mesh.
 */
class vert: public meshpart {
public:
	

	std::vector<int> edgeind;
	std::vector<double> coord;
	 // reserves 2 as this is the size of the array
	std::vector<int> elmind(const mesh &meshin, int dimOveride=-1) const;

	void disp() const;
	void disptree(const mesh &meshin, int n) const;
	void ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu);
	void PrepareForUse(){};
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wunused-parameter"
	bool isready(bool isInMesh) const {return(true);}
	#pragma GCC diagnostic pop
	void read(FILE * fid);
	void write(FILE * fid) const;
	void TightenConnectivity() {sort(edgeind);unique(edgeind);};
	int OrderEdges(const mesh *meshin);
	std::pair<std::vector<int>,int> OrderEdges(const mesh *meshin) const;
	int OrderEdges(const mesh *meshin, std::vector<int> &edgeIndOut) const;
	int SurroundingCoords(const mesh *meshin, grid::coordlist &coordout,
		bool isOrdered=false, std::vector<int>* edgeIndOutPtr=NULL) const;

	int Normal(const mesh *meshin, grid::coordlist &neighCoord,
		coordvec &normalVec, bool isOrdered=false) const;
	coordvec Normal(const mesh *meshin) const;


	vert(){ // Constructor
		index=0;
		coord.reserve(3); // reserves 2 as this is the size of the array
		coord.assign(3,0);
	}
	vert(const vert& oldEdge){ // Copy-Constructor
		index=oldEdge.index;
		edgeind=oldEdge.edgeind;
		coord=oldEdge.coord;
	}
	~vert(){ // Destructor

		edgeind.clear();
		coord.clear();

	}
	void operator=(const vert* other){
		index=other->index;

		edgeind=other->edgeind;
		coord=other->coord;
	}

	int Key() const {return(index);}
};



/**
 * @brief      Class containing the information needed to trim objects from a 
 * mesh.
 */
class ConnecRemv {
public:
	int keepind;
	int typeobj;
	std::vector<int> rmvind;
	std::vector<int> scopeind;
	void disp();
	
};


/**
 * @brief      Class for connecting meshes.

 * Stores a vector of mesh references for parent and children. Needs to support
 * partial meshes for constraint handling. 
*/
class meshdependence {
protected:
	friend class mesh;
	/// Number of parent meshes
	int nParents = 0;
	/// Indices of the active elements of the owning mesh
	std::vector<int> elemind;
	/// Vector of pointers to the mesh which are coarser (parents).
	std::vector<mesh*> parentmesh;
	/// Vector of pointers to the mesh which are finer (children).
	std::vector<mesh*> childmesh;
	/// parent/to self connectivity, 1 vector element per parent. This is an 
	/// vector with the index of each parent element stored at the location of 
	/// each self element.
	std::vector<HashedVectorSafe<int,int>> parentconn;
	// These methods are protected to avoid broken/uni-directional
	// connectivities being generated
	
	int AddParent(mesh* meshin);
	int AddChild(mesh* meshin);
	void AddParent(mesh* meshin, std::vector<int> &parentind);
	void RemoveChild(mesh* meshin);
	void RemoveParent(mesh* meshin);

public:
	int NumberOfParents() const {return this->nParents;}
	const mesh* ParentPointer(int a) const {return this->parentmesh.at(a);}
	std::vector<int> ChildIndices(int parent, int parentVoluIndex) const {
		return this->parentconn.at(parent).findall(parentVoluIndex);}
};


/**
 * @brief      Class for mesh handling.
 * 
 * Class implementing the functionality of this file. The mesh class allow, the
 * robust evolution of a grid. Element connectivity is stored bi-directionnaly.
 * This allows no connectivity to need to be infered and allows very fast and
 * robust traversing of the mesh by using hashed lists of the indices of the
 * mesh components.
 */
class mesh {
private:
	bool borderIsSet=false;
	bool meshDepIsSet=false;
	bool facesAreOriented=false;
	bool edgesLengthsAreSet=false;
	int meshDim=0;
	void SetLastIndex();
	
	void OrientSurfaceVolume();
	void OrientEdgeSurface();
	int OrientRelativeSurfaceVolume(std::vector<int> &surfOrient);
	void ArraysAreHashed();
	void _LinearTransformGeneration(const grid::transformation &transform,
		std::vector<mesh*> meshdependence::*mp);
	friend class snake;
public:
	vertarray verts;
	edgearray edges;
	surfarray surfs;
	voluarray volus;

	meshdependence meshtree;
	// Mesh Lineage
	void RemoveFromFamily();
	void AddChild(mesh* meshin);
	void AddParent(mesh* meshin);
	void AddParent(mesh* meshin, std::vector<int> &parentind);
	void AddChild(mesh* meshin, std::vector<int> &parentind);
	void SetMeshDepElm();
	// Method needed to robustly maintain lineage through the family.
	void MaintainLineage();
	int CountParents() const;
	int SurfInParent(int surfind) const;
	void SurfInParent(std::vector<int> &listInParent) const;
	void ElmOnParentBound(std::vector<int> &listInParent, std::vector<int> &voluInd,
		bool isBorderBound=true,
		bool outerVolume=true) const;
	void SurfOnParentBound(std::vector<int> &listInParent, std::vector<int> &voluInd,
		bool isBorderBound,
		bool outerVolume) const;
	void EdgeOnParentBound(std::vector<int> &listInParent, std::vector<int> &voluInd,
		bool isBorderBound,
		bool outerVolume) const;
	int CountVoluParent() const ;
	void ReturnParentMap(std::vector<int> &currind, std::vector<int> &parentpos,
		std::vector<std::pair<int,int>> &parentcases, std::vector<double> &voluVals) const;
	void MapVolu2Parent(const std::vector<double> &fillIn,
		const std::vector<std::pair<int,int>> &parentcases, double volu::*mp=&volu::fill);
	void MapVolu2Self(const std::vector<double> &fillIn, 
		const std::vector<int> &elms, double volu::*mp=&volu::fill);
	void VoluValuesofParents(int elmInd, std::vector<double> &vals,
		int volType=0) const;
	void VoluValuesofParents(int elmInd, std::vector<double> &vals,
		double volu::*mp) const;
	void SurfValuesofParents(int elmInd, std::vector<double> &vals,
		int volType=0) const;
	void SurfValuesofParents(int elmInd, std::vector<double> &vals,
		double surf::*mp) const;
	int ParentElementIndex(int childElmInd, int parentInd=0) const;
	// Mesh property
	int WhatDim() const {return(meshDim);}
	// basic operations grouped from each field
	void HashArray();
	void SetMaxIndex();
	void GetMaxIndex(int *nVert,int *nEdge,int *nSurf,int *nVolu) const;
	void Init(int nVe,int nE, int nS, int nVo);
	void size(int &nVe,int &nE, int &nS, int &nVo) const;
	void reserve(int nVe,int nE, int nS, int nVo);
	void PrepareForUse(bool needOrder=true);
	void SetEdgeLengths();
	void InvalidateEdgeLength(int iEdge);
	void disp() const;
	void displight() const;
	void Concatenate(const mesh &other);
	bool isready() const;
	void PopulateIndices();
	void TightenConnectivity();
	int TestConnectivity(const char *strRoot="") const;
	int TestConnectivityBiDir(const char *strRoot="",
		bool emptyIsErr=true) const;
	// File I/o
	void write(FILE *fid) const;
	void read(FILE *fid);
	int write(const char *str) const;
	int read(const char *str);
	// Mesh merging
	void MakeCompatible_inplace(mesh &other) const;
	mesh MakeCompatible(mesh other) const;
	void ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu);
	void SwitchIndex(int typeInd, int oldInd, int newInd,
		const std::vector<int> &scopeInd={0});
	void RemoveIndex(int typeInd, int oldInd);
	int ConnectedVertex(std::vector<int> &vertBlock) const;
	int ConnectedVolumes(std::vector<int> &volBlock, 
		const std::vector<bool> &boundaryFaces={}) const;
	void ForceCloseContainers();
	void RemoveSingularConnectors(const std::vector<int> &rmvVertInds={},
		bool voidError=true);
	std::vector<int> MergeGroupedVertices(HashedVector<int, int> &closeVert,
		bool delVerts=true);
	// Mesh Quality
	std::vector<int> OrderEdges();
	void SetBorders();
	void OrientFaces();
	int OrderVertexEdges(int vertIndex);
	// Mesh component comparison
	void GetOffBorderVert(std::vector<int> &vertList, std::vector<int> &voluInd,
		int outerVolume=-1);
	void GetOffBorderVert(std::vector<int> &vertList, std::vector<int> &voluInd,
		int outerVolume=-1) const;
	void GetOffBorderVert3D(std::vector<int> &vertList, std::vector<int> &voluInd,
		int outerVolume=-1) const;
	void GetOffBorderVert2D(std::vector<int> &vertInd, std::vector<int> &surfind,
		int outerVolume=-1) const;
	// Mesh calculations
	coordvec CalcCentreVolu(int ind) const;
	coordvec CalcPseudoNormalSurf(int ind) const;
 	std::vector<int> VertexInVolume(const std::vector<double> testVertices,
 		int sizeVert=3) const; 
	// Mesh size and position
	grid::transformation Scale();
	grid::transformation Scale(const grid::limits &domain);
	void LinearTransform(const grid::transformation &transform);
	void LinearTransformFamily(const grid::transformation &transform);
	void LoadTargetFill(const std::string &fileName);
	grid::limits BoundingBox() const;
	void ReturnBoundingBox(std::array<double,3> &lowerB, 
		std::array<double,3> &upperB) const ;
	// Mesh Splitting and cropping
	void Crop(std::vector<int> indList, int indType=1);
	std::vector<int> AddBoundary(const std::vector<double> &lb, const std::vector<double> &ub);
	void CropAtBoundary(const std::vector<double> &lb, const std::vector<double> &ub);
	// Mesh traversal convenience functions
	int EdgeFromVerts(int v1, int v2) const;
	int SurfFromEdges(int e1, int e2, int repetitionBehaviour=-1) const;
	int VertFromVertEdge(int v, int e) const;
	void VerticesVector(int v1, int v2, coordvec &vec) const;
	void EdgeVector(int e, coordvec &vec) const;
	// Destructor
	~mesh(){
		RemoveFromFamily();
	}
};

// Function declarations

void ConnVertFromConnEdge(const mesh &meshin, const std::vector<int> &edgeind,
	std::vector<int> &vertind);
std::pair<int, int> OrderMatchLists(const std::vector<int> &vec1, int p1, int p2);
int OrderMatchLists(const std::vector<int> &vec1, const std::vector<int> &vec2,
	int p1, int p2);
void CropMeshGreedy(mesh &meshin, const std::vector<double> &lb,
		const std::vector<double> &ub);
int OrderEdgeList(std::vector<int> &edgeind, const mesh &meshin, bool warn=true,
	bool errout=true, const std::vector<int>* edgeIndOrigPtr=NULL,
	const surf* surfin=NULL);
int OrderList(std::vector<int> &edgeind, const std::vector<int> &edge2Vert, bool warn=true,
	bool errout=true, const std::vector<int>* edgeIndOrigPtr=NULL);
void DiffPointsFromCentre(const std::vector<double> &centre, 
	const std::vector<double> &planeVert2,
	const std::vector<double> &planeVert3,
	coordvec &normal, coordvec &temp1);
void DiffPoints(const std::vector<double> &vert1,
	const std::vector<double> &vert2, coordvec &diffVerts);
double Angle3Points(const std::vector<double> &centre, 
	const std::vector<double> &planeVert2,
	const std::vector<double> &planeVert3,
	coordvec &vec1, coordvec &vec2);
double VertexDistanceToPlane(const std::vector<double> &planeVert1, 
	const std::vector<double> &planeVert2,
	const std::vector<double> &planeVert3,
	const std::vector<double> &testVertex,
	coordvec &temp1, 
	coordvec &temp2);
std::vector<double> VerticesDistanceToPlane(const std::vector<double> &planeVert1, 
	const std::vector<double> &planeVert2,
	const std::vector<double> &planeVert3,
	const std::vector<double> &testVertices,
	coordvec &temp1, 
	coordvec &temp2);
double VertexDistanceToPlane(const std::vector<double> &planeVert1, 
	const std::vector<double> &planeVert2,
	const std::vector<double> &planeVert3,
	const std::vector<double> &testVertex);
std::vector<double> VerticesDistanceToPlane(const std::vector<double> &planeVert1, 
	const std::vector<double> &planeVert2,
	const std::vector<double> &planeVert3,
	const std::vector<double> &testVertices);
mesh Points2Mesh(const std::vector<double> &vecPts, int nProp=3);
double PlanesDotProduct(const std::vector<double> &planeVert1, 
	const std::vector<double> &planeVert2,
	const std::vector<double> &planeVert3,
	const std::vector<double> &planeVert4, 
	const std::vector<double> &planeVert5,
	const std::vector<double> &planeVert6, 
	bool normalize=true);
void PlaneNormal(const std::vector<double> &planeVert1, 
	const std::vector<double> &planeVert2,
	const std::vector<double> &planeVert3,
	coordvec &normal, coordvec &temp1);
double PlaneNormalAndAngle(const std::vector<double> &planeVert1, 
	const std::vector<double> &planeVert2,
	const std::vector<double> &planeVert3,
	coordvec &normal, coordvec &temp1);
std::tuple<coordvec,double> VertexNormal(const std::vector<double>& centre, 
	const grid::coordlist &vecPts);

namespace meshhelp {
	template<class T, class V, class W>
	double ProjectRay(int count, const W &&boundBox,
		const T &dir, const V &orig, double minDist=0.0);

	void PlaceBorderVertex(const std::vector<double> &coordIn, 
		const std::vector<double> &coordOut,
		const std::vector<double> &lb,
		const std::vector<double> &ub, std::vector<double> &coordTarg);

	void SplitBorderSurfaceEdgeind(const mesh &meshin,
		const std::vector<bool> &edgeOut, 
		std::vector<int> &vecconnIn, std::vector<int> &vecconnOut);
	void SplitBorderVolumeSurfind(const mesh &meshin, 
		const std::vector<bool> &edgeOut, 
		std::vector<int> &vecconnIn, std::vector<int> &vecconnOut);
	void HandleMultiSurfaceSplit(mesh &meshin, 
		std::vector<int> &edgeindOld, std::vector<int> &edgeindNew,
		std::vector<int> &vertindNew);
	std::vector<int> FindVertInFromEdgeOut(const mesh &meshin, 
		const std::vector<bool> &vertOut,
		const std::vector<int> &edgeList, 
		const std::vector<int> &edgeListCheck);
	std::vector<int> FindEdgeInFromSurfOut(const mesh &meshin, 
		const std::vector<bool> &edgeOut,
		std::vector<int> surfList);
	double VerticesDistanceSquared(const mesh &meshin,
		const std::vector<int> &vertind);
	double VerticesDistance(const mesh &meshin,
		const std::vector<int> &vertind);
	bool IsVerticesDistance0(const mesh &meshin,
		const std::vector<int> &vertind, double eps=__DBL_EPSILON__);
	int VertexInVolume(const mesh &meshin, 
		const std::vector<double> testCoord,
		bool needFlip=false);
	int Get3PointsInSurface(const mesh &meshin, int surfCurr, 
		std::array<int, 3> &surfacePoints);
	int NormalShouldFlip(const std::vector<int> orderedList, int elm1, int elm2,
		const std::vector<int> & voluind, bool innerComparison);
}
//test functions
int Test_ArrayStructures();
int Test_Volu();
int Test_Surf();
int Test_Vert();
int Test_Edge();
int Test_Mesh();
int Test_Crop();

  
template<class T, class V, class W>
double meshhelp::ProjectRay(int count, const W &&boundBox,
	const T &dir, const V &orig, double minDist){
	/*
	Calculates the distance to project a single ray to reach the bounding box

	It is templated to accept any types of containers with `count` elements, 
	accessible by operator[];

	W is an object of size [2][count] accessible by operators [][]
	*/

	double l=-INFINITY;

	for (int i = 0; i < count; ++i)
	{
		l = std::max(l,
				std::min(
					(boundBox[0][i]-orig[i])/dir[i] ,
				(boundBox[1][i]-orig[i])/dir[i] ));
	}
	return std::min(l, minDist);
}

#endif // MESH_H_INCLUDED
