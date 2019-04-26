/**
 * Provides the core restricted surface snake container.
 * 
 * This container allows efficient and robust geometry and topology evolution.
 * 
 *@file
 */

//===============================================
// Include Guards
#ifndef SNAKSTRUCT_H_INCLUDED
#define SNAKSTRUCT_H_INCLUDED

//===============================================
// Levels of debuging Guards
#ifdef DEBUGLVL2 // All Debugging calls
#define DEBUGLVL1

#define TEST_SNAKSTRUCT

#endif

#ifdef DEBUGLVL1 // Debugging of new features.
#endif

//=================================
// forward declared dependencies
// 		class foo; //when you only need a pointer not the actual object
// 		and to avoid circular dependencies


//=================================
// included dependencies
#include <iostream>
#include <vector>
#include <cmath>
#include <cfloat>
#include <unordered_map>

#include "arraystructures.hpp"
#include "mesh.hpp"

//==================================
// Code
// NOTE: function in a class definition are IMPLICITELY INLINED 
//       ie replaced by their code at compile time
using namespace std;

class snax;
class snaxedge;
class snaxsurf;
class snaxarray;

typedef SnakStruct<snaxedge> snaxedgearray;
typedef SnakStruct<snaxsurf> snaxsurfarray;

class snaxarray : public SnakStruct<snax> 
{
protected:
	using SnakStruct<snax>::readyforuse;
	int isOrderedOnEdge=0;

public: 
	friend class snake;
	friend class snax;
	friend void SpawnArrivedSnaxelsDir(snake &fullsnake, snake &partSnake,
		const vector<int> &isImpact, int dir);

	void ReorderOnEdge();
	void OrderOnEdge();
	void CalculateTimeStepOnEdge(vector<double> &dt, vector<bool> &isSnaxDone,
		int edgeInd);
	void DetectImpactOnEdge(vector<int> &isImpact, vector<bool> &isSnaxDone,
		int edgeInd);
	// Functions that need modification
	bool checkready();
	void ForceArrayReady();
	void PrepareForUse();
	void Concatenate(const snaxarray &other);
	snax& operator[](const int a){ 
		isOrderedOnEdge=0;
		return(this->SnakStruct<snax>::operator[](a));
	}

};



class snake  {
private:

	bool snaxDistanceLimit_conserveShape = true;

	bool is3D=true;
	void SetLastIndex(); // Baaaad function do not use if you're not sure.
	bool isFlipped=false;
	void OrientSurfaceVolume();
	void OrientEdgeSurface();
public:
	// Handling of data specific to snake
	snaxarray snaxs; 			// properties associated with snakconn verts
	snaxedgearray snaxedges; 	// properties associated with snakconn edges
	snaxsurfarray snaxsurfs; 	// properties associated with snakconn surfs
	// Using the mesh container to store connectivity
	mesh snakeconn;
	// pointer to snaking mesh
	mesh *snakemesh=NULL;
	vector<bool> isMeshVertIn;
	// basic operations grouped from each field
	void disp() const;

	void displight() const;
	bool isready() const ;
	void PrepareForUse(bool needOrder=true);
	void Init(mesh *snakemesh,int nSnax, int nEdge, int nSurf, int nVolu);
	void reserve(int nSnax, int nEdge, int nSurf, int nVolu);
	inline void GetMaxIndex(int *nVert,int *nEdge,int *nSurf,int *nVolu) const;
	void HashArray(); // Not really needed as handled by PrepareForUse
	void HashArrayNM(); // Not really needed as handled by PrepareForUse
	void HashParent();
	void SetMaxIndex(); // Not really needed as handled by PrepareForUse
	void SetMaxIndexNM(); // SetMaxIndex no mesh
	void Concatenate(const snake &other, int isInternal=0);
	bool Check3D() const {return(is3D);}
	// Snake merging
	void MakeCompatible_inplace(snake &other) const;
	snake MakeCompatible(snake other) const;
	void ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu);
	void ChangeIndicesSnakeMesh(int nVert,int nEdge,int nSurf,int nVolu);
	void ForceCloseContainers();
	// Snake Movement
	void UpdateDistance(double dt,  double maxDstep=1.0);
	void UpdateDistance(const vector<double> &dt,  double maxDstep=1.0);
	void CalculateTimeStep(vector<double> &dt, double dtDefault, 
		double distDefault=1.0);
	void SnaxImpactDetection(vector<int> &isImpact);
	void SnaxAlmostImpactDetection(vector<int> &isImpact, double dDlim);
	void UpdateCoord();
	void Flip(); // reverses snake directions
	grid::limits Scale(const grid::limits &newSize);
	// Snake connectivity operations
	void OrderEdges();
	void SetSnaxSurfs() {}
	void OrientFaces();
	int FindBlockSnakeMeshVerts(vector<int> &vertBlock) const;
	void AssignInternalVerts();
	void CheckConnectivity() const; 
	void TakeSpawnStep(int minIndex, double stepLength);

	void VertIsIn(int vertInd, bool isIn=true);
	void VertIsIn(vector<int> vertInd, bool isIn=true);
	bool ReturnFlip()const{return(isFlipped);}
	// io of snake
	void read(FILE *fid);
	void write(FILE *fid) const;
	int read(const char *str);
	int write(const char *str) const;
	void SetSnaxDistanceLimit_conserveShape(bool in){
		this->snaxDistanceLimit_conserveShape = in;
	}
};

 
class snax : public meshpart , public snakpart {	
public:
	
	double d=0.0;
	double v=0.0;
	int fromvert=0; 	// root vertex of *snakemesh
	int tovert=0; 	// destination vertex of *snakemesh
	int edgeind=0; 	// edge of snakemesh
	int isfreeze=0; 	// freeze status 
	int orderedge=0;	// order on the edge

	// interface functions
	void disp() const;
	void disptree(const mesh &meshin, int n) const;
	void disptree(const snake &snakein, int n) const;
	int Key() const {return (index);};
	int KeyParent() const {return (edgeind);};
	void ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu);
	void ChangeIndicesSnakeMesh(int nVert,int nEdge,int nSurf,int nVolu);
	void PrepareForUse(){};
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wunused-parameter"
	bool isready(bool isInMesh) const {return(true);}
	#pragma GCC diagnostic pop
	void read(FILE * fid);
	void write(FILE * fid) const;
	inline void set(int index, double d,double v,int fromvert,
		int tovert,int edgeind,int isfreeze,int orderedge);
	void SwitchIndex(int typeInd, int oldInd, int newInd);
	void TightenConnectivity(){}
	void TakeSpawnStep(snake &snakein, double stepLength);
};


class snaxedge : public meshpart , public snakpart {
public:
	int surfind=0;
	coordvec normvector;

	void PrepareForUse();
	void disp() const;
	void disptree(const mesh &meshin, int n) const;
	void disptree(const snake &snakein, int n) const;
	int Key() const {return (index);};
	int KeyParent() const {return (surfind);}; 
	void ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu);
	void ChangeIndicesSnakeMesh(int nVert,int nEdge,int nSurf,int nVolu);
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wunused-parameter"
	bool isready(bool isInMesh) const {return(normvector.isready());}
	#pragma GCC diagnostic pop
	void read(FILE * fid);
	void write(FILE *fid)const;
	void SwitchIndex(int typeInd, int oldInd, int newInd);
	void TightenConnectivity(){}

};

class snaxsurf : public meshpart , public snakpart {
public: 
	
	int voluind=0;
	coordvec normvector;
	void PrepareForUse(); 
	void disp() const;
	void disptree(const mesh &meshin, int n) const;
	void disptree(const snake &snakein, int n) const;
	int Key() const {return (index);};
	int KeyParent() const {return (voluind);};
	void ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu);
	void ChangeIndicesSnakeMesh(int nVert,int nEdge,int nSurf,int nVolu);
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wunused-parameter"
	bool isready(bool isInMesh) const {return(normvector.isready());}
	#pragma GCC diagnostic pop
	void read(FILE * fid);
	void write(FILE * fid)const;
	void SwitchIndex(int typeInd, int oldInd, int newInd);
	void TightenConnectivity(){}


};

// Function prototypes
double SnaxImpactDt(const snax &snax1,const snax &snax2);
inline bool IsAproxEqual(double d1,double d2, double tol=DBL_EPSILON)
	{return(fabs(d1-d2)<tol);}
int CompareSnakeInternalStatus(const vector<bool> & thisVec,bool thisFlipped,
	 const vector<bool> & otherVec, bool otherFlipped);
// Test Function prototypes
int Test_SnakeStructures();
int Test_coordvec();
int Test_snax();
int Test_snaxedge();
int Test_snake();
int Test_snakeinit();
int Test_snakeinit_MC();
int Test_snakeOrderEdges();
int Test_snakeinitflat();
void Test_stepalgo(snake &testSnake,  vector<int> &isImpact);
void Test_stepalgo_mergeclean(snake &testSnake,  vector<int> &isImpact);
// Functions needed at Compile time

// set constructors (used to avoid a variable being unknowingly forgotten)

inline void snax::set(int indexin, double din, double vin, int fromvertin,
	int tovertin, int edgeindin, int isfreezein, int orderedgein)
{
	index=indexin;
	d=din;
	v=vin;
	fromvert=fromvertin; 	// root vertex of *snakemesh
	tovert=tovertin; 	// destination vertex of *snakemesh
	edgeind=edgeindin; 	// edge of snakemesh
	isfreeze=isfreezein; 	// freeze status 
	orderedge=orderedgein;
}


#endif //SNAKSTRUCT_H_INCLUDED