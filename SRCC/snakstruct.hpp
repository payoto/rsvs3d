
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

#include "arraystructures.hpp"
#include "postprocessing.hpp"

//==================================
// Code
// NOTE: function in a class definition are IMPLICITELY INLINED 
//       ie replaced by their code at compile time
using namespace std;

template <class T> class SnakStruct; 
class snax;
class snaxedge;
class snaxsurf;
class coordvec;
class snaxarray;

typedef SnakStruct<snaxedge> snaxedgearray;
typedef SnakStruct<snaxsurf> snaxsurfarray;


template <class T> 
class SnakStruct : public ArrayStruct<T> 
{
protected:
	using ArrayStruct<T>::elems;
    using ArrayStruct<T>::readyforuse;

	unordered_multimap<int,int> hashParent;
	int isHashParent=0;

public: 
	friend class snake;

	inline int KeyParent(int a) const ;
	int findparent(int key) const; 
	void findsiblings(int key, vector<int> &siblings) const{
		siblings=ReturnDataEqualRange(key, hashParent);}
	int countparent(int key) const {return(hashParent.count(key));};
	void HashParent();

	// Functions that need modification
	bool checkready();
	void ForceArrayReady();
	void PrepareForUse();
	void Concatenate(const SnakStruct<T> &other);
	T& operator[](const int a){ 
		isHashParent=0;
		return(ArrayStruct<T>::operator[](a));
	}

};

class snaxarray : public SnakStruct<snax> 
{
protected:
	int isOrderedOnEdge=0;

public: 
	friend class snake;


	void ReorderOnEdge();
	void OrderOnEdge();
	void CalculateTimeStepOnEdge(vector<double> &dt, 	vector<bool> &isSnaxDone, int edgeInd);
	void DetectImpactOnEdge(vector<int> &isImpact, vector<bool> &isSnaxDone, int edgeInd);
	// Functions that need modification
	bool checkready();
	void ForceArrayReady();
	void PrepareForUse();
	void Concatenate(const snaxarray &other);
	snax& operator[](const int a){ 
		isOrderedOnEdge=0;

		return(SnakStruct<snax>::operator[](a));
	}

};

class coordvec {
	// Handles the use and norm of a vector for which the norm and the unit 
	// value might be needed
protected:
	vector<double> elems;
	double norm;
	int isuptodate;
	

public:
	double CalcNorm();
	double GetNorm();
	void PrepareForUse();
	coordvec Unit() const ;
	double Unit(const int a) const;
	void assign(double a, double b, double c);
	double& operator[](int a);
	double operator()(int a) const;
	void disp() const;
	bool isready() const {return(bool(isuptodate));};

	coordvec(){
		elems.reserve(3); // reserves 3 as this is the size of the array
		elems.assign(3,0);
		norm=0;
		isuptodate=0;
		#ifdef TEST_SNAKSTRUCT
		cout << "constructor called for coordvec" << endl;
		#endif
	}
};


class snake  {
private:
	bool is3D=true;
	void SetLastIndex(); // Baaaad function do not use if you're not sure.
public:
	// Handling of data specific to snake
	snaxarray snaxs; 			// properties associated with snakconn verts
	snaxedgearray snaxedges; 	// properties associated with snakconn edges
	snaxsurfarray snaxsurfs; 	// properties associated with snakconn surfs
	// Using the mesh container to store connectivity
	mesh snakeconn;
	// pointer to snaking mesh
	mesh *snakemesh=NULL;

	// basic operations grouped from each field
	void disp() const;
	void displight() const;
	bool isready() const ;
	void PrepareForUse();
	void Init(mesh *snakemesh,int nSnax, int nEdge, int nSurf, int nVolu);
	void reserve(int nSnax, int nEdge, int nSurf, int nVolu);
	inline void GetMaxIndex(int *nVert,int *nEdge,int *nSurf,int *nVolu) const;
	void HashArray(); // Not really needed as handled by PrepareForUse
	void SetMaxIndex(); // Not really needed as handled by PrepareForUse
	void SetMaxIndexNM(); // SetMaxIndex no mesh
	void Concatenate(const snake &other);
	bool Check3D(){return(is3D);}
	// Snake merging
	void MakeCompatible_inplace(snake &other) const;
	snake MakeCompatible(snake other) const;
	void ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu);
	void ChangeIndicesSnakeMesh(int nVert,int nEdge,int nSurf,int nVolu);

	// Snake Movement
	void UpdateDistance(double dt);
	void UpdateDistance(const vector<double> &dt);
	void CalculateTimeStep(vector<double> &dt, double dtDefault);
	void SnaxImpactDetection(vector<int> &isImpact);
	void UpdateCoord();
	void Flip(); // reverses snake directions

};

class snakpart { // required functions for parts of snake
public: 
	virtual int KeyParent() const =0;
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

};


class snaxedge : public meshpart , public snakpart {
public:
	int surfind=0;
	coordvec normvector;

	void PrepareForUse();
	void disp() const;
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
inline bool IsAproxEqual(double d1,double d2) {return(fabs(d1-d2)<DBL_EPSILON);}
// Test Function prototypes
int Test_SnakeStructures();
int Test_coordvec();
int Test_snax();
int Test_snaxedge();
int Test_snake();
int Test_snakeinit();
int Test_snakeOrderEdges();
int Test_snakeinitflat();
void Test_stepalgo(snake &testSnake, vector<double> dt, vector<int> isImpact, tecplotfile &outSnake);

// Functions needed at Compile time

// set constructors (used to avoid a variable being unknowingly forgotten)


#include "snakstruct_incl.cpp"

#endif //SNAKSTRUCT_H_INCLUDED