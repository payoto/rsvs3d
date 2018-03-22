
//===============================================
// Include Guards
#ifndef SNAKEENGINE_H_INCLUDED
#define SNAKEENGINE_H_INCLUDED

//===============================================
// Levels of debuging Guards
#ifdef DEBUGLVL2 // All Debugging calls
#define DEBUGLVL1

#define 
#endif

#ifdef DEBUGLVL1 // Debugging of new features.
#define TEST_SNAKSTRUCT
#endif

//=================================
// forward declared dependencies
// 		class foo; //when you only need a pointer not the actual object
// 		and to avoid circular dependencies


//=================================
// included dependencies
#include <iostream>
#include <vector>

#include "arraystructures.hpp"

//==================================
// Code
// NOTE: function in a class definition are IMPLICITELY INLINED 
//       ie replaced by their code at compile time
using namespace std;


class snax;
class snaxedge;
class snaxsurf;
class coordvec;
class snaxarray;

typedef ArrayStruct<snaxedge> snaxedgearray;
typedef ArrayStruct<snaxsurf> snaxsurfarray;




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

class snaxarray : public ArrayStruct<snax> 
{
protected:
	unordered_multimap<int,int> hashEdge;
	int isHashEdge=0;
	int isOrderedOnEdge=0;

public: 
	inline int KeyEdge(int a) const ;
	int findedge(int key) const; 
	void HashArrayEdge();
	void ReorderOnEdge();
	// Functions that need modification
	bool checkready();
	void PrepareForUse();
	snax& operator[](const int a){ 
		isOrderedOnEdge=0;
		isHashEdge=0;
		return(ArrayStruct<snax>::operator[](a));
	}
};

class snake  {
private:
	bool is3D=true;
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
	inline void GetMaxIndex(int *nVert,int *nEdge,int *nSurf,int *nVolu) const;
	void HashArray(); // Not really needed as handled by PrepareForUse
	void SetMaxIndex(); // Not really needed as handled by PrepareForUse	
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
	void UpdateCoord();

};


class snax : public meshpart {	
public:
	int index=0;
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


};


class snaxedge : public meshpart {
public:
	int index=0;
	int surfind=0;
	coordvec normvector;
	void PrepareForUse();
	void disp() const;
	int Key() const {return (index);}; 
	void ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu);
	void ChangeIndicesSnakeMesh(int nVert,int nEdge,int nSurf,int nVolu);
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wunused-parameter"
	bool isready(bool isInMesh) const {return(normvector.isready());}
	#pragma GCC diagnostic pop
	void read(FILE * fid);
	void write(FILE *fid)const;


};

class snaxsurf : public meshpart {
public: 
	int index=0;
	int voluind=0;
	coordvec normvector;
	void PrepareForUse(); 
	void disp() const;
	int Key() const {return (index);};
	void ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu);
	void ChangeIndicesSnakeMesh(int nVert,int nEdge,int nSurf,int nVolu);
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wunused-parameter"
	bool isready(bool isInMesh) const {return(normvector.isready());}
	#pragma GCC diagnostic pop
	void read(FILE * fid);
	void write(FILE * fid)const;


};

// Function prototypes


// Test Function prototypes
int Test_SnakeStructures();
int Test_coordvec();
int Test_snax();
int Test_snaxedge();
int Test_snake();
int Test_snakeinit();
int Test_snakeOrderEdges();
int Test_snakeinitflat();

// Functions needed at Compile time

// set constructors (used to avoid a variable being unknowingly forgotten)
inline void snax::set(int indexin, double din,double vin,int fromvertin,int tovertin,
	int edgeindin,int isfreezein,int orderedgein)
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

inline int snaxarray::KeyEdge(int a) const 
{
	return((*this)(a)->edgeind);
}

#endif //SNAKEENGINE_H_INCLUDED