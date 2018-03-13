
//===============================================
// Include Guards
#ifndef SNAKSTRUCT_H_INCLUDED
#define SNAKSTRUCT_H_INCLUDED

//===============================================
// Levels of debuging Guards
#ifdef DEBUGLVL2 // All Debugging calls
#define DEBUGLVL1

#define 
#endif

#ifdef DEBUGLVL1 // Debugging of new features.
#define TEST_SNAKSTRUCT
#define TEST_RANGE
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

typedef ArrayStruct<snax> snaxarray;
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
	inline double GetNorm();
	coordvec Unit() const ;
	double Unit(const int a) const;
	void assign(double a, double b, double c);
	double& operator[](int a);
	double operator()(int a) const;
	void disp() const;

	coordvec(){
		elems.reserve(3); // reserves 3 as this is the size of the array
		elems.assign(3,0);
		norm=0;
		isuptodate=0;
	}
};


class snake : public mesh {
private:

public:
	// Handling of data specific to snake
	snaxarray snaxs; 			// properties associated with snakconn verts
	snaxedgearray snaxedges; 	// properties associated with snakconn edges
	snaxsurfarray snaxsurfs; 	// properties associated with snakconn surfs
	// Using the mesh container to store connectivity
	mesh snakeconn;
	// pointer to snaking mesh
	mesh *snakemesh;

};


class snax : public meshpart {	
public:
	int index;
	double d;
	double v;
	int fromvert; 	// root vertex of *snakemesh
	int tovert; 	// destination vertex of *snakemesh
	int edgeind; 	// edge of snakemesh
	int isfreeze; 	// freeze status 
	int orderedge;	// order on the edge

	// interface functions
	void disp() const;
	int Key() const {return (index);};
	void ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu){index+=nVert;};


};

class snaxedge : public meshpart {
public:
	int index;
	coordvec normvector;
	void disp() const;
	int Key() const {return (index);}; 
	void ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu){index+=nEdge;};

};

class snaxsurf : public meshpart {
public: 
	int index;
	coordvec normvector;
	void disp() const;
	int Key() const {return (index);};
	void ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu){index+=nSurf;};

};

// Function prototypes
int Test_SnakeStructures();
int Test_coordvec();

#endif //SNAKSTRUCT_H_INCLUDED