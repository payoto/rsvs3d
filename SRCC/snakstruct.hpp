
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
	mesh *snakemesh=NULL;

	// methods
	void disp() const;
	void displight() const;
	bool isready() const ;
	void PrepareForUse();
};


class snax : public meshpart {	
public:
	int index=0;
	double d=0;
	double v=0;
	int fromvert=0; 	// root vertex of *snakemesh
	int tovert=0; 	// destination vertex of *snakemesh
	int edgeind=0; 	// edge of snakemesh
	int isfreeze=0; 	// freeze status 
	int orderedge=0;	// order on the edge

	// interface functions
	void disp() const;
	int Key() const {return (index);};
	void ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu){index+=nVert;};
	void PrepareForUse(){};
	bool isready() const {return(true);}
	/*snax(){
	 index=0;
	 d=0;
	 v=0;
	 fromvert=0; 	// root vertex of *snakemesh
	 tovert=0; 		// destination vertex of *snakemesh
	 edgeind=0; 	// edge of snakemesh
	 isfreeze=0; 	// freeze status 
	 orderedge=0;	// order on the edge
	}*/

};

class snaxedge : public meshpart {
public:
	int index=0;
	int surfind=0;
	coordvec normvector;
	void PrepareForUse();
	void disp() const;
	int Key() const {return (index);}; 
	void ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu){index+=nEdge;};
	bool isready() const {return(normvector.isready());}


};

class snaxsurf : public meshpart {
public: 
	int index=0;
	int voluind=0;
	coordvec normvector;
	void PrepareForUse(); 
	void disp() const;
	int Key() const {return (index);};
	void ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu){index+=nSurf;};
	bool isready() const {return(normvector.isready());}


};

// Function prototypes
int Test_SnakeStructures();
int Test_coordvec();
int Test_snax();
int Test_snaxedge();
int Test_snaxsurf();
int Test_snake();

#endif //SNAKSTRUCT_H_INCLUDED