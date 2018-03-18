
//===============================================
// Include Guards
#ifndef ARRAYSTRUCTS_H_INCLUDED
#define ARRAYSTRUCTS_H_INCLUDED

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
#ifdef DBG_MEMLEAK
#define _CRTDBG_MAP_ALLOC  
#include <stdlib.h>  
#include <crtdbg.h>
#endif //DBG_MEMLEAK

#include <iostream>
#include <vector>
#include <cstdlib>
#include <string>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <functional>

//==================================
// Code
// NOTE: function in a class definition are IMPLICITELY INLINED 
//       ie replaced by their code at compile time
using namespace std;


class meshpart;
class mesh;

class volu;
class surf;
class vert;
class edge;

template <class T> class ArrayStruct;
typedef unsigned int unsigned_int;

typedef ArrayStruct<volu> voluarray;
typedef ArrayStruct<surf> surfarray;
typedef ArrayStruct<edge> edgearray;
typedef ArrayStruct<vert> vertarray;

// Templates
template <class T> 
class ArrayStruct { 
protected: 
	int maxIndex;
	int isHash=0;
	int isSetMI=0;
	bool readyforuse=false;
	vector<T> elems;    // vector of elements (structures) 
	unordered_map<int,int> hashTable; // Hash Table of indexed elements

public: 
	void disp() const;
	int find(int key) const ;
	vector<int> find_list(vector<int> key) const ;
	inline int GetMaxIndex() const;
	inline void Init(int n);
	bool isready() const {return(readyforuse);};
	bool checkready() ;
	void Concatenate(const ArrayStruct<T> &other);
	void PopulateIndices();
	void SetMaxIndex();
	void HashArray();
	void PrepareForUse();
	void ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu);
	// methods needed from vector
	inline int size() const;
	inline int capacity() const;
	inline void assign(int n, T& newelem);
	inline void push_back(T& newelem);
	inline void reserve(int n);
	// Operators
	const T* operator()(const int a) const{ 
	// () Operator returns a constant pointer to the corresponding elems.
	// Cannot be used on the left hand side and can't be used to edit data in elems
		#ifdef SAFE_ACCESS // adds a check in debug mode
		if ((unsigned_int(a)>=elems.size()) | (0>a)){
			throw range_error ("Index is out of range");
		}
		#endif //SAFE_ACCESS
		return(&(elems[a]));
	}
	T& operator[](const int a){ 
	// [] Operator returns a reference to the corresponding elems.
		#ifdef SAFE_ACCESS // adds a check in debug mode
		if ((unsigned_int(a)>=elems.size()) | (0>a)){
			throw range_error ("Index is out of range");
		}
		#endif //SAFE_ACCESS
		isHash=0;
		isSetMI=0;
		readyforuse=false;
		return(elems[a]);
	}

}; 
  


// Base class

class mesh {
private:

public:
	vertarray verts;
	edgearray edges;
	surfarray surfs;
	voluarray volus;
// basic operations grouped from each field
	void HashArray();
	void SetMaxIndex();
	void GetMaxIndex(int *nVert,int *nEdge,int *nSurf,int *nVolu) const;
	void Init(int nVe,int nE, int nS, int nVo);
	void PrepareForUse();
	void disp() const;
	void displight() const;
	void Concatenate(const mesh &other);
	bool isready() const;
	void PopulateIndices();
// Mesh merging
	void MakeCompatible_inplace(mesh &other) const;
	mesh MakeCompatible(mesh other) const;
	void ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu);
};


class meshpart{ // Abstract class to ensure interface is correct
	public : 
	int index=0;
	virtual void disp() const =0 ;
	virtual int Key() const =0 ; 
	virtual void ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu)=0 ;
	virtual void PrepareForUse()=0 ;
	virtual bool isready() const=0 ;
	//virtual operator=( meshpart* other)=0 ;

};

// Derived Classes
class volu: public meshpart {
public:
	int index;
	double fill,target,error;
	vector<int> surfind;

	void ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu);
	void disp() const;
	void PrepareForUse(){};
	bool isready() const {return(true);}

	volu(){ // Constructor
		index=0;
		fill=0;
		target=1;
		error=1;

		#ifdef TEST_ARRAYSTRUCT
		cout << "volu #" << index << " Was created " << surfind.size() << endl;
		#endif
	}
	volu(const volu& oldVolu){ // Copy-Constructor
		index=oldVolu.index;
		fill=oldVolu.fill;
		target=oldVolu.target;
		error=oldVolu.error;
		surfind=oldVolu.surfind;

		#ifdef TEST_ARRAYSTRUCT
		cout << "copyvolu #" << index << " Was created " << surfind.size() << endl;
		#endif
	}
	~volu(){ // Destructor
		surfind.clear();

		#ifdef TEST_ARRAYSTRUCT
		cout << "volu #" << index << " Was deleted " << surfind.size() << endl;
		#endif

	}
	void operator=(const volu* other){
		index=other->index;
		fill=other->fill;
		target=other->target;
		error=other->error;
		surfind=other->surfind;

		#ifdef TEST_ARRAYSTRUCT
		cout << "OTHER: " ; other->disp();
		#endif
	}

	int Key() const {return(index);}
};



class surf: public meshpart {
public:
	int index;
	double fill,target,error;
	vector<int> edgeind;
	vector<int> voluind;
	 // reserves 2 as this is the size of the array

	void disp() const;
	void ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu);
	void PrepareForUse(){};	
	bool isready() const {return(true);}

	surf(){ // Constructor
		index=0;
		fill=1;
		target=1;
		error=1;
		voluind.reserve(2); // reserves 2 as this is the size of the array
		voluind.assign(2,0);
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
	}
	void operator=(const surf* other){
		index=other->index;
		fill=other->fill;
		error=other->error;
		target=other->target;
		edgeind=other->edgeind;
		voluind=other->voluind;
	}

	int Key() const {return(index);}

};


class edge: public meshpart {
public:
	int index;

	vector<int> vertind;
	vector<int> surfind;
	 // reserves 2 as this is the size of the array

	void ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu);
	void disp() const;
	void PrepareForUse(){};
	bool isready() const {return(true);}

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

class vert: public meshpart {
public:
	int index;

	vector<int> edgeind;
	vector<double> coord;
	 // reserves 2 as this is the size of the array

	void disp() const;
	void ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu);
	void PrepareForUse(){};
	bool isready() const {return(true);}

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

// functions
template <class T> bool CompareDisp(T *mesh1,T *mesh2);
template<class T> int TestReadyness(T &stackT, const char* txt, bool errTarg);
bool CompareFuncOut(function<void()> func1, function<void()> func2);
//test functions
int Test_ArrayStructures();
int Test_Volu();
int Test_Surf();
int Test_Vert();
int Test_Edge();
int Test_Mesh();
void PopulateIndices(mesh *meshin);
//template <class T> bool CompareDisp(T *mesh1,T *mesh2);
//bool CompareFuncOut(function<void()> mesh1, function<void()> mesh2);

#include "arraystructures_incl.cpp"

#endif // ARRAYSTRUCTS_H_INCLUDED