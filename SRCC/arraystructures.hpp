
//===============================================
// Include Guards
#ifndef ARRAYSTRUCTS_H_INCLUDED
#define ARRAYSTRUCTS_H_INCLUDED

//===============================================
// Levels of debuging Guards
#ifdef DEBUGLVL2 // All Debugging calls
#define DEBUGLVL1
#define TEST_ARRAYSTRUCT
#define 
#endif

#ifdef DEBUGLVL1 // Debugging of new features.
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
#include <cstdlib>
#include <string>
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
public: 
	vector<T> elems;    // vector of elements (structures) 
	unordered_map<int,int> hashTable; // Hash Table of indexed elements
	void HashArray();
	void disp();
	void Concatenate(ArrayStruct<T> *other);
	void SetMaxIndex();
	void PopulateIndices();
	inline int GetMaxIndex();

	T* operator[](const int& a){ // [] Operator returns a pointer to the corresponding elems.
		#ifdef TEST_RANGE // adds a check in debug mode
		if ((unsigned_int(a)>=elems.size()) | (0>a)){
			throw range_error ("Index is out of range");
		}
		#endif
		return(&(elems[a]));
	}

	void Init(int n){
		T sT;
		elems.reserve(n);
		elems.assign(n,sT);
	}
	int size(){
		return(int(elems.size()));
	}
	int find(int key){
		auto search=hashTable.find(key);
		if (search==hashTable.end()){
			return(-1);
		}
		return(search->second);
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
	void GetMaxIndex(int *nVert,int *nEdge,int *nSurf,int *nVolu);
	void Init(int nVe,int nE, int nS, int nVo);
	void disp();
	void Concatenate(mesh *other);
// Mesh merging
	void MakeCompatible_inplace(mesh *other);
	mesh MakeCompatible(mesh other);
	void ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu);
};


class meshpart{ // Abstract class to ensure interface is correct
	public : 
	int index;
	virtual void disp() =0 ;
	virtual int Key() =0 ; 
	virtual void ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu)=0 ;

};

// Derived Classes
class volu: public meshpart {
public:
	int index;
	double fill,target,error;
	vector<int> surfind;

	void ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu);
	void disp();

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
	void operator=( volu* other){
		index=other->index;
		fill=other->fill;
		target=other->target;
		error=other->error;
		surfind=other->surfind;

		#ifdef TEST_ARRAYSTRUCT
		cout << "OTHER: " ; other->disp();
		#endif
	}

	int Key() {return(index);}
};



class surf: public meshpart {
public:
	int index;
	double fill,target,error;
	vector<int> edgeind;
	vector<int> voluind;
	 // reserves 2 as this is the size of the array

	void disp();
	void ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu);

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
	void operator=( surf* other){
		index=other->index;
		fill=other->fill;
		error=other->error;
		target=other->target;
		edgeind=other->edgeind;
		voluind=other->voluind;
	}

	int Key() {return(index);}

};


class edge: public meshpart {
public:
	int index;

	vector<int> vertind;
	vector<int> surfind;
	 // reserves 2 as this is the size of the array

	void ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu);
	void disp();

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
	void operator=( edge* other){
		index=other->index;

		vertind=other->vertind;
		surfind=other->surfind;
	}

	int Key() {return(index);}
};

class vert: public meshpart {
public:
	int index;

	vector<int> edgeind;
	vector<double> coord;
	 // reserves 2 as this is the size of the array

	void disp();
	void ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu);

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
	void operator=( vert* other){
		index=other->index;

		edgeind=other->edgeind;
		coord=other->coord;
	}

	int Key() {return(index);}
};

// functions

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

// member function definition template <class T> : "ArrayStruct"
template<class T> void ArrayStruct <T>::HashArray(){
hashTable.reserve(elems.size());
for (int i = 0; i < int(elems.size()); ++i)
{
	hashTable.emplace(elems[i].Key(),i);
}
   //cout << "Array Struct Succesfully Hashed" << endl;
}

template<class T> void ArrayStruct <T>::disp(){
for (int ii=0 ; unsigned_int(ii)<elems.size();ii++){
	cout << "Array " << ii << " " ;
	elems[ii].disp();
}
}

template<class T> void ArrayStruct <T>::Concatenate(ArrayStruct <T> *other){
int nCurr,nNew,nTot,ii;
nCurr=this->size();
nNew=other->size();
nTot=nCurr+nNew;
elems.reserve(nTot);

for (ii=0; ii<nNew;ii++){
	elems.push_back(other->elems[ii]);
}
this->HashArray();
}


template<class T> void ArrayStruct <T>::SetMaxIndex(){
int n=elems.size();
maxIndex=-1;
for(int ii=0;ii<n;++ii){
	maxIndex= (maxIndex<this->elems[ii].index) ? elems[ii].index : maxIndex;
}
}

template<class T> inline int ArrayStruct <T>::GetMaxIndex(){
return(maxIndex);
}

template<class T> void ArrayStruct <T>::PopulateIndices(){
int n=elems.size();
maxIndex=-1;
for(int ii=0;ii<n;++ii){
	elems[ii].index=ii+1;
}
maxIndex=n;
}
#endif // ARRAYSTRUCTS_H_INCLUDED