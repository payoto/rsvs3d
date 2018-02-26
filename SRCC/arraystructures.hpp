
//===============================================
// Include Guards
#ifndef ARRAYSTRUCTS_H_INCLUDED
#define ARRAYSTRUCTS_H_INCLUDED

//===============================================
// Levels of debuging Guards
#ifdef DEBUGLVL2 // All Debugging calls
#define DEBUGLVL1
#endif

#ifdef DEBUGLVL1 // Debugging of new features.
#define TEST_ARRAYSTRUCT
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

//==================================
// Code
// NOTE: function in a class definition are IMPLICITELY INLINED 
//       ie replaced by their code at compile time
using namespace std;

class volu;
template <class T> class ArrayStruct;
typedef unsigned int unsigned_int;

typedef ArrayStruct<volu> voluarray;

// Templates
template <class T> // Template for array of structures
class ArrayStruct { 
    public: 
		vector<T> elems;    // vector of elements (structures) 
		unordered_map<int,int> hashTable; // Hash Table of indexed elements
		void HashArray();
		void disp2();
		

		T* operator[](const int& a){ // [] Operator returns a pointer to the corresponding elems.
			#ifdef TEST_ARRAYSTRUCT // adds a check in debug mode
				if ((unsigned_int(a)>=elems.size()) | (0>a)){
					throw range_error ("Index is out of range");
				}
			#endif
			return(&(elems[a]));
		}



}; 


// Classes
class volu {
public:
	int index;
	double fill;
	vector<int> surfind;

	volu(){ // Constructor
		index=0;
		fill=0;

		#ifdef TEST_ARRAYSTRUCT
			cout << "volu #" << index << " Was created " << surfind.size() << endl;
		#endif
	}
	volu(const volu& oldVolu){ // Copy-Constructor
		index=oldVolu.index;
		fill=oldVolu.fill;
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
		surfind=other->surfind;

		#ifdef TEST_ARRAYSTRUCT
			cout << "OTHER: " ; other->disp();
		#endif
	}

	void disp(){
		cout << "volu : index " << index << " | fill " << fill << " | surfind " << surfind.size();
		for (int i=0; unsigned_int(i)<surfind.size();i++){
			cout << "-" << surfind[i];
		}
		cout << endl;
	}
	int Key()    {return(index);}


};


class surf {
public:
	int index;
	double fill;
	vector<int> edgeind;
	vector<int> voluind;

	surf(){ // Constructor
		index=0;
		fill=0;
		voluind.reserve(2); // reserves 2 as this is the size of the array
	}
	~surf(){ // Destructor

		edgeind.clear();
		voluind.clear();
	
	}

	void disp(){cout << "surf #" << index << " Fill " << fill ;}
	int Key()    {return(index);}

	void operator=( surf* other){
		index=other->index;
		fill=other->fill;
		edgeind=other->edgeind;
		voluind=other->voluind;
	}

};
// functions
int test_arraystructures();
int test_volu();
#endif