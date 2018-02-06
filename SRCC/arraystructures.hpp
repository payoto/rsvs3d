
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
typedef unsigned int unsigned_int;
class cell;
//template <class T> class ArrayStruct;

typedef ArrayStruct<cell> cellarray;

// Templates
template <class T> // Template for array of structures
class ArrayStruct { 
    public: 
		vector<T> elems;    // vector of elements (structures) 
		unordered_map<int,int> hashTable; // Hash Table of indexed elements
		void HashArray();
		

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
class cell {
public:
	int index;
	double fill;

	cell(){ // Constructor
		index=0;
		fill=0;
		cout << "Cell #" << index << " Was created" << endl;
	}
	~cell(){ // Destructor
		cout << "Cell #" << index << " Was deleted" << endl;
	}

	void disp(){cout << "Cell #" << index << " Fill " << fill << endl;}
	int Key()    {return(index);}

	void operator=( cell* otherCell){
		index=otherCell->index;
		fill=otherCell->fill;
		cout << "HELLO"  << " Other #" << otherCell->index << " Fill " << otherCell->fill << endl;
	}

};

// functions
int test_arraystructures();

#endif