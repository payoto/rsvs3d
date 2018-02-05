#include <iostream>
#include <vector>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <unordered_map>

#ifndef ARRAYSTRUCTS_H_INCLUDED
#define ARRAYSTRUCTS_H_INCLUDED

// Defines levels of debugging
#ifdef DEBUGLVL2 // All Debugging calls
#define DEBUGLVL1
#endif

#ifdef DEBUGLVL1 // Debugging of new features.
#define TEST_ARRAYSTRUCT
#endif


using namespace std;
template <class T> // Template for array of structures
class ArrayStruct { 
    public: 
		vector<T> elems;    // vector of elements (structures) 
		unordered_map<int,int> hashTable; // Hash Table of indexed elements
		void HashArray();
		

		T* operator[](const int& a){ // [] Operator returns a pointer to the corresponding elems.
			#ifdef TEST_ARRAYSTRUCT // adds a check in debug mode
				if (a>=elems.size() | 0>elems.size()){
					throw range_error ("Index is out of range");
				}
			#endif
			return(&(elems[a]));
		}
}; 
// member function definition
template<class T> void ArrayStruct <T>::HashArray(){
	hashTable.reserve(elems.size());
	for (int i = 0; i < elems.size(); ++i)
	{
		hashTable.emplace(elems[i].Key(),i);
	}
	cout << "Array Struct Succesfully Hashed" << endl;
}


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
	disp(){
		cout << "Cell #" << index << " Fill " << fill << endl;
	}
	int Key(){
		return(index);
	}
	operator=( cell* otherCell){
		index=otherCell->index;
		fill=otherCell->fill;
		cout << "HELLO"  << " Other #" << otherCell->index << " Fill " << otherCell->fill << endl;
		//return(*this);
	}
};



#endif