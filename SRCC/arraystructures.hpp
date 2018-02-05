#include <iostream>
#include <vector>
#include <cstdlib>
#include <string>
#include <stdexcept>

#ifndef ARRAYSTRUCTS_H_INCLUDED
#define ARRAYSTRUCTS_H_INCLUDED

using namespace std;

template <class T> // Template for array of structures
class ArrayStruct { 
   public: 
      vector<T> elems;    // vector of elements (structures) 

      // Hash Table
}; 

class ArrayStruct2 {
public:
	ArrayStruct2(){ // Constructor

		
	}
	~ArrayStruct2(){ // Destructor

		
	}


};

class cell {
public:
	int index;
	double fill;
	cell(){ // Constructor
		index=0;
		fill=0;
	}
	~cell(){ // Destructor
		cout << "Cell #" << index << " Was deleted" << endl;
		
	}
	disp(){
		cout << "Cell #" << index << " Fill " << fill << endl;
	}
};

typedef struct {
	int index;
	int cellind[2];
	int vertex[2];
	int orientation; /* 1 is vertical 0 is Horizontal */
} edge;

typedef struct {
	int index;
	double coord[2]; /* needs to be changed wiht dim */
	
} vertex;


#endif