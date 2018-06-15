
//===============================================
// Include Guards
#ifndef VECTORARRAY_H_INCLUDED
#define VECTORARRAY_H_INCLUDED

//===============================================
// Levels of debuging Guards
#ifdef DEBUGLVL2 // All Debugging calls
#define DEBUGLVL1
#define TEST_ALL
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


//==================================
// Code
// NOTE: function in a class definition are IMPLICITELY INLINED 
//       ie replaced by their code at compile time
using namespace std;

// Template Class

template<class T> class ArrayVec{
protected:
	vector< vector<T> > elems;
	vector<int> dim;
public:
	void assign(int nR,int nC, T& newelem);

	vector<T>& operator[](const int a){ 
	// [] Operator returns a reference to the corresponding elems.
		#ifdef SAFE_ACCESS // adds a check in debug mode
		if (((a)>=int(elems.size())) | (0>a)){
			cerr << "Error in " << __PRETTY_FUNCTION__ << endl;
			throw range_error (" : Index is out of range");
		}
		#endif //SAFE_ACCESS
		return(elems[a]);
	}

};

#include "vectorarray_incl.cpp"

#endif // SNAKEVEL_H_INCLUDED

