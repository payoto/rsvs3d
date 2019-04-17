/**
 * Provides a 2D std::vector based container.
 * 
 *@file
 */

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
#include <fstream>

#include "warning.hpp"

//==================================
// Code
// NOTE: function in a class definition are IMPLICITELY INLINED 
//       ie replaced by their code at compile time
using namespace std;
 
// Template Class

/**
 * @brief      Template class for vector of vectors (matrix).
 * 
 * This is designed to be rectangular.
 *
 * @tparam     T     Type of the vector elements.
 */
template<class T> class ArrayVec{
protected:
	vector< vector<T> > elems;
	vector<int> dim;
public:
	void assign(int nR,int nC, T newelem);
	void size(int &nR, int &nC) const {nR=elems.size();if(nR>0){nC=elems[0].size();}};
	int size() const {return(elems.size());}
	void clear(){
		for (int ii=0; ii<int(elems.size());ii++){
			elems[ii].clear();
		}
		elems.clear();
	}
	void write(ofstream &streamout, const char* sep=", ") const;
	vector<T>& operator[](const int a){ 
	// [] Operator returns a reference to the corresponding elems.
		#ifdef SAFE_ACCESS // adds a check in debug mode
		if (((a)>=int(elems.size())) | (0>a)){
			cerr << "Error in " << __PRETTY_FUNCTION__ << endl;
			RSVS3D_ERROR_RANGE(" : Index is out of range");
		}
		#endif //SAFE_ACCESS
		return(elems[a]);
	}

	const vector<T>& operator[](const int a) const { 
	// [] Operator returns a reference to the corresponding elems.
		#ifdef SAFE_ACCESS // adds a check in debug mode
		if (((a)>=int(elems.size())) | (0>a)){
			cerr << "Error in " << __PRETTY_FUNCTION__ << endl;
			RSVS3D_ERROR_RANGE(" : Index is out of range");
		}
		#endif //SAFE_ACCESS
		return(elems[a]);
	}

	~ArrayVec<T>(){
		clear();
	}

};

#include "vectorarray_incl.cpp"

#endif // SNAKEVEL_H_INCLUDED

