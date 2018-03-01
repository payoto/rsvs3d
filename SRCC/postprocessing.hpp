
//===============================================
// Include Guards
#ifndef POSTPROCESSING_H_INCLUDED
#define POSTPROCESSING_H_INCLUDED

//===============================================
// Levels of debuging Guards
#ifdef DEBUGLVL2 // All Debugging calls
#define DEBUGLVL1

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
#include "arraystructures.hpp"


//==================================
// Code
// NOTE: function in a class definition are IMPLICITELY INLINED 
//       ie replaced by their code at compile time
using namespace std;


// Base classes

class tecplotfile {
	public:
		FILE *fid=NULL;
		int OpenFile();
		int PrintMesh(mesh *meshout);

}

// Derived Classes

// functions



// member function definition template <class T> 



#endif // POSTPROCESSING_H_INCLUDED

