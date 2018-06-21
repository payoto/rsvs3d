//===============================================
// Include Guards
#ifndef MESHREFINEMENT_H_INCLUDED
#define MESHREFINEMENT_H_INCLUDED

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

#include "mesh.hpp"

//==================================
// Code
// NOTE: function in a class definition are IMPLICITELY INLINED 
//       ie replaced by their code at compile time
using namespace std;




// Forward declared templated functions

// Base class


//test functions

#endif // MESH_H_INCLUDED