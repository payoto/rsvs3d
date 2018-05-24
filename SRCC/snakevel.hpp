
//===============================================
// Include Guards
#ifndef SNAKEVEL_H_INCLUDED
#define SNAKEVEL_H_INCLUDED

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
#include "arraystructures.hpp"
#include "snakstruct.hpp" 


//==================================
// Code
// NOTE: function in a class definition are IMPLICITELY INLINED 
//       ie replaced by their code at compile time
using namespace std;


// Base classes


void CalculateSnakeVel(snake &snakein);


#endif // SNAKEVEL_H_INCLUDED

