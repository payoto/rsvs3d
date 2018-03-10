
//===============================================
// Include Guards
#ifndef SNAKSTRUCT_H_INCLUDED
#define SNAKSTRUCT_H_INCLUDED

//===============================================
// Levels of debuging Guards
#ifdef DEBUGLVL2 // All Debugging calls
#define DEBUGLVL1
#define TEST_SNAKSTRUCT
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
#include "arraystruct.hpp"

//==================================
// Code
// NOTE: function in a class definition are IMPLICITELY INLINED 
//       ie replaced by their code at compile time
using namespace std;

typedef ArrayStruct<snax> snaxarray;



class snake : public mesh {
private:

public:
	snaxarray snaxs;

};


class snax : public meshpart {


};