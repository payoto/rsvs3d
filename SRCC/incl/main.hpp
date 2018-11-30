#ifndef MAIN_H_INCLUDED 
#define MAIN_H_INCLUDED 



//=================================
// forward declared dependencies
// 		class foo; //when you only need a pointer not the actual object
// 		and to avoid circular dependencies

//=================================
// included dependencies


//=================================
// Code
// 
// This file defines parameters to be used in other parts of the 
// RSVS snaking framework. Default values are defined in "parameters.cpp"
//
// Substructure names are all 4-5 letters

void CommandLineParser(int argc, char* argv[]);

#endif