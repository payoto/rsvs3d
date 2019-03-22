#ifndef MAIN_H_INCLUDED 
#define MAIN_H_INCLUDED 



//=================================
// forward declared dependencies
// 		class foo; //when you only need a pointer not the actual object
// 		and to avoid circular dependencies

namespace param {
	class parameters;
}

//=================================
// included dependencies

#include <string>

//=================================
// Code
// 
// This file defines parameters to be used in other parts of the 
// RSVS snaking framework. Default values are defined in "parameters.cpp"
//
// Substructure names are all 4-5 letters

int RSVSExecution(int argc, char* argv[]);
void NoExecution(int execFlow, param::parameters &paramconf);
namespace parse {
	int CommandLineParser(int argc, char* argv[], param::parameters &paramconf);

	namespace config{
		void useconfig(const std::string &confCase, param::parameters &paramconf);
		void loadconfig(const std::string &confCase, param::parameters &paramconf);

	}
}
#endif