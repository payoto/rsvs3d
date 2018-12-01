
#ifndef WARNING_H_INCLUDED
#define WARNING_H_INCLUDED

//===============================================
// Levels of debuging Guards
#ifdef DEBUGLVL2 // All Debugging calls
#define DEBUGLVL1

#endif

#ifdef DEBUGLVL1 // Debugging of new features.
#define TEST_WARNING
#endif

//=================================
// forward declared dependencies
// 		class foo; //when you only need a pointer not the actual object
// 		and to avoid circular dependencies



//=================================
// included dependencies
#include <iostream>
#include <stdarg.h>
#include <stdexcept>
#include <fstream>

//==================================
// Code
// NOTE: function in a class definition are IMPLICITELY INLINED 
//       ie replaced by their code at compile time

using namespace std;

// class userwarning {
// private:
// 	const char* message;

// public:
// 	void SetMessage(const char* message)=0;
// 	void ThrowME(){
// 		cerr << message <<endl;
// 	};
// };


void ThrowWarning(const char * message);
template<class T>
void CheckFStream(const T &file, const char* callerID, 
	const std::string &fileName){
	if (!file.is_open()){
		std::string errstr;
		errstr = "Parameter file failed to open ";
		errstr += "in " ;
		errstr += callerID;
		errstr += " \n:  " + fileName;
		throw std::invalid_argument(errstr.c_str());
	}
}
#endif // WARNING_H_INCLUDED