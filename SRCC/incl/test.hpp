
//===============================================
// Include Guards
#ifndef TEST_H_INCLUDED
#define TEST_H_INCLUDED

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
#include <stdarg.h>
#include <ctime>
#include <fstream>

//==================================
// Code
// NOTE: function in a class definition are IMPLICITELY INLINED 
//       ie replaced by their code at compile time
using namespace std;


// Base classes

class customtest {
	private:
		int testCount;
		int errFlag;
		int errCount;
		int unhandledError;
		int prevTime;
		int runTotal;
		int lastRunTime;
	public:
		customtest(){
			testCount=0;
			errFlag=0;
			errCount=0;
			unhandledError=0;
			runTotal=0;
			lastRunTime=0;
			cout << "-------------------------------------------------------------"
				"---------------------------" << endl;
			cout << "-------------------------------------------------------------"
				"---------------------------" << endl;
			cout << "      Start testing " << endl;
			cout << "-------------------------------------------------------------"
				"---------------------------" << endl;
		}
		~customtest(){
			this->PrintSummary();
		}

		int Run(function<int()> test, const char * funcName);
		int RunSilent(function<int()> test, const char * funcName);
		void PrintSummary();


};


// member function definition template <class T> 

#endif // POSTPROCESSING_H_INCLUDED

