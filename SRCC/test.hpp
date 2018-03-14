
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
	public:
		customtest(){
			testCount=0;
			errFlag=0;
			errCount=0;
		}

		void Run(function<int()> test, const char * funcName){
			cout << "----------------------------------------------------------------------------------------" << endl;
			cout << "----------------------------------------------------------------------------------------" << endl;
			cout << "      Start testing " << funcName << endl;
			cout << "----------------------------------------------------------------------------------------" << endl;
			errFlag=test();
			++testCount;
			cout << "----------------------------------------------------------------------------------------" << endl;
			if (errFlag!=0){
				++errCount;
				cout << "Finished testing " << funcName << endl;
				cout << "                      - Caught Error: " << errFlag << endl;
			} else {
				cout << "Finished testing "<< funcName << endl;
				cout << "                      - No Error" << endl;
			}
			cout << "----------------------------------------------------------------------------------------" << endl;
			cout << "----------------------------------------------------------------------------------------" << endl;
			cout << endl;
		}
		void PrintSummary(){
			cout << endl;

			cout << "----------------------------------------------------------------------------------------" << endl;
			cout << "----------------------------------------------------------------------------------------" << endl;
			cout << "      Summary of Tests:" << endl;
			cout << "         " << testCount << " tests completed" << endl;
			cout << "         " << errCount << "  detected errors  " << endl;
			cout << "----------------------------------------------------------------------------------------" << endl;
			cout << "----------------------------------------------------------------------------------------" << endl;
		}


};


// member function definition template <class T> 

#endif // POSTPROCESSING_H_INCLUDED

