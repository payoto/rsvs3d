
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
		int prevTime;
		int runTotal;
	public:
		customtest(){
			testCount=0;
			errFlag=0;
			errCount=0;
			runTotal=0;
		}

		void Run(function<int()> test, const char * funcName){
			cout << "-------------------------------------------------------------"
				"---------------------------" << endl;
			cout << "-------------------------------------------------------------"
				"---------------------------" << endl;
			cout << "      Start testing " << funcName << endl;
			cout << "-------------------------------------------------------------"
				"---------------------------" << endl;
			this->prevTime=clock();
			errFlag=test();
			int runTime = clock() - this->prevTime;
			runTotal += runTime;
			++testCount;
			cout << "-------------------------------------------------------------"
				"---------------------------" << endl;
			if (errFlag!=0){
				++errCount;
				cout << "Finished testing " << funcName << endl;
				cout << "                      - Caught Error: " << errFlag << endl;
			} else {
				cout << "Finished testing "<< funcName << endl;
				cout << "                      - No Error" << endl;
			}
			cout << "                Execution time : " << 
				double(runTime)/double(CLOCKS_PER_SEC)*1000 << " ms" << endl; 
			cout << "-------------------------------------------------------------"
				"---------------------------" << endl;
			cout << "-------------------------------------------------------------"
				"---------------------------" << endl;
			cout << endl;
		}
		void PrintSummary(){
			cout << endl;

			cout << "-------------------------------------------------------------"
				"---------------------------" << endl;
			cout << "-------------------------------------------------------------"
				"---------------------------" << endl;
			cout << "      Summary of Tests:" << endl;
			cout << "         " << testCount << " tests completed" << endl;
			cout << "         " << errCount << "  detected errors  " << endl;
			cout << "      Total run time:" << endl;
			cout << "         " << double(this->runTotal)/double(CLOCKS_PER_SEC) 
				<< "  seconds  " << endl;
			cout << "-------------------------------------------------------------"
				"---------------------------" << endl;
			cout << "-------------------------------------------------------------"
				"---------------------------" << endl;
			#ifdef DBG_MEMLEAK
			_CrtDumpMemoryLeaks();
			#endif
		}


};


// member function definition template <class T> 

#endif // POSTPROCESSING_H_INCLUDED

