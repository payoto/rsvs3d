/**
 * Provides the custom testing system used by the RSVS3D project.
 *  
 *@file
 */
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
#include <string>

//==================================
// Code
// NOTE: function in a class definition are IMPLICITELY INLINED 
//       ie replaced by their code at compile time


/**
 * @brief      Namespace for rsvs tests.
 */
namespace rsvstest {

		/**
	 * @brief      Class for customtest.
	 */
	class customtest {
		private:
			int testCount;
			int errFlag;
			int errCount;
			int unhandledError;
			int prevTime;
			int runTotal;
			int lastRunTime;
			std::string testName;
		public:
			customtest(const char * testNameIn=""){
				testCount=0;
				errFlag=0;
				errCount=0;
				unhandledError=0;
				runTotal=0;
				lastRunTime=0;
				testName=testNameIn;
				std::cout << "-------------------------------------------------------------"
					"---------------------------" << std::endl;
				std::cout << "-------------------------------------------------------------"
					"---------------------------" << std::endl;
				std::cout << "      Start testing " << testName << std::endl;
				std::cout << "-------------------------------------------------------------"
					"---------------------------" << std::endl;
			}
			int ReturnErrCount(){return errCount;}
			~customtest(){
				this->PrintSummary();
			}

			int Run(std::function<int()> test, const char * funcName,
				int expectedTime=-1);
			int RunSilent(std::function<int()> test, const char * funcName, 
				int expectedTime=-1);
			void PrintSummary();
	};

	int maintest();
	int shorttest();
	int newtest();
	int arraystructtemplates();
	int meshprocesses();
	int snakeprocesses();
	int longevolution();
	int RSVSprocesses();
	int RSVS2Dprocesses();
	int tetgenprocesses();
	int JSONprocesses();
	int integrationprocesses();
}

// member function definition template <class T> 

#endif // POSTPROCESSING_H_INCLUDED

