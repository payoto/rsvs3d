#include <iostream>
#include <cstdlib>
#include <string>

#include "arraystructures.hpp"


int main(){
	int errCount,testCount,errFlag;

	errCount=0;
	testCount=0;

	cout << "--------------------------------------------" << endl;
	cout << "--------------------------------------------" << endl;
	cout << "      Start testing arraystructures" << endl;
	cout << "--------------------------------------------" << endl;
	errFlag=test_arraystructures();
	++testCount;
	cout << "--------------------------------------------" << endl;
	if (errFlag!=0){
		++errCount;
		cout << "Finished testing arraystructures" << endl;
		cout << "                      - Caught Error: " << errFlag << endl;
	} else {
		cout << "Finished testing arraystructures" << endl;
		cout << "                      - No Error" << endl;
	}
	cout << "--------------------------------------------" << endl;
	cout << "--------------------------------------------" << endl;
	cout << endl;




	
	
	cout << endl;

	cout << "--------------------------------------------" << endl;
	cout << "--------------------------------------------" << endl;
	cout << "      Summary of Tests:" << endl;
	cout << "         " << testCount << " tests completed" << endl;
	cout << "         " << errCount << "  detected errors  " << endl;
	cout << "--------------------------------------------" << endl;
	cout << "--------------------------------------------" << endl;


}



