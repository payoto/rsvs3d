
#include <iostream>
#include <fstream>
#include <ctime>
#include "warning.hpp"


using namespace std;

#ifndef TIME_EXEC
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored  "-Wunused-parameter"
#endif
int rsvs3d::TimeStamp(const char* str,int start_s){
	int stop_s=clock();
	#ifdef TIME_EXEC
	cout << str << " " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << "ms; ";
	#endif
	return(stop_s);
}
#ifndef TIME_EXEC
	#pragma GCC diagnostic pop
#endif

void ThrowWarning(const char * message){
	cerr << message << endl;
}

