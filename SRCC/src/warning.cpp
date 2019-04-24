
#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include "warning.hpp"


using namespace std;
double rsvs3d::Clock2ms(int clockCycles)
{
	return(double(clockCycles)/double(CLOCKS_PER_SEC)*1000.0);
}

#ifndef TIME_EXEC
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored  "-Wunused-parameter"
#endif
int rsvs3d::TimeStamp(const char* str,int start_s){
	int stop_s=clock();
	#ifdef TIME_EXEC
	if(str!=NULL){
		cout << str << " " << Clock2ms(stop_s-start_s) << "ms; ";
	}
	#endif
	return(stop_s);
}
#ifndef TIME_EXEC
	#pragma GCC diagnostic pop
#endif

void ThrowWarning(const char * message){
	cerr << message << endl;
}

double rsvs3d::SignedLogScale(double in){

	double out = sign(in);
	double logeps = -15; // aprox log10(__DBL_EPSILON__)
	if(out == 0){
		// Do nothing
	} else if (fabs(in)>__DBL_EPSILON__){
		out = out * (log10(fabs(in))-logeps+1.0);
	} else {
		out = -out * 1.0/(log10(fabs(in))-logeps+1.0);
	}

	return out;
}