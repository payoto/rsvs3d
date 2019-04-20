
#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
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

double rsvs3d::SignedLogScale(double in){

	double out = sign(in);

	if(out == 0){
		// Do nothing
	} else if (fabs(in)>__DBL_EPSILON__){
		out = out * (log10(fabs(in))-log10(__DBL_EPSILON__)+1.0);
	} else {
		out = -out * 1.0/(log10(fabs(in))-log10(__DBL_EPSILON__)+1.0);
	}

	return out;
}