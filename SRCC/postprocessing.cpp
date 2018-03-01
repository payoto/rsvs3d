#include <iostream>
#include "arraystructures.hpp"
#include "postprocessing.hpp"


// Functions


// Class function Implementation

int tecplotfile::OpenFile (char *str){

	fid=fopen(str,"w");
	if (fid==NULL){
		return(-1);
	}
	return(0);
}

// Test Code