#include <cstdlib>

#include "parameters.hpp"


using namespace param;


//===========================================
// Voxel class method definitions
//===========================================

voxel::voxel(){
	for(int i=0; i<3; ++i){
		this->domain[i].lb=0.0;
		this->domain[i].ub=1.0;
	}

	this->gridsizebackground = {1, 1, 1};
	this->gridsizesnake = {6, 6, 6};
}

voxel::~voxel(){
}


//===========================================
// snaking class method definitions
//===========================================

snaking::snaking(){
	this->arrivaltolerance = 1e-7;
	this->multiarrivaltolerance = 1e-2; 
	this->initboundary = 1;
}

snaking::~snaking(){
}

//===========================================
// RSVS class method definitions
//===========================================

rsvs::rsvs(){
	this->solveralgorithm = 0;
}

rsvs::~rsvs(){
}