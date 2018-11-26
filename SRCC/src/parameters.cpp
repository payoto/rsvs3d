#include <cstdlib>
#include <iostream>
#include <cassert>

#include "json.hpp"
#include "parameters.hpp"


using nlohmann::json;


//===========================================
// Bounds template class method definitions
//===========================================

template <class T>
void param::to_json(json& j, const bounds<T>& p){
	j = json{
		{"lb", p.lb},
		{"ub", p.ub},
	};

}
template <class T>
void param::from_json(const json& j, bounds<T>& p){
	j.at("lb").get_to(p.lb);
	j.at("ub").get_to(p.ub);
}

//===========================================
// Voxel class method definitions
//===========================================

param::voxel::voxel(){
	for(int i=0; i<3; ++i){
		this->domain[i].lb=0.0;
		this->domain[i].ub=1.0;
	}

	this->gridsizebackground = {1, 1, 1};
	this->gridsizesnake = {6, 6, 6};

}

param::voxel::~voxel(){
}
void param::to_json(json& j, const voxel& p){
	j = json{
		{"domain", p.domain},
		{"gridsizebackground", p.gridsizebackground},
		{"gridsizesnake", p.gridsizesnake},
	};
}
void param::from_json(const json& j, voxel& p){
	j.at("domain").get_to(p.domain);
	j.at("gridsizebackground").get_to(p.gridsizebackground);
	j.at("gridsizesnake").get_to(p.gridsizesnake);
}


//===========================================
// snaking class method definitions
//===========================================

param::snaking::snaking(){
	this->arrivaltolerance = 1e-7;
	this->multiarrivaltolerance = 1e-2; 
	this->initboundary = 1;
	this->maxsteps = 50;
}

param::snaking::~snaking(){
}
void param::to_json(json& j, const snaking& p){
	j = json{
		{"arrivaltolerance", p.arrivaltolerance},
		{"multiarrivaltolerance", p.multiarrivaltolerance},
		{"initboundary", p.initboundary},
		{"maxsteps", p.maxsteps},
	};
}
void param::from_json(const json& j, snaking& p){
	j.at("arrivaltolerance").get_to(p.arrivaltolerance);
	j.at("multiarrivaltolerance").get_to(p.multiarrivaltolerance);
	j.at("initboundary").get_to(p.initboundary);
	j.at("maxsteps").get_to(p.maxsteps);
}

//===========================================
// RSVS class method definitions
//===========================================

param::rsvs::rsvs(){
	this->solveralgorithm = 0;
}

param::rsvs::~rsvs(){
}
void param::to_json(json& j, const rsvs& p){
	j = json{
		{"solveralgorithm", p.solveralgorithm},
	};
}
void param::from_json(const json& j, rsvs& p){
	j.at("solveralgorithm").get_to(p.solveralgorithm);
}

//===========================================
// grid class method definitions
//===========================================

void param::to_json(json& j, const grid& p){
	j = json{
		{"voxel", p.voxel},
	};
}
void param::from_json(const json& j, grid& p){
	j.at("voxel").get_to(p.voxel);
}

//===========================================
// parameters class method definitions
//===========================================
void param::to_json(json& j, const parameters& p){
	j = json{
		{"rsvs", p.rsvs},
		{"snak", p.snak},
		{"grid", p.grid},
	};
}
void param::from_json(const json& j, parameters& p){
	j.at("rsvs").get_to(p.rsvs);
	j.at("snak").get_to(p.snak);
	j.at("grid").get_to(p.grid);
}


int param::test(){

	param::parameters params, params2;
	json j, j2; 

	j = params;

	std::cout << "json object" << std::endl;
	std::cout << j.dump(2) << std::endl;
	std::cout << "flattened json object" << std::endl;
	std::cout << j.flatten().dump(2) << std::endl;

	params2 = j;

	j2 = params2;

	if (j!=j2){
		std::cerr << "Error: Parameter conversion to JSON "
			<<" is not symmetrical" << std::endl;
		std::cerr << "In: " << __PRETTY_FUNCTION__ << std::endl;
		return (1);
	};
	return 0;
}