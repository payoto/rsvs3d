#include <cstdlib>
#include <iostream>
#include <cassert>
#include <fstream>
#include <string>

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

//================================
// io
//================================

void param::io::read(const std::string &fileName, parameters &p){
	std::ifstream file;
	json j, j_fin;

	file.open(fileName);
	j_fin = p; 
	file >> j;
	file.close();

	// Insert values read into the parameter structure
	j_fin.update(j);
	p=j_fin;
}

void param::io::write(const std::string &fileName, const parameters &p){
	std::ofstream file;
	json j;

	file.open(fileName);

	j = p; 
	file << j.dump(2);
	file.close();
}

void param::io::readflat(const std::string &fileName, parameters &p){
	std::ifstream file;
	json j, j_fin;

	file.open(fileName);
	j_fin = p; 
	file >> j;
	file.close();
	j = j.unflatten();
	// Insert values read into the parameter structure
	j_fin.update(j);
	p=j_fin;
}
void param::io::writeflat(const std::string &fileName, const parameters &p){
	std::ofstream file;
	json j;

	file.open(fileName);

	j = p; 
	file << j.flatten().dump(2);
	file.close();
}
void param::io::defaultconf(){
	param::parameters params;
	std::string fileName="config\\defaultconf.json";
	std::string fileNameFlat="config\\defaultconfflat.json";

	param::io::write(fileName, params);
	param::io::writeflat(fileNameFlat, params);

}
//================================
// Tests
//================================
int param::test::base(){
	/*
	
	*/
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
	param::io::defaultconf();
	return(0);
}

int param::test::io(){
	param::parameters params, params2;
	std::string fileName="..\\TESTOUT\\testioparam.json";
	json j1, j2;

	params.snak.arrivaltolerance = 1;
	param::io::write(fileName, params);

	param::io::read(fileName, params2);

	j1 = params; 
	j2 = params2;

	if (j1!=j2){
		std::cerr << "Error: Parameter read/write "
			<<" is not symmetrical" << std::endl;
		std::cerr << "In: " << __PRETTY_FUNCTION__ << std::endl;
		return (1);
	};

	return(0);
}

int param::test::ioflat(){
	param::parameters params, params2;
	std::string fileName="..\\TESTOUT\\testioflatparam.json";
	json j1, j2;

	params.snak.arrivaltolerance = 1;
	param::io::writeflat(fileName, params);

	param::io::readflat(fileName, params2);

	j1 = params; 
	j2 = params2;

	if (j1!=j2){
		std::cerr << "Error: Parameter read/write "
			<<" is not symmetrical" << std::endl;
		std::cerr << "In: " << __PRETTY_FUNCTION__ << std::endl;
		return (1);
	};

	return(0);
}

int param::test::ipartialread(){
	param::parameters params, params2;
	std::string fileName="config\\partialconfflat.json";
	std::string fileName2="config\\partialconfflat_out.json";
	json j1, j2;

	param::io::readflat(fileName, params2);
	param::io::writeflat(fileName2, params2);

	j1 = params; 
	j2 = params2;

	if (j1==j2){
		std::cerr << "Error: Parameter read/write "
			<<" is not symmetrical" << std::endl;
		std::cerr << "In: " << __PRETTY_FUNCTION__ << std::endl;
		return (1);
	} else {
		std::cout << "Partial read succesful, outputs are different" << std::endl;
	}

	return(0);
}
