#include <cstdlib>
#include <iostream>
#include <cassert>
#include <fstream>
#include <string>


#include <iomanip>
#include <ctime>
#include <sstream>

#include "json.hpp"
#include "parameters.hpp"
#include "parameters2json.hpp"


using nlohmann::json;


//===========================================
// Bounds template class method definitions
//===========================================

// template <class T>
// void param::to_json(json& j, const bounds<T>& p){
// 	j = json{
// 		{"lb", p.lb},
// 		{"ub", p.ub},
// 	};
// }
// template <class T>
// void param::from_json(const json& j, bounds<T>& p){
// 	j.at("lb").get_to(p.lb);
// 	j.at("ub").get_to(p.ub);
// }

//===========================================
// voxel class method definitions
//===========================================

param::voxel::voxel(){
	for(int i=0; i<3; ++i){
		this->domain[i][0]=0.0;
		this->domain[i][1]=1.0;
	}

	this->gridsizebackground = {1, 1, 1};
	this->gridsizesnake = {6, 6, 6};

}

param::voxel::~voxel(){
}
void param::voxel::PrepareForUse(){
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
void param::snaking::PrepareForUse(){
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
void param::rsvs::PrepareForUse(){
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
void param::grid::PrepareForUse(){
	this->voxel.PrepareForUse();

}

//===========================================
// File and io classes method definitions
//===========================================


param::ioin::ioin(){
	this->snakemeshname = "";
	this->volumeshname = "";
	this->targetfill = "";
	this->casename = "";
}

void param::to_json(json& j, const ioin& p){
	j = json{
		{"snakemeshname", p.snakemeshname},
		{"volumeshname", p.volumeshname},
		{"targetfill", p.targetfill},
		{"casename", p.casename},
	};
}
void param::from_json(const json& j, ioin& p){
	j.at("snakemeshname").get_to(p.snakemeshname);
	j.at("volumeshname").get_to(p.volumeshname);
	j.at("targetfill").get_to(p.targetfill);
	j.at("casename").get_to(p.casename);
}
void param::ioin::PrepareForUse(){

}

param::ioout::ioout(){
	this->pathoutdir = "../out";
	this->pathpattern = "Archive_%Y_%m/Day_%y-%m-%d";
	this->basenamepattern = "%y%m%dT%H%M%S_";
	this->basenameoutdir="rsvs3d_";
	this->outdir="";

}
void param::ioout::PrepareForUse(){

	auto t = std::time(nullptr);
	auto tm = *std::localtime(&t);
	std::ostringstream oss;
	if (this->pathoutdir.size()==0){
		this->pathoutdir += ".";
	}
	if (this->pathpattern.size()>0){
		oss << std::put_time(&tm, this->pathpattern.c_str());
		this->pathoutdir += "/";
		this->pathoutdir += oss.str();
		oss.str(std::string());
	}
	if (this->basenamepattern.size()>0){
		oss.str(std::string());
		oss << std::put_time(&tm, this->basenamepattern.c_str());
		this->basenameoutdir += oss.str();
	}
	if(this->outdir.size()==0){
		this->outdir = this->pathoutdir + "/" + this->basenameoutdir;
	}
}
void param::to_json(json& j, const ioout& p){
	j = json{
		{"pathoutdir", p.pathoutdir},
		{"pathpattern", p.pathpattern},
		{"basenamepattern", p.basenamepattern},
		{"basenameoutdir", p.basenameoutdir},
		{"outdir", p.outdir},
	};
}
void param::from_json(const json& j, ioout& p){
	j.at("pathoutdir").get_to(p.pathoutdir);
	j.at("pathpattern").get_to(p.pathpattern);
	j.at("basenamepattern").get_to(p.basenamepattern);
	j.at("basenameoutdir").get_to(p.basenameoutdir);
	j.at("outdir").get_to(p.outdir);
}

param::files::files(){
	this->appcasename2outdir=true;
}
void param::files::PrepareForUse(){
	this->ioout.PrepareForUse();
	this->ioin.PrepareForUse();

	if(this->appcasename2outdir){
		this->ioout.outdir += this->ioin.casename;
		this->appcasename2outdir=false;
	}
}
void param::to_json(json& j, const files& p){
	j = json{
		{"ioin", p.ioin},
		{"ioout", p.ioout},
		{"appcasename2outdir", p.appcasename2outdir},
	};
}
void param::from_json(const json& j, files& p){
	j.at("ioin").get_to(p.ioin);
	j.at("ioout").get_to(p.ioout);
	j.at("appcasename2outdir").get_to(p.appcasename2outdir);
}

//===========================================
// parameters class method definitions
//===========================================
void param::parameters::PrepareForUse(){

	this->rsvs.PrepareForUse();
	this->snak.PrepareForUse();
	this->grid.PrepareForUse();
	this->files.PrepareForUse();

}
void param::to_json(json& j, const parameters& p){
	j = json{
		{"rsvs", p.rsvs},
		{"snak", p.snak},
		{"grid", p.grid},
		{"files", p.files}
	};
}
void param::from_json(const json& j, parameters& p){
	j.at("rsvs").get_to(p.rsvs);
	j.at("snak").get_to(p.snak);
	j.at("grid").get_to(p.grid);
	j.at("files").get_to(p.files);
}

//================================
// io
//================================

void param::flatupdate(json& jfin, json& jnew,
	bool isFlatFin, bool isFlatNew){
	/*
	Allows recursing updae into sub fields
	*/
	// std::cout << "file read " << std::endl;
	if(!isFlatNew){
		jnew = jnew.flatten();
	}
	if(!isFlatFin){
		jfin = jfin.flatten();
	}
		// std::cout << "j unflattened " << std::endl;
		// std::cout << jnew.dump(1) << std::endl;
		// std::cout << jfin.dump(1) << std::endl;
	// Insert values read into the parameter structure
	jfin.update(jnew);
		// std::cout << "jfin updated " << std::endl;
		// std::cout << jfin.dump(1) << std::endl;
	jfin = jfin.unflatten();
	jnew = jnew.unflatten();
}


void param::io::read(const std::string &fileName, parameters &p){
	std::ifstream file;
	json j, j_fin;

	file.open(fileName);
	j_fin = p; 
	file >> j;
	file.close();

	// Insert values read into the parameter structure
	param::flatupdate(j_fin,j,false, false);
	
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
	json jnew, jfin;

	file.open(fileName);
	jfin = p; 
	// std::cout << "jfin assigned " << std::endl;
	file >> jnew;
	file.close();
	param::flatupdate(jfin,jnew,false, true);
	p=jfin;
	// std::cout << "p assigned " << std::endl;
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
	std::cout << "Start read" << std::endl;
	param::io::readflat(fileName, params2);
	std::cout << "succesful read" << std::endl;
	param::io::writeflat(fileName2, params2);
	std::cout << "succesful write" << std::endl;

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


int param::test::prepareforuse(){
	param::parameters params, params2;
	std::string fileName="config\\confprepared.json";
	std::string fileName2="config\\confprepared2.json";
	json j1, j2;
	std::cout << "Start read" << std::endl;
	param::io::readflat(fileName, params2);
	std::cout << "succesful read" << std::endl;
	param::io::writeflat(fileName2, params2);
	std::cout << "succesful write" << std::endl;

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
