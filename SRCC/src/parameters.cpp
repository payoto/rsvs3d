#include <cstdlib>
#include <iostream>
#include <cassert>
#include <fstream>
#include <string>


#include <iomanip>
#include <ctime>
#include <sstream>

#include "parameters.hpp"
#include "rsvsjson.hpp"
#include "warning.hpp"


//===========================================
// Bounds template class method definitions
//===========================================
using namespace std;
using json = rsvsjson::json;

template<class T>
void param::to_json(json& j, const filltype<T>& p){
	j = json{
		{"active", p.active},
		{"fill", p.fill},
	};
}
template<class T>
void param::from_json(const json& j, filltype<T>& p){
	j.at("active").get_to(p.active);
	j.at("fill").get_to(p.fill);

}

//===========================================
// voxel class method definitions
//===========================================

param::voxel::voxel(){

	this->gridsizebackground = {1, 1, 1};
	this->gridsizesnake = {6, 6, 6};
}

param::voxel::~voxel(){
}
void param::voxel::PrepareForUse(){
}
void param::to_json(json& j, const voxel& p){
	j = json{
		{"gridsizebackground", p.gridsizebackground},
		{"gridsizesnake", p.gridsizesnake},
	};
}
void param::from_json(const json& j, voxel& p){
	j.at("gridsizebackground").get_to(p.gridsizebackground);
	j.at("gridsizesnake").get_to(p.gridsizesnake);
}

//===========================================
// Voronoi class method definitions
//===========================================

param::voronoi::voronoi(){

	this->inputpoints={0.0};
	this->pointfile = "";
	this->distancebox = 0.51;
	this->snakecoarseness = 0.0;
	this->vorosnakelayers = 2;
}

param::voronoi::~voronoi(){
}
void param::voronoi::PrepareForUse(){
}
void param::to_json(json& j, const voronoi& p){
	j = json{
		{"inputpoints", p.inputpoints},
		{"distancebox", p.distancebox},
		{"pointfile", p.pointfile},
		{"snakecoarseness", p.snakecoarseness},
		{"vorosnakelayers", p.vorosnakelayers},
	};
}
void param::from_json(const json& j, voronoi& p){
	p.inputpoints=j.at("inputpoints").get<std::vector<double>>();
	// j.at("inputpoints").get_to(p.inputpoints);
	j.at("distancebox").get_to(p.distancebox);
	j.at("pointfile").get_to(p.pointfile);
	j.at("snakecoarseness").get_to(p.snakecoarseness);
	j.at("vorosnakelayers").get_to(p.vorosnakelayers);
}
void param::voronoi::ReadPoints(){
	std::ifstream pointstream;

	pointstream.open(this->pointfile, std::ios::in);
	CheckFStream(pointstream, __PRETTY_FUNCTION__, this->pointfile);

	this->inputpoints.clear();

	double temp;
	while(pointstream >> temp){
		this->inputpoints.push_back(temp);
		std::cout << temp << " ";
	}
	pointstream.close();
	std::cout << std::endl;

}

//===========================================
// snaking class method definitions
//===========================================

param::snaking::snaking(){
	this->arrivaltolerance = 1e-7;
	this->multiarrivaltolerance = 1e-2; 
	this->snaxtimestep=0.9;
	this->snaxdiststep=0.9;
	this->initboundary = 1;
	this->maxsteps = 50;
	this->spawnposition = 1e-5;
}

param::snaking::~snaking(){
}
void param::to_json(json& j, const snaking& p){
	j = json{
		{"arrivaltolerance", p.arrivaltolerance},
		{"multiarrivaltolerance", p.multiarrivaltolerance},
		{"snaxtimestep", p.snaxtimestep},
		{"snaxdiststep", p.snaxdiststep},
		{"initboundary", p.initboundary},
		{"maxsteps", p.maxsteps},
		{"spawnposition", p.spawnposition},
	};
}
void param::from_json(const json& j, snaking& p){
	j.at("arrivaltolerance").get_to(p.arrivaltolerance);
	j.at("multiarrivaltolerance").get_to(p.multiarrivaltolerance);
	j.at("snaxtimestep").get_to(p.snaxtimestep);
	j.at("snaxdiststep").get_to(p.snaxdiststep);
	j.at("initboundary").get_to(p.initboundary);
	j.at("maxsteps").get_to(p.maxsteps);
	j.at("spawnposition").get_to(p.spawnposition);
}
void param::snaking::PrepareForUse(){
}

//===========================================
// RSVS class method definitions
//===========================================

param::rsvs::rsvs(){
	this->solveralgorithm = 0;

	this->cstfill.active=false;
	this->cstfill.fill=0.5;

	this->filefill.active=false;
	this->filefill.fill="";
	
	this->makefill.active=false;
	this->makefill.fill="";
}

param::rsvs::~rsvs(){
}
void param::to_json(json& j, const rsvs& p){
	j = json{
		{"solveralgorithm", p.solveralgorithm},
		{"cstfill", p.cstfill},
		{"filefill", p.filefill},
		{"makefill", p.makefill},
	};
}
void param::from_json(const json& j, rsvs& p){
	j.at("solveralgorithm").get_to(p.solveralgorithm);
	j.at("cstfill").get_to(p.cstfill);
	j.at("filefill").get_to(p.filefill);
	j.at("makefill").get_to(p.makefill);
}
void param::rsvs::PrepareForUse(){
}

//===========================================
// grid class method definitions
//===========================================

param::grid::grid(){
	for(int i=0; i<3; ++i){
		this->domain[i][0]=0.0;
		this->domain[i][1]=1.0;
		this->physdomain[i][0]=0.0;
		this->physdomain[i][1]=1.0;
	}

	this->activegrid = "voxel";
}
void param::to_json(json& j, const grid& p){
	j = json{
		{"domain", p.domain},
		{"voxel", p.voxel},
		{"voronoi", p.voronoi},
		{"physdomain", p.physdomain},
		{"activegrid", p.activegrid},
	};
}
void param::from_json(const json& j, grid& p){
	j.at("domain").get_to(p.domain);
	j.at("voxel").get_to(p.voxel);
	j.at("voronoi").get_to(p.voronoi);
	j.at("physdomain").get_to(p.physdomain);
	j.at("activegrid").get_to(p.activegrid);
}
void param::grid::PrepareForUse(){
	if(this->activegrid.compare("voronoi")==0) {
		// Load data into this->voronoi
		if(this->voronoi.inputpoints.size()<=1){
			this->voronoi.ReadPoints();
		}
	}
	this->voxel.PrepareForUse();
	this->voronoi.PrepareForUse();
}

//===========================================
// File and io classes method definitions
//===========================================


param::ioin::ioin(){
	this->snakemeshname = "";
	this->volumeshname = "";
	this->snakefile = "";
	this->casename = "";
}

void param::to_json(json& j, const ioin& p){
	j = json{
		{"snakemeshname", p.snakemeshname},
		{"volumeshname", p.volumeshname},
		{"snakefile", p.snakefile},
		{"casename", p.casename},
	};
}
void param::from_json(const json& j, ioin& p){
	j.at("snakemeshname").get_to(p.snakemeshname);
	j.at("volumeshname").get_to(p.volumeshname);
	j.at("snakefile").get_to(p.snakefile);
	j.at("casename").get_to(p.casename);
}
void param::ioin::PrepareForUse(){
	if(!this->snakefile.empty() 
		&& (this->volumeshname.empty() || this->snakemeshname.empty())){
		RSVS3D_ERROR_ARGUMENT("Incompatible ioin parameters: mesh files must "
			"be defined to load a 'snakefile'.");
	}
}

param::ioout::ioout(){
	this->pathoutdir = "../out";
	this->pathpattern = "Archive_%Y_%m/Day_%y-%m-%d";
	this->basenamepattern = "%y%m%dT%H%M%S_";
	this->basenameoutdir="rsvs3d_";
	this->outdir="";
	this->pattern="";

	this->redirectcout=false;
	this->redirectcerr=false;

	this->logginglvl=2;
	this->outputlvl=2;
}

void param::ioout::PrepareForUse(){

	std::time_t t_ptr;
	auto t = std::time(&t_ptr);
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
		this->pattern = oss.str();
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
		{"pattern", p.pattern},
		{"redirectcout", p.redirectcout},
		{"redirectcerr", p.redirectcerr},
		{"logginglvl", p.logginglvl},
		{"outputlvl", p.outputlvl},
	};
}
void param::from_json(const json& j, ioout& p){
	j.at("pathoutdir").get_to(p.pathoutdir);
	j.at("pathpattern").get_to(p.pathpattern);
	j.at("basenamepattern").get_to(p.basenamepattern);
	j.at("basenameoutdir").get_to(p.basenameoutdir);
	j.at("outdir").get_to(p.outdir);
	j.at("pattern").get_to(p.pattern);
	j.at("redirectcout").get_to(p.redirectcout);
	j.at("redirectcerr").get_to(p.redirectcerr);
	j.at("logginglvl").get_to(p.logginglvl);
	j.at("outputlvl").get_to(p.outputlvl);
}

param::files::files(){
	this->appcasename2outdir=true;
	this->exportconfig = {{"", ""}};
}
void param::files::PrepareForUse(){
	this->ioout.PrepareForUse();
	this->ioin.PrepareForUse();

	if(this->appcasename2outdir){
		this->ioout.outdir += this->ioin.casename;
		this->ioout.pattern += this->ioin.casename;
		this->appcasename2outdir=false;
	}
	std::cout << "Output folder: " << this->ioout.outdir << std::endl;
}
void param::to_json(json& j, const files& p){
	j = json{
		{"ioin", p.ioin},
		{"ioout", p.ioout},
		{"appcasename2outdir", p.appcasename2outdir},
		{"exportconfig", p.exportconfig},
	};
}
void param::from_json(const json& j, files& p){
	j.at("ioin").get_to(p.ioin);
	j.at("ioout").get_to(p.ioout);
	j.at("appcasename2outdir").get_to(p.appcasename2outdir);
	p.exportconfig=j.at("exportconfig").get<param::exports>();
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


void param::io::read(const std::string &fileName, parameters &p){
	std::ifstream file;
	json j, j_fin;

	file.open(fileName);
	CheckFStream(file, __PRETTY_FUNCTION__, fileName);

	j_fin = p; 
	file >> j;
	file.close();

	try{
		j = j.unflatten();
	} catch (std::exception const& ex) {
		// if j is already not flat catch the exception and move on
		// TODO check the correct exception is being thrown (one day)
	}

	// Insert values read into the parameter structure
	rsvsjson::flatupdate(j_fin,j,false, false);
	
	param::from_json(j_fin,p);
}

void param::io::write(const std::string &fileName, const parameters &p){
	std::ofstream file;
	json j;

	file.open(fileName);

	CheckFStream(file, __PRETTY_FUNCTION__, fileName);
	j = p; 
	file << j.dump(2);
	file.flush();
	file.close();
}

void param::io::readflat(const std::string &fileName, parameters &p){
	std::ifstream file;
	json jnew, jfin;

	file.open(fileName);
	CheckFStream(file, __PRETTY_FUNCTION__, fileName);

	jfin = p; 
	// std::cout << "jfin assigned " << std::endl;
	file >> jnew;
	file.close();
	rsvsjson::flatupdate(jfin,jnew,false, true);
	// p=jfin;
	param::from_json(jfin,p);

	// std::cout << "p assigned " << std::endl;
}
void param::io::writeflat(const std::string &fileName, const parameters &p){
	std::ofstream file;
	json j;

	file.open(fileName);
	CheckFStream(file, __PRETTY_FUNCTION__, fileName);

	j = p; 
	file << j.flatten().dump(2);
	file.flush();
	file.close();
}
void param::io::defaultconf(){
	param::parameters params;
	std::string fileName("config/defaultconf.json");
	std::string fileNameFlat("config/defaultconfflat.json");

	param::io::write(fileName, params);
	param::io::writeflat(fileNameFlat, params);

	json j;
	j = params;
	std::cout << j.flatten().dump(2);
}

int param::io::updatefromstring(const std::vector<std::string> &flatjsonKeyVal,
	parameters &p, const std::string&& sep){
	/*Parses json key value pairs returning the number of failures.*/
	int numFail=0;
	json j = p;
	std::size_t pos;
	std::string key, val;
	j=j.flatten();

	for (auto keyVal : flatjsonKeyVal){
		pos = keyVal.find(sep);
		if(pos==std::string::npos){ // if separator was found
			std::cerr << "Warning in: " << std::endl << "  "
				 << __PRETTY_FUNCTION__ << std::endl;
			std::cerr << "  Separator " << sep 
				<< " Was not found in JSON key/val: " 
				<< keyVal << std::endl << std::endl;
			++numFail;
			continue;
		}
		// split the string and parse it
		key = keyVal.substr(0, pos);
		val = keyVal.substr(pos+1);
		if(j[key]==nullptr){
			std::cerr << "Warning in: " << std::endl << "  "
				 << __PRETTY_FUNCTION__ << std::endl;
			std::cerr << "  key " << key << " Was not found in JSON object " 
				<< std::endl << "  The following keyval will not be used: "
				<<  keyVal << std::endl << std::endl;
			++numFail;
			continue;
		}

		if (!j[key].is_number() && !j[key].is_string() && !j[key].is_boolean()){
			std::cerr << "Warning in: " << std::endl << "  "
				 << __PRETTY_FUNCTION__ << std::endl;
			std::cerr << "  Value of j[" << key << "] is neither a string or number, " 
				<< " parsing not supported." << std::endl << std::endl;
			++numFail;
			continue;
		}

		if(j[key].is_string()){
			j[key]=val;
		} else if (j[key].is_number_integer()){
			j[key]=std::stoi(val);
		} else if (j[key].is_number_float()){
			j[key]=std::stod(val);
		} else if (j[key].is_boolean()){
			std::istringstream is(val);
			is >> std::boolalpha >> j[key];
		}
	}
	j=j.unflatten();
	param::from_json(j,p);
	return(numFail);
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

	param::from_json(j,params2);

	j2 = params2;
	std::cout << "j " << j.flatten().dump(2);
	std::cout << "j2 " << j2.flatten().dump(2);

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
	std::string fileName="../TESTOUT/testioparam.json";
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
	std::string fileName="../TESTOUT/testioflatparam.json";
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
	std::string fileName="config/partialconfflat.json";
	std::string fileName2="config/partialconfflat_out.json";
	json j1, j2;
	std::cout << "Start read" << std::endl;
	param::io::readflat(fileName, params2);
	std::cout << "succesful read" << std::endl;
	param::io::writeflat(fileName2, params2);
	std::cout << "succesful write" << std::endl;

	j1 = params; 
	j2 = params2;
	std::cout << j1.flatten().dump(2);
	std::cout << j2.flatten().dump(2);
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
	std::string fileName="config/confprepared.json";
	std::string fileName2="config/confprepared2.json";
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

int param::test::autoflat(){
	param::parameters params, params2;

	json j1, j2;
	
	j1 = params;
	j2 = params2;

	j1 = j1.flatten();
	std::cout << "Dump j1 flattened" << std::endl << j1.dump(1);
	j1 = j1.flatten();
	std::cout << "Dump j1 flattened twice" << std::endl << j1.dump(1);

	try{
		j2 = j2.unflatten();
		std::cout << "Dump j2 unflattened" << std::endl << j2.dump(1);
	} catch (std::exception const& ex) {
		std::cout << "Cannot be unflattened and throws exception" << std::endl;
		std::cout << ex.what() << std::endl;

		std::cout << "Dump j2 unflattened" << std::endl << j2.dump(1);
	}


	return(0);
}

int param::test::symmetry(){
	json j1, j2;
	param::parameters v1;
	// std::vector<double> v1, v2;

	// v2 = v1;
	v1.grid.voronoi.inputpoints.push_back(1.0);
	j1=v1;
	// v2 = j1
	j1.get_to(v1);
	j2=v1;
	
	std::cout << j1 << endl;
	std::cout << j2 << endl;

	if (j1!=j2){
		std::cout << j1 << endl;
		std::cout << j2 << endl;
		std::cerr << "Error: Parameter conversion to JSON "
			<<" is not symmetrical" << std::endl;
		std::cerr << "In: " << __PRETTY_FUNCTION__ << std::endl;
		return (1);
	};

	return(0);
}

int param::test::symmetry_makefillactive(){
	json j1, j2;
	param::parameters v1;
	// std::vector<double> v1, v2;
	if(v1.rsvs.makefill.active){
		RSVS3D_ERROR_NOTHROW("makefile has been turned active at "
			"construction.");
		return (3);
	}

	j1=v1;
	// v2 = j1
	j1.get_to(v1);
	if(v1.rsvs.makefill.active){
		RSVS3D_ERROR_NOTHROW("makefile has been turned active when returned "
			"from json.");
		return (2);
	}
	j2=v1;
	
	std::cout << j1 << endl;
	std::cout << j2 << endl;

	if (j1!=j2){
		std::cout << j1 << endl;
		std::cout << j2 << endl;
		std::cerr << "Error: Parameter conversion to JSON "
			<<" is not symmetrical" << std::endl;
		std::cerr << "In: " << __PRETTY_FUNCTION__ << std::endl;
		return (1);
	};

	return(0);
}
