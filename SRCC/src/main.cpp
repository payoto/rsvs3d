#include <iostream>
#include <string>

#include "main.hpp"
#include "cxxopts.hpp"
#include "parameters.hpp"
#include "RSVSclass.hpp"

int main(int argc, char* argv[]){
	int execFlow=0;
	param::parameters paramconf;

	execFlow = parse::CommandLineParser(argc, argv, paramconf);
	if (execFlow>0){
		integrate::RSVSclass RSVSobj;
	} else {
		// Output parameter file to the directory (NO OVERWRITE?)

		if (execFlow==-2){
			cerr << "Error while parsing the arguments. Check generated"
				"  'conf.json' file";
				exit(-2);
		}
	}
}


int parse::CommandLineParser(int argc, char* argv[], param::parameters &paramconf){
	/*
	Takes in the command line arguments and returns an integer specifying the
	execution flow of the rest of the program.
	*/

	int execFlow=0;
	std::vector<std::string> triggerExec,configParam, configFiles, configPredef;
	// options that will trigger execution of RSVS
	triggerExec = {"use-config","load-config","param"};
	/*
	execFlow controls the following execution of the flow. It is usually only
	modified if it is 0 (by config parameters) but can be overriden by the exec or
	noexec commands.
	execFlow values and corresponding actions:
		0   does nothing
		1   Run RSVS
		-1  No exec stops it
		-2  Error parsing inputs 
	*/

	std::string strDescription;
	strDescription = "\nProgram for the execution of the Restricted-Surface";
	strDescription += " Volume of Solid in 3D";
	cxxopts::Options options(argv[0], strDescription);
	options
	.positional_help("[optional args]")
	.show_positional_help();
	
	options.add_options("")
	("h,help", "Print help")
	;

	options.add_options("Parameter configuration")
	("u,use-config","Use one of the predefined configurations stored in the code.",
		cxxopts::value(configPredef), "STRING")
	("l,load-config", 
		std::string("Load configuration file in JSON format to set parameter structure.")
		+ std::string(" Multiple files can be specified and will be processed")
		+ std::string(" in order of appearance."), cxxopts::value(configFiles), "FILES")
	("p,param", "Define a parameter manually on the command line."
		" The format must be a flat key into the JSON configuration:"
		" (e.g.: '/snak/maxsteps:50' will set 'param.snak.maxsteps=50')",
		 cxxopts::value(configParam),
		"KEY:VAL")
	("default-config","Output the default configuration to a file.",
		cxxopts::value<std::string>()->implicit_value("default_config"), "FILE")
	;

	options.add_options("Execution control")
	("n,noexec",std::string("Do not execute RSVS process, ")
		+ std::string("will only parse the inputs and output ")
		+ std::string("the resulting configuration file to 'arg'"),
		cxxopts::value<std::string>()->implicit_value("./noexec_config.json"))
	("e,exec", "Execute RSVS. With no command line argument the program does nothing.")
	;


	auto result = options.parse(argc, argv);
	// ""
	if (result.count("help")){
		std::cout << options.help({"", "Execution control",
			"Parameter configuration"}) << std::endl;
		exit(0);
		execFlow = -1;
	}
	// "Execution control"
	if(result.count("noexec")>0 && result.count("exec")>0){
		std::cerr << std::endl << " Invalid specification of -e (exec) and -n (noexec)";
		std::cerr << " on the command line." << std::endl;
		std::cerr << " see --help for more info" << std::endl;
		exit(-1);
	}
	if (result.count("exec")){
		execFlow = 1;
	}
	if (result.count("noexec")){
		execFlow = -1;
	}

	// Parameter configuration
	// if one of the triggers is found specify that execution
	// should take place
	if (execFlow==0){
		for (auto i : triggerExec){
			if(result.count(i)){
				execFlow = 1;
				break;
			}
		}
	}
	// Parse a predefined config:
	if (result.count("use-config")){
		for (auto confCase : configPredef) {
			parse::config::useconfig(confCase, paramconf);
		}
	}
	if (result.count("load-config")){
		for (auto confCase : configFiles) {
			parse::config::loadconfig(confCase, paramconf);
		}
	}
	if (result.count("param")){
		if(param::io::updatefromstring(configParam, paramconf)>0){
			execFlow = -2;
		}
	}

	return(execFlow);
}

void parse::config::useconfig(const std::string &confCase,
	param::parameters &paramconf){

	std::cerr << "No predefined cases yet." << std::endl;

}
void parse::config::loadconfig(const std::string &confCase,
	param::parameters &paramconf){

	param::io::read(confCase, paramconf);

}
void parse::config::param(const std::string &confCase,
	param::parameters &paramconf){

	std::cerr << "No predefined cases yet." << std::endl;

}