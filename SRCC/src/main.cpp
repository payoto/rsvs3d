#include <iostream>

#include "main.hpp"
#include "cxxopts.hpp"

int main(int argc, char* argv[])
{

	CommandLineParser(argc, argv);
}


void CommandLineParser(int argc, char* argv[])
{


	cxxopts::Options options(argv[0], " - example command line options");
    options
      .positional_help("[optional args]")
      .show_positional_help();
	
}