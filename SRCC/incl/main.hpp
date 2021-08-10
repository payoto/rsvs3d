/**
 * file containing the main functions and the command line parser.
 * @file
 */

#ifndef MAIN_H_INCLUDED
#define MAIN_H_INCLUDED

//=================================
// forward declared dependencies
// 		class foo; //when you only need a pointer not the actual object
// 		and to avoid circular dependencies

namespace param
{
class parameters;
}
namespace parse
{
class ParserOutput;
}
//=================================
// included dependencies

#include <string>
#include <vector>
//=================================
// Code
//
// This file defines parameters to be used in other parts of the
// RSVS snaking framework. Default values are defined in "parameters.cpp"
//
// Substructure names are all 4-5 letters

int RSVSExecution(int argc, char *argv[]);
void NoExecution(parse::ParserOutput &parseOut, param::parameters &paramconf);
void ExecuteTests(const parse::ParserOutput &parseOut);
#ifdef LIB_RSVS
int main_rsvs3d(int argc, char *argv[]);
#endif
namespace parse
{

class ParserOutput
{
  public:
    int execFlow = 0;
    std::string paramFileOut = "noexec_config.json";
    std::string testCase = "";
};

ParserOutput CommandLineParser(int argc, char *argv[], param::parameters &paramconf);
ParserOutput StringParser(std::vector<std::string> &commands, param::parameters &paramconf);
namespace config
{
void useconfig(const std::string &confCase, param::parameters &paramconf);
void loadconfig(const std::string &confCase, param::parameters &paramconf);

} // namespace config
} // namespace parse
#endif
