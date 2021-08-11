#include "main.hpp"

#include <iostream>
#include <string>

#include "RSVSclass.hpp"
#include "RSVSintegration.hpp"
#include "cxxopts.hpp"
#include "makeontargetchange.h"
#include "parameters.hpp"

#ifdef RSVSTEST
#include "test.hpp"
#endif // RSVSTEST
#ifndef RSVS_NOTESTS
#include "test.hpp"
#endif

using namespace std;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"

#ifndef LIB_RSVS
int main(int argc, char *argv[])
#else
int main_rsvs3d(int argc, char *argv[])
#endif
{
// Execution paths for different compilation flags allows building test
// executables
#ifndef RSVSTEST
    return RSVSExecution(argc, argv);
#elif defined(TEST_ALL)
    return rsvstest::maintest();
#elif defined(TEST_SHORT)
    return rsvstest::shorttest();
#else
    return rsvstest::newtest();
#endif // RSVSTEST
}
#pragma GCC diagnostic pop

int RSVSExecution(int argc, char *argv[])
{
    param::parameters paramconf;

    auto parseOut = parse::CommandLineParser(argc, argv, paramconf);
    if (parseOut.execFlow > 0)
    {
        integrate::RSVSclass RSVSobj(parseOut.isHeadless);
        RSVSobj.paramconf = paramconf;
        if (parseOut.execFlow == 1)
        {
            integrate::execute::All(RSVSobj);
        }
        else if (parseOut.execFlow == 2)
        {
            integrate::execute::Interactive(RSVSobj);
        }
        else
        {
            RSVS3D_ERROR_ARGUMENT("Invalid execution flow");
        }
    }
    else if (parseOut.execFlow == -3)
    {
#ifndef RSVS_NOTESTS
        ExecuteTests(parseOut);
#else
        RSVS3D_ERROR("Tests not compiled; use `make all` to produce an executable"
                     "with tests;\n or `make notest` for an executable without");
#endif
    }
    else
    {
        // Output parameter file to the directory (NO OVERWRITE?)
        NoExecution(parseOut, paramconf);
    }
    return (0);
}

void ExecuteTests(const parse::ParserOutput &parseOut)
{
    int testNum = 0;

#ifndef RSVS_HIDETESTS
    if (parseOut.testCase.compare("all") == 0)
    {
        testNum = rsvstest::maintest();
    }
    else if (parseOut.testCase.compare("short") == 0)
    {
        testNum = rsvstest::shorttest();
    }
    else if (parseOut.testCase.compare("new") == 0)
    {
        testNum = rsvstest::newtest();
    }
    else
    {
        std::stringstream errStr;
        errStr << "The test case '" << parseOut.testCase << "'is not known, valid options are:" << std::endl
               << "'all', 'short', or 'new'.";
        RSVS3D_ERROR_ARGUMENT(errStr.str().c_str());
    }
#endif // RSVS_HIDETESTS

    exit(testNum);
}

void NoExecution(parse::ParserOutput &parseOut, param::parameters &paramconf)
{
    if (parseOut.execFlow == -2)
    {
        parseOut.paramFileOut = "failexec_rsvsconfig.json";
    }

    param::io::writeflat(parseOut.paramFileOut + "flat", paramconf);
    param::io::write(parseOut.paramFileOut, paramconf);

    if (parseOut.execFlow == -2)
    {
        std::cerr << "Error while parsing the arguments. Check generated '" << parseOut.paramFileOut << "' file";
        exit(-2);
    }
}

/**
Takes in the command line arguments and returns an integer specifying the
execution flow of the rest of the program.

@param[in]  argc       The number of arguments
@param      argv       the arguments
@param      paramconf  The parameter struncture

@return     Number indicating the state of the program after parsing: 0   does
            nothing; 1   Run RSVS; -1  No exec stops it; -2  Error parsing
            inputs.
*/
parse::ParserOutput parse::CommandLineParser(int argc, char *argv[], param::parameters &paramconf)
{
    parse::ParserOutput parseOut;
    parseOut.execFlow = 0;
    std::vector<std::string> triggerExec, configParam, configFiles, configPredef, noexecStr, testString;
    // options that will trigger execution of RSVS
    triggerExec = {"use-config", "load-config", "param"};

    std::string strDescription;
    strDescription = "\nProgram for the execution of the Restricted-Surface";
    strDescription += " Volume of Solid in 3D";
    cxxopts::Options options(argv[0], strDescription);
    options.positional_help("[optional args]").show_positional_help();
    // clang-format off
    options.add_options("")("h,help", "Print help");

    options.add_options("Parameter configuration")
        ("u,use-config", "Use one of the predefined configurations stored in the code.", cxxopts::value(configPredef),"STRING")
        ("l,load-config",std::string("Load configuration file in JSON format to set parameter "
                              "structure. Multiple files can be specified and will be processed"
                              " in order of appearance."),cxxopts::value(configFiles), "FILES")
        ("p,param", "Define a parameter manually on the command line."
                    " The format must be a flat key into the JSON configuration:"
                    " (e.g.: '/snak/maxsteps:50' will set 'param.snak.maxsteps=50')",
                    cxxopts::value(configParam),
                    "KEY:VAL")
        ("default-config", "Output the default configuration to a file.",
                            cxxopts::value<std::string>()->implicit_value("default_config"), "FILE");

    options.add_options("Execution control")
        ("n,noexec", std::string("Do not execute RSVS process, will only parse the inputs and output "
                                 "the resulting configuration file to 'arg'"),
                    cxxopts::value(noexecStr)->implicit_value("./noexec_config.json"), "STRING")
        ("e,exec", "Execute RSVS. With no command line argument the program does nothing.")
        ("i,interactive", "Execute the RSVS in interactive mode using polyscope")
        ("no-gui", "Disable all GUI calls and OpenGL functionality.")
#ifndef RSVS_HIDETESTS
        ("test", std::string("Executes specified tests. requires compilation without flag RSVS_NOTESTS."),
            cxxopts::value(testString)->implicit_value("short"), "STRING")
#endif
        ;
    // clang-format on
    auto result = options.parse(argc, argv);
    // ""
    if (result.count("help"))
    {
        std::cout << options.help({"", "Execution control", "Parameter configuration"}) << std::endl;
        exit(0);
    }
    // "Execution control"
    if (result.count("noexec") > 0 && result.count("exec") > 0)
    {
        std::cerr << std::endl << " Invalid specification of -e (exec) and -n (noexec)";
        std::cerr << " on the command line." << std::endl;
        std::cerr << " see --help for more info" << std::endl;
        exit(-1);
    }
    else if ((result.count("exec") || result.count("interactive")) && result.count("test"))
    {
        std::cerr << std::endl << " Invalid specification of -e/i (exec/interactive) and -t (test)";
        std::cerr << " on the command line." << std::endl;
        std::cerr << " see --help for more info" << std::endl;
        exit(-1);
    }
    else if (result.count("noexec") && result.count("test"))
    {
        std::cerr << std::endl << " Invalid specification of -n (noexec) and -t (test)";
        std::cerr << " on the command line." << std::endl;
        std::cerr << " see --help for more info" << std::endl;
        exit(-1);
    }
    else if (result.count("interactive"))
    {
        parseOut.execFlow = 2;
    }
    else if (result.count("exec"))
    {
        parseOut.execFlow = 1;
    }
    else if (result.count("test"))
    {
        parseOut.execFlow = -3;
        for (auto testCase : testString)
        {
            parseOut.testCase = testCase;
            std::cout << parseOut.testCase << " " << testCase << std::endl;
        }
    }
    else if (result.count("noexec"))
    {
        parseOut.execFlow = -1;
        for (auto confCase : noexecStr)
        {
            parseOut.paramFileOut = confCase;
            break;
        }
        std::cout << parseOut.paramFileOut << std::endl;
    }
    if (result.count("no-gui"))
    {
        parseOut.isHeadless = true;
    }
    // Parameter configuration
    // if one of the triggers is found specify that execution
    // should take place
    if (parseOut.execFlow == 0)
    {
        for (auto i : triggerExec)
        {
            if (result.count(i))
            {
                parseOut.execFlow = 1;
                break;
            }
        }
    }
    // Parse a predefined config:
    if (result.count("use-config"))
    {
        for (auto confCase : configPredef)
        {
            parse::config::useconfig(confCase, paramconf);
        }
    }
    if (result.count("load-config"))
    {
        for (auto confCase : configFiles)
        {
            parse::config::loadconfig(confCase, paramconf);
        }
    }
    if (result.count("param"))
    {
        if (param::io::updatefromstring(configParam, paramconf) > 0)
        {
            parseOut.execFlow = -2;
        }
    }

    return (parseOut);
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
void parse::config::useconfig(const std::string &confCase, param::parameters &paramconf)
{
    std::cerr << "No predefined cases yet." << std::endl;
}
#pragma GCC diagnostic pop
void parse::config::loadconfig(const std::string &confCase, param::parameters &paramconf)
{
    param::io::read(confCase, paramconf);
}

parse::ParserOutput parse::StringParser(std::vector<std::string> &commands, param::parameters &paramconf)
{
    int argc = commands.size();
    char **argv;
    argv = new char *[argc];
    for (int i = 0; i < argc; ++i)
    {
        argv[i] = new char[commands[i].length() + 1];
        strcpy(argv[i], commands[i].c_str());
    }

    auto parseOut = parse::CommandLineParser(argc, argv, paramconf);

    for (int i = 0; i < argc; ++i)
    {
        delete[] argv[i];
    }
    return (parseOut);
}
