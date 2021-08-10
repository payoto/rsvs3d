/**
 * Simple class containing all the information needed for the 3D-RSVS execution.
 *
 *@file
 */

#ifndef RSVSCLASS_H_INCLUDED
#define RSVSCLASS_H_INCLUDED

//=================================
// forward declared dependencies
// 		class foo; //when you only need a pointer not the actual object
// 		and to avoid circular dependencies

//=================================
// included dependencies
#include <fstream>
#include <iostream>

#include "RSVScalc.hpp"
#include "mesh.hpp"
#include "parameters.hpp"
#include "polyscopersvs.hpp"
#include "postprocessing.hpp"
#include "snake.hpp"
#include "triangulate.hpp"

//==================================
// Code
// This file defines a god class for conveniently doing the integration
// of the various compoments into a cohesive framework
//
namespace integrate
{
class RSVStecplot
{
  public:
    tecplotfile outSnake;
    tecplotfile outgradientsnake;
    tecplotfile outvectorsnake;
};
class RSVSclass
{
  public:
    param::parameters paramconf;
    tecplotfile outSnake;
    tecplotfile outgradientsnake;
    tecplotfile outvectorsnake;
    polyscopersvs::PolyScopeRSVS viewer;
    snake rsvsSnake;
    mesh snakeMesh;
    mesh stabilityMesh;
    mesh voluMesh;
    triangulation rsvsTri;
    RSVScalc calcObj;
    std::ofstream logFile;
    std::ofstream coutFile;
    std::ofstream cerrFile;
};
} // namespace integrate

#endif
