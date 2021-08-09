/**
 * @file polyscopersvs.hpp
 * @author Alexandre Payot
 * @brief Provide an interface between the RSVS process and polyscope.
 * @version 0.1
 * @date 2021-08-08
 *
 * @copyright Copyright (c) 2021
 *
 *
 * This header includes classes and functions to manage the displaying of RSVS
 * elements through polyscope.
 *
 *
 */

#ifndef POLYSCOPERSVS_H_INCLUDED
#define POLYSCOPERSVS_H_INCLUDED

namespace integrate
{
class RSVSclass;
}

#include "mesh.hpp"
/**
 * @brief Namespace containing interfaces to polyscope for RSVS objects
 *
 */
namespace polyscopersvs
{
/// Flag controlling wether polyscope has been initialised or not
extern bool POLYSCOPE_INITIALISED;

/**
 * @brief A structure containing the information about the polyscope display
 * and the RSVS elements to display.
 *
 */
class PolyScopeRSVS
{
  public:
    PolyScopeRSVS();
    void show(size_t forFrames = 18446744073709551615);
    void addMesh(std::string name, const mesh &meshIn);
    void setInteractiveCallback(integrate::RSVSclass &RSVSobj);
};

/// Test functions for polyscope
namespace test
{
int init();
int show();
int meshShow();
} // namespace test

} // namespace polyscopersvs
#endif
