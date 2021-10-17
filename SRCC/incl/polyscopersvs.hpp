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
 * elements through polyscope. The key class for this interaction is the
 * `PolyScopeRSVS`.
 *
 * To use the functionality of this module call:
 * ```
 * RSVS3D -i
 * ```
 */

#ifndef POLYSCOPERSVS_H_INCLUDED
#define POLYSCOPERSVS_H_INCLUDED

namespace integrate
{
class RSVSclass;
}

class mesh;
class snake;
#include <string>
#include <vector>

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
 * The normal lifetime of that class is to call its `show` method when you want to
 * display the polyscope window. To use the RSVS in interactive mode you can call the
 * `setInteractiveCallback` method which enables interaction with the RSVS's configuration
 * and execution flow.
 *
 */
class PolyScopeRSVS
{
  private:
    /// Captures whether the RSVS window is being displayed or not.
    bool isHeadless;

  public:
    PolyScopeRSVS();
    /**
     * @brief Construct a new Poly Scope RSVS object
     *
     * @param isHeadless Whether to use the headless mode or not, when headless the OpenGL mocked back end is used.
     */
    PolyScopeRSVS(bool isHeadless);
    /**
     * @brief Show the polyscope window.
     *
     * @param forFrames The number of frames for which to display the window.
     *                  This can be used to make show non blocking.
     */
    void show(size_t forFrames = SIZE_MAX);
    /**
     * @brief Add a mesh with a given name to the polyscope window
     *
     * This method plots the given mesh as a surface mesh in polyscope.
     *
     * @param name The name (in polyscope) to give to that mesh, if it is the same as
     *             an existing one the previous one will be overwritten.
     * @param meshIn The mesh object to plot.
     */
    void addMesh(std::string name, const mesh &meshIn);
    /**
     * @brief PLot a snake with it's velocity in the given polyscope window
     *
     * @param name The name of the mesh which will be used to display the snake.
     * @param snakeIn The snake object which will be plotted.
     */
    void addSnake(std::string name, const snake &snakeIn);
    /**
     * @brief Plot specified volume elements into polyscope
     *
     * This method can be used to highlight specific volume elements in the mesh.
     *
     * @param name The name to give to the polyscope object.
     * @param meshIn The mesh from which the elements need to be extracted.
     * @param cellIndices The indicies of the volume elements to be displayed.
     * @param isIndex Wether an index (true) or the position in the array (false) is used to identify volume
     * elements.
     * @return float The value of the target volume fraction in the last volume cell that was accessed.
     */
    float addCells(std::string name, const mesh &meshIn, const std::vector<int> &&cellIndices, bool isIndex = true);
    /**
     * @brief Set polyscope's user callback to a lambda which allows control of the RSVSobj.
     *
     * This method is called to run the RSVS process in interactive mode. Once you have initialised the RSVS process you
     * are ready to call this method. The callback it will set captures the RSVS object and allows an interactive
     * configuration, control and visualization.
     *
     * @param RSVSobj The main object of the RSVS process. A reference is kept in the callback to allow iteration.
     */
    void setInteractiveCallback(integrate::RSVSclass &RSVSobj);
};

/// Test functions for polyscope-RSVS integration
namespace test
{
int init();
int show();
int meshShow();
} // namespace test

} // namespace polyscopersvs
#endif
