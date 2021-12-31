/**
 * Interface for the Snaking algorithm
 *
 * This file provides the base class which defines the interface for
 * a velocity calculation to be compatible with the RSVSclass and the RSVS
 * integration process.
 *
 *@file
 */

#ifndef SNAKEVELOCITYCALCULATOR_H_INCLUDED
#define SNAKEVELOCITYCALCULATOR_H_INCLUDED

namespace param
{
namespace dev
{
class devparam;
}
} // namespace param
class triangulation;

#include <fstream>

namespace rsvs3d
{
/**
 * @brief A base class which needs to be inherited from to implement a new
 * velocity calculation.
 */
class SnakeVelocityCalculator
{
  public:
    /**
     * @brief Set the Development parameters of the Calculator object.
     */
    virtual void setDevParameters(const param::dev::devparam &) = 0;

    /**
     * @brief Calculate volumes and derivatives based on a triangulation
     */
    virtual void CalculateTriangulation(const triangulation &, int = 0) = 0;

    /**
     * @brief Calculate the step required as specified by the algorithm
     */
    virtual void CheckAndCompute(int, bool = false) = 0;

    /**
     * @brief Return constraints (areas, volumes) to a mesh.
     *
     * @param triRSVS  The triangulation object associated with the snake object.
     */
    virtual void ReturnConstrToMesh(triangulation &triRSVS) const = 0;

    /**
     * @brief Returns velocities to the snake in the triangulation object.
     *
     * @param triRSVS  The triangulation object associated with the snake object.
     */
    virtual void ReturnVelocities(triangulation &triRSVS) = 0;

    /**
     * @brief Print convergence information to a log stream.
     *
     * @param file A file stream already opened
     * @param logLevel The level of logging, lower is more logs.
     */
    virtual void ConvergenceLog(std::ofstream &, int = 3) const = 0;
};
} // namespace rsvs3d
#endif // SNAKEVELOCITYCALCULATOR_H_INCLUDED
