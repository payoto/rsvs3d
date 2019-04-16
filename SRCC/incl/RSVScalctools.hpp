/**
 * Provides internal tools for RSVS calc.
 * 
 *@file
 */

//===============================================
// Include Guards
#ifndef RSVSCALCTOOLS_HPP
#define RSVSCALCTOOLS_HPP

//===============================================
// Levels of debuging Guards
#ifdef DEBUGLVL2 // All Debugging calls
#define DEBUGLVL1
#define TEST_ALL
#endif

#ifdef DEBUGLVL1 // Debugging of new features.

#endif

// forward declaration

class triangulation;
class triangle;

// Includes

#include <vector>
#include <Eigen>

#include "arraystructures.hpp"
#include "vectorarray.hpp"



// Declarations 

/**
 * Get all the coordinates needed to define a triangle.
 *
 * @param[in]  triIn    The triangle object to be analysed.
 * @param[in]  triRSVS  The triangulation object containing the triangle.
 *
 * @return     A vector of constant coordinate vectors, these allow the data to
 *             be read but not edited.
 */
vector<vector<double> const *> TrianglePointerCoordinates(const triangle &triIn,
	const triangulation& triRSVS);

/**
 * Get the active design variables defined by a triangle object.
 *
 * @param[in]  triIn     The triangle object.
 * @param[in]  triRSVS   The triangulation containing the triangle.
 * @param[in]  objDvMap  The calculation (RSVScalc) object's design variable
 *                       map.
 *
 * @return     A hashed vector of the active design variables from the triangle
 *             object.
 */
HashedVector<int,int> TriangleActiveDesignVariables(const triangle &triIn,
	const triangulation& triRSVS, const HashedVector<int, int> &objDvMap);

/**
 * Calculate the positional derivatives of a triangle object.
 *
 * This computes the partial derivatives $\frac{\partial coord}{\partial d}$
 * derivatives of snaxels and pseudo-centroids.
 *
 * @param[in]  triIn      The triangle for which the quantities need to be
 *                        computed.
 * @param[in]  triRSVS    The triangulation object from which this object is
 *                        part.
 * @param[in]  dvListMap  A hashed vector containing a list of active design
 *                        variables in the otriangle object as defined by
 *                        function TriangeActiveDesignVariables().
 * @param      dPos       The jacobian of the position derivatives.
 * @param      HPos       The Hessian of the position derivatives.
 */
void TrianglePositionDerivatives(const triangle &triIn, 
	const triangulation &triRSVS, const HashedVector<int,int> &dvListMap,
	Eigen::MatrixXd &dPos, Eigen::MatrixXd &HPos);



#endif // RSVSCALCTOOLS_HPP