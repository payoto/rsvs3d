/**
 * Provides the infrastructure for calculation of the RSVS equations.
 *
 *@file
 */

#ifndef SQPSTEP_H_INCLUDED
#define SQPSTEP_H_INCLUDED

//=================================
// forward declared dependencies
// 		class foo; //when you only need a pointer not the actual object
// 		and to avoid circular dependencies
#include <Eigen>
#include <fstream>
#include <iostream>
#include <vector>

#include "RSVScalc.hpp"
#include "arraystructures.hpp"
#include "mesh.hpp"
#include "snake.hpp"
#include "triangulate.hpp"
#include "warning.hpp"

//=================================
// included dependencies
//==================================
// Code
// NOTE: function in a class definition are IMPLICITELY INLINED
//       ie replaced by their code at compile time

/**
 * Resizes the lagrangian multiplier LagMultAct based on whether any of its
 * values are nan or too large.

 This uses the RSVScalc object to guide the resizing operation if it is needed.

 @param[in]     calcobj     The calculation object.
 @param[in,out] lagMultAct  The std::vector of active lagrangian multipliers.
 @param[out]    isLarge     Returns if lagMultAct is too large.
 @param[out]    isNan       Returns if lagMultAct has Nan values.
*/
void ResizeLagrangianMultiplier(const RSVScalc &calcobj, Eigen::VectorXd &lagMultAct, bool &isLarge, bool &isNan);

/**
 * Template for calculation of an SQP step.
 *
 * This template cannot be deduced and needs the developer to pass the required
 * solver class when it is called.
 *
 * Instantiation options: Eigen::HouseholderQR Eigen::ColPivHouseholderQR
 * Eigen::LLT<Eigen::MatrixXd> (*) <- needs a full type to be defined (see
 * below) Eigen::PartialPivLU
 *
 * For stability info
 * https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html
 * https://eigen.tuxfamily.org/dox/group__DenseDecompositionBenchmark.html
 *
 * @param[in]  calcobj            The calculation object
 * @param[in]  dConstrAct         Active constraint jacobian dh/dx
 * @param[in]  dObjAct            Active objective jacobian dJ/dx
 * @param[in]  constrAct          Active constraint vector
 * @param      lagMultAct         The active lagrangian multipliers
 * @param      deltaDVAct         The active SQP step to take
 * @param[out] isNan              Indicates if lagMult is nan
 * @param[out] isLarge            Indicates if lagMult is large
 * @param[in]  attemptConstrOnly  Should the step algorithm attempt using only
 *                                the constraint to step.
 *
 * @tparam     T                  The Eigen object template type to use. A full
 *                                type will be defined using T<Eigen::MatrixXd>.
 *
 * @return     (isLarge || isNan), if true some form of failure was detected.
 */
template <class T>
bool SQPstep(const RSVScalc &calcobj, const Eigen::MatrixXd &dConstrAct, const Eigen::RowVectorXd &dObjAct,
             const Eigen::VectorXd &constrAct, Eigen::VectorXd &lagMultAct, Eigen::VectorXd &deltaDVAct, bool &isNan,
             bool &isLarge, bool attemptConstrOnly = true);

/**
 * Template for calculation of an SQP step.
 *
 * This template cannot be deduced and needs the developer to pass the required
 * solver class when it is called.
 *
 * Instantiation options: Eigen::HouseholderQR<Eigen::MatrixXd>
 * Eigen::ColPivHouseholderQR<Eigen::MatrixXd> Eigen::LLT<Eigen::MatrixXd>
 * Eigen::PartialPivLU<Eigen::MatrixXd>
 *
 * For stability info
 * https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html
 * https://eigen.tuxfamily.org/dox/group__DenseDecompositionBenchmark.html
 *
 * @param[in]  calcobj            The calculation object
 * @param[in]  dConstrAct         Active constraint jacobian dh/dx
 * @param[in]  dObjAct            Active objective jacobian dJ/dx
 * @param[in]  constrAct          Active constraint vector
 * @param      lagMultAct         The active lagrangian multipliers
 * @param      deltaDVAct         The active SQP step to take
 * @param[out] isNan              Indicates if lagMult is nan
 * @param[out] isLarge            Indicates if lagMult is large
 * @param[in]  attemptConstrOnly  Should the step algorithm attempt using only
 *                                the constraint to step.
 *
 * @tparam     T                  The Eigen object type to use. Should take a
 *                                RSVScalc::HLag as a constructor and support a
 *                                solve method.
 *
 * @return     (isLarge || isNan), if true some form of failure was detected.
 */
template <template <typename> class T>
bool SQPstep(const RSVScalc &calcobj, const Eigen::MatrixXd &dConstrAct, const Eigen::RowVectorXd &dObjAct,
             const Eigen::VectorXd &constrAct, Eigen::VectorXd &lagMultAct, Eigen::VectorXd &deltaDVAct, bool &isNan,
             bool &isLarge, bool attemptConstrOnly = true);

/**
 * Calculation for the sparse SQP sensitivity.
 *
 * @param      sensMult  The sensitivity multiplier.
 * @param[in]  sensInv   The Matrix to invert.
 * @param      sensRes   The sensitivity result.
 *
 * @tparam     T         The Eigen object type to use. Should take a
 *                       RSVScalc::HLag as a constructor and support a solve
 *                       method.
 *
 * @return     true.
 */
template <class T>
bool SQPsens_sparse(const Eigen::MatrixXd &sensMult, const MatrixXd_sparse &sensInv, Eigen::MatrixXd &sensRes);

template <class T>
bool SQPsens(const Eigen::MatrixXd &sensMult, const Eigen::MatrixXd &sensInv, Eigen::MatrixXd &sensRes);

template <template <typename> class T>
bool SQPsens(const Eigen::MatrixXd &sensMult, const Eigen::MatrixXd &sensInv, Eigen::MatrixXd &sensRes);

// Code to be included as templated functions

template <template <typename> class T>
bool SQPstep(const RSVScalc &calcobj, const Eigen::MatrixXd &dConstrAct, const Eigen::RowVectorXd &dObjAct,
             const Eigen::VectorXd &constrAct, Eigen::VectorXd &lagMultAct, Eigen::VectorXd &deltaDVAct, bool &isNan,
             bool &isLarge, bool attemptConstrOnly)
{
    return (SQPstep<T<Eigen::MatrixXd>>(calcobj, dConstrAct, dObjAct, constrAct, lagMultAct, deltaDVAct, isNan, isLarge,
                                        attemptConstrOnly));
}

template <class T>
bool SQPstep(const RSVScalc &calcobj, const Eigen::MatrixXd &dConstrAct, const Eigen::RowVectorXd &dObjAct,
             const Eigen::VectorXd &constrAct, Eigen::VectorXd &lagMultAct, Eigen::VectorXd &deltaDVAct, bool &isNan,
             bool &isLarge, bool attemptConstrOnly)
{
    Eigen::LDLT<Eigen::MatrixXd> condsys(calcobj.HLag);
    Eigen::MatrixXd temp1, temp2;
    T HLagSystem(calcobj.HLag);

    std::cout << " (rcond) " << condsys.rcond();

    temp1 = HLagSystem.solve(dConstrAct.transpose());
    temp2 = HLagSystem.solve(dObjAct.transpose());

    T LagMultSystem(dConstrAct * (temp1));

    lagMultAct = LagMultSystem.solve(constrAct - (dConstrAct * (temp2)));

    ResizeLagrangianMultiplier(calcobj, lagMultAct, isLarge, isNan);
    isLarge = false;
    if (!attemptConstrOnly && (isLarge || isNan))
    {
        std::cout << "(early sqp return) ";
        return (isLarge || isNan);
    }
    if (isLarge || isNan)
    {
        // Use a least squared solver if only using the constraint
        std::cout << "(constrmov) ";
        deltaDVAct = -dConstrAct.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(constrAct);
    }
    else
    {
        deltaDVAct = -(HLagSystem.solve(dObjAct.transpose() + dConstrAct.transpose() * lagMultAct));
    }
    return (isLarge || isNan);
}

template <class T>
bool SQPstep_sparse(const RSVScalc &calcobj, const MatrixXd_sparse &dConstrAct, const Eigen::RowVectorXd &dObjAct,
                    const Eigen::VectorXd &constrAct, Eigen::VectorXd &lagMultAct, Eigen::VectorXd &deltaDVAct,
                    bool &isNan, bool &isLarge, bool attemptConstrOnly)
{
    // Eigen::LDLT<MatrixXd_sparse> condsys(calcobj.HLag_sparse);
    T HLagSystem;
    HLagSystem.compute(calcobj.HLag_sparse);
    // std::cout << " (rcond) " << condsys.rcond();
    if (HLagSystem.info() != Eigen::Success)
    {
        // decomposition failed
        RSVS3D_ERROR_NOTHROW("Failed to decompose Hessian of lagrangian.");
        isNan = true;
        return (isLarge || isNan);
    }
    T LagMultSystem;
    MatrixXd_sparse temp1 = dConstrAct.transpose();
    temp1.makeCompressed();
    MatrixXd_sparse temp2 = dConstrAct * HLagSystem.solve(temp1);
    temp2.makeCompressed();
    LagMultSystem.compute(temp2);
    if (LagMultSystem.info() != Eigen::Success)
    {
        // decomposition failed
        RSVS3D_ERROR_NOTHROW("Failed to decompose Lagrangian multiplier system.");
        isNan = true;
        return (isLarge || isNan);
    }
    lagMultAct = LagMultSystem.solve(constrAct - (dConstrAct * HLagSystem.solve(dObjAct.transpose())));

    ResizeLagrangianMultiplier(calcobj, lagMultAct, isLarge, isNan);
    isLarge = false;
    if (!attemptConstrOnly && (isLarge || isNan))
    {
        std::cout << "(early sqp return) ";
        return (isLarge || isNan);
    }
    if (isLarge || isNan)
    {
        // Use a least squared solver if only using the constraint
        std::cout << "(constrmov) ";
        Eigen::SparseQR<MatrixXd_sparse, Eigen::COLAMDOrdering<int>> gradOnlySolver;
        gradOnlySolver.compute(dConstrAct);
        deltaDVAct = -gradOnlySolver.solve(constrAct);
    }
    else
    {
        deltaDVAct = -(HLagSystem.solve(dObjAct.transpose() + dConstrAct.transpose() * lagMultAct));
    }
    return (isLarge || isNan);
}

template <template <typename> class T>
bool SQPsens(const Eigen::MatrixXd &sensMult, const Eigen::MatrixXd &sensInv, Eigen::MatrixXd &sensRes)
{
    return SQPsens<T<Eigen::MatrixXd>>(sensMult, sensInv, sensRes);
}

template <class T>
bool SQPsens(const Eigen::MatrixXd &sensMult, const Eigen::MatrixXd &sensInv, Eigen::MatrixXd &sensRes)
{
    T HLagSystem(sensInv);

    sensRes = -HLagSystem.solve(sensMult);
    return (true);
}

template <class T>
bool SQPsens_sparse(const Eigen::MatrixXd &sensMult, const MatrixXd_sparse &sensInv, Eigen::MatrixXd &sensRes)
{
    T HLagSystem;
    HLagSystem.compute(sensInv);
    if (HLagSystem.info() != Eigen::Success)
    {
        // decomposition failed
        RSVS3D_ERROR_NOTHROW("Failed to decompose sensitivity system.");
        return (false);
    }
    sensRes = -(HLagSystem.solve(sensMult));
    return (true);
}

#endif
