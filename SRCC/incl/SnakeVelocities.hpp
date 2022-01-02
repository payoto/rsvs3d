/**
 * This file contains declarations of the valid velocity entries.
 *
 * @file SnakeVelocities.hpp
 * @author Alexandre Payot (alexandre.payotli@gmail.com)
 * @brief
 * @version 0.1
 * @date 2021-12-31
 *
 * @copyright Copyright (c) 2021
 *
 */

#ifndef SNAKEVELOCITIES_H_INCLUDED
#define SNAKEVELOCITIES_H_INCLUDED

class snake;

#include <functional>

#include "SnakeVelocityCalculator.hpp"
namespace rsvs3d
{
typedef std::function<void(snake &)> simplevelocity;
class VelocityFunction : public rsvs3d::SnakeVelocityCalculator
{
  private:
    simplevelocity activeFunction;

  public:
    void setDevParameters(const param::dev::devparam &) override;
    void ConvergenceLog(std::ofstream &, int = 3) const override;
    void CalculateVelocities(triangulation &triRSVS, int calculationMethod = 0, bool calculateDerivatives = true,
                             int derivativeMethod = 1) override;
    VelocityFunction(simplevelocity activeFunction) : activeFunction(activeFunction)
    {
    }
};
} // namespace rsvs3d
#endif // SNAKEVELOCITIES_H_INCLUDED
