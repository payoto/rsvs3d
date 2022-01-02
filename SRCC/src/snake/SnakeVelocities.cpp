#include "SnakeVelocities.hpp"

#include "snake.hpp"
#include "triangulate.hpp"

using namespace rsvs3d;

void VelocityFunction::setDevParameters(const param::dev::devparam &devparam)
{
    // For a function this is a no-op
}

void VelocityFunction::ConvergenceLog(std::ofstream &, int) const
{
    // For a function this is a no-op
}

void VelocityFunction::CalculateVelocities(triangulation &triRSVS, int calculationMethod, bool calculateDerivatives,
                                           int derivativeMethod)
{
    this->activeFunction(*triRSVS.snakeDep);
}
