#include "SQPstep.hpp"

#include <Eigen>
#include <cmath>

using namespace Eigen;

void ResizeLagrangianMultiplier(const RSVScalc &calcobj, VectorXd &lagMultAct, bool &isLarge, bool &isNan)
{
    /*
    Resizes the lagrangian multiplier LagMultAct based on whether any of its
    values are nan or too large.

    This uses the RSVScalc object to guide it.
    */
    int ii, ni;

    isLarge = false;
    isNan = false;
    ni = lagMultAct.size();
    for (ii = 0; ii < ni; ++ii)
    {
        if (lagMultAct[ii] < -calcobj.limLag)
        {
            lagMultAct[ii] = -calcobj.limLag;
            isLarge = true;
        }
        else if (lagMultAct[ii] > calcobj.limLag)
        {
            lagMultAct[ii] = calcobj.limLag;
            isLarge = true;
        }
        else if (std::isnan(lagMultAct[ii]))
        {
            lagMultAct[ii] = 0.0;
            isNan = true;
        }
    }
}
