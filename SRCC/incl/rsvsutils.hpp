/**
 * Provides various utility functions for the RSVS operation
 *  
 *@file
 */

//===============================================
// Include Guards
#ifndef RSVSUTILS_H_INCLUDED
#define RSVSUTILS_H_INCLUDED

#include <cmath>
#include <cfloat>

namespace rsvs3d {
	namespace utils {
		inline bool IsAproxEqual(double d1,double d2, 
			double tol=DBL_EPSILON)	{return(fabs(d1-d2)<tol);}
	}
}


#endif // RSVSUTILS_H_INCLUDED
