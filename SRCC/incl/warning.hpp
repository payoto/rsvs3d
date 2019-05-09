/**
 * Provides the error and warning system used by the RSVS3D project.
 *  
 *@file
 */
#ifndef WARNING_H_INCLUDED
#define WARNING_H_INCLUDED

//===============================================
// Levels of debuging Guards
#ifdef DEBUGLVL2 // All Debugging calls
#define DEBUGLVL1

#endif

#ifdef DEBUGLVL1 // Debugging of new features.
#define TEST_WARNING
#endif

//=================================
// forward declared dependencies
// 		class foo; //when you only need a pointer not the actual object
// 		and to avoid circular dependencies



//=================================
// included dependencies
#include <iostream>
#include <stdarg.h>
#include <stdexcept>
#include <fstream>

//==================================
// Code
// NOTE: function in a class definition are IMPLICITELY INLINED 
//       ie replaced by their code at compile time

/**
 * @brief      Namespace for general purpose tools of the RSVS project.
 */
namespace rsvs3d {

	/**
	 * @brief      Exception for signaling rsvs errors.
	 */
	class rsvs_exception : public std::logic_error {
		using std::logic_error::logic_error;
	};

	/**
	 * Custom error function.
	 *
	 * Displays the name of the caller function and throw an exception type
	 * object with the message specified. can be turned off by setting the last
	 * parameter to false.
	 *
	 * @param[in]  message     Error message
	 * @param[in]  caller      Caller function
	 * @param[in]  file        The file in which the caller is.
	 * @param[in]  line        The line at which the caller is.
	 * @param[in]  throwError  should the error be thrown (True) or a warning
	 *                         (False)?
	 *
	 * @tparam     E           Exception type to throw
	 *
	 * Convenience macros are also provided to use this function without typing
	 * all the file, line and caller function macro names:
	 *  - RSVS3D_ERROR(M) : throws the default exception type (std::exception);
	 *  - RSVS3D_ERROR_LOGIC(M) : throws std::logic_error;
	 *  - RSVS3D_ERROR_ARGUMENT(M) : throws std::invalid_argument;
	 *  - RSVS3D_ERROR_TYPE(M, T) : throws T(M);
	 */
	template <class E=rsvs_exception>
	void error(const char* message="", const char* caller="",
		const char *file="", int line=0, bool throwError=true);
	/**
	 * @brief      Returns the sign of a type comparable to 0.
	 *
	 * @param[in]  val   	The value in.
	 *
	 * @tparam     T     	Any type comparable and with a 0 type.
	 *
	 * @return     	-1 for negative, 0 for 0, 1 for positive.
	 */
	template <typename T> int sign(T val);
	int TimeStamp(const char* str,int start_s);
	/**
	 * Returns a signed logscale useful for plotting data.
	 *
	 * This mathematical function is geared towards plotting of data with wildly
	 * varying orders of magnitude and has the following properties:
	 *  + Order preserving
	 *  + 0 preserving
	 *  + logarithmic for R not including [-DBL_EPSILON, DBL_EPSILON]
	 *  + 1/log10 in that range
	 *
	 * @param[in]  in    A double to scale.
	 *
	 * @return     Scaled value of the double.
	 */
	double SignedLogScale(double in);
	double Clock2ms(int clockCycles);
}
/**
 * @brief      Throw generic rsvs errors.
 *
 * @param      M     Message of the error (const char*).
 *
 * @throw     rsvs3d::rsvs_exception
 */
#define RSVS3D_ERROR(M) (rsvs3d::error(M, __PRETTY_FUNCTION__, __FILE__, __LINE__, true))

/**
 * @brief      Generic rsvs warning.
 *
 * @param      M     Message of the warning (const char*).
 */
#define RSVS3D_ERROR_NOTHROW(M) (rsvs3d::error(M, __PRETTY_FUNCTION__, __FILE__, __LINE__, false))

/**
 * @brief      Throw a specific error type.
 *
 * @param      M     Message of the warning (const char*).
 * @tparam      T     Type of the exception to throw.
 *
 * @throw     T
 */
#define RSVS3D_ERROR_TYPE(M, T) (rsvs3d::error<T>(M, __PRETTY_FUNCTION__, __FILE__, __LINE__, true))

/**
 * @brief      Throw a logic_error.
 *
 * @param      M     Message of the error (const char*).
 *
 * @throw     std::logic_error
 */
#define RSVS3D_ERROR_LOGIC(M) (rsvs3d::error<std::logic_error>(M, __PRETTY_FUNCTION__, __FILE__, __LINE__, true))

/**
 * @brief      Throw a invalid_argument.
 *
 * @param      M     Message of the error (const char*).
 *
 * @throw     std::invalid_argument
 */
#define RSVS3D_ERROR_ARGUMENT(M) (rsvs3d::error<std::invalid_argument>(M, __PRETTY_FUNCTION__, __FILE__, __LINE__, true))

#ifndef RSVS_NO_ARGCHECK
#define RSVS3D_ARGCHECK(E,M) if(!(E)){RSVS3D_ERROR_ARGUMENT(M);}
#define RSVS3D_ARGCHECK_WARN(E,M) if(!(E)){RSVS3D_ERROR_NOTHROW(M);}
#else
#define RSVS3D_ARGCHECK(E,M)
#define RSVS3D_ARGCHECK_WARN(E,M)
#endif
/**
 * @brief      Throw a range_error.
 *
 * @param      M     Message of the error (const char*).
 *
 * @throw     std::range_error
 */
#define RSVS3D_ERROR_RANGE(M) (rsvs3d::error<std::range_error>(M, __PRETTY_FUNCTION__, __FILE__, __LINE__, true))

void ThrowWarning(const char * message);
/**
 * @brief      Checks a file stream
 *
 * @param[in]  file      input or output file stream
 * @param[in]  callerID  the name of the caller function as given by pretty
 *                       function
 * @param[in]  fileName  The name of the file opened in the stream.
 *
 * @tparam     T         either ifstream or ofstream, needs to support method
 *                       `T::is_open()`
 */
template<class T>
void CheckFStream(const T &file, const char* callerID, 
	const std::string &fileName){
	if (!file.is_open()){
		std::string errstr;
		errstr = "Parameter file failed to open ";
		errstr += "in " ;
		errstr += callerID;
		errstr += " \n:  " + fileName;
		RSVS3D_ERROR_ARGUMENT(errstr.c_str());
	}
}

namespace rsvs3d {

	template <class E>
	void error(const char* message, const char* caller,
		const char *file, int line, bool throwError){
		std::cerr << std::endl << "Error at: " << file << ":" << line 
			<< std::endl << "    " <<  caller << std::endl;
		if (throwError){
			throw E(message);
		} else {
			std::cerr << "Exception (not thrown)" << message << std::endl;
		}
	}
	template <typename T> int sign(T val) {
	    return (T(0) < val) - (val < T(0));
	}

}
#endif // WARNING_H_INCLUDED