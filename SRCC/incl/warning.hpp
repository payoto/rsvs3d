
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

using namespace std;
namespace rsvs3d {
	class rsvs_exception : public std::logic_error {
		using std::logic_error::logic_error;
	};

	template <class E=rsvs_exception>
	void error(const char* message="", const char* caller="",
		const char *file="", int line=0, bool throwError=true);
}
#define RSVS3D_ERROR(M) (rsvs3d::error(M, __PRETTY_FUNCTION__, __FILE__, __LINE__, true))

#define RSVS3D_ERROR_TYPE(M, T) (rsvs3d::error<T>(M, __PRETTY_FUNCTION__, __FILE__, __LINE__, true))
#define RSVS3D_ERROR_LOGIC(M) (rsvs3d::error<std::logic_error>(M, __PRETTY_FUNCTION__, __FILE__, __LINE__, true))
#define RSVS3D_ERROR_ARGUMENT(M) (rsvs3d::error<std::invalid_argument>(M, __PRETTY_FUNCTION__, __FILE__, __LINE__, true))
#define RSVS3D_ERROR_RANGE(M) (rsvs3d::error<std::range_error>(M, __PRETTY_FUNCTION__, __FILE__, __LINE__, true))

void ThrowWarning(const char * message);
/**
 * @brief      Checks a file stream
 *
 * @param[in]  file      input or output file stream
 * @param[in]  callerID  the name of the caller function as given by pretty
 * function
 * @param[in]  fileName  The name of the file opened in the stream.
 *
 * @tparam     T         either ifstream or ofstream, needs to support method
 * `T::is_open()`
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

	/**
	 * Custom error function. 
	 * 
	 * Displays the name of the caller
	 * function and throw an exception type object with the message
	 * specified. can be turned off by setting the last parameter to
	 * false.
	 *
	 * @param[in]  message     Error message
	 * @param[in]  caller      Caller function
	 * @param[in]  file        The file
	 * @param[in]  line        The line
	 * @param[in]  throwError  should the error be thrown (True) or a warning
	 *                         (False)?
	 *
	 * @tparam     E           Exception type to throw
	 *
	 * Convenience macros are also provided to use this function
	 * without typing all the file, line and caller function macro
	 * names:
	 *  - RSVS3D_ERROR(M) : throws the default exception type (std::exception);
	 *  - RSVS3D_ERROR_LOGIC(M) : throws std::logic_error;
	 *  - RSVS3D_ERROR_ARGUMENT(M) : throws std::invalid_argument;
	 *  - RSVS3D_ERROR_TYPE(M, T) : throws T(M);
	 */
	template <class E>
	void error(const char* message, const char* caller,
		const char *file, int line, bool throwError){
		cerr << endl << "Error at: " << file << ":" << line 
			<< endl << "    " <<  caller << endl;
		if (throwError){
			throw E(message);
		} else {
			cerr << "Exception (not thrown)" << message << endl;
		}
	}


	//(?s)throw\s+([a-zA-Z:_]*)\s*\(([^;]*)\);
}
#endif // WARNING_H_INCLUDED