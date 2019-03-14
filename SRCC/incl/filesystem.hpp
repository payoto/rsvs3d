// Custom filesystem header Faff about with filesystem depending on version
// To give a readable compile time error if incompatible things are attempted.
// Basically a workaround gcc on windows hating me.

// test for c++17 which has its own filesystem header
#if __cplusplus == 201703L
	#if defined(__GNUC__) && (defined(__MINGW32__) || defined(__MINGW64__))
		#if !(defined(USE_BOOST) || defined(USE_CSTD_FILESYSTEM)) 
			#if (__GNUC__==8)
				#error gcc v8.1 does not currently support <filesystem> header  windows (BUG) use <boost/filesystem.hpp>. Alternatively force use of the header by adding flag "-DUSE_CSTD_FILESYSTEM" to your compile command.
			#elif (__GNUC__==7)
				#define USE_CSTD_FILESYSTEM_EXPERIMENTAL
			#else
				// Not even sure we can get there	
				#error gcc <= v6 does not  support <filesystem> header: an addition of c++17. Please install and use boost <boost/filesystem.hpp>.
			#endif
		#endif
	#else
		#define USE_CSTD_FILESYSTEM
	#endif
#endif 
#ifdef USE_CSTD_FILESYSTEM
	#include <filesystem>
#elif defined(USE_CSTD_FILESYSTEM_EXPERIMENTAL)
	#include <experimental/filesystem>
#else
	// else use boost
	#include <boost/filesystem.hpp>
#endif

#ifdef USE_CSTD_FILESYSTEM
	namespace filesystem=std::filesystem;
#else
	// else use boost
	namespace filesystem=boost::filesystem;
#endif