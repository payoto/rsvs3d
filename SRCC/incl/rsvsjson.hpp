/**
 * Interface between the RSVS project and the 
 * [JSON for Modern C++ library](https://github.com/nlohmann/json).
 *  
 *@file
 */



#ifndef PARAMETERS2JSON_H_INCLUDED 
#define PARAMETERS2JSON_H_INCLUDED 

//=================================
// forward declared dependencies
// 		class foo; //when you only need a pointer not the actual object
// 		and to avoid circular dependencies
namespace param{
	template <class T> struct filltype;
	class outputtemplate;
	class rsvs;
	class snaking;
	class voxel;
	class voronoi;
	class grid;
	class parameters;
	class ioin;
	class ioout;
	class files;
}

namespace tetgen{
	class apiparam;
}

//=================================
// included dependencies

#include <cstdlib>
#include <array>
#include <string>
#include "json.hpp"

//==================================
// Code
// 
// This file defines the functions required to convert the
// various classes to json

namespace rsvsjson{
	using json = nlohmann::json;
}
namespace param {

	template<class T>
	void to_json(rsvsjson::json& j, const filltype<T>& p);
	template<class T>
	void from_json(const rsvsjson::json& j, filltype<T>& p);

	void to_json(rsvsjson::json& j, const outputtemplate& p);
	void from_json(const rsvsjson::json& j, outputtemplate& p); 	 

	void to_json(rsvsjson::json& j, const rsvs& p);
	void from_json(const rsvsjson::json& j, rsvs& p); 	

	void to_json(rsvsjson::json& j, const snaking& p);
	void from_json(const rsvsjson::json& j, snaking& p); 

	void to_json(rsvsjson::json& j, const voxel& p);
	void from_json(const rsvsjson::json& j, voxel& p);

	void to_json(rsvsjson::json& j, const voronoi& p);
	void from_json(const rsvsjson::json& j, voronoi& p);

	void to_json(rsvsjson::json& j, const grid& p);
	void from_json(const rsvsjson::json& j, grid& p); 
	
	void to_json(rsvsjson::json& j, const parameters& p);
	void from_json(const rsvsjson::json& j, parameters& p); 

	void to_json(rsvsjson::json& j, const ioin& p);
	void from_json(const rsvsjson::json& j, ioin& p); 

	void to_json(rsvsjson::json& j, const ioout& p);
	void from_json(const rsvsjson::json& j, ioout& p); 

	void to_json(rsvsjson::json& j, const files& p);
	void from_json(const rsvsjson::json& j, files& p); 
}
namespace rsvsjson {
	void flatupdate(json& jfin, json& jnew,
		bool isFlatFin, bool isFlatNew);
}
namespace tetgen {

	void to_json(rsvsjson::json& j, const apiparam& p);
	void from_json(const rsvsjson::json& j, apiparam& p); 
}


#endif