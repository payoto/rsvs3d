#ifndef PARAMETERS2JSON_H_INCLUDED 
#define PARAMETERS2JSON_H_INCLUDED 



//=================================
// forward declared dependencies
// 		class foo; //when you only need a pointer not the actual object
// 		and to avoid circular dependencies
namespace param{
	template <class T> struct filltype;
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
// parameter structure to json

using nlohmann::json;
namespace param {
	template<class T>
	void to_json(json& j, const filltype<T>& p);
	template<class T>
	void from_json(const json& j, filltype<T>& p); 

	void to_json(json& j, const rsvs& p);
	void from_json(const json& j, rsvs& p); 	

	void to_json(json& j, const snaking& p);
	void from_json(const json& j, snaking& p); 

	void to_json(json& j, const voxel& p);
	void from_json(const json& j, voxel& p);

	void to_json(json& j, const voronoi& p);
	void from_json(const json& j, voronoi& p);

	void to_json(json& j, const grid& p);
	void from_json(const json& j, grid& p); 
	
	void to_json(json& j, const parameters& p);
	void from_json(const json& j, parameters& p); 

	void to_json(json& j, const ioin& p);
	void from_json(const json& j, ioin& p); 

	void to_json(json& j, const ioout& p);
	void from_json(const json& j, ioout& p); 

	void to_json(json& j, const files& p);
	void from_json(const json& j, files& p); 

	void flatupdate(json& jfin, json& jnew,
		bool isFlatFin, bool isFlatNew);
}

#endif