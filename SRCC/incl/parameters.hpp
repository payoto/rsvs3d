#ifndef PARAMETERS_H_INCLUDED 
#define PARAMETERS_H_INCLUDED 



//=================================
// forward declared dependencies
// 		class foo; //when you only need a pointer not the actual object
// 		and to avoid circular dependencies

//=================================
// included dependencies

#include <cstdlib>
#include <array>
#include "json.hpp"

//==================================
// Code
// 
// This file defines parameters to be used in other parts of the 
// RSVS snaking framework. Default values are defined in "parameters.cpp"
//
// Substructure names are all 4-5 letters
using nlohmann::json;
namespace param {

	template <class T>
	struct bounds{
		T lb;
		T ub;
	};
	template <class T>
	void to_json(json& j, const bounds<T>& p);
	template <class T>
	void from_json(const json& j, bounds<T>& p); 

	/*
	Parameters related to the Velocity calculation and VOS steps
	*/
	class rsvs
	{
	public:
		// Algorithm used by Eigen to solve the SQP system.
		int solveralgorithm;

		rsvs();
		~rsvs();
		
	};
	void to_json(json& j, const rsvs& p);
	void from_json(const json& j, rsvs& p); 

	/*
	Parameters controlling tuning parameters for the stepping of
	the restricted surface.
	*/
	class snaking
	{
	public:
		// Distance along edge at which a vertex is considered
		// arrived regardless of "d" and "v"
		double arrivaltolerance;
		// Distance along edge at which converging snaxels are
		// considered arrived
		double multiarrivaltolerance;
		// Initialisation boundary (either 0 or 1)
		int initboundary;
		// maximum number of steps
		int maxsteps;

		snaking();
		~snaking();
		
	};
	void to_json(json& j, const snaking& p);
	void from_json(const json& j, snaking& p); 

	/*
	Paramters controlling grid properties
	*/
	class voxel
	{
	public:
		std::array<bounds<double>, 3>  domain;
		std::array<int, 3>  gridsizebackground;
		std::array<int, 3> gridsizesnake;

		voxel();
		~voxel();
		
	};
	void to_json(json& j, const voxel& p);
	void from_json(const json& j, voxel& p); 

	class grid
	{
	public:
		voxel voxel;
		// grid();
		// ~grid();
		
	};
	void to_json(json& j, const grid& p);
	void from_json(const json& j, grid& p); 

	class parameters
	{
	public:
		rsvs rsvs;
		snaking snak;
		grid grid;
		// parameters();
		// ~parameters();
		
	};
	void to_json(json& j, const parameters& p);
	void from_json(const json& j, parameters& p); 

	int test();
}

#endif