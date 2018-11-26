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

//==================================
// Code
// NOTE: function in a class definition are IMPLICITELY INLINED 
//       ie replaced by their code at compile time

namespace param {

	template <class T>
	struct bounds{
		T lb;
		T ub;
	};

	class RSVS
	{
	public:
		RSVS();
		~RSVS();
		
	};

	class snake
	{
	public:
		snake();
		~snake();
		
	};
	class voxel
	{
	public:
		std::array<bounds<double>, 3>  domain;
		std::array<int, 3>  gridsizebackground;
		std::array<int, 3> gridsizesnake;

		voxel();
		~voxel();
		
	};

	class parameters
	{
	public:
		parameters();
		~parameters();
		
	};

}

#endif