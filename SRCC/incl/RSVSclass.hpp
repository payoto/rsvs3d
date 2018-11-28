#ifndef RSVSCLASS_H_INCLUDED 
#define RSVSCLASS_H_INCLUDED 


//=================================
// forward declared dependencies
// 		class foo; //when you only need a pointer not the actual object
// 		and to avoid circular dependencies



//=================================
// included dependencies

#include "mesh.hpp"
#include "snake.hpp"
#include "postprocessing.hpp"
#include "parameters.hpp"
#include "snakevel.hpp"

//==================================
// Code
// This file defines a god class for conveniently doing the integration
// of the various compoments into a cohesive framework
//
namespace integrate {
	class RSVSclass {
	public:
		param::parameters paramconf;
		tecplotfile outSnake;
		snake rsvsSnake;
		mesh snakeMesh;
		mesh voluMesh;
		triangulation rsvsTri;
	};
}

#endif