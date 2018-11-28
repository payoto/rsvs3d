#ifndef RSVSINTEGRATION_H_INCLUDED 
#define RSVSINTEGRATION_H_INCLUDED 


//=================================
// forward declared dependencies
// 		class foo; //when you only need a pointer not the actual object
// 		and to avoid circular dependencies
namespace integrate {
	class RSVSclass;
}
class triangulation;
class snake;
class mesh;
class tecplotfile;
namespace param {
	class grid;
	class snaking;
	class parameters;
}

//=================================
// included dependencies
#include <vector>

// ================================
// declarations

using namespace std; 


void SnakeConnectivityUpdate(snake &testSnake,  vector<int> &isImpact);
void SnakeConnectivityUpdate_2D(snake &testSnake,  vector<int> &isImpact);
void SnakeConnectivityUpdate_legacy(snake &snakein,  vector<int> &isImpact);
void SnakeConnectivityUpdate_robust(snake &snakein,  vector<int> &isImpact);
int TimeStamp(const char* str,int start_s);


namespace integrate {
	void Prepare(RSVSclass RSVSobj);
	namespace prepare {
		void Mesh(
			const param::grid &gridconf,
			mesh &snakeMesh,
			mesh &voluMesh
		);
		void Snake(
			const param::snaking &snakconf, 
			mesh &snakeMesh,
			snake &rsvsSnake
		);
		void Triangulation(
			mesh &snakeMesh,
			snake &rsvsSnake,
			triangulation &rsvsTri
			);
		void PostProcessing(
			const param::parameters &paramconf,
			const param::parameters &origcong,
			tecplotfile &outSnake
			);
	}

	namespace execute{
		void PostProcessing();
		void RSVSiterate();
	}

	namespace test {
		int Prepare();
	}
}

//==================================
// Code
// NOTE: function in a class definition are IMPLICITELY INLINED 
//       ie replaced by their code at compile time
	

#endif //RSVSINTEGRATION_H_INCLUDED 