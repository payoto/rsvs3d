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
class RSVScalc;
namespace param {
	class grid;
	class snaking;
	class parameters;
}

//=================================
// included dependencies
#include <vector>
#include <fstream>

// ================================
// declarations

using namespace std; 


void SnakeConnectivityUpdate(snake &testSnake,  vector<int> &isImpact,
	double impactAlmostRange=0.2);
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
		void Output(
			const param::parameters &paramconf,
			const param::parameters &origcong,
			tecplotfile &outSnake,
			std::ofstream &logFile,
			std::ofstream &coutFile,
			std::ofstream &cerrFile
			);
	}

	namespace execute{
		void RSVSiterate(RSVSclass &RSVSobj);
		void Logging(RSVSclass &RSVSobj,
			double totT, int nVoluZone, int stepNum);
		void PostProcessing(RSVSclass &RSVSobj,
			double totT, int nVoluZone, int stepNum);

		namespace logging{
			// Log only convergence data
			void Log(
				std::ofstream &logFile,
				RSVScalc &calcObj,
				int loglvl);
			// Store the Snake in a file
			void Snake(
				tecplotfile &outSnake,
				snake &rsvsSnake,
				mesh &voluMesh,
				double totT, int nVoluZone);
			// Store All the available information
			void FullTecplot(
				tecplotfile &outSnake, snake &rsvsSnake,
				triangulation &rsvsTri, mesh &voluMesh,
				double totT, int nVoluZone, int stepNum);
		}

		namespace postprocess{
			// Log only convergence data
			// Log only convergence data
			void Log(
				std::ofstream &logFile,
				RSVScalc &calcObj,
				int loglvl);
			// Store the Snake in a file
			void Snake(
				tecplotfile &outSnake, snake &rsvsSnake,
				mesh &voluMesh, double totT, int nVoluZone,
				param::parameters &paramconf);
			// Store All the available information
			void FullTecplot(
				tecplotfile &outSnake, snake &rsvsSnake,
				triangulation &rsvsTri, mesh &voluMesh,
				double totT, int nVoluZone, int stepNum);
		}
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