/**
 * Integration into the full 3 dimensional Restricted Snake Volume of Solid
 * method.
 *  
 *@file
 */


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
	class rsvs;
	class voxel;
	class voronoi;
	class ioin;
	class ioout;
	class tecplottemplate;
}

//=================================
// included dependencies
#include <vector>
#include <fstream>
#include <tuple>

// ================================
// declarations

using namespace std; 


void SnakeConnectivityUpdate(snake &testSnake,  vector<int> &isImpact,
	double impactAlmostRange=0.2);
void SnakeConnectivityUpdate_2D(snake &testSnake,  vector<int> &isImpact);
void SnakeConnectivityUpdate_legacy(snake &snakein,  vector<int> &isImpact);
void SnakeConnectivityUpdate_robust(snake &snakein,  vector<int> &isImpact);
double SnakePositionUpdate(snake &rsvsSnake, std::vector<double> &dt,
	double snaxtimestep, double snaxdiststep);
namespace integrate {
	namespace constants {
		static const std::string tecplotsnake = "rsvs3D_";
		static const std::string tecplotgradient = "rsvsgradients3D_";
	}
	class iteratereturns {
	public:
		int nVoluZone=0;
		int stepNum=0;
		double timeT=0.0;

		iteratereturns(int n, int s, double t){
			this->nVoluZone=n;
			this->stepNum=s;
			this->timeT=t;
		}

	};

	void Prepare(RSVSclass &RSVSobj);
	namespace prepare {
		void Mesh(
			const param::grid &gridconf,
			const param::ioin &ioinconf,
			mesh &snakeMesh,
			mesh &voluMesh
		);
		void Snake(
			const param::snaking &snakconf, 
			const param::rsvs &rsvsconf, 
			const param::ioin &ioinconf, 
			mesh &snakeMesh,
			mesh &voluMesh,
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
			tecplotfile &outgradientsnake,
			std::ofstream &logFile,
			std::ofstream &coutFile,
			std::ofstream &cerrFile
			);
		namespace grid {
			void Voxel(
				const param::grid &gridconf,
				mesh &snakeMesh,
				mesh &voluMesh
			);
			void Voronoi(
				const param::grid &gridconf,
				mesh &snakeMesh,
				mesh &voluMesh
			);
			void Load(
				const param::ioin &ioinconf,
				mesh &snakeMesh,
				mesh &voluMesh
			);
		}
	}

	namespace execute{

		void All(integrate::RSVSclass &RSVSobj);
		iteratereturns RSVSiterate(RSVSclass &RSVSobj);
		void Logging(RSVSclass &RSVSobj,
			double totT, int nVoluZone, int stepNum);
		void PostProcessing(RSVSclass &RSVSobj,
			double totT, int nVoluZone, int stepNum);
		void Exporting(RSVSclass &RSVSobj);

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
			void Gradients(const RSVScalc &calcObj, 
				const triangulation &rsvsTri,
				tecplotfile &outgradientsnake, 
				double totT);
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
				snake &rsvsSnake, mesh &voluMesh,
				param::parameters &paramconf);
			// Store All the available information
			void FullTecplot(
				tecplotfile &outSnake, snake &rsvsSnake,
				triangulation &rsvsTri, mesh &voluMesh,
				double totT, int nVoluZone, int stepNum);
			void Gradients(const RSVScalc &calcObj, 
				const triangulation &rsvsTri,
				tecplotfile &outgradientsnake, 
				double totT);
		}

		namespace exporting {
			void SU2(
				std::string exportStr,
				snake &rsvsSnake,
				param::parameters &paramconf);

		}
	}
	namespace utils {

		/**
		 * @brief      Convenience function to generate file names for RSVS
		 *
		 * @param[in]  paramconf  The current parameter configuration
		 * @param[in]  fileName   The start of the file name
		 * @param[in]  extension  The file extension
		 *
		 * @return     Returns a file path by adding the root directory and the
		 *             file pattern
		 */
		std::string OutputFileName(const param::parameters &paramconf, 
			std::string fileName, std::string extension);
		/**
		 * @brief      Convenience function to generate file names for RSVS
		 *
		 * @param[in]  rootDirectory  The root directory in which the file will
		 *                            be stored.
		 * @param[in]  filePattern    The file diferentiating pattern
		 * @param[in]  fileName       The core file name to be placed at teh
		 *                            start of the basename
		 * @param[in]  extension      The file extension
		 *
		 * @return     A string with a file path for RSVS output.
		 */
		std::string OutputFileName(const std::string rootDirectory, 
			const std::string &filePattern,
			std::string fileName, std::string extension);


		void WriteModifiedTemplate(const std::string &fileIn, 
			const std::string &fileOut, const std::string &oldLine, 
			const std::string newLine);

		void SpecialiseTemplateFiles(const param::parameters &paramconf);
		void SpecialiseTemplateFile(const param::tecplottemplate &tecconfig, int logLvl,
			const param::ioout &ioout, std::string fileName);
	}
	namespace test {
		int Prepare();
		int All();
		int CompareSurfCentreDerivatives();
		int CompareDerivativesSpike();
		int CompareDerivativesSpikeNoDPos();
		int StudyDerivatives();
	}
}


//==================================
// Code
// NOTE: function in a class definition are IMPLICITELY INLINED 
//       ie replaced by their code at compile time
	

#endif //RSVSINTEGRATION_H_INCLUDED 