/**
 * Parameters for the integrated 3DRSVS.
 *  
 *@file
 */

#ifndef PARAMETERS_H_INCLUDED 
#define PARAMETERS_H_INCLUDED 




//=================================
// included dependencies

#include <cstdlib>
#include <array>
#include <string>
#include <vector>
#include <ctime>

//=================================
// Code
// 
// This file defines parameters to be used in other parts of the 
// RSVS snaking framework. Default values are defined in "parameters.cpp"
//
// Substructure names are all 4-5 letters

/**
 * Namespace containing the parameter classes used to control execution of the
 * 3D-RSVS program..
 */
namespace param {
	/// Collects a lower and an upper bound.
	typedef std::array<double,2> realbounds;
	/// Collects the export settings which is a vector of pairs of strings.
	/// Each pair is: ["valid export type", "export config string"]
	typedef std::vector<std::pair<std::string,std::string>> exports;
	/// The input type of fill information.
	template <class T>
	struct filltype {
		bool active=false;
		T fill;
	};
	/**
	Parameters related to the Velocity calculation and VOS steps
	*/
	class rsvs
	{
	public:
		/**
		 * Algorithm used by Eigen to solve the SQP system.
		 * 
		 * See RSVScalc::SQPStep for details of the valid options.
		 */
		int solveralgorithm;
		/// Fill the VOS values with a constant value
		filltype<double> cstfill;
		/// Fill the VOS values from file filefill.fill
		filltype<std::string> filefill;
		/// Fill the VOS values from a run time function accessible from
		/// makefill.fill
		filltype<std::string> makefill;
		rsvs();
		~rsvs();
		void PrepareForUse();
	};

	/**
	Parameters controlling tuning parameters for the stepping of
	the restricted surface.
	*/
	class snaking
	{
	public:
		/// Distance along edge at which a vertex is considered
		/// arrived regardless of "d" and "v"
		double arrivaltolerance;
		/// Distance along edge at which converging snaxels are
		/// considered arrived
		double multiarrivaltolerance;
		/// maximum snake time step length
		double snaxtimestep;
		/// maximum snaxel distance movement
		double snaxdiststep;
		/// Initialisation boundary (either 0 or 1)
		int initboundary;
		/// maximum number of steps
		int maxsteps;

		snaking();
		~snaking();
		void PrepareForUse();
		
	};

	
	/// Parameters controlling cartesian grid properties
	class voxel
	{
	public:
		/// Size of the Background grid on which the VOS is defined.
		std::array<int, 3> gridsizebackground;
		/// Size of the Snaking grid on which the snake is defined.
		/// final size = gridsizebackground*gridsizesnake
		std::array<int, 3> gridsizesnake;

		voxel();
		~voxel();
		void PrepareForUse();
	};

	/// Class for handling of voronoi VOS meshing parameters.
	class voronoi
	{
	public:
		/// Vector of input points, 4 datums per point.
		std::vector<double> inputpoints;
		/// Distance at which to build the bounding box of the mesh.
		double distancebox;
		/// A string pointing to a file containing the set of inputpoints.
		std::string pointfile;
		/// The coarseneness level of the snaking mesh that will be generated.
		/// 1 -> same as VOS, 0.1 -> 1 tenth the edge length of the VOS.
		double snakecoarseness;

		voronoi();
		~voronoi();
		void PrepareForUse();
		void ReadPoints();
	};



	/// Class for parameters of the grid generation.
	class grid
	{
	public:
		voxel voxel;
		voronoi voronoi;
		/// Domain size in internal coordinates.
		std::array<realbounds, 3> domain;
		/// Physical domain size for export.
		std::array<realbounds, 3> physdomain;
		/// The type of grid to use either "voxel" or "voronoi" .
		std::string activegrid;
		
		grid();
		void PrepareForUse();
	};

	/**
	Class containing the input configuration these are files to load.
	*/
	class ioin{
	public:
		std::string snakemeshname;
		std::string volumeshname;
		std::string targetfill;
		std::string casename;

		ioin();
		void PrepareForUse();
	};
	/**
	Class containing the output configuration these are files to store and where
	to store them.
	
	This automatically parses output directory patterns to produce archival
	folders with time stamps and logical numbering.
	*/
	class ioout{
	private:

	public:
		std::string pathoutdir;
		std::string pathpattern;
		std::string basenamepattern;
		std::string basenameoutdir;
		std::string outdir;
		std::string pattern;

		bool redirectcout;
		bool redirectcerr;

		int logginglvl;
		int outputlvl;

		ioout();
		void PrepareForUse();
	};


	/**
	 * @brief      Class containing all parameter settings for file operations.
	 */
	class files
	{
	public:
		bool appcasename2outdir;
		ioin ioin;
		ioout ioout;
		exports exportconfig;
		files();
		void PrepareForUse();
	};


	/**
	 * @brief      Root class for all the parameters.
	 */
	class parameters
	{
	public:
		rsvs rsvs;
		snaking snak;
		grid grid;
		files files;

		void PrepareForUse();
	};

	/**
	 * Provide functions for reading and writing of the parameter structure.
	 */
	namespace io {
		void read(const std::string &fileName, parameters &p);
		void readflat(const std::string &fileName, parameters &p);
		void write(const std::string &fileName, const parameters &p);
		void writeflat(const std::string &fileName, const parameters &p);
		int updatefromstring(const std::vector<std::string> &flatjsonKeyVal,
			parameters &p, const std::string&& sep=std::string(":"));

		void defaultconf();
	}

	/**
	 * Tests for the parameter implementation.
	 */
	namespace test {

		int base();
		int io();
		int ioflat();
		int ipartialread();
		int prepareforuse();
		int autoflat();
		int symmetry();
	}
}



#endif