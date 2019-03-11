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


namespace param {
	typedef std::array<double,2> realbounds; 
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
		// Algorithm used by Eigen to solve the SQP system.
		int solveralgorithm;

		filltype<double> cstfill;
		filltype<std::string> filefill;
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
		// Distance along edge at which a vertex is considered
		// arrived regardless of "d" and "v"
		double arrivaltolerance;
		// Distance along edge at which converging snaxels are
		// considered arrived
		double multiarrivaltolerance;
		// maximum snake time step length
		double snaxtimestep;
		// maximum snaxel distance movement
		double snaxdiststep;
		// Initialisation boundary (either 0 or 1)
		int initboundary;
		// maximum number of steps
		int maxsteps;

		snaking();
		~snaking();
		void PrepareForUse();
		
	};

	/**
	Parameters controlling grid properties
	*/
	class voxel
	{
	public:
		std::array<int, 3> gridsizebackground;
		std::array<int, 3> gridsizesnake;

		voxel();
		~voxel();
		void PrepareForUse();
	};
	class voronoi
	{
	public:
		std::vector<double> inputpoints;
		double distancebox;
		std::string pointfile;

		voronoi();
		~voronoi();
		void PrepareForUse();
		void ReadPoints();
	};


	class grid
	{
	public:
		voxel voxel;
		voronoi voronoi;
		// Domain size
		std::array<realbounds, 3> domain;
		// Internal stretching of the meshes
		std::array<double, 3> stretch;
		std::string activegrid;
		
		grid();
		void PrepareForUse();
	};

	/**
	Class containing the input configuration these are files to load etc
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
	to store them
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
	class files
	{
	public:
		bool appcasename2outdir;
		ioin ioin;
		ioout ioout;
		files();
		void PrepareForUse();
	};

	class parameters
	{
	public:
		rsvs rsvs;
		snaking snak;
		grid grid;
		files files;

		void PrepareForUse();
	};


	namespace io {
		void read(const std::string &fileName, parameters &p);
		void readflat(const std::string &fileName, parameters &p);
		void write(const std::string &fileName, const parameters &p);
		void writeflat(const std::string &fileName, const parameters &p);
		int updatefromstring(const std::vector<std::string> &flatjsonKeyVal,
			parameters &p, const std::string&& sep=std::string(":"));

		void defaultconf();
	}

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