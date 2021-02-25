# 3D-RSVS #

[![Linux build](https://github.com/payoto/rsvs3d/actions/workflows/linux.yml/badge.svg)](https://github.com/payoto/rsvs3d/actions/workflows/linux.yml)
[![Windows build](https://github.com/payoto/rsvs3d/actions/workflows/windows.yml/badge.svg)](https://github.com/payoto/rsvs3d/actions/workflows/windows.yml)

[![3D Restricted Snakes Volume of Solid Parameterisation](https://raw.githubusercontent.com/payoto/rsvs3d/master/SRCC/docs/3DRSVS_LR1.gif)](https://github.com/payoto/rsvs3d/releases)

The 3D-RSVS is a geometry generation tool using volume specification to build
smooth surfaces. Papers describing the method are available on 
[researchgate](https://www.researchgate.net/publication/330197375_Parametric_Surfaces_with_Volume_of_Solid_Control_for_Optimisation_of_Three_Dimensional_Aerodynamic_Topologies) 
and [ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0045793018304730).


Last updated 12/11/2019
[Full documentation](https://payoto.github.io/rsvs3d_docs/)
available.

# What is this repository for? #

This repository is the C++ implementation of the 3D R-Snake Volume of Solid (RSVS) parameterisation.
It includes C++ main utility and Matlab support codes. The C++ code is here to do the heavy lifting,
the MATLAB code is here as it was used to prototype and test ideas. 

Relevant publications for the 2D RSVS are at the end of this readme.

The compiled binary is available for download for Windows 64bits and Linux 64bits.

# Install notes #

## Get the code: ##

This project uses git submodules to handle dependencies:

	git clone --recurse-submodules https://github.com/payoto/rsvs3d

or

	git clone https://github.com/payoto/rsvs3d
	git submodule update --init --recursive


## Required tools for compilation: ##

For this code to work necessary programs:

 + Standalone c++11 compiler for the compilation of console programs 
   (GCC/G++ v7.1 used for development)
 + `make` to build the `RSVS3D` executable.

## Wrappers: ##

Execution wrappers:

 + MATLAB installed (2016a or later)  
 + c++ compiler compatible with MATLAB for the compilation of mex files  
 + fortran (90+) compiler for compilation of flow solvers   

## External libraries ##

Required 3rd party open source libraries for compilation 
(dependencies available at [`github/payoto/rsvs3d-dependencies`](https://github.com/payoto/rsvs3d-externals)):

 + [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page): Library for linear algebra
   (templated, header only).
 + [boost/filesytem](https://www.boost.org/): Use some filesystem command for interface 
   (needs to be compiled).
 + [cxxopts](https://github.com/jarro2783/cxxopts/releases): Handling of command line arguments
   (header only).
 + [JSON for Modern C++](https://github.com/nlohmann/json/releases): JSON handling for c++. Used
   for the parameter handling of the RSVS3D framework (single include header).

Optional 3rd party open source library:
 + [Tetgen](http://tetgen.org) : A Quality Tetrahedral Mesh Generator and a 3D Delaunay
 Triangulator. Download my modified version for this project 
 [payoto/tetgen](https://github.com/payoto/tetgen).

All dependencies are intended to be handled via git submodules, to install them using submodules, run `git submodule update --init --recursive`.

# License #

Any code using the Tetgen interface and functionalities is under the GNU Affero GPL license
which is more restrictive than the LGPL of this project:
 + LGPL : If you write an interface to this project and use it as a library, you can distribute
 closed source version.
 + AGPL : Regardless how you use it and distribute a program you have to open the source.

Refer to the full License terms for more information.

# How do I get set up? #

## Matlab ##
Before being able to call the Matlab codes (these are not necessary for the geometry tool
execution).
	
		>> Init3DMatlab
		>> Include_3DCheckGrid


## C++ ##
To Compile and test the C++ code:

		cd ./SRCC
		make testall      
		testall_RSVS3D

If the tests run (they should), to see some example usages;

		make
		RSVS3D.exe -h 


## Note on using git to update a private copy of the code. ##

Begginning to use git? Follow these [5mn ELI5 explainer](https://dev.to/sublimegeek/git-staging-area-explained-like-im-five-1anh) which 
will help understand what the lingo means, [git begginer guide](https://www.atlassian.com/git/tutorials) which should get you up and
running, or the full git documentation if you're trying to do something [git documentation](https://git-scm.com/doc).

Very minimal guide I wrote a while back specific to the repository:
Updating your files to be up to date with the master branch can be done using git very efficiently. With a few steps.  

 +  Add all your local changes `git add -u` then `git add *.m` then `git add Active_Build*png`   
 +  Commit all your local changes `git commit  -m "Add comment about what was done"`
 +  Switch to the master branch `git checkout master`  
 +  Pull the latest version from the remote repository: `git pull`  
 +  If there are any merge issues resolve them using a text editor (if 
 	there are you will need to run `git add -u` and `git commit` before the next step)   
 +  Switch to your local branch `git checkout <your branch name>`   
 +  Merge the new master with your local branch `git merge master`   
 +  If there are any merge issues resolve them using a text editor (if there are you 
 	will need to run `git add -u` and `git commit`)

## Getting help ##

Use the issues board.

# Using the 3D-RSVS #

Generating a geoemtry using the 3D-RSVS method only requires the executable `RSVS3D.exe`.
For basic usage information from the command line use:  

	RSVS3D --help

Warning:  
  Running `RSVS3D` with no command line arguments does nothing.

## Command line options ##

Below are all the possible commad line options for the RSVS3D program. These can be assembled 
in arbitrary ways to run a specific config. The long name is shown (called with prefixed `--`
on the command line) followed by 

 + help (`-h`): Display command line help; 
 + noexec (`-n`): Do not run the RSVS process and output the configuration file; 
 + exec (`-e`): Execute the RSVS3D for the default case;
 + use-config (`-u`): Use system configuration `STRING` (none specified yet);
 + load-config (`-l`): Load a configuration file from `FILE`;
 + param (`-p`): Overwrite a specific parameter specified by `KEY:VAL`. "key" is the name of that 
 paramaeter as it appears in the flattened JSON parameter files, "val" is the value of that
 patameter;
 + default-config: Outputs a configuration file with all the default value assigned to the
 parameters.
 + test: Runs the specified test suite one of `new`, `short` or `all`. May be disabled in releases.

## Parameter control ##

Internally parameters are controlled by a single structure defined in [parameters.h](SRCC/incl/parameters.h).
Externally parameters are handled using [_JSON_](https://www.json.org/) files. These provide a good balance 
of human and machine readable format. And support intricate tree structures and nesting. The _JSON_ interaction
is handled by an external library [JSON for Modern C++](https://github.com/nlohmann/json/releases). This library
allows two types of _JSON_ files: normal and flat. Default parameter configuration files showing all the parameters
and there default options in [default_config](SRCC/config/default_conf.json) and 
[default_configflat](SRCC/config/default_confflat.json). Below are two _JSON_ examples.

Example normal _JSON_:
```json
{
  "files": {
    "appcasename2outdir": true,
    "ioin": {
      "casename": "",
      "snakemeshname": "",
      "targetfill": "",
      "volumeshname": ""
    },
  },
  "grid": {
  	"voxel": {
  		"gridsizebackground": [1,1,1],
  	},
  },
}
```

Equivalent flat _JSON_:

```json
{
  "/files/appcasename2outdir": true,
  "/files/ioin/casename": "",
  "/files/ioin/snakemeshname": "",
  "/files/ioin/targetfill": "",
  "/files/ioin/volumeshname": "",
  "/grid/voxel/gridsizebackground/0": 1,
  "/grid/voxel/gridsizebackground/1": 1,
  "/grid/voxel/gridsizebackground/2": 1,
}
```

Three command line options are currently available for parameter control:
the _use-config_, _load-config_ and _param_ give control over the execution of the RSVS. It permits
the control of execution flow, output level and output location as well as specific mesh and volume
information.
The _use-config_, _load-config_ and _param_ options can be combined to get the desired set of parameters,
multiple ones of each can be called. The program throws an error if a parameter is not recognised
or not correctly read from a file. These 3 types of options are parsed in order of their appearance
in the help and this (readme):  _use-config_ then _load-config_ and finally _param_. Inputs of the same
type are then parsed in their order of appearance. 

_Load-config_ can load an incomplete set of parameters overwriting only parameters that are specified.

## Non exhaustive parameter list ##

For up to date parameter list check the default [configuration files](SRCC/config).

### files ###

Controls the file interaction of the program including the naming of output folders.

```json
	"/files/appcasename2outdir": true,
	"/files/ioin/casename": "",
	"/files/ioin/snakemeshname": "",
	"/files/ioin/targetfill": "",
	"/files/ioin/volumeshname": "",
	"/files/ioout/basenameoutdir": "rsvs3d_",
	"/files/ioout/basenamepattern": "%y%m%dT%H%M%S_",
	"/files/ioout/logginglvl": 2,
	"/files/ioout/outdir": "",
	"/files/ioout/outputlvl": 2,
	"/files/ioout/pathoutdir": "../out",
	"/files/ioout/pathpattern": "Archive_%Y_%m/Day_%y-%m-%d",
	"/files/ioout/pattern": "",
	"/files/ioout/redirectcerr": false,
	"/files/ioout/redirectcout": false,
```
#### ioin ####

 + `appcasename2outdir`: append casename to output dir path?
 + `casename`: Name of the case.
 + `snakemeshname`: Mesh file to load.
 + `targetfill`: Unused (see [rsvs](#rsvs))

#### ioout ####

 + `basenameoutdir`: Name of the output directory.
 + `basenamepattern`: time format string added to the basenameoutdir.
 + `logginglvl`: Depth of data logging 0-minimal, 1-Logs only, 2-Snake history, 3-All data.
 + `outdir`: Leave empty to use the automatic archive directory 
 trees, otherwise the output directory.
 + `outputlvl`: Depth of final data output.
 + `pathoutdir`: Root directory (relative or absolute) for the archiving tree.
 + `pathpattern`: Directory stub to use as a time format which will be assembled to generate
 an archiving output folder pattern.
 + `pattern`: Used internally to store the pattern generated by `basenamepattern`.
 + `redirectcerr`: redirection of standard error to a file.
 + `redirectcout`: redirection of standard output to a file.

### grid ###

Control the underlying grid if it is generated. It can also be loaded if `"/files/ioin/snakemeshname"`
is specified.

 + `activegrid` : The type of grid to build (`"voxel"`, `"voronoi"` or `"load"`).
 + `domain`: Domain dimensions, each of x, y and z are represented by a lower and upprt bound.
 + `gridsizebackground`: Design grid size on which the volume fractions are specified.
 + `gridsizesnake`: Snaking mesh as a refinement of the background mesh.
 + `distancebox` : for a voronoi VOS mesh the distance outside `domain` at which the
 bounding points will be placed.
 + `inputpoints` : A vector of data containing coordinates used for the Voronoi process.
 + `pointfile` : The file from which these are loaded.

Examples:
	`gridsizebackground=[2, 3, 4]` and `gridsizesnake=[4, 4, 4]` will leed to an actual snaking mesh
	of `[8, 12, 16]`.

Parameters:
```json
    "/grid/activegrid": "voxel",
	"/grid/domain/0/0": 0.0,
	"/grid/domain/0/1": 1.0,
	"/grid/domain/1/0": 0.0,
	"/grid/domain/1/1": 1.0,
	"/grid/domain/2/0": 0.0,
	"/grid/domain/2/1": 1.0,
  	"/grid/stretch/0": 1.0,
 	"/grid/stretch/1": 1.0,
 	"/grid/stretch/2": 1.0,
    "/grid/voronoi/distancebox": 0.1,
    "/grid/voronoi/inputpoints/0": 0.0,
    "/grid/voronoi/pointfile": "",
	"/grid/voxel/gridsizebackground/0": 1,
	"/grid/voxel/gridsizebackground/1": 1,
	"/grid/voxel/gridsizebackground/2": 1,
	"/grid/voxel/gridsizesnake/0": 6,
	"/grid/voxel/gridsizesnake/1": 6,
	"/grid/voxel/gridsizesnake/2": 6,
```
### rsvs ###

RSVS process control. Includes the selection of which volume fraction the 3D-RSVS needs
to match.

 + `cstfill`: constant fill in all the volume cells.
 + `filefill`: Fill is specified in a file (space delimited data).
 + `makefill`: Programmaticaly defined fill information.
 + `solveralgorithm`: Chooses the solution process for the Quadratic Problem of the RSVS.

Only one of `filefill`, `makefill` or `cstfill` is taken into account if they are all set to _active_. 
The order of precendence is:

 1. `filefill`
 2. `makefill`
 3. `cstfill`

Parameters:
```json
	"/rsvs/cstfill/active": false,
	"/rsvs/cstfill/fill": 0.5,
	"/rsvs/filefill/active": false,
	"/rsvs/filefill/fill": "",
	"/rsvs/makefill/active": true,
	"/rsvs/makefill/fill": "",
	"/rsvs/solveralgorithm": 0,
```
### snak ###

Control the restricted snaking process, can have a large impact on the speed and quality 
of the convergence of the RSVS process.

 + `arrivaltolerance`: Distance from a vertex at which a snaxel is considered "arrived".
 + `initboundary`: Initialisation boundary (1 or 0). Which volume fraction boundary should 
 the surface be started at. 
 + `maxsteps`: Maximum number of snake steps.
 + `multiarrivaltolerance`: When two snaxels converge on a vertex what is the radius at which
 the arrival procedure is triggered.
 + `snaxdiststep`: Maximum non-dimensional distance which can be covered by a snaxel in 1 step.
 + `snaxtimestep`: Maximum time step (used for damping of the SQP).

Parameters
```json
	"/snak/arrivaltolerance": 1e-07,
	"/snak/initboundary": 1,
	"/snak/maxsteps": 50,
	"/snak/multiarrivaltolerance": 0.01,
	"/snak/snaxdiststep": 0.9,
	"/snak/snaxtimestep": 0.9
```

# I don't get it what does this ACTUALLY do and who do I talk to? #

For more information about what the code does (i.e. the science of it)  
[Restricted Snakes: a Flexible Topology Parameterisation Method for Aerodynamic Optimisation](https://arc.aiaa.org/doi/pdf/10.2514/6.2017-1410)  
[Mixing and Refinement of Design Variables for Geometry and Topology Optimization in Aerodynamics](https://arc.aiaa.org/doi/pdfplus/10.2514/6.2017-3577)  

(Also available on research gate)

Alexandre Payot - a.payot@bristol.ac.uk  
[ResearchGate profile](https://www.researchgate.net/profile/Alexandre_Payot)  
[Google Scholar profile](https://scholar.google.co.uk/citations?user=JX_AmkwAAAAJ&hl=en)  
[personal GitHub/payoto](https://www.github.com/payoto)  
[Research group GitHub/farg-bristol](https://www.github.com/farg-bristol)  
