# README #
Last updated 13/12/2018

# What is this repository for? #

This repository is the C++ implementation of the 3D R-Snake Volume of Solid (RSVS) parameterisation.
It includes C++ main utility and Matlab support codes. The C++ code is here to do the heavy lifting,
the MATLAB code is here as it was used to prototype and test ideas. 

Relevant publications for the 2D RSVS are at the end of this readme.

The compiled binary is available for download for Windows 64bits and Linux 64bits.

# Pre-requisites #

For this code to work necessary programs:

 + MATLAB installed (2015a or later) - including parallel toolbox  
 + c++ compiler compatible with MATLAB for the compilation of mex files  
 + Standalone c++11 compiler for the compilation of console programs 
   (GCC/G++ v7.1 used for development)  
 + `make` to build the `RSVS3D` executable.  
 + fortran (90+) compiler for compilation of flow solvers   

Required 3rd party open source libraries for compilation:  

 + [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page): Library for linear algebra
   (templated, header only).
 + [boost/filesytem](https://www.boost.org/): Use some filesystem command for interface 
   (needs to be compiled).
 + [cxxopts](https://github.com/jarro2783/cxxopts/releases): Handling of command line arguments
   (header only).
 + [JSON for Modern C++](https://github.com/nlohmann/json/releases): JSON handling for c++. Used
   for the parameter handling of the RSVS3D framework (single include header).


# How do I get set up? #

## Matlab ##
Before being able to call the Matlab codes 
	
		>> Init3DMatlab
		>> Include_3DCheckGrid


## C++ ##
To Compile and test the C++ code:

		cd ./SRCC
		make testall      
		bin/testall_RSVS3D

If the tests run;

		make
		bin/RSVS3D.exe


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

 + `domain`: Domain dimensions, each of x, y and z are represented by a lower and upprt bound.
 + `gridsizebackground`: Design grid size on which the volume fractions are specified.
 + `gridsizesnake`: Snaking mesh as a refinement of the background mesh.

Examples:
	`gridsizebackground=[2, 3, 4]` and `gridsizesnake=[4, 4, 4]` will leed to an actual snaking mesh
	of `[8, 12, 16]`.

Parameters:
```json
	"/grid/voxel/domain/0/0": 0.0,
	"/grid/voxel/domain/0/1": 1.0,
	"/grid/voxel/domain/1/0": 0.0,
	"/grid/voxel/domain/1/1": 1.0,
	"/grid/voxel/domain/2/0": 0.0,
	"/grid/voxel/domain/2/1": 1.0,
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
# Contribution guidelines #

To add new features a few things must be considered (roughly in order of importance):

+ Maintaining the git archive's integrity and minimising merge conflicts for yourself and other users.
+ Execution by other users.
+ Simplicity of deployment.
+ Centralised control of execution flow, i.e. single parameter structure controls 
  all features of the execution and the case being run
+ Maintaining consistent coding style.

## Managing git ##

To manage the git repo's integrity new features should be developed in new folders 
contained in `./SRC<language of dvp>/` Make sure that names of functions are logical 
and express what they do. If large numbers of files are needed, make sure each contain
information about the creator, the feature they are part of and the intended purpose.

Use a separate git branch for development of new experimental features that may break
the main code (and not just the new feature).

Large new features should add their own parameter defaults.

## Execution by other users ##

Execution by other users relies on maintaining a consitent way of calling cases. 
This means that parameters for new features should use the parameter structure syntax already
in place.

## Simple deployment ##

Provide scripts to deploy any new utilities required and make sure they are executed by the `deploylinux.sh` code.
Ideally compilation of C++ code should be done from `SRCC/makefile`  

## Coding Style ##

Lines of code shall not be longer than 80 characters.

camelCase is used:
```Matlab
	MyNewFunctionName 
	myNewVariableName 
	mynewstructure.withitsfield
	% Unfortunately I was not always super consistent for structures you might encounter:
	myNewStructure.withitsfield
	myNewStructure.withItsField
	% Curse past-me silently and get on with your life.
```
The aproach followed in the code is close to functional programming. The idea is to have short functions (20 lines)
which perform a single simple process. The functions should be named according to their purpose and should include 
just below a short (1 sentence) description of what they are supposed to do. Inputs should have clear meaningful 
names and a brief description of what they are should be in the function if necessary.

Objects may be used where they make sense, i.e. to protect the integrity of co-dependant data structures. The example
is meshes that need to be edited: rather than deleting elements "manually" use class functions that ensure connectivity
is maintained safely. Encapsulation of these features is a work in progress.

The code is geared to have robust predictable behaviour rather than speed. It is built to be modular and facilitate the
future implementation of unintended features. Any new feature should be easy to enable from file `StructOptimParam.m`.

All parameters are added to a single structure at the start in `StructOptimParam.m`. These can then be accessed using:
```Matlab
[var1,var2,var3,...,varN]=ExtractVariables({cell array of variable names},paramstructure);
paramstructure=SetVariables({cell array of variable names},{var1,var2,var3,...,varN},paramstructure);
```
# Gosh, how the hell does this parameter thing work? #

## Motivation ##

It seems a bit complicated but the idea is to reduce the number of inputs and outputs functions have and make sure all
functions have all the inputs they need even when they are changed without having to edit an entire stack of upstream
functions. At the most basic level setting parameters has been centralised: it all happens in `StructOptimParam` and
`structInputVar`.

All access to parameters, to extract and more rarely to change the value of a parameter, is done through two functions `ExtractVariables` and `SetVariables`.
An example call is shown below
```Matlab
[var1,var2,var3,...,varN]=ExtractVariables({cell array of variable names},paramstructure);
paramstructure=SetVariables({cell array of variable names},{var1,var2,var3,...,varN},paramstructure);
```
The structure `paramoptim` (or `paroptim` oops I've been inconsistant) is then passed around. This greatly reduces 
the ease with which new features can be tested and implemented. And no globals are needed for passing parameters. 
And it scales really well.

## Initialising the parameter structures ##

All parameters are stored and passed using a single structure: `paramoptim`  
This structure has got fields which reflect the purpose of the parameters:  
`paramoptim.obj.flow` for the aerodynamic objective functions  
`paramoptim.constraint` for optimisation constraints  

All these parameters are set in file `StructOptimParam.m`. All of these fields is set by calling local function:
`DefaultOptim`. This function then calls a set of subfunctions assigning each of the results to a single field. 
The calls for the two examples above are shown here.

	paramoptim.obj.flow=DefaultCutCell_Flow();
    paroptim.constraint=DefaultConstraint();

This populates paramoptim with a set of default values which can then be changed in additional local functions added to `StructOptimParam.m`.

One field is an exception:
`paramoptim.parameterisation` is for parameters needed for the RSVS and is built independently using file `structInputVar.m`.
`structInputVar.m` works exactly like `StructOptimParam.m` but is for any parameters needed in function `Snakes.m` (the RSVS engine).


## Creating a parameter set to run a case ##

To create a new run case a local function must be added at the end of `StructOptimParam.m`. This process is shown for a
function called `NewCaseOfUserX`. This set of parameters once added can be called from the MATLAB console using
`ExecuteOptimisation('NewCaseOfUserX')`. 

`StructOptimParam` works by calling `[paroptim]=eval(caseStr);` where `caseStr` is the first input to `ExecuteOptimisation`.

### Simple case ###

There are three steps to building a new parameter set to that will run a new optimisation case:  
 1.  Building the entire parameter structure (using a callable function with a single output)  
 2.  Aplying standard modifications (using a callable function with a valid structure as an input and output)  
 3.  Applying bespoke parameter changes (using . structure assignements or `SetVariables method`).  

Below is an example of this process:

	function [paroptim]=NewCaseOfUserX()
		% 1. Build parameter structure with default settings
		[paroptim]=DefaultOptim(); 
		
		% 2. Apply a Standard modification to default settings: for example change all the 
		% settings required to change the optimiser for aero, differential evolution cases:
		paroptim=CutCellObjective(paroptim);
		paroptim=OptimDE(paroptim);
		
		% 3. Single changes to specific parameters:
		paroptim.general.maxIter=5; % maximum number of optimisation iteration
		paramoptim=SetVariables({'worker'},{2},paramoptim); % number of parallel workers
		% Previous line equivalent to: paroptim.general.worker=2;
	end


### I liked my previous case and want to keep it but want to change one parameter... I'll just copy paste it, right? ###

No. Don't do that. Part of the reason the current code looks terrifying is because I did this 2 years ago. We never recovered.
Since then I've found better solution(s). Just call the parent case `NewCaseOfUserX` instead of `DefaultOptim()`. Your `NewCaseOfUserX_v2`
will look something like this:

	function [paroptim]=NewCaseOfUserX_v2()
		% 1. Build parameter structure with default settings
		[paroptim]=NewCaseOfUserX(); 
		
		
		% 3. Single changes to specific parameters:
		paroptim.general.maxIter=30; % maximum number of optimisation iteration
	end

`NewCaseOfUserX` is still callable and if you make any changes to it they will be duplicated to case `NewCaseOfUserX_v2`. This is 
very useful when you want to make sure a bunch of cases are similar.
	
### Building your own standard modifications ###

You've run your first case and you liked a lot of the settings? The new feature you implemented requires the modification 
of the same 10 default parameters every time? You've used a parent case, but now you have children of children of 
children and backtracking through that stack is a pain?
 
This is the time to add a new "standard" modification function. These functions do not define a callable execution case
but modify a few of the input parameter directly in the structure and return the modified structure. For this let us look 
at function `CutCellObjective` which is used whenever the cut-cell flow solver needs to be run.

	function paroptim=CutCellObjective(paroptim)
		
		paroptim.general.objectiveName='CutCellFlow';
		paroptim.general.direction='min';
		paroptim.general.defaultVal=1e3;
	end

### Advanced - Making parameter sweeps ###

These bespoke case functions can be adapted to accept inputs.

	function [paroptim]=NewCaseOfUserX_varin(n)
		% 1. Build parameter structure with default settings
		[paroptim]=NewCaseOfUserX(); 
		
		
		% 3. Single changes to specific parameters:
		paroptim.general.maxIter=n; % maximum number of optimisation iteration
	end

To call this function from the console is done with the following command:

	ExecuteOptimisation('NewCaseOfUserX_varin(5)')

this will set the value of `maxIter` to be 5.

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