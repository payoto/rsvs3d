# README #
Last updated 27/02/2017 (for the 2D Version)

# What is this repository for? #

This repository is the C++ implementation of the 3D R-Snake Volume of Solid (RSVS) parameterisation. It includes C++, Matlab and Python code. The C++ code is here to do the heavy lifting, the MATLAB code is here as it was used to prototype and test ideas. Eventually python will be used to build interfaces to use the C++ code.

Relevant publications for the 2D RSVS are at the end of this readme.


# Pre-requisites #

For this code to work necessary programs:   
 +  MATLAB installed (2015a or later) - including parallel toolbox  
 +  c compiler compatible with MATLAB for the compilation of mex files  
 +  Standalone c compiler for the compilation of console programs (GCC/G++ v7.1 used for development) 
 +  `make` for linux and `nmake` on windows should
 +  fortran (90+) compiler for compilation of flow solvers   



# How do I get set up? #

## Matlab ##
To call the Matlab code 
	
		>> Init3DMatlab
		>> Include_3DCheckGrid


## C++ ##
To Compile and test the C++ code:

		cd ./SRCC

		nmake test -a  (windows)
		make test      (linux)

		main_test

If the tests run;

		nmake -a  (windows)
		make      (linux)
		main

## 2D Instructions not yet implemented ##
(Any - At the start of a new MATLAB session) 
	in matlab:
	
		>> InitialiseSnakeFlow
		>> ExecInclude
		
Tests to run:
	Main('SnakesFoilVVSmall4')
	% Tests the Snaking process on a multi-topology case
	ExecuteOptimisation('Test_Rosenbrock')
	% Tests the optimisation framework on an analytical function
	
These tests will save results in:

	../results/Standard_Execution/<Archive_YYYY_MM>/Day_YY_MM_DD/
	
and
	
	../results/Optimisation/<Archive_YYYY_MM>/Day_YY_MM_DD/
	
respectively. All result files and folders are time stamped such that a named sort will return the files in chronological order of creation.
These files and their location are entered into a file  at:
	../results/<archive_name>/Index_<archive_name>.txt

# Note on using git to update a private copy of the code.#

Updating your files to be up to date with the master branch can be done using git very efficiently. With a few steps.  
 +  Add all your local changes `git add -u` then `git add *.m` then `git add Active_Build*png`   
 +  Commit all your local changes `git commit  -m "Add comment about what was done"`
 +  Switch to the master branch `git checkout master`  
 +  Pull the latest version from the remote repository: `git pull`  
 +  If there are any merge issues resolve them using a text editor (if there are you will need to run `git add -u` and `git commit` before the next step)   
 +  Switch to your local branch `git checkout <your branch name>`   
 +  Merge the new master with your local branch `git merge master`   
 +  If there are any merge issues resolve them using a text editor (if there are you will need to run `git add -u` and `git commit`)


# Contribution guidelines #

To add new features a few things must be considered (roughly in order of importance):
+ Maintaining the git archive's integrity and minimising merge conflicts for yourself and other users.
+ Execution by other users.
+ Simplicity of deployment.
+ Centralised control of execution flow, i.e. single parameter structure controls all features of the execution and the case being run
+ Maintaining consistent coding style.

## Managing git ##

To manage the git repo's integrity new features should be developed in new folders contained in `./SRC<language of dvp>/`
Make sure that names of functions are logical and express what they do. If large numbers of files are needed, make sure each contain information
about the creator, the feature they are part of and the intended purpose.

Use a separate git branch for development of new experimental features that may break the main code (and not just the new feature).

Large new features should add their own parameter defaults.

## Execution by other users ##

Execution by other users relies on maintaining a consitent way of calling cases. This means that parameters for new features
should use the parameter structure syntax already in place.

## Simple deployment ##

Provide scripts to deploy any new utilities required and make sure they are executed by the `deploylinux.sh` code.
Ideally compilation of C++ code should be done from `SRCC/makefile`  

## Coding Style ##

Lines of code shall not be longer than 80 characters.

camelCase is used:

	MyNewFunctionName 
	myNewVariableName 
	mynewstructure.withitsfield
	% Unfortunately I was not always super consistent for structures you might encounter:
	myNewStructure.withitsfield
	% Curse past-me silently and get on with your life.
	
The aproach followed in the code is close to functional programming. The idea is to have short functions (20 lines) 
which perform a single simple process. The functions should be named according to their purpose and should include just below a short (1 sentence) description of what they are supposed to do. Inputs should have clear meaningful names and a brief description of what they are should be in the function if necessary.

Objects may be used where they make sense, i.e. to protect the integrity of co-dependant data structures. The example is meshes that need to be edited: rather than deleting elements "manually" use class functions that ensure connectivity is maintained safely. Encapsulation of these features is a work in progress.

The code is geared to have robust predictable behaviour rather than speed. It is built to be modular and facilitate the future implementation of unintended features. Any new feature should be easy to enable from file `StructOptimParam.m`.

All parameters are added to a single structure at the start in `StructOptimParam.m`. These can then be accessed using:

	[var1,var2,var3,...,varN]=ExtractVariables({cell array of variable names},paramstructure);
	paramstructure=SetVariables({cell array of variable names},{var1,var2,var3,...,varN},paramstructure);

# Gosh, how the hell does this parameter thing work? #

## Motivation ##

It seems a bit complicated but the idea is to reduce the number of inputs and outputs functions have and make sure all functions 
have all the inputs they need even when they are changed without having to edit an entire stack of upstream functions. At the most 
basic level setting parameters has been centralised: it all happens in `StructOptimParam` and `structInputVar`.

All access to parameters, to extract and more rarely to change the value of a parameter, is done through two functions `ExtractVariables` and `SetVariables`.
An example call is shown below

	[var1,var2,var3,...,varN]=ExtractVariables({cell array of variable names},paramstructure);
	paramstructure=SetVariables({cell array of variable names},{var1,var2,var3,...,varN},paramstructure);

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

# I don't get it what does this ACTUALLY do and who do I talk to?#

For more information about what the code does (i.e. the science of it)  
[Restricted Snakes: a Flexible Topology Parameterisation Method for Aerodynamic Optimisation](https://arc.aiaa.org/doi/pdf/10.2514/6.2017-1410)  
[Mixing and Refinement of Design Variables for Geometry and Topology Optimization in Aerodynamics](https://arc.aiaa.org/doi/pdfplus/10.2514/6.2017-3577)  

(Also available on research gate)

Alexandre Payot - a.payot@bristol.ac.uk  
[ResearchGate profile](https://www.researchgate.net/profile/Alexandre_Payot)  
[Google Scholar profile](https://scholar.google.co.uk/citations?user=JX_AmkwAAAAJ&hl=en)  