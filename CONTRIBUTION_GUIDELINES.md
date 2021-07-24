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
