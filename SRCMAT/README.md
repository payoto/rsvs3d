# Matlab helper scripts

This folder contains an unfinished MATLAB interface, helper scripts,
and prototypes to test out what needed to be implemented in C++.

:Note: A Matlab license is not necessary to use the 3D-RSVS.

## Contents

### C code generation

The MATLAB symbolic toolbox was used to generate code for the matrix
heavy equations calculating the area and volumes of the 3D geometries.

The MATLAB file doing the code generation is `SRCMAT\PrepareCRSVSSource.m` it was used to generate the files: `SRCC\src\rsvs\RSVSmath_automatic.cpp`
and `SRCC\incl\RSVSmath_automatic.hpp`.

### Interface

A Matlab interface to the program is defined in `RSVS3D_interface.m` unfortunately it is undocumented and in matlab.

### Utilities

A number of files were copied as is from the [RSVS2D repository](https://github.com/payoto/rsvs2d). These are stored as is in the `SRCMAT/2DUTILS/` folder.

### Prototypes and helpers

Some of the data structures and output patterns were prototyped in MATLAB notably:

- `include_3DCheckGrid.m`: prototype the `SRCC/src/grid/mesh.cpp` operations for connectivity exploration.
