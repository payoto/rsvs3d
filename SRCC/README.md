# C++ code structure

The `SRCC` folder contains the C++ source files used to build the executable
for the 3D-RSVS software.

## Contents

The files and folders inside the `SRCC/` folder are:

### Source directories

- `incl`: Contains the "includes" of the program, this encompasses header files
  (`hpp`) and template implementations (`cpp` files that need to be included).
- `src`: Contains the implementation of the software.
- `modules`: Contains git sub-modules used for managing build time dependencies
  of the software.

### Makefiles

- `makefile`: The main makefile defining compilation recipes, important ones
  are:
  - `all`: Build the complete executable with main, and test suite.
  - `distribution`: Builds the executable without any test code and with
    optimisation level of 3.

The main makefile calls on either the default linux or the default windows
configuration depending on the OS:

- `makelinux.mk`: Default configuration for building the program in linux
  (tested on Ubuntu 20.04).
- `makewin.mk`: Default configuration for the building the program on
  windows (tested on windows 10).

### Supporting directories

- `docs`: Doxygen configuration file for building the API docs for the program.
  See the [docs/README](docs/README.md) for instructions on compiling the documentation.
  The compiled documentation can be found at:
  [payoto.github.io/rsvs3d](https://payoto.github.io/rsvs3d/).
- `config`: Template input files for the RSVS process, each of the JSON
  files can be called by doing ``

### Automatically generated directories

- `bin`: The compiled binaries go there.
- `.d`: Dependency management artifact to ensure make rebuilds objects when
  headers change.
- `obj`: Directory for the compiled objects (one per `.cpp` file in `src`).

## Building

In order to build the fully featured executable run:

```bash
make all
```

This will generate the `RSVS3D` executable. When developing a specific part of
the program it can be desirable to build a smaller executable which only compiles
parts of the program being tested. The `testnew` and `testall` recipes in the makefile
fill that role:

```bash
make testnew
```

will create the `test_RSVS3D` executable with debug flags.

## Usage

When the `RSVS3D` executable has been built the `--help` flag can be used to return
the command line flags.

```bash
$ ./RSVS3D --help
Program for the execution of the Restricted-Surface Volume of Solid in 3D
Usage:
  ./RSVS3D [OPTION...]

  -h, --help  Print help

 Execution control options:
  -n, --noexec [=STRING(=./noexec_config.json)]
                                Do not execute RSVS process, will only parse
                                the inputs and output the resulting
                                configuration file to 'arg'
  -e, --exec                    Execute RSVS. With no command line argument
                                the program does nothing.
      --test [=STRING(=short)]  Executes specified tests. requires
                                compilation without flag RSVS_NOTESTS.

 Parameter configuration options:
  -u, --use-config STRING       Use one of the predefined configurations
                                stored in the code.
  -l, --load-config FILES       Load configuration file in JSON format to set
                                parameter structure. Multiple files can be
                                specified and will be processed in order of
                                appearance.
  -p, --param KEY:VAL           Define a parameter manually on the command
                                line. The format must be a flat key into the
                                JSON configuration: (e.g.: '/snak/maxsteps:50'
                                will set 'param.snak.maxsteps=50')
      --default-config [=FILE(=default_config)]
                                Output the default configuration to a file.
```

## Testing

Tests are handled by a custom testing framework developed for this project:

- the testing framework is defined in `incl/test.hpp`;
- calling the framework is done in `src/test/test.cpp`;
- the actual test functions are defined in the `cpp` files they correspond to
  and all start with `Test_`.

The test suite can be executed by building the program with `make all` and
then calling `./RSVS3D --test` which runs a subset of the tests. Running `./RSVS3D --test=all` will run all the tests, which is very slow, and fails (at the moment).
