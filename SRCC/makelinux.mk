
#------------------------------------------------------
#                 TARGET COMPILER
#------------------------------------------------------
# Specify compiler (only g++ is supported)
CC=g++-9
# Specify linker
LINK=g++-9

#------------------------------------------------------
#                 DIRECTORIES and INCLUDES
#------------------------------------------------------
# Executable name
MAINEXECUTABLE=RSVS3D
TESTEXECUTABLE=test_RSVS3D
TESTSHORTEXECUTABLE=testshort_RSVS3D
TESTALLEXECUTABLE=testall_RSVS3D

# source, object and dependency directories
OBJDIR = obj
SRCTREE = src
DEPDIR := .d
BINDIR = bin
INCLPROJ= incl

# Include directories
SHAREINCLROOTDIR = ./modules/rsvs3d-externals
DIREIGEN= $(SHAREINCLROOTDIR)/Eigen/include
DIRJSON= $(SHAREINCLROOTDIR)/json/include
DIRCXXOPTS= $(SHAREINCLROOTDIR)/cxxopts-2.1.1/include
DIRTETGEN= modules/tetgen

DIRPOLYSCOPE= modules/polyscope
INCLPOLYSCOPE= $(DIRPOLYSCOPE)/include $(DIRPOLYSCOPE)/deps/glm
INCLPOLYSCOPE+=$(DIRPOLYSCOPE)/deps/imgui/imgui
LDPOLYSCOPE=$(DIRPOLYSCOPE)/build/src
LDPOLYSCOPEDEPS= $(DIRPOLYSCOPE)/build/deps/stb
LDPOLYSCOPEDEPS+= $(DIRPOLYSCOPE)/build/deps/glad/src
LDPOLYSCOPEDEPS+= $(DIRPOLYSCOPE)/build/deps/glfw/src
LDPOLYSCOPEDEPS+= $(DIRPOLYSCOPE)/build/deps/imgui

ifeq ($(HEADLESS), true)
# glad glfw3 may not be available in headless mode you may also need to
# LIBPOLYSCOPE= polyscope stb imgui pthread
LIBPOLYSCOPE= polyscope stb glad glfw3 imgui dl X11 pthread
else
LIBPOLYSCOPE= polyscope stb glad glfw3 imgui dl X11 pthread
endif
#------------------------------------------------------
#           WARNING AND OPTIMISATION FLAGS
#------------------------------------------------------
# Standard GCC warning, optimisation and c++ standard flags
# WARNFLAGS= -Wall -Wextra -pedantic
WARNFLAGS= -Wall -Wextra -fpermissive
# -Werror
#
OPTIMLVL=3
CCPSTD=17
DEBUG=false

#------------------------------------------------------
#                 USE OF OPTIONAL LIBRARIES
#------------------------------------------------------
# Some libraries are optional in certain conditions
# Specify by "true" or "false" whether they should be used

# Boost is needed for the <boost/filesystem> header
# except under std=c++17 which added a standard filesystem header
# (Broken on GCC v8.1 + MinGW, must still use boost)
USE_BOOST = false
BOOST_VERSION_NUMBER = 1_68
BOOST_BUILD_TOOLCHAIN = mgw63-mt-sd-x64
DIRBOOST= $(SHAREINCLROOTDIR)/boost_$(BOOST_VERSION_NUMBER)/include
# BOOST can be ignored if using c++17 (not on GCC v8.1 on Win (25/01/19))
# boost compiled library file location (for linker)
LDDIRBOOST = $(SHAREINCLROOTDIR)/boost_$(BOOST_VERSION_NUMBER)/lib
# boost libraries to load: boost_<library>
LDFILEBOOST = boost_filesystem boost_system
# boost version appropriate to the system
BOOSTVERSION = -$(BOOST_BUILD_TOOLCHAIN)-$(BOOST_VERSION_NUMBER)

# intel MKL libraries for faster Eigen maths
# Not compatible with GCC on windows
DIRMKL =
# -isystem$(DIRMKL)
INCL=

SRCTREE+= $(DIRTETGEN)
INCLDIR= $(INCLPROJ) $(DIRJSON) $(DIRCXXOPTS) $(DIRTETGEN) $(INCLPOLYSCOPE)
INCLDIRSYS = $(DIREIGEN) $(DIRBOOST)



#------------------------------------------------------
#                 COMPILER FLAGS
#------------------------------------------------------

# Compilation flags
# gdb (debug) flags:
#  -g -ggdb
CCFLAGDBG=

# Flags defined by the RSVS project:
# -DDEBUGLVL1 -DSAFE_ACCESS  -DSAFE_ALGO -DTIME_EXEC
# -DRSVS_DIAGNOSTIC -DRSVS_VERBOSE -DRSVS_DIAGNOSTIC_RESOLVED
# -DRSVS_ACCESS_DEVELOPMENT_PARAMETERS -DTIME_EXEC
# -DRSVS_NO_ARGCHECK
CCFLAGPROJ= -DRSVS_DEBUG -DSAFE_ACCESS -DSAFE_ALGO -DTIME_EXEC
CCFLAGPROJ+= -DRSVS_ACCESS_DEVELOPMENT_PARAMETERS

# boost/backtrace flags:
# -DBOOST_STACKTRACE_USE_BACKTRACE -DUSE_STACKTRACE
CCFLAGBOOST=

# Eigen flags:
# Release: -DEIGEN_NO_DEBUG -DEIGEN_NO_STATIC_ASSERT
# MKL usage: -DEIGEN_USE_MKL_ALL
# License requirements: -DEIGEN_MPL2_ONLY
CCFLAGEIG= -DEIGEN_NO_DEBUG -DEIGEN_NO_STATIC_ASSERT -DEIGEN_MPL2_ONLY

# Tetgen flag to use as library: -DTETLIBRARY
CCFLAGETETGEN= -DTETLIBRARY
CCFLAGUSER=

# Flags to add when the target is "testall"
# -DTEST_ALL_WORKING
CFTESTFLAG= -DTEST_ALL
CF_NORSVSTESTFLAG= -DRSVS_NOTESTS
CF_DISTRIBUTIONFLAG= -DRSVS_HIDETESTS
