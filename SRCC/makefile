
#------------------------------------------------------
#                 TARGET COMPILER
#------------------------------------------------------
# Specify compiler (only g++ is supported)
CC=g++ 
# Specify linker
LINK=g++

#------------------------------------------------------
#                 DIRECTORIES and INCLUDES
#------------------------------------------------------
# Executable name
MAINEXECUTABLE=RSVS3D
TESTEXECUTABLE=test_RSVS3D
TESTALLEXECUTABLE=testall_RSVS3D

# source, object and dependency directories
OBJDIR = obj
SRCTREE = src
DEPDIR := .d
BINDIR = bin
INCLPROJ= incl
#modules/tetgen/

# Sources which contain the main function for all and test
MAINALLSRC = src/main.cpp
MAINTESTSRC = src/test/test.cpp

# NOTNEEDEDNEWTEST = src/rsvs/RSVScalc.cpp

# Include directories
SHAREINCLROOTDIR = C:/ap1949include
DIREIGEN= $(SHAREINCLROOTDIR)/Eigen_all/Eigen
DIRJSON= $(SHAREINCLROOTDIR)/json/single_include/nlohmann
DIRCXXOPTS= $(SHAREINCLROOTDIR)/cxxopts-2.1.1/include
DIRBACTRACE= $(SHAREINCLROOTDIR)/libbacktrace/
DIRTETGEN= modules/tetgen

#------------------------------------------------------
#                 USE OF OPTIONAL LIBRARIES
#------------------------------------------------------
# Some libraries are optional in certain conditions
# Specify by "true" or "false" whether they should be used

# Boost is needed for the <boost/filesystem> header 
# except under std=c++17 which added a standard filesystem header
# (Broken on GCC v8.1 + MinGW, must still use boost)
USE_BOOST = true
DIRBOOST= $(SHAREINCLROOTDIR)/boost_1_67_0
# BOOST can be ignored if using c++17 (not on GCC v8.1 on Win (25/01/19))
# boost compiled library file location (for linker)
LDDIRBOOST = $(DIRBOOST)/stage/lib
# boost libraries to load: boost_<library>
LDFILEBOOST = boost_filesystem boost_system
# boost version appropriate to the system 
BOOSTVERSION = -mgw63-mt-sd-x64-1_67

# intel MKL libraries for faster Eigen maths 
# Not compatible with GCC on windows
DIRMKL = "C:\Program Files (x86)\IntelSWTools\compilers_and_libraries_2019.1.144\windows\mkl\include"
# -isystem$(DIRMKL)
INCL= 

SRCTREE+= $(DIRTETGEN)
INCLDIR= $(INCLPROJ) $(DIRJSON) $(DIRCXXOPTS) $(DIRTETGEN)
INCLDIRSYS = $(DIREIGEN) $(DIRBOOST)

#------------------------------------------------------
#                 COMPILER FLAGS
#------------------------------------------------------
# Standard GCC warning, optimisation and c++ standard flags
WARNFLAGS= -Wall -Wextra -pedantic
# -Werror
#
OPTIMLVL=3
CCPSTD=11

# O3 comes with spurious warnings which is ugly to deal in the code
ifeq ($(OPTIMLVL), 3)
WARNFLAGS += -Wno-maybe-uninitialized
endif

# Other compiler flags
# project flags: # -DDEBUGLVL1 -DSAFE_ACCESS  -DSAFE_ALGO -DTIME_EXEC -DRSVS_DIAGNOSTIC -DRSVS_VERBOSE
# -DRSVS_DIAGNOSTIC_RESOLVED -DRSVS_ACCESS_DEVELOPMENT_PARAMETERS
# gdb (debug) flags: # -g -ggdb 
# boost/backtrace flags: # -DBOOST_STACKTRACE_USE_BACKTRACE -DUSE_STACKTRACE
# eigen flags: # -DEIGEN_NO_DEBUG -DEIGEN_NO_STATIC_ASSERT -march=native -DEIGEN_USE_MKL_ALL
# -funroll-all-loops -ffast-math -ffinite-math-only 
# Tetgen flag to use as library: -DTETLIBRARY
# cxxopts #
CCFLAGDBG= -g -ggdb
# -DDEBUGLVL1 -DSAFE_ACCESS -DSAFE_ALGO -DTIME_EXEC -DRSVS_DIAGNOSTIC
CCFLAGPROJ= -DRSVS_DEBUG -DSAFE_ACCESS -DSAFE_ALGO -DTIME_EXEC -DRSVS_ACCESS_DEVELOPMENT_PARAMETERS
CCFLAGBOOST= 
CCFLAGEIG= -DEIGEN_NO_DEBUG -DEIGEN_NO_STATIC_ASSERT
CCFLAGETETGEN= -DTETLIBRARY 
CCFLAGUSER= 
CCFLAGCUSTOM= $(CCFLAGDBG) $(CCFLAGPROJ) $(CCFLAGBOOST) $(CCFLAGEIG) $(CCFLAGETETGEN)
# CCFLAGETETGEN is needed as Tetgen does not tolerate optimisation 
# beyond O3. It will break using -march and/or -ffast-math and co...
CCFLAGETETGEN += $(CCFLAGDBG)
# Flags to add when the target is "testall"
# -DTEST_ALL_WORKING
CFTESTFLAG= -DTEST_ALL 

#------------------------------------------------------
#    PARSE correct flags to get FILESYSTEM to work
#------------------------------------------------------

ifneq (,$(filter $(CC),gcc g++))
ifneq ($(USE_BOOST),true)
GCCVERSION = $(shell $(CC) --version | grep ^$(CC) | sed 's/^.* //g')
GCCVERSIONMAJOR = $(shell $(CC) --version | grep ^$(CC) | sed 's/^.* //g' | sed 's/^\([0-9]\).*/\1/')

# $(info Using $(CC) $(GCCVERSION) $(GCCVERSIONMAJOR))

ifneq "$(GCCVERSIONMAJOR)" "7"
$(warning Attempting to use GCC != v7 without Boost is not supported on windows)
# space to define the flag specifying the correct file system
CCFLAGCUSTOM += -DUSE_CSTD_FILESYSTEM_EXPERIMENTAL
else
# space to define the flag specifying the correct file system
# CCFLAGCUSTOM += -DUSE_CSTD_FILESYSTEM
endif
endif
endif

#------------------------------------------------------
#                 LINKER FLAGS
#------------------------------------------------------
# LDFILE are the library names, LDDIR are the library 
# directories
LDFILE = stdc++ 
LDDIRBACTRACE= $(SHAREINCLROOTDIR)/libbacktrace-build/lib
# optional boost paths (not needed if using c++17)

ifeq ($(USE_BOOST),true)
LDFILEBOOSTWITHVERSION = $(addsuffix $(BOOSTVERSION),$(LDFILEBOOST))
LDFILE += $(LDFILEBOOSTWITHVERSION)
CCFLAGCUSTOM += -DUSE_BOOST
else
# Leave empty
LDDIRBOOST = 
endif

LDDIR = $(LDDIRBOOST) 
# backtrace linker flags: # -L"C:\ap1949include\libbacktrace-build\lib" -lbacktrace



#------------------------------------------------------
#            PARSE FLAGS BEFORE RECIPES
#------------------------------------------------------

LDFLAGS = 

# Compile list of flags
CFLAGS=-c --std=c++$(CCPSTD) $(WARNFLAGS) -O$(OPTIMLVL) 
ifeq ($(MAKECMDGOALS),testall)
CFLAGS+= $(CFTESTFLAG)
endif
CFLAGS += $(CCFLAGUSER)
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$(notdir $*).Td
# Parse includes
INCL+=$(addprefix -I,$(INCLDIR))
INCL+=$(addprefix -isystem,$(INCLDIRSYS))
# Parse libraries
LDFLAGS+=$(addprefix -L,$(LDDIR))
LDFLAGS+=$(addprefix -l,$(LDFILE))
# Compilation command
COMPILECPP=$(CC) $(DEPFLAGS) $(CFLAGS) $(INCL)
POSTCOMPILE = @mv -f $(DEPDIR)/$(notdir $*).Td $(DEPDIR)/$(notdir $*).d && touch $@


#------------------------------------------------------
#            PARSE FILES and FOLDERS BEFORE RECIPES
#------------------------------------------------------
# Each target requires to choose between main and test.cpp
# and therefore needs its own build statement
# Define a recursive wildcard function, to find sources:
rwildcard=$(wildcard $1$2) $(foreach d,$(wildcard $1*),$(call rwildcard,$d/,$2))

# Make output folders
$(shell mkdir -p $(DEPDIR) >/dev/null)
$(shell mkdir -p $(OBJDIR) >/dev/null)
$(shell mkdir -p $(BINDIR) >/dev/null)
# Object definitions
SOURCES = $(foreach d,$(wildcard $(SRCTREE)*),$(call rwildcard,$d/,*.cpp))
# SOURCES = $(call rwildcard,$(SRCTREE),*.cpp)
# $(info $(SOURCES))
# $(info sources printed)
# Handling target differences
SRCS := $(SOURCES)
# 
ifeq ($(MAKECMDGOALS), )
EXECUTABLE = $(MAINEXECUTABLE)
endif
ifeq ($(MAKECMDGOALS), all)
EXECUTABLE = $(MAINEXECUTABLE)
endif

ifeq ($(MAKECMDGOALS), testnew)
EXECUTABLE = $(TESTEXECUTABLE)
CCFLAGCUSTOM += -DRSVSTEST 
endif

ifeq ($(MAKECMDGOALS), testall)
EXECUTABLE = $(TESTALLEXECUTABLE)
CCFLAGCUSTOM += -DRSVSTEST -DTEST_KNOWN_FAILURES
endif

# $(info $(SRCS))
# $(info sources printed after filter)

# Define objects from sources
OBJECTS=$(addprefix $(OBJDIR)/,$(SRCS:.cpp=.o))

# Allows Conditional build of objects which include "makeontargetchange.h"
PREVIOUSTARGET = $(shell cat $(DEPDIR)/prevtarget.d)
ifneq ($(MAKECMDGOALS), $(PREVIOUSTARGET))
$(shell echo // overwritten when make target changes > $(INCLPROJ)/makeontargetchange.h)
endif



#------------------------------------------------------
#             RECIPES
#------------------------------------------------------

# TARGETS
all : FORCE $(SRCS) $(EXECUTABLE)  
testnew : FORCE $(SRCS) $(EXECUTABLE) 
testall : FORCE $(SRCS) $(EXECUTABLE) 
  
# Linker recipes for the EXECUTABLES
$(EXECUTABLE) : $(BINDIR)/$(EXECUTABLE)
	cp $(BINDIR)/$@ ./
$(BINDIR)/$(EXECUTABLE) : $(OBJECTS) 
	$(LINK)  $(OBJECTS) -o $@ $(LDFLAGS)

# Needed to mirror the dirrectory tree of SRCTREE
# into OBJDIR
.SECONDEXPANSION:

# Compilation recipes to build the object files
# prerequisites : 
#	1. the corresponding .cpp file
#	2. the corresponding .d file (a makefile included in this makefile
#       which specifies the includes on which the object depends)
#	3. the corresponding .f file which marks the existence of the folder
# 	   in which the object file will go.
#
# recipe: 
# 	COMPILECPP (compiler + flags) $< (the source) -o (output flag) $@ (the object)
# 	POSTCOMPILE (Dependency actions and flags)
$(OBJDIR)/%.o : %.cpp
$(OBJDIR)/%.o : %.cpp $(DEPDIR)/%.d $$(@D)/.f 
	$(COMPILECPP) $(CCFLAGCUSTOM) $< -o $@
	$(POSTCOMPILE)

# Special case for this source and object which cannot be compiled without optimisation
# or the object file is too large on windows
$(OBJDIR)/src/rsvs/RSVScalc_SQP.o : src/rsvs/RSVScalc_SQP.cpp
$(OBJDIR)/src/rsvs/RSVScalc_SQP.o : src/rsvs/RSVScalc_SQP.cpp $(DEPDIR)/RSVScalc_SQP.d $$(@D)/.f 
	$(COMPILECPP) $(CCFLAGCUSTOM) -O3 $< -o $@
	$(POSTCOMPILE)

# #  Always on optimisation and disabled warnings in the tetgen module
$(OBJDIR)/modules/tetgen/tetgen.o : modules/tetgen/tetgen.cpp
$(OBJDIR)/modules/tetgen/tetgen.o : modules/tetgen/tetgen.cpp $(DEPDIR)/tetgen.d $$(@D)/.f 
	$(COMPILECPP) $(CCFLAGETETGEN) -O3 -D_CRT_SECURE_NO_WARNINGS $< -o $@
	$(POSTCOMPILE)

$(OBJDIR)/modules/tetgen/predicates.o : modules/tetgen/predicates.cpp
$(OBJDIR)/modules/tetgen/predicates.o : modules/tetgen/predicates.cpp $(DEPDIR)/predicates.d $$(@D)/.f 
	$(COMPILECPP) $(CCFLAGETETGEN) -O3 -D_CRT_SECURE_NO_WARNINGS $< -o $@
	$(POSTCOMPILE)

# Folder depency used to mirror the SRCTREE into OBJDIR
# Called by `$$(@D)/.f` in the object recipes 
%/.f:
	mkdir -p $(dir $@)
	touch $@

.PHONY: clean 
# Clean target
OBJDEL = $(call rwildcard,$(OBJDIR)/,*.o)
clean :
	rm *.o *.exe $(OBJDEL) $(DEPDIR)/*.d 

FORCE: 


$(DEPDIR)/%.d: ;
.PRECIOUS: $(DEPDIR)/%.d  %/.f

$(shell echo -n $(MAKECMDGOALS) > $(DEPDIR)/prevtarget.d)
# $(info $(patsubst %,$(DEPDIR)/%.d,$(notdir $(basename $(SRCS)))))
include $(wildcard $(patsubst %,$(DEPDIR)/%.d,$(notdir $(basename $(SRCS)))))
