
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
EXECUTABLE=main
TESTEXECUTABLE=main_test

# source, object and dependency directories
OBJDIR = obj
SRCTREE = src/
DEPDIR := .d
BINDIR =
INCLDIR= "incl"

# Include directories
SHAREINCLROOTDIR = C:/ap1949include
DIREIGEN= $(SHAREINCLROOTDIR)/Eigen_all/Eigen
# DIRMKL = "C:\Program Files (x86)\IntelSWTools\compilers_and_libraries_2019.1.144\windows\mkl\include"
DIRJSON= $(SHAREINCLROOTDIR)/json/single_include/nlohmann
DIRBOOST= $(SHAREINCLROOTDIR)/boost_1_67_0
DIRBACTRACE= $(SHAREINCLROOTDIR)/libbacktrace/
DIRCXXOPTS= $(SHAREINCLROOTDIR)/cxxopts-2.1.1/include

# -isystem$(DIRMKL)
INCL= 

INCLDIR+= $(DIREIGEN) $(DIRJSON) 
INCLDIRSYS = $(DIRBOOST)

#------------------------------------------------------
#                 COMPILER FLAGS
#------------------------------------------------------
# Standard GCC warning, optimisation and c++ standard flags
WARNFLAGS= -Wall -Wextra -pedantic
OPTIMLVL=3
CCPSTD=c++11

# Other compiler flags
# project flags: # -DDEBUGLVL1 -DSAFE_ACCESS  -DSAFE_ALGO -DTIME_EXEC
# gdb (debug) flags: # -g -ggdb 
# boost/backtrace flags: # -DBOOST_STACKTRACE_USE_BACKTRACE
# eigen flags: # -DEIGEN_NO_DEBUG -DEIGEN_NO_STATIC_ASSERT -march=native -DEIGEN_USE_MKL_ALL
# cxxopts #
CCFLAGDBG= -g -ggdb 
CCFLAGPROJ= -DTIME_EXEC -DDEBUGLVL1 -DSAFE_ACCESS  -DSAFE_ALGO
CCFLAGBOOST=
CCFLAGEIG= 
CCFLAGUSER= 
CCFLAGCUSTOM= $(CCFLAGDBG) $(CCFLAGPROJ) $(CCFLAGBOOST) $(CCFLAGEIG)
# Flags to add when the target is "testall"
CFTESTFLAG= -DTEST_ALL


#------------------------------------------------------
#                 LINKER FLAGS
#------------------------------------------------------
# LDFILE are the library names, LDDIR are the library 
# directories
LDFILE = stdc++

LDDIRBOOST = $(DIRBOOST)/stage/lib
LDFILEBOOST = boost_filesystem boost_system
BOOSTVERSION = -mgw63-mt-sd-x64-1_67
LDDIRBACTRACE= $(SHAREINCLROOTDIR)/libbacktrace-build/lib

LDFILEBOOSTWITHVERSION = $(addsuffix $(BOOSTVERSION),$(LDFILEBOOST))
LDDIR = $(LDDIRBOOST)
LDFILE += $(LDFILEBOOSTWITHVERSION)
# backtrace linker flags: # -L"C:\ap1949include\libbacktrace-build\lib" -lbacktrace


#------------------------------------------------------
#            PARSE FLAGS BEFORE RECIPES
#------------------------------------------------------


# Compile list of flags
CFLAGS=-c -std=$(CCPSTD) $(WARNFLAGS) -O$(OPTIMLVL) $(CCFLAGCUSTOM)
ifeq ($(MAKECMDGOALS),testall)
CFLAGS+= $(CFTESTFLAG)
endif
CFLAGS += $(CCFLAGUSER)
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$(notdir $*).Td
# Parse includes
INCL+=$(addprefix -I,$(INCLDIR))
INCL+=$(addprefix -isystem,$(INCLDIRSYS))
# Parse libraries
LDFLAGS=$(addprefix -L,$(LDDIR))
LDFLAGS+=$(addprefix -l,$(LDFILE))
# Compilation command
COMPILECPP=$(CC) $(DEPFLAGS) $(CFLAGS) $(INCL)
POSTCOMPILE = @mv -f $(DEPDIR)/$(notdir $*).Td $(DEPDIR)/$(notdir $*).d && touch $@


#------------------------------------------------------
#            PARSE FILES and FOLDERS BEFORE RECIPES
#------------------------------------------------------
# Define a recursive wildcard function, to find sources:
rwildcard=$(wildcard $1$2) $(foreach d,$(wildcard $1*),$(call rwildcard,$d/,$2))

# Make output folders
$(shell mkdir -p $(DEPDIR) >/dev/null)
$(shell mkdir -p $(OBJDIR) >/dev/null)
# Object definitions
TESTSOURCES = $(call rwildcard,$(SRCTREE),*.cpp)
SOURCES = $(call rwildcard,$(SRCTREE),*.cpp)
# testall handling
ifeq ($(MAKECMDGOALS),testall)
TESTSOURCES+= $(TESTALLSOURCES)
endif
OBJECTS=$(addprefix $(OBJDIR),$(SOURCES:.cpp=.o))
TESTOBJECTS=$(addprefix $(OBJDIR)/,$(TESTSOURCES:.cpp=.o))

# Handling dependency include:
ifeq ($(MAKECMDGOALS),testall)
SRCS = $(TESTALLSOURCES) $(TESTSOURCES)
endif
ifeq ($(MAKECMDGOALS),testnew)
SRCS =  $(TESTSOURCES)
endif
ifeq ($(MAKECMDGOALS),testall)
SRCS = $(SOURCES)
endif

#------------------------------------------------------
#             RECIPES
#------------------------------------------------------

# TARGETS
all : $(SOURCES) $(EXECUTABLE)  
testnew : $(TESTSOURCES) $(TESTEXECUTABLE) 
testall : $(TESTSOURCES) $(TESTALLSOURCES) $(TESTEXECUTABLE) 
  
# Linker recipes for the EXECUTABLES
$(EXECUTABLE) : $(OBJECTS) 
	$(LINK)  $(OBJECTS) -o $@ $(LDFLAGS)
$(TESTEXECUTABLE) : $(TESTOBJECTS) 
	$(LINK)  $(TESTOBJECTS) -o $@ $(LDFLAGS)

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
	$(COMPILECPP) $< -o $@
	$(POSTCOMPILE)

# Special case for this source and object which cannot be compiled without optimisation
# or the object file is too large on windows
$(OBJDIR)/src/rsvs/RSVScalc_SQP.o : src/rsvs/RSVScalc_SQP.cpp
$(OBJDIR)/src/rsvs/RSVScalc_SQP.o : src/rsvs/RSVScalc_SQP.cpp $(DEPDIR)/RSVScalc_SQP.d $$(@D)/.f 
	$(COMPILECPP) -O3 $< -o $@
	$(POSTCOMPILE)

# Folder depency used to mirror the SRCTREE into OBJDIR
# Called by `$$(@D)/.f` in the object recipes 
%/.f:
	mkdir -p $(dir $@)
	touch $@

# Clean target
OBJDEL = $(call rwildcard,$(OBJDIR)/,*.o)
clean :
	rm *.o *.exe $(OBJDEL) $(DEPDIR)/*.d 

#precompile.hpp : 
#	$(CC) -c -H -DDEBUGLVL1 -ggdb $(INCLDIR) -x c++-header precompile.hpp


$(DEPDIR)/%.d: ;
.PRECIOUS: $(DEPDIR)/%.d  %/.f


# $(info $(patsubst %,$(DEPDIR)/%.d,$(notdir $(basename $(SRCS)))))
include $(wildcard $(patsubst %,$(DEPDIR)/%.d,$(notdir $(basename $(SRCS)))))
