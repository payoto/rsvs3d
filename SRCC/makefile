# The default settings (flags, dependencies, libraries) for linux
# and windows are stored in their own makefiles.
ifeq ($(OS),Windows_NT)
	include makewin.mk
else
	include makelinux.mk
endif

# This flag is meant to be controlled from the command line and governs
# the behaviour of polyscope
HEADLESS=false

#------------------------------------------------------
#          PARSE compiler flags
#------------------------------------------------------
# flags are grouped from the variable in which they are called to the ones
# in which they are called. There can be dependencies on the target and other
# flags

ifeq ($(MAKECMDGOALS), distribution)
	# distribution requires optimisation and No debug and test flags
	OPTIMLVL= 3
	CCFLAGDBG=
	CCFLAGPROJ= -DTIME_EXEC $(CF_NORSVSTESTFLAG) $(CF_DISTRIBUTIONFLAG)
	CCFLAGEIG= -DEIGEN_NO_DEBUG -DEIGEN_NO_STATIC_ASSERT -DEIGEN_MPL2_ONLY
endif

ifeq ($(DEBUG), true)
	OPTIMLVL= 0
	CCFLAGDBG= -g -ggdb
	CCFLAGPROJ= -DRSVS_DEBUG -DSAFE_ACCESS -DSAFE_ALGO -DTIME_EXEC
	CCFLAGPROJ+= -DRSVS_ACCESS_DEVELOPMENT_PARAMETERS
	CCFLAGEIG=
endif

# O3 comes with spurious warnings which is ugly to deal in the code
ifeq ($(OPTIMLVL), 3)
WARNFLAGS += -Wno-maybe-uninitialized
endif

CCFLAGCUSTOM= $(CCFLAGDBG) $(CCFLAGPROJ) $(CCFLAGBOOST) $(CCFLAGEIG) $(CCFLAGETETGEN)
# CCFLAGETETGEN is needed as Tetgen does not tolerate optimisation
# beyond O3. It will break using -march and/or -ffast-math and co...
CCFLAGETETGEN += $(CCFLAGDBG)
#------------------------------------------------------
#    PARSE correct flags to get FILESYSTEM to work
#------------------------------------------------------

ifneq (,$(filter $(CC),gcc g++))
ifneq ($(USE_BOOST),true)
	GCCVERSION = $(shell $(CC) --version | grep ^$(CC) | sed 's/^.* //g')
	GCCVERSIONMAJOR = $(shell $(CC) --version | grep ^$(CC) | sed 's/^.* //g' | sed 's/^\([0-9]\).*/\1/')

	# $(info Using $(CC) $(GCCVERSION) $(GCCVERSIONMAJOR))

	ifneq ("$(GCCVERSIONMAJOR)", "7")
		# $(warning Attempting to use GCC != v7 without Boost is not supported on windows)
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
LDFILE = stdc++ $(LIBPOLYSCOPE)
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

LDDIR = $(LDDIRBOOST) $(LDPOLYSCOPE) $(LDPOLYSCOPEDEPS)
# backtrace linker flags: # -L"C:\ap1949include\libbacktrace-build\lib" -lbacktrace
POLYSCOPELIBFILE= $(LDPOLYSCOPE)/libpolyscope.a


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
# Set default executable
EXECUTABLE = $(MAINEXECUTABLE)

ifeq ($(MAKECMDGOALS), notests)
CCFLAGCUSTOM += $(CF_NORSVSTESTFLAG)
endif
ifeq ($(MAKECMDGOALS), distribution)
CCFLAGCUSTOM += $(CF_NORSVSTESTFLAG) $(CF_DISTRIBUTIONFLAG)
endif

ifeq ($(MAKECMDGOALS), testnew)
EXECUTABLE = $(TESTEXECUTABLE)
CCFLAGCUSTOM += -DRSVSTEST
endif

ifeq ($(MAKECMDGOALS), testshort)
EXECUTABLE = $(TESTEXECUTABLE)
CCFLAGCUSTOM += -DRSVSTEST -DTEST_SHORT
endif

ifeq ($(MAKECMDGOALS), testall)
EXECUTABLE = $(TESTALLEXECUTABLE)
CCFLAGCUSTOM += -DRSVSTEST -DTEST_KNOWN_FAILURES
endif


# make polyscope be built only with the mocked backend
ifeq ($(HEADLESS), true)
# At the moment cannot build with only the mock backend.
# POLYSCOPE_CMAKE_FLAG= -DPOLYSCOPE_BACKEND_OPENGL3_GLFW=OFF
CCFLAGCUSTOM += -DHEADLESS
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
testshort : FORCE $(SRCS) $(EXECUTABLE)
notests : FORCE $(SRCS) $(EXECUTABLE)
distribution : FORCE $(SRCS) $(EXECUTABLE)

# Linker recipes for the EXECUTABLES
$(EXECUTABLE) : $(BINDIR)/$(EXECUTABLE)
	cp $(BINDIR)/$@ ./
$(BINDIR)/$(EXECUTABLE) : $(OBJECTS) $(POLYSCOPELIBFILE)
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

# Build the polyscope library
$(POLYSCOPELIBFILE):
	mkdir -p modules/polyscope/build
	cd modules/polyscope/build && cmake .. -G "Unix Makefiles" $(POLYSCOPE_CMAKE_FLAG)
	cd modules/polyscope/build && $(MAKE)

echo:
	@echo $(POLYSCOPELIBFILE)
	@echo $(LIBPOLYSCOPE)

# Folder depency used to mirror the SRCTREE into OBJDIR
# Called by `$$(@D)/.f` in the object recipes
%/.f:
	mkdir -p $(dir $@)
	touch $@

.PHONY: clean clean-deps
# Clean target
OBJDEL = $(call rwildcard,$(OBJDIR)/,*.o)
clean :
	rm *.o *.exe $(OBJDEL) $(DEPDIR)/*.d

clean-deps:
	rm -r modules/polyscope/build

FORCE:


$(DEPDIR)/%.d: ;
.PRECIOUS: $(DEPDIR)/%.d  %/.f

$(shell echo -n $(MAKECMDGOALS) > $(DEPDIR)/prevtarget.d)
# $(info $(patsubst %,$(DEPDIR)/%.d,$(notdir $(basename $(SRCS)))))
include $(wildcard $(patsubst %,$(DEPDIR)/%.d,$(notdir $(basename $(SRCS)))))
