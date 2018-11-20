
#------------------------------------------------------
#                 INPUTS
#------------------------------------------------------
# Specify compiler
CC=g++ 
# Specify linker
LINK=g++
# Compiler flags
# Warnings all and Define DEBUGLVL2 to include all check code

CFLAGS=-c -std=c++11 -Wall -Wextra -pedantic -DTIME_EXEC -g -ggdb -DDEBUGLVL1 -DSAFE_ACCESS  -DSAFE_ALGO
#-g -ggdb -DDEBUGLVL1 -DSAFE_ACCESS  -DSAFE_ALGO
CUSERFLAGS= 
# -DEIGEN_NO_DEBUG -DEIGEN_NO_STATIC_ASSERT 
# -DDEBUGLVL1 -DSAFE_ACCESS  -DSAFE_ALGO
#-DSAFE_ALGO
CFTESTFLAG= -DTEST_ALL
LDFLAGS=-lstdc++
INCLDIR= -I"C:\ap1949include\Eigen_all\Eigen" -I"incl"
OBJDIR=obj
# make all options
SOURCES=arraystructures.cpp
EXECUTABLE=main
#make testnew/testall options
SRCTREE = src/

TESTEXECUTABLE=main_test

# Dependency handling
DEPDIR := .d


# Make does not offer a recursive wildcard function, so here's one:
rwildcard=$(wildcard $1$2) $(foreach d,$(wildcard $1*),$(call rwildcard,$d/,$2))

TESTSOURCES = $(call rwildcard,$(SRCTREE),*.cpp)
SOURCES = $(call rwildcard,$(SRCTREE),*.cpp)

# TESTSOURCES=test.cpp mesh.cpp snake.cpp arraystructures.cpp snakstruct_test.cpp postprocessing.cpp
# TESTSOURCES+= snakeengine.cpp snakevel.cpp RSVSmath_automatic.cpp RSVSmath.cpp meshrefinement.cpp
# TESTSOURCES+= RSVSinterface.cpp RSVSalgorithm.cpp RSVSintegration.cpp voxel.cpp
# TESTALLSOURCES= voxel.cpp arraystructures_test.cpp

# Dependency handling
$(shell mkdir -p $(DEPDIR) >/dev/null)
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$(notdir $*).Td


#------------------------------------------------------
#            CODE BEFORE RECIPES
#------------------------------------------------------

# testall handling
ifeq ($(MAKECMDGOALS),testall)
TESTSOURCES+= $(TESTALLSOURCES)
CFLAGS+= $(CFTESTFLAG)
endif
CFLAGS += $(CUSERFLAGS)
# Object definitions

OBJECTS=$(addprefix $(OBJDIR),$(SOURCES:.cpp=.o))
TESTOBJECTS=$(addprefix $(OBJDIR)/,$(TESTSOURCES:.cpp=.o))
$(shell mkdir -p $(OBJDIR) >/dev/null)

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

# Compilation command

COMPILECPP=$(CC) $(DEPFLAGS) $(CFLAGS) $(INCLDIR)
POSTCOMPILE = @mv -f $(DEPDIR)/$(notdir $*).Td $(DEPDIR)/$(notdir $*).d && touch $@
#------------------------------------------------------
#             RECIPES
#------------------------------------------------------
#.PHONY : all
all : $(SOURCES) $(EXECUTABLE)
  
testnew : $(TESTSOURCES) $(TESTEXECUTABLE) 
  
testall : $(TESTSOURCES) $(TESTALLSOURCES) $(TESTEXECUTABLE) 
  
# Compile the source files into object files
 

main : $(OBJECTS) 
	$(LINK) $(LDFLAGS) $(OBJECTS) -o $@ 

main_test : $(TESTOBJECTS) 
	$(LINK) $(LDFLAGS) $(TESTOBJECTS) -o $@ 

.SECONDEXPANSION:

$(OBJDIR)/%.o : %.cpp
$(OBJDIR)/%.o : %.cpp $(DEPDIR)/%.d $$(@D)/.f 
	$(COMPILECPP) $< -o $@
	$(POSTCOMPILE)

$(OBJDIR)/src/rsvs/RSVScalc_SQP.o : src/rsvs/RSVScalc_SQP.cpp
$(OBJDIR)/src/rsvs/RSVScalc_SQP.o : src/rsvs/RSVScalc_SQP.cpp $(DEPDIR)/RSVScalc_SQP.d $$(@D)/.f 
	$(COMPILECPP) -O3 $< -o $@
	$(POSTCOMPILE)

%/.f:
	mkdir -p $(dir $@)
	touch $@

#compileflagsalltest :
#CFLAGS=$(CFLAGS) $(CFTESTFLAG)
#TESTOBJECTS=$(TESTOBJECTS) $(TESTALL:.cpp=.o)

# Clean target
#.PHONY : clean
OBJDEL = $(call rwildcard,$(OBJDIR)/,*.o)
clean :
	rm *.o *.exe $(OBJDEL) $(DEPDIR)/*.d 

#precompile.hpp : 
#	$(CC) -c -H -DDEBUGLVL1 -ggdb $(INCLDIR) -x c++-header precompile.hpp


$(DEPDIR)/%.d: ;
.PRECIOUS: $(DEPDIR)/%.d  %/.f


# $(info $(patsubst %,$(DEPDIR)/%.d,$(notdir $(basename $(SRCS)))))
include $(wildcard $(patsubst %,$(DEPDIR)/%.d,$(notdir $(basename $(SRCS)))))
