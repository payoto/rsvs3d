
#------------------------------------------------------
#                 INPUTS
#------------------------------------------------------
# Specify compiler
CC=g++ 
# Specify linker
LINK=g++
# Compiler flags
# Warnings all and Define DEBUGLVL2 to include all check code
CFLAGS=-c -std=c++11 -Wall -Wextra -pedantic -ggdb
#CUSERFLAGS= 
CUSERFLAGS= -DDEBUGLVL1 -DSAFE_ACCESS 
#-DSAFE_ALGO
CFTESTFLAG= -DTEST_ALL
LDFLAGS=-lstdc++
INCLDIR= -I"C:\ap1949include\Eigen_all\Eigen"
OBJDIR=obj
# make all options
SOURCES=arraystructures.cpp
EXECUTABLE=main
#make testnew/testall options
TESTSOURCES=test.cpp arraystructures.cpp  snakstruct.cpp snakstruct_test.cpp postprocessing.cpp snakeengine.cpp
TESTALLSOURCES= voxel.cpp arraystructures_test.cpp
TESTEXECUTABLE=main_test

# Dependency handling
DEPDIR := .d
$(shell mkdir -p $(DEPDIR) >/dev/null)
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.Td


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
POSTCOMPILE = @mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d && touch $@
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


$(OBJDIR)/%.o : %.cpp
$(OBJDIR)/%.o : %.cpp $(DEPDIR)/%.d
	$(COMPILECPP) $< -o $@
	$(POSTCOMPILE)

#compileflagsalltest :
#CFLAGS=$(CFLAGS) $(CFTESTFLAG)
#TESTOBJECTS=$(TESTOBJECTS) $(TESTALL:.cpp=.o)

# Clean target
#.PHONY : clean
clean :
	rm *.o *.exe $(OBJDIR)/*.o $(DEPDIR)/*.d 

#precompile.hpp : 
#	$(CC) -c -H -DDEBUGLVL1 -ggdb $(INCLDIR) -x c++-header precompile.hpp


$(DEPDIR)/%.d: ;
.PRECIOUS: $(DEPDIR)/%.d

include $(wildcard $(patsubst %,$(DEPDIR)/%.d,$(basename $(SRCS))))
