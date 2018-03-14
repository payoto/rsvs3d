# Specify compiler
CC=g++ 
# Specify linker
LINK=g++
# Compiler flags
# Warnings all and Define DEBUGLVL2 to include all check code
CFLAGS=-c -Wall -DDEBUGLVL1 -ggdb
CFTESTFLAG= -DTEST_ALL
LDFLAGS=-lstdc++
INCLDIR= -I"C:\ap1949include\Eigen_all\Eigen"

SOURCES=arraystructures.cpp
EXECUTABLE=main
TESTSOURCES=test.cpp arraystructures.cpp  snakstruct.cpp snakstruct_test.cpp 
TESTALLSOURCES= voxel.cpp postprocessing.cpp arraystructures_test.cpp
TESTEXECUTABLE=main_test

OBJECTS=$(SOURCES:.cpp=.o)

TESTOBJECTS=$(TESTSOURCES:.cpp=.o)

ifeq ($(MAKECMDGOALS),testall)
TESTOBJECTS += $(TESTALLSOURCES:.cpp=.o)
CFLAGS+= $(CFTESTFLAG)
endif




#.PHONY : all
all : $(SOURCES) $(EXECUTABLE)
  
testnew : $(TESTSOURCES) $(TESTEXECUTABLE) 
  
testall : $(TESTSOURCES) $(TESTALLSOURCES) $(TESTEXECUTABLE) 
  
# Compile the source files into object files
 

main : $(OBJECTS) 
	$(LINK) $(LDFLAGS) $(OBJECTS) -o $@ 

main_test : $(TESTOBJECTS) 
	$(LINK) $(LDFLAGS) $(TESTOBJECTS) -o $@ 


.cpp.o :
	$(CC) $(CFLAGS) $(INCLDIR) $< -o $@

#compileflagsalltest :
#CFLAGS=$(CFLAGS) $(CFTESTFLAG)
#TESTOBJECTS=$(TESTOBJECTS) $(TESTALL:.cpp=.o)

# Clean target
#.PHONY : clean
clean :
	rm *.o *.exe

#precompile.hpp : 
#	$(CC) -c -H -DDEBUGLVL1 -ggdb $(INCLDIR) -x c++-header precompile.hpp
