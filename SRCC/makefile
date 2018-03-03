# Specify compiler
CC=g++
# Specify linker
LINK=g++
# Compiler flags
# Warnings all and Define DEBUGLVL2 to include all check code
CFLAGS=-c -Wall -DDEBUGLVL1 -ggdb
LDFLAGS=-lstdc++
#INCLDIR=-include Eigen.gch
#-include other/Eigen.gch
INCLDIR= -I"gch" -I"C:\ap1949include\Eigen_all\Eigen"

SOURCES=arraystructures.cpp
EXECUTABLE=main
TESTSOURCES=test.cpp arraystructures.cpp voxel.cpp
TESTEXECUTABLE=main_test

OBJECTS=$(SOURCES:.cpp=.o)
TESTOBJECTS=$(TESTSOURCES:.cpp=.o)




#.PHONY : all
all : precompile  $(SOURCES) $(EXECUTABLE)

test : precompile $(TESTSOURCES) $(TESTEXECUTABLE) 
  
# Compile the source files into object files
 

main : $(OBJECTS) 
    $(LINK) $(LDFLAGS) $(OBJECTS) -o $@ 

main_test : $(TESTOBJECTS) 
    $(LINK) $(LDFLAGS) $(TESTOBJECTS) -o $@ 

precompile : 
	$(CC) $(CFLAGS) $(INCLDIR) -x c++-header precompile.hpp

.cpp.o :
	$(CC) $(CFLAGS) $(INCLDIR) $< -o $@

# Clean target
#.PHONY : clean
clean :
	rm *.o *.exe
