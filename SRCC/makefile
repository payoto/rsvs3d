# Specify compiler
CC=g++
# Specify linker
LINK=g++
# Compiler flags
# Warnings all and Define DEBUGLVL2 to include all check code
CFLAGS=-c -Wall -DDEBUGLVL1 -O3
LDFLAGS=-lstdc++
INCLDIR= -I"C:\ap1949include\Eigen_all\Eigen"

SOURCES=arraystructures.cpp
EXECUTABLE=main
TESTSOURCES=test.cpp arraystructures.cpp voxel.cpp postprocessing.cpp
TESTEXECUTABLE=main_test

OBJECTS=$(SOURCES:.cpp=.o)
TESTOBJECTS=$(TESTSOURCES:.cpp=.o)




#.PHONY : all
all : $(SOURCES) $(EXECUTABLE)

test : $(TESTSOURCES) $(TESTEXECUTABLE) 
  
# Compile the source files into object files
 

main : $(OBJECTS) 
    $(LINK) $(LDFLAGS) $(OBJECTS) -o $@ 

main_test : $(TESTOBJECTS) 
    $(LINK) $(LDFLAGS) $(TESTOBJECTS) -o $@ 


.cpp.o :
	$(CC) $(CFLAGS) $(INCLDIR) $< -o $@

# Clean target
#.PHONY : clean
clean :
	rm *.o *.exe

#precompile.hpp : 
#	$(CC) -c -H -DDEBUGLVL1 -ggdb $(INCLDIR) -x c++-header precompile.hpp
