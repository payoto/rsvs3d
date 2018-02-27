# Specify compiler
CC=g++
# Specify linker
LINK=g++
# Compiler flags
# Warnings all and Define DEBUGLVL2 to include all check code
CFLAGS=-c -Wall -DDEBUGLVL2
LDFLAGS=-lstdc++

SOURCES=arraystructures.cpp
EXECUTABLE=main
TESTSOURCES=test.cpp arraystructures.cpp
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
	$(CC) $(CFLAGS) $< -o $@


# Clean target
#.PHONY : clean
clean :
	rm *.o *.exe
