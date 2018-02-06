# Specify compiler
CC=g++
# Specify linker
LINK=g++
# Compiler flags
# Warnings all and Define DEBUGLVL2 to include all check code
CFLAGS=-c -Wall -DDEBUGLVL2
LDFLAGS=
SOURCES=test_arraystructures.cpp
TESTSOURCES=test_arraystructures.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=make


#.PHONY : all
all : app
  
# Link the object files into a binary
app : main.o
    $(LINK) -o app main.o -lstdc++
  
# Compile the source files into object files
main.o : test_arraystructures.cpp
    $(CC) -c test_arraystructures.cpp -o main.o
  
# Clean target
#.PHONY : clean
clean :
    rm main.o app 


# CC=g++
# CFLAGS=-c -Wall
# LDFLAGS=
# SOURCES=main.cpp hello.cpp factorial.cpp
# OBJECTS=$(SOURCES:.cpp=.o)
# EXECUTABLE=hello

# all: $(SOURCES) $(EXECUTABLE)
    
# $(EXECUTABLE): $(OBJECTS) 
#     $(CC) $(LDFLAGS) $(OBJECTS) -o $@

# .cpp.o:
#     $(CC) $(CFLAGS) $< -o $@