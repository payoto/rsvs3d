#include <iostream>
#include <vector>
#include <cstdlib>
#include <string>
#include <stdexcept>

#include "arraystructures.hpp"



int main() { 
   try {
      ArrayStruct<int>         intStack;  // stack of ints 
      ArrayStruct<cell>         cellStack;  // stack of ints 
      cell singleCell;
      cellStack.elems.push_back(singleCell); 
      singleCell.index=1;
      cellStack.elems.push_back(singleCell); 

      // manipulate int stack 
      intStack.elems.push_back(7); 

      cout << intStack.elems[0] <<endl; 

      cout << cellStack.elems[0].disp() <<endl; 
      cout << cellStack.elems[1].disp() <<endl; 

   } catch (exception const& ex) { 
      cerr << "Exception: " << ex.what() <<endl; 
      return -1;
   } 
} 