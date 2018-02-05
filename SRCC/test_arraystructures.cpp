#include <iostream>
#include <vector>
#include <cstdlib>
#include <string>
#include <stdexcept>

#include "arraystructures.hpp"



int main() { 
	ArrayStruct<int>         intStack;  // stack of ints 
    ArrayStruct<cell>         cellStack;  // stack of ints 
   int i;
   try {

      
      cell singleCell;

      cout << "reserve" << endl;
      //cellStack.elems.reserve(2);
      cout << "assign" << endl;
      cellStack.elems.assign(3,singleCell); 
      cout << "=" << endl;
      singleCell.index=1;
      cout << "push_back" << endl;
      cellStack.elems.push_back(singleCell); 

      // manipulate int stack 
      intStack.elems.push_back(7); 

      cout << "Calls to display" << endl;
      cout << intStack.elems[0] <<endl; 
      for ( i = 0; i < 4; ++i)
      {
      	cellStack.elems[i].index=i;
      	cellStack[i]->disp() ;
      }
      cout << "GAP" << endl;
      singleCell.index=0;
      cellStack[0]->disp();
      cellStack.elems[0]=cellStack[i-1]; 

      for ( i = 0; i < 4; ++i)
      {
      	cellStack[i]->disp() ;
      }
      cout << cellStack.elems.capacity() << endl;

      cellStack.HashArray();
      


   } catch (exception const& ex) { 
      cerr << "Exception: " << ex.what() <<endl; 
      return -1;
   } 
} 