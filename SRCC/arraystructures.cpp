#include <iostream>
#include <stdexcept>

#include "arraystructures.hpp"


int test_arraystructures() { 
   // Test the functionality provided by arraystructures
   int i;

	ArrayStruct<int>         intStack;  // stack of ints 
   voluarray         cellStack;  // stack of ints 
   volu singleCell;

   try {

      

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
   return(0);
} 

// member function definition template <class T> : "ArrayStruct"
template<class T> void ArrayStruct <T>::HashArray(){
   hashTable.reserve(elems.size());
   for (int i = 0; i < int(elems.size()); ++i)
   {
      hashTable.emplace(elems[i].Key(),i);
   }
   cout << "Array Struct Succesfully Hashed" << endl;
}