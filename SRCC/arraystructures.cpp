#include <iostream>
#include <stdexcept>

#include "arraystructures.hpp"


int test_arraystructures() { 
   // Test the functionality provided by arraystructures
   int errFlag,errTest;

   errFlag=0;

   errTest=test_volu();
   errFlag= errFlag | (errTest!=0);

   return(0);
} 

int test_volu(){
   int i;

   voluarray cellStack;  // stack of ints 
   volu singleCell;

   try {

      
      singleCell.surfind.push_back(3);
      singleCell.surfind.push_back(4);

      cout << "assign 3 in cellStack" << endl;
      cellStack.elems.assign(3,singleCell); 
      cout << "=" << endl;
      singleCell.index=1;
      cout << "push_back of a volu in cellStack" << endl;
      cellStack.elems.push_back(singleCell); 

      // manipulate int stack 

      cout << "Calls to display" << endl;
      cout << "Display each volu in stack" << endl;

      for ( i = 0; i < 4; ++i)
      {
         cellStack.elems[i].index=i;
         cellStack.elems[i].surfind[0]=i;
         cellStack[i]->disp() ;
      }

      cout << "Assign last volu to first volu" << endl;
      cellStack.elems[0]=cellStack.elems[i-1]; 

      cout << "Display each volu in stack" << endl;
      for ( i = 0; i < 4; ++i)
      {
         cellStack[i]->disp() ;
         cellStack.elems[i].disp() ;
      }
      cellStack.disp2();
      cout << cellStack.elems.capacity() << endl;

      cout << "HashCellArrayStack" << endl;
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

template<class T> void ArrayStruct <T>::disp2(){
   for (int ii ; unsigned_int(ii)<elems.size();ii++){
      cout << "D2:" ;
      elems[ii].disp();
   }
}