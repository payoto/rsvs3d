#include <iostream>
#include <stdexcept>


#include "arraystructures.hpp"

// Implementation of features




// Class function definitions

void volu::disp(){
   int i;
   cout << "volu : index " << index << " | fill " << fill << ", " << 
   target << ", "<< error << " | surfind " << surfind.size();
   for (i=0; unsigned_int(i)<surfind.size();i++){
      cout << "-" << surfind[i];
   }
   cout << endl;
}

void surf::disp(){
   int i;
   cout << "surf : index " << index << " | fill " << fill << ", " << 
   target << ", "<< error << " | voluind " << voluind.size();
   for (i=0; unsigned_int(i)<voluind.size();i++){
      cout << "-" << voluind[i];
   }
   cout << " | edgeind " << edgeind.size();
   for (i=0; unsigned_int(i)<edgeind.size();i++){
      cout << "-" << edgeind[i];
   }
   cout << endl;
}

void edge::disp(){
   int i;
   cout << "edge : index " << index <<  " | vertind " << vertind.size();
   for (i=0; unsigned_int(i)<vertind.size();i++){
      cout << "-" << vertind[i];
   }
   cout << " | surfind " << surfind.size();
   for (i=0; unsigned_int(i)<surfind.size();i++){
      cout << "-" << surfind[i];
   }
   cout << endl;
}

void vert::disp(){
   int i;
   cout << "vert : index " << index <<  " | edgeind " << edgeind.size();
   for (i=0; unsigned_int(i)<edgeind.size();i++){
      cout << "-" << edgeind[i];
   }
   cout << " | coord " << coord.size();
   for (i=0; unsigned_int(i)<coord.size();i++){
      cout << "-" << coord[i];
   }
   cout << endl;
}

void mesh::HashArray(){
   verts.HashArray();
   edges.HashArray();
   surfs.HashArray();
   volus.HashArray();
}

void mesh::Init(int nVe,int nE, int nS, int nVo){

   verts.Init(nVe);
   edges.Init(nE);
   surfs.Init(nS);
   volus.Init(nVo);
   
   cout << "Mesh Correctly Assigned!" << endl;

}
// Test code
int Test_ArrayStructures() { 
   // Test the functionality provided by arraystructures

   int errFlag,errTest;


   errFlag=0;

   cout << "--------------------------------------------" << endl;
   cout << "      testing volu" << endl;
   cout << "--------------------------------------------" << endl;
   errTest=Test_Volu();
   errFlag= errFlag | (errTest!=0);

   cout << "--------------------------------------------" << endl;
   cout << "      testing surf" << endl;
   cout << "--------------------------------------------" << endl;
   errTest=Test_Surf();
   errFlag= errFlag | (errTest!=0);

   cout << "--------------------------------------------" << endl;
   cout << "      testing vert" << endl;
   cout << "--------------------------------------------" << endl;
   errTest=Test_Vert();
   errFlag= errFlag | (errTest!=0);

   cout << "--------------------------------------------" << endl;
   cout << "      testing edge" << endl;
   cout << "--------------------------------------------" << endl;
   errTest=Test_Edge();
   errFlag= errFlag | (errTest!=0);

   errFlag= errFlag | (errTest!=0);

   return(errFlag);
} 

int Test_Volu(){
   int i;

   voluarray cellStack;  // stack of ints 
   volu singleCell;

   try {

      
      singleCell.surfind.push_back(3);
      singleCell.surfind.push_back(4);

      cout << "assign 3 in cellStack" << endl;
      cellStack.elems.assign(3,singleCell); 
      
      singleCell.index=1;
      cout << "push_back in cellStack" << endl;
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
      
      cellStack.disp();
      
      cout << cellStack.elems.capacity() << endl;

      cout << "HashCellArrayStack" << endl;
      cellStack.HashArray();
      


   } catch (exception const& ex) { 
      cerr << "Exception: " << ex.what() <<endl; 
      return -1;
   } 
   return(0);

}

int Test_Surf(){
   int i;

   surfarray cellStack;  // stack of ints 
   surf singleCell;

   try {

      for (i=0;i<5;i++){
         singleCell.edgeind.push_back(rand()%100+1);
      }
      
      for (i=0;i<2;i++){
         singleCell.voluind.push_back(rand()%100+1);
      }
      

      cout << "assign 3 in cellStack" << endl;
      cellStack.elems.assign(3,singleCell); 
      
      singleCell.index=1;
      cout << "push_back in cellStack" << endl;
      cellStack.elems.push_back(singleCell); 

      // manipulate int stack 

      cout << "Calls to display" << endl;
      cout << "Display each volu in stack" << endl;

      for ( i = 0; i < 4; ++i)
      {
         cellStack.elems[i].index=i;
         cellStack.elems[i].edgeind[0]=i;
         cellStack[i]->disp() ;
      }

      cout << "Assign last volu to first volu" << endl;
      cellStack.elems[0]=cellStack[i-1]; 

      cout << "Display each volu in stack" << endl;
      
      cellStack.disp();
      

      cout << "HashCellArrayStack" << endl;
      cellStack.HashArray();
      


   } catch (exception const& ex) { 
      cerr << "Exception: " << ex.what() <<endl; 
      return -1;
   } 
   return(0);


}
int Test_Edge(){
   int i;

   edgearray cellStack;  // stack of ints 
   edge singleCell;

   try {

      for (i=0;i<5;i++){
         singleCell.surfind.push_back(rand()%100+1);
      }
      
      for (i=0;i<2;i++){
         singleCell.vertind.push_back(rand()%100+1);
      }
      

      cout << "assign 3 in cellStack" << endl;
      cellStack.elems.assign(3,singleCell); 
      
      singleCell.index=1;
      cout << "push_back in cellStack" << endl;
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
      cellStack.elems[0]=cellStack[i-1]; 

      cout << "Display each volu in stack" << endl;
      
      cellStack.disp();
      

      cout << "HashCellArrayStack" << endl;
      cellStack.HashArray();
      


   } catch (exception const& ex) { 
      cerr << "Exception: " << ex.what() <<endl; 
      return -1;
   } 
   return(0);

}

int Test_Vert(){
   int i;

   vertarray cellStack;  // stack of ints 
   vert singleCell;

   try {

      for (i=0;i<5;i++){
         singleCell.edgeind.push_back(rand()%100+1);
      }
      
      for (i=0;i<2;i++){
         singleCell.coord.push_back(double(rand()%100+1)/100.0);
      }
      

      cout << "assign 3 in cellStack" << endl;
      cellStack.elems.assign(3,singleCell); 
      singleCell.index=1;
      cout << "push_back in cellStack" << endl;
      cellStack.elems.push_back(singleCell); 

      // manipulate int stack 

      cout << "Calls to display" << endl;
      cout << "Display each volu in stack" << endl;

      for ( i = 0; i < 4; ++i)
      {
         cellStack.elems[i].index=i;
         cellStack.elems[i].edgeind[0]=i;
         cellStack[i]->disp() ;
      }

      cout << "Assign last volu to first volu" << endl;
      cellStack.elems[0]=cellStack[i-1]; 

      cout << "Display each volu in stack" << endl;
      
      cellStack.disp();
      

      cout << "HashCellArrayStack" << endl;
      cellStack.HashArray();
      
      cout << "Test Find " << cellStack.find(1) << endl;


   } catch (exception const& ex) { 
      cerr << "Exception: " << ex.what() <<endl; 
      return -1;
   } 
   return(0);

}


