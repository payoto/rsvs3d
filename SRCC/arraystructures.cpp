#include <iostream>
#include <stdexcept>
#include <sstream>
#include <functional>

#include "arraystructures.hpp"

// Implementation of features


template <class T> bool CompareDisp(T *mesh1,T *mesh2){
   bool compFlag;
   stringstream ss1,ss2;
   auto old_buf = cout.rdbuf(ss1.rdbuf()); 

   mesh1->disp();
   cout.rdbuf(ss2.rdbuf()); 
   mesh2->disp();
   std::cout.rdbuf(old_buf);

   compFlag=ss1.str().compare(ss2.str())==0;
   return(compFlag);
}


bool CompareFuncOut(function<void()> func1, function<void()> func2){
   bool compFlag;
   stringstream ss1,ss2;
   auto old_buf = cout.rdbuf(ss1.rdbuf()); 

   func1();
   cout.rdbuf(ss2.rdbuf()); 
   func2();
   std::cout.rdbuf(old_buf);

   compFlag=ss1.str().compare(ss2.str())==0;
   return(compFlag);
}



// Class function definitions
// Methods of meshpart : volu surf edge vert
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

void volu::ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu){
   int i;
   index+=nVolu;
   for (i=0; unsigned_int(i)<surfind.size();i++){
      surfind[i]=surfind[i]+nSurf;
   }
}
void surf::ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu){
   int i;
   index+=nSurf;
   for (i=0; unsigned_int(i)<voluind.size();i++){
      voluind[i]=voluind[i]=nVolu;
   }

   for (i=0; unsigned_int(i)<edgeind.size();i++){
      edgeind[i]=edgeind[i]+nEdge;
   }
}
void edge::ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu){
   int i;
   index+=nEdge;
   for (i=0; unsigned_int(i)<vertind.size();i++){
      vertind[i]=vertind[i]+nVert;
   }

   for (i=0; unsigned_int(i)<surfind.size();i++){
      surfind[i]=surfind[i]+nSurf;
   }
}
void vert::ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu){
   int i;
   index+=nVert;
   for (i=0; unsigned_int(i)<edgeind.size();i++){
      edgeind[i]=edgeind[i]+nEdge;
   }
}
// methods for mesh
void mesh::HashArray(){
   verts.HashArray();
   edges.HashArray();
   surfs.HashArray();
   volus.HashArray();
}
void mesh::SetMaxIndex(){
   verts.SetMaxIndex();
   edges.SetMaxIndex();
   surfs.SetMaxIndex();
   volus.SetMaxIndex();
}
void mesh::GetMaxIndex(int *nVert,int *nEdge,int *nSurf,int *nVolu){
   *nVert=verts.GetMaxIndex();
   *nEdge=edges.GetMaxIndex();
   *nSurf=surfs.GetMaxIndex();
   *nVolu=volus.GetMaxIndex();
}
void mesh::disp(){
   verts.disp();
   edges.disp();
   surfs.disp();
   volus.disp();
}

void mesh::Init(int nVe,int nE, int nS, int nVo){

   verts.Init(nVe);
   edges.Init(nE);
   surfs.Init(nS);
   volus.Init(nVo);

   #ifdef TEST_ARRAYSTRUCTURES
   cout << "Mesh Correctly Assigned!" << endl;
   #endif // TEST_ARRAYSTRUCTURES
}

void mesh::MakeCompatible_inplace(mesh *other){
   // Makes other mesh compatible with this to be 
   // merged without index crashes

   int nVert,nEdge,nVolu,nSurf;

   // Define Max indices in current mesh
   this->GetMaxIndex(&nVert,&nEdge,&nVolu,&nSurf);
   other->ChangeIndices(nVert,nEdge,nVolu,nSurf);
}

void mesh::ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu){
   int ii;
   // Volumes
   for(ii=0;ii<volus.size();++ii){
      volus.elems[ii].ChangeIndices(nVert,nEdge,nSurf,nVolu);
   }
   // Surfaces
   for(ii=0;ii<surfs.size();++ii){
      surfs.elems[ii].ChangeIndices(nVert,nEdge,nSurf,nVolu);
   }
   // Volumes
   for(ii=0;ii<edges.size();++ii){
      edges.elems[ii].ChangeIndices(nVert,nEdge,nSurf,nVolu);
   }
   // Surfaces
   for(ii=0;ii<verts.size();++ii){
      verts.elems[ii].ChangeIndices(nVert,nEdge,nSurf,nVolu);
   }
}

mesh mesh::MakeCompatible(mesh other){
   MakeCompatible_inplace(&other);
   return(other);
}

void mesh::Concatenate(mesh *other){

   this->volus.Concatenate(&other->volus);
   this->edges.Concatenate(&other->edges);
   this->verts.Concatenate(&other->verts);
   this->surfs.Concatenate(&other->surfs);
}

void PopulateIndices(mesh *meshin){
   
   meshin->volus.PopulateIndices();
   meshin->edges.PopulateIndices();
   meshin->verts.PopulateIndices();
   meshin->surfs.PopulateIndices();
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

   cout << "--------------------------------------------" << endl;
   cout << "      testing Mesh" << endl;
   cout << "--------------------------------------------" << endl;
   errTest=Test_Mesh();
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
      cellStack.SetMaxIndex();
      cout << "Max Index in Array: " << cellStack.GetMaxIndex() <<  endl;
      
      cellStack.Concatenate(&cellStack);
      cellStack.disp();

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

int Test_Mesh(){
   mesh mesh1,mesh2,mesh3;
   int errFlag=0;
   bool test1;

   try {

      mesh1.Init(8,12,6,1);
      mesh2.Init(14,20,11,2);

      PopulateIndices(&mesh1);
      PopulateIndices(&mesh2);
      mesh3=mesh1.MakeCompatible(mesh2);
#ifdef TEST_ARRAYSTRUCT_MESH
      cout << "mesh1 " << endl;
      mesh1.disp() ;
      cout << "Compatible mesh #3  generated displaying ";
      cout << "mesh2 " << endl;
      mesh2.disp() ;
#endif

      mesh1.MakeCompatible_inplace(&mesh2);

#ifdef TEST_ARRAYSTRUCT_MESH
      cout << "mesh2 made compatible inplace: ";
      cout << "mesh2 " << endl;
      
      cout << "mesh3 " << endl;
      mesh2.disp() ;
      mesh3.disp() ;
#endif
      test1=CompareDisp(&mesh2,&mesh3);
      errFlag=errFlag || (!test1);
      cout << "Result of Comparison 2&3: " << test1 << endl;
      test1=CompareDisp(&mesh1,&mesh2);
      cout << "Result of Comparison 1&2: " << test1 << endl; 
      errFlag=errFlag || (test1);

      cout << "Result of Comparison Test_Volu&Test_Volu: " << CompareFuncOut(Test_Volu,Test_Volu) << endl;

      cout << "Concatenate mesh 1 and 2" << endl;
      mesh1.Concatenate(&mesh2);
      mesh1.disp();


   } catch (exception const& ex) { 
      cerr << "Exception: " << ex.what() <<endl; 
      errFlag= -1;
   } 
   return(errFlag);
}


