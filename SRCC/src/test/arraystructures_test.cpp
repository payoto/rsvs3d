#include <iostream>
#include <stdexcept>
#include <sstream>
#include <functional>

#include "mesh.hpp"
#include "postprocessing.hpp"
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
   const volu *cellPtr;

   try {


   	singleCell.surfind.push_back(3);
   	singleCell.surfind.push_back(4);

   	cout << "assign 3 in cellStack" << endl;
   	cellStack.assign(3,singleCell); 

   	singleCell.index=1;
   	cout << "push_back in cellStack" << endl;
   	cellStack.push_back(singleCell); 

      // manipulate int stack 

   	cout << "Calls to display" << endl;
   	cout << "Display each volu in stack" << endl;


   	for ( i = 0; i < 4; ++i)
   	{
   		cellStack[i].index=i;
   		cellStack[i].surfind[0]=i;
   		cellStack[i].disp() ;
   	}

   	cout << "Assign last volu to first volu" << endl;
   	cellStack[0]=cellStack[i-1]; 
   	cellStack.disp();
   	cellStack[0]=cellStack(i-2); 
   	cellPtr=cellStack(i-2);
   	cellStack.disp();
   	cout << "output cellptr" << endl;
   	cellPtr->disp();
   	cout << "Display each volu in stack" << endl;



   	cout << cellStack.capacity() << endl;

   	cout << "HashCellArrayStack" << endl;
   	cellStack.HashArray();
   	cellStack.SetMaxIndex();
   	cout << "Max Index in Array: " << cellStack.GetMaxIndex() <<  endl;

   	cellStack.Concatenate(cellStack);
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
   	cellStack.assign(3,singleCell); 

   	singleCell.index=1;
   	cout << "push_back in cellStack" << endl;
   	cellStack.push_back(singleCell); 

      // manipulate int stack 

   	cout << "Calls to display" << endl;
   	cout << "Display each volu in stack" << endl;

   	for ( i = 0; i < 4; ++i)
   	{
   		cellStack[i].index=i;
   		cellStack[i].edgeind[0]=i;
   		cellStack[i].disp() ;
   	}

   	cout << "Assign last volu to first volu" << endl;
   	cellStack[0]=cellStack[i-1]; 

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
   	cellStack.assign(3,singleCell); 

   	singleCell.index=1;
   	cout << "push_back in cellStack" << endl;
   	cellStack.push_back(singleCell); 

      // manipulate int stack 

   	cout << "Calls to display" << endl;
   	cout << "Display each volu in stack" << endl;

   	for ( i = 0; i < 4; ++i)
   	{
   		cellStack[i].index=i;
   		cellStack[i].surfind[0]=i;
   		cellStack[i].disp() ;
   	}

   	cout << "Assign last volu to first volu" << endl;
   	cellStack[0]=cellStack[i-1]; 

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
   	cellStack.assign(3,singleCell); 
   	singleCell.index=1;
   	cout << "push_back in cellStack" << endl;
   	cellStack.push_back(singleCell); 

      // manipulate int stack 

   	cout << "Calls to display" << endl;
   	cout << "Display each volu in stack" << endl;

   	for ( i = 0; i < 4; ++i)
   	{
   		cellStack[i].index=i;
   		cellStack[i].edgeind[0]=i;
   		cellStack[i].disp() ;
   	}

   	cout << "Assign last volu to first volu" << endl;
   	cellStack[0]=cellStack[i-1]; 

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

		mesh1.PopulateIndices();
		mesh2.PopulateIndices();
		mesh3=mesh1.MakeCompatible(mesh2);
    #ifdef TEST_ARRAYSTRUCT_MESH
	cout << "mesh1 " << endl;
	mesh1.disp() ;
	cout << "Compatible mesh #3  generated displaying ";
	cout << "mesh2 " << endl;
	mesh2.disp() ;
    #endif

	mesh1.MakeCompatible_inplace(mesh2);

	#ifdef TEST_ARRAYSTRUCT_MESH
	cout << "mesh2 made compatible inplace: ";
	cout << "mesh2 " << endl;

	cout << "mesh3 " << endl;
	mesh2.disp() ;
	mesh3.disp() ;
	#endif
	test1=CompareDisp(mesh2,mesh3);
	errFlag=errFlag || (!test1);
	cout << "Result of Comparison 2&3: " << test1 << endl;
	test1=CompareDisp(mesh1,mesh2);
	cout << "Result of Comparison 1&2: " << test1 << endl; 
	errFlag=errFlag || (test1);

	cout << "Result of Comparison Test_Volu&Test_Volu: " 
		<< CompareFuncOut(Test_Volu,Test_Volu) << endl;

	cout << "Concatenate mesh 1 and 2" << endl;
	mesh1.Concatenate(mesh2);
	mesh1.disp();


	} catch (exception const& ex) { 
		cerr << "Exception: " << ex.what() <<endl; 
		errFlag= -1;
	} 
	return(errFlag);
}

int Test_Crop(){
	mesh meshin;
	tecplotfile tecout;
	std::vector<double> lb, ub;

   try {
      meshin.read("../TESTOUT/mesh234.dat");
      tecout.OpenFile("../TESTOUT/meshcrop.plt");

      lb = {0.1,0.1,0.1};
      ub = {0.9,0.9,0.9};
      meshin.PrepareForUse();
      tecout.PrintMesh(meshin);
      tecout.PrintMesh(meshin,0,0,rsvs3d::constants::tecplot::line);
      auto verts = meshin.AddBoundary(lb, ub);
      tecout.PrintMesh(meshin,0,0,rsvs3d::constants::tecplot::line);
      tecout.PrintMesh(meshin);
      meshin.Crop(verts, 1);
      tecout.PrintMesh(meshin);
   } catch (exception const& ex) { 
      cerr << "Exception: " << ex.what() <<endl; 
      return(-1);
   } 


	return(0);
}
