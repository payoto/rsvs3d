#include <iostream>
#include <stdexcept>
#include <Eigen/Dense>

#include "voxel.hpp"
#include "arraystructures.hpp"
// Implementation of features

int BuildBlockGrid(RowVector3i dimGrid, mesh* blockGrid) {

   int nVolu, nVert, nSurf,nEdge;
   RowVector3i nSurfDim,nEdgeDim;
   Matrix3i surfProp, edgeProp;
   //int ii;

   // Calculate Number of Elements

   nVolu=dimGrid.prod();
   nVert=(dimGrid.array()+1).prod();

   surfProp=(MatrixXi::Identity(3,3)).rowwise().reverse();
   edgeProp=(1-MatrixXi::Identity(3,3).array());

   nSurfDim=((dimGrid.colwise().replicate(3).array())+surfProp.array()).rowwise().prod();
   nEdgeDim=((dimGrid.colwise().replicate(3).array())+edgeProp.array()).rowwise().prod();

   nSurf=nSurfDim.sum();
   nEdge=nEdgeDim.sum();

   blockGrid->Init(nVert,nEdge,nSurf,nVolu);

   #ifdef TEST_VOXEL
   cout << "nVert: " << nVert << " nVolu: " << nVolu 
      << " nEdge: " << nEdge << " nSurf: " << nSurf << endl; 
   cout << "surfProp: " << endl << surfProp << endl;
   cout << "edgeProp: " << endl << edgeProp << endl;
   cout << "nSurfDim: " << endl << nSurfDim << endl;
   cout << "nEdgeDim: " << endl << nEdgeDim << endl;
   cout << endl;
   
   #endif //TEST_VOXEL

   #ifdef TEST_EIGENEXT
   cout << "cumprod(nSurfDim): " << endl << cumprod(nSurfDim,0) << endl;
   cout << "nSurfDim: " << endl << nSurfDim << endl;
   cout << "cumsum(nSurfDim): " << endl << cumsum(nSurfDim,0) << endl;
   cout << "nSurfDim: " << endl << nSurfDim << endl;
   cout << "cumsum(surfProp) : " << endl << cumsum(surfProp,0) << endl;
   cout << "cumprod(edgeProp): " << endl << cumprod(edgeProp,0) << endl;
   cout << "edgeProp: " << endl << edgeProp << endl;
   cout << "cumsum(edgeProp')'' : " << endl << cumsum(edgeProp,1)<< endl;
   cout << "edgeProp: " << endl << edgeProp << endl;
   cout << "cumprod(edgeProp')': " << endl << cumprod(edgeProp,1) << endl;
   cout << "edgeProp: " << endl << edgeProp << endl;
   cout << endl << "------------------------" << endl << endl;
   #endif // TEST_EIGENEXT

   BuildBlockVolu(dimGrid, nVolu, blockGrid, nSurfDim, surfProp);

   return(0);
}

int BuildBlockVolu(RowVector3i dimGrid, int nVolu,  mesh* blockGrid,
	RowVector3i nSurfDim, Matrix3i surfProp){

	RowVector3i incrSurf, pos;
	Matrix3i matSurf;
	Matrix<int,3,3> incPos;
	Matrix<int,6,1> incrSurf2, tempSurfInd;
	int jj;

	matSurf=(dimGrid.colwise().replicate(3).array())+surfProp.array();

	incrSurf << 0 , nSurfDim.head<2>() ;
	incrSurf=cumsum(incrSurf,0);
	incrSurf2 << incrSurf.rowwise().replicate(2).transpose();

	incPos << Matrix<int,3,1>::Ones() ,  matSurf.leftCols<2>();
	incPos=cumprod(incPos,0);

	cout << "Size of volus: " << blockGrid->volus.elems.capacity() << endl;

	for (int ii=0;ii<nVolu;++ii){
		blockGrid->volus.elems[ii].index=ii+1;
		blockGrid->volus.elems[ii].fill=double(rand()%1000+1)/1000.0;

		pos(0)=ii%(dimGrid(0));
		pos(1)=int(floor(float(ii)/float(dimGrid(0))))%(dimGrid(1));
		pos(2)=int(floor(double(ii)/double(dimGrid(0)*dimGrid(1))))%dimGrid(2);

		tempSurfInd.head<3>()= pos.colwise().replicate(3).cwiseProduct(incPos).rowwise().sum();
		tempSurfInd.tail<3>()=(pos.colwise().replicate(3)+surfProp).cwiseProduct(incPos).rowwise().sum();
		tempSurfInd=(tempSurfInd+incrSurf2).array()+1;


		blockGrid->volus.elems[ii].surfind.assign(tempSurfInd.size(),0);
		for (jj=0;jj<tempSurfInd.size();jj++){
			blockGrid->volus.elems[ii].surfind[jj]=tempSurfInd(jj);

		}
		cout << endl;

	}

	#ifdef TEST_VOXEL
	cout << "matSurf: " << endl << matSurf << endl;
	cout << "incrSurf: " << endl << incrSurf << endl;
	cout << "incrSurf2: " << endl << incrSurf2 << endl;
	cout << "incPos: " << endl << incPos << endl;
	cout << "tempSurfInd: " << endl << tempSurfInd << endl;
	blockGrid->volus.disp();
	#endif //TEST_VOXEL

	return(0);

}


// Test code
int Test_BuildBlockGrid() { 
   // Test the functionality provided by arraystructures

   int errFlag,errTest;
   RowVector3i dimGrid(2,3,4);
   mesh blockGrid;

   errFlag=0;

   cout << "--------------------------------------------" << endl;
   cout << "      testing BuildBlockGrid" << endl;
   cout << "--------------------------------------------" << endl;
   errTest=BuildBlockGrid(dimGrid,&blockGrid);
   errFlag= errFlag | (errTest!=0);

   return(errFlag);
} 