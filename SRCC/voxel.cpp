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
   BuildBlockSurf( dimGrid,  nSurf , blockGrid ,
	 surfProp,  edgeProp,  nSurfDim,  nEdgeDim, nEdge);
   BuildBlockEdge(dimGrid,  blockGrid,  nEdge , nEdgeDim,nSurfDim,
     edgeProp, surfProp);
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

	}

	#ifdef TEST_VOXEL_VOLU
	cout << "--------------------------------------------" << endl << "Test BuildBlockVolu" << endl;
	cout << "matSurf: " << endl << matSurf << endl;
	cout << "incrSurf: " << endl << incrSurf << endl;
	cout << "incrSurf2: " << endl << incrSurf2 << endl;
	cout << "incPos: " << endl << incPos << endl;
	cout << "tempSurfInd: " << endl << tempSurfInd << endl;
	blockGrid->volus.disp();
	#endif //TEST_VOXEL

	return(0);

}

int BuildBlockSurf(RowVector3i dimGrid, int nSurf ,mesh* blockGrid ,
	Matrix3i surfProp, Matrix3i edgeProp, RowVector3i nSurfDim, RowVector3i nEdgeDim,
	 int nEdge){
	/*Builds the surface information for the generation of block grids of size dimGrid*/

	RowVector3i incrSurf, incrEdge, cumSurfDim,dimGridCur, pos, cumProdDimGrid, ind;
	Matrix<bool, 1,3> maskVec;
	RowVector2i currDim;
	Matrix3i matSurf, matEdge;
	Matrix<int,3,3> incPos;
	Matrix<int,6,1> incrSurf2;
	Matrix<int,4,3> mask,incPosTemp;
	Matrix<int,4,1> tempEdgeInd,incrEdge2;
	int  isFill,jPlane, boundaryFlag, jj, kk;

	// Calculate arrays which can be precalculated
	matSurf=(dimGrid.colwise().replicate(3).array())+surfProp.array();
	matEdge=(dimGrid.colwise().replicate(3).array())+edgeProp.array();

	incrSurf << 0 , nSurfDim.head<2>() ;
	incrSurf=cumsum(incrSurf,0);
	incrSurf2 << incrSurf.rowwise().replicate(2).transpose();


	incrEdge << 0 , nEdgeDim.head<2>() ;
	incrEdge=cumsum(incrEdge,0);
	

	incPos << Matrix<int,3,1>::Ones() ,  matEdge.leftCols<2>();
	incPos=cumprod(incPos,0);

	cumSurfDim=cumsum(nSurfDim,0);

	cumProdDimGrid << 1 , dimGrid.head(2);
	cumProdDimGrid=cumprod(cumProdDimGrid,0);
	isFill=1-(dimGrid.all());
	ind << 0,1,2;
	// Assign content to Array
	for (int ii=0;ii<nSurf;++ii){

		blockGrid->surfs.elems[ii].index=ii+1;
		blockGrid->surfs.elems[ii].fill=isFill*double(rand()%1000+1)/1000.0;

		/*jPlane=0;
		for (jj=0;jj<3;++jj){
			jPlane=jPlane+(ii>=cumSurfDim[jj]);
		}*/
		jPlane=(ii>=cumSurfDim.array()).cast<int>().sum();
		

		dimGridCur=matSurf.row(jPlane);

		pos(0)=(ii-incrSurf[jPlane])%(dimGridCur[0]);
		pos(1)=int(floor(double(ii-incrSurf[jPlane])/double(dimGridCur[0])))%dimGridCur[1];
		pos(2)=int(floor(double(ii-incrSurf[jPlane])/double(dimGridCur[0]*dimGridCur[1])))%dimGridCur[2];

		// define voluind
		boundaryFlag=!((pos.array()>=dimGrid.array()).any());
		blockGrid->surfs.elems[ii].voluind[0]=boundaryFlag*((pos*cumProdDimGrid.transpose())+1);
		boundaryFlag=!(((pos-surfProp.row(jPlane)).array()<0).any());
		blockGrid->surfs.elems[ii].voluind[1]=boundaryFlag*
			(((pos-surfProp.row(jPlane))*cumProdDimGrid.transpose())+1);

		//define edgeind
		maskVec=(ind.reverse().array()!=jPlane);
		kk=0;
		for(jj=0;jj<3;jj++){
			currDim(kk)=ind(jj);
			kk=kk+maskVec(jj);
			if (kk==2) {break;}
		}

		mask=mask.setZero();

		for (jj=0;jj<2;jj++){ 
			mask(jj+1,currDim(jj))=1;// assign identity in mask
			incPosTemp.row(jj)=incPos.row(currDim(jj));// extract needed parts of incPos into Temp
			incPosTemp.row(jj+2)=incPos.row(currDim(jj));
			incrEdge2(jj)=incrEdge(currDim(jj));
			incrEdge2(jj+2)=incrEdge(currDim(jj));
		}
		tempEdgeInd=((pos.colwise().replicate<4>()+mask).cwiseProduct(incPosTemp)
			.rowwise().sum()+incrEdge2).array()+1;
		

		blockGrid->surfs.elems[ii].edgeind.reserve(tempEdgeInd.size());
		blockGrid->surfs.elems[ii].edgeind.assign(tempEdgeInd.size(),0);
		for (jj=0;jj<tempEdgeInd.size();jj++){
			blockGrid->surfs.elems[ii].edgeind[jj]=tempEdgeInd[jj];
		}
	}

	#ifdef TEST_VOXEL_SURF
	cout << "--------------------------------------------" << endl << "Test BuildBlockSurf" << endl;
	cout << "matSurf: " << endl << matSurf << endl;
	cout << "incrSurf: " << endl << incrSurf << endl;
	cout << "incrSurf2: " << endl << incrSurf2 << endl;
	cout << "incrEdge: " << endl << incrEdge << endl;
	cout << "incrEdge2: " << endl << incrEdge2 << endl;
	cout << "incPos: " << endl << incPos << endl;
	cout << "tempEdgeInd: " << endl << tempEdgeInd << endl;
	blockGrid->surfs.disp();
	#endif //TEST_VOXEL_SURF


	return(0);
}


int BuildBlockEdge(RowVector3i dimGrid, mesh* blockGrid, int nEdge ,RowVector3i nEdgeDim,
  RowVector3i nSurfDim, Matrix3i edgeProp,Matrix3i surfProp ){
	/*Function to build the Edges of a simple cartesian block grid*/

	RowVector3i incrSurf, incrEdge, cumSurfDim,dimGridCur, pos, cumProdDimGrid, ind;
	RowVector2i currDim;
	Matrix3i matSurf, matEdge;
	Matrix<bool, 1,3> maskVec;
	Matrix<int,3,3> incPos;
	Matrix<int,6,1> incrSurf2;
	Matrix<int,4,3> mask,incPosTemp;
	Matrix<int,4,1> tempEdgeInd,incrEdge2;
	Matrix<int,2,1> vertIndTemp;
	Matrix<int,2,3> posMatTemp;

	int  isFill,jPlane, boundaryFlag, jj, kk;

	// Calculate arrays which can be precalculated
	matSurf=(dimGrid.colwise().replicate(3).array())+surfProp.array();
	matEdge=(dimGrid.colwise().replicate(3).array())+edgeProp.array();

	incrSurf << 0 , nSurfDim.head<2>() ;
	incrSurf=cumsum(incrSurf,0);
	incrSurf2 << incrSurf.rowwise().replicate(2).transpose();


	incrEdge << 0 , nEdgeDim.head<2>() ;
	incrEdge=cumsum(incrEdge,0);
	

	incPos << Matrix<int,3,1>::Ones() ,  matSurf.leftCols<2>();
	incPos=cumprod(incPos,0);

	cumSurfDim=cumsum(nSurfDim,0);

	cumProdDimGrid << 1 , dimGrid.head(2).array()+1;
	cumProdDimGrid=cumprod(cumProdDimGrid,0);
	isFill=1-(dimGrid.all());
	ind << 0,1,2;

	
	for (int ii=0; ii<nEdge; ii++){
		blockGrid->surfs.elems[ii].index=ii+1;
		jPlane=(incrEdge.array()<ii).cast<int>().sum();
		dimGridCur=matEdge.row(jPlane);

		pos(0)=(ii-incrEdge[jPlane])%(dimGridCur[0]);
		pos(1)=int(floor(double(ii-incrEdge[jPlane])/double(dimGridCur[0])))%dimGridCur[1];
		pos(2)=int(floor(double(ii-incrEdge[jPlane])/double(dimGridCur[0]*dimGridCur[1])))%dimGridCur[2];

		posMatTemp << pos , (pos.array()+1)-edgeProp.row(jPlane).array();
		vertIndTemp=(posMatTemp*cumProdDimGrid.transpose()).array()+1;
	}

	#ifdef TEST_VOXEL_EDGE
	cout << "--------------------------------------------" << endl << "Test BuildBlockSurf" << endl;
	cout << "matSurf: " << endl << matSurf << endl;
	cout << "incrSurf: " << endl << incrSurf << endl;
	cout << "incrSurf2: " << endl << incrSurf2 << endl;
	cout << "incrEdge: " << endl << incrEdge << endl;
	cout << "incrEdge2: " << endl << incrEdge2 << endl;
	cout << "incPos: " << endl << incPos << endl;
	cout << "tempEdgeInd: " << endl << tempEdgeInd << endl;
	blockGrid->surfs.disp();
	#endif //TEST_VOXEL_EDGE

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