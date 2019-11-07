/*#include <iostream>
#include <stdexcept>
#include <Eigen>

#include "voxel.hpp"
#include "arraystructures.hpp"
*/
//#pragma GCC diagnostic ignored "-Wignored-attributes" // Added because the precompiled voxel error throws this a lot

#include <iostream>
#include <numeric>      // std::partial_sum
#include <Eigen>
#include <ctime>

#include "mesh.hpp"
#include "voxel.hpp"
#include "postprocessing.hpp"

// Namespaces
using namespace std;

// Main code
int BuildBlockGrid(std::array<int, 3> &dimGrid, mesh& blockGrid){

	int out; 
	Eigen::RowVector3i dimGridEig;
		//  = Eigen::Map<Eigen::RowVector3i, 
		// Eigen::Unaligned>(dimGrid.data(), dimGrid.size());

	for (int i = 0; i < 3; ++i)
	{
		dimGridEig[i] = dimGrid[i];
	}

	out = BuildBlockGrid(dimGridEig,blockGrid);
	return(out);
}

// Implementation of Cartesian Block grids
int BuildBlockGrid(const Eigen::RowVector3i  &dimGrid, mesh& blockGrid) {

	int nVolu, nVert, nSurf,nEdge;
	Eigen::RowVector3i nSurfDim,nEdgeDim;
	Eigen::Matrix3i surfProp, edgeProp;
   //int ii;

   // Calculate Number of Elements

	nVolu=dimGrid.prod();
	nVert=(dimGrid.array()+1).prod();

	surfProp=(Eigen::MatrixXi::Identity(3,3)).rowwise().reverse();
	edgeProp=(1-Eigen::MatrixXi::Identity(3,3).array());

	nSurfDim=((dimGrid.colwise().replicate(3).array())
		+surfProp.array()).rowwise().prod();
	nEdgeDim=((dimGrid.colwise().replicate(3).array())
		+edgeProp.array()).rowwise().prod();

	nSurf=nSurfDim.sum();
	nEdge=nEdgeDim.sum();

	blockGrid.Init(nVert,nEdge,nSurf,nVolu);

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
		surfProp,  edgeProp,  nSurfDim,  nEdgeDim);
	BuildBlockEdge(dimGrid,  blockGrid,  nEdge , nEdgeDim,nSurfDim,
		edgeProp, surfProp);
	BuildBlockVert( dimGrid, blockGrid, nVert, edgeProp, nEdgeDim);
	
	blockGrid.PrepareForUse();
	blockGrid.TightenConnectivity();
	blockGrid.OrderEdges(); 
	blockGrid.OrientFaces();
	return(0);
}

int BuildBlockVolu(const Eigen::RowVector3i  &dimGrid, int nVolu,  mesh& blockGrid,
	const Eigen::RowVector3i  &nSurfDim, const Eigen::Matrix3i  &surfProp){

	Eigen::RowVector3i incrSurf, pos;
	Eigen::Matrix3i matSurf;
	Eigen::Matrix<int,3,3> incPos;
	Eigen::Matrix<int,6,1> incrSurf2, tempSurfInd;
	int jj;

	matSurf=(dimGrid.colwise().replicate(3).array())+surfProp.array();

	incrSurf << 0 , nSurfDim.head<2>() ;
	incrSurf=cumsum(incrSurf,0);
	incrSurf2 << incrSurf.rowwise().replicate(2).transpose();

	incPos << Eigen::Matrix<int,3,1>::Ones() ,  matSurf.leftCols<2>();
	incPos=cumprod(incPos,0);

	//cout << "Size of volus: " << blockGrid.volus.capacity() << endl;

	for (int ii=0;ii<nVolu;++ii){
		blockGrid.volus[ii].index=ii+1;
		blockGrid.volus[ii].fill=double(rand()%1000+1)/1000.0;

		pos(0)=ii%(dimGrid(0));
		pos(1)=int(floor(float(ii)/float(dimGrid(0))))%(dimGrid(1));
		pos(2)=int(floor(double(ii)/double(dimGrid(0)*dimGrid(1))))%dimGrid(2);

		tempSurfInd.head<3>()= pos.colwise().replicate(3).cwiseProduct(incPos).rowwise().sum();
		tempSurfInd.tail<3>()=(pos.colwise().replicate(3)+surfProp).cwiseProduct(incPos).rowwise().sum();
		tempSurfInd=(tempSurfInd+incrSurf2).array()+1;


		blockGrid.volus[ii].surfind.assign(tempSurfInd.size(),0);
		for (jj=0;jj<tempSurfInd.size();jj++){
			blockGrid.volus[ii].surfind[jj]=tempSurfInd(jj);

		}

	}

	#ifdef TEST_VOXEL_VOLU
	cout << "--------------------------------------------" << endl << "Test BuildBlockVolu" << endl;
	cout << "matSurf: " << endl << matSurf << endl;
	cout << "incrSurf: " << endl << incrSurf << endl;
	cout << "incrSurf2: " << endl << incrSurf2 << endl;
	cout << "incPos: " << endl << incPos << endl;
	cout << "tempSurfInd: " << endl << tempSurfInd << endl;
	blockGrid.volus.disp();
	#endif //TEST_VOXEL_VOLU

	return(0);
}

int BuildBlockSurf(const Eigen::RowVector3i  &dimGrid, int nSurf ,mesh& blockGrid ,
	const Eigen::Matrix3i  &surfProp, const Eigen::Matrix3i  &edgeProp, const Eigen::RowVector3i  &nSurfDim, const Eigen::RowVector3i  &nEdgeDim){
	/*Builds the surface information for the generation of block grids of size dimGrid*/

	Eigen::RowVector3i incrSurf, incrEdge, cumSurfDim,dimGridCur, pos, 
		cumProdDimGrid, ind;
	Eigen::Matrix<bool, 1,3> maskVec;
	Eigen::RowVector2i currDim;
	Eigen::Matrix3i matSurf, matEdge;
	Eigen::Matrix<int,3,3> incPos;
	Eigen::Matrix<int,6,1> incrSurf2;
	Eigen::Matrix<int,4,3> mask,incPosTemp;
	Eigen::Matrix<int,4,1> tempEdgeInd,incrEdge2;
	int  isFill,jPlane, boundaryFlag, jj, kk;

	// Calculate arrays which can be precalculated
	matSurf=(dimGrid.colwise().replicate(3).array())+surfProp.array();
	matEdge=(dimGrid.colwise().replicate(3).array())+edgeProp.array();

	incrSurf << 0 , nSurfDim.head<2>() ;
	incrSurf=cumsum(incrSurf,0);
	incrSurf2 << incrSurf.rowwise().replicate(2).transpose();


	incrEdge << 0 , nEdgeDim.head<2>() ;
	incrEdge=cumsum(incrEdge,0);
	

	incPos << Eigen::Matrix<int,3,1>::Ones() ,  matEdge.leftCols<2>();
	incPos=cumprod(incPos,0);

	cumSurfDim=cumsum(nSurfDim,0);

	cumProdDimGrid << 1 , dimGrid.head(2);
	cumProdDimGrid=cumprod(cumProdDimGrid,0);
	isFill=1-(dimGrid.all());
	ind << 0,1,2;
	// Assign content to Array
	for (int ii=0;ii<nSurf;++ii){

		blockGrid.surfs[ii].index=ii+1;
		blockGrid.surfs[ii].fill=isFill*double(rand()%1000+1)/1000.0;

		jPlane=(ii>=cumSurfDim.array()).cast<int>().sum();
		
		dimGridCur=matSurf.row(jPlane);

		pos(0)=(ii-incrSurf[jPlane])%(dimGridCur[0]);
		pos(1)=int(floor(double(ii-incrSurf[jPlane])/double(dimGridCur[0])))%dimGridCur[1];
		pos(2)=int(floor(double(ii-incrSurf[jPlane])/double(dimGridCur[0]*dimGridCur[1])))%dimGridCur[2];

		// define voluind
		boundaryFlag=!((pos.array()>=dimGrid.array()).any());
		blockGrid.surfs[ii].voluind[0]=boundaryFlag*((pos*cumProdDimGrid.transpose())+1);
		boundaryFlag=!(((pos-surfProp.row(jPlane)).array()<0).any());
		blockGrid.surfs[ii].voluind[1]=boundaryFlag*
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
		

		blockGrid.surfs[ii].edgeind.reserve(tempEdgeInd.size());
		blockGrid.surfs[ii].edgeind.assign(tempEdgeInd.size(),0);
		for (jj=0;jj<tempEdgeInd.size();jj++){
			blockGrid.surfs[ii].edgeind[jj]=tempEdgeInd[jj];
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
	blockGrid.surfs.disp();
	#endif //TEST_VOXEL_SURF


	return(0);
}

int BuildBlockEdge(const Eigen::RowVector3i  &dimGrid, mesh& blockGrid, 
	int nEdge ,const Eigen::RowVector3i  &nEdgeDim,
	const Eigen::RowVector3i  &nSurfDim, const Eigen::Matrix3i  &edgeProp,
	const Eigen::Matrix3i  &surfProp ){
	/*Function to build the Edges of a simple cartesian block grid*/

	Eigen::RowVector3i incrSurf, incrEdge, cumSurfDim,dimGridCur, pos, cumProdDimGrid, ind;
	Eigen::RowVector2i currDim;
	Eigen::Matrix3i matSurf, matEdge;
	Eigen::Matrix<bool, 1,3> maskVec;
	Eigen::Matrix<bool, 6,1> surfLog;
	Eigen::Matrix<int,3,3> incPos;
	Eigen::Matrix<int,6,1> incrSurf2,surfIndTemp;
	Eigen::Matrix<int,6,3> mask, surfLogTemp;
	Eigen::Matrix<int,4,3> incPosTemp;
	Eigen::Matrix<int,4,1> tempEdgeInd,incrEdge2;
	Eigen::Matrix<int,2,1> vertIndTemp,indTemp1,indTemp2;
	Eigen::Matrix<int,2,3> posMatTemp;

	int  jPlane, jj, kk, nSurfEdge;

	// Calculate arrays which can be precalculated
	matSurf=(dimGrid.colwise().replicate(3).array())+surfProp.array();
	matEdge=(dimGrid.colwise().replicate(3).array())+edgeProp.array();

	incrSurf << 0 , nSurfDim.head<2>() ;
	incrSurf=cumsum(incrSurf,0);
	incrSurf2 << incrSurf.rowwise().replicate(2).transpose();


	incrEdge << 0 , nEdgeDim.head<2>() ;
	incrEdge=cumsum(incrEdge,0);
	

	incPos << Eigen::Matrix<int,3,1>::Ones() ,  matSurf.leftCols<2>();
	incPos=cumprod(incPos,0);

	cumSurfDim=cumsum(nSurfDim,0);

	cumProdDimGrid << 1 , dimGrid.head(2).array()+1; // Warning ! Not the same as for Surf
	cumProdDimGrid=cumprod(cumProdDimGrid,0);

	
	ind << 0,1,2;

	
	for (int ii=0; ii<nEdge; ii++){

		blockGrid.edges[ii].index=ii+1;
		jPlane=(incrEdge.array()<=ii).cast<int>().sum()-1;
		dimGridCur=matEdge.row(jPlane);

		pos(0)=(ii-incrEdge[jPlane])%(dimGridCur[0]);
		pos(1)=int(floor(double(ii-incrEdge[jPlane])/double(dimGridCur[0])))%dimGridCur[1];
		pos(2)=int(floor(double(ii-incrEdge[jPlane])/double(dimGridCur[0]*dimGridCur[1])))%dimGridCur[2];

		//cout << ii << " | " << pos << " | " << jPlane << endl;
		// Assign vertind
		posMatTemp.row(0)=pos ;
		posMatTemp.row(1)=(pos.array()+1)-edgeProp.row(jPlane).array();
		vertIndTemp=(posMatTemp*cumProdDimGrid.transpose()).array()+1;

		blockGrid.edges[ii].vertind.reserve(vertIndTemp.size());
		blockGrid.edges[ii].vertind.assign(vertIndTemp.size(),0);
		for (jj=0;jj<vertIndTemp.size();jj++){
			blockGrid.edges[ii].vertind[jj]=vertIndTemp[jj];
		}

		// Assign surfind values
		mask=mask.setZero();
		maskVec=(ind.reverse().array()!=jPlane);
		kk=0;
		for(jj=0;jj<3;jj++){
			indTemp1(kk)=ind(jj);
			kk=kk+maskVec(jj);
			if (kk==2) {break;}
		}
		maskVec=(ind.array()!=jPlane);
		kk=0;
		for(jj=0;jj<3;jj++){
			indTemp2(kk)=ind(jj);
			kk=kk+maskVec(jj);
			if (kk==2) {break;}
		}
		for (jj=0;jj<2;jj++){
			mask(indTemp1(jj)+3,indTemp2(jj))=-1;
		}

		surfLogTemp=pos.colwise().replicate(6)+mask;

		surfIndTemp=(surfLogTemp.cwiseProduct(incPos.colwise().replicate(2)).rowwise().sum()
			+incrSurf2).array()+1;

		surfLog=!(((surfLogTemp.array()<0) || 
			(surfLogTemp.array()>(matSurf.colwise().replicate(2).array()-1))).rowwise().any()
		|| (ind.rowwise().replicate(2).reverse().array()==jPlane).transpose());
		nSurfEdge=max(surfLog.cast<int>().sum(),2);
		blockGrid.edges[ii].surfind.reserve(nSurfEdge);
		blockGrid.edges[ii].surfind.assign(nSurfEdge,0);
		kk=0;
		for (jj=0;jj<surfLog.size();jj++){
			if (surfLog(jj)){
				blockGrid.edges[ii].surfind[kk]=surfIndTemp[jj];
				kk++;
			}
		}
		//if (kk==1){blockGrid.edges[ii].surfind[kk]=0;} // not needed as already initialised
	}

	#ifdef TEST_VOXEL_EDGE
	cout << "--------------------------------------------" << endl << "Test BuildBlockEdge" << endl;
	cout << "matSurf: " << endl << matSurf << endl;
	cout << "incrSurf: " << endl << incrSurf << endl;
	cout << "incrEdge: " << endl << incrEdge << endl;
	cout << "incPos: " << endl << incPos << endl;
	cout << "vertIndTemp: " << endl << vertIndTemp << endl;
	blockGrid.edges.disp();
	#endif //TEST_VOXEL_EDGE

	return(0);
}

int BuildBlockVert(const Eigen::RowVector3i  &dimGrid, mesh& blockGrid, int nVert, 
	const Eigen::Matrix3i  &edgeProp, const Eigen::RowVector3i  &nEdgeDim){
	// Builds vertices for block grid
	Eigen::RowVector3i  cumSurfDim,dimGridVert,dimGridAct, pos, cumProdDimGrid,
		incrEdge;
	Eigen::RowVector3d coordTemp;
	Eigen::RowVector2i currDim;
	Eigen::Matrix3i  matEdge, incPos;
	Eigen::Matrix<int, 6,1> edgeIndTemp;
	Eigen::Matrix<bool, 6,1> edgeLog;
	Eigen::Matrix<int, 6,3> edgeLogTemp;


	int ii,jj,kk;

	incrEdge << 0 , nEdgeDim.head<2>() ;
	incrEdge=cumsum(incrEdge,0);

	dimGridVert=dimGrid.array()+1;
	dimGridAct=dimGrid.cwiseMax(1);

	matEdge=(dimGrid.colwise().replicate(3).array())+edgeProp.array();
	Eigen::Matrix3i edgeProp2 = 1-edgeProp.array();
	incPos << Eigen::Matrix<int,3,1>::Ones() ,  matEdge.leftCols<2>();
	incPos=cumprod(incPos,0);


	for(ii=0;ii<nVert;ii++){
		blockGrid.verts[ii].index=ii+1;

		pos(0)=ii%(dimGridVert[0]);
		pos(1)=int(floor(float(ii)/float(dimGridVert[0])))%(dimGridVert[1]);
		pos(2)=int(floor(double(ii)/double(dimGridVert[0]*dimGridVert[1])))%dimGridVert[2];

		coordTemp=pos.cast<double>().array()/dimGridAct.cast<double>().array();
		for (jj=0;jj<3;jj++){
			blockGrid.verts[ii].coord[jj]=coordTemp[jj];
		}

		edgeLogTemp << pos.colwise().replicate(3), (pos.colwise().replicate(3)-edgeProp2);
		edgeIndTemp=(edgeLogTemp.cwiseProduct(incPos.colwise().replicate(2))
			.rowwise().sum()+(incrEdge.rowwise().replicate(2).transpose())).array()+1;

		edgeLog=!(((edgeLogTemp.array()<0) 
			|| (edgeLogTemp.array()>(matEdge.colwise().replicate(2).array()-1))).rowwise().any());

		blockGrid.verts[ii].edgeind.assign(int(edgeLog.cast<int>().sum()),0);
		kk=0;

		for (jj=0;jj<6;jj++){
			if(edgeLog[jj]){
				blockGrid.verts[ii].edgeind[kk]=edgeIndTemp[jj];
				++kk;
			}
		}

	}

	#ifdef TEST_VOXEL_VERT
	cout << "--------------------------------------------" << endl << "Test BuildBlockVert" << endl;
	cout << "matEdge: " << endl << matEdge << endl;
	cout << "incPos: " << endl << incPos << endl;

	blockGrid.verts.disp();
	#endif
	return(0);
}

// Implementation of Refinement of block grids

// Test code
int Test_BuildBlockGrid_noout() { 
   // Test the functionality provided by arraystructures

	int errFlag,errTest;
	Eigen::RowVector3i dimGrid(2,3,4);
	mesh blockGrid;

	errFlag=0;

	errTest=BuildBlockGrid(dimGrid,blockGrid);
	errFlag= errFlag | (errTest!=0);

	return(errFlag);
} 



int Test_MeshOut(){
	int errFlag,errTest, start_s,stop_s;
	Eigen::RowVector3i dimGrid1(6,6,12), dimGrid2(100,100,0), 
		dimGrid3(20,30,10), dimGrid4(2,3,4);
	mesh blockGrid;
	const char *fileToOpen;
	
	tecplotfile outmesh1, outmesh3, outmesh2, outmesh4;

	errFlag=0;
	errTest=BuildBlockGrid(dimGrid1,blockGrid);
	errFlag+= (errTest!=0);
	fileToOpen="../TESTOUT/tecout6612.plt";
	errTest=outmesh1.OpenFile(fileToOpen);
	errFlag+= (errTest!=0);
	errTest=outmesh1.PrintMesh(blockGrid);
	errFlag+= (errTest!=0);
	fileToOpen="../TESTOUT/mesh6612.dat";
	errFlag+= TestCompareReadWrite(fileToOpen, blockGrid, outmesh1);

	errTest+=BuildBlockGrid(dimGrid2,blockGrid);
	errFlag+= (errTest!=0);
	fileToOpen="../TESTOUT/tecout100100.plt";
	errTest=outmesh2.OpenFile(fileToOpen);
	errFlag+= (errTest!=0);
	fileToOpen="../TESTOUT/tecout100100.dat";
	errFlag+= TestCompareReadWrite(fileToOpen, blockGrid, outmesh2);

	errTest=outmesh2.PrintMesh(blockGrid);
	errFlag+= (errTest!=0);
	//scanf("%i %i %i",&dimGrid3[0],&dimGrid3[1],&dimGrid3[2]);
	start_s=clock();
	errTest=BuildBlockGrid(dimGrid3,blockGrid);
	errFlag+= (errTest!=0);
	// the code you wish to time goes here
	stop_s=clock();
	cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << "ms" << endl;
	fileToOpen="../TESTOUT/tecout203010.plt";
	errTest=outmesh3.OpenFile(fileToOpen);
	errFlag+= (errTest!=0); 
	errTest=outmesh3.PrintMesh(blockGrid);
	errFlag+= (errTest!=0);

	fileToOpen="../TESTOUT/mesh203010.dat";
	errFlag+= TestCompareReadWrite(fileToOpen, blockGrid, outmesh3);

	return(errFlag);
}


int Test_MeshOut2(){
	int errFlag;
	Eigen::RowVector3i dimGrid1(6,6,12), dimGrid2(100,100,0), 
		dimGrid3(20,30,10), dimGrid4(2,3,4);


	errFlag=0;
	errFlag+=Test_MeshOut_Size(dimGrid1, false);
	errFlag+=Test_MeshOut_Size(dimGrid2, false);
	errFlag+=Test_MeshOut_Size(dimGrid3, false);
	errFlag+=Test_MeshOut_Size(dimGrid4, false);


	return(errFlag);
}


int Test_MeshOut_Size(Eigen::RowVector3i dimGrid, bool runReOrder){
	mesh blockGrid;
	std::string fileToOpen, gridSize;
	int errFlag, errTest, start_s;
	tecplotfile outmesh;

	gridSize = std::to_string(dimGrid[0]) + std::to_string(dimGrid[1])
		+ std::to_string(dimGrid[2]);

	if (runReOrder)
	{
		gridSize += "_reordered";
	}

	errFlag = 0;
	// Build grid
	start_s=clock();
	errTest=BuildBlockGrid(dimGrid,blockGrid);
	if (runReOrder)	{
		blockGrid.ReOrder();
	}
	errFlag+= (errTest!=0);
	rsvs3d::TimeStamp((std::string("Grid of size ") + gridSize + " : ").c_str(),
		start_s);

	// Print to PLT
	fileToOpen=std::string("../TESTOUT/tecout") + gridSize + ".plt";
	errTest=outmesh.OpenFile(fileToOpen.c_str());
	errFlag+= (errTest!=0); 
	errTest=outmesh.PrintMesh(blockGrid);
	errFlag+= (errTest!=0);

	// print to .msh file
	fileToOpen=std::string("../TESTOUT/mesh") + gridSize + ".dat";
	errFlag+= TestCompareReadWrite(fileToOpen.c_str(), blockGrid, outmesh);

	std::cout << std::endl;
	return (errFlag);
}

int Test_MeshReOrder(){
	int errFlag;
	Eigen::RowVector3i dimGrid1(6,6,12), dimGrid2(100,100,0), 
		dimGrid3(20,30,10), dimGrid4(2,3,4);


	errFlag=0;
	errFlag+=Test_MeshOut_Size(dimGrid1, true);
	errFlag+=Test_MeshOut_Size(dimGrid2, true);
	errFlag+=Test_MeshOut_Size(dimGrid3, true);
	errFlag+=Test_MeshOut_Size(dimGrid4, true);


	return(errFlag);
}
