#include <iostream>
#include "arraystructures.hpp"
#include "postprocessing.hpp"


// Functions
void ExtractMeshData(const mesh &grid,int *nVert, int *nEdge, int *nVolu, int *nSurf, int *totNumFaceNode){
	// Extracts Data needed to write out a mesh to tecplot
	int ii;

	*nVert=grid.verts.size();
	*nVolu=grid.volus.size();
	*nSurf=grid.surfs.size();
	*nEdge=grid.edges.size();
	*totNumFaceNode=0;
	for (ii=0;ii<*nSurf;++ii){
		*totNumFaceNode=*totNumFaceNode+grid.surfs(ii)->edgeind.size();
	}
}

int tecplotfile::VolDataBlock(const mesh& meshout,int nVert,int nVolu, int nVertDat){
	// Prints the Coord and Fill Data blocks to the tecplot file

	int ii,jj;
	// Print vertex Data
	for (jj=0;jj<int(meshout.verts(0)->coord.size());++jj){
		for ( ii = 0; ii<nVert; ++ii){
			this->Print("%.16lf ",meshout.verts(ii)->coord[jj]);
		}
		fprintf(fid,"\n");this->ResetLine();
	}
	for (jj=int(meshout.verts(0)->coord.size());jj<nVertDat;++jj){
		for ( ii = 0; ii<nVert; ++ii){
			this->Print("%lf ",0);
		}
		fprintf(fid,"\n");this->ResetLine();
	}
	// Print Cell Data
	for ( ii = 0; ii<nVolu; ++ii){
		this->Print("%.16lf ",meshout.volus(ii)->fill);
	}
	fprintf(fid,"\n");this->ResetLine();
	for ( ii = 0; ii<nVolu; ++ii){
		this->Print("%.16lf ",meshout.volus(ii)->target);
	}
	fprintf(fid,"\n");this->ResetLine();
	for ( ii = 0; ii<nVolu; ++ii){
		this->Print("%.16lf ",meshout.volus(ii)->error);
	}
	fprintf(fid,"\n");this->ResetLine();
	return(0);
}

int tecplotfile::SurfDataBlock(const mesh &meshout,int nVert,int nSurf, int nVertDat){
	// Prints the Coord and Fill Data blocks to the tecplot file

	int ii,jj;
	// Print vertex Data
	for (jj=0;jj<int(meshout.verts(0)->coord.size());++jj){
		for ( ii = 0; ii<nVert; ++ii){
			this->Print("%.16lf ",meshout.verts(ii)->coord[jj]);
		}
		fprintf(fid,"\n");this->ResetLine();
	}
	for (jj=int(meshout.verts(0)->coord.size());jj<nVertDat;++jj){
		for ( ii = 0; ii<nVert; ++ii){
			this->Print("%lf ",0);
		}
		fprintf(fid,"\n");this->ResetLine();
	}
	// Print Cell Data
	for ( ii = 0; ii<nSurf; ++ii){
		this->Print("%.16lf ",meshout.surfs(ii)->fill);
	}
	fprintf(fid,"\n");this->ResetLine();
	for ( ii = 0; ii<nSurf; ++ii){
		this->Print("%.16lf ",meshout.surfs(ii)->target);
	}
	fprintf(fid,"\n");this->ResetLine();
	for ( ii = 0; ii<nSurf; ++ii){
		this->Print("%.16lf ",meshout.surfs(ii)->error);
	}
	fprintf(fid,"\n");this->ResetLine();
	return(0);
}

int tecplotfile::LineDataBlock(const mesh &meshout,int nVert,int nEdge, int nVertDat,int nCellDat){
	// Prints the Coord and Fill Data blocks to the tecplot file

	int ii,jj;
	// Print vertex Data
	for (jj=0;jj<int(meshout.verts(0)->coord.size());++jj){
		for ( ii = 0; ii<nVert; ++ii){
			this->Print("%.16lf ",meshout.verts(ii)->coord[jj]);
		}
		fprintf(fid,"\n");this->ResetLine();
	}
	for (jj=int(meshout.verts(0)->coord.size());jj<nVertDat;++jj){
		for ( ii = 0; ii<nVert; ++ii){
			this->Print("%lf ",0);
		}
		fprintf(fid,"\n");this->ResetLine();
	}
	// Print Cell Data
	for(jj=0;jj<nCellDat;++jj){
		for ( ii = 0; ii<nEdge; ++ii){
			this->Print("%i ",0);
		}
		fprintf(fid,"\n");this->ResetLine();
	}
	return(0);
}
int tecplotfile::SurfFaceMap(const mesh &meshout,int nEdge){
	int ii,jj,actVert;
	int verts[2];

	for (ii=0;ii<nEdge;++ii){ // Print Number of vertices per face
		verts[0]=meshout.edges(ii)->vertind[0];
		verts[1]=meshout.edges(ii)->vertind[1];
		this->Print("%i %i\n",meshout.verts.find(verts[0])+1,meshout.verts.find(verts[1])+1);
		this->ResetLine();
	}
	


	for (jj=0;jj<2;++jj){// print index of left and right facing volumes
		for (ii=0;ii<nEdge;++ii){
			//cout << ii << "," << jj << endl ;
			actVert=meshout.edges(ii)->surfind[jj];
			this->Print("%i ",meshout.surfs.find(actVert)+1);
		}
		fprintf(fid,"\n");this->ResetLine();
	}

	return(0);
}

int tecplotfile::LineFaceMap(const mesh &meshout,int nEdge){
	int ii;
	int verts[2];

	for (ii=0;ii<nEdge;++ii){ // Print Number of vertices per face
		verts[0]=meshout.edges(ii)->vertind[0];
		verts[1]=meshout.edges(ii)->vertind[1];
		this->Print("%i %i\n",meshout.verts.find(verts[0])+1,meshout.verts.find(verts[1])+1);
		this->ResetLine();
	}

	return(0);
}

int tecplotfile::VolFaceMap(const mesh &meshout,int nSurf){
	int ii,jj,actVert,edgeCurr;
	int verts[2],vertsPast[2];
	for (ii=0;ii<nSurf;++ii){ // Print Number of vertices per face
		this->Print("%i ",meshout.surfs(ii)->edgeind.size());
	}
	fprintf(fid,"\n");this->ResetLine();

	for (ii=0;ii<nSurf;++ii){// print ordered  list of vertices in face
		jj=int(meshout.surfs(ii)->edgeind.size())-1;

		edgeCurr=meshout.edges.find(meshout.surfs(ii)->edgeind[jj]);
		verts[0]=meshout.verts.find(meshout.edges(edgeCurr)->vertind[0]);
		verts[1]=meshout.verts.find(meshout.edges(edgeCurr)->vertind[1]);
		vertsPast[0]=verts[0];
		vertsPast[1]=verts[1];

		for(jj=0;jj<int(meshout.surfs(ii)->edgeind.size());++jj){
			edgeCurr=meshout.edges.find(meshout.surfs(ii)->edgeind[jj]);

			verts[0]=meshout.verts.find(meshout.edges(edgeCurr)->vertind[0]);
			verts[1]=meshout.verts.find(meshout.edges(edgeCurr)->vertind[1]);

			if ((verts[0]==vertsPast[0]) || (verts[1]==vertsPast[0])){
				actVert=0;
			} 
			#ifdef TEST_POSTPROCESSING
			else if ((verts[0]==vertsPast[1]) || (verts[1]==vertsPast[1])) {
				actVert=1;
			}
			#endif //TEST_POSTPROCESSING
			else {
				actVert=1;
				#ifdef TEST_POSTPROCESSING
				cerr << "Warning: postprocessing.cpp:tecplotfile::VolFaceMap"<< endl ;
				cerr << "		Mesh Output failed in output of facemap data" << endl;
				cerr << "		Surface is not ordered " << endl;
				return(-1);
				#endif
			}

			this->Print("%i ",vertsPast[actVert]+1);
			vertsPast[0]=verts[0];
			vertsPast[1]=verts[1];

		}
		fprintf(fid,"\n");this->ResetLine();
		

	}
	for (jj=0;jj<2;++jj){// print index of left and right facing volumes
		for (ii=0;ii<nSurf;++ii){
			//cout << ii << "," << jj << endl ;
			actVert=meshout.surfs(ii)->voluind[jj];
			this->Print("%i ",meshout.volus.find(actVert)+1);
		}
		fprintf(fid,"\n");this->ResetLine();
	}

	return(0);
}
// Class function Implementation
int tecplotfile::PrintMesh(const mesh& meshout,int strandID, double timeStep, int forceOutType){

	int nVert,nEdge,nVolu,nSurf,totNumFaceNode,nVertDat,nCellDat;

	if(nZones==0){
		fprintf(fid, "VARIABLES = \"X\" ,\"Y\" , \"Z\" ,\"v1\" ,\"v2\", \"v3\"\n" );
	}

	this->NewZone();

	if(strandID>0){
		this->StrandTime(strandID, timeStep);
	}
	ExtractMeshData(meshout,&nVert,&nEdge,&nVolu, &nSurf, &totNumFaceNode);
	// Fixed by the dimensionality of the mesh
	nVertDat=3;
	nCellDat=3;

	if (forceOutType==0){
		if(nVolu>0){
			forceOutType=1; // output as volume data (FEPOLYHEDRON)
		} else if (nSurf>0){
			forceOutType=2;// output as Surface data (FEPOLYGON)
		} else {
			forceOutType=3; // output as line data (FELINESEG)
		}
	}


	if (forceOutType==1){
		this->ZoneHeaderPolyhedron(nVert,nVolu,nSurf,totNumFaceNode,nVertDat,nCellDat);
		this->VolDataBlock(meshout,nVert,nVolu, nVertDat);
		this->VolFaceMap(meshout,nSurf);
	}
	else if (forceOutType==2){
		this->ZoneHeaderPolygon(nVert, nEdge,nSurf,nVertDat,nCellDat);
		this->SurfDataBlock(meshout,nVert,nSurf, nVertDat);
		this->SurfFaceMap(meshout,nEdge);
	} else if (forceOutType==3){
		this->ZoneHeaderFelineseg(nVert, nEdge,nVertDat,nCellDat);
		this->LineDataBlock(meshout,nVert,nEdge, nVertDat,nCellDat);
		this->LineFaceMap(meshout,nEdge);
	}
	return(0);
}




void tecplotfile::ZoneHeaderPolyhedron(int nVert, int nVolu, int nSurf, int totNumFaceNode,
	int nVertDat, int nCellDat){

	fprintf(fid, "ZONETYPE = FEPOLYHEDRON\n");
	fprintf(fid, "NODES = %i\n", nVert);
	fprintf(fid, "ELEMENTS = %i\n",nVolu);
	this->Print( "FACES = %i\n",nSurf);
	fprintf(fid, "TOTALNUMFACENODES = %i\n",totNumFaceNode);
	fprintf(fid, "NUMCONNECTEDBOUNDARYFACES = 0\n");
	fprintf(fid, "TOTALNUMBOUNDARYCONNECTIONS = 0\n");
	fprintf(fid, "VARLOCATION=([%i-%i]=NODAL ,[%i-%i]=CELLCENTERED)\n",1,nVertDat, nVertDat+1,nVertDat+nCellDat);
	fprintf(fid, "DATAPACKING=BLOCK\n");

}

void tecplotfile::ZoneHeaderPolygon(int nVert,int nEdge,  int nSurf, int nVertDat, int nCellDat){

	fprintf(fid, "ZONETYPE = FEPOLYGON\n");
	fprintf(fid, "NODES = %i\n", nVert);
	fprintf(fid, "ELEMENTS = %i\n",nSurf);
	this->Print( "FACES = %i\n",nEdge);
	fprintf(fid, "TOTALNUMFACENODES = %i\n",2*nEdge);
	fprintf(fid, "NUMCONNECTEDBOUNDARYFACES = 0\n");
	fprintf(fid, "TOTALNUMBOUNDARYCONNECTIONS = 0\n");
	fprintf(fid, "VARLOCATION=([%i-%i]=NODAL ,[%i-%i]=CELLCENTERED)\n",1,nVertDat, nVertDat+1,nVertDat+nCellDat);
	fprintf(fid, "DATAPACKING=BLOCK\n");

}

void tecplotfile::ZoneHeaderFelineseg(int nVert,int nEdge,  int nVertDat, int nCellDat){

	fprintf(fid, "ZONETYPE = FELINESEG\n");
	fprintf(fid, "NODES = %i\n", nVert);
	fprintf(fid, "ELEMENTS = %i\n",nEdge);
	this->Print( "FACES = %i\n",nEdge);
	fprintf(fid, "VARLOCATION=([%i-%i]=NODAL ,[%i-%i]=CELLCENTERED)\n",1,nVertDat, nVertDat+1,nVertDat+nCellDat);
	fprintf(fid, "DATAPACKING=BLOCK\n");

}

int tecplotfile::OpenFile (const char *str){
	fclose(fid);
	fid=fopen(str,"w");
	if (fid==NULL){
		cout << "File '" << str << "' failed to open" << endl; 
		return(-1);
	}
	return(0);
}
void tecplotfile::CloseFile (){

	fclose(fid);
	
}

// Test Code

int Test_tecplotfile(){

	const char *fileToOpen;
	tecplotfile outmesh;
	int errFlag;

	fileToOpen="..\\TESTOUT\\tecout.plt";

	errFlag=outmesh.OpenFile(fileToOpen);
	if (errFlag!=0){
		return(errFlag);
	}



	return(0);
}

int TestCompareReadWrite(const char* fileToOpen, mesh &blockGrid, tecplotfile &outmesh1){
	
	int errFlag=0;
	int errTest=0;
	mesh blockGrid2;
	FILE *fid;

	fid=fopen(fileToOpen,"w");
	if(fid!=NULL){
		blockGrid.PrepareForUse();
		blockGrid.SetBorders();
		blockGrid.write(fid);
		fclose(fid);
		fid=fopen(fileToOpen,"r");
		if(fid!=NULL){
			blockGrid2.read(fid);
			fclose(fid);
			blockGrid.PrepareForUse();
			blockGrid2.PrepareForUse();

			outmesh1.PrintMesh(blockGrid2);
			errTest=CompareDisp(blockGrid,blockGrid2);
			if (!errTest){
				//blockGrid.disp();
				//blockGrid2.disp();
				cerr << "Error: Displays were not the same after write read" << endl;
				errFlag++;
			}
		} else {
			cout << "File for mesh out failed to open" << endl; 
		}
	} else {
		cout << "File for mesh out failed to open" << endl; 
	}
	return(errFlag);
}