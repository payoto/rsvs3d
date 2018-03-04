#include <iostream>
#include "arraystructures.hpp"
#include "postprocessing.hpp"


// Functions
void ExtractMeshData(mesh *grid,int *nVert, int *nVolu, int *nSurf, int *totNumFaceNode){
	// Extracts Data needed to write out a mesh to tecplot
	int ii;

	*nVert=grid->verts.size();
	*nVolu=grid->volus.size();
	*nSurf=grid->surfs.size();
	*totNumFaceNode=0;
	for (ii=0;ii<*nSurf;++ii){
		*totNumFaceNode=*totNumFaceNode+grid->surfs[ii]->edgeind.size();
	}
}

int tecplotfile::MeshDataBlock(mesh *meshout,int nVert,int nVolu, int nVertDat, int nCellDat){
	// Prints the Coord and Fill Data blocks to the tecplot file

	int ii,jj;
	// Print vertex Data
	for (jj=0;jj<int(meshout->verts[0]->coord.size());++jj){
		for ( ii = 0; ii<nVert; ++ii){
			this->Print("%.16lf ",meshout->verts[ii]->coord[jj]);
		}
		fprintf(fid,"\n");this->ResetLine();
	}
	for (jj=int(meshout->verts[0]->coord.size());jj<nVertDat;++jj){
		for ( ii = 0; ii<nVert; ++ii){
			this->Print("%lf ",0);
		}
		fprintf(fid,"\n");this->ResetLine();
	}
	// Print Cell Data
	for ( ii = 0; ii<nVolu; ++ii){
		this->Print("%.16lf ",meshout->volus[ii]->fill);
	}
	fprintf(fid,"\n");this->ResetLine();
	for ( ii = 0; ii<nVolu; ++ii){
		this->Print("%.16lf ",meshout->volus[ii]->target);
	}
	fprintf(fid,"\n");this->ResetLine();
	for ( ii = 0; ii<nVolu; ++ii){
		this->Print("%.16lf ",meshout->volus[ii]->error);
	}
	fprintf(fid,"\n");this->ResetLine();
	return(0);
}

int tecplotfile::MeshFaceMap(mesh *meshout,int nVert,int nSurf,int nVolu){
	int ii,jj,actVert,edgeCurr;
	int verts[2],vertsPast[2];
	for (ii=0;ii<nSurf;++ii){ // Print Number of vertices per face
		this->Print("%i ",meshout->surfs[ii]->edgeind.size());
	}
	fprintf(fid,"\n");this->ResetLine();

	for (ii=0;ii<nSurf;++ii){// print ordered  list of vertices in face
		jj=int(meshout->surfs[ii]->edgeind.size())-1;

		edgeCurr=meshout->edges.find(meshout->surfs[ii]->edgeind[jj]);
		verts[0]=meshout->verts.find(meshout->edges[edgeCurr]->vertind[0]);
		verts[1]=meshout->verts.find(meshout->edges[edgeCurr]->vertind[1]);
		vertsPast[0]=verts[0];
		vertsPast[1]=verts[1];

		for(jj=0;jj<int(meshout->surfs[ii]->edgeind.size());++jj){
			edgeCurr=meshout->edges.find(meshout->surfs[ii]->edgeind[jj]);

			verts[0]=meshout->verts.find(meshout->edges[edgeCurr]->vertind[0]);
			verts[1]=meshout->verts.find(meshout->edges[edgeCurr]->vertind[1]);

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
			actVert=meshout->surfs[ii]->voluind[jj];
			this->Print("%i ",meshout->volus.find(actVert)+1);
		}
		fprintf(fid,"\n");this->ResetLine();
	}

	return(0);
}

// Class function Implementation

int tecplotfile::PrintMesh(mesh *meshout){

	int nVert,nVolu,nSurf,totNumFaceNode,nVertDat,nCellDat;

	fprintf(fid, "VARIABLES = \"X\" ,\"Y\" , \"Z\" ,\"v1\" ,\"v2\", \"v3\"\n" );
	this->NewZone();
	ExtractMeshData(meshout,&nVert, &nVolu, &nSurf, &totNumFaceNode);
	// Fixed by the dimensionality of the mesh
	nVertDat=3;
	nCellDat=3;
	this->ZoneHeaderPolyHedron(nVert,nVolu,nSurf,totNumFaceNode,nVertDat,nCellDat);
	this->MeshDataBlock(meshout,nVert,nVolu, nVertDat, nCellDat);
	this->MeshFaceMap(meshout,nVert,nSurf,nVolu);
	return(0);
}


void tecplotfile::ZoneHeaderPolyHedron(int nVert, int nVolu, int nSurf, int totNumFaceNode,
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

int tecplotfile::OpenFile (const char *str){

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