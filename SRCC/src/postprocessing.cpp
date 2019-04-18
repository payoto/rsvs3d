#include <iostream>

#include "mesh.hpp"
#include "postprocessing.hpp" 
#include "triangulate.hpp"
#include "warning.hpp"
#include "RSVScalc.hpp"

using namespace rsvs3d::constants;

// Functions
void ExtractMeshData(const mesh &grid,int *nVert, int *nEdge,
	int *nVolu, int *nSurf, int *totNumFaceNode){
	// Extracts Data needed to write out a mesh to tecplot
	int ii;

	*nVert=grid.verts.size();
	*nVolu=grid.volus.size();
	*nSurf=grid.surfs.size();
	*nEdge=grid.edges.size();
	*totNumFaceNode=0;
	for (ii=0;ii<*nSurf;++ii){
		*totNumFaceNode=*totNumFaceNode+int(grid.surfs(ii)->edgeind.size());
	}
}

int tecplotfile::VolDataBlock(const mesh& meshout,int nVert,int nVolu,
		int nVertDat, const std::vector<int> &voluList,
		const std::vector<int> &vertList){
	// Prints the Coord and Fill Data blocks to the tecplot file

	int ii,jj,nCoord;
	// Print vertex Data
	nCoord=int(meshout.verts(0)->coord.size());
	nCoord= nCoord > nVertDat ? nVertDat : nCoord;
	if(vertList.size()==0){
		for (jj=0;jj<nCoord;++jj){
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
	} else {
		for (jj=0;jj<nCoord;++jj){
			for(auto ii : vertList){

				this->Print("%.16lf ",meshout.verts.isearch(ii)->coord[jj]);
			}
			fprintf(fid,"\n");this->ResetLine();
		}
		nVert = vertList.size();
		for (jj=int(meshout.verts(0)->coord.size());jj<nVertDat;++jj){

			for(ii = 0; ii<nVert; ++ii){
				this->Print("%lf ",0);
			}
			fprintf(fid,"\n");this->ResetLine();
		}
	}
 	// Print Cell Data
	if(voluList.size()==0){
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
	} else {
		for(auto ii : voluList){
			this->Print("%.16lf ",meshout.volus.isearch(ii)->fill);
		}
		fprintf(fid,"\n");this->ResetLine();
		
		for(auto ii : voluList){
			this->Print("%.16lf ",meshout.volus.isearch(ii)->target);
		}
		fprintf(fid,"\n");this->ResetLine();
		for(auto ii : voluList){
			this->Print("%.16lf ",meshout.volus.isearch(ii)->error);
		}
		fprintf(fid,"\n");this->ResetLine();
	}
	return(0);
}

int tecplotfile::SnakeDataBlock(const snake& snakeout,int nVert, int nVertDat){
	// Prints the Coord and Fill Data blocks to the tecplot file

	int ii,jj,nCoord;
	// Print vertex Data
	nCoord=int(snakeout.snakeconn.verts(0)->coord.size());
	nCoord= nCoord > nVertDat ? nVertDat : nCoord;
	for (jj=0;jj<nCoord;++jj){
		for ( ii = 0; ii<nVert; ++ii){
			this->Print("%.16lf ",snakeout.snakeconn.verts(ii)->coord[jj]);
		}
		fprintf(fid,"\n");this->ResetLine();
	}
	for (jj=int(snakeout.snakeconn.verts(0)->coord.size());jj<nVertDat;++jj){
		for ( ii = 0; ii<nVert; ++ii){
			this->Print("%lf ",0);
		}
		fprintf(fid,"\n");this->ResetLine();
	}
	// Print Cell Data
	for ( ii = 0; ii<nVert; ++ii){
		this->Print("%.16lf ",snakeout.snaxs(ii)->d);
	}
	fprintf(fid,"\n");this->ResetLine();
	for ( ii = 0; ii<nVert; ++ii){
		this->Print("%.16lf ",snakeout.snaxs(ii)->v);
	}
	fprintf(fid,"\n");this->ResetLine();
	for ( ii = 0; ii<nVert; ++ii){
		auto temp = snakeout.snakeconn.verts(ii)->elmind(snakeout.snakeconn,2);
		int maxSize = 0;
		for (auto surfSub : temp){
			int tempSize = snakeout.snakeconn.surfs.isearch(surfSub)->edgeind.size();
			maxSize = maxSize < tempSize ? tempSize : maxSize;
		}
		// this->Print("%i ",snakeout.snaxs(ii)->isfreeze);
		this->Print("%i ", maxSize);
	}
	fprintf(fid,"\n");this->ResetLine();
	return(0);
}

int tecplotfile::SurfDataBlock(const mesh &meshout,int nVert,int nSurf, 
	int nVertDat){
	// Prints the Coord and Fill Data blocks to the tecplot file

	int ii,jj, nCoord;
	// Print vertex Data
	nCoord=int(meshout.verts(0)->coord.size());
	nCoord= nCoord > nVertDat ? nVertDat : nCoord;
	for (jj=0;jj<nCoord;++jj){
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

int tecplotfile::LineDataBlock(const mesh &meshout,int nVert,int nEdge, 
	int nVertDat,int nCellDat){
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

int tecplotfile::VertDataBlock(const mesh &meshout,int nVert, int nVertDat,
	int nCellDat, const vector<int> &vertList){
	// Prints the Coord and Fill Data blocks to the tecplot file

	int ii,jj, nCoord;
	nCoord=int(meshout.verts(0)->coord.size());
	// Print vertex Data
	if(vertList.size()>0){ // if vertList is a list of index
		for (jj=0;jj<nCoord;++jj){
			for ( ii = 0; ii<nVert; ++ii){
				this->Print("%.16lf ",meshout.verts.isearch(vertList[ii])->coord[jj]);
			}
			fprintf(fid,"\n");this->ResetLine();
		}
	} else if(int(vertList.size())==meshout.verts.size()) { // vertList is a boolean
		for (jj=0;jj<nCoord;++jj){
			for ( ii = 0; ii<nVert; ++ii){
				if(vertList[ii]){
					this->Print("%.16lf ",meshout.verts(ii)->coord[jj]);
				}
			}
			fprintf(fid,"\n");this->ResetLine();
		}
	} else {
		for (jj=0;jj<nCoord;++jj){
			for ( ii = 0; ii<nVert; ++ii){
				this->Print("%.16lf ",meshout.verts(ii)->coord[jj]);
			}
			fprintf(fid,"\n");this->ResetLine();
		}
	}
	for (jj=nCoord;jj<nVertDat;++jj){
		for ( ii = 0; ii<nVert; ++ii){
			this->Print("%lf ",0);
		}
		fprintf(fid,"\n");this->ResetLine();
	}
	// Print Cell Data
	for(jj=0;jj<nCellDat;++jj){
		for ( ii = 0; ii<nVert; ++ii){
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
				//meshout.surfs(ii)->disptree(meshout,4);
				
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
int tecplotfile::VolFaceMap(const mesh &meshout, 
	const std::vector<int> &surfList, const std::vector<int> &voluList,
	const std::vector<int> &vertList){
	int jj,actVert,edgeCurr;
	int verts[2],vertsPast[2];
	for (auto ii : meshout.surfs.find_list(surfList)){ // Print Number of vertices per face
		this->Print("%i ",meshout.surfs(ii)->edgeind.size());
	}
	fprintf(fid,"\n");this->ResetLine();
	HashedVector<int, int> verthash;
	verthash.vec=vertList;
	verthash.GenerateHash();
	for (auto ii : meshout.surfs.find_list(surfList))
	{// print ordered  list of vertices in face
		jj=int(meshout.surfs(ii)->edgeind.size())-1;

		edgeCurr=meshout.edges.find(meshout.surfs(ii)->edgeind[jj]);
		verts[0]=verthash.find(meshout.edges(edgeCurr)->vertind[0]);
		verts[1]=verthash.find(meshout.edges(edgeCurr)->vertind[1]);
		vertsPast[0]=verts[0];
		vertsPast[1]=verts[1];

		for(jj=0;jj<int(meshout.surfs(ii)->edgeind.size());++jj){
			edgeCurr=meshout.edges.find(meshout.surfs(ii)->edgeind[jj]);

			verts[0]=verthash.find(meshout.edges(edgeCurr)->vertind[0]);
			verts[1]=verthash.find(meshout.edges(edgeCurr)->vertind[1]);

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
				//meshout.surfs(ii)->disptree(meshout,4);
				
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
		for (auto ii : meshout.surfs.find_list(surfList)){
			//cout << ii << "," << jj << endl ;
			actVert=meshout.surfs(ii)->voluind[jj];
			int k=0;
			bool notFoundCell=true;
			for(auto j : voluList){
				k++;
				if(j==actVert){notFoundCell=false;break;}
			}
			if (notFoundCell){k=0;}
			this->Print("%i ",k);
		}
		fprintf(fid,"\n");this->ResetLine();
	}

	return(0);
}
// Class function Implementation
int tecplotfile::PrintMesh(const mesh& meshout,int strandID, double timeStep,
	int forceOutType, const vector<int> &vertList){

	int nVert,nEdge,nVolu,nSurf,totNumFaceNode,nVertDat,nCellDat;

	ExtractMeshData(meshout,&nVert,&nEdge,&nVolu, &nSurf, &totNumFaceNode);
	if(nVert==0){ // Don't print mesh
		return(1);
	}
	if(nZones==0){
		fprintf(fid, "VARIABLES = \"X\" ,\"Y\" , \"Z\" ,\"v1\" "
			",\"v2\", \"v3\"\n" );
	}

	this->NewZone();

	if(strandID>0){
		this->StrandTime(strandID, timeStep);
	}
	// Fixed by the dimensionality of the mesh
	nVertDat=3;
	nCellDat=3;

	if (forceOutType==tecplot::autoselect){
		if(nVolu>0){
			forceOutType=tecplot::polyhedron; // output as volume data (FEPOLYHEDRON)
		} else if (nSurf>0){
			forceOutType=tecplot::polygon;// output as Surface data (FEPOLYGON)
		} else if (nEdge>0){
			forceOutType=tecplot::line; // output as line data (FELINESEG)
		} else {
			forceOutType=tecplot::point;
		}
	}


	if (forceOutType==tecplot::polyhedron){
		this->ZoneHeaderPolyhedron(nVert,nVolu,nSurf,totNumFaceNode,
			nVertDat,nCellDat);
		this->VolDataBlock(meshout,nVert,nVolu, nVertDat);
		this->VolFaceMap(meshout,nSurf);
	}
	else if (forceOutType==tecplot::polygon){
		this->ZoneHeaderPolygon(nVert, nEdge,nSurf,nVertDat,nCellDat);
		this->SurfDataBlock(meshout,nVert,nSurf, nVertDat);
		this->SurfFaceMap(meshout,nEdge);
	} else if (forceOutType==tecplot::line){
		this->ZoneHeaderFelineseg(nVert, nEdge,nVertDat,nCellDat);
		this->LineDataBlock(meshout,nVert,nEdge, nVertDat,nCellDat);
		this->LineFaceMap(meshout,nEdge);
	} else if (forceOutType==tecplot::point){
		if(int(vertList.size())==nVert){
			nVert=0;
			for (int ii=0; ii< int(vertList.size());++ii){
				nVert += int(vertList[ii]); 
			}
		} else if (vertList.size()>0){
			nVert=int(vertList.size());
		}
		this->ZoneHeaderOrdered(nVert,nVertDat,nCellDat);
		this->VertDataBlock(meshout,nVert, nVertDat,nCellDat,vertList);
		// No map just points
	} else if (forceOutType==5){
		if (vertList.size()>0){
			nVolu=int(vertList.size());
		} else {
			RSVS3D_ERROR_ARGUMENT("voluList with forceOutType == 5 must at"
				" least be length 1");
		}
		nSurf=0; totNumFaceNode=0;
		auto voluSubList = meshout.volus.find_list(vertList);
		auto surfList = ConcatenateVectorField(meshout.volus,
			&volu::surfind, voluSubList);
		sort(surfList);
		unique(surfList);
		nSurf = surfList.size();
		auto surfSubList = meshout.surfs.find_list(surfList);
		auto edgeList = ConcatenateVectorField(meshout.surfs,
			&surf::edgeind, surfSubList);
		auto edgeSubList = meshout.edges.find_list(edgeList);
		auto actualVert = ConcatenateVectorField(meshout.edges,
			&edge::vertind, edgeSubList);
		sort(actualVert);
		unique(actualVert);
		nVert = actualVert.size();

		for(auto j : surfList){
			totNumFaceNode += meshout.surfs.isearch(j)->edgeind.size();
		}

		this->ZoneHeaderPolyhedron(nVert,nVolu,nSurf,totNumFaceNode,
			nVertDat,nCellDat);
		this->VolDataBlock(meshout,nVert,nVolu, nVertDat,vertList, actualVert);
		this->VolFaceMap(meshout, surfList,vertList, actualVert);
	}
	return(0);
}

// Class function Implementation
int tecplotfile::PrintSnake(const snake& snakeout,int strandID, double timeStep, 
	int forceOutType, const vector<int> &vertList){

	int nVert,nEdge,nVolu,nSurf,totNumFaceNode,nVertDat,nCellDat;

	ExtractMeshData(snakeout.snakeconn,&nVert,&nEdge,&nVolu, &nSurf, &totNumFaceNode);
	if(nVert==0){ // Don't print mesh
		return(1);
	}
	
	if(nZones==0){
		fprintf(fid, "VARIABLES = \"X\" ,\"Y\" , \"Z\" ,\"v1\" ,\"v2\", \"v3\"\n" );
	}

	this->NewZone();

	if(strandID>0){
		this->StrandTime(strandID, timeStep);
	}
	// Fixed by the dimensionality of the mesh
	nVertDat=3;
	nCellDat=3;

	if (forceOutType==tecplot::autoselect){
		if(nVolu>0){
			forceOutType=tecplot::polyhedron; // output as volume data (FEPOLYHEDRON)
		} else if (nSurf>0){
			forceOutType=tecplot::polygon;// output as Surface data (FEPOLYGON)
		} else if (nEdge>0){
			forceOutType=tecplot::line; // output as line data (FELINESEG)
		} else {
			forceOutType=tecplot::point;
		}
	}


	if (forceOutType==tecplot::polyhedron){
		this->ZoneHeaderPolyhedronSnake(nVert,nVolu,nSurf,totNumFaceNode,
			nVertDat,nCellDat);
		this->SnakeDataBlock(snakeout,nVert, nVertDat);
		this->VolFaceMap(snakeout.snakeconn,nSurf);
	}
	else if (forceOutType==tecplot::polygon){
		this->ZoneHeaderPolygonSnake(nVert, nEdge,nSurf,nVertDat,nCellDat);
		this->SnakeDataBlock(snakeout,nVert, nVertDat);
		this->SurfFaceMap(snakeout.snakeconn,nEdge);
	} else if (forceOutType==tecplot::line){
		this->ZoneHeaderFelinesegSnake(nVert, nEdge,nVertDat,nCellDat);
		this->SnakeDataBlock(snakeout,nVert, nVertDat);
		this->LineFaceMap(snakeout.snakeconn,nEdge);
	} else if (forceOutType==tecplot::point){
		if(int(vertList.size())==nVert){
			nVert=0;
			for (int ii=0; ii< int(vertList.size());++ii){
				nVert += int(vertList[ii]); 
			}
		} else if (vertList.size()>0){
			nVert=int(vertList.size());
		}
		this->ZoneHeaderOrdered(nVert,nVertDat,nCellDat);
		this->SnakeDataBlock(snakeout,nVert, nVertDat);
		// No map just points
	}
	return(0);
}


int tecplotfile::PrintVolumeDat(const mesh &meshout, int shareZone, 
	int strandID, double timeStep){

	int nVert,nEdge,nVolu,nSurf,totNumFaceNode,nVertDat,nCellDat;
	int forceOutType;
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
	// forceOutType=tecplot::autoselect;
	// if (forceOutType==tecplot::autoselect){
	if(nVolu>0){
		forceOutType=tecplot::polyhedron; // output as volume data (FEPOLYHEDRON)
	} else if (nSurf>0){
		forceOutType=tecplot::polygon;// output as Surface data (FEPOLYGON)
	} else if (nEdge>0){
		forceOutType=tecplot::line; // output as line data (FELINESEG)
	} else {
		forceOutType=tecplot::point;
	}
	// }


	if (forceOutType==tecplot::polyhedron){
		this->ZoneHeaderPolyhedron(nVert,nVolu,nSurf,totNumFaceNode,nVertDat,nCellDat);
		this->DefShareZoneVolume(shareZone, nVertDat);
		this->VolDataBlock(meshout,nVert,nVolu, 0);
		//this->VolFaceMap(meshout,nSurf);
	}
	else if (forceOutType==tecplot::polygon){
		this->ZoneHeaderPolygon(nVert, nEdge,nSurf,nVertDat,nCellDat);
		this->DefShareZoneVolume(shareZone, nVertDat);
		this->SurfDataBlock(meshout,nVert,nSurf, 0);
		//this->SurfFaceMap(meshout,nEdge);

	} else if (forceOutType==tecplot::line){
		RSVS3D_ERROR_ARGUMENT("Cannot output volume of line");
	} else if (forceOutType==tecplot::point){
		RSVS3D_ERROR_ARGUMENT("Cannot output volume of point");
	}
	return(0);
	//CONNECTIVITYSHAREZONE=<zone>
	//VARSHARELIST=([set-of-vars]=<zone>, [set-ofvars]=<zone>)

}
int tecplotfile::DefShareZoneVolume(int shareZone, int nVertDat){

	fprintf(fid, "CONNECTIVITYSHAREZONE=%i\n",shareZone);
	fprintf(fid, "VARSHARELIST=([%i-%i]=%i )\n",1, nVertDat, shareZone);

	return(0);
}

int tecplotfile::PrintSnakeInternalPts(const snake &snakein,int strandID, 
	double timeStep){
	int jj;
	std::vector<int> vertList;
	vertList.clear();
		for(jj=0;jj<int(snakein.isMeshVertIn.size()); ++jj){
			if(snakein.isMeshVertIn[jj]){
				vertList.push_back(snakein.snakemesh->verts(jj)->index);
			}
		}
		if(int(snakein.isMeshVertIn.size())==0){
			vertList.push_back(snakein.snakemesh->verts(0)->index);
		}
		jj=PrintMesh(*(snakein.snakemesh),strandID,timeStep,4,vertList);
		return(jj);
}

// Functions triarray
void ExtractTriData(const triangulation &triout, triarray triangulation::*mp,
	int *nVert, int *nEdge, int *nVolu, int *nSurf, int *totNumFaceNode){
	// Extracts Data needed to write out a mesh to tecplot


	*nVert=(triout.*mp).size()*3;
	*nVolu=1;
	*nSurf=(triout.*mp).size();
	*nEdge=(triout.*mp).size()*3;
	*totNumFaceNode=(triout.*mp).size()*3;
}

void ExtractTriData(int nTriangles, int *nVert, int *nEdge, int *nVolu, 
	int *nSurf, int *totNumFaceNode){
	// Extracts Data needed to write out a mesh to tecplot


	*nVert=nTriangles*3;
	*nVolu=1;
	*nSurf=nTriangles;
	*nEdge=nTriangles*3;
	*totNumFaceNode=nTriangles*3;
}

int tecplotfile::VolDataBlock(const triangulation &triout, 
	triarray triangulation::*mp,int nVert,int nVolu, int nVertDat){
	// Prints the Coord and Fill Data blocks to the tecplot file

	int ii,jj,kk;
	int nTri,currInd,currType,nCoord;
	nTri=(triout.*mp).size();
	nCoord=int(triout.meshDep->verts(0)->coord.size());
	// Print vertex Data
	for (jj=0;jj<nCoord;++jj){
		for ( ii = 0; ii<nTri; ++ii){
			for ( kk = 0; kk < (3); ++kk){
				currType=(triout.*mp)(ii)->pointtype[kk];
				currInd=(triout.*mp)(ii)->pointind[kk];
				if(currType==1) {
					this->Print("%.16lf ",triout.meshDep->verts.isearch(currInd)->coord[jj]);
				} else if (currType==2) {
					this->Print("%.16lf ",triout.snakeDep->snakeconn.verts.isearch(currInd)->coord[jj]);
				} else if (currType==3) {
					this->Print("%.16lf ",triout.trivert.isearch(currInd)->coord(jj));
				}
			}
		}
		fprintf(fid,"\n");this->ResetLine();
	}
	for (jj=nCoord;jj<nVertDat;++jj){
		for ( ii = 0; ii<nVert; ++ii){
			this->Print("%lf ",0);
		}
		fprintf(fid,"\n");this->ResetLine();
	}

	// Print Cell Data
	for ( ii = 0; ii<nVolu; ++ii){
		this->Print("%.16lf ",0);
	}
	fprintf(fid,"\n");this->ResetLine();
	for ( ii = 0; ii<nVolu; ++ii){
		this->Print("%.16lf ",0);
	}
	fprintf(fid,"\n");this->ResetLine();
	for ( ii = 0; ii<nVolu; ++ii){
		this->Print("%.16lf ",0);
	}
	fprintf(fid,"\n");this->ResetLine();
	return(0);
}

int tecplotfile::SurfDataBlock(const triangulation &triout, 
	triarray triangulation::*mp,int nVert,int nSurf, int nVertDat){
	// Prints the Coord and Fill Data blocks to the tecplot file

	int ii,jj,kk;
	int nTri,currInd,currType,nCoord;
	nTri=(triout.*mp).size();
	nCoord=int(triout.meshDep->verts(0)->coord.size());
	// Print vertex Data
	for (jj=0;jj<nCoord;++jj){
		for ( ii = 0; ii<nTri; ++ii){
			for ( kk = 0; kk < (3); ++kk){
				currType=(triout.*mp)(ii)->pointtype[kk];
				currInd=(triout.*mp)(ii)->pointind[kk];
				if(currType==1) {
					this->Print("%.16lf ",triout.meshDep->verts.isearch(currInd)->coord[jj]);
				} else if (currType==2) {
					this->Print("%.16lf ",triout.snakeDep->snakeconn.verts.isearch(currInd)->coord[jj]);
				} else if (currType==3) {
					this->Print("%.16lf ",triout.trivert.isearch(currInd)->coord(jj));
				}
			}
		}
		fprintf(fid,"\n");this->ResetLine();
	}
	for (jj=nCoord;jj<nVertDat;++jj){
		for ( ii = 0; ii<nVert; ++ii){
			this->Print("%lf ",0);
		}
		fprintf(fid,"\n");this->ResetLine();
	}

	// Print Cell Data 
	for ( ii = 0; ii<nSurf; ++ii){
		this->Print("%.16lf ",0);
	}
	fprintf(fid,"\n");this->ResetLine();
	for ( ii = 0; ii<nSurf; ++ii){
		this->Print("%.16lf ",0);
	}
	fprintf(fid,"\n");this->ResetLine();
	for ( ii = 0; ii<nSurf; ++ii){
		this->Print("%.16lf ",0);
	}
	fprintf(fid,"\n");this->ResetLine();
	return(0);
}

int tecplotfile::LineDataBlock(const triangulation &triout, 
	triarray triangulation::*mp,int nVert,int nEdge, int nVertDat,int nCellDat){
	// Prints the Coord and Fill Data blocks to the tecplot file

	int ii,jj,kk;
	int nTri,currInd,currType,nCoord;
	nTri=(triout.*mp).size();
	nCoord=int(triout.meshDep->verts(0)->coord.size());
	// Print vertex Data
	for (jj=0;jj<nCoord;++jj){
		for ( ii = 0; ii<nTri; ++ii){
			for ( kk = 0; kk < (3); ++kk){
				currType=(triout.*mp)(ii)->pointtype[kk];
				currInd=(triout.*mp)(ii)->pointind[kk];
				if(currType==1) {
					this->Print("%.16lf ",triout.meshDep->verts.isearch(currInd)->coord[jj]);
				} else if (currType==2) {
					this->Print("%.16lf ",triout.snakeDep->snakeconn.verts.isearch(currInd)->coord[jj]);
				} else if (currType==3) {
					this->Print("%.16lf ",triout.trivert.isearch(currInd)->coord(jj));
				} else {
					RSVS3D_ERROR_ARGUMENT("unknown point type");
				}
			}
		}
		fprintf(fid,"\n");this->ResetLine();
	}
	for (jj=nCoord;jj<nVertDat;++jj){
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


int tecplotfile::LineDataBlock(const triangulation &triout, 
	triarray triangulation::*mp,int nVert,int nEdge, int nVertDat,int nCellDat,
	const vector<int> &triList){
	// Prints the Coord and Fill Data blocks to the tecplot file

	int ii,jj,kk;
	int nTri,currInd,currType,nCoord;
	nTri=int(triList.size());
	nCoord=int(triout.meshDep->verts(0)->coord.size());
	// Print vertex Data
	for (jj=0;jj<nCoord;++jj){
		for ( ii = 0; ii<nTri; ++ii){
			for ( kk = 0; kk < (3); ++kk){
				currType=(triout.*mp).isearch(triList[ii])->pointtype[kk];
				currInd=(triout.*mp).isearch(triList[ii])->pointind[kk];
				if(currType==1) {
					this->Print("%.16lf ",triout.meshDep->verts.isearch(currInd)->coord[jj]);
				} else if (currType==2) {
					this->Print("%.16lf ",triout.snakeDep->snakeconn.verts.isearch(currInd)->coord[jj]);
				} else if (currType==3) {
					this->Print("%.16lf ",triout.trivert.isearch(currInd)->coord(jj));
				} else {
					RSVS3D_ERROR_ARGUMENT("unknown point type");
				}
			}
		}
		fprintf(fid,"\n");this->ResetLine();
	}
	for (jj=nCoord;jj<nVertDat;++jj){
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

int tecplotfile::SurfFaceMap(const triangulation &triout, 
	triarray triangulation::*mp){
	int ii,jj,kk;

	int nTri;
	nTri=(triout.*mp).size();
	//RSVS3D_ERROR_ARGUMENT("Surface Map not supported for triangulation")
	kk=1;
	for (ii=0;ii<nTri;++ii){ // Print Number of vertices per face
		for(jj=0;jj<3;++jj){
			if (jj==2){
				this->Print("%i %i\n",kk-2,kk);
			} else {
				this->Print("%i %i\n",kk,kk+1);
			}
			this->ResetLine();
			++kk;
		}
	}
	

	for (ii=0;ii<nTri;++ii){ // Print connectedFace number
		for(jj=0;jj<3;++jj){
			this->Print("%i ",ii+1);
		}
	}
	fprintf(fid,"\n");this->ResetLine();
	for (ii=0;ii<3*nTri;++ii){ // Print 0
			this->Print("%i ",0);
	}
	fprintf(fid,"\n");this->ResetLine();
	return(0);
}

int tecplotfile::LineFaceMap(const triangulation &triout, 
	triarray triangulation::*mp){
	int ii,jj,kk;

	int nTri;
	nTri=(triout.*mp).size();
	//RSVS3D_ERROR_ARGUMENT("Surface Map not supported for triangulation")
	kk=1;
	for (ii=0;ii<nTri;++ii){ // Print Number of vertices per face
		for(jj=0;jj<3;++jj){
			if (jj==2){
				this->Print("%i %i\n",kk-2,kk);
			} else {
				this->Print("%i %i\n",kk,kk+1);
			}
			this->ResetLine();
			++kk;
		}
	}

	return(0);
}

int tecplotfile::LineFaceMap(const vector<int> &triList){
	int ii,jj,kk;

	int nTri;
	nTri=int(triList.size());
	//RSVS3D_ERROR_ARGUMENT("Surface Map not supported for triangulation")
	kk=1;
	for (ii=0;ii<nTri;++ii){ // Print Number of vertices per face
		for(jj=0;jj<3;++jj){
			if (jj==2){
				this->Print("%i %i\n",kk-2,kk);
			} else {
				this->Print("%i %i\n",kk,kk+1);
			}
			this->ResetLine();
			++kk;
		}
	}

	return(0);
}

int tecplotfile::VolFaceMap(const triangulation &triout, 
	triarray triangulation::*mp,int nSurf){
	int ii,jj,kk,n;
	n=(triout.*mp)(0)->pointind.size();
	for (ii=0;ii<nSurf;++ii){ // Print Number of vertices per face
		this->Print("%i ",n);
	}
	fprintf(fid,"\n");this->ResetLine();
	kk=1;
	for (ii=0;ii<nSurf;++ii){// print ordered  list of vertices in face
		for(jj=0;jj<3;++jj){
			this->Print("%i ",kk);
			kk++;
		}
		fprintf(fid,"\n");this->ResetLine();
	}
	for (ii=0;ii<nSurf;++ii){
		this->Print("%i ",1);
	}
	fprintf(fid,"\n");this->ResetLine();
	for (ii=0;ii<nSurf;++ii){
		this->Print("%i ",0);
	}

	return(0);
}
// Class function Implementation
int tecplotfile::PrintTriangulation(const triangulation &triout, 
	triarray triangulation::*mp,int strandID, double timeStep, 
	int forceOutType, const vector<int> &triList){

	int nVert,nEdge,nVolu,nSurf,totNumFaceNode,nVertDat,nCellDat;

	if(triList.size()==0 || forceOutType!=3){
		ExtractTriData(triout,mp,&nVert,&nEdge,&nVolu, &nSurf, &totNumFaceNode);
	} else {
		ExtractTriData(triList.size(),&nVert,&nEdge,&nVolu, &nSurf, &totNumFaceNode);
	}
	if(nVert>0){
		if(nZones==0){
			fprintf(fid, "VARIABLES = \"X\" ,\"Y\" , \"Z\" ,\"v1\" ,\"v2\", \"v3\"\n" );
		}

		this->NewZone();

		if(strandID>0){
			this->StrandTime(strandID, timeStep);
		}
	// Fixed by the dimensionality of the mesh
		nVertDat=3;
		nCellDat=3;

		if (forceOutType==tecplot::autoselect){
			if(nVolu>0 && triout.snakeDep->Check3D()){
				forceOutType=tecplot::polyhedron; // output as volume data (FEPOLYHEDRON)
			} else if (nSurf>0){
				forceOutType=tecplot::polygon;// output as Surface data (FEPOLYGON)
			} else {
				forceOutType=tecplot::line; // output as line data (FELINESEG)
			}
		}


		if (forceOutType==tecplot::polyhedron){
			this->ZoneHeaderPolyhedron(nVert,nVolu,nSurf,totNumFaceNode,nVertDat,nCellDat);
			this->VolDataBlock(triout,mp,nVert,nVolu, nVertDat);
			this->VolFaceMap(triout,mp,nSurf);
		}
		else if (forceOutType==tecplot::polygon){
			this->ZoneHeaderPolygon(nVert, nEdge,nSurf,nVertDat,nCellDat);
			this->SurfDataBlock(triout,mp,nVert,nSurf, nVertDat);
			this->SurfFaceMap(triout,mp);
		} else if (forceOutType==tecplot::line){
			this->ZoneHeaderFelineseg(nVert, nEdge,nVertDat,nCellDat);
			if(triList.size()==0){
				this->LineDataBlock(triout,mp,nVert,nEdge, nVertDat,nCellDat);
				this->LineFaceMap(triout,mp);
			} else {
				this->LineDataBlock(triout,mp,nVert,nEdge, nVertDat,nCellDat,triList);
				this->LineFaceMap(triList);
			}
			
		}
	}
	return(0);
}


// Functions
void ExtractTriData(const triangulation &triout,
	trisurfarray triangulation::*mp, int *nVert, int *nEdge, int *nVolu, 
	int *nSurf, int *totNumFaceNode){
	// Extracts Data needed to write out a mesh to tecplot
	int ii;

	*nSurf=(triout.*mp).size();
	*nVolu=1;
	*nVert=0;
	*nEdge=0;
	for (ii=0; ii<*nSurf;ii++){
		*nVert+=(triout.*mp)(ii)->indvert.size();
	}

	*nEdge=*nVert;
	*totNumFaceNode=*nVert;
}

int tecplotfile::VolDataBlock(const triangulation &triout,
	trisurfarray triangulation::*mp,int nVert,int nVolu, int nVertDat){
	// Prints the Coord and Fill Data blocks to the tecplot file

	int ii,jj,kk;
	int nTri,currInd,currType,nCoord;
	nTri=(triout.*mp).size();
	nCoord=int(triout.meshDep->verts(0)->coord.size());
	// Print vertex Data
	for (jj=0;jj<nCoord;++jj){
		for ( ii = 0; ii<nTri; ++ii){
			for ( kk = 0; kk < int((triout.*mp)(ii)->typevert.size()); ++kk){
				currType=(triout.*mp)(ii)->typevert[kk];
				currInd=(triout.*mp)(ii)->indvert[kk];
				if(currType==1) {
					this->Print("%.16lf ",triout.meshDep->verts.isearch(currInd)->coord[jj]);
				} else if (currType==2) {
					this->Print("%.16lf ",triout.snakeDep->snakeconn.verts.isearch(currInd)->coord[jj]);
				} else if (currType==3) {
					this->Print("%.16lf ",triout.trivert.isearch(currInd)->coord(jj));
				}
			}
		}
		fprintf(fid,"\n");this->ResetLine();
	}
	for (jj=nCoord;jj<nVertDat;++jj){
		for ( ii = 0; ii<nVert; ++ii){
			this->Print("%lf ",0);
		}
		fprintf(fid,"\n");this->ResetLine();
	}

	// Print Cell Data
	for ( ii = 0; ii<nVolu; ++ii){
		this->Print("%.16lf ",0);
	}
	fprintf(fid,"\n");this->ResetLine();
	for ( ii = 0; ii<nVolu; ++ii){
		this->Print("%.16lf ",0);
	}
	fprintf(fid,"\n");this->ResetLine();
	for ( ii = 0; ii<nVolu; ++ii){
		this->Print("%.16lf ",0);
	}
	fprintf(fid,"\n");this->ResetLine();
	return(0);
}
 
int tecplotfile::SurfDataBlock(const triangulation &triout,
	trisurfarray triangulation::*mp,int nVert,int nSurf, int nVertDat){
	// Prints the Coord and Fill Data blocks to the tecplot file

	int ii,jj,kk;
	int nTri,currInd,currType,nCoord;
	nTri=(triout.*mp).size();
	nCoord=int(triout.meshDep->verts(0)->coord.size());
	// Print vertex Data
	for (jj=0;jj<nCoord;++jj){
		for ( ii = 0; ii<nTri; ++ii){
			for ( kk = 0; kk < int((triout.*mp)(ii)->typevert.size()); ++kk) {
				currType=(triout.*mp)(ii)->typevert[kk];
				currInd=(triout.*mp)(ii)->indvert[kk];
				if(currType==1) {
					this->Print("%.16lf ",triout.meshDep->verts.isearch(currInd)->coord[jj]);
				} else if (currType==2) {
					this->Print("%.16lf ",triout.snakeDep->snakeconn.verts.isearch(currInd)->coord[jj]);
				} else if (currType==3) {
					this->Print("%.16lf ",triout.trivert.isearch(currInd)->coord(jj));
				}
			}
		} 
		fprintf(fid,"\n");this->ResetLine();
	}
	for (jj=nCoord;jj<nVertDat;++jj){
		for ( ii = 0; ii<nVert; ++ii){
			this->Print("%lf ",0);
		}
		fprintf(fid,"\n");this->ResetLine();
	}

	// Print Cell Data 
	for ( ii = 0; ii<nSurf; ++ii){
		this->Print("%.16lf ",0);
	}
	fprintf(fid,"\n");this->ResetLine();
	for ( ii = 0; ii<nSurf; ++ii){
		this->Print("%.16lf ",0);
	}
	fprintf(fid,"\n");this->ResetLine();
	for ( ii = 0; ii<nSurf; ++ii){
		this->Print("%.16lf ",0);
	}
	fprintf(fid,"\n");this->ResetLine();
	return(0);
}

int tecplotfile::LineDataBlock(const triangulation &triout,
	trisurfarray triangulation::*mp,int nVert,int nEdge, int nVertDat,
	int nCellDat){
	// Prints the Coord and Fill Data blocks to the tecplot file

	int ii,jj,kk;
	int nTri,currInd,currType,nCoord;
	nTri=(triout.*mp).size();
	nCoord=int(triout.meshDep->verts(0)->coord.size());
	// Print vertex Data
	for (jj=0;jj<nCoord;++jj){
		for ( ii = 0; ii<nTri; ++ii){
			for ( kk = 0; kk < int((triout.*mp)(ii)->typevert.size()); ++kk){
				currType=(triout.*mp)(ii)->typevert[kk];
				currInd=(triout.*mp)(ii)->indvert[kk];
				if(currType==1) {
					this->Print("%.16lf ",triout.meshDep->verts.isearch(currInd)->coord[jj]);
				} else if (currType==2) {
					this->Print("%.16lf ",triout.snakeDep->snakeconn.verts.isearch(currInd)->coord[jj]);
				} else if (currType==3) {
					this->Print("%.16lf ",triout.trivert.isearch(currInd)->coord(jj));
				} else {
					RSVS3D_ERROR_ARGUMENT("unknown point type");
				}
			}
		}
		fprintf(fid,"\n");this->ResetLine();
	}
	for (jj=nCoord;jj<nVertDat;++jj){
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
int tecplotfile::SurfFaceMap(const triangulation &triout,
	trisurfarray triangulation::*mp){
	int ii,jj,kk;

	int nTri;
	nTri=(triout.*mp).size();
	//RSVS3D_ERROR_ARGUMENT("Surface Map not supported for triangulation")
	kk=1;
	for (ii=0;ii<nTri;++ii){ // Print Number of vertices per face
		for(jj=0;jj<int((triout.*mp)(ii)->typevert.size());++jj){
			if (jj==(int((triout.*mp)(ii)->typevert.size())-1)){
				this->Print("%i %i\n",kk-(int((triout.*mp)(ii)->typevert.size())-1),kk);
			} else {
				this->Print("%i %i\n",kk,kk+1);
			}
			this->ResetLine();
			++kk;
		}
	}
	

	for (ii=0;ii<nTri;++ii){ // Print connectedFace number
		for(jj=0;jj<3;++jj){
			this->Print("%i ",ii+1);
		}
	}
	fprintf(fid,"\n");this->ResetLine();
	for (ii=0;ii<3*nTri;++ii){ // Print 0
			this->Print("%i ",0);
	}
	fprintf(fid,"\n");this->ResetLine();
	return(0);
}

int tecplotfile::LineFaceMap(const triangulation &triout,
	trisurfarray triangulation::*mp){
	int ii,jj,kk;

	int nTri;
	nTri=(triout.*mp).size();
	//RSVS3D_ERROR_ARGUMENT("Surface Map not supported for triangulation")
	kk=1;
	for (ii=0;ii<nTri;++ii){ // Print Number of vertices per face
		for(jj=0;jj<int((triout.*mp)(ii)->typevert.size());++jj){
			if (jj==(int((triout.*mp)(ii)->typevert.size())-1)){
				this->Print("%i %i\n",kk-(int((triout.*mp)(ii)->typevert.size())-1),kk);
			} else {
				this->Print("%i %i\n",kk,kk+1);
			}
			this->ResetLine();
			++kk;
		}
	}

	return(0);
}

int tecplotfile::VolFaceMap(const triangulation &triout,
	trisurfarray triangulation::*mp,int nSurf){
	int ii,jj,kk;

	for (ii=0;ii<nSurf;++ii){ // Print Number of vertices per face
		this->Print("%i ",int((triout.*mp)(ii)->typevert.size()));
	}
	fprintf(fid,"\n");this->ResetLine();
	kk=1;
	for (ii=0;ii<nSurf;++ii){// print ordered  list of vertices in face
		for(jj=0;jj<int((triout.*mp)(ii)->typevert.size());++jj){
			this->Print("%i ",kk);
			kk++;
		}
		fprintf(fid,"\n");this->ResetLine();
	}
	for (ii=0;ii<nSurf;++ii){
		this->Print("%i ",1);
	}
	fprintf(fid,"\n");this->ResetLine();
	for (ii=0;ii<nSurf;++ii){
		this->Print("%i ",0);
	}

	return(0);
}
// Class function Implementation
int tecplotfile::PrintTriangulation(const triangulation &triout,
	trisurfarray triangulation::*mp,int strandID, double timeStep,
	int forceOutType){

	int nVert,nEdge,nVolu,nSurf,totNumFaceNode,nVertDat,nCellDat;

	
	ExtractTriData(triout,mp,&nVert,&nEdge,&nVolu, &nSurf, &totNumFaceNode);
	if(nVert>0){
		if(nZones==0){
			fprintf(fid, "VARIABLES = \"X\" ,\"Y\" , \"Z\" ,\"v1\" ,\"v2\", \"v3\"\n" );
		}

		this->NewZone();

		if(strandID>0){
			this->StrandTime(strandID, timeStep);
		}
	// Fixed by the dimensionality of the mesh
		nVertDat=3;
		nCellDat=3;

		if (forceOutType==tecplot::autoselect){
			if(nVolu>0 && triout.snakeDep->Check3D()){
				forceOutType=tecplot::polyhedron; // output as volume data (FEPOLYHEDRON)
			} else if (nSurf>0){
				forceOutType=tecplot::polygon;// output as Surface data (FEPOLYGON)
			} else {
				forceOutType=tecplot::line; // output as line data (FELINESEG)
			}
		}


		if (forceOutType==tecplot::polyhedron){
			this->ZoneHeaderPolyhedron(nVert,nVolu,nSurf,totNumFaceNode,nVertDat,nCellDat);
			this->VolDataBlock(triout,mp,nVert,nVolu, nVertDat);
			this->VolFaceMap(triout,mp,nSurf);
		}
		else if (forceOutType==tecplot::polygon){
			this->ZoneHeaderPolygon(nVert, nEdge,nSurf,nVertDat,nCellDat);
			this->SurfDataBlock(triout,mp,nVert,nSurf, nVertDat);
			this->SurfFaceMap(triout,mp);
		} else if (forceOutType==tecplot::line){
			this->ZoneHeaderFelineseg(nVert, nEdge,nVertDat,nCellDat);
			this->LineDataBlock(triout,mp,nVert,nEdge, nVertDat,nCellDat);
			this->LineFaceMap(triout,mp);
		}
	}
	return(0);
}



void tecplotfile::ZoneHeaderPolyhedron(int nVert, int nVolu, int nSurf,
	int totNumFaceNode,
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

void tecplotfile::ZoneHeaderPolygon(int nVert,int nEdge,  int nSurf,
	int nVertDat, int nCellDat){

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

void tecplotfile::ZoneHeaderFelineseg(int nVert,int nEdge,  int nVertDat, 
	int nCellDat){

	fprintf(fid, "ZONETYPE = FELINESEG\n");
	fprintf(fid, "NODES = %i\n", nVert);
	fprintf(fid, "ELEMENTS = %i\n",nEdge);
	this->Print( "FACES = %i\n",nEdge);
	fprintf(fid, "VARLOCATION=([%i-%i]=NODAL ,[%i-%i]=CELLCENTERED)\n",1,nVertDat, nVertDat+1,nVertDat+nCellDat);
	fprintf(fid, "DATAPACKING=BLOCK\n");
}


void tecplotfile::ZoneHeaderPolyhedronSnake(int nVert, int nVolu, int nSurf,
	int totNumFaceNode,	int nVertDat, int nCellDat, int nSensDat){

	fprintf(fid, "ZONETYPE = FEPOLYHEDRON\n");
	fprintf(fid, "NODES = %i\n", nVert);
	fprintf(fid, "ELEMENTS = %i\n",nVolu);
	this->Print( "FACES = %i\n",nSurf);
	fprintf(fid, "TOTALNUMFACENODES = %i\n",totNumFaceNode);
	fprintf(fid, "NUMCONNECTEDBOUNDARYFACES = 0\n");
	fprintf(fid, "TOTALNUMBOUNDARYCONNECTIONS = 0\n");
	fprintf(fid, "VARLOCATION=([%i-%i]=NODAL)\n",1,nVertDat+nCellDat+nSensDat);

	fprintf(fid, "DATAPACKING=BLOCK\n");
}

void tecplotfile::ZoneHeaderPolygonSnake(int nVert,int nEdge,  int nSurf,
	int nVertDat, int nCellDat, int nSensDat){

	fprintf(fid, "ZONETYPE = FEPOLYGON\n");
	fprintf(fid, "NODES = %i\n", nVert);
	fprintf(fid, "ELEMENTS = %i\n",nSurf);
	this->Print( "FACES = %i\n",nEdge);
	fprintf(fid, "TOTALNUMFACENODES = %i\n",2*nEdge);
	fprintf(fid, "NUMCONNECTEDBOUNDARYFACES = 0\n");
	fprintf(fid, "TOTALNUMBOUNDARYCONNECTIONS = 0\n");
	fprintf(fid, "VARLOCATION=([%i-%i]=NODAL)\n",1,nVertDat+nCellDat+nSensDat);


	fprintf(fid, "DATAPACKING=BLOCK\n");
}

void tecplotfile::ZoneHeaderFelinesegSnake(int nVert,int nEdge,  int nVertDat,
	int nCellDat, int nSensDat){

	fprintf(fid, "ZONETYPE = FELINESEG\n");
	fprintf(fid, "NODES = %i\n", nVert);
	fprintf(fid, "ELEMENTS = %i\n",nEdge);
	this->Print( "FACES = %i\n",nEdge);
	fprintf(fid, "VARLOCATION=([%i-%i]=NODAL)\n",1,nVertDat+nCellDat+nSensDat);

	fprintf(fid, "DATAPACKING=BLOCK\n");
}

void tecplotfile::ZoneHeaderOrdered(int nVert,  int nVertDat, int nCellDat,
	int nSensDat){

	fprintf(fid, "ZONETYPE = ORDERED\n");
	fprintf(fid, "I = %i\n", nVert);
	fprintf(fid, "VARLOCATION=([%i-%i]=NODAL)\n",1,nVertDat+nCellDat+nSensDat);
	fprintf(fid, "DATAPACKING=BLOCK\n");
}


// Print snake with sensitivity
int tecplotfile::PrintSnakeSensitivity(const triangulation& triRSVS,
	const RSVScalc &calcObj, int strandID, double timeStep, int forceOutType, 
	const vector<int> &vertList){

	int nVert,nEdge,nVolu,nSurf,totNumFaceNode,nVertDat,nCellDat, nSensDat;
	const snake& snakeout = *triRSVS.snakeDep;
	ExtractMeshData(snakeout.snakeconn,&nVert,&nEdge,&nVolu, &nSurf, &totNumFaceNode);
	if(nVert==0){ // Don't print mesh
		return(1);
	}
	
	if(nZones==0){
		fprintf(fid, "VARIABLES = \"X\" ,\"Y\" , \"Z\" ,\"v1\" ,\"v2\", \"v3\"" );
		for (int i = 0; i < calcObj.numConstr(); ++i)
		{
			fprintf(fid, ", \"sens_%i\"",i);
		}
		fprintf(fid, "\n");
	}

	this->NewZone();

	if(strandID>0){
		this->StrandTime(strandID, timeStep);
	}
	// Fixed by the dimensionality of the mesh
	nVertDat=3;
	nCellDat=3;
	nSensDat = calcObj.numConstr();

	if (forceOutType==tecplot::autoselect){
		if(nVolu>0){
			forceOutType=tecplot::polyhedron; // output as volume data (FEPOLYHEDRON)
		} else if (nSurf>0){
			forceOutType=tecplot::polygon;// output as Surface data (FEPOLYGON)
		} else if (nEdge>0){
			forceOutType=tecplot::line; // output as line data (FELINESEG)
		} else {
			forceOutType=tecplot::point;
		}
	}


	if (forceOutType==tecplot::polyhedron){
		this->ZoneHeaderPolyhedronSnake(nVert,nVolu,nSurf,totNumFaceNode,
			nVertDat,nCellDat,nSensDat);
		this->SnakeDataBlock(snakeout,nVert, nVertDat);
		this->RSVScalcDataBlock(triRSVS, calcObj, nVert, nSensDat);
		this->VolFaceMap(snakeout.snakeconn,nSurf);
	}
	else if (forceOutType==tecplot::polygon){
		this->ZoneHeaderPolygonSnake(nVert, nEdge,nSurf,nVertDat,nCellDat,
			nSensDat);
		this->SnakeDataBlock(snakeout,nVert, nVertDat);
		this->RSVScalcDataBlock(triRSVS, calcObj, nVert, nSensDat);
		this->SurfFaceMap(snakeout.snakeconn,nEdge);
	} else if (forceOutType==tecplot::line){
		this->ZoneHeaderFelinesegSnake(nVert, nEdge,nVertDat,nCellDat,
			nSensDat);
		this->SnakeDataBlock(snakeout,nVert, nVertDat);
		this->RSVScalcDataBlock(triRSVS, calcObj, nVert, nSensDat);
		this->LineFaceMap(snakeout.snakeconn,nEdge);
	} else if (forceOutType==tecplot::point){
		if(int(vertList.size())==nVert){
			nVert=0;
			for (int ii=0; ii< int(vertList.size());++ii){
				nVert += int(vertList[ii]); 
			}
		} else if (vertList.size()>0){
			nVert=int(vertList.size());
		}
		this->ZoneHeaderOrdered(nVert,nVertDat,nCellDat,nSensDat);
		this->SnakeDataBlock(snakeout,nVert, nVertDat);
		this->RSVScalcDataBlock(triRSVS,	calcObj, nVert, nSensDat);
		// No map just points
	}
	return(0);
}

int tecplotfile::PrintSnakeGradients(const triangulation& triRSVS,
	const RSVScalc &calcObj, int strandID, double timeStep, int forceOutType, 
	const vector<int> &vertList){

	int nVert,nEdge,nVolu,nSurf,totNumFaceNode,nVertDat,nCellDat, nSensDat;
	const snake& snakeout = *triRSVS.snakeDep;
	ExtractMeshData(snakeout.snakeconn,&nVert,&nEdge,&nVolu, &nSurf, &totNumFaceNode);
	if(nVert==0){ // Don't print mesh
		return(1);
	}
	int nPreConstr =5;
	if(nZones==0){
		fprintf(fid, "VARIABLES = \"X\" ,\"Y\" , \"Z\" ,\"v1\" ,\"v2\", \"v3\"" );
		fprintf(fid, ", \"dObj\"");
		fprintf(fid, ", \"dLag\"");
		fprintf(fid, ", \"HObj\"");
		fprintf(fid, ", \"HConstr\"");
		fprintf(fid, ", \"HLag\"");
		for (int i = 0; i < calcObj.numConstr(); ++i)
		{
			fprintf(fid, ", \"dConstr_%i\"",i);
		}
		fprintf(fid, "\n");
	}

	this->NewZone();

	if(strandID>0){
		this->StrandTime(strandID, timeStep);
	}
	// Fixed by the dimensionality of the mesh
	nVertDat=3;
	nCellDat=3;
	nSensDat = calcObj.numConstr()+nPreConstr;

	if (forceOutType==tecplot::autoselect){
		if(nVolu>0){
			forceOutType=tecplot::polyhedron; // output as volume data (FEPOLYHEDRON)
		} else if (nSurf>0){
			forceOutType=tecplot::polygon;// output as Surface data (FEPOLYGON)
		} else if (nEdge>0){
			forceOutType=tecplot::line; // output as line data (FELINESEG)
		} else {
			forceOutType=tecplot::point;
		}
	}


	if (forceOutType==tecplot::polyhedron){
		this->ZoneHeaderPolyhedronSnake(nVert,nVolu,nSurf,totNumFaceNode,
			nVertDat,nCellDat,nSensDat);
		this->SnakeDataBlock(snakeout,nVert, nVertDat);
		this->RSVScalcDataBlock(triRSVS, calcObj, nVert, nSensDat-nPreConstr,-nPreConstr,2);
		this->VolFaceMap(snakeout.snakeconn,nSurf);
	}
	else if (forceOutType==tecplot::polygon){
		this->ZoneHeaderPolygonSnake(nVert, nEdge,nSurf,nVertDat,nCellDat,
			nSensDat);
		this->SnakeDataBlock(snakeout,nVert, nVertDat);
		this->RSVScalcDataBlock(triRSVS, calcObj, nVert, nSensDat-nPreConstr,-nPreConstr,2);
		this->SurfFaceMap(snakeout.snakeconn,nEdge);
	} else if (forceOutType==tecplot::line){
		this->ZoneHeaderFelinesegSnake(nVert, nEdge,nVertDat,nCellDat,
			nSensDat);
		this->SnakeDataBlock(snakeout,nVert, nVertDat);
		this->RSVScalcDataBlock(triRSVS, calcObj, nVert, nSensDat-nPreConstr,-nPreConstr,2);
		this->LineFaceMap(snakeout.snakeconn,nEdge);
	} else if (forceOutType==tecplot::point){
		if(int(vertList.size())==nVert){
			nVert=0;
			for (int ii=0; ii< int(vertList.size());++ii){
				nVert += int(vertList[ii]); 
			}
		} else if (vertList.size()>0){
			nVert=int(vertList.size());
		}
		this->ZoneHeaderOrdered(nVert,nVertDat,nCellDat,nSensDat);
		this->SnakeDataBlock(snakeout,nVert, nVertDat);
		this->RSVScalcDataBlock(triRSVS,	calcObj, nVert, nSensDat-nPreConstr,-nPreConstr,2);
		// No map just points
	}
	return(0);
}

int tecplotfile::RSVScalcDataBlock(const triangulation& triRSVS, 
	const RSVScalc &calcObj, int nVert, int nSensDat, int sensStart,
	int methodProcess){
	// Prints the Coord and Fill Data blocks to the tecplot file

	int ii,jj;
	// Print vertex Data
	std::vector<double> sensTemp;
	void (RSVScalc::*mp) (const triangulation&, std::vector<double>&, int) const;
	if (methodProcess==1){
		mp = &RSVScalc::ReturnSensitivities;
	} else if(methodProcess==2){
		mp = &RSVScalc::ReturnGradient;
	} else {
		RSVS3D_ERROR_ARGUMENT("Unknown methodProcess value accepted are 1 "
			"(sensitivities) and 2 (derivatives).");
	}

	for (jj=sensStart;jj<nSensDat;++jj){
		(calcObj.*mp)(triRSVS, sensTemp, jj);
		for ( ii = 0; ii<nVert; ++ii){
			this->Print("%.16lf ",sensTemp[ii]);
		}
		fprintf(fid,"\n");this->ResetLine();
	}
	return 0;
}


// Print snake with sensitivity
int tecplotfile::PrintSnakeSensitivityTime(const triangulation& triRSVS, 
	const RSVScalc &calcObj, int strandID, double timeStep, int forceOutType, 
	const vector<int> &vertList){

	int nVert,nEdge,nVolu,nSurf,totNumFaceNode,nVertDat,nCellDat, nSensDat, nSensStep;
	const snake& snakeout = *triRSVS.snakeDep;
	ExtractMeshData(snakeout.snakeconn,&nVert,&nEdge,&nVolu, &nSurf, &totNumFaceNode);
	if(nVert==0){ // Don't print mesh
		return(1);
	}
	nVertDat=3;
	nCellDat=3;
	nSensDat = 1;
	nSensStep = calcObj.numConstr();
	
	if(nZones==0){
		fprintf(fid, "VARIABLES = \"X\" ,\"Y\" , \"Z\" ,\"v1\" ,\"v2\", \"v3\"" );
		for (int i = 0; i < nSensDat; ++i)
		{
			fprintf(fid, ", \"sens_%i\"",i);
		}
		fprintf(fid, "\n");
	}

	this->NewZone();
	if(strandID<=0){
		strandID = rand();
	}
	if(strandID>0){
		this->StrandTime(strandID, timeStep);
	}
	// Fixed by the dimensionality of the mesh

	if (forceOutType==tecplot::autoselect){
		if(nVolu>0){
			forceOutType=tecplot::polyhedron; // output as volume data (FEPOLYHEDRON)
		} else if (nSurf>0){
			forceOutType=tecplot::polygon;// output as Surface data (FEPOLYGON)
		} else if (nEdge>0){
			forceOutType=tecplot::line; // output as line data (FELINESEG)
		} else {
			forceOutType=tecplot::point;
		}
	}


	if (forceOutType==tecplot::polyhedron){
		this->ZoneHeaderPolyhedronSnake(nVert,nVolu,nSurf,totNumFaceNode,
			nVertDat,nCellDat,nSensDat);
		this->SnakeDataBlock(snakeout,nVert, nVertDat);
		this->RSVScalcDataBlock(triRSVS,	calcObj, nVert, nSensDat);
		this->VolFaceMap(snakeout.snakeconn,nSurf);
	}
	else if (forceOutType==tecplot::polygon){
		this->ZoneHeaderPolygonSnake(nVert, nEdge,nSurf,nVertDat,nCellDat,
			nSensDat);
		this->SnakeDataBlock(snakeout,nVert, nVertDat);
		this->RSVScalcDataBlock(triRSVS,	calcObj, nVert, nSensDat);
		this->SurfFaceMap(snakeout.snakeconn,nEdge);
	} else if (forceOutType==tecplot::line){
		this->ZoneHeaderFelinesegSnake(nVert, nEdge,nVertDat,nCellDat,
			nSensDat);
		this->SnakeDataBlock(snakeout,nVert, nVertDat);
		this->RSVScalcDataBlock(triRSVS,	calcObj, nVert, nSensDat);
		this->LineFaceMap(snakeout.snakeconn,nEdge);
	} else if (forceOutType==tecplot::point){
		if(int(vertList.size())==nVert){
			nVert=0;
			for (int ii=0; ii< int(vertList.size());++ii){
				nVert += int(vertList[ii]); 
			}
		} else if (vertList.size()>0){
			nVert=int(vertList.size());
		}
		this->ZoneHeaderOrdered(nVert,nVertDat,nCellDat,nSensDat);
		this->SnakeDataBlock(snakeout,nVert, nVertDat);
		this->RSVScalcDataBlock(triRSVS,	calcObj, nVert, nSensDat);
		// No map just points
	}
	int shareZone = this->ZoneNum();
	double tStepMultiplier = pow(10,- ceil(log10(double(nSensStep))));
	for (int i = nSensDat; i < nSensStep; ++i)
	{
		this->NewZone();
		this->StrandTime(strandID, timeStep+tStepMultiplier*i);
		if (forceOutType==tecplot::polyhedron){
			this->ZoneHeaderPolyhedronSnake(nVert,nVolu,nSurf,totNumFaceNode,
				nVertDat,nCellDat,nSensDat);
		} else if (forceOutType==tecplot::polygon){
			this->ZoneHeaderPolygonSnake(nVert, nEdge,nSurf,nVertDat,nCellDat,
				nSensDat);
		} else if (forceOutType==tecplot::line){
			this->ZoneHeaderFelinesegSnake(nVert, nEdge,nVertDat,nCellDat,
				nSensDat);
		} else if (forceOutType==tecplot::point){
			this->ZoneHeaderOrdered(nVert,nVertDat,nCellDat,nSensDat);
		}
		this->DefShareZoneVolume(shareZone, nVertDat+nCellDat);
		this->RSVScalcDataBlock(triRSVS,	calcObj, nVert, i+1, i);
	}

	return(0);
}

// tecplotfile operations

int tecplotfile::OpenFile (const char *str, const char *mode){
	if (fid!=NULL){
		fclose(fid);
	}
	fid=fopen(str,mode);
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

	fileToOpen="../TESTOUT/tecout.plt";

	errFlag=outmesh.OpenFile(fileToOpen);
	if (errFlag!=0){
		return(errFlag);
	}



	return(0);
}

int TestCompareReadWrite(const char* fileToOpen, mesh &blockGrid,
	tecplotfile &outmesh1){
	
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