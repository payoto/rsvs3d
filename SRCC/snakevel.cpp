
#include <iostream>
#include <cstdlib>
#include "arraystructures.hpp"
#include "mesh.hpp"
#include "snake.hpp" 
#include "snakevel.hpp"
#include "RSVSmath.hpp"

using namespace std;	

void CalculateSnakeVel(snake &snakein){

	int ii=0;

	for (ii=0;ii<int(snakein.snaxs.size());++ii){
		if (snakein.snaxs(ii)->isfreeze==1){
			snakein.snaxs[ii].v=(0.5-snakein.snaxs[ii].d)*0.3;
			snakein.snaxs[ii].isfreeze=0;
		}
		//snakein.snaxs[ii].v=(double(rand()%1001)/1000.0+0.5)*snakein.snaxs[ii].v;
		snakein.snaxs[ii].v=(0.4*(double(rand()%1001)/1000.0)+0.8)*snakein.snaxs[ii].v;
	}
}

void TriangulateMesh(mesh& meshin, triangulation &triangleRSVS){

	TriangulateContainer(meshin,triangleRSVS , 1); 
	triangleRSVS.meshDep=&meshin;

}

void TriangulateSnake(snake& snakein, triangulation &triangleRSVS){

	TriangulateContainer(snakein.snakeconn, triangleRSVS , 2);
	triangleRSVS.snakeDep=&snakein;
}

void MaintainTriangulateSnake(triangulation &triangleRSVS){

	vector<int> surfReTriangulate;
	int ii,n;

	if(triangleRSVS.snakeDep!=NULL){
		triangleRSVS.snakeDep->snakeconn.surfs.ReturnModifInd(surfReTriangulate);
		surfReTriangulate=triangleRSVS.snakeDep->snakeconn.surfs.find_list(surfReTriangulate);
		triangleRSVS.CleanDynaTri();
		triangleRSVS.trivert.SetMaxIndex();
		TriangulateContainer(triangleRSVS.snakeDep->snakeconn, triangleRSVS , 2,surfReTriangulate);

		BuildTriSurfGridSnakeIntersect(triangleRSVS);
		surfReTriangulate.clear();
		TriangulateContainer(triangleRSVS.snakeDep->snakeconn, triangleRSVS , 3,surfReTriangulate);

		n=triangleRSVS.trivert.size();
		// Still need to recompute coordinates
		for (ii=0;ii<n;++ii){
			triangleRSVS.CalcTriVertPosDyna(ii);
		}

		triangleRSVS.snakeDep->snakeconn.surfs.SetNoModif();
		triangleRSVS.snakeDep->snakeconn.edges.SetNoModif();
	}
	triangleRSVS.PrepareForUse();
}

void TriangulateContainer(const mesh& meshin, triangulation &triangleRSVS , 
	const int typeMesh, const vector<int> &subList){
	int ii,n,nTriS,nTriE,maxIndVert, maxIndTriangle;
	triarray triangulation::*mp;

	if (typeMesh==1){
		mp=&triangulation::stattri;
	} else if (typeMesh==2) { 
		mp=&triangulation::dynatri;
	} else if (typeMesh==3) {
		mp=&triangulation::intertri; 
	}

	maxIndVert=triangleRSVS.trivert.GetMaxIndex();
	maxIndTriangle=int((triangleRSVS.*mp).GetMaxIndex());
	
	nTriS=int((triangleRSVS.*mp).size());
	if (int(subList.size())==0){ // if there is no list of surfaces to triangulate
		
		if (typeMesh!=3) {
			n=meshin.surfs.size();
			for (ii=0; ii<n; ++ii){
				TriangulateSurface(*(meshin.surfs(ii)),meshin,triangleRSVS.*mp, 
					triangleRSVS.trivert, typeMesh, maxIndVert+ii+1);
			}
		} else { // trisurf triangulation
			n=triangleRSVS.trisurf.size();
			for (ii=0; ii<n; ++ii){
				TriangulateTriSurface(*(triangleRSVS.trisurf(ii)),meshin,triangleRSVS.*mp, 
					triangleRSVS.trivert, typeMesh, maxIndVert+ii+1);
			}

		}
	} else { // if there is a list of specific surfaces to triangulate
		n=subList.size();
		if (typeMesh!=3) {
			for (ii=0; ii<n; ++ii){
				TriangulateSurface(*(meshin.surfs(subList[ii])),meshin,triangleRSVS.*mp, 
					triangleRSVS.trivert, typeMesh, maxIndVert+ii+1);
			}
		} else { // trisurf triangulation
			for (ii=0; ii<n; ++ii){ 
				TriangulateTriSurface(*(triangleRSVS.trisurf(subList[ii])),meshin,triangleRSVS.*mp, 
					triangleRSVS.trivert, typeMesh, maxIndVert+ii+1);
			}
		}

	}
	nTriE=int((triangleRSVS.*mp).size());
	for (ii=0; ii<(nTriE-nTriS); ++ii){
		(triangleRSVS.*mp)[nTriS+ii].index=ii+maxIndTriangle+1;
	}
}

void TriangulateSurface(const surf &surfin,const mesh& meshin, 
	triarray &triangul, tripointarray& trivert, const int typeMesh, int trivertMaxInd){
	// Generates the triangulation parts 
	// typeMesh=1 is a static mesh, type 2 is a dynamic one (snake)

	int ii,n;
	coordvec edgeCentre;
	//double surfLength,edgeLength;
	trianglepoint surfCentre;
	triangle triangleEdge;
	vector<int> orderVert;

	// Loop around edges calculating centre and length of each edge
	//surfLength=0;
	//edgeLength=0;
	n=int(surfin.edgeind.size());
	surfin.OrderedVerts(&meshin,orderVert);
	triangleEdge.SetPointType(typeMesh,typeMesh,3);
	triangleEdge.pointind[2]=trivertMaxInd+1;
	triangleEdge.parentsurf=surfin.index;

	triangleEdge.connec.celltarg=surfin.voluind;
	triangleEdge.connec.constrinfluence={-1.0, 1.0};

	if (n>3){
		for(ii=0; ii<n; ++ii){
			//meshin.edges.isearch(surfin.edgeind[ii])->
			//	GeometricProperties(&meshin,edgeCentre,edgeLength);
			//edgeCentre.mult(edgeLength);
			//surfCentre.coord.add(edgeCentre.usedata());
			//surfLength+=edgeLength;

			// triangleEdge.pointind[0]=meshin.edges.isearch(surfin.edgeind[ii])->vertind[0];
			// triangleEdge.pointind[1]=meshin.edges.isearch(surfin.edgeind[ii])->vertind[1];
			triangleEdge.pointind[0]=orderVert[ii];
			triangleEdge.pointind[1]=orderVert[(ii+1)%n];
			triangul.push_back(triangleEdge);
		}

		//surfCentre.coord.div(surfLength);
		surfCentre.index=trivertMaxInd+1;
		surfCentre.parentsurf=surfin.index;
		surfCentre.parentType=typeMesh;
		if(typeMesh==2){
			surfCentre.nInfluences=n;
		}
		trivert.push_back(surfCentre);
	} else {

		triangleEdge.SetPointType(typeMesh,typeMesh,typeMesh);

		triangleEdge.pointind=orderVert;
		triangleEdge.parentsurf=surfin.index;
		triangul.push_back(triangleEdge);
	}
}

void TriangulateTriSurface(const trianglesurf &surfin,const mesh& meshin, 
	triarray &triangul, tripointarray& trivert, const int typeMesh, int trivertMaxInd){
	// Generates the triangulation parts 
	// typeMesh=1 is a static mesh, type 2 is a dynamic one (snake)

	int ii,n,kk;
	coordvec edgeCentre;
	trianglepoint surfCentre;
	triangle triangleEdge;

	// Loop around edges calculating centre and length of each edge


	n=int(surfin.indvert.size());

	triangleEdge.SetPointType(typeMesh,typeMesh,3);
	triangleEdge.pointind[2]=trivertMaxInd+1;
	triangleEdge.parentsurf=surfin.index; 
	triangleEdge.connec.celltarg=surfin.voluind;
	triangleEdge.connec.constrinfluence={-1.0, 1.0};
	kk=0;
	if (n>3){
		for(ii=0; ii<n; ++ii){

			triangleEdge.SetPointType(surfin.typevert[ii],surfin.typevert[(ii+1)%n],3);

			triangleEdge.pointind[0]=surfin.indvert[ii];
			triangleEdge.pointind[1]=surfin.indvert[(ii+1)%n];
			triangul.push_back(triangleEdge);
			kk+=(surfin.typevert[ii]==2);
		}

		surfCentre.index=trivertMaxInd+1;
		surfCentre.parentsurf=surfin.index;
		surfCentre.parentType=typeMesh;
		if(typeMesh==2){
			surfCentre.nInfluences=kk;
		}
		trivert.push_back(surfCentre);
	} else if (n==3) {

		triangleEdge.SetPointType(surfin.typevert[0],surfin.typevert[1],surfin.typevert[2]);

		triangleEdge.pointind=surfin.indvert;
		triangleEdge.parentsurf=surfin.index;
		triangul.push_back(triangleEdge);
	}
}

void SurfaceCentroid_fun2(coordvec &coord,const surf &surfin, const mesh& meshin){
	int ii,n;
	coordvec edgeCentre;
	double edgeLength,surfLength=0;
	coord.assign(0,0,0);
	n=int(surfin.edgeind.size());
	for(ii=0; ii<n; ++ii){
		meshin.edges.isearch(surfin.edgeind[ii])->GeometricProperties(&meshin,edgeCentre,edgeLength);
		edgeCentre.mult(edgeLength);
		coord.add(edgeCentre.usedata());
		surfLength+=edgeLength;
	}

	coord.div(surfLength);
}


void SnakeSurfaceCentroid_fun(coordvec &coord,const surf &surfin, const mesh& meshin){ 
	int ii,n;
	vector<int> vertind;
	vector<vector<double> const *> veccoord;
	SurfCentroid tempCalc;
	ArrayVec<double> tempCoord,jac,hes;

	coord.assign(0,0,0);
	n=int(surfin.edgeind.size());

	veccoord.reserve(n);
	ConnVertFromConnEdge(meshin, surfin.edgeind,vertind);

	for(ii=0; ii<n; ++ii){
		veccoord.push_back(&(meshin.verts.isearch(vertind[ii])->coord));
	}

	tempCalc.assign(veccoord);
	tempCalc.Calc();

	tempCalc.ReturnDat(tempCoord,jac,hes);
	coord.assign(tempCoord[0][0],tempCoord[1][0],tempCoord[2][0]);
}

void HybridSurfaceCentroid_fun(coordvec &coord,const trianglesurf &surfin, const mesh& meshin,
	const mesh& snakeconn){ 
	int ii,n;
	vector<int> vertind;
	vector<vector<double> const *> veccoord;
	SurfCentroid tempCalc;
	ArrayVec<double> tempCoord,jac,hes;

	coord.assign(0,0,0);

	n=surfin.indvert.size();
	for(ii=0; ii<n; ++ii){
		if(surfin.typevert[ii]==1){
			veccoord.push_back(&(meshin.verts.isearch(surfin.indvert[ii])->coord));
		} else if (surfin.typevert[ii]==2){
			veccoord.push_back(&(snakeconn.verts.isearch(surfin.indvert[ii])->coord));
		}
	}

	tempCalc.assign(veccoord);
	tempCalc.Calc();

	tempCalc.ReturnDat(tempCoord,jac,hes);
	coord.assign(tempCoord[0][0],tempCoord[1][0],tempCoord[2][0]);
}


void MeshTriangulation(mesh &meshout,const mesh& meshin,triarray &triangul, tripointarray& trivert){
	// Adds a triarray and corresponding 
	int ii,jj,kk,ll,n, nSub,subSurf,nSurfInd,mm;
	int nNewVert, nNewEdge, nNewSurf,maxIndVert, maxIndEdge, maxIndSurf,vertSub;
	vector<bool> isTriDone;
	vector<int> tempSub,tempType,tempVert,tempEdge;
	mesh tempMesh;
	vert buildVert;edge buildEdge;surf buildSurf;
	bool flag;
	ConnecRemv tempConnec,tempConnec2;
	vector<ConnecRemv> conn;

	meshout=meshin;

	meshout.PrepareForUse();
	triangul.PrepareForUse();
	trivert.PrepareForUse();
	meshout.SetMaxIndex();

	n=int(triangul.size());
	isTriDone.assign(n,false);

	nNewVert=0; nNewEdge=0; nNewSurf=0;
	maxIndVert=meshout.verts.GetMaxIndex(); 
	maxIndEdge=meshout.edges.GetMaxIndex();
	maxIndSurf=meshout.surfs.GetMaxIndex();
	//cout << " " << maxIndVert << " " << maxIndEdge << " " << maxIndSurf << endl;

	buildEdge.vertind.assign(2,0);
	buildEdge.surfind.assign(2,0);

	for (ii=0 ; ii < n ; ++ii){
		if(!isTriDone[ii]){
			tempType=triangul(ii)->pointtype;
			sort(tempType);
			unique(tempType);
			if (int(tempType.size())>1){
				triangul.findsiblings(triangul(ii)->KeyParent(),tempSub);
				nSub=tempSub.size();
				//cout << " " << nSub ;
				for(jj=0 ; jj< nSub;++jj){
					isTriDone[tempSub[jj]]=true;
				}

				subSurf=meshin.surfs.find(triangul(ii)->KeyParent());
				tempVert=ConcatenateVectorField(meshin.edges,&edge::vertind,meshin.edges.find_list(meshin.surfs(subSurf)->edgeind));
				sort(tempVert);
				unique(tempVert);

				tempEdge=meshin.edges.find_list(meshin.surfs(subSurf)->edgeind);
				// Build nSub edges
				// 1 vertex
				// nSub-1 surfaces


				// Add the vertex
				for(jj=0;jj<3;++jj){
					if(triangul(ii)->pointtype[jj]==3){
						buildVert.index=maxIndVert+nNewVert+1;
						nNewVert++;
						buildVert.coord=trivert.isearch(triangul(ii)->pointind[jj])->coord.usedata(); 
						buildVert.edgeind.clear();
						for(kk=0; kk<int(tempVert.size());++kk){
							buildVert.edgeind.push_back(maxIndEdge+nNewEdge+1+kk);
						}
						tempMesh.verts.push_back(buildVert);
						break;
					}
				}
				// Add the edges 
				ll=0;
				kk=0;
				flag=(meshin.edges(tempEdge[kk])->vertind[ll]
					==meshin.edges(tempEdge[kk+1])->vertind[0]) | 
				(meshin.edges(tempEdge[kk])->vertind[ll]
					==meshin.edges(tempEdge[kk+1])->vertind[1]);
				if(!flag){ll=((ll+1)%2);}
				

				nSub=int(tempEdge.size());
				for(kk=0; kk<nSub;++kk){
					buildEdge.index=(maxIndEdge+nNewEdge+1+kk);
					buildEdge.vertind[0]=buildVert.index;
					buildEdge.vertind[1]=meshin.edges(tempEdge[kk])->vertind[ll];
					vertSub=meshin.verts.find(meshin.edges(tempEdge[kk])->vertind[ll]);
					meshout.verts[vertSub].edgeind.push_back(maxIndEdge+nNewEdge+1+kk);

					buildEdge.surfind[0]=(kk==0)*(meshin.surfs(subSurf)->index)
					+(kk!=0)*(maxIndSurf+nNewSurf+kk);
					buildEdge.surfind[1]=(kk==(nSub-1))*(meshin.surfs(subSurf)->index)
					+(kk!=(nSub-1))*(maxIndSurf+nNewSurf+kk+1);

					tempMesh.edges.push_back(buildEdge);

					if(kk==0){

						meshout.surfs[subSurf].edgeind.clear();
						meshout.surfs[subSurf].edgeind.push_back(maxIndEdge+nNewEdge+1);
						meshout.surfs[subSurf].edgeind.push_back(maxIndEdge+nNewEdge+nSub);
						meshout.surfs[subSurf].edgeind.push_back(meshin.edges(tempEdge[kk])->index);
						
					} else {

						buildSurf=*(meshin.surfs(subSurf));
						buildSurf.index=maxIndSurf+nNewSurf+kk;
						buildSurf.edgeind.clear();
						buildSurf.edgeind.push_back(maxIndEdge+nNewEdge+kk);
						buildSurf.edgeind.push_back(maxIndEdge+nNewEdge+kk+1);
						buildSurf.edgeind.push_back(meshin.edges(tempEdge[kk])->index);

						tempMesh.surfs.push_back(buildSurf);

						nSurfInd=meshout.edges[tempEdge[kk]].surfind.size();
						for (mm=0;mm<nSurfInd; ++mm){
							if(meshout.edges[tempEdge[kk]].surfind[mm]==meshin.surfs(subSurf)->index){
								meshout.edges[tempEdge[kk]].surfind[mm]=maxIndSurf+nNewSurf+kk;
							}
						}
						nSurfInd=buildSurf.voluind.size();
						for (mm=0;mm<nSurfInd; ++mm){
							if(buildSurf.voluind[mm]!=0){
								meshout.volus[meshin.volus.find(buildSurf.voluind[mm])].surfind.push_back(buildSurf.index);
							}
						}
					}


					if(kk<int(tempEdge.size()-1)){
						
						flag=(meshin.edges(tempEdge[kk])->vertind[0]
							==meshin.edges(tempEdge[kk+1])->vertind[ll]) | 
						(meshin.edges(tempEdge[kk])->vertind[1]
							==meshin.edges(tempEdge[kk+1])->vertind[ll]);
						ll=((ll+1)%2)*(flag)+ll*(!flag);
					}
				}
				nNewEdge=nNewEdge+nSub;
				nNewSurf=nNewSurf+nSub-1;

			}
			isTriDone[ii]=true;
		}
	}

	meshout.Concatenate(tempMesh);
	meshout.PrepareForUse();
	//meshout.TestConnectivityBiDir();
}


bool FollowSnaxEdgeConnection(int actSnax, int actSurf,int followSnaxEdge,  const snake &snakeRSVS, vector<bool> &isSnaxEdgeDone, int & returnIndex){
	// Snaxel Operation:
	// *- follow appropriate connection* < HERE
	// - reverse onto edge
	// -if its the final snaxel on the edge -> Go to Vertex
	// - else find the correct snaxel -> Has it been explored?
	//    - No: got to snaxel
	//    - Yes: finish 
	int snaxSub, snaxedgeSub,ii,kk,nE;
	bool flagIter,isRepeat;
	//Build edge index sets to intersect

	snaxSub=snakeRSVS.snakeconn.verts.find(actSnax);
	nE=snakeRSVS.snakeconn.verts(snaxSub)->edgeind.size();
	// Define snaxel edge to follow
	ii=0;isRepeat=0;
	if (followSnaxEdge==0){
		followSnaxEdge=0;
		flagIter=false;
	} else {
		flagIter=true;
	}

	while (ii<nE && !flagIter) {
		followSnaxEdge=snakeRSVS.snakeconn.verts(snaxSub)->edgeind[ii];
		flagIter=snakeRSVS.snaxedges.isearch(followSnaxEdge)->surfind==actSurf;
		ii++;
	}

	if(!flagIter){
		cerr << "Error : Cannot build grid/snake intersection a suitable edge in the surface" << endl;
		cerr << "        was not found." << endl; 
		cerr << "Error in " << __PRETTY_FUNCTION__ << endl;
		throw invalid_argument ("edge not found in surface");
	}
	
	snaxedgeSub=snakeRSVS.snaxedges.find(followSnaxEdge);

	isRepeat=(isSnaxEdgeDone[snaxedgeSub]);
	isSnaxEdgeDone[snaxedgeSub]=true;
	// Find the snaxel to look from
	kk=1;	
	kk=snakeRSVS.snakeconn.edges(snaxedgeSub)->vertind[kk]!=actSnax;

	actSnax=snakeRSVS.snakeconn.edges(snaxedgeSub)->vertind[kk];
	returnIndex=actSnax;
	return(isRepeat);
}

int FollowSnaxelDirection(int actSnax,const snake &snakeRSVS, int &returnIndex, int &returnType, int &actEdge){
	// - reverse onto edge
	// -if its the final snaxel on the edge -> Go to Vertex
	// - else find the correct snaxel -> Has it been explored?
	//    - No: got to snaxel
	//    - Yes: finish 
	// returnType 0 (unassigned) 1 (vertex) 2 ()


	bool dirSnax; // 0 reverse, 1 forward;
	int snaxSub,nSib;
	int ii;
	int currSnaxOrd,nextSnaxOrd, nextSnaxPos, testOrder;
	vector<int> snaxSiblings;

	snaxSub=snakeRSVS.snaxs.find(actSnax);
	dirSnax=snakeRSVS.snaxs(snaxSub)->fromvert < snakeRSVS.snaxs(snaxSub)->tovert;

	snakeRSVS.snaxs.findsiblings(snakeRSVS.snaxs(snaxSub)->edgeind,snaxSiblings);
	actEdge=snakeRSVS.snaxs(snaxSub)->edgeind;
	nSib=snaxSiblings.size();
	if(nSib==1){ // if there is a single snaxel return only that one
		returnType=1;
		returnIndex=snakeRSVS.snaxs(snaxSub)->fromvert;
	} else {
		currSnaxOrd=snakeRSVS.snaxs(snaxSub)->orderedge;
		nextSnaxOrd=0;
		nextSnaxPos=-1;

		if (dirSnax){
			for(ii=0 ; ii<nSib; ++ii){
				testOrder=snakeRSVS.snaxs(snaxSiblings[ii])->orderedge;

				if(testOrder<currSnaxOrd && (testOrder>nextSnaxOrd || nextSnaxPos==-1)){
					nextSnaxPos=ii;
					nextSnaxOrd=testOrder;
				}
			}
		}else{
			for(ii=0 ; ii<nSib; ++ii){
				testOrder=snakeRSVS.snaxs(snaxSiblings[ii])->orderedge;

				if(testOrder>currSnaxOrd && (testOrder<nextSnaxOrd || nextSnaxPos==-1)){
					nextSnaxPos=ii;
					nextSnaxOrd=testOrder;
				}
			}
		}
		if (nextSnaxPos==-1){
			returnType=1;
			returnIndex=snakeRSVS.snaxs(snaxSub)->fromvert;
		} else {
			returnType=2;
			returnIndex=snakeRSVS.snaxs(snaxSiblings[nextSnaxPos])->index;
		}

	}
	return(0);

}

int FollowVertexConnection(int actVert, int prevEdge, const HashedVector<int,int> &edgeSurfInd,
	const HashedVector<int,int> &vertSurfInd, const snake &snakeRSVS, const mesh &meshRSVS, int &returnIndex,
	int &returnType, int &nextEdge){
	// Vertex Operation
	// - Find next edge
	// - if there is a snaxel on that edge 
	//    -> Goto the first snaxel.
	//    -> else go to the other vertex.
	vector<int> tempVert,snaxTemp;
	int ii,jj,kk, nSnax, currOrd, currSnax;
	bool isDir;

	tempVert=vertSurfInd.findall(actVert);

	if(tempVert.size()!=2){
		throw invalid_argument("Edge and Vertex list out of surface are broken");
	}

	jj=1;	
	jj= edgeSurfInd.vec[tempVert[jj]/2]!=prevEdge;
	nextEdge=edgeSurfInd.vec[tempVert[jj]/2];

	snakeRSVS.snaxs.findsiblings(nextEdge,snaxTemp);
	nSnax=snaxTemp.size();
	if(nSnax==0){
		returnType=1;
		returnIndex=vertSurfInd.vec[tempVert[jj]+1-((tempVert[jj]%2)*2)];
	} else if (nSnax==1){
		returnType=2;
		returnIndex=snakeRSVS.snaxs(snaxTemp[0])->index;
	}	else {
		kk=tempVert[jj]%2;

		isDir=meshRSVS.edges.isearch(nextEdge)->vertind[kk]<actVert;
		currOrd=snakeRSVS.snaxs(snaxTemp[0])->orderedge;
		currSnax=snakeRSVS.snaxs(snaxTemp[0])->index;
		if(isDir){
			for (ii=1; ii< nSnax; ++ii) {
				if(snakeRSVS.snaxs(snaxTemp[ii])->orderedge>currOrd){
					currOrd=snakeRSVS.snaxs(snaxTemp[ii])->orderedge;
					currSnax=snakeRSVS.snaxs(snaxTemp[0])->index;
				}

			}
		} else {
			for (ii=1; ii< nSnax; ++ii) {
				if(snakeRSVS.snaxs(snaxTemp[ii])->orderedge<currOrd){
					currOrd=snakeRSVS.snaxs(snaxTemp[ii])->orderedge;
					currSnax=snakeRSVS.snaxs(snaxTemp[0])->index;
				}
			}
		}

		returnType=2;
		returnIndex=currSnax;

	}

	return(0);
}

void BuildTriSurfGridSnakeIntersect(triangulation &triangleRSVS){

	vector<bool> isSnaxEdgeDone;
	int ii,n2,nVert,nSnax;
	int actIndex, actSurf,actEdge,actSurfSub,maxNEdge,isFlip;
	int returnType, returnIndex, returnEdge, actType;
	bool flagDone;
	HashedVector<int,int> hashedEdgeInd, vertSurfList;
	vector<int> edgeSub,pseudoEdgeInd;
	trianglesurf newTrisSurf;

	n2=triangleRSVS.snakeDep->snaxedges.size();
	isSnaxEdgeDone.assign(n2,false);
	pseudoEdgeInd.assign(4,0);
	triangleRSVS.trisurf.clear();
	newTrisSurf.index=0;
	for(ii=0;ii<n2;ii++){
		if(!isSnaxEdgeDone[ii]){

			actIndex=0;
			actEdge=triangleRSVS.snakeDep->snaxedges(ii)->index;
			actSurf=triangleRSVS.snakeDep->snaxedges(ii)->surfind;
			actSurfSub=triangleRSVS.meshDep->surfs.find(actSurf);
			if(actSurfSub>0){
				newTrisSurf.parentsurfmesh=actSurf;
				newTrisSurf.indvert.clear();
				newTrisSurf.typevert.clear();
				newTrisSurf.indvert.reserve(8);
				newTrisSurf.typevert.reserve(8);
				newTrisSurf.index++;
				
				hashedEdgeInd.vec=triangleRSVS.meshDep->surfs(actSurfSub)->edgeind;
				edgeSub=triangleRSVS.meshDep->edges.find_list(hashedEdgeInd.vec);
				hashedEdgeInd.GenerateHash();
				vertSurfList.vec=ConcatenateVectorField(triangleRSVS.meshDep->edges, &edge::vertind, 
					edgeSub);
				vertSurfList.GenerateHash();

				actIndex=triangleRSVS.snakeDep->snakeconn.edges(ii)->vertind[0];
				actType=2;
				flagDone=false;
				nVert=0;nSnax=0;
				maxNEdge=hashedEdgeInd.vec.size();
				// Prepare edge lists and vertlists
				while(!flagDone && nVert<maxNEdge+2){

					if (actType==1){ // Type vert
						newTrisSurf.indvert.push_back(actIndex);
						newTrisSurf.typevert.push_back(actType);
						FollowVertexConnection(actIndex, actEdge, hashedEdgeInd, vertSurfList, *(triangleRSVS.snakeDep),
							*(triangleRSVS.meshDep), returnIndex, returnType, returnEdge);
						pseudoEdgeInd[1]=actEdge;
						pseudoEdgeInd[2]=returnEdge;
						actEdge=returnEdge;
						actType=returnType;
						actIndex=returnIndex;
						nVert++;
					} else if (actType==2) { // Snaxel operations
						flagDone=FollowSnaxEdgeConnection(actIndex, actSurf,actEdge, 
							*(triangleRSVS.snakeDep), isSnaxEdgeDone,returnIndex);
						nSnax++;
						if(!flagDone){
							newTrisSurf.indvert.push_back(actIndex);
							newTrisSurf.typevert.push_back(actType);
							nSnax++;
							actIndex=returnIndex;

							newTrisSurf.indvert.push_back(actIndex);
							newTrisSurf.typevert.push_back(actType);
							returnIndex=-1;returnType=-1;
							FollowSnaxelDirection(actIndex,*(triangleRSVS.snakeDep), returnIndex, returnType,actEdge);
							actType=returnType;
							actIndex=returnIndex;
						} else if (nVert==0){ // Define the direction of the trisurf if no vertex was explored
							pseudoEdgeInd[1]=triangleRSVS.snakeDep->snaxs.isearch(actIndex)->edgeind;
							FollowVertexConnection(triangleRSVS.snakeDep->snaxs.isearch(actIndex)->tovert, 
								pseudoEdgeInd[1], hashedEdgeInd, vertSurfList, *(triangleRSVS.snakeDep),
								*(triangleRSVS.meshDep), returnIndex, returnType, returnEdge);
							pseudoEdgeInd[2]=returnEdge;
						}
					}
					if (actType==2) {
						actEdge=0;
					}
					//cout << endl << "nVert= "<< nVert << "  nSnax=" << nSnax ;
					if (nVert>maxNEdge) {
						cout << "Error";
					}
				}
				
				newTrisSurf.voluind=triangleRSVS.meshDep->surfs(actSurfSub)->voluind;
				isFlip=OrderMatchLists(pseudoEdgeInd, triangleRSVS.meshDep->surfs(actSurfSub)->edgeind, 
							 pseudoEdgeInd[1],pseudoEdgeInd[2]);
				if(isFlip==-1){
					isFlip=newTrisSurf.voluind[0];
					newTrisSurf.voluind[0]=newTrisSurf.voluind[1];
					newTrisSurf.voluind[1]=isFlip;
				}
				triangleRSVS.trisurf.push_back(newTrisSurf);
			}
			isSnaxEdgeDone[ii]=true;
		}
	}


	// Snaxel Operation:
	// - follow appropriate connection
	// - reverse onto edge
	// -if its the final snaxel on the edge -> Go to Vertex
	// - else find the correct snaxel -> Has it been explored?
	//    - No: got to snaxel
	//    - Yes: finish 

	// Vertex Operation
	// - Find next edge
	// - if there is a snaxel on that edge 
	//    -> Goto the first snaxel.
	//    -> else go to the other vertex.

}



// Triangulation class Methods

void triangulation::disp() const{

	cout << "This is a triangulation object" << endl;
} 
void triangulation::PrepareForUse() {
	stattri.PrepareForUse();
	dynatri.PrepareForUse();
	intertri.PrepareForUse();
	trivert.PrepareForUse();
	trisurf.PrepareForUse();
	
} 

void triangulation::CleanDynaTri(){
	vector<int> triDel,pDEl;
	int ii,jj,n,n2;


	n=dynatri.size();
	// remove the triangles of which the surface is gone or modif
	for (ii=0;ii<n;++ii){ 
		if(snakeDep->snakeconn.surfs.find(dynatri(ii)->KeyParent())==-1){
			triDel.push_back(dynatri(ii)->index);
		} else if (snakeDep->snakeconn.surfs.isearch(dynatri(ii)->KeyParent())->returnIsModif()) { 
			triDel.push_back(dynatri(ii)->index);
		}
	}
	n=triDel.size();
	// Remove the surface centroid points that were in the removed triangles
	for (ii=0;ii<n;++ii){ 
		n2=3;//dynatri.isearch(triDel[ii])->pointtype.size();
		for(jj=0;jj<n2;++jj){
			if(dynatri.isearch(triDel[ii])->pointtype[jj]==3){
				pDEl.push_back(dynatri.isearch(triDel[ii])->pointind[jj]);
			}
		}
	}

	n=trivert.size();
	// Remove the surface centroid points that are in the intertri
	for (ii=0;ii<n;++ii){ 
		if(trivert(ii)->parentType==3){
			pDEl.push_back(trivert(ii)->index);
		}
	}

	dynatri.remove(triDel);
	trivert.remove(pDEl);
	intertri.clear();
	trisurf.clear();

	PrepareForUse();
}
void triangulation::CalcTriVertPos(){
	int ii,n;
	n=trivert.size();
	for (ii=0 ; ii<n ; ++ii) {
		CalcTriVertPos(ii);
	}
}
void triangulation::CalcTriVertPos(int ii){ 
	if(trivert(ii)->parentType==1 ){
		SnakeSurfaceCentroid_fun(trivert[ii].coord,
			*(meshDep->surfs.isearch(
				trivert(ii)->parentsurf)),*meshDep); 
	} else if(trivert(ii)->parentType==2 ){
		SnakeSurfaceCentroid_fun(trivert[ii].coord,
			*(snakeDep->snakeconn.surfs.isearch(
				trivert(ii)->parentsurf)),snakeDep->snakeconn); 
	} else if(trivert(ii)->parentType==3){
		HybridSurfaceCentroid_fun(trivert[ii].coord,*(trisurf.isearch(
			trivert(ii)->parentsurf)),*(meshDep),
			snakeDep->snakeconn);  
	}
}
void triangulation::CalcTriVertPosDyna(){
	int ii,n;
	n=trivert.size();
	for (ii=0 ; ii<n ; ++ii) {
		CalcTriVertPosDyna(ii);
	}
}
void triangulation::CalcTriVertPosDyna(int ii){ 
	if(trivert(ii)->parentType==2){
		SnakeSurfaceCentroid_fun(trivert[ii].coord,
			*(snakeDep->snakeconn.surfs.isearch(
				trivert(ii)->parentsurf)),snakeDep->snakeconn); 
	} else if(trivert(ii)->parentType==3){
		HybridSurfaceCentroid_fun(trivert[ii].coord,*(trisurf.isearch(
			trivert(ii)->parentsurf)),*(meshDep),
			snakeDep->snakeconn);  
	}
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
void triangle::disp() const{} 
void triangle::disptree(mesh const&, int) const {}
void triangle::ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu){}
void triangle::read(FILE * fid){}
void triangle::write(FILE * fid) const {}

void trianglepoint::disp() const{} 
void trianglepoint::disptree(mesh const&, int) const {}
void trianglepoint::ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu){}
void trianglepoint::ChangeIndicesSnakeMesh(int nVert,int nEdge,int nSurf,int nVolu){}
void trianglepoint::read(FILE * fid){}
void trianglepoint::write(FILE * fid) const {}

void trianglesurf::disp() const{
	cout << "trianglesurf " << index << " parent " << parentsurfmesh << endl;
	cout << "indvert: " ;DisplayVector(indvert);
	cout << endl;
	cout << "typevert: " ;DisplayVector(typevert);
	cout << endl;
} 
void trianglesurf::disptree(mesh const&, int) const {}
void trianglesurf::ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu){}
void trianglesurf::ChangeIndicesSnakeMesh(int nVert,int nEdge,int nSurf,int nVolu){}
void trianglesurf::read(FILE * fid){}
void trianglesurf::write(FILE * fid) const {}
#pragma GCC diagnostic pop