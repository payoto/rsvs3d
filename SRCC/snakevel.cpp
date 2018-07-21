
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
		n=triangleRSVS.trivert.size();
		// Still need to recompute coordinates
		for (ii=0;ii<n;++ii){
			if(triangleRSVS.trivert(ii)->parentType==2){
				SurfaceCentroid_fun(triangleRSVS.trivert[ii].coord,
					*(triangleRSVS.snakeDep->snakeconn.surfs.isearch(triangleRSVS.trivert(ii)->parentsurf)),
					 triangleRSVS.snakeDep->snakeconn); 
			}
		}
		triangleRSVS.trivert.SetMaxIndex();
		TriangulateContainer(triangleRSVS.snakeDep->snakeconn, triangleRSVS , 2,surfReTriangulate);
		

		triangleRSVS.snakeDep->snakeconn.surfs.SetNoModif();
	}
	triangleRSVS.PrepareForUse();
}

void TriangulateContainer(const mesh& meshin, triangulation &triangleRSVS , const int typeMesh, const vector<int> &subList){
	int ii,n,nTriS,nTriE,maxIndVert, maxIndTriangle;
	triarray triangulation::*mp;

	if (typeMesh==1){
		mp=&triangulation::stattri;
	} else {
		mp=&triangulation::dynatri;
	}

	maxIndVert=triangleRSVS.trivert.GetMaxIndex();
	maxIndTriangle=int((triangleRSVS.*mp).GetMaxIndex());
	
	nTriS=int((triangleRSVS.*mp).size());
	if (int(subList.size())==0){
		n=meshin.surfs.size();
		for (ii=0; ii<n; ++ii){
			TriangulateSurface(*(meshin.surfs(ii)),meshin,triangleRSVS.*mp, 
				triangleRSVS.trivert, typeMesh, maxIndVert+ii+1);
		}
	} else {
		n=subList.size();
		for (ii=0; ii<n; ++ii){
			TriangulateSurface(*(meshin.surfs(subList[ii])),meshin,triangleRSVS.*mp, 
				triangleRSVS.trivert, typeMesh, maxIndVert+ii+1);
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
	double surfLength,edgeLength;
	trianglepoint surfCentre;
	triangle triangleEdge;

	// Loop around edges calculating centre and length of each edge
	surfLength=0;
	edgeLength=0;
	n=int(surfin.edgeind.size());

	triangleEdge.SetPointType(typeMesh,typeMesh,3);
	triangleEdge.pointind[2]=trivertMaxInd+1;
	triangleEdge.parentsurf=surfin.index;
	if (n>3){
		for(ii=0; ii<n; ++ii){
			meshin.edges.isearch(surfin.edgeind[ii])->GeometricProperties(&meshin,edgeCentre,edgeLength);
			edgeCentre.mult(edgeLength);
			surfCentre.coord.add(edgeCentre.usedata());
			surfLength+=edgeLength;

			triangleEdge.pointind[0]=meshin.edges.isearch(surfin.edgeind[ii])->vertind[0];
			triangleEdge.pointind[1]=meshin.edges.isearch(surfin.edgeind[ii])->vertind[1];
			triangul.push_back(triangleEdge);
		}

		surfCentre.coord.div(surfLength);
		surfCentre.index=trivertMaxInd+1;
		surfCentre.parentsurf=surfin.index;
		surfCentre.parentType=typeMesh;
		trivert.push_back(surfCentre);
	} else {

		triangleEdge.SetPointType(typeMesh,typeMesh,typeMesh);

		triangleEdge.pointind=ConcatenateVectorField(meshin.edges,&edge::vertind,meshin.edges.find_list(surfin.edgeind));
		sort(triangleEdge.pointind);
		unique(triangleEdge.pointind);
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


void SurfaceCentroid_fun(coordvec &coord,const surf &surfin, const mesh& meshin){ 
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


void MeshTriangulation(mesh &meshout,const mesh& meshin,triarray &triangul, tripointarray& trivert){

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

void TriangulateGridSnakeIntersect(triangulation &triangleRSVS){

	// Snaxel Operation:
	// - 

}

// Triangulation class Methods

void triangulation::disp() const{

	cout << "This is a triangulation object" << endl;
} 
void triangulation::PrepareForUse() {
	stattri.PrepareForUse();
	dynatri.PrepareForUse();
	trivert.PrepareForUse();
	
} 

void triangulation::CleanDynaTri(){
	vector<int> triDel,pDEl;
	int ii,jj,n,n2;


	n=dynatri.size();
	for (ii=0;ii<n;++ii){
		if(snakeDep->snakeconn.surfs.find(dynatri(ii)->KeyParent())==-1){
			triDel.push_back(dynatri(ii)->index);
		} else if (snakeDep->snakeconn.surfs.isearch(dynatri(ii)->KeyParent())->returnIsModif()) { 
			triDel.push_back(dynatri(ii)->index);
		}
	}
	n=triDel.size();
	for (ii=0;ii<n;++ii){
		n2=3;//dynatri.isearch(triDel[ii])->pointtype.size();
		for(jj=0;jj<n2;++jj){
			if(dynatri.isearch(triDel[ii])->pointtype[jj]==3){
				pDEl.push_back(dynatri.isearch(triDel[ii])->pointind[jj]);
			}
		}
	}

	dynatri.remove(triDel);
	trivert.remove(pDEl);

	PrepareForUse();
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
#pragma GCC diagnostic pop