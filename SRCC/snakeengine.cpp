#include <iostream>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <ctime>

#include "snakstruct.hpp"
#include "snakeengine.hpp"
#include "arraystructures.hpp"

using namespace std;



// Snake Spawning at Vertex.
void SpawnAtVertex(snake& snakein,int indVert){

	snake newsnake;
	int subVert,nVert, nEdge,nSurf,nVolu;
	bool is3D;
	vector<int> vertInds,edgeInds,surfInds,voluInds,edgeInds2;
	vector<int> vertSubs,edgeSubs,surfSubs,voluSubs;
	vector<int> vertSubsTemp,edgeSubsTemp,surfSubsTemp,voluSubsTemp;
	vector<int>::iterator itVecInt; 
	unordered_multimap<int,int> hashEdgeInds,hashVoluInds,hashSurfInds,hashVertInds;


	is3D=snakein.snakemesh->volus.size()>0;
	// Extract Data corresponding to vertex from Mesh
	subVert=snakein.snakemesh->verts.find(indVert);

	
	edgeInds=snakein.snakemesh->verts(subVert)->edgeind;
	edgeSubs=snakein.snakemesh->edges.find_list(edgeInds);
	surfInds=ConcatenateVectorField(snakein.snakemesh->edges, &edge::surfind, edgeSubs);
	sort(surfInds);
	unique(surfInds);

	surfSubs=snakein.snakemesh->surfs.find_list(surfInds);
	voluInds=ConcatenateVectorField(snakein.snakemesh->surfs, &surf::voluind, surfSubs);
	sort(voluInds);
	unique(voluInds);
	if(is3D){
		voluSubs=snakein.snakemesh->volus.find_list(voluInds);
	} else {

		voluSubs=snakein.snakemesh->volus.find_list(voluInds);
	}
	//OperArrayStructMethod(snakein.snakemesh->surfs, surfSubs, &surf::isready, 
	//	ii, std::logical_and<bool>());
	nVert=edgeInds.size();
	nEdge=surfInds.size();
	nSurf=voluInds.size();
	nVolu=int(is3D);

	newsnake.Init(snakein.snakemesh,nVert,nEdge,nSurf,nVolu);
	// Generates snaxels and vertices
	SpawnAtVertexVert(newsnake,nVert, indVert,subVert, surfInds,edgeInds,
		edgeSubs,hashSurfInds);
	// Generate snake edges
	SpawnAtVertexEdge(newsnake, nEdge ,surfInds, edgeInds,	voluInds,
		surfSubs, hashEdgeInds,hashVoluInds);
	// Generate Snake surfaces
	if (is3D){
		SpawnAtVertexSurf3D(newsnake,nSurf,surfInds ,voluInds,voluSubs,hashSurfInds);
	// Generate Volume
		SpawnAtVertexVolu(newsnake,nSurf);
	} else {
		SpawnAtVertexSurf2D( newsnake, nEdge, voluInds);
	}
	
	snakein.SetMaxIndexNM();

	snakein.MakeCompatible_inplace(newsnake);

	// DO NOT RUN TO MAITAIN orederedge newsnake.PrepareForUse();
	snakein.Concatenate(newsnake);

}

void SpawnAtVertexVert(snake& newsnake, int nVert,int indVert, int subVert, const vector<int> &surfInds,
	const vector<int> &edgeInds,const vector<int> &edgeSubs, unordered_multimap<int,int> &hashSurfInds){
	int ii,jj;
	vector<int> edgeSubsTemp;

	newsnake.snakeconn.verts.PopulateIndices();
	newsnake.snaxs.PopulateIndices();
	for (ii=0;ii<nVert;++ii){
		// Finds the to vertex
		jj=int(newsnake.snakemesh->edges(edgeSubs[ii])->vertind[0]==indVert);
		newsnake.snaxs[ii].set(newsnake.snaxs(ii)->index,0.0,0.5,indVert,
			newsnake.snakemesh->edges(edgeSubs[ii])->vertind[jj],edgeInds[ii],0,-1);

		edgeSubsTemp=FindSubList(newsnake.snakemesh->edges(edgeSubs[ii])->surfind,
			surfInds,hashSurfInds);
		newsnake.snakeconn.verts[ii].edgeind=edgeSubsTemp;
		newsnake.snakeconn.verts[ii].coord=newsnake.snakemesh->verts(subVert)->coord;
	}
	newsnake.snakeconn.verts.ChangeIndices(0,1,0,0);
}

void SpawnAtVertexEdge(snake& newsnake,int nEdge,const vector<int> &surfInds,const vector<int> &edgeInds,
	const vector<int> &voluInds,const vector<int> &surfSubs,unordered_multimap<int,int> &hashEdgeInds, unordered_multimap<int,int> &hashVoluInds){
	int ii,jj,kk;
	vector<int> surfSubsTemp,vertSubsTemp;

	newsnake.snakeconn.edges.PopulateIndices();
	newsnake.snaxedges.PopulateIndices();
	for (ii=0;ii<nEdge;++ii){
		newsnake.snaxedges[ii].surfind=surfInds[ii];

		surfSubsTemp=FindSubList(newsnake.snakemesh->surfs(surfSubs[ii])->voluind,
			voluInds,hashVoluInds);
		newsnake.snakeconn.edges[ii].surfind=surfSubsTemp;

		// Assign vertind (can be done WAY more efficiently the other way round)
		// But liek this we can check the logic
		vertSubsTemp=FindSubList(newsnake.snakemesh->surfs(surfSubs[ii])->edgeind,
			edgeInds,hashEdgeInds);
		kk=0;
		for(jj=0;jj<int(vertSubsTemp.size());++jj){
			if (vertSubsTemp[jj]>=0){
				newsnake.snakeconn.edges[ii].vertind[kk]=vertSubsTemp[jj];
				kk++;
			}
		}
	}
	newsnake.snakeconn.edges.ChangeIndices(1,0,1,0);
	if(!newsnake.Check3D()){
		for (ii=0;ii<nEdge;++ii){
			newsnake.snakeconn.edges[ii].surfind[0]=1;
		}
	}
}
void SpawnAtVertexSurf3D(snake& newsnake,int nSurf,const vector<int> &surfInds, const vector<int> &voluInds,
	const vector<int> &voluSubs,unordered_multimap<int,int> &hashSurfInds){

	int ii,jj;
	vector<int> surfSubsTemp;

	newsnake.snakeconn.surfs.PopulateIndices();
	newsnake.snaxsurfs.PopulateIndices();
	for(ii=0;ii<nSurf;++ii){
		newsnake.snaxsurfs[ii].voluind=voluInds[ii];
		newsnake.snakeconn.surfs[ii].voluind[0]=1;
		// Assign edgeind (can be done WAY more efficiently the other way round)
		// But liek this we can check the logic
		surfSubsTemp=FindSubList(newsnake.snakemesh->volus(voluSubs[ii])->surfind,
			surfInds,hashSurfInds);
		// Needs to be modified to work with 2D (surfSubsTemps does not come out right) 
		for(jj=0;jj<int(surfSubsTemp.size());++jj){
			if (surfSubsTemp[jj]>=0){
				newsnake.snakeconn.surfs[ii].edgeind.push_back(surfSubsTemp[jj]);
				
			}
		}
	}
	newsnake.snakeconn.surfs.ChangeIndices(0,1,0,0);
}

void SpawnAtVertexSurf2D(snake& newsnake,int nEdge, const vector<int> &voluInds){

	int ii,jj;
	vector<int> surfSubsTemp;

	newsnake.snakeconn.surfs.PopulateIndices();
	newsnake.snaxsurfs.PopulateIndices();
	ii=0;
	newsnake.snaxsurfs[ii].voluind=voluInds[ii];
	newsnake.snakeconn.surfs[ii].voluind[0]=0;

	for(jj=0;jj<int(nEdge);++jj){
		newsnake.snakeconn.surfs[ii].edgeind.push_back(jj+1);
	}
	
	//newsnake.snakeconn.surfs.ChangeIndices(0,1,0,0);
}

void SpawnAtVertexVolu(snake& newsnake, int nSurf){
	int ii;

	newsnake.snakeconn.volus.PopulateIndices();
	newsnake.snakeconn.volus[0].surfind.reserve(nSurf);
	for (ii=0;ii<nSurf;++ii){
		newsnake.snakeconn.volus[0].surfind.push_back(ii+1);
	}
}

// Merge vertices in contact

void MergeAllContactVertices(snake &fullsnake, vector<int> &isImpact){

	// in isImpact needs to be hashed to rapidly check
	int ii,jj,nImpacts;
	vector<int>  snaxToRemove, vertSameSub,subVelTo0;
	vector<bool> isImpactDone;
	HashedVector<int,snax> impactInd, impactTarg;
	
	nImpacts=isImpact.size()/2;
	impactInd.vec.reserve(nImpacts);
	impactTarg.vec.reserve(nImpacts);
	snaxToRemove.reserve(nImpacts);
	isImpactDone.assign(nImpacts,false);

	for(ii=0;ii<int(isImpact.size());ii=ii+2){
		impactInd.vec.push_back(isImpact[ii]);
	}
	for(ii=1;ii<int(isImpact.size());ii=ii+2){
		impactTarg.vec.push_back(isImpact[ii]);
	}
	impactInd.GenerateHash();
	impactTarg.GenerateHash();


	for(ii=0; ii<nImpacts; ++ii){
		if(!isImpactDone[ii]){
			isImpactDone[ii]=true;
			if(impactTarg.vec[ii]>0){
				fullsnake.snakeconn.SwitchIndex(1,impactInd.vec[ii],impactTarg.vec[ii]);
				subVelTo0.push_back(fullsnake.snaxs.find(impactTarg.vec[ii])); 

				snaxToRemove.push_back(impactInd.vec[ii]);
				vertSameSub=ReturnDataEqualRange(impactTarg.vec[ii], impactInd.hashTable);
				
				for(jj=0;jj< int(vertSameSub.size()) ; jj++){
					isImpactDone[vertSameSub[jj]]=true;
				}
				vertSameSub=ReturnDataEqualRange(impactTarg.vec[ii], impactTarg.hashTable);
				for(jj=0;jj< int(vertSameSub.size()) ; jj++){
					isImpactDone[vertSameSub[jj]]=true;
				}
			}
		}
	}
	for (ii=0;ii<int(subVelTo0.size());ii++){
		fullsnake.snaxs[subVelTo0[ii]].v=0.0; 
	}
	fullsnake.snaxs.remove(snaxToRemove);
	fullsnake.snakeconn.verts.remove(snaxToRemove);
	
}

void SpawnArrivedSnaxels(snake &fullsnake,const vector<int> &isImpact){

	snake fwdSnake, bwdSnake;
	int start_s,stop_s;

	fwdSnake.Init(fullsnake.snakemesh,0,0,0,0);
	bwdSnake.Init(fullsnake.snakemesh,0,0,0,0);

	// Generate fwd spawn
	
	start_s=clock();
	SpawnArrivedSnaxelsDir(fullsnake,bwdSnake,isImpact,-1);
	SpawnArrivedSnaxelsDir(fullsnake,fwdSnake,isImpact,-2);

	stop_s=clock();
	cout << "Spawn: " << double(stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << "ms  " ;
	start_s=clock();

	bwdSnake.Flip();

	fwdSnake.SetMaxIndexNM();
	fwdSnake.MakeCompatible_inplace(bwdSnake);
	// DO NOT RUN TO MAITAIN orederedge newsnake.PrepareForUse();
	fwdSnake.Concatenate(bwdSnake);

	fullsnake.SetMaxIndexNM();
	fullsnake.MakeCompatible_inplace(fwdSnake);
	// DO NOT RUN TO MAITAIN orederedge newsnake.PrepareForUse();
	fullsnake.Concatenate(fwdSnake);
	stop_s=clock();
	cout << "Concatenate: " << double(stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << "ms  " ;

	cout << " nSnax " << fwdSnake.snaxs.size() << "  ";
	fullsnake.PrepareForUse();

}

void SpawnArrivedSnaxelsDir(const snake &fullsnake,snake &partSnake,const vector<int> &isImpact,int dir){

	int nVert, nEdge, nSurf, nVolu,ii,jj,kk;
	vector<int> vertSpawn,subList;
	nVert=0; nEdge=0; nSurf=0; nVolu=0;jj=-1;
	if(dir==-1){
		for(ii=0;ii<int(isImpact.size());ii=ii+2){
			if(isImpact[ii+1]==dir){
				if(!fullsnake.snakemesh->verts.isearch(fullsnake.snaxs(jj)->fromvert)->isBorder){
					jj=fullsnake.snaxs.find(isImpact[ii]);
					vertSpawn.push_back(fullsnake.snaxs(jj)->fromvert);
					nVolu++;
					nVert=nVert+fullsnake.snakemesh->verts(jj)->edgeind.size();


					subList=fullsnake.snakemesh->edges.find_list(fullsnake.snakemesh->verts(jj)->edgeind);

					for(kk=0;kk<int(subList.size());kk++){
						nEdge=nEdge+fullsnake.snakemesh->edges(subList[kk])->surfind.size();
					}
				}

			}
		}
	}

	if(dir==-2){
		for(ii=0;ii<int(isImpact.size());ii=ii+2){
			if(isImpact[ii+1]==dir){
				jj=fullsnake.snaxs.find(isImpact[ii]);
				if(!fullsnake.snakemesh->verts.isearch(fullsnake.snaxs(jj)->tovert)->isBorder){
					vertSpawn.push_back(fullsnake.snaxs(jj)->tovert);
					nVolu++;
					nVert=nVert+fullsnake.snakemesh->verts(jj)->edgeind.size();


					subList=fullsnake.snakemesh->edges.find_list(fullsnake.snakemesh->verts(jj)->edgeind);

					for(kk=0;kk<int(subList.size());kk++){
						nEdge=nEdge+fullsnake.snakemesh->edges(subList[kk])->surfind.size();
					}
				}

			}
		}
	}

	nSurf=nEdge/2;
	partSnake.reserve(nVert, nEdge, nSurf, nVolu);

	sort(vertSpawn);
	unique(vertSpawn);

	for (ii=0;ii<int(vertSpawn.size());ii++){
		SpawnAtVertex(partSnake,vertSpawn[ii]);
	}

}
void IdentifyMergeConnec(snake &snakein, vector<ConnecRemv> &connecEdit);
void CleanupSnakeConnec(snake &snakein){
	vector<ConnecRemv> connecEdit;

	IdentifyMergeConnec(snakein, connecEdit);
}


void IdentifyMergeConnec(snake &snakein, vector<ConnecRemv> &connecEdit){

	vector<bool> isObjDone;
	//vector<int> objSub;
	int nSnax,nSnaxEdge,nSnaxSurf,ii,nParent;
	ConnecRemv tempConnec;

	nSnax=snakein.snaxs.size();
	nSnaxEdge=snakein.snaxedges.size();
	nSnaxSurf=snakein.snaxsurfs.size();

	isObjDone.reserve(nSnaxEdge);

	isObjDone.assign(nSnax,false);
	for(ii=0; ii<nSnax ; ++ii){
		if(!isObjDone[ii]){
			isObjDone[ii]=true;
			nParent=snakein.snaxs.countparent(snakein.snaxs(ii)->KeyParent());
			if(nParent>1){
				snakein.snaxs.findsiblings(snakein.snaxs(ii)->KeyParent(),tempConnec.rmvind);
			}
		}
	}

}