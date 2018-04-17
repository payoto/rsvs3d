#include <iostream>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <ctime>

#include "snakstruct.hpp"
#include "snakeengine.hpp"
#include "arraystructures.hpp"

using namespace std;


void ConnecRemv::disp() {
	cout << "connrmv:  ki " << keepind << " | tobj " << typeobj << " | rmvind ";
	DisplayVector(rmvind);
	cout << endl;
}

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
		for(jj=0;jj<int(surfSubsTemp.size());++jj){

			newsnake.snakeconn.edges[ii].surfind[jj]++;

		}
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
	newsnake.snakeconn.edges.ChangeIndices(1,0,0,0);
	if(!newsnake.Check3D()){
		for (ii=0;ii<nEdge;++ii){
			newsnake.snakeconn.edges[ii].surfind[0]=1;
			newsnake.snakeconn.edges[ii].surfind[1]=0;
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
	HashedVector<int,int> vertNoSpawn;

	fwdSnake.Init(fullsnake.snakemesh,0,0,0,0);
	bwdSnake.Init(fullsnake.snakemesh,0,0,0,0);

	// Generate fwd spawn
	
	start_s=clock();
	vertNoSpawn.vec.push_back(0);
	vertNoSpawn.GenerateHash();
	SpawnArrivedSnaxelsDir(fullsnake,bwdSnake,isImpact,-1,vertNoSpawn);
	SpawnArrivedSnaxelsDir(fullsnake,fwdSnake,isImpact,-2,vertNoSpawn);

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

void SpawnArrivedSnaxelsDir(snake &fullsnake,snake &partSnake,const vector<int> &isImpact,int dir,
	HashedVector<int,int> &vertNoSpawn){

	int nVert, nEdge, nSurf, nVolu,ii,jj,kk;
	vector<int> vertSpawn,subList;
	int snax::*mp;

	nVert=0; nEdge=0; nSurf=0; nVolu=0;jj=-1;

	if(dir==-1){
		mp=&snax::fromvert;
	} else if (dir==-2){
		mp=&snax::tovert;

	}

	for(ii=0;ii<int(isImpact.size());ii=ii+2){
		if(isImpact[ii+1]==dir){
			jj=fullsnake.snaxs.find(isImpact[ii]);
			if(!fullsnake.snakemesh->verts.isearch(fullsnake.snaxs(jj)->*mp)->isBorder){
				jj=fullsnake.snaxs.find(isImpact[ii]);
				vertSpawn.push_back(fullsnake.snaxs(jj)->*mp);
				nVolu++;
				nVert=nVert+fullsnake.snakemesh->verts.isearch(fullsnake.snaxs(jj)->*mp)->edgeind.size();


				subList=fullsnake.snakemesh->edges.find_list(fullsnake.snakemesh->verts.isearch(fullsnake.snaxs(jj)->*mp)->edgeind);

				for(kk=0;kk<int(subList.size());kk++){
					nEdge=nEdge+fullsnake.snakemesh->edges(subList[kk])->surfind.size();
				}
			} else {
				fullsnake.snaxs[jj].isfreeze=1;
				fullsnake.snaxs.ForceArrayReady();
			}

		}
	}
	

	nSurf=nEdge/2;
	partSnake.reserve(nVert, nEdge, nSurf, nVolu);

	sort(vertSpawn);
	unique(vertSpawn);
	/*
	cout << endl; 
	DisplayVector(vertSpawn);
	cout << endl;
*/
	kk=int(vertSpawn.size());

	for (ii=0;ii<kk;++ii){
		if (vertNoSpawn.find(vertSpawn[ii])==-1){
			SpawnAtVertex(partSnake,vertSpawn[ii]);
		} else {
			cout << endl << "vertex ignored for spawn " <<endl;
		}
	}
	vertNoSpawn.vec.insert(vertNoSpawn.vec.end(),vertSpawn.begin(),vertSpawn.end());
	vertNoSpawn.GenerateHash();
}

void dispconnrmv(vector<ConnecRemv> conn){
	for (int i = 0; i < int(conn.size()); ++i)
	{
		conn[i].disp();
	} 
}


void SnaxEdgeConnecDetection(const snake &snakein, vector<ConnecRemv> &connecEdit){

	int ii,snaxSub1,snaxSub2,nSnaxEdge;
	ConnecRemv tempConnec,tempConnec2;
	nSnaxEdge=snakein.snaxedges.size();

	tempConnec.typeobj=1;
	tempConnec2.typeobj=2;

	for(ii=0;ii<nSnaxEdge;++ii){
		snaxSub1=snakein.snaxs.find(snakein.snakeconn.edges(ii)->vertind[0]);
		snaxSub2=snakein.snaxs.find(snakein.snakeconn.edges(ii)->vertind[1]);
		if(snakein.snaxs(snaxSub1)->edgeind==snakein.snaxs(snaxSub2)->edgeind)
		{
			tempConnec.rmvind.clear();
			tempConnec.keepind=snakein.snaxs(snaxSub1)->index;
			tempConnec.rmvind.push_back(snakein.snaxs(snaxSub2)->index);
			tempConnec2.rmvind.clear();
			tempConnec2.keepind=snakein.snaxedges(ii)->index;
			tempConnec2.rmvind.push_back(snakein.snaxedges(ii)->index);

			connecEdit.push_back(tempConnec);
			connecEdit.push_back(tempConnec2);
		}
	}

}


void SnaxNoConnecDetection(const snake &snakein, vector<ConnecRemv> &connecEdit){

	int ii,nSnax;
	ConnecRemv tempConnec;
	nSnax=snakein.snaxs.size();

	tempConnec.typeobj=1;

	for(ii=0;ii<nSnax;++ii){
		if(snakein.snakeconn.verts(ii)->edgeind.size()==0)
		{
			tempConnec.rmvind.clear();
			tempConnec.keepind=snakein.snaxs(ii)->index;
			tempConnec.rmvind.push_back(snakein.snaxs(ii)->index);

			connecEdit.push_back(tempConnec);
		}
	}

}

void CleanupSnakeConnec(snake &snakein){

	vector<ConnecRemv> connecEdit;

	vector<int> indRmvVert,indRmvEdge,indRmvSurf,indRmvVolu,indDelSurf;
	vector<int> tempSub,isImpact;
	int ii,jj,kk,nEdgeConn,nSurfConn,nEdgeSurfConn,nVertConn,nSnaxConn;
	bool flag, iterFlag;
	HashedVector<int,int> indDelEdge;
	iterFlag=true;
	//indRmvEdge.reserve(snakein.snakeconn.edges.size());
	//indRmvSurf.reserve(snakein.snakeconn.surfs.size());
	auto itVert=indRmvVert.begin();
	auto itEdge=indRmvEdge.begin();
	auto itSurf=indRmvSurf.begin();
	auto itVolu=indRmvVolu.begin();


	while(iterFlag){
		indRmvVert.clear();
		indRmvEdge.clear();
		indRmvSurf.clear();
		indRmvVolu.clear();
		connecEdit.clear();
		indDelSurf.clear();
		indDelEdge.vec.clear();

		snakein.PrepareForUse();
		itVert=indRmvVert.begin();
		itEdge=indRmvEdge.begin();
		itSurf=indRmvSurf.begin();
		itVolu=indRmvVolu.begin();
		// Identify invalid vertex connections
		SnaxNoConnecDetection(snakein, connecEdit);

		nSnaxConn=int(connecEdit.size());
		SnaxEdgeConnecDetection(snakein, connecEdit);
		nVertConn=int(connecEdit.size());
		for(ii=nSnaxConn; ii < nVertConn;ii=ii+2){
			// Skipping the edges which are marked here for removal.
			snakein.snakeconn.SwitchIndex(connecEdit[ii].typeobj,connecEdit[ii].rmvind[0],
				connecEdit[ii].keepind,connecEdit[ii].scopeind);
		}

		// Identify Edge Connections
		IdentifyMergEdgeConnec(snakein, connecEdit);

		nEdgeConn=int(connecEdit.size());
		iterFlag=int(nEdgeConn)>0;
		if(iterFlag){	

			for(ii=nVertConn; ii < nEdgeConn;++ii){
				for(jj=0; jj < int(connecEdit[ii].rmvind.size());++jj){
					snakein.snakeconn.SwitchIndex(connecEdit[ii].typeobj,connecEdit[ii].rmvind[jj],
						connecEdit[ii].keepind,connecEdit[ii].scopeind);
				}
			}

			// Identify surface connections
			if (snakein.Check3D()){
				IdentifyMergSurfConnec(snakein, connecEdit);

				nSurfConn=int(connecEdit.size());
				for(ii=nEdgeConn; ii < nSurfConn;++ii){
					for(jj=0; jj < int(connecEdit[ii].rmvind.size());++jj){
						snakein.snakeconn.SwitchIndex(connecEdit[ii].typeobj,connecEdit[ii].rmvind[jj],
							connecEdit[ii].keepind,connecEdit[ii].scopeind);
					}
				}
			} 

			// Build removal ind of objects
			for(ii=0; ii<nEdgeConn;++ii){
				if (connecEdit[ii].typeobj==2){
					indRmvEdge.insert(itEdge, connecEdit[ii].rmvind.begin(),connecEdit[ii].rmvind.end());
					itEdge=indRmvEdge.end();

				} else if (connecEdit[ii].typeobj==3){
					cerr<< "Should not be here "<< endl;
				} else if (connecEdit[ii].typeobj==5 || connecEdit[ii].typeobj==1){
					indRmvVert.insert(itVert, connecEdit[ii].rmvind.begin(),connecEdit[ii].rmvind.end());
					itVert=indRmvVert.end();
				} 
			}
			sort(indRmvVert);
			sort(indRmvEdge);
			unique(indRmvVert);
			unique(indRmvEdge);

			// Identify Volumes from vertices
			if (snakein.Check3D()){
				ModifyMergVoluConnec(snakein, connecEdit, indRmvVert);

			} else {
				ModifyMergSurf2DConnec(snakein, connecEdit,indRmvVert);
			}


			nEdgeSurfConn=int(connecEdit.size());
			for(ii=nEdgeConn; ii<nEdgeSurfConn;++ii){

				if (connecEdit[ii].typeobj==3){
					indRmvSurf.insert(itSurf, connecEdit[ii].rmvind.begin(),connecEdit[ii].rmvind.end());
					itSurf=indRmvSurf.end();
				} else if (connecEdit[ii].typeobj==4){
					indRmvVolu.insert(itVolu, connecEdit[ii].rmvind.begin(),connecEdit[ii].rmvind.end());
					itVolu=indRmvVolu.end();

				} else if (connecEdit[ii].typeobj==2){
					cerr << "Should not be Here " << endl;

				} else  if (connecEdit[ii].typeobj==5){

					cerr << "Should not be Here " << endl;
				} 
			}

			sort(indRmvSurf);
			unique(indRmvSurf);
			sort(indRmvVolu);
			unique(indRmvVolu);
			// Identify collapsed edges from vertind
			// Remove these like the switch index
			// Look for edges attached to a deleted vertex and delete them as well
			indDelEdge.vec.reserve(nEdgeConn);
			for (ii=0;ii<nEdgeConn;++ii){
				if (connecEdit[ii].typeobj==2){
					if (snakein.snakeconn.edges.isearch(connecEdit[ii].keepind)->vertind[0]==
						snakein.snakeconn.edges.isearch(connecEdit[ii].keepind)->vertind[1])
					{
						indDelEdge.vec.push_back(connecEdit[ii].keepind);
					}
				}
			}
			indDelEdge.GenerateHash();

			// Identify collapsed surfaces if all edges are collapsed in it
			/*
			indDelSurf.reserve(nEdgeSurfConn-nEdgeConn);
			for (ii=nEdgeConn;ii<nEdgeSurfConn;++ii){
				if (connecEdit[ii].typeobj==3){
					flag=true;
					jj=0;
					tempSub=snakein.snakeconn.surfs.isearch(connecEdit[ii].keepind)->edgeind;
					kk=tempSub.size();
					while(flag && jj<kk){
						flag=indDelEdge.find(tempSub[jj])!=-1;
						++jj;
					}
					if(flag){
						indDelSurf.push_back(connecEdit[ii].keepind);
					}
				}
			}
			*/
			kk=indDelEdge.vec.size();
			for (ii=0;ii<kk;++ii){
				snakein.snakeconn.RemoveIndex(2,indDelEdge.vec[ii]);
			}

			kk=int(snakein.snakeconn.surfs.size());
			indDelSurf.reserve(nEdgeSurfConn-nEdgeConn);
			for (ii=0;ii<kk;++ii){
				flag=int(snakein.snakeconn.surfs(ii)->edgeind.size())==0;
				if(flag){
					indDelSurf.push_back(snakein.snakeconn.surfs(ii)->index);
				}
			}

			kk=indDelSurf.size();
			for (ii=0;ii<kk;++ii){
				snakein.snakeconn.RemoveIndex(3,indDelSurf[ii]);
			}
			kk=indRmvSurf.size()+indDelSurf.size()+2;
			//indRmvSurf.reserve(kk);
			indRmvSurf.insert(indRmvSurf.end(), indDelSurf.begin(),indDelSurf.end());

			kk=indRmvEdge.size()+indDelEdge.vec.size()+2;
			//indRmvEdge.reserve(kk);
			indRmvEdge.insert(indRmvEdge.end(), indDelEdge.vec.begin(),indDelEdge.vec.end());

			sort(indRmvEdge);
			sort(indRmvSurf);
			unique(indRmvEdge);
			unique(indRmvSurf);

			snakein.snakeconn.surfs.remove(indRmvSurf);
			snakein.snakeconn.edges.remove(indRmvEdge);
			snakein.snakeconn.verts.remove(indRmvVert);
			snakein.snakeconn.volus.remove(indRmvVolu);
			snakein.snaxs.remove(indRmvVert);
			snakein.snaxedges.remove(indRmvEdge);
			snakein.snaxsurfs.remove(indRmvSurf);

			snakein.snakeconn.TightenConnectivity();

			snakein.PrepareForUse();
			#ifdef SAFE_ALGO
			if (snakein.Check3D()){
				snakein.snakeconn.TestConnectivity();
			}
			#endif
			//tecout.PrintMesh(snakein.snakeconn,2,ttt);
			//tecout.PrintMesh(snakein.snakeconn,3,ttt,3);
			//ttt++;

		}
	}
}


void IdentifyMergEdgeConnec(const snake &snakein, vector<ConnecRemv> &connecEdit){

	vector<bool> isObjDone;
	vector<int> tempSub,tempSub2, tempCount;
	HashedVector<int,int> tempIndHash; 
	//vector<int> objSub;
	int nSnaxEdge, ii,nParent; //nSnax, nSnaxSurf,
	ConnecRemv tempConnec, tempConnec2;

	//nSnax=snakein.snaxs.size();
	nSnaxEdge=snakein.snaxedges.size();
	//nSnaxSurf=snakein.snaxsurfs.size();

	isObjDone.reserve(nSnaxEdge);


	isObjDone.assign(nSnaxEdge,false);
	for(ii=0; ii<nSnaxEdge ; ++ii){
		if(!isObjDone[ii]){
			nParent=snakein.snaxedges.countparent(snakein.snaxedges(ii)->KeyParent());
			if(nParent>1){
				snakein.snaxedges.findsiblings(snakein.snaxedges(ii)->KeyParent(),tempSub);
				IdentifyMergeEdgeGeneral(snakein, isObjDone,connecEdit, tempConnec,  tempConnec2,tempSub,tempSub2, tempCount,tempIndHash);
			}
			isObjDone[ii]=true;
		}
	}

}

void IdentifyMergeEdgeGeneral(const snake &snakein, vector<bool> &isObjDone,vector<ConnecRemv> &connecEdit, ConnecRemv &tempConnec,  ConnecRemv &tempConnec2,vector<int> &tempSub,vector<int> &tempSub2, vector<int> &tempCount, HashedVector<int,int> &tempIndHash) 
{

	int jj,jjStart, nTemp;

	
	// check if the edges are connected
	tempIndHash.vec.clear();
	tempCount.clear();
	tempIndHash.vec=ConcatenateVectorField(snakein.snakeconn.edges, &edge::vertind,tempSub);
	tempIndHash.GenerateHash();
	tempCount=tempIndHash.count(tempIndHash.vec);
	tempCount.push_back(0);
	nTemp=tempCount.size()-1;

	tempConnec2.scopeind.clear();
	for (jj=0;jj<int(tempSub.size());++jj){
		tempConnec2.scopeind.push_back(snakein.snaxedges(tempSub[jj])->index);

	}
	// perform all chains starting at a 1
	jjStart=0;
	while (tempCount[jjStart]!=1 && jjStart<nTemp){jjStart++;}
	while (jjStart<nTemp){ 
		IdentifyMergeEdgeGeneralChain(snakein, isObjDone,connecEdit, tempConnec,  tempConnec2,tempSub,tempSub2, tempCount, tempIndHash,    jjStart);

		jjStart=0;
		while (tempCount[jjStart]!=1 && jjStart<nTemp){jjStart++;}
	}
	// perform all loops starting at a 2
	jjStart=0;
	while (tempCount[jjStart]!=2 && jjStart<nTemp){jjStart++;}
	while (jjStart<nTemp){ 
		IdentifyMergeEdgeGeneralChain(snakein, isObjDone,connecEdit, tempConnec,  tempConnec2,tempSub,tempSub2, tempCount, tempIndHash,    jjStart);
		/*cout << endl;
		tempConnec.disp();
		tempConnec2.disp();*/
		jjStart=0;
		while (tempCount[jjStart]!=2 && jjStart<nTemp){jjStart++;}
	}
}

void IdentifyMergeEdgeGeneralChain(const snake &snakein, vector<bool> &isObjDone,vector<ConnecRemv> &connecEdit, ConnecRemv &tempConnec,  ConnecRemv &tempConnec2,vector<int> &tempSub,vector<int> &tempSub2, vector<int> &tempCount, HashedVector<int,int> &tempIndHash, int jjStart) {

	int jj, jjNext;
	bool flagMoved;
	jj=jjStart;
	tempConnec.rmvind.clear();
	tempConnec2.rmvind.clear();
	tempConnec.typeobj=2;
	tempConnec2.typeobj=5;

	tempCount[jjStart]=0;
	tempConnec.keepind=snakein.snaxedges(tempSub[jjStart/2])->index;
	isObjDone[tempSub[jjStart/2]]=true;
	// if second part of an edge check the other part
	jjNext=jj+(1-((jj%2)*2)); // equivalend of jj+ (jj%2 ? -1 : 1) 
	flagMoved=false;
	while(tempCount[jjNext]>1){ 
			// Builds one group
		flagMoved=true;
		#ifdef SAFE_ALGO
		if (tempCount[jjNext]>2){
			cerr << "Error: Algorithm not conceived for this case "<< endl;
			cerr << " snake has more than 2 edges connected to the same snaxel inside the same surface "<< endl;
			cerr << "	in function:" <<  __PRETTY_FUNCTION__ << endl;
			throw invalid_argument ("Unexpected algorithmic behaviour");
		}
		#endif // SAFE_ALGO

		tempConnec2.rmvind.push_back(tempIndHash.vec[jjNext]);
		tempCount[jjNext]=0;
		tempSub2=tempIndHash.findall(tempIndHash.vec[jjNext]);
		jj=0;
		while(tempSub2[jj]==jjNext && jj<4){++jj;}

		#ifdef SAFE_ALGO
		if (jj>2 || jj<0){
			cerr << "Error: Algorithm not conceived for this case "<< endl;
			cerr << " jj>3 Unsafe read has happened "<< endl;
			cerr << "	in function:" <<  __PRETTY_FUNCTION__ << endl;
			throw invalid_argument ("Unexpected algorithmic behaviour");
		}
		#endif // SAFE_ALGO

		jj=tempSub2[jj];
		tempCount[jj]=0;
		tempConnec.rmvind.push_back(snakein.snaxedges(tempSub[jj/2])->index);
		isObjDone[tempSub[jj/2]]=true;

		jjNext=jj+(1-((jj%2)*2));

	}
	tempConnec2.keepind=tempIndHash.vec[jjNext];
	if (flagMoved){
		connecEdit.push_back(tempConnec);
		connecEdit.push_back(tempConnec2);
	}

}

/*
void IdentifyMergeEdgeGeneralOLD(const snake &snakein, vector<bool> &isObjDone,vector<ConnecRemv> &connecEdit, ConnecRemv &tempConnec,  ConnecRemv &tempConnec2,vector<int> &tempSub,vector<int> &tempSub2, vector<int> &tempCount, HashedVector<int,int> &tempIndHash) 
{

	int jj,jjNext,jjStart, nTemp;


	// check if the edges are connected
	tempIndHash.vec.clear();
	tempCount.clear();
	tempIndHash.vec=ConcatenateVectorField(snakein.snakeconn.edges, &edge::vertind,tempSub);
	tempIndHash.GenerateHash();
	tempCount=tempIndHash.count(tempIndHash.vec);
	nTemp=tempCount.size();

	tempConnec2.scopeind.clear();
	for (jj=0;jj<int(tempSub.size());++jj){
		tempConnec2.scopeind.push_back(snakein.snaxedges(tempSub[jj])->index);

	}

	jjStart=0;
	while (tempCount[jjStart]!=1 && jjStart<nTemp){jjStart++;}
	if (jjStart>=nTemp){ 
	// if all 2s

		tempConnec.rmvind.clear();
		tempConnec.keepind=snakein.snaxedges(tempSub[0])->index;
		tempConnec2.rmvind.clear();
		tempConnec2.keepind=snakein.snakeconn.edges(tempSub[0])->vertind[0];
		tempConnec.typeobj=2;
		tempConnec2.typeobj=5;
		for (jj=0;jj<int(tempSub.size());++jj){
			tempConnec.rmvind.push_back(snakein.snaxedges(tempSub[jj])->index);
			isObjDone[tempSub[jj]]=true;
			tempConnec2.rmvind.push_back(snakein.snakeconn.edges(tempSub[jj])->vertind[0]);
			tempConnec2.rmvind.push_back(snakein.snakeconn.edges(tempSub[jj])->vertind[1]);
			#ifdef SAFE_ALGO
			if (tempCount[jj*2]!=2 && tempCount[jj*2+1]!=2){
				cerr << "Error: Unexpected  behaviour "<< endl;
				cerr << " jjStart not found but vertex does not have 2 connections "<< endl;
				cerr << "	in function:" <<  __PRETTY_FUNCTION__ << endl;
				throw invalid_argument ("Unexpected algorithmic behaviour"); 
			}
			#endif //SAFE_ALGO
		}
		connecEdit.push_back(tempConnec);
		sort(tempConnec2.rmvind);
		unique(tempConnec2.rmvind);
		connecEdit.push_back(tempConnec2);

	} else {
		do{ 
		// group edges by connection groups
			jj=jjStart;
			tempConnec.rmvind.clear();
			tempConnec2.rmvind.clear();
			tempConnec.typeobj=2;
			tempConnec2.typeobj=5;

			tempCount[jjStart]=0;
			tempConnec.keepind=snakein.snaxedges(tempSub[jjStart/2])->index;
			isObjDone[tempSub[jjStart/2]]=true;
						// if second part of an edge check the other part
			jjNext=jj+(1-((jj%2)*2)); // equivalend of jj+ (jj%2 ? -1 : 1) 
			while(tempCount[jjNext]>1){ 
			// Builds one group

				#ifdef SAFE_ALGO
				if (tempCount[jjNext]>2){
					cerr << "Error: Algorithm not conceived for this case "<< endl;
					cerr << " snake has more than 2 edges connected to the same snaxel inside the same surface "<< endl;
					cerr << "	in function:" <<  __PRETTY_FUNCTION__ << endl;
					throw invalid_argument ("Unexpected algorithmic behaviour");
				}
				#endif // SAFE_ALGO

				tempConnec2.rmvind.push_back(tempIndHash.vec[jjNext]);
				tempCount[jjNext]=0;
				tempSub2=tempIndHash.findall(tempIndHash.vec[jjNext]);
				jj=0;
				while(tempSub2[jj]==jjNext && jj<4){++jj;}

				#ifdef SAFE_ALGO
				if (jj>2 || jj<0){
					cerr << "Error: Algorithm not conceived for this case "<< endl;
					cerr << " jj>3 Unsafe read has happened "<< endl;
					cerr << "	in function:" <<  __PRETTY_FUNCTION__ << endl;
					throw invalid_argument ("Unexpected algorithmic behaviour");
				}
				#endif // SAFE_ALGO

				jj=tempSub2[jj];
				tempCount[jj]=0;
				tempConnec.rmvind.push_back(snakein.snaxedges(tempSub[jj/2])->index);
				isObjDone[tempSub[jj/2]]=true;

				jjNext=jj+(1-((jj%2)*2));

			}
			tempConnec2.keepind=tempIndHash.vec[jjNext];
			if (jj!=jjStart){
				connecEdit.push_back(tempConnec);
				connecEdit.push_back(tempConnec2);
			}

			jjStart=0;
			while (tempCount[jjStart]!=1 && jjStart<nTemp){jjStart++;}
			if(jjStart>=nTemp){
				jjStart=0;
				while (tempCount[jjStart]!=2 && jjStart<nTemp){jjStart++;}
			}
		} while (jjStart<nTemp);
	}
}
*/

void IdentifyMergSurfConnec( const snake &snakein, vector<ConnecRemv> &connecEdit){

	vector<bool> isObjDone;
	vector<int> tempSub,tempSub2, tempCount;
	HashedVector<int,int> tempIndHash,edge2Surf; 
	//vector<int> objSub;
	int nSnaxSurf, ii,nParent; //nSnax, nSnaxSurf,
	ConnecRemv tempConnec, tempConnec2;

	//nSnax=snakein.snaxs.size();
	nSnaxSurf=snakein.snaxsurfs.size();
	//nSnaxSurf=snakein.snaxsurfs.size();

	isObjDone.reserve(nSnaxSurf);


	isObjDone.assign(nSnaxSurf,false);
	for(ii=0; ii<nSnaxSurf ; ++ii){
		if(!isObjDone[ii]){
			nParent=snakein.snaxsurfs.countparent(snakein.snaxsurfs(ii)->KeyParent());
			if(nParent>1){
				tempSub.clear();
				snakein.snaxsurfs.findsiblings(snakein.snaxsurfs(ii)->KeyParent(),tempSub);
				IdentifyMergeSurfGeneral(snakein, isObjDone,connecEdit, tempConnec,tempSub,tempSub2, tempCount,edge2Surf,tempIndHash);
			}
			isObjDone[ii]=true;
		}
	}

}



void IdentifyMergeSurfGeneral(const snake &snakein, vector<bool> &isObjDone,vector<ConnecRemv> &connecEdit, 
	ConnecRemv &tempConnec,vector<int> &tempSub,vector<int> &tempSub2,
	vector<int> &tempCount,HashedVector<int,int> &edge2Surf, HashedVector<int,int> &tempIndHash) 
{
	// tempSub is the sub of surfaces in snakeconn that are in the same volume
	int ii,jj,jjStart, nTemp;

	
	edge2Surf.vec.clear();
	tempIndHash.vec=ConcatenateVectorField(snakein.snakeconn.surfs, &surf::edgeind,tempSub);
	// tempIndHash is a hashed vector of concatenate (surfs(tempSub).edgeind)
	for(ii=0; ii<int(tempSub.size());++ii){
		for(jj=0; jj <int(snakein.snakeconn.surfs(tempSub[ii])->edgeind.size()); ++jj){
			edge2Surf.vec.push_back(ii);
		}
	}
	// edge2Surf is a hashed vector of the subscripts into tempSub of the surf matching the 
	// edges in tempHashInd

	tempIndHash.GenerateHash();
	edge2Surf.GenerateHash();
	tempCount=tempIndHash.count(tempIndHash.vec);
	tempCount.push_back(0);
	// tempCount is the vector counting the number of occurences of each edge at each edges location
	nTemp=tempCount.size()-1;




	jjStart=0;
	while (tempCount[jjStart]<=1 && jjStart<nTemp){jjStart++;}
	// jjStart must start at a point were tempCount[jjStart] > 1 otherwise there 
	// is no merging needed in the cell

	if(jjStart<nTemp){ 
		do 
		{
			// if can't find a count above 1 we're done
			tempConnec.rmvind.clear();

			tempSub2=tempIndHash.findall(tempIndHash.vec[jjStart]);
			// tempSub2 is the position of edges matching that detected by jjStart
			tempConnec.typeobj=3;
			isObjDone[tempSub[edge2Surf.vec[jjStart]]]=true;
			tempConnec.keepind=snakein.snakeconn.surfs(tempSub[edge2Surf.vec[jjStart]])->index;
			// Kept index is the last surf to be detected as having the edge of jjStart
			tempCount[jjStart]=0; // To ensure this edge is not set again set tempCount to 0

			IdentifyMergeSurfRecursive( snakein,isObjDone, tempCount,edge2Surf, tempIndHash, 
				tempConnec, tempSub, tempSub2, jjStart);
			if (tempConnec.rmvind.size()>0){
				sort(tempConnec.rmvind);
				unique(tempConnec.rmvind);
				connecEdit.push_back(tempConnec);

				//cout << " " << tempConnec.rmvind.size() << " " << tempSub.size() << " ; ";
			}


			while (tempCount[jjStart]<=1 && jjStart<nTemp)
				{jjStart++;}
		} while (jjStart<nTemp);
	}
	// Note:
	// Check for surface collapse
	// if all of the edges are "collapsed edges" the surfaces need to be assembled in a single surface
	// and made a collapsed surface ie marked for deletion.
}

void IdentifyMergeSurfRecursive(const snake &snakein, vector<bool> &isObjDone, vector<int> &tempCount,const HashedVector<int,int> &edge2Surf, const HashedVector<int,int> &tempIndHash, ConnecRemv &tempConnec, const vector<int> &tempSub, const vector<int> &tempSub2, int excludeSub){

	// tempSub2 is the position of edges matching that detected by jjStart
	// excludeSub is used to not recurse into the caller edge
	int ii, jj ;
	vector<int> tempSurf, tempRecur;

	for(ii = 0 ; ii< int(tempSub2.size()); ++ii){
		if(tempSub2[ii]!=excludeSub){
			if(!isObjDone[tempSub[edge2Surf.vec[tempSub2[ii]]]]){
				tempConnec.rmvind.push_back(snakein.snakeconn.surfs(
					tempSub[edge2Surf.vec[tempSub2[ii]]])->index);

				isObjDone[tempSub[edge2Surf.vec[tempSub2[ii]]]]=true;
				tempCount[tempSub2[ii]]=0; 
				// this edge is explored set to 0;
				// add all edges detected on the same edge as merge targets
				tempSurf=edge2Surf.findall(edge2Surf.vec[tempSub2[ii]]);

				// find all the occurences of that surf
				for (jj=0; jj<int(tempSurf.size()); ++jj){
					if(tempCount[tempSurf[jj]]>1){
					// for each edge of that cell which is tempCount>1 recurs
						tempCount[tempSurf[jj]]=0;
						tempRecur=tempIndHash.findall(tempIndHash.vec[tempSurf[jj]]);
						IdentifyMergeSurfRecursive( snakein,isObjDone, tempCount,edge2Surf, tempIndHash, 
							tempConnec, tempSub, tempRecur, tempSurf[jj]);
					}
				}
			}
		}
	}
}



void ModifyMergVoluConnec(snake &snakein, vector<ConnecRemv> &connecEdit, const vector<int> &indRmvVert){

	
	vector<int> tempSub,tempSub2;	//vector<int> objSub;
	int  ii,jj,kk; //nSnax, nSnaxSurf,
	ConnecRemv tempConnec;

	for(ii=0; ii<int(indRmvVert.size()) ; ++ii){
		// find edges connected to vertex
		tempSub=snakein.snakeconn.edges.find_list(
			snakein.snakeconn.verts.isearch(indRmvVert[ii])->edgeind);

		tempSub2=ConcatenateVectorField(snakein.snakeconn.edges, &edge::surfind, tempSub);
		tempSub=ConcatenateVectorField(snakein.snakeconn.surfs, &surf::voluind, 
			snakein.snakeconn.surfs.find_list(tempSub2));

		sort(tempSub);
		unique(tempSub);

		if (int(tempSub.size())>1){
			tempConnec.typeobj=4;
			kk=0;
			if (tempSub[kk]==0){
				kk++;
			}
			if (int(tempSub.size())>(kk+1)){
				tempConnec.keepind=tempSub[kk];
				tempConnec.rmvind.clear();
				for(jj=kk+1; jj<int(tempSub.size());++jj){
					if(tempSub[jj]>0){
						tempConnec.rmvind.push_back(tempSub[jj]);
					}
				}
				if(int(tempConnec.rmvind.size())>0){
					for (jj=0; jj<int(tempConnec.rmvind.size());jj++){
						snakein.snakeconn.SwitchIndex(tempConnec.typeobj,tempConnec.rmvind[jj],
							tempConnec.keepind,tempConnec.scopeind);
					}
					connecEdit.push_back(tempConnec);
				}
			}

		}

	}

}


void ModifyMergSurf2DConnec(snake &snakein, vector<ConnecRemv> &connecEdit, const vector<int> &indRmvVert){

	
	vector<int> tempSub,tempSub2;	//vector<int> objSub;
	int  ii,jj,kk; //nSnax, nSnaxSurf,
	ConnecRemv tempConnec;

	for(ii=0; ii<int(indRmvVert.size()) ; ++ii){
		// find edges connected to vertex
		tempSub=snakein.snakeconn.edges.find_list(
			snakein.snakeconn.verts.isearch(indRmvVert[ii])->edgeind);

		tempSub2=ConcatenateVectorField(snakein.snakeconn.edges, &edge::surfind, tempSub);
		
		sort(tempSub2);
		unique(tempSub2);

		if (int(tempSub2.size())>1){
			tempConnec.typeobj=3;
			kk=0;
			if (tempSub2[kk]==0){
				kk++;
			}
			if (int(tempSub2.size())>(kk+1)){
				tempConnec.keepind=tempSub2[kk];
				tempConnec.rmvind.clear();
				for(jj=kk+1; jj<int(tempSub2.size());++jj){
					if(tempSub2[jj]>0){
						tempConnec.rmvind.push_back(tempSub2[jj]);
					}
				}
				if(int(tempConnec.rmvind.size())>0){
					for (jj=0; jj<int(tempConnec.rmvind.size());jj++){
						snakein.snakeconn.SwitchIndex(tempConnec.typeobj,tempConnec.rmvind[jj],
							tempConnec.keepind,tempConnec.scopeind);
					}
					connecEdit.push_back(tempConnec);
				}
			}

		}

	}

}