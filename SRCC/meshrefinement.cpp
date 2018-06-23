#include <iostream>
#include <cmath>

#include "mesh.hpp"
#include "meshrefinement.hpp"

using namespace std; 

void CoarsenMesh(const mesh &meshchild, mesh &newparent, const vector<int> &elmMapping){

	int ii,n,nDim;
	bool flag;
	HashedVector<int,int> hashedMapping;
	vector<int> indRmvVert,indRmvEdge,indRmvSurf,indRmvVolu;

	newparent=meshchild;
	nDim=newparent.WhatDim();
	// Check Validity of elmMapping (that all the volumes already exist)

	switch (newparent.WhatDim()){
		case 0: // verts
		throw invalid_argument ("unhandled dimension");
		case 1: // edges
		throw invalid_argument ("unhandled dimension");
		case 2: // surfs
		n=int(newparent.surfs.size());
		flag=false;
		flag=n!=int(elmMapping.size());
		if(!flag){
			for (ii=0;ii < n ; ii++){
				if (newparent.surfs.find(elmMapping[ii])==-1){flag=true;break;}
			}
		}
		case 3: // volumes
		n=int(newparent.volus.size());
		flag=false;
		flag=n!=int(elmMapping.size());
		if(!flag){
			for (ii=0;ii < n ; ii++){
				if (newparent.volus.find(elmMapping[ii])==-1){flag=true;break;}
			}
		}
	}

	if(flag){
		throw invalid_argument ("Element mapping and mesh did not match");
	}

	// Switch element indices according to elmMapping
	switch (nDim){
		case 0:
		for(ii=0;ii< n ;ii++){newparent.SwitchIndex(nDim+1,newparent.verts(ii)->index,elmMapping[ii]);}
			case 1:
		for(ii=0;ii< n ;ii++){newparent.SwitchIndex(nDim+1,newparent.edges(ii)->index,elmMapping[ii]);}
			case 2:
		for(ii=0;ii< n ;ii++){newparent.SwitchIndex(nDim+1,newparent.surfs(ii)->index,elmMapping[ii]);}
			case 3:
		for(ii=0;ii< n ;ii++){newparent.SwitchIndex(nDim+1,newparent.volus(ii)->index,elmMapping[ii]);}
	}
	newparent.TightenConnectivity();

	// Identify Surfs, edges and verts to remove after merging of the higher order elements
	//// This section could be a member function of mesh that removes invalid connections in any mesh
	if(nDim>=3){
		for (ii=0; ii< newparent.surfs.size();++ii){
			if(newparent.surfs(ii)->voluind[0]==newparent.surfs(ii)->voluind[1]){
				newparent.RemoveIndex(3,newparent.surfs(ii)->index);
				indRmvSurf.push_back(newparent.surfs(ii)->index);
			}
		}
	}
	if(nDim>=2){
		for (ii=0; ii< newparent.edges.size();++ii){
			if(newparent.edges(ii)->surfind.size()==0){
				newparent.RemoveIndex(2,newparent.edges(ii)->index);
				indRmvEdge.push_back(newparent.edges(ii)->index);
			}
		}
	}
	if(nDim>=1){
		for (ii=0; ii< newparent.verts.size();++ii){
			if(newparent.verts(ii)->edgeind.size()==0){
				indRmvVert.push_back(newparent.verts(ii)->index);
			}
		}
	}
	//// This section could be a member function of mesh << to here

	// Build list of elements to remove
	hashedMapping.vec=elmMapping;
	hashedMapping.GenerateHash();

	switch (nDim){
		case 0:
		case 1:
		case 2:
		for (ii=0; ii<n ;++ii){
			if(hashedMapping.find(newparent.surfs(ii)->index)==-1){
				indRmvSurf.push_back(newparent.surfs(ii)->index);
			}
		}

		case 3:
		for (ii=0; ii<n ;++ii){
			if(hashedMapping.find(newparent.volus(ii)->index)==-1){
				indRmvVolu.push_back(newparent.volus(ii)->index);
			}
		}
	}

	
	// Remove elements
	sort(indRmvVert);
	sort(indRmvEdge);
	sort(indRmvSurf);
	sort(indRmvVolu);
	unique(indRmvVert);
	unique(indRmvEdge);
	unique(indRmvSurf);
	unique(indRmvVolu);
	newparent.surfs.remove(indRmvSurf);
	newparent.edges.remove(indRmvEdge);
	newparent.verts.remove(indRmvVert);
	newparent.volus.remove(indRmvVolu);

	newparent.PrepareForUse();
}


void CartesianMapping(const mesh& meshin, vector<int> &elmMapping, vector<int> &dims){

	int ii, jj, kk, n,nBox,sub;
	coordvec minCoord,maxCoord,cellCoord,deltaCoord;
	vector<int> boxMap,vertList,edgeList;


	n=meshin.verts.size();
	// Define bounds;
	minCoord=meshin.verts(0)->coord;
	maxCoord=meshin.verts(0)->coord;
	for(ii=1;ii<n;ii++){
		minCoord.min(meshin.verts(ii)->coord);
		maxCoord.max(meshin.verts(ii)->coord);
	}
	deltaCoord=maxCoord;
	deltaCoord.substract(minCoord.usedata());
	// integer vector for each possible locations based on dim split
	nBox=1;
	for(ii=0;ii<3;++ii){
		dims[ii]=dims[ii]>0?dims[ii]:1;
		nBox=nBox*dims[ii];
	}
	boxMap.assign(nBox,0);
	elmMapping.clear();
	elmMapping.assign(meshin.volus.size(),0);
	// Go through volumes calculating their location find which box it fits in
	n=meshin.volus.size();
	for(ii=0;ii<n;ii++){

		edgeList=ConcatenateVectorField(meshin.surfs,&surf::edgeind,
			meshin.surfs.find_list(meshin.volus(ii)->surfind));
		vertList=ConcatenateVectorField(meshin.edges,&edge::vertind,
			meshin.edges.find_list(edgeList));
		vertList=meshin.verts.find_list(vertList);
		sort(vertList);
		unique(vertList);	
		kk=int(vertList.size());
		cellCoord.assign(0.0, 0.0 , 0.0);
		for (jj=0;jj<kk;++jj){
			cellCoord.add(meshin.verts(vertList[jj])->coord);
		}
		cellCoord.div(double(kk));

		cellCoord.substract(minCoord.usedata());
		cellCoord.div(deltaCoord.usedata());
		for (jj=0;jj<3;++jj){
			cellCoord[jj]=cellCoord[jj]*double(dims[jj]);
			cellCoord[jj]=floor(cellCoord[jj]);
			cellCoord[jj]=(cellCoord[jj]>=double(dims[jj]))?double(dims[jj]-1):cellCoord[jj];
		}
		sub=int(cellCoord[0])+(dims[0])*int(cellCoord[1])
			+(dims[0])*(dims[1])*int(cellCoord[2]);
		if (sub>=nBox){
			throw invalid_argument("sub was larger than available size");
		}
		if(boxMap[sub]==0){
		// if the box as never been explored put the index in it
			boxMap[sub]=meshin.volus(ii)->index;
			elmMapping[ii]=boxMap[sub];
		} else {
		// put the number of the box in the elmMapping
			elmMapping[ii]=boxMap[sub];
		}
	}


}