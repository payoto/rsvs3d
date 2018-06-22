#include <iostream>

#include "meshrefinement.hpp"
#include "mesh.hpp"

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