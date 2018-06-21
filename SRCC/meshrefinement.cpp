#include <iostream>

#include "meshrefinement.hpp"
#include "mesh.hpp"

using namespace std; 

void CoarsenMesh(const mesh &meshchild, mesh &newparent, const vector<int> &elmMapping){

	int ii,n;
	bool flag;

	newparent=meshchild;
	n=int(newparent.volus.size());
	// Check Validity of elmMapping (that all the volumes already exist)
	flag=false;
	flag=n==int(elmMapping.size());
	if(!flag){
		for (ii=0;ii < n ; ii++){
			if (newparent.volus.find(elmMapping[ii])==-1){
				flag=true;
				break;
			}
		}
	}	

	// Switch element indices according to elmMapping

	// Build list of elements to remove

	// Remove elements
}