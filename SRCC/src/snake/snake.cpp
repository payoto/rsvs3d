#define _USE_MATH_DEFINES

#include <iostream>
#include <stdexcept>
#include <vector>
#include <cmath>

#include "snake.hpp"
#include "warning.hpp"

using namespace std;

// Implementation of features


// Class Method implementaton 


// --------------------------------------------------------------------------
// Implementatation of snax - snaxedge - snaxsurf 
// --------------------------------------------------------------------------
void snax::disp() const{
	cout << "snax #" << index << "; d " << d  << "; v " << v  
		<< "; edgeind " << edgeind << "; fromvert " << fromvert  
		<< "; tovert " << tovert  << "; isfreeze " << isfreeze 
		<< "; orderedge " << orderedge << " | isBorder " << isBorder << endl;
}

void snaxedge::disp() const{
	cout << "snaxedge #" << index << " | surfind " << surfind 
		<< " " << " | isBorder " << isBorder << " " ;
	normvector.disp();
}

void snaxsurf::disp() const{
		cout << "snaxsurf #" << index << index << " | voluind "
		<< voluind << " " << " | isBorder " << isBorder << " ";
	normvector.disp();
}

// file i/o
void snax::write(FILE *fid) const{

	fprintf(fid, "%i %lf %lf %i %i %i %i %i \n",index,d,v,edgeind, fromvert,
		tovert,isfreeze, orderedge);
}

void snaxedge::write(FILE *fid) const{
	fprintf(fid, "%i %i %lf %lf %lf \n",index,surfind,normvector(0),
		normvector(1),normvector(2));
}

void snaxsurf::write(FILE *fid) const{
	fprintf(fid, "%i %i %lf %lf %lf \n",index,voluind,normvector(0),
		normvector(1),normvector(2));
}
void snax::read(FILE *fid) {

	fscanf(fid,"%i %lf %lf %i %i %i %i %i ", &index, &d, &v, &edgeind, 
		&fromvert, &tovert,	&isfreeze, &orderedge);
}

void snaxedge::read(FILE *fid) {
	fscanf(fid,"%i %i %lf %lf %lf ",&index,&surfind,&normvector[0],
		&normvector[1],&normvector[2]);
}

void snaxsurf::read(FILE *fid) {
	fscanf(fid,"%i %i %lf %lf %lf ",&index,&voluind,&normvector[0],
		&normvector[1],&normvector[2]);
}

void snaxedge::PrepareForUse() {
	normvector.PrepareForUse();
}

void snaxsurf::PrepareForUse() {
	normvector.PrepareForUse();
}


#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
void snax::disptree(const mesh &meshin, int n) const{
   
   disp();
   meshin.verts.isearch(index)->disptree(meshin, n);
}
void snaxedge::disptree(const mesh &meshin, int n) const{
   
   disp();
   meshin.edges.isearch(index)->disptree(meshin, n);
}
void snaxsurf::disptree(const mesh &meshin, int n) const{
   
   disp();
   meshin.surfs.isearch(index)->disptree(meshin, n);

}
void snax::disptree(const snake &snakein, int n) const{
   
   int i;
	for(i=0;i<n;i++){cout<< "  ";}
	disp();
	for(i=0;i<n;i++){cout<< "  ";}
	snakein.snakeconn.verts.isearch(index)->disp();
	if(n>0){
		for (i=0; unsigned_int(i)<
			snakein.snakeconn.verts.isearch(index)->edgeind.size();i++){
			snakein.snaxedges.isearch(
				snakein.snakeconn.verts.isearch(index)->edgeind[i]
				)->disptree(snakein,n-1);
		}
	}

}
void snaxedge::disptree(const snake &snakein, int n) const{
   
   int i;
	for(i=0;i<n;i++){cout<< "  ";}
	disp();
	for(i=0;i<n;i++){cout<< "  ";}
	snakein.snakeconn.edges.isearch(index)->disp();
	if(n>0){
		for (i=0; unsigned_int(i)<
			snakein.snakeconn.edges.isearch(index)->vertind.size();i++){
			snakein.snaxs.isearch( 
				snakein.snakeconn.edges.isearch(index)->vertind[i]
				)->disptree(snakein,n-1);
		}
		for (i=0; unsigned_int(i)<
			snakein.snakeconn.edges.isearch(index)->surfind.size();i++){
			snakein.snaxsurfs.isearch(
				snakein.snakeconn.edges.isearch(index)->surfind[i]
				)->disptree(snakein,n-1);
		}
	}

}
void snaxsurf::disptree(const snake &snakein, int n) const{
   
	int i;
	for(i=0;i<n;i++){cout<< "  ";}
	disp();
	for(i=0;i<n;i++){cout<< "  ";}
	snakein.snakeconn.surfs.isearch(index)->disp();
	if(n>0){
		for (i=0; unsigned_int(i)<
			snakein.snakeconn.surfs.isearch(index)->edgeind.size();i++){
			snakein.snaxedges.isearch(
				snakein.snakeconn.surfs.isearch(index)->edgeind[i]
				)->disptree(snakein,n-1);
		}
	}
}

void snax::ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu){
	index+=nVert;
}
void snaxedge::ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu){
	index+=nEdge;
}
void snaxsurf::ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu){
	index+=nSurf;
}
void snax::ChangeIndicesSnakeMesh(int nVert,int nEdge,int nSurf,int nVolu){
	fromvert+=nVert;
	tovert+=nVert;
	edgeind+=nEdge;
}
void snaxedge::ChangeIndicesSnakeMesh(int nVert,int nEdge,int nSurf,int nVolu){
	surfind+=nSurf;
}
void snaxsurf::ChangeIndicesSnakeMesh(int nVert,int nEdge,int nSurf,int nVolu){
	voluind+=nVolu;
}
#pragma GCC diagnostic pop

// ---------------------------------------------------------------------------
// Extensions of SnakStruct for snaxarray
// ---------------------------------------------------------------------------


bool snaxarray::checkready() 
{	
	this->readyforuse=SnakStruct<snax>::checkready();
	this->readyforuse=(this->readyforuse && this->isOrderedOnEdge);
	return(this->readyforuse);
}
void snaxarray::PrepareForUse()
{
	SnakStruct<snax>::PrepareForUse();
	if (isOrderedOnEdge==0){
		this->ReorderOnEdge();
	}
	readyforuse=true;
}
void snaxarray::Concatenate(const snaxarray &other)
{
	SnakStruct<snax>::Concatenate(other);
	isOrderedOnEdge=0;
}
void snaxarray::ForceArrayReady()
{
	SnakStruct<snax>::ForceArrayReady();
	isOrderedOnEdge=1;
}

// -------------------------------------------------------------------------------------------
// Implementatation of snake
// -------------------------------------------------------------------------------------------

void snake::disp() const{
	snaxs.disp();
	snaxedges.disp();
	snaxsurfs.disp();
	snakeconn.disp();
	cout << "Snaking Mesh with  ";
	snakemesh->displight();
}
void snake::displight() const{

	cout << "snake: snax " << snaxs.size();
	cout << "; snaxedges " << snaxedges.size();
	cout << "; snaxsurfs " << snaxsurfs.size() << endl;

	cout << "Snake with connectivity  ";
	snakeconn.displight();

	cout << "Snaking Mesh with  ";
	snakemesh->displight();
}

void snake::read(FILE *fid){
	snaxs.read(fid);
	snaxedges.read(fid);
	snaxsurfs.read(fid);

	snakeconn.read(fid);
	
}
void snake::write(FILE *fid) const {
   //fprintf(fid,"%i \n",int(borderIsSet));
	snaxs.write(fid);
	snaxedges.write(fid);
	snaxsurfs.write(fid);

	snakeconn.write(fid);
}
int snake::read(const char *str) {
   // Convenience read function taking a single char argument
	FILE *fid;
	fid = fopen( str, "r");
	if (fid!=NULL){
		this->read(fid);
		fclose(fid);
	} else {
		cout << "File " << str << "Could not be opened to read" << endl;
		return(1);
	}
	return(0);
}
int snake::write(const char *str) const {
	FILE *fid;
	fid = fopen( str, "w");
	if (fid!=NULL){
		this->write(fid);
		fclose(fid);
	} else {
		cout << "File " << str << "Could not be opened to write" << endl;
		return(1);
	}
	return(0);
}

bool snake::isready() const{

	bool readyforuse=true;
	readyforuse=readyforuse & snaxs.isready();
	readyforuse=readyforuse & snaxedges.isready();
	readyforuse=readyforuse & snaxsurfs.isready();
	readyforuse=readyforuse & snakeconn.isready();
	readyforuse=readyforuse & snakemesh->isready();

	return(readyforuse);
}

inline void snake::GetMaxIndex(int *nVert,int *nEdge,int *nSurf,int *nVolu) const
{
	snakeconn.GetMaxIndex(nVert,nEdge,nSurf,nVolu);
}

void snake::PrepareForUse(bool needOrder) {

	snaxs.isInMesh=true;
	snaxedges.isInMesh=true;
	snaxsurfs.isInMesh=true;
	
	snaxs.PrepareForUse();
	snaxedges.PrepareForUse();
	snaxsurfs.PrepareForUse();
	snakeconn.PrepareForUse(needOrder);
	if(this->snakemesh!=NULL){
		snakemesh->PrepareForUse();
	} else {
		cerr << endl << "Warning: Mesh containing snake unset in snake"
			" object.\n Set attribute `snakemesh` to a mesh pointer" << endl;
	}

	is3D=int(snakemesh->volus.size())>0;

	if (int(isMeshVertIn.size()) != int(snakemesh->verts.size())){
		isMeshVertIn.assign(snakemesh->verts.size(),false);
	}
}
void snake::HashArray(){
	// prefer use of PrepareForUse in final code as it includes checks to 
	// avoid unnecessary repetiion of work
	// Use HashArray when debugging to be sure that fields are correctly hashed
	snaxs.HashArray();
	snaxedges.HashArray();
	snaxsurfs.HashArray();
	snakeconn.HashArray();
	snakemesh->HashArray();

}
void snake::HashParent(){
	// prefer use of PrepareForUse in final code as it includes checks to 
	// avoid unnecessary repetiion of work
	// Use HashArray when debugging to be sure that fields are correctly hashed
	snaxs.HashParent();
	snaxedges.HashParent();
	snaxsurfs.HashParent();

}
void snake::HashArrayNM(){
	// prefer use of PrepareForUse in final code as it includes checks to 
	// avoid unnecessary repetiion of work
	// Use HashArray when debugging to be sure that fields are correctly hashed
	snaxs.HashArray();
	snaxedges.HashArray();
	snaxsurfs.HashArray();
	snakeconn.HashArray();

}
void snake::SetMaxIndex(){
	// prefer use of PrepareForUse in final code as it includes checks to 
	// avoid unnecessary repetiion of work

	snaxs.SetMaxIndex();
	snaxedges.SetMaxIndex();
	snaxsurfs.SetMaxIndex();
	snakeconn.SetMaxIndex();
	snakemesh->SetMaxIndex();

}
void snake::SetMaxIndexNM(){
	// prefer use of PrepareForUse in final code as it includes checks to 
	// avoid unnecessary repetiion of work

	snaxs.SetMaxIndex();
	snaxedges.SetMaxIndex();
	snaxsurfs.SetMaxIndex();
	snakeconn.SetMaxIndex();

}
void snake::SetLastIndex(){
	// prefer use of PrepareForUse in final code as it includes checks to 
	// avoid unnecessary repetiion of work

	snaxs.SetLastIndex();
	snaxedges.SetLastIndex();
	snaxsurfs.SetLastIndex();
	snakeconn.SetLastIndex();
	snakemesh->SetLastIndex();

}

void snake::Init(mesh *snakemeshin,int nSnax, int nEdge, int nSurf, int nVolu){
// Initialise a 
	is3D=snakemeshin->volus.size()>0;

	snaxs.Init(nSnax);
	snaxedges.Init(nEdge);
	snaxsurfs.Init(nSurf);

	snaxs.isInMesh=true;
	snaxedges.isInMesh=true;
	snaxsurfs.isInMesh=true;

	snakeconn.Init(nSnax, nEdge, nSurf, nVolu);

	snakemesh=snakemeshin;
	isMeshVertIn.assign(snakemeshin->verts.size(),false);
}

void snake::reserve(int nSnax, int nEdge, int nSurf, int nVolu){
// Initialise a 

	snaxs.reserve(nSnax);
	snaxedges.reserve(nEdge);
	snaxsurfs.reserve(nSurf);

	snaxs.isInMesh=true;
	snaxedges.isInMesh=true;
	snaxsurfs.isInMesh=true;

	snakeconn.reserve(nSnax, nEdge, nSurf, nVolu);


}

void snake::ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu){

	snaxs.ChangeIndices(nVert,nEdge,nSurf,nVolu); 
	snaxedges.ChangeIndices(nVert,nEdge,nSurf,nVolu);
	snaxsurfs.ChangeIndices(nVert,nEdge,nSurf,nVolu);
	snakeconn.ChangeIndices(nVert,nEdge,nSurf,nVolu);

}
void snake::ChangeIndicesSnakeMesh(int nVert,int nEdge,int nSurf,int nVolu){
	int ii=0;
	for(ii=0;ii<snaxs.size();++ii){
		snaxs[ii].ChangeIndicesSnakeMesh(nVert,nEdge,nSurf,nVolu); 
	}
	for(ii=0;ii<snaxedges.size();++ii){
		snaxedges[ii].ChangeIndicesSnakeMesh(nVert,nEdge,nSurf,nVolu);
	}
	for(ii=0;ii<snaxsurfs.size();++ii){
		snaxsurfs[ii].ChangeIndicesSnakeMesh(nVert,nEdge,nSurf,nVolu);
	}
	// Nothing to do for snakeconn
	//snakeconn.ChangeIndicesSnakeMesh(nVert,nEdge,nSurf,nVolu);

	snakemesh->ChangeIndices(nVert,nEdge,nSurf,nVolu); //Maybe this one is a bad idea

}

int CompareSnakeInternalStatus(const vector<bool> & thisVec,bool thisFlipped,
	 const vector<bool> & otherVec, bool otherFlipped){
	/*
	Checks if a snake is inside another using information from the points inside
	the snaking mesh.
	Works by finding counter examples. 
	This is not a complete test and assumes a relatively simple relationship between
	the two snakes. of either one completely internal or the two side by side.
	*/
	bool isIndep;
	int ii, ni;
	ni=thisVec.size();
	isIndep=true;
	for (ii=0; ii<ni; ++ii){
		if((thisVec[ii]^thisFlipped) && (otherVec[ii]^otherFlipped)){
			isIndep	= false;
			break;
		}
	}

	return(int(isIndep));
}
void snake::Concatenate(const snake &other, int isInternal){
	/*
	isInternal=0 status unknown
	isInternal=1 other snake is inside this
	isInternal=2 this snake is inside other
	isInternal=-1 other and this are independant
	*/
	int i,ni, isIndep;
	if(this->snakemesh==other.snakemesh){
		this->snaxs.Concatenate(other.snaxs);
		this->snaxedges.Concatenate(other.snaxedges);
		this->snaxsurfs.Concatenate(other.snaxsurfs);
		this->snakeconn.Concatenate(other.snakeconn);
		// Maintain the internal status
		ni=isMeshVertIn.size();

		if(isInternal==0){
			isIndep=CompareSnakeInternalStatus(isMeshVertIn,ReturnFlip(),
				other.isMeshVertIn, other.ReturnFlip());
		} else {
			isIndep=isInternal<=0;
		}
		//cout << "Indep status " << isIndep << isFlipped << other.isFlipped << " " ;
		if(isIndep){
			for (i=0;i<ni;++i){
				isMeshVertIn[i]= ((isMeshVertIn[i] ^ isFlipped)
					|| (other.isMeshVertIn[i] ^ other.isFlipped)) ^ isFlipped;
			}
		} else {
			for (i=0;i<ni;++i){
				isMeshVertIn[i] = ((isMeshVertIn[i] ^ isFlipped)
					&& (other.isMeshVertIn[i] ^ other.isFlipped)) ^ isFlipped;
			}
		}

		// if(!isFlipped && !other.isFlipped){
		// 	for (i=0;i<ni;++i){
		// 		isMeshVertIn[i]= isMeshVertIn[i] || other.isMeshVertIn[i];
		// 	}
		// } else if (isFlipped && other.isFlipped) {
		// 	for (i=0;i<ni;++i){
		// 		isMeshVertIn[i]= isMeshVertIn[i] || other.isMeshVertIn[i];
		// 	}
		// } else if (!isFlipped && other.isFlipped) {
		// 	for (i=0;i<ni;++i){
		// 		isMeshVertIn[i]= isMeshVertIn[i] ^ other.isMeshVertIn[i];
		// 	} //xor
		// } else if (isFlipped && !other.isFlipped) {
		// 	for (i=0;i<ni;++i){
		// 		isMeshVertIn[i]= isMeshVertIn[i] ^ other.isMeshVertIn[i];
		// 	} //xor
		// }
	} else {
		RSVS3D_ERROR_ARGUMENT("Concatenation of snakes with different base meshes not supported");
	}

}


void snake::MakeCompatible_inplace(snake &other) const{
   // Makes other mesh compatible with this to be 
   // merged without index crashes

	int nVert,nEdge,nVolu,nSurf,ii;

   // Define Max indices in current mesh
	this->GetMaxIndex(&nVert,&nEdge,&nSurf,&nVolu);
	other.ChangeIndices(nVert,nEdge,nSurf,nVolu);
	for(ii=0;ii<other.snaxs.size();++ii){
		other.snaxs[ii].orderedge=-1;
	}
}
snake snake::MakeCompatible(snake other) const{
	MakeCompatible_inplace(other);
	return(other);
}


void snake::VertIsIn(int vertInd , bool isIn){

	if (snakemesh!=NULL){
		int i=snakemesh->verts.find(vertInd);
		if (i!=-1){
			isMeshVertIn[i]=isIn;
		}
	} 
}
void snake::VertIsIn(vector<int> vertInd, bool isIn){
	if (snakemesh!=NULL){
		int nj=vertInd.size();
		for (int j=0; j<nj; j++){
			int i=snakemesh->verts.find(vertInd[j]);
			if (i!=-1){
				isMeshVertIn[i]=isIn;
			}
		}
	}
}

void snake::CheckConnectivity() const {

	int ii,jj,ni,kk, ll,errCount, errTot;
	vector<int> testSub;
	bool isErr;

	testSub.reserve(10);

	errCount=0;
	errTot=0;
	ni=int(snaxs.size());
	for (ii=0; ii<ni;++ii){
		if (snaxs.countparent(snaxs(ii)->KeyParent())>1){
			snaxs.findsiblings(snaxs(ii)->KeyParent(), testSub);
			for(jj=0; jj< int(testSub.size());++jj){
				isErr=false;
				if(snaxedges(ii)->index!=snaxedges(testSub[jj])->index){
					for(kk=0; kk< int(snakeconn.verts(ii)->edgeind.size());kk++){
						if((snakeconn.edges.isearch(
							snakeconn.verts(ii)->edgeind[kk]
							)->vertind[0]==snaxs(testSub[jj])->index)
							|| (snakeconn.edges.isearch(
							snakeconn.verts(ii)->edgeind[kk]
							)->vertind[1]==snaxs(testSub[jj])->index)){
							isErr=true;
							break;
						}
					}
				}
				if (isErr){
					cerr << endl << " Test Snake Connectivity Error :" 
						<< errCount << " snaxel " << snakeconn.verts(ii)->index
						<< " is connected to snaxel " << snaxs(testSub[jj])->index 
						<< " list: " << endl; 
					snakeconn.verts(ii)->disp();
					snakeconn.verts(testSub[jj])->disp();
					errCount++;
				}
			}
		}
	}
	if (errCount>0){
		cerr << "Test Connectivity snax (same edge connectivity) Errors :" 
			<< errCount << endl;
	}
	errTot+=errCount;
	errCount=0;
	ni=int(snaxs.size());
	/*Snaxel does not have the same number of edges as its parent edge
	has surface connections*/
	for (ii=0; ii<ni;++ii){
		if (snakeconn.verts(ii)->edgeind.size()
				!=snakemesh->edges.isearch(
					snaxs(ii)->edgeind
					)->surfind.size()){
			cerr << endl << " Test snax number of edges error :" << errCount 
				 << " snaxel-vertex " << snakeconn.verts(ii)->index
				 << " is connected to snaxel " << snaxs(ii)->edgeind 
				 << " list: " << endl; 
			snakeconn.verts(ii)->disp();
			snakemesh->edges.isearch(snaxs(ii)->edgeind)->disp();
			errCount++;
		}
	}
	if (errCount>0){
		cerr << "Test Connectivity snax (edges in surfs) Errors :"
			<< errCount << endl;
	}
	errTot+=errCount;


	errCount=0;
	errTot=0;
	ni=int(snaxedges.size());
	for (ii=0; ii<ni;++ii){
		if (snaxedges.countparent(snaxedges(ii)->KeyParent())>1){
			snaxedges.findsiblings(snaxedges(ii)->KeyParent(), testSub);
			for(jj=0; jj< int(testSub.size());++jj){
				isErr=false;
				if(snaxedges(ii)->index!=snaxedges(testSub[jj])->index){
					for(kk=0; kk< int(snakeconn.edges(ii)->vertind.size());kk++){
						if((snakeconn.edges(ii)->vertind[kk]
								==snakeconn.edges(testSub[jj])->vertind[0])
							|| 
							(snakeconn.edges(ii)->vertind[kk]
								==snakeconn.edges(testSub[jj])->vertind[1])){
							isErr=true;
							break;
						}
					}
				}
				if (isErr){
					cerr << endl << " Test SnakeEdges Connectivity Error :" << errCount << " snaxedge " << snakeconn.edges(ii)->index
					<< " is connected to snaxedge " << snaxedges(testSub[jj])->index << " list: " << endl; 
					snakeconn.edges(ii)->disp();
					snakeconn.edges(testSub[jj])->disp();
					errCount++;
				}
			}
		}
	}
	if (errCount>0){
		cerr << "Test Connectivity snaxedge (same surf connectivity) Errors :" << errCount << endl;
	}
	errTot+=errCount;


	errCount=0;
	errTot=0;
	ni=int(snaxsurfs.size());
	for (ii=0; ii<ni;++ii){
		if (snaxsurfs.countparent(snaxsurfs(ii)->KeyParent())>1){
			snaxsurfs.findsiblings(snaxsurfs(ii)->KeyParent(), testSub);
			for(jj=0; jj< int(testSub.size());++jj){
				isErr=false;
				if(snaxsurfs(ii)->index!=snaxsurfs(testSub[jj])->index){
					for(kk=0; kk< int(snakeconn.surfs(ii)->edgeind.size());kk++){
						for(ll=0; ll<int(snakeconn.surfs(testSub[jj])->edgeind.size());++ll){
							if(snakeconn.surfs(testSub[jj])->edgeind[ll]==snakeconn.surfs(ii)->edgeind[kk]){
								isErr=true;
								break;
							}
						}
						if(isErr){
							break;
						}
					}
				}
				if (isErr){
					cerr << endl << " Test SnakeSurf Connectivity Error :" << errCount << " snaxsurf " << snakeconn.surfs(ii)->index
					<< " is connected to snaxsurf " << snaxsurfs(testSub[jj])->index << " list: " << endl; 
					snakeconn.surfs(ii)->disp();
					snakeconn.surfs(testSub[jj])->disp();
					errCount++;
				}
			}
		}
	}
	if (errCount>0){
		cerr << "Test Connectivity snaxsurf (same volu connectivity) Errors :" << errCount << endl;
	}
	errTot+=errCount;

	if (errTot>0){
		cerr << errTot << "  Total errors were detected in the connectivity list" <<endl;
	}

}

void snake::OrderEdges(){

	int nNewSurfs;
	snaxsurf newSurf;
	vector<int> newSurfPairs;


	newSurfPairs=snakeconn.OrderEdges();
	nNewSurfs = int(newSurfPairs.size()/2);
	for (int i = 0; i < nNewSurfs; ++i)
	{
		newSurf=*(snaxsurfs.isearch(newSurfPairs[i*2]));
		newSurf.index = newSurfPairs[i*2+1];
		snaxsurfs.push_back(newSurf);
	}


}

// Snake Movement
void snake::UpdateDistance(double dt, double maxDstep){
	int ii;
	double dStep;
	for (ii = 0; ii < int(snaxs.size()); ++ii)
	{
		if(snaxs[ii].v*dt>maxDstep){
			dStep = maxDstep;
		} else if (snaxs[ii].v*dt<-maxDstep){
			dStep = -maxDstep;
		} else{
			dStep = snaxs[ii].v*dt; 
		}
		snaxs[ii].d=snaxs[ii].d+dStep;
		snaxs[ii].d=snaxs[ii].d<0.0 ? 0.0 : (snaxs[ii].d>1.0 ? 1.0 : snaxs[ii].d);
		
	}
}
void snake::UpdateDistance(const vector<double> &dt,  double maxDstep){
	int ii;
	double dStep;

	if (int(dt.size())==snaxs.size()){
		for (ii = 0; ii < int(snaxs.size()); ++ii)
		{
			if(snaxs[ii].v*dt[ii]>maxDstep){
			dStep = maxDstep;
		} else if (snaxs[ii].v*dt[ii]<-maxDstep){
			dStep = -maxDstep;
		} else{
			dStep = snaxs[ii].v*dt[ii]; 
		}
		snaxs[ii].d=snaxs[ii].d+dStep;
			snaxs[ii].d=snaxs[ii].d<0.0 ? 0.0 : (snaxs[ii].d>1.0 ? 1.0 : snaxs[ii].d);

		}
	} else {
		cerr << "Error : Mismatching sizes passed to Uptdate distance of snake"  << endl;
		cerr << "Error in " << __PRETTY_FUNCTION__ << endl;
		throw length_error("Snake and Dt had mismatched sizes at snake::UpdateDistance");
	}
}

void snake::UpdateCoord(){
	int fromVertSub,toVertSub;
	int ii,jj;
	for (ii = 0; ii < int(snaxs.size()); ++ii)
	{
		fromVertSub=snakemesh->verts.find(snaxs(ii)->fromvert);
		toVertSub=snakemesh->verts.find(snaxs(ii)->tovert);

		for(jj=0;jj<3;++jj){
			snakeconn.verts.elems[ii].coord[jj]=(snakemesh->verts(toVertSub)->coord[jj]
				-snakemesh->verts(fromVertSub)->coord[jj])*snaxs(ii)->d
			+snakemesh->verts(fromVertSub)->coord[jj];
		}
	}
}

// snaxarray extension


void snaxarray::OrderOnEdge(){

	vector<int> edgeInds,snaxSubs,orderSubs;
	vector<double> orderSnax,unorderSnax;
	unordered_multimap<double,int> hashOrder;
	int isBwd;
	int ii, jj, nEdge;

	edgeInds=ConcatenateScalarField(*this,&snax::edgeind,0,elems.size());
	sort(edgeInds);
	unique(edgeInds);
	for(ii=0;ii<int(edgeInds.size());++ii){
		nEdge=hashParent.count(edgeInds[ii]);

		if (nEdge>1){
			snaxSubs=ReturnDataEqualRange(edgeInds[ii], hashParent);
			orderSnax.resize(snaxSubs.size());
			for(jj=0;jj<int(snaxSubs.size());++jj){
				isBwd=elems[snaxSubs[jj]].fromvert>elems[snaxSubs[jj]].tovert;
				orderSnax[jj]= isBwd ? (1.0-elems[snaxSubs[jj]].d) : elems[snaxSubs[jj]].d ;
			}
			unorderSnax=orderSnax;
			// This will break with snakes in same place
			sort(orderSnax);

			orderSubs=FindSubList(orderSnax,unorderSnax,hashOrder);
			for(jj=0;jj<int(snaxSubs.size());++jj){
				elems[snaxSubs[orderSubs[jj]]].orderedge=jj+1;
			}
		}

	}

	isOrderedOnEdge=1;
}


void snaxarray::ReorderOnEdge()
{

	vector<int> edgeInds,snaxSubs;
	vector<double> orderSnax;
	int isBwd;
	int ii, jj, nEdge, maxOrd=0;
	bool needIncrement;

	edgeInds=ConcatenateScalarField(*this,&snax::edgeind,0,elems.size());
	sort(edgeInds);
	unique(edgeInds);
	for(ii=0;ii<int(edgeInds.size());++ii){

		nEdge=hashParent.count(edgeInds[ii]);

		if (nEdge>1){
			snaxSubs=ReturnDataEqualRange(edgeInds[ii], hashParent);
			orderSnax.resize(snaxSubs.size());
			needIncrement=false;
			maxOrd=0;
			for(jj=0;jj<int(snaxSubs.size());++jj){

				maxOrd=(maxOrd<=elems[snaxSubs[jj]].orderedge) ? elems[snaxSubs[jj]].orderedge : maxOrd;
			}
			for(jj=0;jj<int(snaxSubs.size());++jj){
				if (elems[snaxSubs[jj]].orderedge==-1){
					isBwd=elems[snaxSubs[jj]].fromvert>elems[snaxSubs[jj]].tovert;

					elems[snaxSubs[jj]].orderedge= (isBwd && (elems[snaxSubs[jj]].d==0.0)) ? maxOrd+1 
					:  ((!isBwd && (elems[snaxSubs[jj]].d==0.0)) ? 0 
						: ((isBwd && (elems[snaxSubs[jj]].d==1.0)) ? 0 
							:  ((!isBwd && (elems[snaxSubs[jj]].d==1.0)) ? maxOrd+1 
								: elems[snaxSubs[jj]].orderedge ))); 
					if (elems[snaxSubs[jj]].orderedge==-1){
						cerr << "Error Order not correctly assigned" << endl;
					}
					if (elems[snaxSubs[jj]].orderedge==0){
						needIncrement=true;
					}
				}
			}
			if(needIncrement){
				for(jj=0;jj<int(snaxSubs.size());++jj){
					elems[snaxSubs[jj]].orderedge++;
				}
			}
		} else {
			snaxSubs=ReturnDataEqualRange(edgeInds[ii], hashParent);
			elems[snaxSubs[0]].orderedge=1;
		}

	}

	isOrderedOnEdge=1;
}

void snaxarray::CalculateTimeStepOnEdge(vector<double> &dt, 
	vector<bool> &isSnaxDone, int edgeInd){

	int nSnax,ii,jj;
	vector<int> snaxSubs;
	double impactTime;

	snaxSubs=ReturnDataEqualRange(edgeInd, hashParent);
	nSnax=snaxSubs.size();
	for(ii=0;ii<nSnax;++ii){
		isSnaxDone[snaxSubs[ii]]=true;
		for(jj=ii+1;jj<nSnax;++jj){
			impactTime=SnaxImpactDt(elems[snaxSubs[ii]],elems[snaxSubs[jj]]);
			if (impactTime>=0.0){
				dt[snaxSubs[ii]]= (dt[snaxSubs[ii]]>impactTime) ?
					impactTime : dt[snaxSubs[ii]];
				dt[snaxSubs[jj]]= (dt[snaxSubs[jj]]>impactTime) ?
					impactTime : dt[snaxSubs[jj]];
			}
		}
	}
}

double ValidateEquilibrateD(double newD, double oldD, double stepLength){

	double maxLim = pow(stepLength,1.5);
	double minLim = pow(stepLength,0.5);
	if (oldD<0.5){
		if (!isfinite(newD)){
			newD = oldD;
		} else if (newD<minLim){
			newD = minLim;
		} else if (newD>maxLim) {
			newD = maxLim;
		}
	} else {

		if (!isfinite(newD)){
			newD = oldD;
		} else if (newD<(1.0-maxLim)){
			newD = (1.0-maxLim);
		} else if (newD>(1.0-minLim)) {
			newD = (1.0-minLim);
		}
	}

	return newD;
}

std::pair<double,double> EquilibrateEdge(snake &snakein,
	int snaxA, int snaxB, int vertA, int vertB){

	// calculate centre of edge.
	coordvec c, del1, vec;
	double dg0;
	c = snakein.snakeconn.verts.isearch(snaxA)->coord;
	c.add(snakein.snakeconn.verts.isearch(snaxB)->coord);
	c.div(2.0);

	del1 = snakein.snakeconn.verts.isearch(vertA)->coord;
	c.substract(snakein.snakeconn.verts.isearch(vertB)->coord);
	del1.substract(snakein.snakeconn.verts.isearch(vertB)->coord);
	vec = del1.cross(c.usedata());

	// c is the plane normal at this point;
	double d = vec.dot(snakein.snakeconn.verts.isearch(vertA)->coord);

	// Vertex A
	del1 = snakein.snakemesh->verts.isearch(snakein.snaxs.isearch(snaxA)->tovert)->coord;
	del1.substract(snakein.snakemesh->verts.isearch(
		snakein.snaxs.isearch(snaxA)->fromvert)->coord);
	dg0 = vec.dot(snakein.snakemesh->verts.isearch(
		snakein.snaxs.isearch(snaxA)->fromvert)->coord);
	double dga = vec.dot(del1.usedata());
	double da = (d-dg0)/dga;



	// Vertex B
	dg0 = vec.dot(snakein.snakemesh->verts.isearch(
		snakein.snaxs.isearch(snaxB)->fromvert)->coord);
	del1 = snakein.snakemesh->verts.isearch(snakein.snaxs.isearch(snaxB)->tovert)->coord;
	del1.substract(snakein.snakemesh->verts.isearch(
		snakein.snaxs.isearch(snaxB)->fromvert)->coord);
	double dgb = vec.dot(del1.usedata());
	double db = (d-dg0)/dgb;
	


	return {da, db};
}
// template <class T=int, int n=2>
void DisplayVectorArray(const std::vector< std::array<int,2>> &in){
	int count = in.size();
	std::cout << count << ",";
	if(!count){
		return;
	}
	int count2 = in[0].size();
	std::cout << count2 << " : ";
	for (int i = 0; i < count; ++i)
	{
		for (int j = 0; j < count2; ++j)
		{
			std::cout << in[i][j] << " ";
		}
		std::cout << " : ";
	}
	std::cout << " | ";
}

int SmoothStep3(int spawnVert, snake &snakein, double spawnDist){

	// identify snaxels

	HashedVector<int, int> snaxInds;
	double spawnDistTol = spawnDist * 1.1;
	snaxInds.reserve(20); 
	for (auto parentEdge : snakein.snakemesh->verts.isearch(spawnVert)->edgeind)
	{
		std::vector<int> snaxIndsTemp;
		snaxIndsTemp.reserve(20); 
		snakein.snaxs.findsiblings(parentEdge,snaxIndsTemp);

		for (auto snaxCandidate : snaxIndsTemp){
			auto currSnax = snakein.snaxs(snaxCandidate);
			if(currSnax->CloseToVertex()==spawnVert 
				&& snaxInds.find(currSnax->index)==rsvs3d::constants::__notfound
				&& (currSnax->d<spawnDistTol || (1-spawnDistTol)<currSnax->d)){
				snaxInds.push_back(currSnax->index);
			}
		}
	}
	if (snaxInds.size()<2){
		return 0;
	}
	HashedVector<int,int> borderSurfs;
	for (auto iSnax : snaxInds.vec){
		for(auto iSurf : snakein.snakeconn.verts.isearch(iSnax)
			->elmind(snakein.snakeconn,2)){
			bool isSpawned = false;
			for (auto iVert : snakein.snakeconn.surfs.isearch(iSurf)
				->vertind(snakein.snakeconn)){
				if(iVert==rsvs3d::constants::__notfound){
					isSpawned = true; 
					break;
				}
			}
			if (!isSpawned 
				&& borderSurfs.find(iSurf)==rsvs3d::constants::__notfound){
				borderSurfs.push_back(iSurf);
			}
		}
	}


	// Consistancy of this list can be checked making sure all the snaxels come
	// from the same vertex or go to the same vertex.
	HashedVector<int,int> externalPts;
	std::vector<std::array<int,2>> externalLinks;
	externalPts.reserve(snaxInds.size());
	externalLinks.reserve(snaxInds.size());
	for (auto iSnax : snaxInds.vec){
		for (auto iEdge : snakein.snakeconn.verts.isearch(iSnax)->edgeind){
			for (auto iVert : snakein.snakeconn.edges.isearch(iEdge)->vertind){
				if(snaxInds.find(iVert)==rsvs3d::constants::__notfound){
					externalPts.push_back(iSnax);
					externalLinks.push_back({iEdge,iVert});
				}
			}
		}
	}
	// Compute the two distances
	std::vector<int> edgeProcess;
	edgeProcess.reserve(borderSurfs.size());
	for (auto border : borderSurfs.vec){
		bool success = false;
		for (auto iEdge : snakein.snakeconn.surfs.isearch(border)->edgeind){
			auto& vertTemp = snakein.snakeconn.edges.isearch(iEdge)->vertind;
			if(snaxInds.find(vertTemp[0])!=rsvs3d::constants::__notfound
				&& snaxInds.find(vertTemp[1])!=rsvs3d::constants::__notfound
				&& success){
				success = false;
				edgeProcess.pop_back();
				break;
			} else if(snaxInds.find(vertTemp[0])!=rsvs3d::constants::__notfound
				&& snaxInds.find(vertTemp[1])!=rsvs3d::constants::__notfound){
				success = true;
				edgeProcess.push_back(iEdge);
			}
		}
		if(!success){
			edgeProcess.push_back(rsvs3d::constants::__notfound);
			// RSVS3D_ERROR_LOGIC("A border edge was not found when looking for it.");
		}
	}

	std::vector<double>  newD;
	std::vector<int>  countAccess;
	newD.assign(snaxInds.size(), 0.0);
	countAccess.assign(snaxInds.size(),0);
	int count1 = edgeProcess.size();
	for(int i=0; i< count1; ++i){
		int actSurf = borderSurfs[i];
		int iEdge = edgeProcess[i];
		if(iEdge==rsvs3d::constants::__notfound){
			continue;
		}
		// pass the edge index, the two snaxels, the two vertices
		int snaxA = snakein.snakeconn.edges.isearch(iEdge)->vertind[0];
		int snaxB = snakein.snakeconn.edges.isearch(iEdge)->vertind[1];

		int posA = rsvs3d::constants::__notfound;
		int posB = rsvs3d::constants::__notfound;
		bool success = false;
		
		for (auto iA : externalPts.findall(snaxA)){
			// Find the one where they all share the same surface
			posA = iA;
			auto testEdge = snakein.snakeconn.edges.isearch(
				externalLinks[posA][0]);
			if(actSurf == testEdge->surfind[0] || actSurf == testEdge->surfind[1]){
				success = true;
				break;
			}
		}
		if(!success){
			std::cerr << snaxA << " " << posA << " | ";
			DisplayVector(externalPts.vec);
			DisplayVector(borderSurfs.vec);
			DisplayVector(edgeProcess);
			DisplayVectorArray(externalLinks);
			std::cerr << std::endl;
			snakein.snakeconn.edges.isearch(iEdge)->disp();
			snakein.snakeconn.surfs.isearch(actSurf)->disp();
			snakein.snakeconn.verts.isearch(snaxA)->disp();
			RSVS3D_ERROR_LOGIC("Pos A could not be found");
		}
		success = false;
		for (auto iA : externalPts.findall(snaxB)){
			// Find the one where they all share the same surface
			posB = iA;
			if (borderSurfs.find(actSurf)!=rsvs3d::constants::__notfound){
				auto testEdge = snakein.snakeconn.edges.isearch(
					externalLinks[posA][0]);
				if(actSurf == testEdge->surfind[0] 
					|| actSurf == testEdge->surfind[1]){
					success = true;
					break;
				}
			}
		}
		if(!success){
			RSVS3D_ERROR_LOGIC("Pos B could not be found");
		}

		int vertA = externalLinks[posA][1];
		int vertB = externalLinks[posB][1];
		auto tempD = EquilibrateEdge(snakein, snaxA, snaxB, vertA, vertB);
		posA = snaxInds.find(snaxA);
		posB = snaxInds.find(snaxB);
		ValidateEquilibrateD(tempD.first, snakein.snaxs.isearch(snaxA)->d, spawnDist);
		// newD[posA] += tempD.first;
		newD[posA] += ValidateEquilibrateD(tempD.first, 
			snakein.snaxs.isearch(snaxA)->d, spawnDist);
		countAccess[posA]++;
		// newD[posB] += tempD.second;
		newD[posB] += ValidateEquilibrateD(tempD.second, 
			snakein.snaxs.isearch(snaxB)->d, spawnDist);
		countAccess[posB]++;
	}
	// How to set all the values not set by the border algorithm?



	// Set the values
	int count = snaxInds.size();
	int noissue=0;
	for (int i = 0; i < count; ++i)
	{
		if(countAccess[i]>0){
			if((newD[i]<0.0 || countAccess[i]!=2 || newD[i]/countAccess[i]>1.0)
				&& false){
				DisplayVector(newD);
				DisplayVector(countAccess);
				std::cout << endl;
				DisplayVector(externalPts.vec);
				DisplayVector(snaxInds.vec);
				DisplayVector(borderSurfs.vec);
				DisplayVector(edgeProcess);

				RSVS3D_ERROR_NOTHROW("Negative d was calculated");
			} else {
				noissue++;
				snakein.snaxs[snakein.snaxs.find(snaxInds[i], true)].d 
					= newD[i]/double(countAccess[i]);
			}
		}
	}
	return(noissue);
}


/**
 * @brief      Calculates the vertex normal weighted by surface angle partitions
 *
 * @param[in]  centre  The coordinate at which the normal needs to be evaluated.
 * @param[in]  vecPts  The points used to compute the normal.
 *
 * @return     A tuple with the normal vector and the total angle between the 
 * surfaces. This angle is a measure of the local curvature.
 */
std::tuple<coordvec,double> VertexNormal(const std::vector<double>& centre, 
	const std::vector<const std::vector<double>*> &vecPts){

	// Calculates a normal for each face
	// and an angle
	// Need to make sure it is inward pointing
	// 

	if(vecPts.size()==0){
		RSVS3D_ERROR_ARGUMENT("Attempted to define a smooth step for an empty"
			"vector of points.");
	}

	coordvec planeCurr, planePrev, normal, temp;
	double totalNormalAngle=0.0;
	double totalTangentAngle=0.0;
	double currAngle;
	int count = vecPts.size();
	currAngle = PlaneNormalAndAngle(centre, *vecPts.back(), *vecPts[0], planeCurr, temp);
	planeCurr.Normalize();
	for (int i = 0; i < count; ++i)
	{
		planePrev.swap(planeCurr);
		currAngle = PlaneNormalAndAngle(centre, *vecPts[i], *vecPts[(i+1)%count],
			planeCurr, temp);
		planeCurr.Normalize();
		totalNormalAngle += currAngle;
		planeCurr.mult(currAngle);
		normal.add(planeCurr.usedata());
		planeCurr.Normalize();
		totalTangentAngle += planeCurr.angle(planePrev);
	}
	normal.div(totalNormalAngle);
	return std::make_tuple(normal,totalTangentAngle);
}


/**
 * @brief      Calculates the parametric position along a line at which it
 *             intersects a sphere
 *
 * The return value of this function when multiplied by the line direction and
 * added to the reference point gives the coordinate of the intersection points.
 *
 * @param[in]  lineVec       The vector direction followed by the line
 * @param[in]  offset        The vector going from the sphere centre to a point
 *                           on the line.
 * @param[in]  sphereRadius  The radius of the sphere.
 *
 * @return     pair of quadratic roots. Returns infinity if there is not intersection.
 */
std::array<double, 2> IntersectLineSphere(const coordvec &lineVec, 
	const coordvec &offset, double sphereRadius){

	double a, b, c;
	a = lineVec.dot(lineVec.usedata());
	b = 2 * lineVec.dot(offset.usedata());
	c = offset.dot(offset.usedata()) - (sphereRadius*sphereRadius);

	double det = b*b - (4 * a*c);
	if (det<0.0){
		return {-INFINITY, INFINITY};
	}
	if(IsAproxEqual(det, 0.0)){
		double sol = -b/(2*a);
		return {sol, sol};
	}
	det = sqrt(det);
	return {(-b-det)/(2*a), (-b+det)/(2*a)};
}


std::array<double, 2> IntersectLineSphere(const coordvec &lineVec,
	const vector<double> &linePoint, const coordvec &sphereCentre,
	double sphereRadius){

	static coordvec offset;
	offset = linePoint;
	offset.substract(sphereCentre.usedata());
	return IntersectLineSphere(lineVec, offset,sphereRadius);
}

/**
 * @brief      Smooths the step away from a spawn vertex
 * 
 * THis function smooths the spawn step by placing the vertices on a sphere
 * corresponding to a small progression of the inset sphere at this point.
 *
 * @param[in]  spawnVert  The spawn vertex in the snaking mesh
 * @param      snakein    The snake which needs smoothing
 * @param[in]  spawnDist  The spawn step length
 *
 * @return     The number of succesfully smoothed vertex steps.
 */
int SmoothStep(int spawnVert, snake &snakein, double spawnDist){

	// identify snaxels

	if(IsAproxEqual(spawnDist,0.0)){
		return  0;
	}

	HashedVector<int, int> snaxInds;
	double spawnDistTol = spawnDist * 1.1;
	snaxInds.reserve(20); 
	for (auto parentEdge : snakein.snakemesh->verts.isearch(spawnVert)->edgeind)
	{
		std::vector<int> snaxIndsTemp;
		snaxIndsTemp.reserve(20); 
		snakein.snaxs.findsiblings(parentEdge,snaxIndsTemp);

		for (auto snaxCandidate : snaxIndsTemp){
			auto currSnax = snakein.snaxs(snaxCandidate);
			if(currSnax->CloseToVertex()==spawnVert 
				&& snaxInds.find(currSnax->index)==rsvs3d::constants::__notfound
				&& (currSnax->d<spawnDistTol || (1-spawnDistTol)<currSnax->d)){
				snaxInds.push_back(currSnax->index);
			}
		}
	}
	if (snaxInds.size()<2){
		return 0;
	}

	// Consistancy of this list can be checked making sure all the snaxels come
	// from the same vertex or go to the same vertex.

	std::vector<int> externalEdges, externalPts, externalSurfsPairs, 
		externalEdgesOrder, internalPts, internalPtsOrdered;
	externalPts.reserve(snaxInds.size());
	externalEdges.reserve(snaxInds.size());
	internalPts.reserve(snaxInds.size());
	externalSurfsPairs.reserve(snaxInds.size()*2);
	int k = 0;
	for (auto iSnax : snaxInds.vec){
		for (auto iEdge : snakein.snakeconn.verts.isearch(iSnax)->edgeind){
			auto currEdge = snakein.snakeconn.edges.isearch(iEdge);
			for (auto iVert : currEdge->vertind){
				if(snaxInds.find(iVert)==rsvs3d::constants::__notfound){
					internalPts.push_back(iSnax);
					externalEdges.push_back(iEdge);
					externalPts.push_back(iVert);
					externalSurfsPairs.push_back(currEdge->surfind[0]);
					externalSurfsPairs.push_back(currEdge->surfind[1]);
					externalEdgesOrder.push_back(k++);
				}
			}
		}
	}
	int retValue = OrderList(externalEdgesOrder, externalSurfsPairs, false, true);
	int nVerts = externalEdgesOrder.size();
	if(retValue != rsvs3d::constants::ordering::ordered){
		RSVS3D_ERROR_LOGIC("Failed to order the list of edges.");
	}
	if(nVerts<2){
		RSVS3D_ERROR_LOGIC("Only one point was left, cannot define direction.");
	}
	std::vector<const std::vector<double>*> vecPts;
	vecPts.reserve(nVerts);
	internalPtsOrdered.reserve(nVerts);
	for (int i = 0; i < nVerts; ++i)
	{
		vecPts.push_back(
			&(snakein.snakeconn.verts.isearch(
				externalPts[externalEdgesOrder[i]])->coord)
			);
	}
	auto& centre = snakein.snakemesh->verts.isearch(spawnVert)->coord;
	auto centreData = VertexNormal(centre, vecPts);
	coordvec& vertNormal = get<0>(centreData);
	// Figure out if the normal is inward facing or outward facing.
	// Order the internal edges
	internalPtsOrdered.reserve(nVerts);
	int v1 = rsvs3d::constants::__notfound;
	int v2 = rsvs3d::constants::__notfound;
	for (int i = 0; i < nVerts; ++i){
		internalPtsOrdered.push_back(internalPts[externalEdgesOrder[i]]);
	}
	v1 = internalPtsOrdered[0];
	for (int i = 1; i < nVerts; ++i)
	{
		if (v1 != internalPtsOrdered[i]){
			v2 = internalPtsOrdered[i];
			break;
		}
	}
	if(v2==rsvs3d::constants::__notfound){
		RSVS3D_ERROR_LOGIC("Only one point was found, cannot define direction.");
	}
	// Extract a surface that is between two points (on an edge)
	int edgeInd = snakein.snakeconn.EdgeFromVerts(v1, v2);
	if(edgeInd==rsvs3d::constants::__notfound){
		RSVS3D_ERROR_LOGIC("The two points do not share an edge.");
	}
	// Compare the list orders
	unique(internalPtsOrdered);
	if(internalPtsOrdered.back() == internalPtsOrdered[0]){
		internalPtsOrdered.pop_back();
	}
	if(internalPtsOrdered.size()<2){
		RSVS3D_ERROR_LOGIC("Only one point was left, cannot define direction.");
	}
	int surfInd = snakein.snakeconn.edges.isearch(edgeInd)->surfind[0];
	auto& edgeindSurf = snakein.snakeconn.surfs.isearch(surfInd)->edgeind;
	int nEdges =  edgeindSurf.size();
	int correctSurf = 0;
	for (int i = 0; i < nEdges; ++i)
 	{
 		bool breakOut = false;
		for (auto iVert : snakein.snakeconn.edges.isearch(edgeindSurf[i])->vertind){
			if (iVert!=v1 && iVert!=v2 
				&& snaxInds.find(iVert)==rsvs3d::constants::__notfound){
				correctSurf++;
				breakOut = true;
				break;
			}
		}
		if (breakOut){break;}
	}
	surfInd = snakein.snakeconn.edges.isearch(edgeInd)->surfind[correctSurf];
	snakein.snakeconn.surfs.isearch(surfInd)->OrderedVerts(&snakein.snakeconn, internalPts);
	int isSameOrder = OrderMatchLists(internalPtsOrdered, internalPts, v1, v2);
	// compares the list vec1 and vec2 returning 
	// 1 if indices p1 and p2 appear in the same order 
	// -1 if indices p1 and p2 appear in opposite orders
	int flipMultiplier = 0;
	if ((snakein.snakeconn.surfs.isearch(surfInd)->voluind[1]==0)
		&& (isSameOrder==-1)){ // case 1 right way if pointing through
		flipMultiplier = 1;
	} else if ((snakein.snakeconn.surfs.isearch(surfInd)->voluind[0]==0)
		&& (isSameOrder==1)) { // case 2 right way
		flipMultiplier = 1;
	} else if ((snakein.snakeconn.surfs.isearch(surfInd)->voluind[1]==0)
		&& (isSameOrder==1)) { // case 3 wrong way
		flipMultiplier = -1;
	} else if ((snakein.snakeconn.surfs.isearch(surfInd)->voluind[0]==0)
		&& (isSameOrder==-1)) { // case 4 wrong way
		flipMultiplier = -1;
	}
	if (flipMultiplier==0){
		std::cerr << std::endl;
		snakein.snakeconn.surfs.isearch(surfInd)->disp();
		std::cerr << " isSameOrder " << isSameOrder << std::endl;
		RSVS3D_ERROR_LOGIC("Flip multiplier was not set in the cases.");
	}

	if (flipMultiplier==1){ // Vector needs to be inward facing
		vertNormal.flipsign();
	}
	vertNormal.Normalize();

	/**
	At this point we have an outwards facing normal and a "measure" of 
	curvature. From this we define a sphere 
	if the angle is 0 the sphere is infinitely big
	if the angle is 360deg the sphere is of radius 0

	Bound these metrics 
		max offset: longest edge
		min offset: 0
	This will position the sphere centre as
		centre = spawnVert - (offset*normal)
	The sphere will have radius
		radius = offset+(spawnStep*shortestEdge(snakeMesh)). 
	*/
	int radiusStep = 1;
	if (!snakein.isMeshVertIn[snakein.snakemesh->verts.find(spawnVert)]){
		radiusStep = -1;
	}
	double maxOffset = 0.0;
	for (auto iEdge : externalEdges){
		double maxOffsetTemp = snakein.snakeconn.edges.isearch(iEdge)
			->Length(snakein.snakeconn);
		maxOffset = max(maxOffsetTemp, maxOffset);
	}
	double minEdgeLength = 0.0;
	for (auto iSnax : snaxInds.vec){
		int iEdge = snakein.snaxs.isearch(iSnax)->edgeind;
		double minTemp = snakein.snakemesh->edges.isearch(iEdge)
			->Length(*snakein.snakemesh);
		minEdgeLength = max(minTemp, minEdgeLength);
	}

	double offsetMultiplier = 1.0 - (get<1>(centreData)/ (2.0 * M_PI));
	offsetMultiplier = max(offsetMultiplier, 0.0);
	double offset = offsetMultiplier*maxOffset; 
	vertNormal.mult(offset);
	// vertNormal.add(centre); // This is now the sphere centre
	/**
	 * For each vertex in snaxInds compute the equation below, lineVec needs to 
	 * be defined from the spawnVert. Then some ds will need to be flipped.
	 * 
	 * the offset  
	 */
	// std::array<double, 2> IntersectLineSphere(const coordvec &lineVec, 
	// 	vertNormal, sphereRadius)
	double sphereRadius = offset + radiusStep*minEdgeLength * spawnDist;
	std::vector<double> newD;
	coordvec lineVec;
	newD.reserve(snaxInds.vec.size());
	for (auto iSnax : snaxInds.vec){
		lineVec = centre;
		int farVert = snakein.snaxs.isearch(iSnax)->CloseToVertex(true);
		lineVec.substractfrom(snakein.snakemesh->verts.isearch(farVert)->coord);
		auto ds = IntersectLineSphere(lineVec, vertNormal, sphereRadius);
		#ifdef RSVS_DIAGNOSTIC
		std::cerr << "ds : " << ds[0] << "  " << ds[1] << std::endl;
		#endif
		if (radiusStep==-1 && rsvs3d::sign(ds[0])==1
			&& rsvs3d::sign(ds[1])==1){
			newD.push_back(min(ds[0], ds[1]));
		} else if ((rsvs3d::sign(ds[0])==rsvs3d::sign(ds[1]))){
			std::cerr << "ds : " << ds[0] << "  " << ds[1] << std::endl;
			RSVS3D_ERROR_LOGIC("The intersection should return a positive and"
				"negative number.");
		} else {
			int dUse = (rsvs3d::sign(ds[1])+1)/2;
			newD.push_back(ds[dUse]);
		}
	}
	#ifdef RSVS_DIAGNOSTIC
	std::cerr << radiusStep << std::endl;
	DisplayVector(newD);
	std::cerr  << std::endl;
	#endif
	int noissue=0;
	int count = snaxInds.vec.size();
	for (int i = 0; i < count; ++i)
	{
		noissue++;
		int j = snakein.snaxs.find(snaxInds[i], true);
		snakein.snaxs[j].d = snakein.snaxs(j)->fromvert==spawnVert ?
			newD[i] : 1.0-newD[i];
	}

	return(noissue);
}

void SmoothStep2(int spawnVert, snake &snakein, double spawnDist){

	// identify snaxels

	HashedVector<int, int> snaxInds;
	double spawnDistTol = spawnDist * 1.1;
	snaxInds.reserve(20); 
	for (auto parentEdge : snakein.snakemesh->verts.isearch(spawnVert)->edgeind)
	{
		std::vector<int> snaxIndsTemp;
		snaxIndsTemp.reserve(20); 
		snakein.snaxs.findsiblings(parentEdge,snaxIndsTemp);

		for (auto snaxCandidate : snaxIndsTemp){
			auto currSnax = snakein.snaxs(snaxCandidate);
			if(currSnax->CloseToVertex()==spawnVert 
				&& snaxInds.find(currSnax->index)==rsvs3d::constants::__notfound
				&& (currSnax->d<spawnDistTol || (1-spawnDistTol)<currSnax->d)){
				snaxInds.push_back(currSnax->index);
			}
		}
	}
	if (snaxInds.size()<2){
		return;
	}
	// Consistancy of this list can be checked making sure all the snaxels come
	// from the same vertex or go to the same vertex.
	HashedVector<int,int> externalPts;
	std::vector<std::array<int,2>> externalLinks;
	externalPts.reserve(snaxInds.size());
	externalLinks.reserve(snaxInds.size());
	for (auto iSnax : snaxInds.vec){
		for (auto iEdge : snakein.snakeconn.verts.isearch(iSnax)->edgeind){
			for (auto iVert : snakein.snakeconn.edges.isearch(iEdge)->vertind){
				if(snaxInds.find(iVert)==rsvs3d::constants::__notfound){
					externalPts.push_back(iSnax);
					externalLinks.push_back({iEdge,iVert});
				}
			}
		}
	}
	// Compute the two distances
	std::vector<int> edgeProcess;
	edgeProcess.reserve(externalPts.size());
	for (auto iSnax : snaxInds.vec){
		for (auto iEdge : snakein.snakeconn.verts.isearch(iSnax)->edgeind){
			int iVert1 = snakein.snakeconn.edges.isearch(iEdge)->vertind[0];
			int iVert2 = snakein.snakeconn.edges.isearch(iEdge)->vertind[1];
			if(externalPts.find(iVert1)!=rsvs3d::constants::__notfound
				&& externalPts.find(iVert2)!=rsvs3d::constants::__notfound){
				edgeProcess.push_back(iEdge);
			}
		}
	}

	std::vector<double>  newD;
	newD.assign(externalPts.size(), 0.0);

	for(auto iEdge : edgeProcess){
		// pass the edge index, the two snaxels, the two vertices
		int snaxA = snakein.snakeconn.edges.isearch(iEdge)->vertind[0];
		int snaxB = snakein.snakeconn.edges.isearch(iEdge)->vertind[1];
		int vertA = externalLinks[externalPts.find(snaxA)][1];
		int vertB = externalLinks[externalPts.find(snaxB)][1];
		auto tempD = EquilibrateEdge(snakein, snaxA, snaxB,
			vertA, vertB);
		newD[externalPts.find(snaxA)] +=tempD.first;
		newD[externalPts.find(snaxB)] += tempD.second;
	}
	// How to set all the values not set by the border algorithm?



	// Set the values
	int count = externalPts.size();
	for (int i = 0; i < count; ++i)
	{
		if(newD[i]<__DBL_EPSILON__){
			DisplayVector(newD);
			DisplayVector(externalPts.vec);
			DisplayVector(edgeProcess);

			RSVS3D_ERROR_NOTHROW("Negative d was calculated");
		}
		snakein.snaxs[snakein.snaxs.find(externalPts[i], true)].d = newD[i]/2.0;
	}
}

void snake::TakeSpawnStep(int minIndex, double stepLength){

	int nSnax = this->snaxs.size();
	std::vector<int> spawnVerts;
	spawnVerts.reserve(nSnax);
	for (int i = 0; i < nSnax; ++i)
	{
		if(this->snaxs(i)->index>minIndex){
			this->snaxs.elems[i].TakeSpawnStep(*this, stepLength);
		}
	}
}

void snake::TakeSmoothSpawnStep(int minIndex, double stepLength){

	int nSnax = this->snaxs.size();
	std::vector<int> spawnVerts;
	spawnVerts.reserve(nSnax);
	for (int i = 0; i < nSnax; ++i)
	{
		if(this->snaxs(i)->index>minIndex){
			spawnVerts.push_back(this->snaxs(i)->CloseToVertex());
		}
	}
	sort(spawnVerts);
	unique(spawnVerts);
	int numSuccess= 0;
	for (auto spawnVert : spawnVerts){
		numSuccess += SmoothStep(spawnVert, *this, stepLength);
		this->snaxs.ForceArrayReady();
	}
	std::cout << std::endl << "Number of successful smoothed steps : " 
		<< numSuccess << std::endl;  
}

void snax::TakeSpawnStep(snake &snakein, double stepLength){

	if(IsAproxEqual(stepLength,0.0)){
		return;
	}
	if(stepLength<-__DBL_EPSILON__ || stepLength>=0.5){
		RSVS3D_ERROR_ARGUMENT("The spawn step must be 0.5>=x>=0");
	}
	int closeVert = this->fromvert;
	if(this->d>0.5){ // if spawned in reverse
		closeVert = this->tovert;
		stepLength = -stepLength;
	}

	double minEdgeLength = INFINITY;
	for (auto edgeInd : snakein.snakemesh->verts.isearch(closeVert)->edgeind)
	{
		double testLength = snakein.snakemesh->edges.isearch(edgeInd)
			->Length(*snakein.snakemesh);
		minEdgeLength = testLength < minEdgeLength ?
			testLength : minEdgeLength;
	}
	stepLength = stepLength * minEdgeLength / snakein.snakemesh->edges.isearch(this->edgeind)
			->Length(*snakein.snakemesh);
	int nEdge=snakein.snaxs.countparent(this->edgeind);
	this->d += stepLength;
	if (this->d<0 || this->d>1.0){
		RSVS3D_ERROR_LOGIC("Distance should be between 0 and 1");
	}
	if(nEdge==1){
		return;
	} else if (nEdge>1){
		int nChanges=0;
		std::vector<int> snaxSubs;
		snakein.snaxs.findsiblings(this->edgeind,snaxSubs);
		for (auto snaxtest : snaxSubs)
		{
			auto tempSnax = snakein.snaxs(snaxtest);
			if (this->index!=tempSnax->index)
			{
				if(!((tempSnax->d<0.5) ^ (this->d<0.5))){
					if (fabs(tempSnax->d-0.5)>fabs(this->d))
					{
						this->d = tempSnax->d;
						nChanges++;
					}
				}
			} 
		}
		if (nChanges>0){
			std::cout << " " << nChanges << "," << this->d << ",";
		}
	} else {
		RSVS3D_ERROR_RANGE("Parent hashing returned a number <=0");
	}
}

// Snake operations

void snake::ForceCloseContainers(){

	this->snakeconn.ForceCloseContainers();
	if(!Check3D()){
		snaxsurfs.elems.clear();
		snaxsurfs.Init(snakeconn.surfs.size());
		snaxsurfs.PopulateIndices();
		snaxsurfs.HashArray();
		this->SetSnaxSurfs();
	}
}

void snake::CalculateTimeStep(vector<double> &dt, double dtDefault,
	double distDefault){

	int nSnax,ii,nEdge;
	vector<bool> isSnaxDone;
	double newMindDt=dtDefault;
	nSnax=snaxs.size();
	isSnaxDone.assign(nSnax,false);

	if(!snaxs.checkready()){
		snaxs.PrepareForUse();
	}
	if(this->snaxDistanceLimit_conserveShape){
		for (int i = 0; i < nSnax; ++i)
		{
			if(fabs(snaxs(i)->v >__DBL_EPSILON__)){
				newMindDt = distDefault/fabs(snaxs(i)->v);
				dtDefault =newMindDt<dtDefault ? newMindDt : dtDefault;
			}
		}
	}

	dt.assign(nSnax,dtDefault);
	dt.resize(nSnax);
	

	for(ii=0;ii<nSnax;++ii){
		if(!isSnaxDone[ii]){
			nEdge=snaxs.countparent(snaxs(ii)->edgeind);
			if (nEdge==1){
				isSnaxDone[ii]=true;
			} else if (nEdge>1){
				snaxs.CalculateTimeStepOnEdge(dt,isSnaxDone,snaxs(ii)->edgeind);
			} else {
				cerr << "Error: hashParent was not up to date" << endl;
				cerr << "Error in " << __PRETTY_FUNCTION__ << endl;
				RSVS3D_ERROR_RANGE("Incorrect hash table provided");
			}
		}
	}
}

void snake::Flip(){
	int temp;
	for (int ii=0;ii<snaxs.size();ii++){
		temp=snaxs[ii].fromvert;
		snaxs[ii].fromvert=snaxs(ii)->tovert;
		snaxs[ii].tovert=temp;

		snaxs[ii].d=1.0-snaxs[ii].d;
		snaxs[ii].v=-snaxs[ii].v;
	}
	isFlipped= !isFlipped;
	snaxs.ForceArrayReady();
}

double SnaxImpactDt(const snax &snax1,const snax &snax2){

	int isSameDir;
	double dD, dV;

	isSameDir=snax1.fromvert==snax2.fromvert;

	dD=((1.0*!isSameDir)+(1.0+(-2.0)*!isSameDir)*snax2.d)-snax1.d;
	dV=(1.0+(-2.0)*!isSameDir)*snax2.v-snax1.v;

	if (IsAproxEqual(dD,0.0)){
		return(0.0);
	}
	if (IsAproxEqual(dV,0.0)){
		return(-1.0);
	}

	return(-dD/dV);
}

void snake::SnaxImpactDetection(vector<int> &isImpact){
	// TODO make an approximate arrival 
	int nSnax,ii,nEdge;
	vector<bool> isSnaxDone;

	nSnax=snaxs.size();
	isSnaxDone.assign(nSnax,false);
	isImpact.clear();
	isImpact.reserve(nSnax*2);

	if(!snaxs.checkready()){
		snaxs.PrepareForUse();
	}

	for(ii=0;ii<nSnax;++ii){
		if(!isSnaxDone[ii]){
			nEdge=snaxs.countparent(snaxs(ii)->edgeind);
			if (nEdge==1){
				if(IsAproxEqual(snaxs(ii)->d,0.0) && (snaxs(ii)->v<=0.0)) {
					isImpact.push_back(snaxs(ii)->index);
					isImpact.push_back(-1);

				} else if (IsAproxEqual(snaxs(ii)->d,1.0) && (snaxs(ii)->v>=0.0)){
					isImpact.push_back(snaxs(ii)->index);
					isImpact.push_back(-2);
				}
				isSnaxDone[ii]=true;

			} else if (nEdge>1){
				snaxs.DetectImpactOnEdge(isImpact,isSnaxDone,snaxs(ii)->edgeind);
			} else {
				cerr << "Error: hashParent was not up to date" << endl;
				cerr << "Error in " << __PRETTY_FUNCTION__ << endl;
				RSVS3D_ERROR_RANGE("Incorrect hash table provided");
			}
		}
	}
}
void snake::SnaxAlmostImpactDetection(vector<int> &isImpact, double dDlim){

	int i;
	int nVerts, nSnax, vertInd;
	vector<int> vertSnaxCloseCnt, vertSnaxMinInd;
	vector<double> vertSnaxMinD;
	double d;
	bool flag;

	nVerts = snakemesh->verts.size();
	nSnax = snaxs.size();

	// isImpact.clear();
	// isImpact.reserve(nSnax*2);

	vertSnaxCloseCnt.assign(nVerts,0);
	vertSnaxMinInd.assign(nVerts,0);
	vertSnaxMinD.assign(nVerts,1.0);

	/*Find the closest snaxel to each vertex and the number of close
	snaxels*/
	for (i = 0; i < nSnax; ++i)
	{
		flag = false;
		if((snaxs(i)->d<dDlim) && (snaxs(i)->v<=0.0) 
			&& !IsAproxEqual(snaxs(i)->d,0.0)){
			flag = true;
			vertInd = snaxs(i)->fromvert;
			d = snaxs(i)->d;
		} else if (((1.0-snaxs(i)->d)<dDlim) && (snaxs(i)->v>=0.0) 
			&& !IsAproxEqual(snaxs(i)->d,1.0)){
			flag = true;
			vertInd = snaxs(i)->tovert;
			d = 1-snaxs(i)->d;
		}
		if(flag){
			vertInd = snakemesh->verts.find(vertInd);
			vertSnaxCloseCnt[vertInd]++;

			if(vertSnaxMinD[vertInd]>=d){
				vertSnaxMinD[vertInd]=d;
				vertSnaxMinInd[vertInd]=i;
			}
		}
	}
	// TODO: Need to add a check for the order of the snaxel on the edge.
	// if it is not 1 or the largest this action will cause crossovers. 
	for (i = 0; i < nVerts; ++i){
		if(vertSnaxCloseCnt[i]>=2){
			isImpact.push_back(snaxs(vertSnaxMinInd[i])->index);
			if(snaxs(vertSnaxMinInd[i])->d<dDlim){
				// Collapse on the from vertex
				snaxs[vertSnaxMinInd[i]].d=0.0;
				isImpact.push_back(-1);
			} else {
				// Collapse on the to vertex
				snaxs[vertSnaxMinInd[i]].d=1.0;
				isImpact.push_back(-2);
			}
		}
	}
	snaxs.PrepareForUse();
}

void snaxarray::DetectImpactOnEdge(vector<int> &isImpact, 
	vector<bool> &isSnaxDone, int edgeInd){
	int nSnax,ii,jj, dOrd;
	vector<int> snaxSubs;
	double impactTime;
	snaxSubs=ReturnDataEqualRange(edgeInd, hashParent);
	nSnax=snaxSubs.size();
	for(ii=0;ii<nSnax;++ii){

		isSnaxDone[snaxSubs[ii]]=true;

		for(jj=ii+1; jj<nSnax; ++jj){
			
			impactTime=SnaxImpactDt(elems[snaxSubs[ii]],elems[snaxSubs[jj]]);
			dOrd=abs(elems[snaxSubs[ii]].orderedge-elems[snaxSubs[jj]].orderedge);
			if(dOrd==1 && IsAproxEqual(impactTime,0.0)){

				isImpact.push_back((elems[snaxSubs[ii]].index));
				isImpact.push_back((elems[snaxSubs[jj]].index));


				isImpact.push_back((elems[snaxSubs[jj]].index));
				isImpact.push_back((elems[snaxSubs[ii]].index));

			} 
			
		}
		
		if(IsAproxEqual(elems[snaxSubs[ii]].d,0.0) && (elems[snaxSubs[ii]].v<=0.0) 
			&& elems[snaxSubs[ii]].orderedge==1) {
			isImpact.push_back(elems[snaxSubs[ii]].index);
			isImpact.push_back(-1);


		} else if (IsAproxEqual(elems[snaxSubs[ii]].d,1.0) && (elems[snaxSubs[ii]].v>=0.0)
			&& elems[snaxSubs[ii]].orderedge==nSnax){
			isImpact.push_back(elems[snaxSubs[ii]].index);
			isImpact.push_back(-2);

		}
	}
}

void snake::OrientFaces(){
	/*
	Orients either surfaces or edges depending. 
	*/
	if (this->Check3D()){
		this->OrientSurfaceVolume();
	} else {
		this->OrientEdgeSurface();
	}
}
void snake::OrientEdgeSurface(){
	cerr << "Warning: not coded yet in " << __PRETTY_FUNCTION__ << endl;
}
void snake::OrientSurfaceVolume(){ 
	// Orders the surf.voluind [c0 c1] such that the surface normal vector points
	// from c0 to c1
	// This is done by using the surface normals and checking they go towards
	// the centre of the cell

	int nBlocks,ii,jj, ni,nj,kk,ll;
	vector<int> surfOrient;
	vector<bool> isFlip;
	double dotProd;
	coordvec centreVolu, normalVec;


	nBlocks=snakeconn.OrientRelativeSurfaceVolume(surfOrient);
	isFlip.assign(nBlocks,false);
	//========================================
	//  Select direction using coordinate geometry
	//     use a surface 

	for (ii=1; ii<= nBlocks; ii++){
		jj=-1; nj=surfOrient.size();
		do{
			jj++;
			while(jj<nj && ii!=abs(surfOrient[jj]))
				{jj++;}
			if(jj==nj){ // if the orientation cannot be defined
				dotProd=1.0;
				kk=0;
				// cerr << "Warning: Cell orientations could not be computed " << endl;
				// cerr << "			in " << __PRETTY_FUNCTION__ << endl;
				break;
			}
			kk=snakeconn.surfs(jj)->voluind[0]==0;
			ll=snaxs.isearch(snakeconn.edges.isearch(snakeconn.surfs(jj)->edgeind[0])->vertind[0])->tovert;
			centreVolu=snakemesh->verts.isearch(ll)->coord;
			ll=snaxs.isearch(snakeconn.edges.isearch(snakeconn.surfs(jj)->edgeind[0])->vertind[0])->fromvert;
			centreVolu.substract(snakemesh->verts.isearch(ll)->coord);
			
			normalVec=snakeconn.CalcPseudoNormalSurf(snakeconn.surfs(jj)->index);

			dotProd=centreVolu.dot(normalVec.usedata());
		} while (!isfinite(dotProd) || (fabs(dotProd)<numeric_limits<double>::epsilon()));


		isFlip[ii-1]= (((dotProd<0.0) && (kk==0)) || ((dotProd>0.0) && (kk==1)));
	}
	if(nBlocks>0){
		ni=surfOrient.size();
		for(ii=0; ii< ni; ++ii){
			if(isFlip[abs(surfOrient[ii])-1]){
				snakeconn.surfs.elems[ii].FlipVolus();
			}
		}
	}
}
void snake::AssignInternalVerts() {
	/*
	Assigns internal vertices to the snake based on the snakemesh connectivity.
	*/
	int ii, ni;
	vector<int> vertBlock;

	FindBlockSnakeMeshVerts(vertBlock);
	this->isFlipped=false;
	ni=vertBlock.size();
	for(ii=0; ii < ni; ++ii){
		isMeshVertIn[ii] = vertBlock[ii]>0;
	}
}
int snake::FindBlockSnakeMeshVerts(vector<int> &vertBlock) const{
	// int mesh::ConnectedVertex(vector<int> &vertBlock) const{
	// Fills a vector with a number for each vertex corresponding to a
    // group of connected edges it is part of , can be used close surfaces in 2D or volumes
    // in 3D.
    // Uses a flood fill with queue method


	int nVertExplored,nVerts,nBlocks,nCurr,nEdgesCurr,ii,jj,kk,nSnaxs,ll,nl;
   vector<bool> vertStatus, snaxStatus; // 1 explored 0 not explored

   vector<int> currQueue, nextQueue, tempSnax; // Current and next queues of indices

   // Preparation of the arrays;
   nVerts=snakemesh->verts.size();
   nSnaxs=snaxs.size();
   nBlocks=0;
   nVertExplored=0;

   snaxStatus.assign(nSnaxs,false);
   vertStatus.assign(nVerts,false);
   vertBlock.assign(nVerts,0);
   currQueue.reserve(nVerts/2);
   nextQueue.reserve(nVerts/2);

   
   // While Loop, while not all vertices explored
	while(nVertExplored<nSnaxs){

		  // if currQueue is empty start new block
		if(currQueue.size()<1){

			 //cout << "Block " << nBlocks << " - " << nVertExplored << " - " << nVerts << endl;
			ii=0;
			while(ii<nSnaxs && snaxStatus[ii]){
				ii++;
			}
			if (snaxStatus[ii]){
				cerr << "Error starting point for loop not found despite max number of vertex not reached" <<endl;
				cerr << "Error in " << __PRETTY_FUNCTION__ << endl;
				RSVS3D_ERROR_RANGE(" : Starting point for block not found");
			}
			currQueue.push_back(snakemesh->verts.find(snaxs(ii)->fromvert));
			nBlocks++;

		}
		  // Explore current queue
		nCurr=currQueue.size();
		// cout << "nCurr " << nCurr << " nVertExplored " << nVertExplored << 
		// 	" nSnax " << nSnaxs << " " << nBlocks << endl;
		for (ii = 0; ii < nCurr; ++ii){
			if (!vertStatus[currQueue[ii]]){
				vertBlock[currQueue[ii]]=nBlocks;
				nEdgesCurr=snakemesh->verts(currQueue[ii])->edgeind.size();
				for(jj=0;jj<nEdgesCurr;++jj){
					// only follow edge if no snaxel exists on it
					if(snaxs.countparent(snakemesh->verts(currQueue[ii])->edgeind[jj])<1){
						kk=int(snakemesh->edges.isearch(
							snakemesh->verts(currQueue[ii])->edgeind[jj])->vertind[0]
							==snakemesh->verts(currQueue[ii])->index);
						nextQueue.push_back(snakemesh->verts.find(
							snakemesh->edges.isearch(
								snakemesh->verts(currQueue[ii])->edgeind[jj])->vertind[kk]));
						#ifdef SAFE_ALGO
						if (snakemesh->verts.find(
							snakemesh->edges.isearch(snakemesh->verts(
								currQueue[ii])->edgeind[jj])->vertind[kk])==-1){
							cerr << "Edge index: " << snakemesh->verts(
								currQueue[ii])->edgeind[jj] << " vertex index:" <<  
							snakemesh->edges.isearch(snakemesh->verts(
								currQueue[ii])->edgeind[jj])->vertind[kk] << endl;
							cerr << "Edge connected to non existant vertex" <<endl;
							cerr << "Error in " << __PRETTY_FUNCTION__ << endl;
							RSVS3D_ERROR_RANGE(" : Vertex not found");
						}
						#endif
					} else {
						snaxs.findsiblings(snakemesh->verts(currQueue[ii])->edgeind[jj],
							tempSnax);
						nl=tempSnax.size();
						for(ll=0; ll<nl; ++ll){
							if(snakemesh->verts(currQueue[ii])->index
								==snaxs(tempSnax[ll])->fromvert){
								nVertExplored+=int(!snaxStatus[tempSnax[ll]]);
								snaxStatus[tempSnax[ll]]=true;
							}
						}
					}
				}
				vertStatus[currQueue[ii]]=true;
			}
	    }

		// Reset current queue and set to next queue
		currQueue.clear();
		currQueue.swap(nextQueue);

	}
	return(nBlocks);
}
grid::limits snake::Scale(const grid::limits &newSize){

	auto boundBox = this->snakemesh->BoundingBox();
	// Scales the snakemesh to the correct size
	auto transformation = this->snakemesh->Scale(newSize);
	// Uses the same transformation on the snakeconn
	this->snakeconn.LinearTransform(transformation);
	// Applies the transformation to the family of the snakemesh
	// ie: the VOS mesh. 
	this->snakemesh->LinearTransformFamily(transformation);

	return boundBox;
}




// -------------------------------------------------------------------------------------------
// TEST CODE
// -------------------------------------------------------------------------------------------

//in snakstruct_test.cpp