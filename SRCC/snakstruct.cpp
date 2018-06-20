#include <iostream>
#include <stdexcept>
#include <vector>
#include <cmath>

#include "mesh.hpp"
#include "snakstruct.hpp"
#include "arraystructures.hpp"

using namespace std;

// Implementation of features


// Class Method implementaton 

// -------------------------------------------------------------------------------------------
// Implementatation of coordvec 
// -------------------------------------------------------------------------------------------
double& coordvec::operator[](int a)
{
	// returns a reference to the value and can be used to set a value
	#ifdef SAFE_ACCESS // adds a check in debug mode
	if ((unsigned_int(a)>=3) | (0>a)){
		throw range_error ("Index is out of range");
	}
	#endif //SAFE_ACCESS
	isuptodate=0;
	return(elems[a]);
}
double coordvec::operator()(int a) const
{
	// returns the value (cannot be used to set data)
	#ifdef SAFE_ACCESS // adds a check in debug mode
	if ((unsigned_int(a)>=3) | (0>a)){
		throw range_error ("Index is out of range");
	}
	#endif // SAFE_ACCESS
	return(elems[a]);
}

coordvec coordvec::Unit() const 
{
	coordvec unitCoordVec;
	for (int ii=0;ii<3;++ii){
		unitCoordVec.elems[ii]=elems[ii]/norm;
	}
	unitCoordVec.norm=1;
	unitCoordVec.isuptodate=1;
	return (unitCoordVec);
}

double coordvec::Unit(const int a) const
{
	#ifdef SAFE_ACCESS // adds a check in debug mode
	if ((unsigned_int(a)>=3) | (0>a)){
		throw range_error ("Index is out of range");
	}
	#endif // SAFE_ACCESS
	if (isuptodate==0){
		cerr << "Warning: NORM of coordvec potentially obsolete " << endl;
		cerr << "          in coordvec::Unit(const int a) const" << endl; 
		cerr << "          To avoid this message perform read operations on coordvec using the () operator" << endl; 
	}
	return(elems[a]/norm);
}

double coordvec::GetNorm(){
	// TEST_RANGE
	if (~isuptodate){
		this->CalcNorm();
	}
	return(norm);
}

double coordvec::CalcNorm(){
	norm=sqrt(pow(elems[0],2)+pow(elems[1],2)+pow(elems[2],2));
	isuptodate=1;
	return(norm);
}

void coordvec::PrepareForUse(){
	this->CalcNorm();
}

void coordvec::assign(double a,double b,double c){
	elems[0]=a;
	elems[1]=b;
	elems[2]=c;
	this->CalcNorm();
}	

void coordvec::disp() const{
	cout << "coordvec [" << elems[0] << ","<< elems[1]<< ","<< elems[2] << "] norm "
	<< norm << " utd " << isuptodate <<endl;
}
// -------------------------------------------------------------------------------------------
// Implementatation of snax - snaxedge - snaxsurf 
// -------------------------------------------------------------------------------------------
void snax::disp() const{
	cout << "snax #" << index << "; d " << d  << "; v " << v  << "; edgeind " << edgeind <<
	"; fromvert " << fromvert  << "; tovert " << tovert  << "; isfreeze " << isfreeze << "; orderedge " << orderedge << " | isBorder " << isBorder << endl;
}

void snaxedge::disp() const{
	cout << "snaxedge #" << index << " | surfind " << surfind << " " << " | isBorder " << isBorder << " " ;
	normvector.disp();
}

void snaxsurf::disp() const{
	cout << "snaxsurf #" << index << index << " | voluind " << voluind << " " << " | isBorder " << isBorder << " ";
	normvector.disp();
}

// file i/o
void snax::write(FILE *fid) const{

	fprintf(fid, "%i %lf %lf %i %i %i %i %i \n",index,d,v,edgeind, fromvert, tovert, 
		isfreeze, orderedge);
}

void snaxedge::write(FILE *fid) const{
	fprintf(fid, "%i %i %lf %lf %lf \n",index,surfind,normvector(0),normvector(1),normvector(2));
}

void snaxsurf::write(FILE *fid) const{
	fprintf(fid, "%i %i %lf %lf %lf \n",index,voluind,normvector(0),normvector(1),normvector(2));
}
void snax::read(FILE *fid) {

	fscanf(fid,"%i %lf %lf %i %i %i %i %i ", &index, &d, &v, &edgeind, &fromvert, &tovert,
		&isfreeze, &orderedge);
}

void snaxedge::read(FILE *fid) {
	fscanf(fid,"%i %i %lf %lf %lf ",&index,&surfind,&normvector[0],&normvector[1],&normvector[2]);
}

void snaxsurf::read(FILE *fid) {
	fscanf(fid,"%i %i %lf %lf %lf ",&index,&voluind,&normvector[0],&normvector[1],&normvector[2]);
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
}
void snaxedge::disptree(const mesh &meshin, int n) const{
   
   disp();
}
void snaxsurf::disptree(const mesh &meshin, int n) const{
   
   disp();
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

// -------------------------------------------------------------------------------------------
// Extensions of SnakStruct for snaxarray
// -------------------------------------------------------------------------------------------


bool snaxarray::checkready() 
{	
	readyforuse=SnakStruct<snax>::checkready();
	
	readyforuse=(readyforuse && isOrderedOnEdge);

	return(readyforuse);
	
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
	snakemesh->PrepareForUse();
	is3D=int(snakemesh->volus.size())>0;

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

void snake::Concatenate(const snake &other){

	if(this->snakemesh==other.snakemesh){
		this->snaxs.Concatenate(other.snaxs);
		this->snaxedges.Concatenate(other.snaxedges);
		this->snaxsurfs.Concatenate(other.snaxsurfs);
		this->snakeconn.Concatenate(other.snakeconn);
	} else {

		throw invalid_argument("Concatenation of snakes with different base meshes not supported");
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



// Snake Movement
void snake::UpdateDistance(double dt){
	int ii;

	for (ii = 0; ii < int(snaxs.size()); ++ii)
	{
		snaxs[ii].d=snaxs[ii].d+snaxs[ii].v*dt;
		snaxs[ii].d=snaxs[ii].d<0.0 ? 0.0 : (snaxs[ii].d>1.0 ? 1.0 : snaxs[ii].d);
		
	}
}
void snake::UpdateDistance(const vector<double> &dt){
	int ii;
	if (int(dt.size())==snaxs.size()){
		for (ii = 0; ii < int(snaxs.size()); ++ii)
		{
			snaxs[ii].d=snaxs[ii].d+snaxs[ii].v*dt[ii];
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
			snakeconn.verts[ii].coord[jj]=(snakemesh->verts(toVertSub)->coord[jj]
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

	vector<int> edgeInds,snaxSubs,orderSubs;
	vector<double> orderSnax,unorderSnax;
	unordered_multimap<double,int> hashOrder;
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

void snaxarray::CalculateTimeStepOnEdge(vector<double> &dt, vector<bool> &isSnaxDone, int edgeInd){
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
				dt[snaxSubs[ii]]= (dt[snaxSubs[ii]]>impactTime) ? impactTime : dt[snaxSubs[ii]];
				dt[snaxSubs[jj]]= (dt[snaxSubs[jj]]>impactTime) ? impactTime : dt[snaxSubs[jj]];
			}
		}
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

void snake::CalculateTimeStep(vector<double> &dt, double dtDefault){

	int nSnax,ii,nEdge;
	vector<bool> isSnaxDone;

	nSnax=snaxs.size();
	isSnaxDone.assign(nSnax,false);
	dt.assign(nSnax,dtDefault);
	dt.resize(nSnax);

	if(!snaxs.checkready()){
		snaxs.PrepareForUse();
	}

	for(ii=0;ii<nSnax;++ii){
		if(!isSnaxDone[ii]){
			nEdge=snaxs.countparent(snaxs(ii)->edgeind);
			if (nEdge==1){
				//dt[ii]=dtDefault;
				isSnaxDone[ii]=true;
			} else if (nEdge>1){
				snaxs.CalculateTimeStepOnEdge(dt,isSnaxDone,snaxs(ii)->edgeind);
			} else {
				cerr << "Error: hashParent was not up to date" << endl;
				cerr << "Error in " << __PRETTY_FUNCTION__ << endl;
				throw range_error ("Incorrect hash table provided");
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
				if(IsAproxEqual(snaxs(ii)->d,0.0) && (snaxs(ii)->v<0.0)) {
					isImpact.push_back(snaxs(ii)->index);
					isImpact.push_back(-1);

				} else if (IsAproxEqual(snaxs(ii)->d,1.0) && (snaxs(ii)->v>0.0)){
					isImpact.push_back(snaxs(ii)->index);
					isImpact.push_back(-2);
				}
				isSnaxDone[ii]=true;

			} else if (nEdge>1){
				snaxs.DetectImpactOnEdge(isImpact,isSnaxDone,snaxs(ii)->edgeind);
			} else {
				cerr << "Error: hashParent was not up to date" << endl;
				cerr << "Error in " << __PRETTY_FUNCTION__ << endl;
				throw range_error ("Incorrect hash table provided");
			}
		}
	}

}

void snaxarray::DetectImpactOnEdge(vector<int> &isImpact, vector<bool> &isSnaxDone, int edgeInd){
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
		
		if(IsAproxEqual(elems[snaxSubs[ii]].d,0.0) && (elems[snaxSubs[ii]].v<0.0) 
			&& elems[snaxSubs[ii]].orderedge==1) {
			isImpact.push_back(elems[snaxSubs[ii]].index);
		isImpact.push_back(-1);


	} else if (IsAproxEqual(elems[snaxSubs[ii]].d,1.0) && (elems[snaxSubs[ii]].v>0.0)
		&& elems[snaxSubs[ii]].orderedge==nSnax){
		isImpact.push_back(elems[snaxSubs[ii]].index);
		isImpact.push_back(-2);

	}
}

}




// -------------------------------------------------------------------------------------------
// TEST CODE
// -------------------------------------------------------------------------------------------

//in snakstruct_test.cpp