#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <unordered_map>

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
	"; fromvert " << fromvert  << "; tovert " << tovert  << "; isfreeze " << isfreeze << "; orderedge " << orderedge << endl;
}

void snaxedge::disp() const{
	cout << "snaxedge #" << index << " ";
	normvector.disp();
}

void snaxsurf::disp() const{
	cout << "snaxsurf #" << index << " ";
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

// set constructors (used to avoid a variable being unknowingly forgotten)
inline void snax::set(int indexin, double din,double vin,int fromvertin,int tovertin,
	int edgeindin,int isfreezein,int orderedgein)
{
	index=indexin;
	d=din;
	v=vin;
	 fromvert=fromvertin; 	// root vertex of *snakemesh
	 tovert=tovertin; 	// destination vertex of *snakemesh
	 edgeind=edgeindin; 	// edge of snakemesh
	 isfreeze=isfreezein; 	// freeze status 
	 orderedge=orderedgein;
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

	void snake::PrepareForUse() {

		snaxs.PrepareForUse();
		snaxedges.PrepareForUse();
		snaxsurfs.PrepareForUse();
		snakeconn.PrepareForUse();
		snakemesh->PrepareForUse();

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
	void snake::SetMaxIndex(){
	// prefer use of PrepareForUse in final code as it includes checks to 
	// avoid unnecessary repetiion of work

		snaxs.SetMaxIndex();
		snaxedges.SetMaxIndex();
		snaxsurfs.SetMaxIndex();
		snakeconn.SetMaxIndex();
		snakemesh->SetMaxIndex();

	}

	void snake::Init(mesh *snakemeshin,int nSnax, int nEdge, int nSurf, int nVolu){
// Initialise a 
		is3D=snakemeshin->volus.size()>0;

		snaxs.Init(nSnax);
		snaxedges.Init(nEdge);
		snaxsurfs.Init(nSurf);

		snakeconn.Init(nSnax, nEdge, nSurf, nVolu);

		snakemesh=snakemeshin;

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

	int nVert,nEdge,nVolu,nSurf;

   // Define Max indices in current mesh
	this->GetMaxIndex(&nVert,&nEdge,&nSurf,&nVolu);
	other.ChangeIndices(nVert,nEdge,nSurf,nVolu);
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
		snaxs[ii].d=snaxs[ii].d<0.0 ? 0.0 : snaxs[ii].d;
		snaxs[ii].d=snaxs[ii].d>1.0 ? 1.0 : snaxs[ii].d;
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

// Snake Spawning at Vertex.
void SpawnAtVertex(snake& snakein,int indVert){

	snake newsnake;
	int subVert,nVert, nEdge,nSurf,nVolu;

	vector<int> vertInds,edgeInds,surfInds,voluInds,edgeInds2;
	vector<int> vertSubs,edgeSubs,surfSubs,voluSubs;
	vector<int> vertSubsTemp,edgeSubsTemp,surfSubsTemp,voluSubsTemp;
	vector<int>::iterator itVecInt; 
	unordered_multimap<int,int> hashEdgeInds,hashVoluInds,hashSurfInds,hashVertInds;


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
	voluSubs=snakein.snakemesh->volus.find_list(voluInds);

	//OperArrayStructMethod(snakein.snakemesh->surfs, surfSubs, &surf::isready, 
	//	ii, std::logical_and<bool>());
	nVert=edgeInds.size();
	nEdge=surfInds.size();
	nSurf=voluInds.size();
	nVolu=1;

	newsnake.Init(snakein.snakemesh,nVert,nEdge,nSurf,nVolu);
	// Generates snaxels and vertices
	SpawnAtVertexVert(newsnake,nVert, indVert,subVert, surfInds,edgeInds,
		edgeSubs,hashSurfInds);
	// Generate snake edges
	SpawnAtVertexEdge(newsnake, nEdge ,surfInds, edgeInds,	voluInds,
		surfSubs, hashEdgeInds,hashVoluInds);
	// Generate Snake surfaces
	SpawnAtVertexSurf(newsnake,nSurf,surfInds ,voluInds,voluSubs,hashSurfInds);
	// Generate Volume
	SpawnAtVertexVolu(newsnake,nSurf);

	

	snakein.PrepareForUse();
	newsnake.PrepareForUse();

	snakein.MakeCompatible_inplace(newsnake);
	newsnake.PrepareForUse();
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
			newsnake.snakemesh->edges(edgeSubs[ii])->vertind[jj],edgeInds[ii],0,1);

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
}
void SpawnAtVertexSurf(snake& newsnake,int nSurf,const vector<int> &surfInds, const vector<int> &voluInds,
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
		
		for(jj=0;jj<int(surfSubsTemp.size());++jj){
			if (surfSubsTemp[jj]>=0){
				newsnake.snakeconn.surfs[ii].edgeind.push_back(surfSubsTemp[jj]);
				
			}
		}
	}
	newsnake.snakeconn.surfs.ChangeIndices(0,1,0,0);
}

void SpawnAtVertexVolu(snake& newsnake, int nSurf){
	int ii;

	newsnake.snakeconn.volus.PopulateIndices();
	newsnake.snakeconn.volus[0].surfind.reserve(nSurf);
	for (ii=0;ii<nSurf;++ii){
		newsnake.snakeconn.volus[0].surfind.push_back(ii+1);
	}
}
// -------------------------------------------------------------------------------------------
// TEST CODE
// -------------------------------------------------------------------------------------------

//in snakstruct_test.cpp