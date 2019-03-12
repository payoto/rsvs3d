
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <fstream>
#include <functional>
#include <cmath>

#include "mesh.hpp"
#include "warning.hpp"


// Class function definitions
// Mesh Linking methods
using namespace std;


// -------------------------------------------------------------------------------------------
// Implementatation of coordvec 
// -------------------------------------------------------------------------------------------
double& coordvec::operator[](int a)
{
	// returns a reference to the value and can be used to set a value
	#ifdef SAFE_ACCESS // adds a check in debug mode
	if ((unsigned_int(a)>=3) | (0>a)){
		RSVS3D_ERROR_RANGE("Index is out of range");
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
		RSVS3D_ERROR_RANGE("Index is out of range");
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
		RSVS3D_ERROR_RANGE("Index is out of range");
	}
	#endif // SAFE_ACCESS
	if (isuptodate==0){
		cerr << "Warning: NORM of coordvec potentially obsolete " << endl;
		cerr << "		  in coordvec::Unit(const int a) const" << endl; 
		cerr << "		  To avoid this message perform read operations on coordvec using the () operator" << endl; 
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
double coordvec::GetNorm() const {
	// TEST_RANGE
	if (~isuptodate){
		RSVS3D_ERROR_ARGUMENT("coordvec is not ready for norm return");
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

// Math implementations

void coordvec::max(const vector<double> &vecin){
	int n=int(elems.size()<=vecin.size()?elems.size():vecin.size());
	for (int ii = 0; ii < n; ++ii)
	{
		elems[ii]=elems[ii]>=vecin[ii] ? elems[ii] : vecin[ii];
	}
}
void coordvec::min(const vector<double> &vecin){
	int n=int(elems.size()<=vecin.size()?elems.size():vecin.size());
	for (int ii = 0; ii < n; ++ii)
	{
		elems[ii]=elems[ii]<=vecin[ii] ? elems[ii] : vecin[ii];
	}
}

void coordvec::add(const vector<double> &vecin){
	int n=int(elems.size()<=vecin.size()?elems.size():vecin.size());
	for (int ii = 0; ii < n; ++ii)
	{
		elems[ii]=elems[ii]+vecin[ii];
	}
}
void coordvec::substract(const vector<double> &vecin){
	int n=int(elems.size()<=vecin.size()?elems.size():vecin.size());
	for (int ii = 0; ii < n; ++ii)
	{
		elems[ii]=elems[ii]-vecin[ii];
	}
}
void coordvec::substractfrom(const vector<double> &vecin){
	int n=int(elems.size()<=vecin.size()?elems.size():vecin.size());
	for (int ii = 0; ii < n; ++ii)
	{
		elems[ii]=vecin[ii]-elems[ii];
	}
}
void coordvec::div(const vector<double> &vecin){
	int n=int(elems.size()<=vecin.size()?elems.size():vecin.size());
	for (int ii = 0; ii < n; ++ii)
	{
		elems[ii]=elems[ii]/vecin[ii];
	}
}
void coordvec::mult(const vector<double> &vecin){
	int n=int(elems.size()<=vecin.size()?elems.size():vecin.size());
	for (int ii = 0; ii < n; ++ii)
	{
		elems[ii]=elems[ii]*vecin[ii];
	}
} 

void coordvec::div(double scalin){
	int n=int(elems.size());
	for (int ii = 0; ii < n; ++ii)
	{
		elems[ii]=elems[ii]/scalin;
	}
}
void coordvec::mult(double scalin){
	int n=int(elems.size());
	for (int ii = 0; ii < n; ++ii)
	{
		elems[ii]=elems[ii]*scalin;
	}
}

vector<double> coordvec::cross(const std::vector<double> &vecin) const {

	vector<double> retVec;
	retVec.assign(3,0.0);

	for(int ii=3; ii<6; ii++){
		retVec[ii%3]=elems[(ii+1)%3]*vecin[(ii+2)%3]
		-elems[(ii-1)%3]*vecin[(ii-2)%3];
	}
	return(retVec);
}

double coordvec::dot(const std::vector<double> &vecin) const {

	double retVec=0.0;
	
	for(int ii=0; ii<3; ii++){
		retVec+=elems[ii]*vecin[ii];
	}
	return(retVec);
}
double coordvec::angle(const coordvec &coordin) const {

	return acos(
		this->dot(coordin.usedata())
		/(coordin.GetNorm()*this->GetNorm())
		);
}

//// ----------------------------------------
// Implementation of point location
//// ----------------------------------------
// Detects which volume (if any a point is in)

/**
 * @brief      Calculates the distance from a vertex to a plane.
 * calculates the distance from a plane to a vertex, with the plane
 * defined by three vertices. 
 * 
 * Two optional arguments can be provided to avoid the need for
 * memory allocation if this is called in a loop. For max speedup 
 * if testing a surface multiple times against many vertices temp2
 * can be reused
 *
 * @param[in]  planeVert1  The plane vertex 1
 * @param[in]  planeVert2  The plane vertex 2
 * @param[in]  planeVert3  The plane vertex 3
 * @param[in]  testVertex  The test vertex
 * @param[in|out]  temp1       The temporary array 1
 * @param[in|out]  temp2       The temporary array 2
 *
 * @return     the distance from the plane to the vertex
 */
double VertexDistanceToPlane(const vector<double> &planeVert1, 
	const vector<double> &planeVert2,
	const vector<double> &planeVert3,
	const vector<double> &testVertex,
	coordvec &temp1, coordvec &temp2){
	
	double planeDistance=0.0;
	if(testVertex.size()!=3){
		RSVS3D_ERROR_ARGUMENT("testVertices.size() must be of size 3.");
	}
	temp1 = planeVert1;
	temp2 = planeVert1;
	temp1.substractfrom(planeVert2);
	temp2.substractfrom(planeVert3);

	temp2 = temp1.cross(temp2.usedata()); // causes allocation

	planeDistance = temp2.dot(testVertex)-temp2.dot(planeVert1);

	return planeDistance;
}
/**
 * Calculates the distance from a set of vertices to a plane.
 * 
 * Calculates the distance from a plane to a vertex, with the plane
 * defined by three vertices. 
 * 
 * Two optional arguments can be provided to avoid the need for
 * memory allocation if this is called in a loop. For max speedup 
 * if testing a surface multiple times against many vertices temp2
 * is reused internally allowing surface properties to only be 
 * calculated once.
 *
 * @param[in]  planeVert1  The plane vertex 1
 * @param[in]  planeVert2  The plane vertex 2
 * @param[in]  planeVert3  The plane vertex 3
 * @param[in]  testVertices  The test vertices
 * @param[in|out]  temp1       The temporary array 1
 * @param[in|out]  temp2       The temporary array 2
 *
 * @return     the distance from the plane to the vertex
 */
vector<double> VerticesDistanceToPlane(const vector<double> &planeVert1, 
	const vector<double> &planeVert2,
	const vector<double> &planeVert3,
	const vector<double> &testVertices,
	coordvec &temp1, coordvec &temp2){
	
	if(testVertices.size()%3 != 0){
		RSVS3D_ERROR_ARGUMENT("testVertices.size() must be a multiple of 3.");
	}

	std::vector<double> planeDistances, testVertex;
	size_t nPts = testVertices.size()/3;
	planeDistances.assign(nPts, 0.0);
	if(nPts==0){
		return planeDistances;
	}
	testVertex.assign(3, 0.0);
	testVertex.assign(testVertices.begin(),	testVertices.begin()+3);
	planeDistances[0]=VertexDistanceToPlane(planeVert1,planeVert2,planeVert3,
		testVertex, temp1, temp2);
	auto planePosition = temp2.dot(planeVert1);
	for (size_t i = 1; i < nPts; ++i)
	{
		testVertex.assign(testVertices.begin()+3*i, testVertices.begin()+3*i);
		planeDistances[i] = temp2.dot(testVertex)-planePosition;
	}

	return planeDistances;
}

double VertexDistanceToPlane(const vector<double> &planeVert1, 
	const vector<double> &planeVert2,
	const vector<double> &planeVert3,
	const vector<double> &testVertex){
	coordvec temp1, temp2;
	return VertexDistanceToPlane(planeVert1, planeVert2, 
		planeVert3, testVertex, temp1, temp2);
}
vector<double> VerticesDistanceToPlane(const vector<double> &planeVert1, 
	const vector<double> &planeVert2,
	const vector<double> &planeVert3,
	const vector<double> &testVertices){
	coordvec temp1, temp2;
	return VerticesDistanceToPlane(planeVert1, planeVert2, 
		planeVert3, testVertices, temp1, temp2);
}

/**
 * @brief      Finds for each vertex, the volume object containing it.
 * 
 * This only works robustly for outside points for convex meshes.
 *
 * @param[in]  testVertices  The test vertices
 * @param[in]  sizeVert      The size of each vertex data
 *
 * @return     returns a list of indices containing the same number of
 * 	values as there are input vertices (testVertices/sizeVert) 
 */
vector<int> mesh::VertexInVolume(const vector<double> testVertices, int sizeVert) const {
	if(sizeVert<3){
		RSVS3D_ERROR_ARGUMENT("sizeVert must be at least 3.");
	}
	if(testVertices.size()%sizeVert != 0){
		RSVS3D_ERROR_ARGUMENT("testVertices.size() must be a multiple of sizeVert.");
	}

	vector<int> vertVolumeIndices;
	size_t nPts = testVertices.size()/sizeVert;
	vector<double> tempCoord;
	vertVolumeIndices.assign(nPts,-1);
	#ifdef RSVS_DIAGNOSTIC_RESOLVED
	cout << endl << " Number of volumes explored before final "
		"candidate out of " << this->volus.size() << endl ;
	#endif
	// Check orientation (OrientFaces is relative which means it can be facing in or out)
	// simply test if the 1st one is 0:
	tempCoord.assign(testVertices.begin(),testVertices.begin()+3);
	bool needFlip = meshhelp::VertexInVolume(*this, tempCoord)==0;
	for (size_t i = 0; i < nPts; ++i)
	{
		tempCoord.assign(testVertices.begin()+i*sizeVert,
			testVertices.begin()+(i*sizeVert+3));
		vertVolumeIndices[i] = meshhelp::VertexInVolume(*this, tempCoord, needFlip);
	}
	#ifdef RSVS_DIAGNOSTIC_RESOLVED 
	cout << endl ;
	#endif
	return vertVolumeIndices;
}
//// ----------------------------------------
// Implementation of mesh dependence
//// ----------------------------------------

int meshdependence::AddParent(mesh* meshin){
	HashedVectorSafe<int,int> temp;

	parentmesh.push_back(meshin);
	parentconn.push_back(temp);
	return(parentmesh.size());
}
int meshdependence::AddChild(mesh* meshin){
	HashedVectorSafe<int,int> temp;

	childmesh.push_back(meshin); 
	return(childmesh.size());
}
void meshdependence::RemoveChild(mesh *meshin){

	for (int i = 0; i < int(childmesh.size()); ++i)
	{
		if(meshin==childmesh[i]){
			childmesh.erase(childmesh.begin()+i);
		} 
	}
}
void meshdependence::RemoveParent(mesh *meshin){

	for (int i = 0; i < int(parentmesh.size()); ++i)
	{
		if(meshin==parentmesh[i]){
			parentmesh.erase(parentmesh.begin()+i);
			parentconn.erase(parentconn.begin()+i);
			break;
		}
	}
	nParents=parentmesh.size();
}

void meshdependence::AddParent(mesh* meshin, vector<int> &parentind){
	/*
	parentind needs to be a list of the parent.volus.index matched to elemind
	that means that running temp.find(parent.volus.index) will return the 
	subscribts of all the child.volus that are contained in that volu

	Alternatively running temp(volu.find(volus[a].index)) will return the 
	corresponding parent.volus index.
	*/
	HashedVectorSafe<int,int> temp;

	if (parentind.size()!=elemind.size()){
		RSVS3D_ERROR_ARGUMENT("parent and child index list must be the same size.");
	}

	temp=parentind;
	temp.GenerateHash(); 

	parentmesh.push_back(meshin);
	parentconn.push_back(temp);
	nParents=parentmesh.size();
}

void mesh::RemoveFromFamily(){
	int jj;

	for (jj = 0; jj<int(meshtree.parentmesh.size()); jj++){ 
		meshtree.parentmesh[jj]->meshtree.RemoveChild(this);
	}
	for (jj = 0; jj<int(meshtree.childmesh.size()); jj++){ 
		meshtree.childmesh[jj]->meshtree.RemoveParent(this);
	}
}

void mesh::AddChild(mesh *meshin){
	meshtree.AddChild(meshin);
	meshin->meshtree.AddParent(this);
}
void mesh::AddParent(mesh *meshin){
	meshtree.AddParent(meshin);
	meshin->meshtree.AddChild(this);
}

void mesh::AddChild(mesh *meshin, vector<int> &parentind){
	if(!meshDepIsSet){
		SetMeshDepElm();
	}
	meshtree.AddChild(meshin);
	meshin->meshtree.AddParent(this,parentind);
}
void mesh::AddParent(mesh *meshin, vector<int> &parentind){
	if(!meshDepIsSet){
		SetMeshDepElm();
	}
	meshtree.AddParent(meshin,parentind);
	meshin->meshtree.AddChild(this);
}

void mesh::SetMeshDepElm(){
	// 
	int ii;
	meshtree.elemind.clear();
	switch (meshDim){
		case 0:
		meshtree.elemind.reserve(verts.size());
		for (ii=0; ii<int(verts.size());ii++){
			meshtree.elemind.push_back(verts(ii)->index);
		}
		break;
		case 1:
		meshtree.elemind.reserve(edges.size());
		for (ii=0; ii<int(edges.size());ii++){
			meshtree.elemind.push_back(edges(ii)->index);
		}
		break;
		case 2:
		meshtree.elemind.reserve(surfs.size());
		for (ii=0; ii<int(surfs.size());ii++){
			meshtree.elemind.push_back(surfs(ii)->index);
		}
		break;
		case 3:
		meshtree.elemind.reserve(volus.size());
		for (ii=0; ii<int(volus.size());ii++){
			meshtree.elemind.push_back(volus(ii)->index);
		}
		break;
	}
	meshDepIsSet=true;

}

void mesh::ReturnParentMap(vector<int> &currind, vector<int> &parentpos,
	vector<pair<int,int>> &parentcases, vector<double> &voluVals) const {

	int ii, ni, jj, nj, nElm, nParCases;
	pair<int,int> tempPair;

	currind.clear();
	parentpos.clear();
	parentcases.clear();

	nParCases=0;
	ni=meshtree.nParents;
	nj=meshtree.elemind.size();
	nElm=meshtree.elemind.size();

	currind.reserve(nElm*meshtree.nParents);
	parentpos.reserve(nElm*meshtree.nParents);
	voluVals.reserve(nElm*meshtree.nParents);
	parentcases.reserve(this->CountVoluParent());

	for (ii=0; ii< ni; ++ii){
		// Build the list of parent cases
		nj=meshtree.parentmesh[ii]->volus.size();
		tempPair.first=ii;
		for (jj = 0; jj < nj; ++jj){
			tempPair.second=meshtree.parentmesh[ii]->volus(jj)->index;
			parentcases.push_back(tempPair);
			voluVals.push_back(meshtree.parentmesh[ii]->volus(jj)->volume
				* meshtree.parentmesh[ii]->volus(jj)->target);
		}
		for (jj = 0; jj < nElm; ++jj){
			currind.push_back(meshtree.elemind[jj]);
			parentpos.push_back(meshtree.parentmesh[ii]->volus.find(
				meshtree.parentconn[ii][jj])+nParCases);
			
		}
		nParCases+=nj;
	}
}

void mesh::MapVolu2Parent(const vector<double> &fillIn, const vector<pair<int,int>> &parentcases,
	double volu::*mp){
	int ii, ni, sub, sub2; 

	ni=parentcases.size();
	// cout << endl << "new fills: ";
	//DisplayVector(fillIn);
	// cout << endl << "replacing: " << ni << " | " ;
	for(ii=0; ii< ni; ii++){
		sub2=this->meshtree.parentmesh[parentcases[ii].first]->volus.isHash;
		sub=this->meshtree.parentmesh[parentcases[ii].first]
		->volus.find(parentcases[ii].second);
		// cout << "(" << this->meshtree.parentmesh[parentcases[ii].first]
		// 	->volus[sub].fill << " , ";
		this->meshtree.parentmesh[parentcases[ii].first]
		->volus[sub].*mp=fillIn[ii];

		// cout  << this->meshtree.parentmesh[parentcases[ii].first]
		// 	->volus[sub].fill << ") ";
		this->meshtree.parentmesh[parentcases[ii].first]->volus.isHash=sub2;
	}
	// cout<< endl;

}
void mesh::MapVolu2Self(const vector<double> &fillIn, const vector<int> &elms,
	double volu::*mp){
	int ii, ni, sub; 

	ni=elms.size();
	sub=volus.isHash;
	for(ii=0; ii< ni; ii++){
		volus[volus.find(elms[ii])].*mp=fillIn[ii];
		volus.isHash=sub;
	}
	// cout<< endl;

}

void mesh::MaintainLineage(){
	// Method not implemented yet, features:
	//	- recognise the modifications needed depending on child and parent dimensionality 
	//	- Modify the parentconn vec 

}
int mesh::CountParents() const {
	return(meshtree.parentmesh.size());
}
int mesh::SurfInParent(int surfind) const{
	int kk1,kk2;
	int jj=-1;
	int ii = surfs.find(surfind);
	bool isParentSurf=false;
	if(ii>=0){
		kk1=volus.find(surfs(ii)->voluind[0]);
		kk2=volus.find(surfs(ii)->voluind[1]);
		while(!isParentSurf && jj<meshtree.nParents-1){
			++jj;
			if((kk1!=-1) ^ (kk2!=-1)){
				isParentSurf=true;
			} else if (kk1!=-1){
				isParentSurf= meshtree.parentconn[jj][kk1] != meshtree.parentconn[jj][kk2];
			}
		}
	}
	if (isParentSurf){
		return(jj);
	} else {
		return(-1);
	}
}

void mesh::SurfInParent(vector<int> &listInParent) const{
	int ii, nSurf = surfs.size();
	listInParent.clear();
	for (ii=0; ii< nSurf; ++ii){
		if (SurfInParent(surfs(ii)->index)>=0){
			listInParent.push_back((ii));
		}
	} 
}

void mesh::SurfValuesofParents(int elmInd, vector<double> &vals, int volType) const{
		/*
	Extracts the surfaces in the parents corresponding to element elmInd
	volType is a selector for which value to extract
	*/
	double surf::*mp;
	mp=NULL;
	switch(volType){
		case 0:
		mp=&surf::area;
		break;
		case 1:
		mp=&surf::fill;
		break;
		case 2:
		mp=&surf::target;
		break;
		case 3:
		mp=&surf::error;
		break;
	}

	this->SurfValuesofParents(elmInd, vals, mp);
}

void mesh::SurfValuesofParents(int elmInd, vector<double> &vals, double surf::*mp) const{
	/*
	Extracts the volumes in the parents corresponding to element elmInd
	volType is a selector for which value to extract
	*/
	int sub, ii;
	vals.clear();


	sub = surfs.find(elmInd);
	for (ii=0; ii < meshtree.nParents; ++ii){
		vals.push_back(meshtree.parentmesh[ii]->surfs.isearch
			(meshtree.parentconn[ii][sub])->*mp);

	}

}

void mesh::VoluValuesofParents(int elmInd, vector<double> &vals, int volType) const{
		/*
	Extracts the volumes in the parents corresponding to element elmInd
	volType is a selector for which value to extract
	*/
	
	double volu::*mp;
	mp=NULL;
	switch(volType){
		case 0:
		mp=&volu::volume;
		break;
		case 1:
		mp=&volu::fill;
		break;
		case 2:
		mp=&volu::target;
		break;
		case 3:
		mp=&volu::error;
		break;
	}

	this->VoluValuesofParents(elmInd, vals, mp);
}

void mesh::VoluValuesofParents(int elmInd, vector<double> &vals, double volu::*mp) const{
	/*
	Extracts the volumes in the parents corresponding to element elmInd
	volType is a selector for which value to extract
	*/
	int sub, ii;
	vals.clear();


	sub = volus.find(elmInd);
	for (ii=0; ii < meshtree.nParents; ++ii){
		vals.push_back(meshtree.parentmesh[ii]->volus.isearch
			(meshtree.parentconn[ii][sub])->*mp);

	}

}

void mesh::ElmOnParentBound(vector<int> &listInParent, vector<int> &voluInd,
	bool isBorderBound,
	bool outerVolume) const{

	if (this->WhatDim()==3){
		this->SurfOnParentBound(listInParent, voluInd,
			isBorderBound,
			outerVolume);
	} else if (this->WhatDim()==2) {
		this->EdgeOnParentBound(listInParent, voluInd,
			isBorderBound,
			outerVolume);
	}
}

void mesh::SurfOnParentBound(vector<int> &listInParent, vector<int> &voluInd,
	bool isBorderBound,
	bool outerVolume) const{
	/*
	Returns a list of surfaces which are on the outer or inner boundaries.
	*/
	int ii, jj, boundVolume, nSurf = surfs.size();
	bool isOnBound;
	vector<double> vals0, vals1;
	listInParent.clear();

	if(outerVolume){
		boundVolume = 0.0;
	} else {
		boundVolume = 1.0;
	}// Mesh component comparison

	vals0.reserve(meshtree.nParents);
	vals1.reserve(meshtree.nParents);
	voluInd.clear();
	voluInd.reserve(volus.size());

	for (ii=0; ii< nSurf; ++ii){
		isOnBound=false;
		if (surfs(ii)->voluind[0]==0 || surfs(ii)->voluind[1]==0){
			isOnBound=isBorderBound;
			if(isOnBound && (surfs(ii)->voluind[0]!=0)){
				voluInd.push_back(surfs(ii)->voluind[0]);
			}
			if(isOnBound && (surfs(ii)->voluind[1]!=0)){
				voluInd.push_back(surfs(ii)->voluind[1]);
			}
		} else if (SurfInParent(surfs(ii)->index)>=0){
			// Check parent volume values
			this->VoluValuesofParents(surfs(ii)->voluind[0],vals0, &volu::target);
			this->VoluValuesofParents(surfs(ii)->voluind[1],vals1, &volu::target);
			jj=0;
			// DisplayVector(vals0);
			// DisplayVector(vals1);
			while(!isOnBound && jj< meshtree.nParents){

				isOnBound = (vals0[jj]!=vals1[jj]) 
				&& ((vals0[jj]==boundVolume) || (vals1[jj]==boundVolume));
				++jj;
			}
			
			if(isOnBound && (vals1[jj-1]==boundVolume)){
				voluInd.push_back(surfs(ii)->voluind[1]);
				if(boundVolume==1.0 && volus.isearch(voluInd.back())->isBorder){
					voluInd.pop_back();
				}
			}
			if(isOnBound && (vals0[jj-1]==boundVolume)){
				voluInd.push_back(surfs(ii)->voluind[0]);
				if(boundVolume==1.0 && volus.isearch(voluInd.back())->isBorder){
					voluInd.pop_back();
				}
			}
			
			//cout << "? " << isOnBound << " || ";
		}
		if (isOnBound){
			listInParent.push_back(surfs(ii)->index);
		}
		
	} 
}


void mesh::EdgeOnParentBound(vector<int> &listInParent, vector<int> &surfInd,
	bool isBorderBound,
	bool outerVolume) const{
	/*
	Returns a list of edges which are on the outer or inner boundaries.
	*/
	int ii, jj, boundVolume, nSurf = edges.size();
	bool isOnBound;
	vector<double> vals0, vals1;
	listInParent.clear();

	if(outerVolume){
		boundVolume = 0.0;
	} else {
		boundVolume = 1.0;
	}// Mesh component comparison

	vals0.reserve(meshtree.nParents);
	vals1.reserve(meshtree.nParents);
	surfInd.clear();
	surfInd.reserve(surfs.size());

	for (ii=0; ii< nSurf; ++ii){
		isOnBound=false;
		if (edges(ii)->surfind[0]==0 || edges(ii)->surfind[1]==0){
			isOnBound=isBorderBound;
			if(isOnBound && (edges(ii)->surfind[0]!=0)){
				surfInd.push_back(edges(ii)->surfind[0]);
			}
			if(isOnBound && (edges(ii)->surfind[1]!=0)){
				surfInd.push_back(edges(ii)->surfind[1]);
			}
		} else if (SurfInParent(edges(ii)->index)>=0){
			// Check parent surfme values
			this->SurfValuesofParents(edges(ii)->surfind[0],vals0, &surf::target);
			this->SurfValuesofParents(edges(ii)->surfind[1],vals1, &surf::target);
			jj=0;
			// DisplayVector(vals0);
			// DisplayVector(vals1);
			while(!isOnBound && jj< meshtree.nParents){

				isOnBound = (vals0[jj]!=vals1[jj]) 
				&& ((vals0[jj]==boundVolume) || (vals1[jj]==boundVolume));
				++jj;
			}
			
			if(isOnBound && (vals1[jj-1]==boundVolume)){
				surfInd.push_back(edges(ii)->surfind[1]);
				if(boundVolume==1.0 && surfs.isearch(surfInd.back())->isBorder){
					surfInd.pop_back();
				}
			}
			if(isOnBound && (vals0[jj-1]==boundVolume)){
				surfInd.push_back(edges(ii)->surfind[0]);
				if(boundVolume==1.0 && surfs.isearch(surfInd.back())->isBorder){
					surfInd.pop_back();
				}
			}
			
			//cout << "? " << isOnBound << " || ";
		}
		if (isOnBound){
			listInParent.push_back(edges(ii)->index);
		}
		
	} 
}

int mesh::CountVoluParent() const {
	int n=0;
	for (int i = 0; i < int(meshtree.parentmesh.size()); ++i)
	{
		n+=meshtree.parentmesh[i]->volus.size();
	}
	return(n);
}

int mesh::ParentElementIndex(int childElmInd, int parentInd) const{
	/*
	 Returns the parent of an element.
	 */
	if(parentInd>=int(this->meshtree.parentconn.size())
		|| parentInd<0){
		RSVS3D_ERROR_ARGUMENT("parentInd is larger than the number of mesh parents");
	}
	if (this->WhatDim()==3){
		return (this->meshtree.parentconn[parentInd][this->volus.find(childElmInd)]);
	} else if (this->WhatDim()==2){
		return (this->meshtree.parentconn[parentInd][this->surfs.find(childElmInd)]);
	} else {
		RSVS3D_ERROR_ARGUMENT("Dimensionality other than 2 or 3 not supported yet");
		return(-1);
	}

}

/// MAth operations in mesh
void edge::GeometricProperties(const mesh *meshin, coordvec &centre, double &length) const {

	int sub;
	centre.assign(0,0,0);
	length=0;
	sub=meshin->verts.find(vertind[0]);
	centre.substract(meshin->verts(sub)->coord);
	centre.add(meshin->verts.isearch(vertind[1])->coord);
	length=centre.CalcNorm(); 
	centre.add(meshin->verts(sub)->coord);
	centre.add(meshin->verts(sub)->coord);
	centre.div(2);
}

/**
 * @brief      Calculate squared edge length
 *
 * @param[in]  meshin  the mesh in which the edge existes
 *
 * @return     the squared length of the edge
 */
double edge::LengthSquared(const mesh &meshin) const {
	double lengthSquared=0.0;

	for (int i = 0; i < meshin.WhatDim(); ++i)
	{
		lengthSquared += pow(
				meshin.verts.isearch(this->vertind[0])->coord[i]
				-meshin.verts.isearch(this->vertind[1])->coord[i]
			,2.0);
	}

	return lengthSquared;
}
/**
 * @brief      Calculate the edge length
 *
 * @param[in]  meshin  the mesh in which the edge existes
 *
 * @return     the length of the edge
 */
double edge::Length(const mesh &meshin) const {
	return sqrt(this->LengthSquared(meshin));
}
/**
 * @brief      Returns
 *
 * @param[in]  meshin  the mesh in which the edge existes
 * @param[in]  eps     Tolerance, number under which the length must 
 * 						be to be considered 0. Defaults to __DBL_EPSILON__.
 *
 * @return     Wether Length squared is below eps squared.
 */
bool edge::IsLength0(const mesh &meshin, double eps) const {
	return this->LengthSquared(meshin)< pow(eps,2.0);
}

// Methods of meshpart : volu surf edge vert
void volu::disp() const{
	int i;
	cout << "volu : index " << index << " | fill " << fill << ", "  <<
	target << ", "<< error << " | isBorder " << isBorder << " | surfind " << surfind.size();
	if (surfind.size()<30){
		for (i=0; unsigned_int(i)<surfind.size();i++){
			cout << "-" << surfind[i];
		}
	}
	cout << endl;
}

void surf::disp() const{
	int i;
	cout << "surf : index " << index << " | fill " << fill  << ", " << 
	target << ", "<< error << " | isBorder " << isBorder << " | voluind " << voluind.size();
	for (i=0; unsigned_int(i)<voluind.size();i++){
		cout << "-" << voluind[i];
	}
	cout << " | edgeind " << edgeind.size();
	for (i=0; unsigned_int(i)<edgeind.size();i++){
		cout << "-" << edgeind[i];
	}
	cout << endl;
}

void edge::disp() const{
	int i;
	cout << "edge : index " << index << " | isBorder " << isBorder << " | vertind " << vertind.size();
	for (i=0; unsigned_int(i)<vertind.size();i++){
		cout << "-" << vertind[i];
	}
	cout << " | surfind " << surfind.size();
	for (i=0; unsigned_int(i)<surfind.size();i++){
		cout << "-" << surfind[i];
	}
	cout << endl;
}

void vert::disp() const{
	int i;
	cout << "vert : index " << index << " | isBorder " << isBorder << " | edgeind " << edgeind.size();
	for (i=0; unsigned_int(i)<edgeind.size();i++){
		cout << "-" << edgeind[i];
	}
	cout << " | coord " << coord.size();
	for (i=0; unsigned_int(i)<coord.size();i++){
		cout << "-" << coord[i];
	}
	cout << endl;
}

// Class function definitions
// Methods of meshpart : volu surf edge vert
void volu::disptree(const mesh &meshin, int n) const{
	int i;
	for(i=0;i<n;i++){cout<< "  ";}
		disp();
	if(n>0 && surfind.size()<8){
		for (i=0; unsigned_int(i)<surfind.size();i++){
			meshin.surfs.isearch(surfind[i])->disptree(meshin,n-1);
		}
	}
}

void surf::disptree(const mesh &meshin, int n) const{
	int i;
	for(i=0;i<n;i++){cout<< "  ";}
		disp();
	if(n>0){
		for (i=0; unsigned_int(i)<edgeind.size();i++){
			meshin.edges.isearch(edgeind[i])->disptree(meshin,n-1);
		}
		for (i=0; unsigned_int(i)<voluind.size();i++){
			if(voluind[i]>0){
				meshin.volus.isearch(voluind[i])->disptree(meshin,n-1);
			}
		}
	}
}

void edge::disptree(const mesh &meshin, int n) const{
	int i;
	for(i=0;i<n;i++){cout<< "  ";}
		disp();
	if(n>0){
		for (i=0; unsigned_int(i)<vertind.size();i++){
			meshin.verts.isearch(vertind[i])->disptree(meshin,n-1);
		}
		for (i=0; unsigned_int(i)<surfind.size();i++){
			meshin.surfs.isearch(surfind[i])->disptree(meshin,n-1);
		}
	}
}

void vert::disptree(const mesh &meshin, int n) const{
	int i;
	for(i=0;i<n;i++){cout<< "  ";}
		disp();
	if(n>0){
		for (i=0; unsigned_int(i)<edgeind.size();i++){
			meshin.edges.isearch(edgeind[i])->disptree(meshin,n-1);
		}
	}
}

// Input and output
void volu::write(FILE *fid) const{

	int i;

	fprintf(fid, "%i %.16lf %.16lf %.16lf %i ",index,fill,target, error,int(isBorder));
	fprintf(fid, "%i ",int(surfind.size()));
	for (i=0; unsigned_int(i)<surfind.size();i++){
		fprintf(fid, "%i ",surfind[i]);
	}
	fprintf(fid,"\n");
}

void surf::write(FILE * fid) const{
	int i;

	fprintf(fid, "%i %.16lf %.16lf %.16lf %i ",index,fill,target, error,int(isBorder));
	fprintf(fid, "%i ",int(voluind.size()));
	for (i=0; unsigned_int(i)<voluind.size();i++){
		fprintf(fid, "%i ",voluind[i]);
	}
	fprintf(fid, "%i ",int(edgeind.size()));
	for (i=0; unsigned_int(i)<edgeind.size();i++){
		fprintf(fid, "%i ",edgeind[i]);
	}
	fprintf(fid,"\n");
}

void edge::write(FILE * fid) const{
	int i;

	fprintf(fid, "%i %i ",index,int(isBorder));
	fprintf(fid, "%i ",int(vertind.size()));
	for (i=0; unsigned_int(i)<vertind.size();i++){
		fprintf(fid, "%i ",vertind[i]);
	}
	fprintf(fid, "%i ",int(surfind.size()));
	for (i=0; unsigned_int(i)<surfind.size();i++){
		fprintf(fid, "%i ",surfind[i]);
	}
	fprintf(fid,"\n");
}

void vert::write(FILE * fid) const{
	int i;

	fprintf(fid, "%i %i ",index,int(isBorder));
	fprintf(fid, "%i ",int(edgeind.size()));
	for (i=0; unsigned_int(i)<edgeind.size();i++){
		fprintf(fid, "%i ",edgeind[i]);
	}
	fprintf(fid, "%i ",int(coord.size()));
	for (i=0; unsigned_int(i)<coord.size();i++){
		fprintf(fid, "%.16lf ",coord[i]);
	}
	fprintf(fid,"\n");
}

void volu::read(FILE * fid) {

	int i,n;

	fscanf(fid, "%i %lf %lf %lf %i ",&index,&fill,&target, &error, &i);
	isBorder=bool(i);
	fscanf(fid, "%i ",&n);
	surfind.assign(n,0);
	for (i=0; unsigned_int(i)<surfind.size();i++){
		fscanf(fid, "%i ",&surfind[i]);
	}

}

void surf::read(FILE * fid) {
	int i,n;

	fscanf(fid, "%i %lf %lf %lf %i ",&index,&fill,&target, &error, &i);
	isBorder=bool(i);
	fscanf(fid, "%i ",&n);
	voluind.assign(n,0);
	for (i=0; unsigned_int(i)<voluind.size();i++){
		fscanf(fid, "%i ",&voluind[i]);
	}
	fscanf(fid, "%i ",&n);
	edgeind.assign(n,0);
	for (i=0; unsigned_int(i)<edgeind.size();i++){
		fscanf(fid, "%i ",&edgeind[i]);
	}

}

void edge::read(FILE * fid) {
	int i,n;

	fscanf(fid, "%i %i ",&index, &i);
	isBorder=bool(i);
	fscanf(fid, "%i ",&n);
	vertind.assign(n,0);
	for (i=0; unsigned_int(i)<vertind.size();i++){
		fscanf(fid, "%i ",&vertind[i]);
	}
	fscanf(fid, "%i ",&n);
	surfind.assign(n,0);
	for (i=0; unsigned_int(i)<surfind.size();i++){
		fscanf(fid, "%i ",&surfind[i]);
	}

}

void vert::read(FILE * fid) {
	int i,n;

	fscanf(fid, "%i %i ",&index, &i);
	isBorder=bool(i);
	fscanf(fid, "%i ",&n);
	edgeind.assign(n,0);
	for (i=0; unsigned_int(i)<edgeind.size();i++){
		fscanf(fid, "%i ",&edgeind[i]);
	}
	fscanf(fid, "%i ",&n);
	coord.assign(n,0);
	for (i=0; unsigned_int(i)<coord.size();i++){
		fscanf(fid, "%lf ",&coord[i]);
	}

}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
void volu::ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu){
	int i;
	index+=nVolu;
	for (i=0; unsigned_int(i)<surfind.size();i++){
		surfind[i]= (surfind[i]>0)? (surfind[i]+nSurf) : surfind[i];
	}
}
void surf::ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu){
	int i;
	index+=nSurf;
	for (i=0; unsigned_int(i)<voluind.size();i++){
		voluind[i]= (voluind[i]>0) ? (voluind[i]+nVolu) : voluind[i];
	}

	for (i=0; unsigned_int(i)<edgeind.size();i++){
		edgeind[i]=edgeind[i]+nEdge;
	}
}
void edge::ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu){
	int i;
	index+=nEdge;
	for (i=0; unsigned_int(i)<vertind.size();i++){
		vertind[i]=vertind[i]+nVert;
	}

	for (i=0; unsigned_int(i)<surfind.size();i++){
		surfind[i]=(surfind[i]>0) ? (surfind[i]+nSurf) : surfind[i];
	}
}
void vert::ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu){
	int i;
	index+=nVert;
	for (i=0; unsigned_int(i)<edgeind.size();i++){
		edgeind[i]=edgeind[i]+nEdge;
	}
}
void mesh::TightenConnectivity(){
	verts.TightenConnectivity();
	surfs.TightenConnectivity();
	edges.TightenConnectivity();
	volus.TightenConnectivity();
}

void mesh::SwitchIndex(int typeInd, int oldInd, int newInd, const vector<int> &scopeInd)
{
	/*Switch the indices of an element for another in all the appropriate
	connectivity lists.

	Args:
		typeInd [1-8]: type of index */
	int ii,jj,kk,newSub,oldSub;
	vector<int> subList;
	HashedVector<int,int> tempSub;
	bool is3DMesh=volus.size()>0;

	if(typeInd==1){ // Switch out a vertex
		newSub=verts.find(newInd);
		oldSub=verts.find(oldInd);
		subList=edges.find_list(verts[oldSub].edgeind);
		for (ii=0;ii<int(subList.size());++ii){ // update vertind
			jj=edges(subList[ii])->vertind[1]==oldInd;
			edges[subList[ii]].vertind[jj]=newInd;
			verts[newSub].edgeind.push_back(edges[subList[ii]].index); // update vertex edgeind
			//cout << " " << edges[subList[ii]].index <<  " ";
			for (jj=0;jj<int(verts(oldSub)->edgeind.size());++jj){
				if(verts(oldSub)->edgeind[jj]==edges[subList[ii]].index){
					verts[oldSub].edgeind.erase(
						verts[oldSub].edgeind.begin()+jj);
					jj--;
				}
			}
		}
		sort(verts[newSub].edgeind);
		unique(verts[newSub].edgeind);
		if(meshDepIsSet && meshDim==0){
			// for (ii=0;ii<int(oldSub.size());++ii){
			// 	meshtree.elemind[oldSub[ii]]=newInd;
			// }
			meshtree.elemind[oldSub]=newInd;
		}
		// Hashing has not been invalidated
		edges.isHash=1;
		verts.isHash=1;

	} else if (typeInd==2){
		newSub=edges.find(newInd);
		oldSub=edges.find(oldInd);
		subList=verts.find_list(edges(edges.find(oldInd))->vertind);
		for (ii=0;ii<int(subList.size());++ii){
			for (jj=0;jj<int(verts(subList[ii])->edgeind.size());++jj){
				if(verts(subList[ii])->edgeind[jj]==oldInd){
			   //cout << " " << verts(subList[ii])->index << " " << endl;DisplayVector(verts[subList[ii]].edgeind);cout << endl;
					verts[subList[ii]].edgeind[jj]=newInd;
			   //DisplayVector(verts[subList[ii]].edgeind);cout << endl;

				}
			}  
		}
	  // Changes the indices of 
		subList=surfs.find_list(edges(edges.find(oldInd))->surfind);
		for (ii=0;ii<int(subList.size());++ii){
			if(subList[ii]!=-1 || is3DMesh){
				for (jj=0;jj<int(surfs(subList[ii])->edgeind.size());++jj){
					if(surfs(subList[ii])->edgeind[jj]==oldInd){
						surfs[subList[ii]].edgeind[jj]=newInd;
						edges[newSub].surfind.push_back(surfs[subList[ii]].index);
						surfs[subList[ii]].isordered=false;
					}
				}
			}  
		}

		if(meshDepIsSet && meshDim==1){
	  	// for (ii=0;ii<int(oldSub.size());++ii){
	  	// 	meshtree.elemind[oldSub[ii]]=newInd;
	  	// }
			meshtree.elemind[oldSub]=newInd;
		}

		edges.isHash=1;
		verts.isHash=1;
		surfs.isHash=1;
		edges.isSetMI=1;
		verts.isSetMI=1;
		surfs.isSetMI=1;


	} else if (typeInd==3){
		newSub=surfs.find(newInd);
		oldSub=surfs.find(oldInd);
		subList=edges.find_list(surfs(surfs.find(oldInd))->edgeind);
		for (ii=0;ii<int(subList.size());++ii){
			for (jj=0;jj<int(edges(subList[ii])->surfind.size());++jj){
				if(edges(subList[ii])->surfind[jj]==oldInd){
					edges[subList[ii]].surfind[jj]=newInd;
					surfs[newSub].edgeind.push_back(edges[subList[ii]].index);
				}
			}  
		}
		surfs.isHash=1;

		subList=volus.find_list(surfs(surfs.find(oldInd))->voluind);
		for (ii=0;ii<int(subList.size());++ii){
			if(subList[ii]!=-1){
				for (jj=0;jj<int(volus(subList[ii])->surfind.size());++jj){
					if(volus(subList[ii])->surfind[jj]==oldInd){
						volus[subList[ii]].surfind[jj]=newInd;
						surfs[newSub].voluind.push_back(volus[subList[ii]].index);
					}
				}
			}  
		}

		if(meshDepIsSet && meshDim==2){
	  	// for (ii=0;ii<int(oldSub.size());++ii){
	  	// 	meshtree.elemind[oldSub[ii]]=newInd;
	  	// }
			meshtree.elemind[oldSub]=newInd;
		}

		surfs[newSub].isordered=false;
		surfs.isHash=1;
		edges.isHash=1;
		volus.isHash=1;
		surfs.isSetMI=1;
		edges.isSetMI=1;
		volus.isSetMI=1;

	} else if (typeInd==4){
		newSub=volus.find(newInd);
		oldSub=volus.find(oldInd);
		subList=surfs.find_list(volus[volus.find(oldInd)].surfind);
	  for (ii=0;ii<int(subList.size());++ii){ // update vertind
		 //jj=surfs(subList[ii])->voluind[1]==oldInd;
		 //surfs[subList[ii]].voluind[jj]=newInd;
	  	for (jj=0;jj<int(surfs(subList[ii])->voluind.size());++jj){
	  		if(surfs[subList[ii]].voluind[jj]==oldInd){
	  			surfs[subList[ii]].voluind[jj]=newInd; 
				  volus[newSub].surfind.push_back(surfs[subList[ii]].index); // update vertex edgeind
				}
			}

		}
		if(meshDepIsSet && meshDim==3){
	  	// for (ii=0;ii<int(oldSub.size());++ii){
	  	// 	meshtree.elemind[oldSub[ii]]=newInd;
	  	// }
			meshtree.elemind[oldSub]=newInd;
		}
	  // Hashing has not been invalidated
		volus.isHash=1;
		surfs.isHash=1;
		volus.isSetMI=1;
		surfs.isSetMI=1;
   } else if (typeInd==5){ // Modify vertex index in scoped mode

		newSub=verts.find(newInd);
		oldSub=verts.find(oldInd);

		subList=edges.find_list(scopeInd);
		tempSub.vec=edges.find_list(verts(oldSub)->edgeind);
		tempSub.GenerateHash();
		for (ii=0;ii<int(subList.size());++ii){
			if(tempSub.find(subList[ii])!=-1){
				for (jj=0;jj<int(edges(subList[ii])->vertind.size());++jj){
					if(edges(subList[ii])->vertind[jj]==oldInd){
						edges[subList[ii]].vertind[jj]=newInd;
						verts[newSub].edgeind.push_back(edges[subList[ii]].index); 

				  //cout << " " << edges[subList[ii]].index <<  " ";
						for (kk=0;kk<int(verts(oldSub)->edgeind.size());++kk){
							if(verts(oldSub)->edgeind[kk]==edges[subList[ii]].index){
								verts[oldSub].edgeind.erase(
									verts[oldSub].edgeind.begin()+kk);
								kk--;
							}
						}
					}
				}
			}  
		}
		if(meshDepIsSet && meshDim==0){
			// for (ii=0;ii<int(oldSub.size());++ii){
			// 	meshtree.elemind[oldSub[ii]]=newInd;
			// }
			meshtree.elemind[oldSub]=newInd;
		}
		edges.isHash=1;
		verts.isHash=1;
		edges.isSetMI=1;
		verts.isSetMI=1;
	} else if (typeInd==6){ // Modify surface index in scoped mode

		newSub=surfs.find(newInd);
		oldSub=surfs.find(oldInd);

		subList=edges.find_list(scopeInd);
		for (ii=0;ii<int(subList.size());++ii){
			for (jj=0;jj<int(edges(subList[ii])->surfind.size());++jj){
				if(edges(subList[ii])->surfind[jj]==oldInd){
					edges[subList[ii]].surfind[jj]=newInd;
					surfs[newSub].edgeind.push_back(edges[subList[ii]].index); 

			   //cout << " " << edges[subList[ii]].index <<  " ";
					for (kk=0;kk<int(surfs(oldSub)->edgeind.size());++kk){
						if(surfs(oldSub)->edgeind[kk]==edges[subList[ii]].index){
							surfs[oldSub].edgeind.erase(
								surfs[oldSub].edgeind.begin()+kk);
							kk--;
						}
					}
				}
			} 
		}
		for(ii=0;ii<int(surfs(newSub)->voluind.size());ii++){
			if(surfs(newSub)->voluind[ii]>0){
				volus.elems[volus.find(surfs(newSub)->voluind[ii])].surfind.push_back(newInd);
			}
		}
		if(meshDepIsSet && meshDim==2){
			// for (ii=0;ii<int(oldSub.size());++ii){
			// 	meshtree.elemind[oldSub[ii]]=newInd;
			// }
			meshtree.elemind[oldSub]=newInd;
		}
		edges.isHash=1;
		surfs.isHash=1;
		edges.isSetMI=1;
		surfs.isSetMI=1;
   	} else if (typeInd==7){ // Modify volume index in scoped mode

		newSub=volus.find(newInd);
		oldSub=volus.find(oldInd);

		subList=surfs.find_list(scopeInd);
		for (ii=0;ii<int(subList.size());++ii){
			for (jj=0;jj<int(surfs(subList[ii])->voluind.size());++jj){
				if(surfs(subList[ii])->voluind[jj]==oldInd){
					surfs[subList[ii]].voluind[jj]=newInd;
					volus[newSub].surfind.push_back(surfs[subList[ii]].index); 

			   //cout << " " << surfs[subList[ii]].index <<  " ";
					for (kk=0;kk<int(volus(oldSub)->surfind.size());++kk){
						if(volus(oldSub)->surfind[kk]==surfs[subList[ii]].index){
							volus[oldSub].surfind.erase(
								volus[oldSub].surfind.begin()+kk);
							kk--;
						}
					}
				}
			} 
		}

		if(meshDepIsSet && meshDim==3){
			// for (ii=0;ii<int(oldSub.size());++ii){
			// 	meshtree.elemind[oldSub[ii]]=newInd;
			// }
			meshtree.elemind[oldSub]=newInd;
		}
		surfs.isHash=1;
		volus.isHash=1;
		surfs.isSetMI=1;
		volus.isSetMI=1;
   	} else {

		cerr << "Error unknown type " << typeInd << " of object for index switching" <<endl;
		cerr << "Error in " << __PRETTY_FUNCTION__ << endl;
		RSVS3D_ERROR_ARGUMENT(" Type is out of range");
	}
}

void mesh::RemoveIndex(int typeInd, int oldInd)
{

	int ii,jj, sub;
	vector<int> subList;

	if(typeInd==1){

		sub = verts.find(oldInd);
		subList=edges.find_list(verts(sub)->edgeind);
		for (ii=0;ii<int(subList.size());++ii){
			for (jj=0;jj<int(edges(subList[ii])->vertind.size());++jj){
				if(edges(subList[ii])->vertind[jj]==oldInd){
					edges[subList[ii]].vertind.erase(
						edges[subList[ii]].vertind.begin()+jj);
					jj--;
				}
			}  
		}
		edges.isHash=1;
		verts.isHash=1;
		edges.isSetMI=1;
		verts.isSetMI=1;

	} else if (typeInd==2){

		sub = edges.find(oldInd);
		subList=verts.find_list(edges(sub)->vertind);
		for (ii=0;ii<int(subList.size());++ii){
			for (jj=0;jj<int(verts(subList[ii])->edgeind.size());++jj){
				if(verts(subList[ii])->edgeind[jj]==oldInd){
					verts[subList[ii]].edgeind.erase(
						verts[subList[ii]].edgeind.begin()+jj);
					jj--;
				}
			}  
		}

		subList=surfs.find_list(edges(sub)->surfind);
		for (ii=0;ii<int(subList.size());++ii){
			if(subList[ii]!=-1){
				for (jj=0;jj<int(surfs(subList[ii])->edgeind.size());++jj){
					if(surfs(subList[ii])->edgeind[jj]==oldInd){
						surfs[subList[ii]].edgeind.erase(
							surfs[subList[ii]].edgeind.begin()+jj);
						surfs[subList[ii]].isordered=false;
						jj--;
					}
				}
			}  
		}
		edges[sub].vertind.clear();
		edges[sub].surfind.clear();

		edges.isHash=1;
		verts.isHash=1;
		surfs.isHash=1;
		edges.isSetMI=1;
		verts.isSetMI=1;
		surfs.isSetMI=1;


	} else if (typeInd==3){

		sub = surfs.find(oldInd);

		subList=edges.find_list(surfs(sub)->edgeind);
		for (ii=0;ii<int(subList.size());++ii){
			for (jj=0;jj<int(edges(subList[ii])->surfind.size());++jj){
				if(edges(subList[ii])->surfind[jj]==oldInd){
					edges[subList[ii]].surfind.erase(
						edges[subList[ii]].surfind.begin()+jj);
					jj--;
				}
			}  
		}

		subList=volus.find_list(surfs(sub)->voluind);
		for (ii=0;ii<int(subList.size());++ii){
			if(subList[ii]!=-1){
				for (jj=0;jj<int(volus(subList[ii])->surfind.size());++jj){
					if(volus(subList[ii])->surfind[jj]==oldInd){
						volus[subList[ii]].surfind.erase(
							volus[subList[ii]].surfind.begin()+jj);
						jj--;
					}
				}
			}  
		}

		surfs[sub].edgeind.clear();
		surfs[sub].voluind.clear();
		surfs.isHash=1;
		edges.isHash=1;
		volus.isHash=1;
		surfs.isSetMI=1;
		edges.isSetMI=1;
		volus.isSetMI=1;

	} else if (typeInd==4){

		sub = volus.find(oldInd);

		subList=surfs.find_list(volus(sub)->surfind);
		for (ii=0;ii<int(subList.size());++ii){
			for (jj=0;jj<int(surfs(subList[ii])->voluind.size());++jj){
				if(surfs(subList[ii])->voluind[jj]==oldInd){
					surfs[subList[ii]].voluind.erase(
						surfs[subList[ii]].voluind.begin()+jj);
					if(surfs[subList[ii]].voluind.size()<2){
						surfs[subList[ii]].voluind.push_back(0);
					}
					jj--;
				}
			}  
		}

		volus[sub].surfind.clear();
		surfs.isHash=1;
		volus.isHash=1;
		surfs.isSetMI=1;
		volus.isSetMI=1;
   } else if (typeInd==5){ // Modify vertex index in scoped mode

   	cerr << "not coded yet" << endl;
   	throw;
   } else {

   	cerr << "Error unknown type of object for index switching" <<endl;
   	cerr << "Error in " << __PRETTY_FUNCTION__ << endl;
   	RSVS3D_ERROR_ARGUMENT(" : Type is out of range");
   }
}

int mesh::TestConnectivity(const char *strRoot) const{
	int ii,jj,kk,kk2,errCount, errTot;
	vector<int> testSub;

	errCount=0;
	errTot=0;
	kk=int(verts.size());
	for (ii=0; ii<kk;++ii){
		if(verts(ii)->edgeind.size()==0){
			errCount++;
			cerr << " Test Connectivity Error :" << errCount << " vertex " << verts(ii)->index
				<< " Has empty connectivity list "; 
			cerr << endl;

		} else {
			#ifdef RSVS_DIAGNOSTIC
			if (verts(ii)->edgeind.size()==2) {
				errCount++;
				cerr << " Test Connectivity Error :" << errCount << " vertex " 
					<< verts(ii)->index << " Has connectivity list of length 2 ";
				DisplayVector(verts(ii)->edgeind); 
				cerr << endl;
			}
			#endif
			testSub=edges.find_list(verts(ii)->edgeind);
			kk2=testSub.size();
			for(jj=0;jj< kk2; ++jj){
				if (testSub[jj]<0 && verts(ii)->edgeind[jj]!=0){
					errCount++;
					cerr << " Test Connectivity Error :" << errCount << " vertex " << verts(ii)->index
					<< " makes unknown reference to edge " << verts(ii)->edgeind[jj] << " list: " ; 
					DisplayVector(verts(ii)->edgeind);
					cerr << endl;
				}
			}
		}
	}
	if (errCount>0){
		cerr << "Test Connectivity vertex (edgeind) Errors :" << errCount << endl;
	}


	errTot+=errCount;
	errCount=0;
	kk=int(edges.size());
	for (ii=0; ii<kk;++ii){
		testSub=verts.find_list(edges(ii)->vertind);
		kk2=testSub.size();
		for(jj=0;jj< kk2; ++jj){
			if (testSub[jj]<0 && edges(ii)->vertind[jj]!=0){
				errCount++;
				cerr << " Test Connectivity Error :" << errCount << " edge " << edges(ii)->index
				<< " makes unknown reference to vertex " << edges(ii)->vertind[jj] << " list: " ; DisplayVector(edges(ii)->vertind); 
				cerr << endl;
			}
		}
	}
	if (errCount>0){
		cerr << "Test Connectivity edges (vertind) Errors :" << errCount << endl;
	}


	errTot+=errCount;
	errCount=0;
	for (ii=0; ii<kk;++ii){
		testSub=surfs.find_list(edges(ii)->surfind);
		kk2=testSub.size();
		for(jj=0;jj< kk2; ++jj){
			if (testSub[jj]<0 && edges(ii)->surfind[jj]!=0){
				errCount++;
				cerr << " Test Connectivity Error :" << errCount << " edge " << edges(ii)->index
				<< " makes unknown reference to surface " << edges(ii)->surfind[jj] << endl;
			}
		}
	}
	if (errCount>0){
		cerr << "Test Connectivity edges (surfind) Errors :" << errCount << endl;
	}



	errTot+=errCount;
	errCount=0;
	kk=int(surfs.size());
	for (ii=0; ii<kk;++ii){
		testSub=edges.find_list(surfs(ii)->edgeind);
		kk2=testSub.size();
		for(jj=0;jj< kk2; ++jj){
			if (testSub[jj]<0){
				errCount++;
				cerr << " Test Connectivity Error :" << errCount << " surf " << surfs(ii)->index
				<< " makes unknown reference to edge " << surfs(ii)->edgeind[jj] << endl;
			}
		}
		if (int(testSub.size())==0){
			errCount++;
			cerr << " Test Connectivity Error :" << errCount << " surf " << surfs(ii)->index
			<< " has empty edgeind " <<  endl;
		}
	}
	if (errCount>0) {
		cerr << "Test Connectivity surfs (edgeind) Errors :" << errCount << endl;
	}

	errTot+=errCount;
	errCount=0;
	kk=int(surfs.size());
	for (ii=0; ii<kk;++ii){
		testSub=volus.find_list(surfs(ii)->voluind);
		kk2=testSub.size();
		for(jj=0;jj< kk2; ++jj){
			if (testSub[jj]<0 && surfs(ii)->voluind[jj]!=0){
				errCount++;
				cerr << " Test Connectivity Error :" << errCount << " surf " << surfs(ii)->index
				<< " makes unknown reference to volu " << surfs(ii)->voluind[jj] << endl;
			}
		}
	}
	if (errCount>0) {
		cerr << "Test Connectivity surfs (voluind) Errors :" << errCount << endl;
	}


	errTot+=errCount;
	errCount=0;
	kk=int(volus.size());
	for (ii=0; ii<kk;++ii){
		testSub=surfs.find_list(volus(ii)->surfind);
		kk2=testSub.size();
		for(jj=0;jj< kk2; ++jj){
			if (testSub[jj]<0 && volus(ii)->surfind[jj]!=0){
				errCount++;
				cerr << " Test Connectivity Error :" << errCount << " volu " << volus(ii)->index
				<< " makes unknown reference to surface " << volus(ii)->surfind[jj] << endl;
			}
		}
	}
	if (errCount>0){
		cerr << "Test Connectivity volus (surfind) Errors :" << errCount << endl;
	}
	errTot+=errCount;
	if (errTot>0){
		cerr << errTot << "  Total errors were detected in the connectivity "
			"list in: " << endl << strRoot <<endl;
	}
	return(errTot);
}

int mesh::TestConnectivityBiDir(const char *strRoot,
	bool emptyIsErr) const{
	int ii,jj,ll,ll2,kk,kk2,errCount, errTot,errCountBiDir,errTotBiDir;
	bool flag;
	vector<int> testSub;

	errCount=0;
	errCountBiDir=0;
	errTot=0;
	errTotBiDir=0;
	kk=int(verts.size());
	for (ii=0; ii<kk;++ii){
		if(emptyIsErr && verts(ii)->edgeind.size()==0){
			errCount++;
			cerr << " Test Connectivity Error :" << errCount << " vertex " << verts(ii)->index
			<< " Has empty connectivity list "; 
			cerr << endl;

		} else {
			testSub=edges.find_list(verts(ii)->edgeind);
			kk2=testSub.size();
			for(jj=0;jj< kk2; ++jj){
				if (testSub[jj]<0 && verts(ii)->edgeind[jj]!=0){
					errCount++;
					cerr << " Test Connectivity Error :" << errCount << " vertex " << verts(ii)->index
					<< " makes unknown reference to edge " << verts(ii)->edgeind[jj] << " list: " ; 
					DisplayVector(verts(ii)->edgeind);
					cerr << endl;
				} else if (verts(ii)->edgeind[jj]!=0){
					ll2=edges(testSub[jj])->vertind.size();
					flag=false;
					for(ll=0;ll<ll2;ll++){
						flag=flag || (edges(testSub[jj])->vertind[ll]==verts(ii)->index);
					}
					if (!flag){
						errCountBiDir++;
						cerr << " Test Connectivity Error :" << errCountBiDir << " vertex " << verts(ii)->index
						<< " makes uni-directional reference to edge " << verts(ii)->edgeind[jj] << " list: " ; 
						DisplayVector(verts(ii)->edgeind); cout << " list (edge.vertind): " ; 
						DisplayVector(edges(jj)->vertind);cerr << endl;
					}
				}
			}
		}
	}
	if (errCount>0){
		cerr << "Test Connectivity vertex (edgeind) Errors :" << errCount << endl;
	}
	if (errCountBiDir>0){
		cerr << "Test Connectivity vertex (edgeind) uni-directional Errors :" << errCountBiDir << endl;
	}


	errTot+=errCount;
	errTotBiDir+=errCountBiDir;
	errCount=0;
	errCountBiDir=0;
	kk=int(edges.size());
	for (ii=0; ii<kk;++ii){
		testSub=verts.find_list(edges(ii)->vertind);
		kk2=testSub.size();
		if(emptyIsErr && testSub.size()!=2){ // Vertind should be of size 2
			errCount++;
			cerr <<  " Test Connectivity Error :"  << errCount << " edge " << edges(ii)->index
				<< " has vertind of length " << testSub.size() << " list (edge::vertind): " ;
			DisplayVector(edges(ii)->vertind); 
			cerr << endl;
		}
		for(jj=0;jj< kk2; ++jj){
			if (testSub[jj]<0 && edges(ii)->vertind[jj]!=0){
				errCount++;
				cerr << " Test Connectivity Error :" << errCount << " edge " << edges(ii)->index
				<< " makes unknown reference to vertex " << edges(ii)->vertind[jj] << " list: " ; DisplayVector(edges(ii)->vertind); 
				cerr << endl;
			} else if (edges(ii)->vertind[jj]!=0){
				ll2=verts(testSub[jj])->edgeind.size();
				flag=false;
				for(ll=0;ll<ll2;ll++){
					flag=flag || (verts(testSub[jj])->edgeind[ll]==edges(ii)->index);
				}
				if (!flag){
					errCountBiDir++;
					cerr << " Test Connectivity Error :" << errCountBiDir << " edge " << edges(ii)->index
					<< " makes uni-directional reference to vertex " << edges(ii)->vertind[jj] << " list: " ; 
					DisplayVector(edges(ii)->vertind); cout << " list (vert.edgeind): " ; 
					DisplayVector(verts(jj)->edgeind);cerr << endl;
				}
			}
		}
	}
	if (errCount>0){
		cerr << "Test Connectivity edges (vertind) Errors :" << errCount << endl;
	}
	if (errCountBiDir>0){
		cerr << "Test Connectivity edges (vertind) uni-directional  Errors :" << errCountBiDir << endl;
	}


	errTot+=errCount;
	errTotBiDir+=errCountBiDir;
	errCount=0;
	errCountBiDir=0;
	for (ii=0; ii<kk;++ii){
		testSub=surfs.find_list(edges(ii)->surfind);
		kk2=testSub.size();
		for(jj=0;jj< kk2; ++jj){
			if (testSub[jj]<0 && edges(ii)->surfind[jj]!=0){
				errCount++;
				cerr << " Test Connectivity Error :" << errCount << " edge " << edges(ii)->index
				<< " makes unknown reference to surface " << edges(ii)->surfind[jj] << endl;
			} else if (edges(ii)->surfind[jj]!=0){
				ll2=surfs(testSub[jj])->edgeind.size();
				flag=false;
				for(ll=0;ll<ll2;ll++){
					flag=flag || (surfs(testSub[jj])->edgeind[ll]==edges(ii)->index);
				}
				if (!flag){
					errCountBiDir++;
					cerr << " Test Connectivity Error :" << errCountBiDir << " edge " << edges(ii)->index
					<< " makes uni-directional reference to surface " << edges(ii)->surfind[jj] << " list: " ; 
					DisplayVector(edges(ii)->surfind); cout << " list (surf.edgeind): " ; 
					DisplayVector(surfs(jj)->edgeind);cerr << endl;
				}
			}
		}
	}
	if (errCount>0){
		cerr << "Test Connectivity edges (surfind) Errors :" << errCount << endl;
	}
	if (errCountBiDir>0){
		cerr << "Test Connectivity edges (surfind) uni-directional  Errors :" << errCountBiDir << endl;
	}



	errTot+=errCount;
	errTotBiDir+=errCountBiDir;
	errCount=0;
	errCountBiDir=0;
	kk=int(surfs.size());
	for (ii=0; ii<kk;++ii){
		testSub=edges.find_list(surfs(ii)->edgeind);
		kk2=testSub.size();
		if (int(testSub.size())==0){
			errCount++;
			cerr << " Test Connectivity Error :" << errCount << " surf " << surfs(ii)->index
			<< " has empty edgeind " <<  endl;
		}
		for(jj=0;jj< kk2; ++jj){
			if (testSub[jj]<0){
				errCount++;
				cerr << " Test Connectivity Error :" << errCount << " surf " << surfs(ii)->index
				<< " makes unknown reference to edge " << surfs(ii)->edgeind[jj] << endl;
			} else if (surfs(ii)->edgeind[jj]!=0){
				ll2=edges(testSub[jj])->surfind.size();
				flag=false;
				for(ll=0;ll<ll2;ll++){
					flag=flag || (edges(testSub[jj])->surfind[ll]==surfs(ii)->index);
				}
				if (!flag){
					errCountBiDir++;
					cerr << " Test Connectivity Error :" << errCountBiDir << " surf " << surfs(ii)->index
					<< " makes uni-directional reference to edge " << surfs(ii)->edgeind[jj] << " list: " ; 
					DisplayVector(surfs(ii)->edgeind); cout << " list (edge.surfind): " ; 
					DisplayVector(edges(jj)->surfind);cerr << endl;
				}
			}
		}
	}
	if (errCount>0) {
		cerr << "Test Connectivity surfs (edgeind) Errors :" << errCount << endl;
	}
	if (errCountBiDir>0) {
		cerr << "Test Connectivity surfs (edgeind) uni-directional Errors :" << errCountBiDir << endl;
	}

	errTot+=errCount;
	errTotBiDir+=errCountBiDir;
	errCount=0;
	errCountBiDir=0;
	kk=int(surfs.size());
	for (ii=0; ii<kk;++ii){
		testSub=volus.find_list(surfs(ii)->voluind);
		kk2=testSub.size();
		if(emptyIsErr && kk2!=2){ // voluind should be of size 2
			errCount++;
			cerr <<  " Test Connectivity Error :"  << errCount << " surf " << surfs(ii)->index
				<< " has voluind of length " <<kk2 << " list (surf::voluind): " ;
			DisplayVector(surfs(ii)->voluind); 
			cerr << endl;
		}
		for(jj=0;jj< kk2; ++jj){
			if (testSub[jj]<0 && surfs(ii)->voluind[jj]!=0){
				errCount++;
				cerr << " Test Connectivity Error :" << errCount << " surf " << surfs(ii)->index
				<< " makes unknown reference to volu " << surfs(ii)->voluind[jj] << endl;
			} else if (surfs(ii)->voluind[jj]!=0){
				ll2=volus(testSub[jj])->surfind.size();
				flag=false;
				for(ll=0;ll<ll2;ll++){
					flag=flag || (volus(testSub[jj])->surfind[ll]==surfs(ii)->index);
				}
				if (!flag){
					errCountBiDir++;
					cerr << " Test Connectivity Error :" << errCountBiDir << " surf " << surfs(ii)->index
					<< " makes uni-directional reference to volume " << surfs(ii)->voluind[jj] << " list: " ; 
					DisplayVector(surfs(ii)->voluind); cout << " list (volu.surfind): " ; 
					DisplayVector(volus(jj)->surfind);cerr << endl;
				}
			}
		}
	}
	if (errCount>0) {
		cerr << "Test Connectivity surfs (voluind) Errors :" << errCount << endl;
	}
	if (errCountBiDir>0) {
		cerr << "Test Connectivity surfs (voluind) uni-directional Errors :" << errCountBiDir << endl;
	}


	errTot+=errCount;
	errTotBiDir+=errCountBiDir;
	errCount=0;
	errCountBiDir=0;
	kk=int(volus.size());
	for (ii=0; ii<kk;++ii){
		testSub=surfs.find_list(volus(ii)->surfind);
		kk2=testSub.size();
		for(jj=0;jj< kk2; ++jj){
			if (testSub[jj]<0 && volus(ii)->surfind[jj]!=0){
				errCount++;
				cerr << " Test Connectivity Error :" << errCount << " volu " << volus(ii)->index
				<< " makes unknown reference to surface " << volus(ii)->surfind[jj] << endl;
			} else if (volus(ii)->surfind[jj]!=0){
				ll2=surfs(testSub[jj])->voluind.size();
				flag=false;
				for(ll=0;ll<ll2;ll++){
					flag=flag || (surfs(testSub[jj])->voluind[ll]==volus(ii)->index);
				}
				if (!flag){
					errCountBiDir++;
					cerr << " Test Connectivity Error :" << errCountBiDir << " volu " << volus(ii)->index
					<< " makes uni-directional reference to surface " << volus(ii)->surfind[jj] << " list: " ; 
					DisplayVector(volus(ii)->surfind); cout << " list (surfs.voluind): " ; 
					DisplayVector(surfs(testSub[jj])->voluind);cerr << endl;
				}
			}
		}
	}
	if (errCount>0){
		cerr << "Test Connectivity volus (surfind) Errors :" << errCount << endl;
	}
	if (errCountBiDir>0){
		cerr << "Test Connectivity volus (surfind) uni-directional Errors :" << errCountBiDir << endl;
	}
	errTot+=errCount;
	errTotBiDir+=errCountBiDir;
	if (errTot>0){
		cerr << errTot << "  Total errors were detected in the connectivity list "
			" in: " << endl << strRoot <<endl;
	}
	if (errTotBiDir>0){
		cerr << errTotBiDir << "  Total errors were detected in the bi-directionality of the connectivity  list in: " << endl << strRoot <<endl;
	}
	return(errTot);
}

#pragma GCC diagnostic pop
// methods for mesh
void mesh::HashArray(){
	verts.HashArray();
	edges.HashArray();
	surfs.HashArray();
	volus.HashArray();
}
void mesh::SetMaxIndex(){
	verts.SetMaxIndex();
	edges.SetMaxIndex();
	surfs.SetMaxIndex();
	volus.SetMaxIndex();
}
void mesh::SetLastIndex(){
	verts.SetLastIndex();
	edges.SetLastIndex();
	surfs.SetLastIndex();
	volus.SetLastIndex();
}

void mesh::ArraysAreHashed(){
	this->verts.isHash=1;
	this->edges.isHash=1;
	this->surfs.isHash=1;
	this->volus.isHash=1;
}

void mesh::PrepareForUse(bool needOrder){

	verts.isInMesh=true;
	edges.isInMesh=true;
	surfs.isInMesh=true;
	volus.isInMesh=true;

	if(meshDim==0){
		meshDim=3;
		if (volus.size()==0){
			meshDim--;
			if (surfs.size()==0){
				meshDim--;
				if (edges.size()==0){
					meshDim--;
				}
			}
		}
	}

	verts.PrepareForUse();
	edges.PrepareForUse();
	surfs.PrepareForUse();
	volus.PrepareForUse();
   // Additional mesh preparation steps
	// if(!meshDepIsSet){
	// 	SetMeshDepElm();
	// }

	if (!borderIsSet){
		this->SetBorders();
	}
	if(needOrder){
		this->OrderEdges();
	}
	verts.ForceArrayReady();
	edges.ForceArrayReady();
	surfs.ForceArrayReady();
	volus.ForceArrayReady();
}

void mesh::GetMaxIndex(int *nVert,int *nEdge,int *nSurf,int *nVolu) const{
	*nVert=verts.GetMaxIndex();
	*nEdge=edges.GetMaxIndex();
	*nSurf=surfs.GetMaxIndex();
	*nVolu=volus.GetMaxIndex();
}
void mesh::disp() const {
	verts.disp();
	edges.disp();
	surfs.disp();
	volus.disp();
}
void mesh::read(FILE *fid) {
	verts.read(fid);
	edges.read(fid);
	surfs.read(fid);
	volus.read(fid);

	verts.isInMesh=true;
	edges.isInMesh=true;
	surfs.isInMesh=true;
	volus.isInMesh=true;
}
void mesh::write(FILE *fid) const {
   //fprintf(fid,"%i \n",int(borderIsSet));
	verts.write(fid);
	edges.write(fid);
	surfs.write(fid);
	volus.write(fid);
}
int mesh::read(const char *str) {
   // Convenience read function taking a single char argument
	FILE *fid;
	fid=fopen(str,"r");
	if (fid!=NULL){
		this->read(fid);
		fclose(fid);
	} else {
		cout << "File " << str << "Could not be opened to read" << endl;
		return(1);
	}
	return(0);
}
int mesh::write(const char *str) const {
	FILE *fid;
	fid=fopen(str,"w");
	if (fid!=NULL){
		this->write(fid);
		fclose(fid);
	} else {
		cout << "File " << str << "Could not be opened to write" << endl;
		return(1);
	}
	return(0);
}
bool mesh::isready() const {
	bool readyforuse=true;
	readyforuse=readyforuse & verts.isready();
	readyforuse=readyforuse & edges.isready();
	readyforuse=readyforuse & surfs.isready();
	readyforuse=readyforuse & volus.isready();

	return(readyforuse);
}

void mesh::displight() const {
	cout << "mesh: vert " << verts.size();
	cout << "; edges " << edges.size();
	cout << "; surfs " << surfs.size();
	cout << "; volus " << volus.size() << endl;
}
 
void mesh::Init(int nVe,int nE, int nS, int nVo)
{
	borderIsSet=false;
	meshDim=0;
	verts.Init(nVe);
	edges.Init(nE);
	surfs.Init(nS);
	volus.Init(nVo);

	verts.isInMesh=true;
	edges.isInMesh=true;
	surfs.isInMesh=true;
	volus.isInMesh=true;

   #ifdef TEST_ARRAYSTRUCTURES
	cout << "Mesh Correctly Assigned!" << endl;
   #endif // TEST_ARRAYSTRUCTURES
}

void mesh::reserve(int nVe,int nE, int nS, int nVo)
{

	verts.reserve(nVe);
	edges.reserve(nE);
	surfs.reserve(nS);
	volus.reserve(nVo);

	verts.isInMesh=true;
	edges.isInMesh=true;
	surfs.isInMesh=true;
	volus.isInMesh=true;
}

void mesh::MakeCompatible_inplace(mesh &other) const{
   // Makes other mesh compatible with this to be 
   // merged without index crashes

	int nVert,nEdge,nSurf,nVolu;

   // Define Max indices in current mesh
	this->GetMaxIndex(&nVert,&nEdge,&nSurf,&nVolu);
	other.ChangeIndices(nVert,nEdge,nSurf,nVolu);
}

void mesh::ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu){
   /*int ii;
   // Volumes
   for(ii=0;ii<volus.size();++ii){
	  volus[ii].ChangeIndices(nVert,nEdge,nSurf,nVolu);
   }
   // Surfaces
   for(ii=0;ii<surfs.size();++ii){
	  surfs[ii].ChangeIndices(nVert,nEdge,nSurf,nVolu);
   }
   // Volumes
   for(ii=0;ii<edges.size();++ii){
	  edges[ii].ChangeIndices(nVert,nEdge,nSurf,nVolu);
   }
   // Surfaces
   for(ii=0;ii<verts.size();++ii){
	  verts[ii].ChangeIndices(nVert,nEdge,nSurf,nVolu);
   }*/
	volus.ChangeIndices(nVert,nEdge,nSurf,nVolu);
	edges.ChangeIndices(nVert,nEdge,nSurf,nVolu);
	surfs.ChangeIndices(nVert,nEdge,nSurf,nVolu);
	verts.ChangeIndices(nVert,nEdge,nSurf,nVolu);
}

mesh mesh::MakeCompatible(mesh other) const{
	MakeCompatible_inplace(other);
	return(other);
}

void mesh::Concatenate(const mesh &other){

	this->volus.Concatenate(other.volus);
	this->edges.Concatenate(other.edges);
	this->verts.Concatenate(other.verts);
	this->surfs.Concatenate(other.surfs);
}

void mesh::PopulateIndices(){

	volus.PopulateIndices();
	edges.PopulateIndices();
	verts.PopulateIndices();
	surfs.PopulateIndices();
	meshDepIsSet=false;
}


// Field Specific operations


/**
 * @brief      Orders a list of edge to be connected. This list of edges is
 *             ordered in place. THis will not work for self crossing lists.
 *
 * @param      edgeind  The edgeind
 * @param[in]  meshin   The meshin
 * @param[in]  warn     The warning
 * @param[in]  errout   The errout
 *
 * @return     the return value is a flag which can be: 
 * 				0 - the edges have been ordered and closed 
 *              1 - the edges have been ordered and closed but
 *             the list was truncated. 
 *              <0 - The edges are ordered but not closed.
 *             (need errout to be false) for edgeind of size
 *             -<return val>
 */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
int OrderEdgeList(vector<int> &edgeind, const mesh &meshin,	bool warn,
	bool errout, const vector<int>* edgeIndOrigPtr,
	const surf* surfin){
#pragma GCC diagnostic pop

	unordered_multimap<int,int> vert2Edge;
	vector<bool> isDone;
	vector<int> edgeIndOrig;
	int vertCurr, edgeCurr, ii, jj, endsExpl=0;

	int return_ordered=0, return_truncated=1, return_open=-1;

	if (edgeind.size()<=0){
		return(return_ordered);
	}
	// edgeIndOrig2=edgeind;
	if(edgeIndOrigPtr==NULL){
		// edgeInOrigPtr is not defined and there is no need
		sort(edgeind);
		unique(edgeind);
		edgeIndOrig = edgeind;
		edgeIndOrigPtr = &edgeIndOrig;
	} 

	isDone.assign(edgeIndOrigPtr->size(),false);
	auto edgeSub=meshin.edges.find_list(edgeind);
	auto edge2Vert=ConcatenateVectorField(meshin.edges,&edge::vertind,edgeSub);

	HashVector(edge2Vert, vert2Edge);

	vertCurr=edge2Vert[0];
	edgeCurr=edgeind[0];
	isDone[0]=true;
	auto it=vert2Edge.end();
	for(ii=1;ii<int(edgeind.size());++ii){
		auto range=vert2Edge.equal_range(vertCurr);
		if (vert2Edge.count(vertCurr)==1){
			if(warn || errout){
				if(errout){
					cerr << "Error";
				} else {
					cerr << "Warning";
				}
				cerr 	<< " in :" << __PRETTY_FUNCTION__ ;
				cerr << " - edgeind is not closed (single vertex "<<
					 vertCurr <<")" << endl;
				if(errout){
					RSVS3D_ERROR_ARGUMENT("edgelist is not closed");
				}
			}
			if(!endsExpl){ // explore the one way chain the other way
				endsExpl++;
				vertCurr=edge2Vert[1];
				edgeCurr=edgeind[0];
				range=vert2Edge.equal_range(vertCurr);
			} else { // Exit

				return return_open*ii;
			}
		}
	 	#ifdef SAFE_ACCESS
		if (range.first==vert2Edge.end()){
			surfin->disptree(meshin,0);

			cerr << ii << " vert " << vertCurr << "  ";
			DisplayVector(edge2Vert);
			DisplayVector(edgeind);
			meshin.verts.isearch(vertCurr)->disp();
			cout << it->second << " " << 1/2 << 2/3 <<  endl;
			cerr << "Error in :" << __PRETTY_FUNCTION__ << endl;
			RSVS3D_ERROR_RANGE("unordered_multimap went beyond its range in OrderEdges");
		}
	 	#ifdef RSVS_DIAGNOSTIC
		if (vert2Edge.count(vertCurr)!=2){
			cerr << "DIAGNOSTIC in " << __PRETTY_FUNCTION__ << endl;
			cerr << " current vertex " << vertCurr << " appears " 
				<< vert2Edge.count(vertCurr) 
				<< " in the surf::edgeind ...->vertind" << endl;
		}
		#endif //RSVS_DIAGNOSTIC
	 	#endif // SAFE_ACCESS
		jj=(*edgeIndOrigPtr)[(range.first->second)/2]==edgeCurr;

		it=range.first;
		if (jj){++it;}
		if (!isDone[(it->second)/2]){
			isDone[(it->second)/2]=true;
			edgeCurr=(*edgeIndOrigPtr)[(it->second)/2];
			// Warning ((x/2)*2) not necessarily equal x
			// Done to round down to the nearest even
			jj=edge2Vert[((it->second)/2)*2]==vertCurr; 
			vertCurr=edge2Vert[((it->second)/2)*2+jj];
			edgeind[ii]=edgeCurr;
		} else if (endsExpl) {
			return return_open*ii;
		} else {
		 	#ifdef RSVS_DIAGNOSTIC_FIXED
			cerr << "DIAGNOSTIC in " << __PRETTY_FUNCTION__ 
				<< " surface about to be truncated, if this should not be the case"
				" this can be due to duplicates in the edgeind list" << endl;
			cerr << "edgeind : "; DisplayVector(edgeind);
			cerr << endl << "edgeIndOrig : "; DisplayVector((*edgeIndOrigPtr));
			cerr << endl << "Original connectivity setup:" << endl;
			auto temp = edgeind;
			edgeind=(*edgeIndOrigPtr);
			surfin->disptree(meshin,1);
			edgeind = temp;
			#endif //RSVS_DIAGNOSTIC

			edgeind.erase(edgeind.begin()+ii, edgeind.end());

		 	#ifdef RSVS_DIAGNOSTIC_FIXED
			cerr << endl << "New connectivity setup:" << endl;
			surfin->disptree(meshin,1);
			#endif //RSVS_DIAGNOSTIC

			return return_truncated;
			break; // Unnecessary
		}
	}
	return return_ordered;
}


/**
 * @brief      { function_description }
 *
 * @param      meshin  The meshin
 *
 * @return     { description_of_the_return_value }
 */
/*
int surf::OrderEdges(mesh *meshin)
{
	unordered_multimap<int,int> vert2Edge;
	vector<int> edge2Vert,edgeSub,edgeIndOrig;
	vector<bool> isDone;
	bool isTruncated;
	int vertCurr,edgeCurr;
	int ii,jj, newSurfInd, kk, prevHash;
	// std::pair<unordered_multimap<int,int>::iterator,unordered_multimap<int,int>::iterator> range;
	unordered_multimap<int,int>::iterator it;
	isTruncated=false;
	newSurfInd=-1; 
	if (edgeind.size()>0){
		// edgeIndOrig2=edgeind;
		sort(edgeind);
		unique(edgeind);

		edgeIndOrig=edgeind;
		isDone.assign(edgeIndOrig.size(),false);
		edgeSub=meshin->edges.find_list(edgeind);
		edge2Vert=ConcatenateVectorField(meshin->edges,&edge::vertind,edgeSub);

		HashVector(edge2Vert, vert2Edge);

		vertCurr=edge2Vert[0];
		edgeCurr=edgeind[0];
		isDone[0]=true;
		it=vert2Edge.end();
		for(ii=1;ii<int(edgeind.size());++ii){
			auto range=vert2Edge.equal_range(vertCurr);
		 	#ifdef SAFE_ACCESS
			if (range.first==vert2Edge.end()){
				disptree((*meshin),0);

				cerr << ii << " vert " << vertCurr << "  ";
				DisplayVector(edge2Vert);
				DisplayVector(edgeind);
				meshin->verts.isearch(vertCurr)->disp();
				cout << it->second << " " << 1/2 << 2/3 <<  endl;
				cerr << "Error in :" << __PRETTY_FUNCTION__ << endl;
				RSVS3D_ERROR_RANGE("unordered_multimap went beyond its range in OrderEdges");
			}
			if (vert2Edge.count(vertCurr)==1){
				// cerr << "ERROR : Surface does not form closed loop" << endl;
				// disptree((*meshin),0);

				// jj=edgeIndOrig[(range.first->second)/2]==edgeCurr;
				// DisplayVector(edge2Vert);
				// DisplayVector(edgeind);
				// cerr <<endl;

				// meshin->verts.isearch(vertCurr)->disp();
				// cerr << "ii is : " << ii << " jj is : " << jj << " count is : " ;
				// jj=vert2Edge.count(vertCurr);
				// cerr << jj  <<  endl; 
				cerr << "Error in :" << __PRETTY_FUNCTION__ ;
				//RSVS3D_ERROR_RANGE("Found a single vertex - surface is not closed");
				cerr << " - surface "<< index <<" is not closed (single vertex "<<
					 vertCurr <<")" << endl;
				this->edgeind = edgeIndOrig;
			 	this->disp();
			 	#ifdef RSVS_DIAGNOSTIC
			 	isordered=true; // uncomment these two to allow run through for diagnostic
				return(newSurfInd);
				#endif //RSVS_DIAGNOSTIC
				RSVS3D_ERROR_ARGUMENT("Surface is not closed");
			}
		 	#ifdef RSVS_DIAGNOSTIC
			if (vert2Edge.count(vertCurr)!=2){
				cerr << "DIAGNOSTIC in " << __PRETTY_FUNCTION__ << endl;
				cerr << " current vertex " << vertCurr << " appears " 
					<< vert2Edge.count(vertCurr) 
					<< " in the surf::edgeind ...->vertind" << endl;
			}
			#endif //RSVS_DIAGNOSTIC
		 	#endif // SAFE_ACCESS
			jj=edgeIndOrig[(range.first->second)/2]==edgeCurr;

			it=range.first;
			if (jj){++it;}
			if (!isDone[(it->second)/2]){
				isDone[(it->second)/2]=true;
				edgeCurr=edgeIndOrig[(it->second)/2];
				// Warning ((x/2)*2) not necessarily equal x
				// Done to round down to the nearest even
				jj=edge2Vert[((it->second)/2)*2]==vertCurr; 
				vertCurr=edge2Vert[((it->second)/2)*2+jj];
				edgeind[ii]=edgeCurr;
			} else {
			 	#ifdef RSVS_DIAGNOSTIC_FIXED
				cerr << "DIAGNOSTIC in " << __PRETTY_FUNCTION__ 
					<< " surface about to be truncated, if this should not be the case"
					" this can be due to duplicates in the edgeind list" << endl;
				cerr << "edgeind : "; DisplayVector(edgeind);
				cerr << endl << "edgeIndOrig : "; DisplayVector(edgeIndOrig);
				cerr << endl << "Original connectivity setup:" << endl;
				auto temp = edgeind;
				edgeind=edgeIndOrig;
				disptree((*meshin),1);
				edgeind = temp;
				#endif //RSVS_DIAGNOSTIC

				edgeind.erase(edgeind.begin()+ii, edgeind.end());

			 	#ifdef RSVS_DIAGNOSTIC_FIXED
				cerr << endl << "New connectivity setup:" << endl;
				disptree((*meshin),1);
				#endif //RSVS_DIAGNOSTIC

				isTruncated=true;
				if ((meshin->surfs.capacity()-meshin->surfs.size())==0){
					// Checks that there is space to add a new face if truncation
					// is in order
			 		#ifdef RSVS_DIAGNOSTIC
			 		cerr << "DIAGNOSTIC in " << __PRETTY_FUNCTION__ 
			 			<< " reallocation necessary."
			 			<< endl; 
					#endif //RSVS_DIAGNOSTIC
					this->edgeind = edgeIndOrig;
					this->isordered = false;
					meshin->surfs.reserve(meshin->surfs.capacity() 
						+ meshin->surfs.capacity() / 2);
					return(newSurfInd);
				}
				break;
			}
		}
		if (isTruncated){ 
			// This adds a second surface if the surface closes early
			vector<int> newSurfedgeind;
			surf newSurf=*this;
			meshin->surfs.SetMaxIndex();
			newSurf.index=meshin->surfs.GetMaxIndex()+1;
			// newSurf.voluind=voluind;
			// gets added to newSurf.edgeind in SwitchIndex
			newSurf.edgeind.clear();
			newSurfedgeind.clear(); 
			for (ii=0;ii<int(edgeIndOrig.size());++ii){
				if(!isDone[ii]){
					newSurfedgeind.push_back(edgeIndOrig[ii]);
				}
			}

			// triggers a reallocation invalidating the this pointer			
			meshin->surfs.push_back(newSurf); 
			// meshin->surfs(meshin->surfs.size()-1)->disp();
			prevHash = meshin->surfs.isHash;
			meshin->surfs.isHash=1;
			meshin->SwitchIndex(6,this->index,newSurf.index,newSurfedgeind);
			meshin->surfs[meshin->surfs.size()-1].OrderEdges(meshin);

			meshin->surfs.isHash=prevHash;
			prevHash = meshin->volus.isHash;
		 	#ifdef SAFE_ACCESS
			if (this->voluind.size()<2){
				cerr << "Error in " << __PRETTY_FUNCTION__ << endl;
				cerr << " voluind is of size " << this->voluind.size() << endl;
				RSVS3D_ERROR_LOGIC("surf::voluind should be size 2");
			}
		 	#endif // SAFE_ACCESS
			for (int i = 0; i < 2; ++i)
			{
				kk=meshin->volus.find(voluind[i]);
				if(kk!=-1){
					meshin->volus[kk].surfind.push_back(newSurf.index);
				}
				meshin->volus.isHash=prevHash;
			}
			
			newSurfInd = newSurf.index;
		}
		isordered=true;
	}
	return(newSurfInd);
}*/

int surf::OrderEdges(mesh *meshin)
{
	unordered_multimap<int,int> vert2Edge;
	vector<int> edgeIndOrig;
	vector<bool> isDone;
	bool isTruncated;
	int  newSurfInd;
	// std::pair<unordered_multimap<int,int>::iterator,unordered_multimap<int,int>::iterator> range;
	unordered_multimap<int,int>::iterator it;
	isTruncated=false;
	newSurfInd=-1; 
	// do nothing if edgeind is empty
	if(meshin==NULL){
		this->isordered=false;
		return(newSurfInd);
	}
	if (edgeind.size()<=0){
		this->isordered=true;
		return(newSurfInd);
	}
	// edgeIndOrig2=edgeind;
	sort(this->edgeind);
	unique(this->edgeind);
	edgeIndOrig = this->edgeind;

	int retFlag = OrderEdgeList(this->edgeind, *meshin,true, true, 
		&edgeIndOrig, this);

	isTruncated = retFlag==1;
	
	if(!isTruncated){
		newSurfInd = -1; // No new surface index
		this->isordered=true;
		return(newSurfInd);
	}

	if ((meshin->surfs.capacity()-meshin->surfs.size())==0){
		// Checks that there is space to add a new face if surface truncation
		// is in order. If there is no space add space and exit
 		#ifdef RSVS_DIAGNOSTIC_FIXED
 		cerr << "DIAGNOSTIC in " << __PRETTY_FUNCTION__ 
 			<< " reallocation necessary."
 			<< endl; 
		#endif //RSVS_DIAGNOSTIC_FIXED
		this->edgeind = edgeIndOrig;
		this->isordered = false;
		meshin->surfs.reserve(meshin->surfs.capacity() 
			+ meshin->surfs.capacity() / 2);
		newSurfInd = -1; // No new surface index
		return(newSurfInd);
	}

	newSurfInd=this->SplitSurface(*meshin, edgeIndOrig);

	this->isordered=true;
	return(newSurfInd);
}

int surf::SplitSurface(mesh &meshin, const vector<int> &fullEdgeInd){

	int newSurfInd , kk;
	// This adds a second surface if the surface closes early
	vector<int> newSurfedgeind;
	surf newSurf=*this;
	meshin.surfs.SetMaxIndex();
	newSurf.index=meshin.surfs.GetMaxIndex()+1;
	// newSurf.voluind=voluind;
	// gets added to newSurf.edgeind in SwitchIndex
	newSurf.edgeind.clear();

	// build new surface edgeind
	for(auto allInd : fullEdgeInd){
		bool flag = false; 
		for(auto currInd : this->edgeind){
			flag = currInd==allInd;
			if(flag){break;}
		}
		if(!flag){newSurf.edgeind.push_back(allInd);}
	}
	#ifdef SAFE_ALGO
	if ((newSurf.edgeind.size()!=(fullEdgeInd.size()-this->edgeind.size())) ||
		newSurf.edgeind.size()==0 || this->edgeind.size()==0){
		cerr << "Error in:  " << __PRETTY_FUNCTION__  << endl;
		cerr << "  New surface edgeind is not of the correct size:" << endl
			<< " newSurf.edgeind (" << newSurf.edgeind.size()
			<< ") !=  fullEdgeInd (" << fullEdgeInd.size()
			<< ")-this::edgeind(" << this->edgeind.size()
			<< ")" << endl;
		RSVS3D_ERROR_ARGUMENT("Size of edgeind of truncated surface is incorrect.");
	}
	#endif
	#ifdef RSVS_DIAGNOSTIC_FIXED
	DisplayVector(this->edgeind);
	DisplayVector(newSurf.edgeind);
	#endif // RSVS_DIAGNOSTIC_FIXED
	// triggers a reallocation invalidating the this pointer			
	meshin.surfs.push_back(newSurf); 
	// meshin.surfs(meshin.surfs.size()-1)->disp();
	auto prevHash = meshin.surfs.isHash;
	meshin.surfs.isHash=1;
	meshin.SwitchIndex(6,this->index,newSurf.index,newSurf.edgeind);
	#ifdef SAFE_ALGO
	if(this->edgeind.size()==0){
		RSVS3D_ERROR_ARGUMENT("surf::edgeind was deleted in "
			"surf::SplitSurface process");
	}
	#endif //SAFE_ALGO
	meshin.surfs[meshin.surfs.size()-1].OrderEdges(&meshin);

	meshin.surfs.isHash=prevHash;
	prevHash = meshin.volus.isHash;
 	#ifdef SAFE_ACCESS
	if (this->voluind.size()<2){
		cerr << "Error in " << __PRETTY_FUNCTION__ << endl;
		cerr << " voluind is of size " << this->voluind.size() << endl;
		RSVS3D_ERROR_LOGIC("surf::voluind should be size 2");
	}
 	#endif // SAFE_ACCESS
	for (int i = 0; i < 2; ++i)
	{
		kk=meshin.volus.find(voluind[i]);
		if(kk!=-1){
			meshin.volus[kk].surfind.push_back(newSurf.index);
		}
		meshin.volus.isHash=prevHash;
	}
	
	newSurfInd = newSurf.index;
	#ifdef RSVS_DIAGNOSTIC_FIXED
	cerr << " Succesful surface split " << this->index << " | "
		<< newSurfInd << endl;
	#endif
	return(newSurfInd);
}

void surf::OrderedVerts(const mesh *meshin, vector<int> &vertList) const{
	int jj,n,actVert,edgeCurr;
	int verts[2],vertsPast[2];
	
	
	n=int(edgeind.size());
	vertList.clear();
	vertList.reserve(n);
	//This comes out in the same order as edgeind
	edgeCurr=meshin->edges.find(edgeind[n-1]);
	verts[0]=(meshin->edges(edgeCurr)->vertind[0]);
	verts[1]=(meshin->edges(edgeCurr)->vertind[1]);
	vertsPast[0]=verts[0];
	vertsPast[1]=verts[1];

	for(jj=0;jj<int(n);++jj){
		edgeCurr=meshin->edges.find(edgeind[jj]);

		verts[0]=(meshin->edges(edgeCurr)->vertind[0]);
		verts[1]=(meshin->edges(edgeCurr)->vertind[1]);

		if ((verts[0]==vertsPast[0]) || (verts[1]==vertsPast[0])){
			actVert=0;
		} 
		#ifdef SAFE_ALGO
		else if ((verts[0]==vertsPast[1]) || (verts[1]==vertsPast[1])) {
			actVert=1;
		}
		#endif //TEST_POSTPROCESSING
		else {
			actVert=1;
			#ifdef SAFE_ALGO
			cerr << "Error: Surface is not ordered " << endl;
			cerr << "	in " << __PRETTY_FUNCTION__ << endl;
			RSVS3D_ERROR_ARGUMENT("Surface not ordered");
			#endif
		}

		vertList.push_back(vertsPast[actVert]);
		vertsPast[0]=verts[0];
		vertsPast[1]=verts[1];
	}
} 

void surf::FlipVolus(){
	int interm;
	interm=voluind[0];
	voluind[0]=voluind[1];
	voluind[1]=interm;
}

bool edge::vertconneq(const edge &other) const{

	return 
		(this->vertind[0]==other.vertind[0] 
			&& this->vertind[1]==other.vertind[1])
		|| (this->vertind[1]==other.vertind[0] 
			&& this->vertind[0]==other.vertind[1]);
}

bool surf::edgeconneq(const surf &other, bool recurse) const{
	/*
	Explores all the edgeind combinations looking for a one to one 
	match with each.

	Boolean recurse controls the recursive behaviour of the function.
	The default is true so that a user calling it will get robust
	equality check between the connectivity lists. THis is then called by
	itself with recurse=false to get the other direction.

	This function is very inneficient at O(2*m*n) but exact, constant, and 
	requires no memory allocation. 
	*/
	bool out, temp;
	int count=this->edgeind.size();
	int count2=other.edgeind.size();
	out =true;
	for (int i = 0; i < count; ++i)
	{
		temp = false;
		for (int j = 0; j < count2; ++j)
		{
			temp = temp || (this->edgeind[i]==other.edgeind[j]);
		}
		out = out && temp;
		if(!out){
			break;
		}
	}
	if (recurse){
		return out && other.edgeconneq(*this, false);
	} else {
		return out;
	}
}

vector<int> mesh::OrderEdges(){
	int ii, kk=0;
	int newInd;
	vector<int> newSurfPairs;
	// Use the address of the first element to ensure that no reallocation
	if(surfs.size()>0){
		auto origSurf = surfs(0);
		auto capPrev = surfs.capacity();
		bool pntrFlag, capFlag;

		do{
			kk++;
			origSurf = surfs(0);
			capPrev = surfs.capacity();
			for (ii = 0; ii < surfs.size(); ++ii)
			{
				if (!surfs(ii)->isready(true)){
					newInd=surfs.elems[ii].OrderEdges(this);
					if(newInd>-1){
						newSurfPairs.push_back(surfs(ii)->index);
						newSurfPairs.push_back(newInd);
					}
				}
			}
		// while the vector has been reallocated.
			pntrFlag =(origSurf != surfs(0)); // test 1st pntr
			capFlag = (capPrev != surfs.capacity());
		} while(pntrFlag || capFlag) ;
		#ifdef SAFE_ALGO
		for (ii = 0; ii < surfs.size(); ++ii)
		{
			if (!surfs(ii)->isready(true)){
				RSVS3D_ERROR_LOGIC("After ordering edges a surface is not ordered.");
			}
		}
		#endif //SAFE_ALGO
	}
	return(newSurfPairs);
}

void mesh::GetOffBorderVert(vector<int> &vertInd,  vector<int> &voluInd,
	int outerVolume){
	/*
	Gets vertices that are in a vlume that is on the edge of the design
	space but off th eedge themselves

	outerVolume indicates the additional condition of the volume needing 
	to be full or empty in the parents.
	-1 return all vertices
	0 return vertices where the target volume is 1.0
	1 return vertices where the target volume is >0.0
	*/

	if(!borderIsSet){
		this->SetBorders();
	}
	this->GetOffBorderVert(vertInd, voluInd, outerVolume);
}

void mesh::GetOffBorderVert(vector<int> &vertInd, vector<int> &voluInd,
	int outerVolume) const {

	if (this->WhatDim()==3){
		this->GetOffBorderVert3D(vertInd, voluInd, outerVolume);
	} else if (this->WhatDim()==2) {
		this->GetOffBorderVert2D(vertInd, voluInd, outerVolume);
	}
}

void mesh::GetOffBorderVert3D(vector<int> &vertInd, vector<int> &voluInd,
	int outerVolume) const{
	/*
	Gets vertices that are in a vlume that is on the edge of the design
	space but off th eedge themselves

	outerVolume indicates the additional condition of the volume needing 
	to be full or empty in the parents.
	-1 return all vertices
	0 return vertices where the target volume is 1.0
	1 return vertices where the target volume is >0.0
	*/
	int ii, ni, jj, nj;
	vector<double> vals;
	bool voluOnBoundary;

	ni=volus.size();
	nj=verts.size();
	vertInd.clear();
	vertInd.reserve(nj);
	voluInd.clear();
	voluInd.reserve(ni);

	for (ii=0; ii<ni; ++ii){
		voluOnBoundary=volus(ii)->isBorder;
		if(voluOnBoundary){
			if(outerVolume==0){
				this->VoluValuesofParents(volus(ii)->index,vals, &volu::target);
				voluOnBoundary=false;
				jj=0;
				while (!voluOnBoundary && jj<meshtree.nParents){
					voluOnBoundary=vals[jj]==1.0;
					++jj;
				}
			} else if(outerVolume==1){
				this->VoluValuesofParents(volus(ii)->index,vals, &volu::target);
				voluOnBoundary=false;
				jj=0;
				while (!voluOnBoundary && jj<meshtree.nParents){
					voluOnBoundary=vals[jj]>0.0;
					++jj;
				}
			} else if(outerVolume!=-1){
				RSVS3D_ERROR_ARGUMENT("Unkownn value of outerVolume -1,0, or 1");
			}
		} 
		// THIS IS DODGY HOW to pick the side to delete can fail
		if(voluOnBoundary){
			// Pick one of the volumes to delete
			if(outerVolume==0){
				for(auto surfAct : volus(ii)->surfind){
					for(auto voluAct : surfs.isearch(surfAct)->voluind){
						if(voluAct!=0 && !(volus.isearch(voluAct)->isBorder)){
							voluInd.push_back(voluAct);
						}
					}
				}
			} else {
				voluInd.push_back(volus(ii)->index);
			}
			for(auto surfAct : volus(ii)->surfind){
				for(auto edgeAct : surfs.isearch(surfAct)->edgeind){
					for (auto vertAct : edges.isearch(edgeAct)->vertind){
						if (!verts.isearch(vertAct)->isBorder){
							vertInd.push_back(vertAct);
						}
					}
				}
			}
		}
	}
}

void mesh::GetOffBorderVert2D(vector<int> &vertInd, vector<int> &surfind,
	int outerVolume) const{
	/*
	Gets vertices that are in a vlume that is on the edge of the design
	space but off th eedge themselves

	outerVolume indicates the additional condition of the volume needing 
	to be full or empty in the parents.
	-1 return all vertices
	0 return vertices where the target volume is 1.0
	1 return vertices where the target volume is >0.0
	*/
	int ii, ni, jj, nj,kk,nk, vertTemp, edgeSub,surfSub;
	vector<double> vals;
	bool surfCond;

	ni=surfs.size();
	nj=verts.size();
	vertInd.clear();
	vertInd.reserve(nj);
	surfind.clear();
	surfind.reserve(ni);

	for (ii=0; ii<ni; ++ii){
		surfCond=surfs(ii)->isBorder;
		if(surfCond){
			if(outerVolume==0){
				this->SurfValuesofParents(surfs(ii)->index,vals, &surf::target);
				surfCond=false;
				jj=0;
				while (!surfCond && jj<meshtree.nParents){
					surfCond=vals[jj]==1.0;
					++jj;
				}
			} else if(outerVolume==1){
				this->SurfValuesofParents(surfs(ii)->index,vals, &surf::target);
				surfCond=false;
				jj=0;
				while (!surfCond && jj<meshtree.nParents){
					surfCond=vals[jj]>0.0;
					++jj;
				}
			} else if(outerVolume!=-1){
				RSVS3D_ERROR_ARGUMENT("Unkownn value of outerVolume -1,0, or 1");
			}
		} 
		// THIS IS DODGY HOW to pick the side to delete can fail
		if(surfCond){
			// Pick one of the surfaces to delete
			if(outerVolume==0){
				nj = surfs(ii)->edgeind.size();
				for(jj=0;jj<nj;++jj){
					surfSub = edges.find(surfs(ii)->edgeind[jj]);
					for(kk=0;kk<2;++kk){
						if(edges(surfSub)->surfind[kk]!=0){
							if(!(surfs.isearch(edges(surfSub)->surfind[kk])->isBorder)){
								surfind.push_back(edges(surfSub)->surfind[kk]);
							}
						}
					}
				}
			} else {
				surfind.push_back(surfs(ii)->index);
			}
			surfSub=ii;
			nj=surfs(surfSub)->edgeind.size();
			for(jj=0; jj<nj; ++jj){
				edgeSub=edges.find(surfs(surfSub)->edgeind[jj]);
				nk=edges(edgeSub)->vertind.size();
				for (kk=0; kk<nk; ++kk){
					vertTemp=edges(edgeSub)->vertind[kk];
					if (!verts.isearch(vertTemp)->isBorder){
						vertInd.push_back(vertTemp);
					}
				}
			}

		}
	}
}

void mesh::SetBorders(){
	int ii,jj,nT;
	if (int(volus.size())>0){
	  // Update border status of edges
		for(ii=0;ii<surfs.size();++ii){
			jj=0;
			nT=surfs(ii)->voluind.size();
			surfs[ii].isBorder=false;
			while(jj<nT && !surfs(ii)->isBorder){
				surfs[ii].isBorder=surfs[ii].voluind[jj]==0;
				jj++;
			}
		}
		surfs.ForceArrayReady();

	  // Update border status of volus
		for(ii=0;ii<volus.size();++ii){
			jj=0;
			nT=volus(ii)->surfind.size();
			volus[ii].isBorder=false;
			while(jj<nT && !volus(ii)->isBorder){
				volus[ii].isBorder=surfs(surfs.find(volus[ii].surfind[jj]))->isBorder;
				jj++;
			}
		}
		volus.ForceArrayReady();

	  // Update border status of edges  
		for(ii=0;ii<edges.size();++ii){
			jj=0;
			nT=edges(ii)->surfind.size();
			edges[ii].isBorder=false;
			while(jj<nT && !edges(ii)->isBorder){
				edges[ii].isBorder=surfs(surfs.find(edges[ii].surfind[jj]))->isBorder;
				jj++;
			}
		}
		edges.ForceArrayReady();

	} else {
		for(ii=0;ii<edges.size();++ii){
			jj=0;
			nT=edges(ii)->surfind.size();
			edges[ii].isBorder=false;
			while(jj<nT && !edges(ii)->isBorder){
				edges[ii].isBorder=(edges[ii].surfind[jj]==0);
				jj++;
			}
		}
		edges.ForceArrayReady();

		for(ii=0;ii<surfs.size();++ii){
			jj=0;
			nT=surfs(ii)->edgeind.size();
			surfs[ii].isBorder=false;
			while(jj<nT && !surfs(ii)->isBorder){
				surfs[ii].isBorder=edges(edges.find(surfs[ii].edgeind[jj]))->isBorder;
				jj++;
			}
		}
		surfs.ForceArrayReady();
	}


   // Update border status of edges
	for(ii=0;ii<verts.size();++ii){
		jj=0;
		nT=verts(ii)->edgeind.size();
		verts[ii].isBorder=false;
		while(jj<nT && !verts(ii)->isBorder){
			verts[ii].isBorder=edges(edges.find(verts[ii].edgeind[jj]))->isBorder;
			jj++;
		}
	}
	verts.ForceArrayReady();

	borderIsSet=true;
}

void mesh::ForceCloseContainers(){

	int ii,jj,iEdge,iSurf,kk;
	int nVert,nEdge,nSurf,nBlocks;
	bool is3DMesh=volus.size()>0;
	vector<int> vertBlock;
	vector<bool> surfsIsDone;

	vertBlock.clear();
	nBlocks=this->ConnectedVertex(vertBlock);
	surfsIsDone.assign(surfs.size(),false);
	nVert=verts.size();
	if (is3DMesh){
	  // reassign volumes
		volus.elems.clear();
		volus.Init(nBlocks);
		volus.PopulateIndices();
		volus.HashArray();
		for(ii=0;ii<nVert;ii++){
			nEdge=verts(ii)->edgeind.size();
			for(jj=0;jj<nEdge;++jj){
				iEdge=edges.find(verts(ii)->edgeind[jj]);
				#ifdef SAFE_ALGO
				if(vertBlock[verts.find(edges(iEdge)->vertind[0])] 
					!= vertBlock[verts.find(edges(iEdge)->vertind[1])]){
					cerr << endl <<"An edge has connections in 2 vertBlocks:" << endl;
							cerr << __PRETTY_FUNCTION__ << endl;
				}
				#endif //SAFE_ALGO
				nSurf=edges(iEdge)->surfind.size();
				for (kk=0;kk<nSurf;++kk){
					iSurf=surfs.find(edges(iEdge)->surfind[kk]);
					if (surfsIsDone[iSurf]){
						#ifdef SAFE_ALGO
						if(surfs.elems[iSurf].voluind[0]!=vertBlock[ii]){
							cerr << endl <<"Surf voluind assignement performed twice during:" << endl;
							cerr << __PRETTY_FUNCTION__ << endl;
						}
						#endif // SAFE_ALGO
					}
					volus.elems[volus.find(vertBlock[ii])].surfind.push_back(edges(iEdge)->surfind[kk]);
					surfs.elems[iSurf].voluind.clear();
					surfs.elems[iSurf].voluind.push_back(vertBlock[ii]);
					surfs.elems[iSurf].voluind.push_back(0);
					surfsIsDone[iSurf] = true;
				}
			}
		}
	} else {
	  // reassign surfaces
		surfs.elems.clear();
		surfs.Init(nBlocks);
		surfs.PopulateIndices();
		surfs.HashArray();

		for(ii=0;ii<nVert;ii++){
			nEdge=verts(ii)->edgeind.size();
			for(jj=0;jj<nEdge;++jj){
				iEdge=edges.find(verts(ii)->edgeind[jj]);
				surfs.elems[surfs.find(vertBlock[ii])].edgeind.push_back(verts(ii)->edgeind[jj]);
				edges.elems[iEdge].surfind.clear();
				edges.elems[iEdge].surfind.push_back(vertBlock[ii]);
				edges.elems[iEdge].surfind.push_back(0);
			}
		}
	}



	verts.ForceArrayReady();
	surfs.ForceArrayReady();
	edges.ForceArrayReady();
	volus.ForceArrayReady();
}

std::vector<int> mesh::MergeGroupedVertices(HashedVector<int, int> &closeVert,
	bool delVerts){
	/*
	Uses a hashed vector to merge points grouped together.

	Deletes points which have been merged out.
	*/
	std::vector<bool> isSnaxDone;
	std::vector<int> sameSnaxs, rmvInds;
	int nSnax, countJ, countRep;

	nSnax = closeVert.vec.size();
	if (!closeVert.isHash){
		closeVert.GenerateHash();
	}
	rmvInds.reserve(nSnax);
	isSnaxDone.assign(nSnax,false);

	// Find matching elements and perform the replacement process
	// DisplayVector(closeVert.vec);
	for (int i = 0; i < nSnax; ++i)
	{
		if((!isSnaxDone[i]) && (closeVert.vec[i]!=-1))
		{
			sameSnaxs.clear();
			sameSnaxs=closeVert.findall(closeVert.vec[i]);
			countJ = sameSnaxs.size();
			for (int j = 1; j < countJ; ++j)
			{
				this->SwitchIndex(1, this->verts(sameSnaxs[j])->index, 
					this->verts(sameSnaxs[0])->index);
				countRep++;
				rmvInds.push_back(this->verts(sameSnaxs[j])->index);
				isSnaxDone[sameSnaxs[j]]=true;
			}
			isSnaxDone[sameSnaxs[0]]=true;

		}  else{
			isSnaxDone[i]=true;
		}
	}

	if(rmvInds.size()>0 && delVerts){
		sort(rmvInds);
		unique(rmvInds);
		this->verts.remove(rmvInds);
		this->verts.PrepareForUse();
	}
	return(rmvInds);
}

void mesh::RemoveSingularConnectors(const std::vector<int> &rmvVertInds, bool voidError){
	/*
	Removes degenerate connectors.

	This removes edges, surfaces and volumes which do not have enough
	children (ie connector elements of lower dimensionality) to be a closed container.
	Edges < 2 vertices
	Surfaces < 3 edges
	Volumes < 4 Surfaces

	Offers the option to delete vertices as well.
	*/

	std::vector<int> rmvInds,rmvInds2;
	int nEdgeSurf;
	int count, countRep;

	if(rmvVertInds.size()>0){
		rmvInds=rmvVertInds;

		sort(rmvInds);
		unique(rmvInds);
		this->verts.remove(rmvInds);
		this->verts.PrepareForUse();
		#ifdef RSVS_VERBOSE
		std::cout << "Number of removed vertices " << rmvInds.size() << std::endl;
		#endif //RSVS_VERBOSE
		rmvInds.clear();
	}
	
	// Remove Edges
	count  = this->edges.size();
	countRep= 0;
	for (int i = 0; i < count; ++i)
	{
		if(this->edges(i)->vertind[0]==this->edges(i)->vertind[1]){
			this->RemoveIndex(2, this->edges(i)->index);
			rmvInds.push_back(this->edges(i)->index);
			countRep++;
		}
	}
	#ifdef RSVS_VERBOSE
	std::cout << "Number of removed edges " << countRep << std::endl;
	#endif //RSVS_VERBOSE
	sort(rmvInds);
	unique(rmvInds);
	this->edges.remove(rmvInds);
	this->edges.PrepareForUse();
	rmvInds.clear();

	// Remove Surfs
	count  = this->surfs.size();
	countRep= 0;
	rmvInds2.clear();
	for (int i = 0; i < count; ++i)
	{
		nEdgeSurf = this->surfs(i)->edgeind.size();
		if (nEdgeSurf==2){
			bool edgeConnEq = this->edges.isearch(this->surfs(i)->edgeind[0])
				->vertconneq(*(this->edges.isearch(this->surfs(i)->edgeind[1])));
			if (edgeConnEq && (this->surfs(i)->edgeind[1]!=this->surfs(i)->edgeind[0])){
				rmvInds2.push_back(this->surfs(i)->edgeind[1]);
				this->SwitchIndex(2, this->surfs(i)->edgeind[1], 
					this->surfs(i)->edgeind[0]);
			}
		}
		if(nEdgeSurf<3){
			this->RemoveIndex(3, this->surfs(i)->index);
			rmvInds.push_back(this->surfs(i)->index);
			countRep++;
		}
	}
	#ifdef RSVS_VERBOSE
	std::cout << "Number of removed surfs " << countRep << std::endl;
	#endif //RSVS_VERBOSE
	sort(rmvInds);
	unique(rmvInds);
	this->surfs.remove(rmvInds);
	this->surfs.PrepareForUse();
	rmvInds.clear();
	sort(rmvInds2);
	unique(rmvInds2);
	#ifdef RSVS_VERBOSE
	std::cout << "Number of removed edges (duplicates) " << rmvInds2.size() << std::endl;
	#endif //RSVS_VERBOSE
	this->edges.remove(rmvInds2);
	this->edges.PrepareForUse();
	rmvInds2.clear();

	// Remove Volus
	count  = this->volus.size();
	countRep= 0;
	// std::cout << this->surfs.isHash <<  "  ";
	for (int i = 0; i < count; ++i)
	{
		nEdgeSurf = this->volus(i)->surfind.size();
		if(nEdgeSurf==2){
			// std::cout << i << " " << this->surfs.isHash <<  " | ";
			bool surfConnEq = this->surfs.isearch(this->volus(i)->surfind[0])
				->edgeconneq(*(this->surfs.isearch(this->volus(i)->surfind[1])));
			// if the surfaces have the same edge connectivities and different indices
			if(surfConnEq && (this->volus(i)->surfind[0]!=this->volus(i)->surfind[1])){
				rmvInds2.push_back(this->volus(i)->surfind[1]);
				this->SwitchIndex(3, this->volus(i)->surfind[1], 
					this->volus(i)->surfind[0]);
				// std::cout << i << " " << this->surfs.isHash <<  " | ";
			}
		} else if (nEdgeSurf==3 && voidError){
			// This case is unhandled as either there is still a degenerate face
			// which should not be the case or a void is being created which should not
			// be the case either.
			RSVS3D_ERROR_LOGIC("Unhandled case - Degenerate face or "
				"unwanted Void creation");
		}
		if(nEdgeSurf<4){
			this->RemoveIndex(4, this->volus(i)->index);
			rmvInds.push_back(this->volus(i)->index);
			countRep++;
		}
	}
	#ifdef RSVS_VERBOSE
	std::cout << "Number of removed volus " << countRep << std::endl;
	#endif //RSVS_VERBOSE
	sort(rmvInds);
	unique(rmvInds);
	this->volus.remove(rmvInds);
	this->volus.PrepareForUse();
	rmvInds.clear();
	sort(rmvInds2);
	unique(rmvInds2);
	#ifdef RSVS_VERBOSE
	std::cout << "Number of removed surfs (duplicates) " << rmvInds2.size() << std::endl;
	#endif //RSVS_VERBOSE
	this->surfs.remove(rmvInds2);
	this->surfs.PrepareForUse();
	rmvInds2.clear();

}

int mesh::ConnectedVertex(vector<int> &vertBlock) const{
   // Fills a vector with a number for each vertex corresponding to a
   // group of connected edges it is part of , can be used close surfaces in 2D or volumes
   // in 3D.
   // Uses a flood fill with queue method


	int nVertExplored,nVerts,nBlocks,nCurr,nEdgesCurr,ii,jj,kk;
	vector<bool> vertStatus; // 1 explored 0 not explored

	vector<int> currQueue, nextQueue; // Current and next queues of indices

	// Preparation of the arrays;
	nVerts=verts.size();
	nBlocks=0;
	nVertExplored=0;

	vertStatus.assign(nVerts,false);

	if(int(vertBlock.size())==0){
		// Standard behaviour of grouping all connected vertices when vertBlock
		// is empty
		vertBlock.clear();
		vertBlock.assign(nVerts,0);
	} else if(int(vertBlock.size())==nVerts){
		// Treat vertices with non zero value as already processed.
		// This allows these vertices to separate blocks of vertices.
		for (int i = 0; i < nVerts; ++i){
			if(vertBlock[i]!=0){
				vertStatus[i] = true;
				nVertExplored++;
			}
		}	
	} else {
		// If vertBlock and nVerts are difference sizes the behaviour is
		// undefined.
		cerr << "Error in " << __PRETTY_FUNCTION__ << ": " << endl
			<< "vertBlock should be empty or the size of mesh.verts";
		RSVS3D_ERROR_ARGUMENT("vertBlock vector incompatible with mesh.");
	}
	currQueue.reserve(nVerts/2);
	nextQueue.reserve(nVerts/2);

   
   	// While Loop, while not all vertices explored
   	while(nVertExplored<nVerts){

	  	// if currQueue is empty start new block
		if(currQueue.size()<1){

			 //cout << "Block " << nBlocks << " - " << nVertExplored << " - " << nVerts << endl;
	   		ii=0;
	   		while(ii<nVerts && vertStatus[ii]){
	   			ii++;
	   		}
	   		if (vertStatus[ii]){
	   			cerr << "Error starting point for loop not found despite max number of vertex not reached" <<endl;
	   			cerr << "Error in " << __PRETTY_FUNCTION__ << endl;
	   			RSVS3D_ERROR_RANGE(" : Starting point for block not found");
	   		}
	   		currQueue.push_back(ii);
	   		nBlocks++;

	   	}
		  // Explore current queue
	   	nCurr=currQueue.size();
	   	for (ii = 0; ii < nCurr; ++ii){
	   		if (!vertStatus[currQueue[ii]]){
	   			vertBlock[currQueue[ii]]=nBlocks;
	   			nEdgesCurr=verts(currQueue[ii])->edgeind.size();
	   			for(jj=0;jj<nEdgesCurr;++jj){
	   				kk=int(edges.isearch(verts(currQueue[ii])->edgeind[jj])->vertind[0]
	   					==verts(currQueue[ii])->index);
	   				nextQueue.push_back(verts.find(
	   					edges.isearch(verts(currQueue[ii])->edgeind[jj])->vertind[kk]));
				   	#ifdef SAFE_ALGO
	   				if (verts.find(
	   					edges.isearch(verts(currQueue[ii])->edgeind[jj])->vertind[kk])==-1){
	   					cerr << "Edge index: " << verts(currQueue[ii])->edgeind[jj] << " vertex index:" <<  
		   				edges.isearch(verts(currQueue[ii])->edgeind[jj])->vertind[kk] << endl;
		   				cerr << "Edge connected to non existant vertex" <<endl;
		   				cerr << "Error in " << __PRETTY_FUNCTION__ << endl;
	   					RSVS3D_ERROR_RANGE(" : Vertex not found");
	   				}
					#endif

	   			}
				vertStatus[currQueue[ii]]=true;
				nVertExplored++;
	   		}
	   	}

		// Reset current queue and set to next queue
		currQueue.clear();
		currQueue.swap(nextQueue);

	}
	return(nBlocks);
}

int mesh::ConnectedVolumes(vector<int> &voluBlock, const vector<bool> &boundaryFaces) const {
	// Fills a vector with a number for each volume corresponding to a
	// group of connected edges it is part of , can be used close surfaces in 2D or volumes
	// in 3D.
	// Uses a flood fill with queue method
	// Boundary faces are faces which stop the the flood exploring through them
	// it is a boolean 

	int nVoluExplored,nVolus, nSurfs,nBlocks,nCurr,nSurfsCurr,ii,jj,kk;
	vector<bool> volStatus; // 1 explored 0 not explored

	vector<int> currQueue, nextQueue; // Current and next queues of indices
	bool boundConstrAct;
	// Preparation of the arrays;
	nVolus=volus.size();
	nSurfs=surfs.size();
	nBlocks=0;
	nVoluExplored=0;

	volStatus.assign(nVolus,false);

	if(int(voluBlock.size())==0){
		// Standard behaviour of grouping all connected volumes when voluBlock
		// is empty
		voluBlock.clear();
		voluBlock.assign(nVolus,0);
	} else if(int(voluBlock.size())==nVolus){
		// Treat volumes with non zero value as already processed.
		// This allows these volumes to separate blocks of volumes.
		for (int i = 0; i < nVolus; ++i){
			if(voluBlock[i]!=0){
				volStatus[i] = true;
				nVoluExplored++;
			}
		}	
	} else {
		// If voluBlock and nVolus are difference sizes the behaviour is
		// undefined.
		cerr << "Error in " << __PRETTY_FUNCTION__ << ": " << endl
			<< "voluBlock should be empty or the size of mesh.volus";
		RSVS3D_ERROR_ARGUMENT("voluBlock vector incompatible with mesh.");
	}
	if(boundaryFaces.size()==0){
		boundConstrAct =false;
	} else if (int(boundaryFaces.size())==nSurfs){
		boundConstrAct =true;
	} else {
		cerr << "Error in " << __PRETTY_FUNCTION__ << ": " << endl
			<< "boundaryFaces should be empty or the size of mesh.surfs";
		RSVS3D_ERROR_ARGUMENT("boundaryFaces vector incompatible with mesh.");
	}
	currQueue.reserve(nVolus/2);
	nextQueue.reserve(nVolus/2);

   
   	// While Loop, while not all volumes explored
   	while(nVoluExplored<nVolus){

	  	// if currQueue is empty start new block
		if(currQueue.size()<1){

			 //cout << "Block " << nBlocks << " - " << nVoluExplored << " - " << nVolus << endl;
	   		ii=0;
	   		while(ii<nVolus && volStatus[ii]){
	   			ii++;
	   		}
	   		if (volStatus[ii]){
	   			cerr << "Error starting point for loop not found despite max number of volume not reached" <<endl;
	   			cerr << "Error in " << __PRETTY_FUNCTION__ << endl;
	   			RSVS3D_ERROR_RANGE(" : Starting point for block not found");
	   		}
	   		currQueue.push_back(ii);
	   		nBlocks++;

	   	}
		  // Explore current queue
	   	nCurr=currQueue.size();
	   	for (ii = 0; ii < nCurr; ++ii){
	   		if (volStatus[currQueue[ii]]){
	   			// Volume already explored
	   			continue;
	   		}
   			voluBlock[currQueue[ii]]=nBlocks;
   			nSurfsCurr=volus(currQueue[ii])->surfind.size();
   			for(jj=0;jj<nSurfsCurr;++jj){
   				int sTempInd = surfs.find(volus(currQueue[ii])->surfind[jj]);
   				if(boundConstrAct && boundaryFaces[sTempInd]){
   					// If a boundary has been provided and this face is one
   					continue;
   				}
   				kk=int(surfs(sTempInd)->voluind[0]==volus(currQueue[ii])->index);
   				if(surfs(sTempInd)->voluind[kk]==0){
   					// if volume index is 0 skip this
   					continue;
   				}
   				nextQueue.push_back(volus.find(surfs(sTempInd)->voluind[kk]));
			   	#ifdef SAFE_ALGO
   				if (volus.find(surfs(sTempInd)->voluind[kk])==-1){
   					cerr << "Surf index: " << volus(currQueue[ii])->surfind[jj] << " volume index:" <<  
	   				surfs(sTempInd)->voluind[kk] << endl;
	   				cerr << "Surf connected to non existant volume" <<endl;
	   				cerr << "Error in " << __PRETTY_FUNCTION__ << endl;
   					RSVS3D_ERROR_RANGE(" : Volume not found");
   				}
				#endif
	   			

   			}
			volStatus[currQueue[ii]]=true;
			nVoluExplored++;
	   		
	   	}

		// Reset current queue and set to next queue
		currQueue.clear();
		currQueue.swap(nextQueue);

	}
	return(nBlocks);
}

coordvec mesh::CalcCentreVolu(int ind) const{
	//takes the position int

	coordvec ret;
	coordvec temp;
	double edgeLength; 
	double voluLength;
	int ii,ni,jj,nj;
	int cSurf,a;
	voluLength=0;
	a=volus.find(ind);
	ni=volus(a)->surfind.size();
	for(ii=0; ii<ni ; ++ii){
		cSurf=surfs.find(volus(a)->surfind[ii]);
		nj=surfs(cSurf)->edgeind.size();
		for(jj=0; jj<nj ; ++jj){
			temp.assign(0,0,0);
			temp.add(verts.isearch(edges.isearch(surfs(cSurf)->edgeind[jj])->vertind[0])->coord);
			temp.substract(verts.isearch(edges.isearch(surfs(cSurf)->edgeind[jj])->vertind[1])->coord);

			edgeLength=temp.CalcNorm();
			voluLength+=edgeLength;

			temp.add(verts.isearch(edges.isearch(surfs(cSurf)->edgeind[jj])->vertind[1])->coord);
			temp.add(verts.isearch(edges.isearch(surfs(cSurf)->edgeind[jj])->vertind[1])->coord);
			temp.mult(edgeLength);

			ret.add(temp.usedata());
		}	
	}
	ret.div(voluLength);
	return(ret);
}


coordvec mesh::CalcPseudoNormalSurf(int ind) const{
	coordvec ret;
	coordvec temp1,temp2,temp3;
	double edgeLength;
	double voluLength;
	vector<int> vertList;
	int jj,nj,cSurf;
	
	voluLength=0;
	cSurf=surfs.find(ind);
	surfs(cSurf)->OrderedVerts(this,vertList);

	nj=vertList.size();
	ret.assign(0,0,0);
	for(jj=0; jj<nj ; ++jj){
		temp1.assign(0,0,0);
		temp2.assign(0,0,0);
		temp1.add(verts.isearch(vertList[(jj+1)%nj])->coord);
		temp2.add(verts.isearch(vertList[(jj+nj-1)%nj])->coord);
		temp1.substract(verts.isearch(vertList[jj])->coord);
		temp2.substract(verts.isearch(vertList[jj])->coord);

		temp3=temp1.cross(temp2.usedata());

		edgeLength=temp3.CalcNorm();
		voluLength+=edgeLength;

		ret.add(temp3.usedata());
	}	
	
	ret.div(voluLength);
	return(ret);
}

void mesh::OrientFaces(){
	/*
	Orients either surfaces or edges depending. 
	*/
	if (this->WhatDim()==3){
		this->OrientSurfaceVolume();
	} else if(this->WhatDim()==2){
		this->OrientEdgeSurface();
	} else {

	}
}

void mesh::OrientEdgeSurface(){
	cerr << "Warning: not coded yet in " << __PRETTY_FUNCTION__ << endl;
}

void mesh::OrientSurfaceVolume(){
	// Orders the surf.voluind [c0 c1] such that the surface normal vector points
	// from c0 to c1
	// This is done by using the surface normals and checking they go towards
	// the centre of the cell

	int nBlocks,ii,jj, ni,nj,kk;
	vector<int> surfOrient;
	vector<bool> isFlip;
	double dotProd;
	coordvec centreVolu, normalVec;


	nBlocks=OrientRelativeSurfaceVolume(surfOrient);
	isFlip.assign(nBlocks,false);
	//========================================
	//  Select direction using coordinate geometry
	//	 use a surface 

	for (ii=1; ii<= nBlocks; ii++){
		jj=-1; nj=surfOrient.size();
		do{
			jj++;
			while(jj<nj && ii!=abs(surfOrient[jj]))
				{jj++;}
			if(jj==nj){ // if the orientation cannot be defined
				dotProd=1.0;
				kk=0;
				cerr << "Warning: Cell orientations could not be computed " << endl;
				cerr << "			in " << __PRETTY_FUNCTION__ << endl;
				break;
			}
			kk=surfs(jj)->voluind[0]==0;
			centreVolu=CalcCentreVolu(surfs(jj)->voluind[kk]);
			normalVec=CalcPseudoNormalSurf(surfs(jj)->index);

			centreVolu.substractfrom(verts.isearch(
				edges.isearch(surfs(jj)->edgeind[0])->vertind[0]
				)->coord);
			dotProd=centreVolu.dot(normalVec.usedata());
		} while (!isfinite(dotProd) || (fabs(dotProd)<numeric_limits<double>::epsilon()));


		isFlip[ii-1]= (((dotProd<0.0) && (kk==0)) || ((dotProd>0.0) && (kk==1)));
	}
	ni=surfOrient.size();
	for(ii=0; ii< ni; ++ii){
		if(isFlip[abs(surfOrient[ii])-1]){
			surfs.elems[ii].FlipVolus();
		}
	}
}

int mesh::OrientRelativeSurfaceVolume(vector<int> &surfOrient){

	int nSurfExplored,nSurfs,nBlocks,nCurr,currEdge,testSurf,relOrient;
	int ii,jj,kk,nj,nk;
	vector<bool> surfStatus; // 1 explored 0 not explored
	vector<vector<int>> orderVert;
	bool isConnec, t0,t1,t3,t4,isFlip;
	vector<int> currQueue, nextQueue, emptVert; // Current and next queues of indices

	// Preparation of the arrays;
	nSurfs=surfs.size();
	

	surfStatus.assign(nSurfs,false);
	surfOrient.assign(nSurfs,0);
	currQueue.reserve(nSurfs/2); // list of positions
	nextQueue.reserve(nSurfs/2);
	orderVert.reserve(nSurfs);
 
	// =======================================
	// Collect surface vertex lists
	emptVert.assign(6,0);
	for(ii=0; ii< nSurfs; ii++){
		orderVert.push_back(emptVert);
		surfs(ii)->OrderedVerts(this,orderVert[ii]);		
	}

	// ========================================
	// Flooding to find relative orientations of surfaces
	// Start from a surf that is in no list, 
	//		look at each its edges 
	//		find adjacent surfaces(checking they share a cell)
	//			-> compute relative orientation (using idea of contra-rotating adjacent surfs) 
	//			-> add surface to queue
	nBlocks=0;
	nSurfExplored=0;
	if(meshDim<3){
		return(nBlocks);
	}
	//cout << " " << nSurfs << " | " ;
	while(nSurfExplored<nSurfs){
		// if currQueue is empty start new block
		//cout << " " << nSurfExplored ;
		if(currQueue.size()<1){
			ii=0;
			while(ii<nSurfs && surfStatus[ii])
				{ ii++; }
			if (ii==nSurfs){
				//cout << " | " << nSurfs << " | " ;
				RSVS3D_ERROR_RANGE(" Start point not found");
			}
			// ii=nSurfs-1;
			// while(ii>=0 && surfStatus[ii])
			// 	{ ii--; }
			// if (ii==-1){
			// 	//cout << " | " << nSurfs << " | " ;
			// 	RSVS3D_ERROR_RANGE(" Start point not found");
			// }
			currQueue.push_back(ii);
			nBlocks++;
			surfStatus[ii]=true;
			nSurfExplored++;
			surfOrient[ii]=nBlocks;

		}
		// Explore current queue
		nCurr=currQueue.size();
		for (ii = 0; ii < nCurr; ++ii){
			nj=surfs(currQueue[ii])->edgeind.size();
			for (jj=0; jj < nj; ++jj){
				currEdge=edges.find(surfs(currQueue[ii])->edgeind[jj]);
				nk=edges(currEdge)->surfind.size();
				for (kk=0; kk < nk; ++kk){
					// tN -> volu[N] is shared 
					testSurf=surfs.find(edges(currEdge)->surfind[kk]);
					t0 = (surfs(testSurf)->voluind[0]==surfs(currQueue[ii])->voluind[0] ||
						surfs(testSurf)->voluind[1]==surfs(currQueue[ii])->voluind[0]) && 
					(surfs(currQueue[ii])->voluind[0]!=0);

					t1 = (surfs(testSurf)->voluind[0]==surfs(currQueue[ii])->voluind[1] ||
						surfs(testSurf)->voluind[1]==surfs(currQueue[ii])->voluind[1]) && 
					(surfs(currQueue[ii])->voluind[1]!=0);
					// if either volume is shared surface is to be flooded
					isConnec = (edges(currEdge)->surfind[kk]!=surfs(currQueue[ii])->index) 
					&& (t0 || t1) && (!surfStatus[testSurf]);

					if(isConnec){
						// Test rotation
						relOrient=OrderMatchLists(orderVert[currQueue[ii]], orderVert[testSurf], 
							edges(currEdge)->vertind[0],edges(currEdge)->vertind[1]);
						// Add to the next queue
						nextQueue.push_back(testSurf);
						nSurfExplored++;
						surfStatus[testSurf]=true;
						surfOrient[testSurf]=-1*relOrient*surfOrient[currQueue[ii]];

						// Flip volumes to match
						t3=(surfs(testSurf)->voluind[0]==surfs(currQueue[ii])->voluind[0]);
						t4=(surfs(testSurf)->voluind[1]==surfs(currQueue[ii])->voluind[1]);
						isFlip=((relOrient==-1) && ((t0 && !t3) || (t1 && !t4))) 
						|| ((relOrient==1) && ((t0 && t3) || (t1 && t4)));
						if(isFlip){
							surfs.elems[testSurf].FlipVolus();
						}
					}
				}
			}
		}
		// Reset current queue and set to next queue
		currQueue.clear();
		currQueue.swap(nextQueue);
	}
	return(nBlocks);
}

void mesh::Scale(){
	std::array<std::array<double, 2>,3> domain;
	for (int i = 0; i < 3; ++i)
	{
		domain[i]={0.0, 1.0};
	}
	this->Scale(domain);
}

void mesh::Scale(const std::array<std::array<double, 2>,3> &domain){
	/*
	Scale the mesh to be in the range domain where domain is a 2D array:
	[
	x [lb ub]
	y [lb ub]
	z [lb ub]
	]
	*/
	std::array<std::array<double, 2>,3> currDomain;
	std::array<std::array<double, 3>,3> transformation;
	int nVerts, ii, jj;


	nVerts = int(this->verts.size());

	// Find the current range of values of the grid
	for (int i = 0; i < 3; ++i)
	{
		currDomain[i]={std::numeric_limits<double>::infinity(), 
			-std::numeric_limits<double>::infinity()};
	}
	for(ii=0; ii<nVerts; ++ii){
		for(jj=0; jj<3; ++jj){
			currDomain[jj][0]= currDomain[jj][0] <= this->verts(ii)->coord[jj] ?
				currDomain[jj][0]:this->verts(ii)->coord[jj];
			currDomain[jj][1]= currDomain[jj][0] >= this->verts(ii)->coord[jj] ?
				currDomain[jj][1]:this->verts(ii)->coord[jj];
		}
	}

	// Calculate modification parameters
	for (int i = 0; i < 3; ++i)
	{
		transformation[i][0]=domain[i][0];
		transformation[i][1]=currDomain[i][0];
		transformation[i][2]=(domain[i][1]-domain[i][0])
			/(currDomain[i][1]-currDomain[i][0]);
		// Handle flat domains by ignoring the dimension
		if(!isfinite(transformation[i][2])){
			transformation[i][2]=1.0;
		}
	}

	// Recalculate grid vertex positions
	for(ii=0; ii<nVerts; ++ii){
		for(jj=0; jj<3; ++jj){
			this->verts.elems[ii].coord[jj]=transformation[jj][0]
				+((this->verts(ii)->coord[jj]-transformation[jj][1])
				*transformation[jj][2]);
		}
	}
}

void mesh::ReturnBoundingBox(std::array<double,3> &lowerB, 
	std::array<double,3> &upperB){
	int count = this->verts.size();

	lowerB.fill(INFINITY);
	upperB.fill(-INFINITY);
	for (int i = 0; i < count; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			lowerB[j] = lowerB[j] < this->verts(i)->coord[j] ?
				lowerB[j] : this->verts(i)->coord[j];
			upperB[j] = upperB[j] > this->verts(i)->coord[j] ?
				upperB[j] : this->verts(i)->coord[j];
		}
	}
}

void mesh::LoadTargetFill(const std::string &fileName){
	
	vector<double> fillVals;
	double ft;
	int nElms, count;
	ifstream file;
	
	file.open(fileName);
	CheckFStream(file, __PRETTY_FUNCTION__, fileName);

	if (this->WhatDim()==3){
		nElms = this->volus.size();
	} else if(this->WhatDim()==2){
		nElms = this->surfs.size();
	} else {
		RSVS3D_ERROR_ARGUMENT("Dimensionality of the target not supported.");
	}

	fillVals.reserve(nElms);
	count=0;
	while(!file.eof() && count<nElms) {
		file >> ft;
		fillVals.push_back(ft);
		++count;
	}
	cout << "fill loaded : " ;
	DisplayVector(fillVals);
	cout << endl; 
	if (this->WhatDim()==3){
		for (int i = 0; i < nElms; ++i) {
			this->volus[i].target = fillVals[i%count];
		}
	} else if(this->WhatDim()==2){
		for (int i = 0; i < nElms; ++i) {
			this->surfs[i].target = fillVals[i%count];
		}
	}
	this->PrepareForUse();
}

/**
 * @brief      Adds boundaries alond max and min xyz planes.
 *	Arguments
 * @param[in]  lb    lower boundary vector of 3 doubles.
 * @param[in]  ub    upper boundary vector of 3 doubles.
 *
 * @return     List of vertex indices in the mesh which lie outside.
 * 
 * Raises:
 * 	- logic_error, 
 * 	-
 * 	-
 * 	
 * 	Process:
 * 	THis method could be readily refactored to allow treatment of more
 *	complex boundaries

	Steps:
		1 - Identify vertices lying outside
		2 - Indentify connectors lying on boundary
			a - edges
			b - surfs
			c - volus
		3 - Introduce boundary vertices (BV)
		4 - Connect those BV to form new boundary edges (BE)
		5 - Assemble BEs inside a volu into a boundary surf (BS)
		(This process is similar to the voronoisation)
 */
std::vector<int> mesh::AddBoundary(const std::vector<double> &lb,
		const std::vector<double> &ub){
	/*
		Adds elements on a boundary defined by upper bounds and lower bounds.

		THis method could be readily refactored to allow treatment of more
		complex boundaries

		Steps:
			1 - Identify vertices lying outside
			2 - Indentify connectors lying on boundary
				a - edges
				b - surfs
				c - volus
			3 - Introduce boundary vertices (BV)
			4 - Connect those BV to form new boundary edges (BE)
			5 - Assemble BEs inside a volu into a boundary surf (BS)
			(This process is similar to the voronoisation)
		
		Args: 
			lb : lower bounds in vector form.
			up : upper bounds in vector form.

		Returns:
			A list of vertex indices which lie outside the boundary.
	*/
	mesh meshAdd; 
	std::vector<int>  vertOutInd; // return
	#ifdef SAFE_ALGO
	bool flagErr1, flagErr2;
	#endif
	std::vector<bool> vertOut;
	std::vector<bool> edgeOut;
	std::vector<int>  edgeBound;
	std::vector<bool> surfOut;
	std::vector<int>  surfBound;
	std::vector<bool> voluOut;
	std::vector<int>  voluBound;

	int vertSize, edgeSize, edgeSize2, surfSize, surfSize2,
		voluSize, voluSize2, count, countPrev;

	vertSize=this->verts.size();
	edgeSize=this->edges.size();
	surfSize=this->surfs.size();
	voluSize=this->volus.size();
	vertOut.assign(vertSize, false);
	edgeBound.reserve(edgeSize);
	surfBound.reserve(surfSize);
	voluBound.reserve(voluSize);


	meshAdd.Init(0,0,0,0);

	// Mark outer vertices
	for (int i = 0; i < vertSize; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			if(this->verts(i)->coord[j]>ub[j] 
				|| this->verts(i)->coord[j]<lb[j]){
				vertOut[i] = true;
				break;
			}
		}
	}
	for (int i = 0; i < edgeSize; ++i)
	{ // for each edge if it's vertices are one in one out.
		if(vertOut[this->verts.find(this->edges(i)->vertind[0])] 
			^ vertOut[this->verts.find(this->edges(i)->vertind[1])]){
			edgeBound.push_back(this->edges(i)->index);
		}
	}
	count = edgeBound.size();
	auto tempVec = this->edges.find_list(edgeBound);
	surfBound = ConcatenateVectorField(this->edges, &edge::surfind, tempVec);
	sort(surfBound);
	unique(surfBound);
	tempVec.clear();
	tempVec = this->surfs.find_list(surfBound);
	voluBound = ConcatenateVectorField(this->surfs, &surf::voluind, tempVec);
	sort(voluBound);
	unique(voluBound);
	if(voluBound.size() && voluBound[0]==0){
		voluBound[0] = voluBound.back();
		voluBound.pop_back();
	}
	// Create all the mesh containers for the new mesh parts
	meshAdd.Init(edgeBound.size(), edgeBound.size()+surfBound.size(),
		surfBound.size()+voluBound.size(), voluBound.size());
	#ifdef RSVS_VERBOSE
	this->displight();
	#endif //RSVS_VERBOSE
	meshAdd.PopulateIndices();
	this->MakeCompatible_inplace(meshAdd);
	// meshAdd.disp();
	this->Concatenate(meshAdd);
	this->HashArray();
	#ifdef RSVS_VERBOSE
	this->displight();
	#endif //RSVS_VERBOSE

	// Handle edges:
	//   1 - Double the edge.
	//   2 - Add the vertex
	//   3 - Connect the vertex
	// 4 - Assign edges as internal or external
	count = edgeBound.size();
	tempVec = this->edges.find_list(edgeBound);
	for (int i = 0; i < count; ++i)
	{
		// Replicate edge (changing index)
		auto tempInd = this->edges[i+edgeSize].index;
		this->edges[i+edgeSize] = this->edges[tempVec[i]];
		this->edges[i+edgeSize].index = tempInd;
		for (auto temp : this->surfs.find_list(this->edges[i+edgeSize].surfind)){
			this->surfs[temp].edgeind.push_back(tempInd);
		}
		// Compute vertex position
		int iIn=0, iOut=0, subIn, subOut;
		iOut += !vertOut[this->verts.find(this->edges[tempVec[i]].vertind[0])];
		iIn += vertOut[this->verts.find(this->edges[tempVec[i]].vertind[0])];
		#ifdef SAFE_ALGO
		if (!(vertOut[this->verts.find(this->edges[tempVec[i]].vertind[0])]
				^ vertOut[this->verts.find(this->edges[tempVec[i]].vertind[1])])){
			RSVS3D_ERROR_LOGIC("Edge Should be connected to one inner and one outer"
				"vertex, this is not the case.");
		}
		flagErr1 = vertOut[this->verts.find(this->edges[tempVec[i]].vertind[iIn])];
		flagErr2 = !vertOut[this->verts.find(this->edges[tempVec[i]].vertind[iOut])];
		if (flagErr1 && flagErr2) {
			RSVS3D_ERROR_LOGIC("Coord In is Out and Coord out is in");
		} else if(flagErr1){
			RSVS3D_ERROR_LOGIC("Coord In is Out");
		} else if(flagErr2){
			RSVS3D_ERROR_LOGIC("Coord out is in");
		}
		#endif // SAFE_ALGO
		subIn = this->verts.find(this->edges[tempVec[i]].vertind[iIn]);
		subOut = this->verts.find(this->edges[tempVec[i]].vertind[iOut]);
		meshhelp::PlaceBorderVertex(
			this->verts[subIn].coord, this->verts[subOut].coord,
			lb, ub, this->verts[i+vertSize].coord
			);

		int count2 = this->verts[subOut].edgeind.size();
		for (int j = 0; j < count2; ++j)
		{
			if(this->verts[subOut].edgeind[j]==this->edges[tempVec[i]].index){
				this->verts[subOut].edgeind[j]=this->edges[i+edgeSize].index;
			}
		}

		this->edges[i+edgeSize].vertind[iIn]=this->verts[i+vertSize].index;
		this->edges[tempVec[i]].vertind[iOut]=this->verts[i+vertSize].index;
		this->verts[i+vertSize].edgeind.push_back(this->edges[tempVec[i]].index);
		this->verts[i+vertSize].edgeind.push_back(this->edges[i+edgeSize].index);
		vertOut.push_back(false); // Add new vertex as an inside vertex
		this->ArraysAreHashed();
	}
	this->verts.HashArray();
	this->edges.HashArray();
	edgeSize2 = edgeSize+count;
	edgeOut.assign(edgeSize2, false);
	for (int i = 0; i < edgeSize2; ++i)
	{ // assigns an edge as outside if either of it's vertices are out
		edgeOut[i] = vertOut[this->verts.find(this->edges(i)->vertind[0])] 
			|| vertOut[this->verts.find(this->edges(i)->vertind[1])] ;
	}
	// Handle surfaces
	//  1 - double the surface
	//  2 - add the new edge
	//  3 - close the two surfaces
	// Need to handle multiple cuts in the same surface
	countPrev=0;
	do {
		count = surfBound.size();
		tempVec.clear(); tempVec = this->surfs.find_list(surfBound);
		//===============================================
		for (int i = countPrev; i < count; ++i)
		{
			// Replicate surf (changing index)
			auto tempInd = this->surfs[i+surfSize].index;
			this->surfs[i+surfSize] = this->surfs[tempVec[i]];
			this->surfs[i+surfSize].index = tempInd;

			meshhelp::SplitBorderSurfaceEdgeind(*this, edgeOut, 
				this->surfs[tempVec[i]].edgeind, this->surfs[i+surfSize].edgeind);

			// Find the vertind connectivity of the new edge.
			this->edges[edgeSize2+i].vertind = 
				meshhelp::FindVertInFromEdgeOut(*this, vertOut, 
				this->surfs(i+surfSize)->edgeind,
				this->surfs(tempVec[i])->edgeind);
			
			if(this->edges[edgeSize2+i].vertind.size()>2){
				this->ArraysAreHashed();
				meshhelp::HandleMultiSurfaceSplit(*this, 
					this->surfs.elems[tempVec[i]].edgeind,
					this->surfs.elems[surfSize+i].edgeind,
					this->edges.elems[edgeSize2+i].vertind);
				surfBound.push_back(this->surfs(tempVec[i])->index);
			}
			// Switch the edge connectivity of certain edges in scoped mode
			this->ArraysAreHashed();
			this->SwitchIndex(6, this->surfs(tempVec[i])->index,
				this->surfs(i+surfSize)->index,this->surfs(i+surfSize)->edgeind);
			// assign new edge.surfind and surf.edgeind connectivity
			this->edges[edgeSize2+i].surfind.push_back(this->surfs[i+surfSize].index);
			this->edges[edgeSize2+i].surfind.push_back(this->surfs[tempVec[i]].index);
			this->surfs[i+surfSize].edgeind.push_back(this->edges[edgeSize2+i].index);
			this->surfs[tempVec[i]].edgeind.push_back(this->edges[edgeSize2+i].index);
			// Assign the vertices edgeind connectivity
			int count2 = this->edges[edgeSize2+i].vertind.size();
			for (int j = 0; j < count2; ++j)
			{
				this->verts.elems[this->verts.find(this->edges[edgeSize2+i].vertind[j])]
					.edgeind.push_back(this->edges[edgeSize2+i].index);
			}
			edgeOut.push_back(false); // New edge is inside
			this->ArraysAreHashed();
		}
		//===============================================
		// if surfBound has grown surfaces require cutting more than once.
		countPrev = count;
		count=surfBound.size();
		if(count!=countPrev){
			meshAdd.Init(0, count-countPrev,
				count-countPrev, 0);
			meshAdd.PopulateIndices();
			this->SetMaxIndex();
			meshAdd.SetMaxIndex();
			this->MakeCompatible_inplace(meshAdd);
			this->Concatenate(meshAdd);
			this->HashArray();
		}
	} while (count != countPrev);
	this->verts.HashArray();
	this->edges.HashArray();
	this->surfs.HashArray();
	surfSize2 = surfSize+count;
	surfOut.assign(surfSize2, false);
	for (int i = 0; i < surfSize2; ++i)
	{ // assigns a surf as outside if any of it's edges are out
		int count2 = this->surfs(i)->edgeind.size();
		for (int j = 0; j < count2; ++j)
		{
			surfOut[i] = surfOut[i] 
				|| edgeOut[this->edges.find(this->surfs(i)->edgeind[j])] ;
		}
	}
	// Handle Volumes
	//  1 - Double the volume
	//  2 - Add the new surface
	//  3 - Close the two volumes
	count = voluBound.size();
	tempVec.clear(); tempVec = this->volus.find_list(voluBound);
	for (int i = 0; i < count; ++i)
	{
		// Replicate volu (changing index)
		auto tempInd = this->volus[i+voluSize].index;
		this->volus[i+voluSize] = this->volus[tempVec[i]];
		this->volus[i+voluSize].index = tempInd;
		
		// Split surfind connectivity
		meshhelp::SplitBorderVolumeSurfind(*this, surfOut, 
			this->volus[tempVec[i]].surfind, this->volus[i+voluSize].surfind);
		// Switch the volumes connectivity of certain surfaces in scoped mode
		this->ArraysAreHashed();
		this->SwitchIndex(7, this->volus(tempVec[i])->index,
			this->volus(i+voluSize)->index,this->volus(i+voluSize)->surfind);
		// ^ triggers a warning which will need to be suppressed
		// Find the edgeind connectivity of the new surf.
		this->surfs[surfSize2+i].edgeind = 
			meshhelp::FindEdgeInFromSurfOut(*this, edgeOut,
			this->volus(i+voluSize)->surfind);
		// Connect the new surface to the two volumes which were split
		this->surfs[surfSize2+i].voluind[0] = this->volus[i+voluSize].index; 
		this->surfs[surfSize2+i].voluind[1] = this->volus[tempVec[i]].index;
		this->volus[tempVec[i]].surfind.push_back(this->surfs[surfSize2+i].index);
		this->volus[i+voluSize].surfind.push_back(this->surfs[surfSize2+i].index);
		// Assign the edges surfind connectivity
		int count2 = this->surfs(surfSize2+i)->edgeind.size();
		this->ArraysAreHashed();
		for (int j = 0; j < count2; ++j)
		{
			this->edges.elems[this->edges.find(this->surfs[surfSize2+i].edgeind[j])]
				.surfind.push_back(this->surfs[surfSize2+i].index);
		}
		surfOut.push_back(false); // New edge is inside
		this->ArraysAreHashed();
	}
	this->HashArray();
	voluSize2 = voluSize+count;
	voluOut.assign(voluSize2, false);
	for (int i = 0; i < voluSize2; ++i)
	{ // assigns a surf as outside if any of it's surfs are out
		for (auto temp : this->volus(i)->surfind)
		{
			voluOut[i] = voluOut[i] 
				|| surfOut[this->surfs.find(temp)] ;
		}
	}
	this->HashArray();
	this->TightenConnectivity();
	this->TestConnectivityBiDir(__PRETTY_FUNCTION__);
	this->PrepareForUse();
	vertOutInd.reserve(this->verts.size());
	count = this->verts.size();
	for (int i = 0; i < count; ++i)
	{
		if(vertOut[i]){
			vertOutInd.push_back(this->verts(i)->index);
		}
	}

	return(vertOutInd);
}

void mesh::Crop(vector<int> indList, int indType){
	/*
	Crops a mesh removing all the elements described by indlist.

	The cropping process is such that it leaves a functional mesh at the end.
	This process is performed in two cycles: one upward removing the elements
	of higher dimensionality connected to the identified objects. Then the
	process goes downwards removing all the orphan elements which have been
	created during the cleanup.

	Args:
		indList : List of indices of a given element type that need to 
			be removed.
		indType : Type of the indices provided in the list. Default is 1
			which corresponds to vertices. vertex (indType=1), edge (2),
			surf (3), volu (4)

	Raises:
		invalid_argument : if the `indType` is not between 1 and 4.

	*/

	// Check validity of input
	if (indType>4 || indType<1){
		RSVS3D_ERROR_ARGUMENT("Type of index not recognised: \n"
			" vertex (indType=1), edge (2), surf (3), volu (4)");
	}

	std::vector<int> vertDel, edgeDel, surfDel, voluDel;

	int i=1;
	if(indType==i++){
		vertDel = indList;
	} else if (indType==i++) {
		edgeDel = indList;
	} else if (indType==i++) {
		surfDel = indList;
	} else if (indType==i++) {
		voluDel = indList;
	}

	switch(indType){
		case 1:{
			auto vertSub = this->verts.find_list(vertDel);
			edgeDel = ConcatenateVectorField(this->verts, &vert::edgeind, vertSub);
			[[gnu::fallthrough]];
		}
		case 2:{
			auto edgeSub = this->edges.find_list(edgeDel);
			surfDel = ConcatenateVectorField(this->edges, &edge::surfind, edgeSub);
			[[gnu::fallthrough]];
		}
		case 3:{
			auto surfSub = this->surfs.find_list(surfDel);
			voluDel = ConcatenateVectorField(this->surfs, &surf::voluind, surfSub);
		}
	}
	// Remove indices from connectivity and propagate downwards from the volumes.
	// Remove volumes
	sort(voluDel);
	unique(voluDel);
	int count = voluDel.size();
	if (count>0){
		if (voluDel[0]==0){voluDel.erase(voluDel.begin());}
		count--;
		for (int i = 0; i < count; ++i)	{this->RemoveIndex(4, voluDel[i]);}
	}
	// Expand surf removal list
	count = this->surfs.size();
	for (int i = 0; i < count; ++i)
	{
		if(this->surfs(i)->voluind.size()==0){
			surfDel.push_back(this->surfs(i)->index);
		} else if (this->surfs(i)->voluind.size()==1){
			this->surfs[i].voluind.push_back(0);
		}
	}
	this->ArraysAreHashed();
	sort(surfDel);
	unique(surfDel);
	count = surfDel.size();
	for (int i = 0; i < count; ++i)	{this->RemoveIndex(3, surfDel[i]);}

	count = this->edges.size();
	for (int i = 0; i < count; ++i)
	{
		if(this->edges(i)->surfind.size()==0){
			edgeDel.push_back(this->edges(i)->index);
		} 
	}
	sort(edgeDel);
	unique(edgeDel);
	count = edgeDel.size();
	for (int i = 0; i < count; ++i)	{this->RemoveIndex(2, edgeDel[i]);}

	count = this->verts.size();
	for (int i = 0; i < count; ++i)
	{
		if(this->verts(i)->edgeind.size()==0){
			vertDel.push_back(this->verts(i)->index);
		} 
	}
	sort(vertDel);
	unique(vertDel);
	count = vertDel.size();
	for (int i = 0; i < count; ++i)	{this->RemoveIndex(1, vertDel[i]);}


	this->verts.remove(vertDel);
	this->edges.remove(edgeDel);
	this->surfs.remove(surfDel);
	this->volus.remove(voluDel);

	this->PrepareForUse();
	this->RemoveSingularConnectors();
	this->PrepareForUse();
}

void mesh::CropAtBoundary(const vector<double> &lb, const vector<double> &ub){
	/*
	Crops a mesh along a specified lower bound and upper bound.

	This function uses `mesh::AddBoundary` and `mesh::Crop` to perform this
	action.
	*/
	auto vecDel = this->AddBoundary(lb, ub);
	this->Crop(vecDel, 1);
}


int OrderMatchLists(const vector<int> &vec1, const vector<int> &vec2, int p1, int p2){
	// compares the list vec1 and vec2 returning 
	// 1 if indices p1 and p2 appear in the same order 
	// -1 if indices p1 and p2 appear in opposite orders
	int ii, n, ord1, ord2, retVal, kk;

	
	ord1=0;kk=0;
	n=vec1.size();
	for(ii=0; ii<n; ++ii){
		if (vec1[ii]==p1) {
			ord1+=ii;
			kk++;
		} else if (vec1[ii]==p2) {
			ord1+=-ii;
			kk++;
		}
	}
	if(ord1>1){ord1=-1;}
	if(ord1<-1){ord1=1;}
	if (kk!=2) {
		cerr << "Error : indices were not found in lists " << endl;
		cerr << " 	p1 and/or p2 did not appear in vec " << endl;
		cerr << " 	in " << __PRETTY_FUNCTION__ << endl;
		RSVS3D_ERROR_ARGUMENT("Incaompatible list and index");
	}

	ord2=0;kk=0;
	n=vec2.size();
	for(ii=0; ii<n; ++ii){
		if (vec2[ii]==p1) {
			ord2+=ii;
			kk++;
		} else if (vec2[ii]==p2) {
			ord2+=-ii;
			kk++;
		}
	}
	if(ord2>1){ord2=-1;}
	if(ord2<-1){ord2=1;}
	if (kk!=2) {
		cerr << "Error : indices were not found in lists " << endl;
		cerr << " 	p1 and/or p2 did not appear in vec " << endl;
		cerr << " 	in " << __PRETTY_FUNCTION__ << endl;
		RSVS3D_ERROR_ARGUMENT("Incaompatible list and index");
	}

	retVal=(ord1==ord2)*2-1;

	return(retVal);
}

void ConnVertFromConnEdge(const mesh &meshin, const vector<int> &edgeind, vector<int> &vertind){
	// Returns a list of connected vertices matching a list of connected edges
	bool flag;
	int kk,ll,nEdge; 

	nEdge=int(edgeind.size());
	vertind.reserve(nEdge);
	ll=0;
	kk=0;
	flag=(meshin.edges.isearch(edgeind[kk])->vertind[ll]
		==meshin.edges.isearch(edgeind[kk+1])->vertind[0]) | 
	(meshin.edges.isearch(edgeind[kk])->vertind[ll]
		==meshin.edges.isearch(edgeind[kk+1])->vertind[1]);
	if(!flag){
		ll=((ll+1)%2);
	}
	for(kk=0; kk<nEdge;++kk){
		vertind.push_back(meshin.edges.isearch(edgeind[kk])->vertind[ll]);
		
		if(kk<(nEdge-1)){
			flag=(meshin.edges.isearch(edgeind[kk])->vertind[0]
				==meshin.edges.isearch(edgeind[kk+1])->vertind[ll]) | 
			(meshin.edges.isearch(edgeind[kk])->vertind[1]
				==meshin.edges.isearch(edgeind[kk+1])->vertind[ll]);
			ll=((ll+1)%2)*(flag)+ll*(!flag);
		}
	}
}


void CropMeshGreedy(mesh &meshin, const std::vector<double> &lb,
		const std::vector<double> &ub){
	/*
	Crops a mesh to only the elements inside the cropBox

	Anything impinging on the cropBox is deleted.
	Steps:
		1 - Find vertices out of the box
		2 - Delete those vertices
		3 - Propagate the deletion to the higher level containers
		4 - Propagate back to the lower level containers removing empty
			connectors.
	*/
	int nVerts;
	std::vector<int> vertDel;
	
	// Primary deletion (upwards from the vertices)
	nVerts = meshin.verts.size();
	for (int i = 0; i < nVerts; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			if(meshin.verts(i)->coord[j]>ub[j] 
				|| meshin.verts(i)->coord[j]<lb[j]){
				vertDel.push_back(meshin.verts(i)->index);
				break;
			}
		}
	}
	meshin.Crop(vertDel, 1);
}

mesh Points2Mesh(const std::vector<double> &vecPts, int nProp){
		/*
		Takes in a set of points and returns a mesh of points ready for
		voronisation.
		*/
		mesh meshpts;

		if((vecPts.size()%nProp)!=0){
			std::cerr << "Error in: " << __PRETTY_FUNCTION__ << std::endl;
			RSVS3D_ERROR_ARGUMENT("Vector of points is an invalid size");
		}
		int nPts = vecPts.size()/nProp;
		meshpts.Init(nPts, 0, 0, 0);
		meshpts.PopulateIndices();
		for (int i = 0; i < nPts; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				meshpts.verts[i].coord[j] = vecPts[i*nProp+j]; 
			}
		}
		meshpts.verts.PrepareForUse();
		return (meshpts);
	}

namespace meshhelp {
	/*
	A namespace of helper function not designed to be propagated throughout
	*/
	void PlaceBorderVertex(const std::vector<double> &coordIn, 
		const std::vector<double> &coordOut,
		const std::vector<double> &lb,
		const std::vector<double> &ub, std::vector<double> &coordTarg){

		for (int i = 0; i < 3; ++i)
		{
			coordTarg[i] = coordIn[i]-coordOut[i];
		}



		double lStep = meshhelp::ProjectRay(3, vector<vector<double>>({lb, ub}), 
				coordTarg, coordIn,0.0);

		#ifdef SAFE_ALGO
		if(lStep>+0.0 || lStep<-1.001){
			cout << "Error in: " << __PRETTY_FUNCTION__ << endl << " coordIn : ";
			DisplayVector(coordIn);
			cout << " coordOut : ";
			DisplayVector(coordOut);
			cout << " coordTarg : ";
			DisplayVector(coordTarg);
			cout << endl << "lStep : " << lStep << endl; 
			RSVS3D_ERROR_LOGIC("The lstep for a vertex should be between 0"
				" and 1 for it to lie on the original edge.");
		}
		#endif // SAFE_ALGO


		for (int i = 0; i < 3; ++i)
		{
			coordTarg[i] = coordTarg[i]*lStep + coordIn[i];
		}

	}

	void SplitBorderSurfaceEdgeind(const mesh &meshin, const std::vector<bool> &edgeOut, 
		std::vector<int> &vecconnIn, std::vector<int> &vecconnOut){
		/*
		Splits the connectivity of a surface lying on the border into 
		those edges lying inside and outside of the new forming boundary

		this needs to be templated
		*/

		vecconnOut.clear();
		int count = vecconnIn.size();
		auto temp = meshin.edges.find_list(vecconnIn, true);
		vecconnIn.clear();
		for (int i = 0; i < count; ++i)
		{
			if(temp[i]==-1){
				std::cerr << "Error in : " << __PRETTY_FUNCTION__ << std::endl;
				RSVS3D_ERROR_ARGUMENT("Edge index in connectivity list not found");
			}
			if(edgeOut[temp[i]]){
				vecconnOut.push_back(meshin.edges(temp[i])->index);
			} else {
				vecconnIn.push_back(meshin.edges(temp[i])->index);
			}
		}
	}

	/**
	 * @brief      Handles case for AddBoundary where a surface more than one split. 
	 *
	 * @param      meshin      The mesh
	 * @param      edgeindOld  the edgeind internal
	 * @param      edgeindNew  the edgeind out of the boundary (that will be split)
	 * @param      vertindNew  the vertex index that must be trimmed to be of size 2
	 */
	void HandleMultiSurfaceSplit(mesh &meshin, 
		vector<int> &edgeindOld, vector<int> &edgeindNew,
		vector<int> &vertindNew){
		sort(edgeindNew); 
		unique(edgeindNew);
		auto edgeOrig = edgeindNew;
		int outFlag=OrderEdgeList(edgeindNew, meshin,
			false, false, &edgeOrig);
		// trim edgeind
		if(outFlag!=1){
			edgeindNew.erase(
				edgeindNew.begin()-outFlag,
				edgeindNew.end()
				);
		}
		for (auto edgeIndAll : edgeOrig){
			bool flag=true;
			for (auto edgeindOut : edgeindNew){
				if(edgeIndAll==edgeindOut){
					flag=false;
					break;
				}
			}
			if(flag){
				edgeindOld.push_back(edgeIndAll);
			}
		}
		if(outFlag>=0){
			RSVS3D_ERROR_LOGIC("The result of OrderEdgeList should be "
				"negative to indicate an open surface");
		}
		std::array<int, 2> inds;
		int k=0;
		int nVertExtra = vertindNew.size();
		for (auto edgei : edgeindNew){
			int temp = meshin.edges.find(edgei);
			for (int j = 0; j < 2; ++j)
			{
				for (int l = 0; l < nVertExtra; ++l)
				{
					if(vertindNew[l]
						==meshin.edges(temp)->vertind[j]){
						inds[k++]=l;
						break;
					}
				}
				if(k==2){break;}
			}
			if(k==2){break;}
		}
		vertindNew = {vertindNew[inds[0]], vertindNew[inds[1]]};

	}

	void SplitBorderVolumeSurfind(const mesh &meshin, const std::vector<bool> &edgeOut, 
		std::vector<int> &vecconnIn, std::vector<int> &vecconnOut){
		/*
		Splits the connectivity of a volume lying on the border into 
		those surfs lying inside and outside of the new forming boundary

		this needs to be templated
		*/

		vecconnOut.clear();
		int count = vecconnIn.size();
		auto temp = meshin.surfs.find_list(vecconnIn, true);
		vecconnIn.clear();
		for (int i = 0; i < count; ++i)
		{
			if(temp[i]==-1){
				std::cerr << "Error in : " << __PRETTY_FUNCTION__ << std::endl;
				RSVS3D_ERROR_ARGUMENT("Edge index in connectivity list not found");
			}
			if(edgeOut[temp[i]]){
				vecconnOut.push_back(meshin.surfs(temp[i])->index);
			} else {
				vecconnIn.push_back(meshin.surfs(temp[i])->index);
			}
		}
	}

	std::vector<int> FindVertInFromEdgeOut(const mesh &meshin, const std::vector<bool> &vertOut,
		const std::vector<int> &edgeList, const std::vector<int> &edgeListCheck){
		/*
		Returns the vertices which are in (false in vertOut) from a list of edges
		 which are outside.

		 This may be larger than two as edgeList are not guaranteed to be part of
		 a single surface. If there are two cuts this would yield 4 vertices.
		*/
		std::vector<int> vertList;
		int count = edgeList.size();
		for (int i = 0; i < count; ++i)
		{
			for (int j = 0; j < 2; ++j)
			{
				int temp = meshin.edges.isearch(edgeList[i])->vertind[j];
				if(!vertOut[meshin.verts.find(temp)]){
					vertList.push_back(temp);
				}
			}
		}
		sort(vertList);
		unique(vertList);
		if(vertList.size()!=2){
			// DisplayVector(vertList);
			std::vector<int> vertListTemp={};
			vertList.swap(vertListTemp);

			for (auto k : edgeListCheck)
			{ // 
				for (auto i : meshin.edges.isearch(k)->vertind){
					for (int j : vertListTemp){
						if(i==j){
							vertList.push_back(i);
							break;
						}
					}
				}
			}
		}
		sort(vertList);
		unique(vertList);

		return(vertList);
	}
	std::vector<int> FindEdgeInFromSurfOut(const mesh &meshin, const std::vector<bool> &edgeOut,
		std::vector<int> surfList){
		/*
		Returns the edgeices which are in (false in edgeOut) from a list of surfs
		 which are outside.

		 Needs to be templated
		*/
		std::vector<int> edgeList;
		int count = surfList.size();
		for (int i = 0; i < count; ++i)
		{
			for (auto temp : meshin.surfs.isearch(surfList[i])->edgeind)
			{	
				if(!edgeOut[meshin.edges.find(temp)]){
					edgeList.push_back(temp);
				}
			}
		}

		return(edgeList);
	}

	int VertexInVolume(const mesh &meshin, const vector<double> testCoord, bool needFlip){
		/*
		 for each point
		 Start at a surface (any)
		 Find which side it is on
		 -> gives me a candidate volume 
		 -> candidate volume is stored and provides a list of 
		 surfaces to test against
		 -> if candidate is 0 and outside of the mesh convex its done. 	
		 */

		if(meshin.volus.size()==0){
			return 0;
		}
		coordvec temp1, temp2;
		std::array<int, 3> surfacePoints = {0,0,0};
		bool isInCandidate=false, candidateIsValid;
		int volumeCandidate=0;
		int pointInVolu, alternateVolume, nVolusExplored=0, nVolus;
		double distToPlane, multiplierDist;
		multiplierDist = needFlip ? -1.0 : 1.0;
		nVolus=meshin.volus.size();
		volumeCandidate=meshin.volus(0)->index;
		do{
			// assume the vertex is out of the candidate
			isInCandidate=true;
			candidateIsValid=false;
			// for each face of the candidate volume
			for(auto surfCurr : meshin.volus.isearch(volumeCandidate)->surfind){
				#ifdef SAFE_ALGO
				// Checks that ordering and edgeind are correct
				if(!meshin.surfs.isearch(surfCurr)->isready(true)){
					RSVS3D_ERROR_ARGUMENT("Surfaces of `meshin` need to be ordered."
						" Run meshin.OrderEdges() before calling");
				}
				if(meshin.surfs.isearch(surfCurr)->edgeind.size()<3){
					std::string strErr;
					strErr = "Surface "+to_string(surfCurr)+" has less than 3 edges.";
					RSVS3D_ERROR_ARGUMENT(strErr.c_str());
				}
				#endif

				// assign surfacePoints the indices of 3 non coincident points appearing
				// in the surface
				int surfacePointsIdentified = 0;
				for (auto edgeCurr : meshin.surfs.isearch(surfCurr)->edgeind)
				{
					if(meshin.edges.isearch(edgeCurr)->IsLength0(meshin)){
						continue;
					}
					if(surfacePointsIdentified==0){
						surfacePoints[0] = meshin.edges.isearch(edgeCurr)->vertind[0];
						surfacePoints[1] = meshin.edges.isearch(edgeCurr)->vertind[1];
						surfacePointsIdentified=2;
					} else {
						for (auto vertCurr : meshin.edges.isearch(edgeCurr)->vertind){
							if(vertCurr!=surfacePoints[0] 
								&& vertCurr!=surfacePoints[1]
								&& !(meshhelp::IsVerticesDistance0(meshin, 
									{vertCurr, surfacePoints[0]}))
								&& !(meshhelp::IsVerticesDistance0(meshin, 
									{vertCurr, surfacePoints[1]}))
							){
								surfacePoints[2] = vertCurr;
								surfacePointsIdentified++;
								break;
							}
						}
					}
					if(surfacePointsIdentified==3){
						break;
					}
				}
				if(surfacePointsIdentified!=3){
					// if the surface has less than 3 non-coincident points it is 
					// degenerate and equivalent to an edge. It does not invalidate or
					// validate the candidate
					int jj=meshin.surfs.isearch(surfCurr)->voluind[0]==volumeCandidate;
					alternateVolume = meshin.surfs.isearch(surfCurr)->voluind[jj];
					continue;
				}
				// candidate is valid if it has one non-degenerate face
				candidateIsValid=true; 
				distToPlane = VertexDistanceToPlane(
					meshin.verts.isearch(surfacePoints[0])->coord,
					meshin.verts.isearch(surfacePoints[1])->coord,
					meshin.verts.isearch(surfacePoints[2])->coord,
					testCoord,temp1, temp2);
				if(fabs(distToPlane)<__DBL_EPSILON__){
					continue;
				}
				pointInVolu = meshin.surfs.isearch(surfCurr)->
					voluind[int((distToPlane*multiplierDist)>0)];
				if(pointInVolu==0){
					volumeCandidate=pointInVolu;
					break;
				} else if(pointInVolu!=volumeCandidate){
					isInCandidate=false;
					volumeCandidate=pointInVolu;
					break;
				}
			}
			if(!candidateIsValid && isInCandidate){
				isInCandidate=false;
				volumeCandidate = alternateVolume;
			}
			nVolusExplored++;
		} while (!isInCandidate && nVolusExplored<nVolus);
		if (!isInCandidate){
			RSVS3D_ERROR_LOGIC("Algorithm did not find the cell in which the vertex lies.");
		}
		#ifdef RSVS_DIAGNOSTIC_RESOLVED
		cout << nVolusExplored << " ";
		#endif
		return volumeCandidate;
	}

	double VerticesDistanceSquared(const mesh &meshin,
		const vector<int> &vertind){
		double lengthSquared=0.0;

		for (int i = 0; i < meshin.WhatDim(); ++i)
		{
			lengthSquared += pow(
					meshin.verts.isearch(vertind[0])->coord[i]
					-meshin.verts.isearch(vertind[1])->coord[i]
				,2.0);
		}

		return lengthSquared;
	}
	double VerticesDistance(const mesh &meshin,
		const vector<int> &vertind){
		return sqrt(VerticesDistanceSquared(meshin, vertind));
	}
	bool IsVerticesDistance0(const mesh &meshin,
		const vector<int> &vertind, double eps){ 
		return VerticesDistanceSquared(meshin,vertind)< pow(eps,2.0);
	}
}
