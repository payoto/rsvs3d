
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <functional>
#include <cmath>

#include "mesh.hpp"


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

vector<double> coordvec::cross(const std::vector<double> &vecin){

	vector<double> retVec;
	retVec.assign(3,0.0);

	for(int ii=3; ii<6; ii++){
		retVec[ii%3]=elems[(ii+1)%3]*vecin[(ii+2)%3]
		-elems[(ii-1)%3]*vecin[(ii-2)%3];
	}
	return(retVec);
}

double coordvec::dot(const std::vector<double> &vecin){

	double retVec=0.0;
	
	for(int ii=0; ii<3; ii++){
		retVec+=elems[ii]*vecin[ii];
	}
	return(retVec);
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
		cerr << "Error in " << __PRETTY_FUNCTION__ << endl; 
		throw invalid_argument("parent and child index list must be the same size.");
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


void mesh::VoluValuesofParents(int elmInd, vector<double> &vals, int volType) const{
		/*
	Extracts the volumes in the parents corresponding to element elmInd
	volType is a selector for which value to extract
	*/
	
	double volu::*mp;

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

int mesh::CountVoluParent() const {
	int n=0;
	for (int i = 0; i < int(meshtree.parentmesh.size()); ++i)
	{
		n+=meshtree.parentmesh[i]->volus.size();
	}
	return(n);
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

void mesh::SwitchIndex(int typeInd, int oldInd, int newInd, vector<int> scopeInd)
{

	int ii,jj,kk,newSub,oldSub;
	vector<int> subList;
	HashedVector<int,int> tempSub;
	bool is3DMesh=volus.size()>0;

	if(typeInd==1){
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
   	if(meshDepIsSet && meshDim==3){
	  	// for (ii=0;ii<int(oldSub.size());++ii){
	  	// 	meshtree.elemind[oldSub[ii]]=newInd;
	  	// }
   		meshtree.elemind[oldSub]=newInd;
   	}
   	edges.isHash=1;
   	surfs.isHash=1;
   	edges.isSetMI=1;
   	surfs.isSetMI=1;
   } else {

   	cerr << "Error unknown type " << typeInd << " of object for index switching" <<endl;
   	cerr << "Error in " << __PRETTY_FUNCTION__ << endl;
   	throw invalid_argument (" Type is out of range");
   }

}

void mesh::RemoveIndex(int typeInd, int oldInd)
{

	int ii,jj;
	vector<int> subList;

	if(typeInd==1){

		cerr << "not coded yet" << endl;
		throw;

	} else if (typeInd==2){


		subList=verts.find_list(edges(edges.find(oldInd))->vertind);
		for (ii=0;ii<int(subList.size());++ii){
			for (jj=0;jj<int(verts(subList[ii])->edgeind.size());++jj){
				if(verts(subList[ii])->edgeind[jj]==oldInd){
					verts[subList[ii]].edgeind.erase(
						verts[subList[ii]].edgeind.begin()+jj);
					jj--;
				}
			}  
		}

		subList=surfs.find_list(edges(edges.find(oldInd))->surfind);
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

		edges.isHash=1;
		verts.isHash=1;
		surfs.isHash=1;
		edges.isSetMI=1;
		verts.isSetMI=1;
		surfs.isSetMI=1;


	} else if (typeInd==3){


		subList=edges.find_list(surfs(surfs.find(oldInd))->edgeind);
		for (ii=0;ii<int(subList.size());++ii){
			for (jj=0;jj<int(edges(subList[ii])->surfind.size());++jj){
				if(edges(subList[ii])->surfind[jj]==oldInd){
					edges[subList[ii]].surfind.erase(
						edges[subList[ii]].surfind.begin()+jj);
					jj--;
				}
			}  
		}

		subList=volus.find_list(surfs(surfs.find(oldInd))->voluind);
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

		surfs.isHash=1;
		edges.isHash=1;
		volus.isHash=1;
		surfs.isSetMI=1;
		edges.isSetMI=1;
		volus.isSetMI=1;

	} else if (typeInd==4){

		cerr << "not coded yet" << endl;
		throw;
   } else if (typeInd==5){ // Modify vertex index in scoped mode

   	cerr << "not coded yet" << endl;
   	throw;
   } else {

   	cerr << "Error unknown type of object for index switching" <<endl;
   	cerr << "Error in " << __PRETTY_FUNCTION__ << endl;
   	throw invalid_argument (" : Type is out of range");
   }

}

void mesh::TestConnectivity(){
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
		cerr << errTot << "  Total errors were detected in the connectivity list" <<endl;
	}
}


void mesh::TestConnectivityBiDir(){
	int ii,jj,ll,ll2,kk,kk2,errCount, errTot,errCountBiDir,errTotBiDir;
	bool flag;
	vector<int> testSub;

	errCount=0;
	errCountBiDir=0;
	errTot=0;
	errTotBiDir=0;
	kk=int(verts.size());
	for (ii=0; ii<kk;++ii){
		if(verts(ii)->edgeind.size()==0){
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
					cerr << " Test Connectivity Error :" << errCountBiDir << " surf " << volus(ii)->index
					<< " makes uni-directional reference to volume " << volus(ii)->surfind[jj] << " list: " ; 
					DisplayVector(volus(ii)->surfind); cout << " list (surfs.voluind): " ; 
					DisplayVector(surfs(jj)->voluind);cerr << endl;
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
		cerr << errTot << "  Total errors were detected in the connectivity list" <<endl;
	}
	if (errTotBiDir>0){
		cerr << errTotBiDir << "  Total errors were detected in the bi-directionality of the connectivity list" << endl;
	}
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
	this->GetMaxIndex(&nVert,&nEdge,&nVolu,&nSurf);
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
void surf::OrderEdges(mesh *meshin)
{
	unordered_multimap<int,int> vert2Edge;
	vector<int> edge2Vert,edgeSub,edgeIndOrig;
	vector<bool> isDone;
	bool isTruncated;
	int vertCurr,edgeCurr;
	int ii,jj;
	std::pair<unordered_multimap<int,int>::iterator,unordered_multimap<int,int>::iterator> range;
	unordered_multimap<int,int>::iterator it;
	isTruncated=false;
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
			range=vert2Edge.equal_range(vertCurr);
		 	#ifdef SAFE_ACCESS
			if (range.first==vert2Edge.end()){
				disptree((*meshin),0);

				cerr << ii << " vert " << vertCurr << "  ";
				DisplayVector(edge2Vert);
				DisplayVector(edgeind);
				meshin->verts.isearch(vertCurr)->disp();
				cout << it->second << " " << 1/2 << 2/3 <<  endl;
				cerr << "Error in :" << __PRETTY_FUNCTION__ << endl;
				throw range_error ("unordered_multimap went beyond its range in OrderEdges");
			}
			if (vert2Edge.count(vertCurr)==1){
				cerr << "ERROR : Surface does not form closed loop" << endl;
				disptree((*meshin),0);

				jj=edgeIndOrig[(range.first->second)/2]==edgeCurr;
				DisplayVector(edge2Vert);
				DisplayVector(edgeind);
				cerr <<endl;

				meshin->verts.isearch(vertCurr)->disp();
				cerr << "ii is : " << ii << " jj is : " << jj << " count is : " ;
				jj=vert2Edge.count(vertCurr);
				cerr << jj  <<  endl; 
				cerr << "Error in :" << __PRETTY_FUNCTION__ << endl;
				//throw range_error ("Found a single vertex - surface is not closed");
				cerr << "Found a single vertex - surface is not closed" << endl;
				return;
			}
		 	#endif // SAFe_ACCESS
			jj=edgeIndOrig[(range.first->second)/2]==edgeCurr;

			it=range.first;
			if (jj){++it;}
			if (!isDone[(it->second)/2]){
				isDone[(it->second)/2]=true;
				edgeCurr=edgeIndOrig[(it->second)/2];
				jj=edge2Vert[((it->second)/2)*2]==vertCurr; // Warnign ((x/2)*2) not necessarily equal x
				vertCurr=edge2Vert[((it->second)/2)*2+jj];
				edgeind[ii]=edgeCurr;
			} else {
				cerr << endl;
				// edgeind = edgeIndOrig;
				edgeind.erase(edgeind.begin()+ii, edgeind.end());
				// DisplayVector(edgeIndOrig);
				disptree((*meshin),1);

				// cerr << "Error in :" << __PRETTY_FUNCTION__ << endl;

				isTruncated=true;
				break;
			}
		}
		if (isTruncated){ 
			// This adds a second surface if the surface closes early

			surf newSurf=*this;
			meshin->surfs.SetMaxIndex();
			newSurf.index=meshin->surfs.GetMaxIndex()+1;
			newSurf.voluind=voluind;
			newSurf.edgeind.clear();
			for (ii=0;ii<int(edgeIndOrig.size());++ii){
				if(!isDone[ii]){
					newSurf.edgeind.push_back(edgeIndOrig[ii]);
				}
			}
			newSurf.OrderEdges(meshin);
			meshin->surfs.push_back(newSurf);
			meshin->SwitchIndex(6,index,newSurf.index,newSurf.edgeind);
			meshin->surfs(meshin->surfs.size()-1)->disptree(*meshin,1);

		}
		isordered=true;
	}
}
/*void surf::OrderEdges(mesh *meshin)
{
	unordered_multimap<int,int> vert2Edge;
	vector<int> edge2Vert,edgeSub,edgeIndOrig,edgeIndOrig2,edgeFinalShrunk;
	vector<bool> isDone;
	bool isTruncated;
	int vertCurr,edgeCurr;
	int ii,jj;
	std::pair<unordered_multimap<int,int>::iterator,unordered_multimap<int,int>::iterator> range;
	unordered_multimap<int,int>::iterator it;
	isTruncated=false;
	if (edgeind.size()>0){
		edgeIndOrig2=edgeind;
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
			range=vert2Edge.equal_range(vertCurr);
		 #ifdef SAFE_ACCESS
			if (range.first==vert2Edge.end()){
				cerr << ii << " vert " << vertCurr << "  ";
				DisplayVector(edge2Vert);
				DisplayVector(edgeind);
				meshin->verts.isearch(vertCurr)->disp();
				cout << it->second << " " << 1/2 << 2/3 <<  endl;
				cerr << "Error in :" << __PRETTY_FUNCTION__ << endl;
				cout<< "throw range_error (" <<"unordered_multimap went beyond its range in OrderEdges" <<")";
			}
			if (vert2Edge.count(vertCurr)==1){
				jj=edgeIndOrig[(range.first->second)/2]==edgeCurr;
				DisplayVector(edge2Vert);
				DisplayVector(edgeind);
				cerr <<endl;

				meshin->verts.isearch(vertCurr)->disp();
				cerr << "Error : Surface does not form closed loop" << endl;
				cerr << "ii is : " << ii << " jj is : " << jj << " count is : " ;
				jj=vert2Edge.count(vertCurr);
				cerr << jj  <<  endl; 
				cerr << "Error in :" << __PRETTY_FUNCTION__ << endl;
				cout << "throw range_error ("<<"Found a single vertex - surface is not closed"<<")";
				cout << endl;
				disp();
			}
		 #endif // SAFe_ACCESS
			jj=edgeIndOrig[(range.first->second)/2]==edgeCurr;

			it=range.first;
			if (jj){++it;}
			if (!isDone[(it->second)/2]){
				isDone[(it->second)/2]=true;
				edgeCurr=edgeIndOrig[(it->second)/2];
			jj=edge2Vert[((it->second)/2)*2]==vertCurr; // Warnign ((x/2)*2) not necessarily equal x
			vertCurr=edge2Vert[((it->second)/2)*2+jj];
			edgeind[ii]=edgeCurr;
		} else {
			edgeind.erase(edgeind.begin()+ii);
			isTruncated=true;
			break;
		}
	}
	if (isTruncated){

		surf newSurf=*this;
		newSurf.index=meshin->surfs.GetMaxIndex()+1;
		newSurf.voluind=voluind;
		newSurf.edgeind.clear();
		for (ii=0;ii<int(edgeIndOrig.size());++ii){
			if(!isDone[ii]){
				newSurf.edgeind.push_back(edgeIndOrig[ii]);
			}
		}
		meshin->surfs.push_back(newSurf);
		meshin->SwitchIndex(6,index,newSurf.index,newSurf.edgeind);
		meshin->surfs(meshin->surfs.size()-1)->disptree(*meshin,1);
	}
	  #ifdef SAFE_ALGO
	edgeFinalShrunk=edgeind;
	sort(edgeFinalShrunk);
	unique(edgeFinalShrunk);
	if(edgeFinalShrunk.size()!=edgeIndOrig.size()){
		cout << endl;
		cout << "edgeind (pre)- ";
		DisplayVector(edgeIndOrig2);
		cout << endl;
		cout << "edgeind (sort)- ";
		DisplayVector(edgeIndOrig);
		cout << endl;
		cout << "edgeind (final)- ";
		DisplayVector(edgeind);
		cout << endl;

		 //disptree((*meshin),2);
		cerr << "Warning: OrderEdges removed some edges from the surface list "<< endl;
		cerr << " The surface had incorrect connectivity, it was self crossing "<< endl;
		cerr << "   in function:" <<  __PRETTY_FUNCTION__ << endl;
		 //throw invalid_argument ("Incorrect surface connectivity");
	}
	  #endif //SAFE_ALGO
	isordered=true;
}
}*/


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
			throw invalid_argument("Surface not ordered");
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

int mesh::OrderEdges(){
	int ii;
	bool kk;

	for (ii = 0; ii < surfs.size(); ++ii)
	{
		if (!surfs(ii)->isready(true)){
			surfs.elems[ii].OrderEdges(this);
			kk=true;
		}
	}
	return(kk);
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
	int ii, ni, jj, nj,kk,nk,ll,nl, vertTemp, edgeSub,surfSub;
	vector<double> vals;
	bool surfCond;

	ni=volus.size();
	nj=verts.size();
	vertInd.clear();
	vertInd.reserve(nj);
	voluInd.clear();
	voluInd.reserve(ni);

	for (ii=0; ii<ni; ++ii){
		surfCond=volus(ii)->isBorder;
		if(surfCond){
			if(outerVolume==0){
				this->VoluValuesofParents(volus(ii)->index,vals, &volu::target);
				surfCond=false;
				jj=0;
				while (!surfCond && jj<meshtree.nParents){
					surfCond=vals[jj]==1.0;
					++jj;
				}
			} else if(outerVolume==1){
				this->VoluValuesofParents(volus(ii)->index,vals, &volu::target);
				surfCond=false;
				jj=0;
				while (!surfCond && jj<meshtree.nParents){
					surfCond=vals[jj]>0.0;
					++jj;
				}
			} else if(outerVolume!=-1){
				throw invalid_argument("Unkownn value of outerVolume -1,0, or 1");
			}
		} 
		// THIS IS DODGY HOW to pick the side to delete can fail
		if(surfCond){
			// Pick one of the volumes to delete
			if(outerVolume==0){
				nj = volus(ii)->surfind.size();
				for(jj=0;jj<nj;++jj){
					surfSub = surfs.find(volus(ii)->surfind[jj]);
					for(kk=0;kk<2;++kk){
						if(surfs(surfSub)->voluind[kk]!=0){
							if(!(volus.isearch(surfs(surfSub)->voluind[kk])->isBorder)){
								voluInd.push_back(surfs(surfSub)->voluind[kk]);
							}
						}
					}
				}
			} else {
				voluInd.push_back(volus(ii)->index);
			}
			nl=volus(ii)->surfind.size();
			for(ll=0; ll<nl; ++ll){
				surfSub=surfs.find(volus(ii)->surfind[ll]);
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
			nT=surfs(ii)->voluind.size();
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


	nBlocks=this->ConnectedVertex(vertBlock);

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
				nSurf=edges(iEdge)->surfind.size();
				for (kk=0;kk<nSurf;++kk){
					iSurf=surfs.find(edges(iEdge)->surfind[kk]);
					volus.elems[volus.find(vertBlock[ii])].surfind.push_back(edges(iEdge)->surfind[kk]);
					surfs.elems[iSurf].voluind.clear();
					surfs.elems[iSurf].voluind.push_back(vertBlock[ii]);
					surfs.elems[iSurf].voluind.push_back(0);
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
   vertBlock.assign(nVerts,0);
   currQueue.reserve(nVerts/2);
   nextQueue.reserve(nVerts/2);

   
   // While Loop, while not all vertices explored
   while(nVertExplored<nVerts){

	  // if currQueue is empty start new block
   	if(currQueue.size()<1){

		 //cout << "Block " << nBlocks << " - " << nVertExplored << " - " << nVerts << endl;
   		ii=0;
   		while(vertStatus[ii] && ii<nVerts){
   			ii++;
   		}
   		if (vertStatus[ii]){
   			cerr << "Error starting point for loop not found despite max number of vertex not reached" <<endl;
   			cerr << "Error in " << __PRETTY_FUNCTION__ << endl;
   			throw range_error (" : Starting point for block not found");
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
   				throw range_error (" : Vertex not found");
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
				throw range_error (" Start point not found");
			}
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
		throw invalid_argument("Incaompatible list and index");
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
		throw invalid_argument("Incaompatible list and index");
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

