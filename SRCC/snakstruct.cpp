#include <iostream>
#include <cmath>

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

void snaxedge::PrepareForUse() {
	normvector.PrepareForUse();
}

void snaxsurf::PrepareForUse() {
	normvector.PrepareForUse();
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

bool snake::isready() {
	
	readyforuse=true;

	readyforuse=readyforuse & snaxs.isready();
	readyforuse=readyforuse & snaxedges.isready();
	readyforuse=readyforuse & snaxsurfs.isready();
	readyforuse=readyforuse & snakeconn.isready();
	readyforuse=readyforuse & snakemesh->isready();

	return(readyforuse);
}

void snake::PrepareForUse() {
	if (!readyforuse){
		snaxs.PrepareForUse();
		snaxedges.PrepareForUse();
		snaxsurfs.PrepareForUse();
		snakeconn.PrepareForUse();
		snakemesh->PrepareForUse();
		readyforuse=true;
	}
}



// -------------------------------------------------------------------------------------------
// TEST CODE
// -------------------------------------------------------------------------------------------

//in snakstruct_test.cpp