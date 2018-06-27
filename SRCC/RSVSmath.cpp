
#include "RSVSmath.hpp"


using namespace std;


bool TriFunc::MakeValidField(vector<double>* TriFunc::*mp){
	// Tries to make a valid array
	// Warning : This operates directly on the data stored in the field of the pointer
	
	int ii,n;
	bool fieldReady; 
	if ((this->*mp)==NULL){
		fieldReady=false;
	} else {
		n=int((this->*mp)->size());
		if (n<nTarg){
			for(ii=0;ii<(nTarg-n);++ii){
				(this->*mp)->push_back(0);
			}
		} else if (n>nTarg){
			for(ii=0;ii<(n-nTarg);++ii){
				(this->*mp)->pop_back();
			}
		}
	}
	return(fieldReady);
}

bool TriFunc::CheckValid(){
	isReady=true;
	if (p0==NULL){isReady=false;} else if (int(p0->size())!=nTarg){isReady=false;}
	if (p1==NULL){isReady=false;} else if (int(p1->size())!=nTarg){isReady=false;}
	if (p2==NULL){isReady=false;} else if (int(p2->size())!=nTarg){isReady=false;}

	return(isReady);
}

bool TriFunc::MakeValid(){
	// Tries to make a valid array
	// Warning : This operates directly on the data of the pointer

	isReady=true;
	bool fieldReady;

	fieldReady=MakeValidField(&TriFunc::p0);isReady=isReady & fieldReady;
	fieldReady=MakeValidField(&TriFunc::p1);isReady=isReady & fieldReady;
	fieldReady=MakeValidField(&TriFunc::p2);isReady=isReady & fieldReady;

	return(isReady);
}

void TriFunc::assign(vector<double> &in0,vector<double> &in1,vector<double> &in2){
	p0=&in0;
	p1=&in1;
	p2=&in2;

	isCalc=false;
	isReady=false;
}

void TriFunc::assign(int pRepI,vector<double> &pRep){
	
	switch(pRepI){
		case 0:
		p0=&pRep;
		case 1: 
		p1=&pRep;
		case 2:
		p2=&pRep;
	}
	
	isCalc=false;
	isReady=false;
}

void TriFunc::PreCalc(){
	
	if (!isReady){
		CheckValid();
	}
	if (!isReady){
		MakeValid();
	}
	if (!isReady){
		cerr << "TriFunc cannot be made to match the required paremeters " << endl;
	}
	
}

 
void Volume::Calc(){
	#ifdef SAFE_ALGO
	PreCalc(); 
	#endif

	Volume_f(*p0,*p1,*p2,fun);
	Volume_df(*p0,*p1,*p2,jac);
	Volume_ddf(*p0,*p1,*p2,hes);
}
 
void Area::Calc(){
	#ifdef SAFE_ALGO
	PreCalc();
	#endif

	Area_f(*p0,*p1,*p2,fun);
	Area_df(*p0,*p1,*p2,jac);
	Area_ddf(*p0,*p1,*p2,hes);
}
 
void LengthEdge::Calc(){ 
	#ifdef SAFE_ALGO
	PreCalc();
	#endif

	LengthEdge_f(*p0,*p1,fun);
	LengthEdge_df(*p0,*p1,jac);
	LengthEdge_ddf(*p0,*p1,hes);
}