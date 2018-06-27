
#include "RSVSmath.hpp"
#include "RSVSmath_automatic.hpp"


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

/// CoordFunc supports the same stuff as tri func but can have any number of points

bool CoordFunc::MakeValidField(vector<double>* mp){
	// Tries to make a valid array
	// Warning : This operates directly on the data stored in the field of the pointer
	
	int ii,n;
	bool fieldReady; 
	if ((mp)==NULL){
		fieldReady=false;
	} else {
		n=int((mp)->size());
		if (n<nDim){
			for(ii=0;ii<(nDim-n);++ii){
				(mp)->push_back(0);
			}
		} else if (n>nDim){
			for(ii=0;ii<(n-nDim);++ii){
				(mp)->pop_back();
			}
		} 
	}
	return(fieldReady);
}

bool CoordFunc::CheckValid(){
	isReady=true;

	for (int ii=0; ii<nCoord; ++ii){
		if (coords[ii]==NULL){
			isReady=false;
		} else if (int(coords[ii]->size())!=nDim){
			isReady=false;
		}
	}
	return(isReady);
} 

bool CoordFunc::MakeValid(){
	// Tries to make a valid array
	// Warning : This operates directly on the data of the pointer

	isReady=true;
	bool fieldReady;
	for (int ii=0; ii<nCoord; ++ii){
		fieldReady=MakeValidField(coords[ii]);isReady=isReady & fieldReady;
	}

	return(isReady);
}

void CoordFunc::assign(vector<vector<double>*> &pRep){
	if (int(pRep.size())==nCoord){
		coords=pRep;
	} else {
		ResetNCoord(int(pRep.size()));
		coords=pRep;
	}
	isCalc=false;
	isReady=false;
}

void CoordFunc::assign(int pRepI,vector<double> &pRep){
	 
	if ((pRepI<nCoord) & (pRepI>=0)){
		coords[pRepI]=&pRep;
	} else {
		throw invalid_argument("pointer index out off bound in CoordFunc::assign");
	}
	
	isCalc=false;
	isReady=false;
}

void CoordFunc::PreCalc(){
	
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

void CoordFunc::InitialiseArrays(){
	fun=0;
	if (nFun>1){
		funA.assign(nFun,1,fun); 
	}
	jac.assign(nFun,nCoord*nDim,fun); 
	hes.assign(nCoord*nDim,nCoord*nDim,fun);
	coords.assign(nCoord,NULL);
	isReady=false;
	isCalc=false;
}

// void CoordFunc::ResetDim(int nDim){ // implicitely inlined
// 	InitialiseArrays();
// }
// void CoordFunc::ResetNCoord(int nCoord){
// 	InitialiseArrays();
// }

 
// Derived classes

void Volume::Calc(){
	#ifdef SAFE_ALGO
	PreCalc(); 
	#endif

	if(!isCalc){
		Volume_f(*p0,*p1,*p2,fun);
		Volume_df(*p0,*p1,*p2,jac);
		Volume_ddf(*p0,*p1,*p2,hes);
	}
}
 
void Area::Calc(){
	#ifdef SAFE_ALGO
	PreCalc();
	#endif

	if(!isCalc){
		Area_f(*p0,*p1,*p2,fun);
		Area_df(*p0,*p1,*p2,jac);
		Area_ddf(*p0,*p1,*p2,hes);
	}
}
 
void LengthEdge::Calc(){  
	#ifdef SAFE_ALGO
	PreCalc();
	#endif
	if(!isCalc){
		LengthEdge_f(*coords[0],*coords[1],fun);
		LengthEdge_df(*coords[0],*coords[1],jac);
		LengthEdge_ddf(*coords[0],*coords[1],hes);
	}
}

void SurfCentroid::Calc(){

	vector<double> x,y,z;
	int ii,n;

	if (!isCalc){
		if  (nCoord<4){
			throw invalid_argument("nCoord <= 3 surfCentroid Not implemented");
		} else if (nCoord<7){
			x.reserve(nCoord);
			y.reserve(nCoord);
			z.reserve(nCoord);
			// Build Full x,y,z vectors
			for (ii=0; ii<nCoord; ++ii){
				x.push_back((*coords[ii])[0]);
				y.push_back((*coords[ii])[1]);
				z.push_back((*coords[ii])[2]);
			}

			switch(nCoord){
				case 4:
					SurfCentroid4_f(x,y,z,funA);
					SurfCentroid4_df(x,y,z,jac);
				case 5:
					SurfCentroid5_f(x,y,z,funA);
					SurfCentroid5_df(x,y,z,jac);
				case 6:
					SurfCentroid6_f(x,y,z,funA);
					SurfCentroid6_df(x,y,z,jac);
			}
		} else {
			x.reserve(6);
			y.reserve(6);
			z.reserve(6);

		}
	}
}