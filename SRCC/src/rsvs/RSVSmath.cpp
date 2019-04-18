#include <iostream>
#include "RSVSmath.hpp"
#include "RSVSmath_automatic.hpp"
#include "arraystructures.hpp" // for use of DisplayVector
#include "warning.hpp"

using namespace std;


bool TriFunc::MakeValidField(vector<double>* TriFunc::*mp){
	// Tries to make a valid array
	// Warning : This operates directly on the data stored in the field of the pointer
	 
	int ii,n;
	bool fieldReady=true; 
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
	//bool fieldReady;
	CheckValid();
	// fieldReady=MakeValidField(&TriFunc::p0);isReady=isReady & fieldReady;
	// fieldReady=MakeValidField(&TriFunc::p1);isReady=isReady & fieldReady;
	// fieldReady=MakeValidField(&TriFunc::p2);isReady=isReady & fieldReady;

	return(isReady);
}

void TriFunc::assign(const vector<double> &in0,const vector<double> &in1,const vector<double> &in2){
	p0=&in0;
	p1=&in1;
	p2=&in2;

	isCalc=false;
	isReady=false;
}
void TriFunc::assign(const vector<double> *in0,const vector<double> *in1,const vector<double> *in2){
	p0=in0;
	p1=in1;
	p2=in2;

	isCalc=false;
	isReady=false;
}
void TriFunc::assign(int pRepI,const vector<double> &pRep){
	
	switch(pRepI){
		case 0:
		p0=&pRep;
		break;
		case 1: 
		p1=&pRep;
		break;
		case 2:
		p2=&pRep;
		break;
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


void TriFunc::ReturnDatPoint(double **a, ArrayVec<double> **b,ArrayVec<double> **c) {
	*a=&fun;
	*b=&jac;
	*c=&hes; 
}


/// CoordFunc supports the same stuff as tri func but can have any number of points

bool CoordFunc::MakeValidField(vector<double> const * mp){
	// Tries to make a valid array
	// Warning : This operates directly on the data stored in the field of the pointer
	
	// int ii,n;
	bool fieldReady; 
	if ((mp)==NULL){
		fieldReady=false;
	} else {
		fieldReady=false;
		// n=int((mp)->size());
		// if (n<nDim){
		// 	for(ii=0;ii<(nDim-n);++ii){
		// 		(mp)->push_back(0);
		// 	}
		// } else if (n>nDim){
		// 	for(ii=0;ii<(n-nDim);++ii){
		// 		(mp)->pop_back();
		// 	}
		// } 
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

void CoordFunc::assign(vector<vector<double> const*> &pRep){
	if (int(pRep.size())==nCoord){
		coords=pRep;
	} else {
		ResetNCoord(int(pRep.size()));
		coords=pRep;
	}
	isCalc=false; 
	isReady=false;
}

void CoordFunc::assign(int pRepI, const vector<double> &pRep){
	 
	if ((pRepI<nCoord) & (pRepI>=0)){
		coords[pRepI]=&pRep;
	} else {
		RSVS3D_ERROR_ARGUMENT("pointer index out off bound in CoordFunc::assign");
	}
	
	isCalc=false;
	isReady=false;
}


void CoordFunc::ReturnDat(double &a, ArrayVec<double> &b,ArrayVec<double> &c){
	a=fun;
	b=jac;
	c=hes; 
}
void CoordFunc::ReturnDat(ArrayVec<double> &a, ArrayVec<double> &b,ArrayVec<double> &c){
	a=funA;
	b=jac;
	c=hes;
}
void CoordFunc::ReturnDatPoint(double **a, ArrayVec<double> **b,ArrayVec<double> **c){
	*a=&fun;
	*b=&jac;
	*c=&hes; 
}
void CoordFunc::ReturnDatPoint(ArrayVec<double> **a, ArrayVec<double> **b,ArrayVec<double> **c){
	*a=&funA;
	*b=&jac;
	*c=&hes;
}

void CoordFunc::PreCalc(){
	
	if (!isReady){
		CheckValid();
	}
	if (!isReady){
		MakeValid();
	}
	if (!isReady){
		cerr << "CoordFunc cannot be made to match the required paremeters " << endl;
	}
	
}

void CoordFunc::InitialiseArrays(){
	fun=0;
	if (nFun>1){
		funA.assign(nFun,1,fun); 
	}
	jac.assign(nFun,nCoord*nDim,fun); 
	hes.assign(nCoord*nDim,nCoord*nDim*nFun,fun);
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
		/* This function calculates Volume Jacobian and Hessian using
matlab automatically generated code

the jacobian is arranged :
	  x0 y0 z0  x1 y1 z1  x2 y2 z2
 V [                               ] 

 The Hessian is arranged:
 					V 	
	   x0 y0 z0  x1 y1 z1  x2 y2 z2
 x0 [                               ]  
 y0 [                               ] 
 z0 [                               ] 
 x1 [                               ] 
 y1 [                               ] 
 z1 [                               ] 
 x2 [                               ] 
 y2 [                               ] 
 z2 [                               ] 

*/
	#ifdef SAFE_ALGO
	PreCalc(); 
	#endif

	if(!isCalc){
		Volume_f(*p0,*p1,*p2,fun);
		Volume_df(*p0,*p1,*p2,jac);
		Volume_ddf(*p0,*p1,*p2,hes);
	}
}

void Volume::CalcFD(){
		/* This function calculates Volume Jacobian and Hessian using
matlab automatically generated code

the jacobian is arranged :
	  x0 y0 z0  x1 y1 z1  x2 y2 z2
 V [                               ] 

 The Hessian is arranged:
 					V 	
	   x0 y0 z0  x1 y1 z1  x2 y2 z2
 x0 [                               ]  
 y0 [                               ] 
 z0 [                               ] 
 x1 [                               ] 
 y1 [                               ] 
 z1 [                               ] 
 x2 [                               ] 
 y2 [                               ] 
 z2 [                               ] 

*/  
	#ifdef SAFE_ALGO
	PreCalc(); 
	#endif
	std::vector<double> v;
	double fdStep=1e-6;
	double tempVal;
	
	vector<vector<double>> vecs;

	v.reserve(3);

	if(!isCalc){
		this->Calc();

		vecs.clear();
		vecs.reserve(3);
		vecs.push_back(*p0);
		vecs.push_back(*p1);
		vecs.push_back(*p2);
		for(int ii=0; ii<3; ++ii){
			for(int jj=0; jj< 3; ++jj){
				vecs[ii][jj] = vecs[ii][jj] + fdStep;
				Volume_f(vecs[0],vecs[1],vecs[2],tempVal);
				jac[0][ii*3+jj] = (tempVal-this->fun)/fdStep;
				vecs[ii][jj] = vecs[ii][jj] - fdStep;
			}

		}

	}
}
 
void Area::Calc(){
		/* This function calculates Area Jacobian and Hessian using
matlab automatically generated code

the jacobian is arranged :
	  x0 y0 z0  x1 y1 z1  x2 y2 z2
 A [                               ] 

 The Hessian is arranged:
 					A 	
	   x0 y0 z0  x1 y1 z1  x2 y2 z2
 x0 [                               ]  
 y0 [                               ] 
 z0 [                               ] 
 x1 [                               ] 
 y1 [                               ] 
 z1 [                               ] 
 x2 [                               ] 
 y2 [                               ] 
 z2 [                               ] 

*/
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


// void SurfCentroid::Calcddf(){
	
// }

void SurfCentroid::CalcFD(){
	
	this->Calc();
	auto baseVal = this->funA;
	std::vector<double> vtemp;
	double fdStep=1e-6;

	this->calcDeriv = false;
	this->jac.assign(3,nCoord*3, 0.0);

	for (int i = 0; i < nCoord; ++i)
	{
		auto saveVec = this->coords[i];
		this->coords[i] = &vtemp;
		for (int j = 0; j < 3; ++j)
		{	
			this->isCalc = false;
			vtemp = *saveVec;
			vtemp[j] = vtemp[j]+fdStep;
			this->Calc();
			for (int k = 0; k < 3; ++k)
			{
				this->jac[k][i+nCoord*j] = (this->funA[k][0] - baseVal[k][0])/fdStep;
			}
		}
		this->coords[i] = saveVec;

	}

}

void SurfCentroid::Calc(){
/* This function calculates centroid Jacobian and Hessian using
matlab automatically generated code

the jacobian is arranged :
	  x0 x1  xn y0 y1  yn  z0 z1  zn
 xc [                               ] 
 yc [                               ]
 zc [                               ]

 The Hessian is arranged:
 					xc 				 					yc 				 				zc 				
	  x0 x1  xn y0 y1  yn  z0 z1  zn	x0 x1  xn y0 y1  yn  z0 z1  zn	    x0 x1  xn y0 y1  yn  z0 z1  zn
 x0 [                               ]  [                               ]  [                               ] 
 x1 [                               ]  [                               ]  [                               ] 
 xn [                               ]  [                               ]  [                               ] 
 y0 [                               ]  [                               ]  [                               ] 
 y1 [                               ]  [                               ]  [                               ] 
 yn [                               ]  [                               ]  [                               ] 
 z0 [                               ]  [                               ]  [                               ] 
 z1 [                               ]  [                               ]  [                               ] 
 zn [                               ]  [                               ]  [                               ] 

*/

	vector<double> x,y,z,centroidTLength;
	ArrayVec<double> temp;
	double currEdge;
	int ii,ii1,ii2,jj,jj2,nR,nC,ind,ind2;

	this->PreCalc();

	
	if (!isCalc){
		// Calculate centroid and 
		centroid.assign(nDim,0);
		edgeLength=0;
		for(jj=0;jj<nCoord;++jj){
			currEdge=0;
			LengthEdge_f(*coords[jj],*coords[(jj+1)%nCoord],currEdge);
			edgeLength+=currEdge;
			for(ii=0; ii<nDim; ++ii){
				centroid[ii]+=currEdge*(
					(*(coords[jj]))[ii]
					+(*(coords[(jj+1)%nCoord]))[ii])/2;
			}
			
		}
		if(edgeLength<0.000000000000001){
			cout << "Warning edgeLength is 0" << endl;
			for(ii=0; ii<nDim; ++ii){
				funA[ii][0]=(*(coords[0]))[ii];
			}
			edgeLength=0.0000001;
		} else {
			for(ii=0; ii<nDim; ++ii){
				funA[ii][0]=centroid[ii]/edgeLength;
			}
		}
		if(!this->calcDeriv){
			return;
		}
		if  (nCoord<4){
			RSVS3D_ERROR_ARGUMENT("nCoord <= 3 surfCentroid Not implemented");
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
					//SurfCentroid4_f(x,y,z,edgeLength,centroid[0],centroid[1],centroid[2],funA);
					SurfCentroid4_df(x,y,z,edgeLength,centroid[0],centroid[1],
						centroid[2],jac);
					break;
				case 5:
					//SurfCentroid5_f(x,y,z,edgeLength,centroid[0],centroid[1],centroid[2],funA);
					SurfCentroid5_df(x,y,z,edgeLength,centroid[0],centroid[1],
						centroid[2],jac);
					break;
				case 6:
					//SurfCentroid6_f(x,y,z,edgeLength,centroid[0],centroid[1],centroid[2],funA);
					SurfCentroid6_df(x,y,z,edgeLength,centroid[0],centroid[1],
						centroid[2],jac);
					break;
			}
		} else {
			/*the jacobian is arranged :
			  x0 x1  xn y0 y1  yn  z0 z1  zn
			 xc [                               ] 
			 yc [                               ]
			 zc [                               ]*/
			x.assign(3,0);
			y.assign(3,0);
			z.assign(3,0);
			temp.assign(nDim,nDim,0); 
			for (jj=0;jj<nCoord;jj++){
				for (ii=0; ii<3 ; ++ii){ 
					x[ii]=((*coords[(jj+ii)%nCoord])[0]);
					y[ii]=((*coords[(jj+ii)%nCoord])[1]);
					z[ii]=((*coords[(jj+ii)%nCoord])[2]);
				}
				SurfCentroidSelf_df(x,y,z,edgeLength,centroid[0],centroid[1],centroid[2],temp);
				ind=(jj+1)%nCoord;
				// WARNING : This needs to be checked to make sure it matches the automaticly generated ones
				for(ii1=0; ii1<nDim ; ii1++){
					for(ii2=0; ii2<nDim ; ii2++){
						jac[ii1][ind+(ii2*nCoord)]=temp[ii1][ii2];
					}
				}
			}
		}

		// Calculate Hessian

		/* The Hessian is arranged:
 					xc 				 					yc 				 				zc 				
	  x0 x1  xn y0 y1  yn  z0 z1  zn	x0 x1  xn y0 y1  yn  z0 z1  zn	    x0 x1  xn y0 y1  yn  z0 z1  zn
 x0 [                               ]  [                               ]  [                               ] 
 x1 [                               ]  [                               ]  [                               ] 
 xn [                               ]  [                               ]  [                               ] 
 y0 [                               ]  [                               ]  [                               ] 
 y1 [                               ]  [                               ]  [                               ] 
 yn [                               ]  [                               ]  [                               ] 
 z0 [                               ]  [                               ]  [                               ] 
 z1 [                               ]  [                               ]  [                               ] 
 zn [                               ]  [                               ]  [                               ] */

		// Squared terms of the Hessian
		hes.size(nR,nC);
		x.assign(3,0);
		y.assign(3,0);
		z.assign(3,0);
		temp.assign(nDim,nDim*nFun,0); 
		for (jj=0;jj<nCoord;jj++){
			for (ii=0; ii<3 ; ++ii){ 
				x[ii]=((*coords[(jj+ii)%nCoord])[0]);
				y[ii]=((*coords[(jj+ii)%nCoord])[1]);
				z[ii]=((*coords[(jj+ii)%nCoord])[2]);
			}
			SurfCentroidSelf_ddf(x,y,z,edgeLength,centroid[0],centroid[1],centroid[2],temp);
			ind=(jj+1)%nCoord;
			// WARNING : This needs to be checked to make sure it matches the automaticly generated ones
			for(ii1=0; ii1<nDim ; ii1++){
				for(ii2=0; ii2<nDim ; ii2++){
					for(ii=0;ii<nFun;ii++){
						hes[ind+(ii1*nCoord)][ind+(ii2*nCoord)+(ii*nDim*nCoord)]=temp[ii1][ii2+(nDim*ii)];
					}
				}
			}
		}

		// Connected terms of the Hessian
		x.assign(6,0);
		y.assign(6,0);
		z.assign(6,0);
		temp.assign(nDim*2,2*nDim*nFun,0); 
		for (jj=0;jj<nCoord;jj++){
			jj2=(jj+1)%nCoord;
			for (ii=0; ii<4 ; ++ii){ 
				x[ii]=((*coords[(jj+ii)%nCoord])[0]);
				y[ii]=((*coords[(jj+ii)%nCoord])[1]);
				z[ii]=((*coords[(jj+ii)%nCoord])[2]);
			}
			SurfCentroidConnec_ddf(x,y,z,edgeLength,centroid[0],centroid[1],centroid[2],temp);
			ind=(jj+1)%nCoord;
			ind2=(jj2+1)%nCoord;
				// WARNING : This needs to be checked to make sure it matches the automaticly generated ones
			for(ii1=0; ii1<nDim ; ii1++){
				for(ii2=0; ii2<nDim ; ii2++){
					for(ii=0;ii<nFun;ii++){
						hes[ind+(ii1*nCoord)][ind2+(ii2*nCoord)+(ii*nDim*nCoord)]=temp[ii1][nDim+ii2+(nDim*ii)];
						hes[ind2+(ii1*nCoord)][ind+(ii2*nCoord)+(ii*nDim*nCoord)]=temp[ii1+nDim][ii2+(nDim*ii)];
					}
				}
			}
			
		}


		// Not Connected terms of the Hessian
		x.assign(8,0);
		y.assign(8,0);
		z.assign(8,0);
		temp.assign(nDim*2,2*nDim*nFun,0); 
		for (jj=0;jj<nCoord;jj++){
			for(jj2=jj+2;jj2<nCoord; ++jj2){
				for (ii=0; ii<3 ; ++ii){ 
					x[ii]=((*coords[(jj+ii)%nCoord])[0]);
					y[ii]=((*coords[(jj+ii)%nCoord])[1]);
					z[ii]=((*coords[(jj+ii)%nCoord])[2]);
				}

				for (ii=0; ii<3 ; ++ii){ 
					x[ii+4]=((*coords[(jj2+ii)%nCoord])[0]);
					y[ii+4]=((*coords[(jj2+ii)%nCoord])[1]);
					z[ii+4]=((*coords[(jj2+ii)%nCoord])[2]);
				}
				SurfCentroidNoConnec_ddf(x,y,z,edgeLength,centroid[0],centroid[1],centroid[2],temp);
				ind=(jj+1)%nCoord;
				ind2=(jj2+1)%nCoord;
				// WARNING : This needs to be checked to make sure it matches the automaticly generated ones
				for(ii1=0; ii1<nDim ; ii1++){
					for(ii2=0; ii2<nDim ; ii2++){
						for(ii=0;ii<nFun;ii++){
							hes[ind+(ii1*nCoord)][ind2+(ii2*nCoord)+(ii*nDim*nCoord)]=temp[ii1][nDim+ii2+(nDim*ii)];
							hes[ind2+(ii1*nCoord)][ind+(ii2*nCoord)+(ii*nDim*nCoord)]=temp[ii1+nDim][ii2+(nDim*ii)];
						}
					}
				}
			}
		}

	}
}


void SurfCentroid::Disp(){
	 for (int ii=0; ii< nCoord ; ii++){
	 	cout << "c " << ii << " ";
	 	DisplayVector(*coords[ii]); 
	 	cout << endl;
	 }
 	cout << "edgeLength" << edgeLength << endl << "centroid ";
 	DisplayVector(centroid); 
 	cout << endl;
}


void Volume2::Calc(){

	#ifdef SAFE_ALGO
	PreCalc(); 
	#endif
	if(!isCalc){
		Volume2_f((*coords[0])[0],(*coords[0])[1],(*coords[0])[2],
			*coords[1],*coords[2],*coords[3],*coords[4],*coords[5],*coords[6],
			fun);
		Volume2_df((*coords[0])[0],(*coords[0])[1],(*coords[0])[2],
			*coords[1],*coords[2],*coords[3],*coords[4],*coords[5],*coords[6],
			jac);
		Volume2_ddf((*coords[0])[0],(*coords[0])[1],(*coords[0])[2],
			*coords[1],*coords[2],*coords[3],*coords[4],*coords[5],*coords[6],
			hes);
	}
}

