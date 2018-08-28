
#include <cstdlib>
#include <iostream> 
#include "snake.hpp"
#include "snakevel.hpp"
#include "vectorarray.hpp"
#include "RSVSinterface.hpp"
#include <Eigen>

using namespace std; 
using namespace Eigen; 



void SQPcalc::CalculateTriangulation(const triangulation &triRSVS){

	int ii,ni;
	int nDv, nConstr;
	vector<int> vecin;
	// prepare the SQP object
	nConstr=triRSVS.meshDep->CountVoluParent(); 
	nDv=triRSVS.snakeDep->snaxs.size();
	BuildMathArrays(nDv, nConstr);

	// TODO this needs to be supported by mapping each volume to the constraint position
	// There can be more than one constraint for each cell.
	BuildConstrMap(triRSVS);

	vecin.clear();
	vecin.reserve(nDv);
	for(ii=0; ii<nDv; ++ii){
		vecin.push_back(triRSVS.snakeDep->snaxs(ii)->index);
	}
	BuildDVMap(vecin);

	// Calculate the SQP object
	ni=triRSVS.dynatri.size();
	for(ii = 0; ii< ni ; ii++){
		CalcTriangle(*(triRSVS.dynatri(ii)), triRSVS);
	} 
	ni=triRSVS.intertri.size();
	for(ii = 0; ii< ni ; ii++){
		CalcTriangle(*(triRSVS.intertri(ii)), triRSVS);
	} 
	ni=triRSVS.acttri.size();
	for(ii = 0; ii< ni ; ii++){
		CalcTriangle(*(triRSVS.stattri.isearch(triRSVS.acttri[ii])), triRSVS);
	} 

	// Output some data to check it makes sense
	ReturnConstrToMesh(triRSVS);
}

void SQPcalc::Print2Screen()const {
	cout << "Math result: obj " << obj << " false access " << falseaccess << endl;
	cout << "Constr " ;
	for (int i = 0; i < nConstr; ++i)
	{
		cout << constr[i] << " ";
	}
	cout << endl;
}

void SQPcalc::ReturnConstrToMesh(const triangulation &triRSVS) const {
	
	int ii, ni;
	vector<double> temp;
	ni=constr.size();
	temp.reserve(ni);

	for(ii=0; ii<ni; ii++){
		temp.push_back(constr[ii]);
	}

	triRSVS.meshDep->MapFill2Parent(temp, this->constrList);


}

void SQPcalc::BuildMathArrays(int nDvIn, int nConstrIn){
	// Builds the target math arrays
	
	nDv=nDvIn;
	nConstr=nConstrIn;

	dConstr.setZero(nConstr,nDv);
	HConstr.setZero(nDv,nDv); 
	HObj.setZero(nDv,nDv);  
	dObj.setZero(nDv);
	constr.setZero(nConstr);
	lagMult.setZero(nConstr);
}


void SQPcalc::BuildConstrMap(const triangulation &triangleRSVS){

	// explore parents of mesh adding 1 by 1 elemind
	// for each parent
	// for each snakemesh.volu
	// Assign to constrMap.targ = the position in parent.volu of the parentconn
	triangleRSVS.meshDep->ReturnParentMap(constrMap.vec,constrMap.targ,constrList);

	constrMap.GenerateHash();
}

void SQPcalc::BuildDVMap(const vector<int> &vecin){
	dvMap.vec=vecin;
	sort(dvMap.vec);
	unique(dvMap.vec);
	dvMap.GenerateHash();
}

void SQPcalc::CalcTriangle(const triangle& triIn, const triangulation &triRSVS){


	int ii,ni,jj,nj,kk,nCellTarg;
	int isCentre,subTemp,subTemp1,subTemp2,subTemp3,nDvAct;
	SurfCentroid centreCalc;
	Volume VolumeCalc;
	Area AreaCalc;
	vector<int> dvList,subTempVec;
	HashedVector<int,int> dvListMap;
	vector<vector<double> const *> veccoord;
	MatrixXd HPos, dPos;
	MatrixXd HVal, dVal;
	MatrixXd dConstrPart,HConstrPart, HObjPart;
	RowVectorXd dObjPart;
	double constrPart, objPart;
	ArrayVec<double>* HValpnt=NULL, *dValpnt=NULL;
	double * retVal;


	veccoord.reserve(3);
	isCentre=0;
	ni=3;
	for(ii=0; ii<ni; ++ii){
		if(triIn.pointtype[ii]==1){
			veccoord.push_back(&(triRSVS.meshDep->verts.isearch(triIn.pointind[ii])->coord));
		} else if (triIn.pointtype[ii]==2){
			veccoord.push_back(&(triRSVS.snakeDep->snakeconn.verts.isearch(triIn.pointind[ii])->coord));
		} else if (triIn.pointtype[ii]==3){
			veccoord.push_back((triRSVS.trivert.isearch(triIn.pointind[ii])->coord.retPtr()));
			isCentre++;
		}
	}

	// Constr and objective
	VolumeCalc.assign(veccoord[0],veccoord[1],veccoord[2]);
	AreaCalc.assign(veccoord[0],veccoord[1],veccoord[2]);

	// Active DV lists

	for(ii=0; ii<ni; ++ii){
		if (triIn.pointtype[ii]==2){
			dvList.push_back(triIn.pointind[ii]);
		} else if (triIn.pointtype[ii]==3 && false){
			subTemp=triRSVS.trivert.find(triIn.pointind[ii]);
			nj=triRSVS.trisurf(subTemp)->indvert.size();
			for(jj=0;jj<nj;++jj){
				if (triRSVS.trisurf(subTemp)->typevert[jj]==2){
					dvList.push_back(triRSVS.trisurf(subTemp)->indvert[jj]);
				}
			}
		}
	}
	dvListMap.vec=dvList;
	sort(dvListMap.vec);
	unique(dvListMap.vec);
	dvListMap.GenerateHash();
	nDvAct=dvListMap.vec.size();

	// Positional Derivatives

	// HERE -> function to calculate SurfCentroid (dc/dd)^T Hm (dc/dd)

	HPos.setZero(nDvAct,nDvAct);
	dPos.setZero(9,nDvAct);
	kk=0;
	for(ii=0; ii<ni; ++ii){
		if (triIn.pointtype[ii]==2){
			subTemp=triRSVS.snakeDep->snaxs.find(triIn.pointind[ii]);
			subTemp1=triRSVS.meshDep->verts.find(triRSVS.snakeDep->snaxs(subTemp)->fromvert);
			subTemp2=triRSVS.meshDep->verts.find(triRSVS.snakeDep->snaxs(subTemp)->tovert);
			subTemp3=dvListMap.find(triIn.pointind[ii]);
			for(jj=0;jj<3;++jj){
				dPos(kk*3+jj,subTemp3)+=triRSVS.meshDep->verts(subTemp2)->coord[jj]
					-triRSVS.meshDep->verts(subTemp1)->coord[jj];
			}
			kk++;
		} else if (triIn.pointtype[ii]==3 && false){

		}
	}
	// Total
	VolumeCalc.Calc();
	AreaCalc.Calc();

	HVal.setZero(9,9);
	dVal.setZero(1,9);
	dConstrPart.setZero(1,nDvAct);
	dObjPart.setZero(1,nDvAct);
	HConstrPart.setZero(nDvAct,nDvAct);
	HObjPart.setZero(nDvAct,nDvAct);

	VolumeCalc.ReturnDatPoint(&retVal, &dValpnt, &HValpnt); 
	ArrayVec2MatrixXd(*HValpnt, HVal);
	ArrayVec2MatrixXd(*dValpnt, dVal);
	Deriv1stChainScalar(dVal, dPos,dConstrPart);
	Deriv2ndChainScalar(dVal,dPos,HVal,HPos,HConstrPart);
	constrPart=*retVal;

	AreaCalc.ReturnDatPoint(&retVal, &dValpnt, &HValpnt); 
	ArrayVec2MatrixXd(*HValpnt, HVal);
	ArrayVec2MatrixXd(*dValpnt, dVal);
	Deriv1stChainScalar(dVal, dPos,dConstrPart);
	Deriv2ndChainScalar(dVal,dPos,HVal,HPos,HConstrPart);
	objPart=*retVal;

	// Assign to main part of the object

	obj+=objPart; // assign objective function

	// assign constraint value
	nCellTarg=triIn.connec.celltarg.size(); 
	for(ii=0; ii< nCellTarg;++ii){
		subTempVec=constrMap.findall(triIn.connec.celltarg[ii]);
		nj=subTempVec.size();
		for(jj=0; jj< nj; ++jj){
			if (subTempVec[jj]!=-1){
				constr[subTempVec[jj]] += triIn.connec.constrinfluence[ii]*constrPart;
			} else {
				falseaccess++;
			}
		}
	}
} 

void ArrayVec2MatrixXd(const ArrayVec<double> &arrayIn, MatrixXd &matOut){

	int nR,nC;
	arrayIn.size(nR,nC);
	matOut.setZero(nR,nC);

	for(int i=0; i< nR; ++i){
		for(int j=0; j< nC; ++j){
			matOut(i,j)=arrayIn[i][j];
		}
	}

}

void VecBy3DimArray(const MatrixXd &vec, const MatrixXd &arr3dim, MatrixXd &retArray){
	/*
	3 Dimensional array is expected as 
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
	vector expected as a row vector
	[xc yc zc]

	will sum the three Hessians according to the values in the vector
	*/

	int nRow,nCol,nVec,nColFin;
	int ii,jj,kk;
	nRow=arr3dim.rows();
	nCol=arr3dim.cols();
	nVec=vec.rows();
	nColFin=nCol/nVec;

	#ifdef SAFE_ALGO
	// size checks
	if ((nVec*nColFin)!=nCol){
		throw invalid_argument("Sizes do not match in 3D array collapse");
	}
	#endif

	retArray.setZero(nRow,nColFin);
	// needs to add checks for matching sizes
	for (ii=0;ii<nRow;ii++){
		for (jj=0;jj<nColFin;jj++){
			for (kk=0;kk<nVec;kk++){
				retArray(ii,jj)+=arr3dim(ii,jj*nVec+kk)*vec(0,kk);
			}	
		}	
	}



}



void Deriv1stChainScalar(const MatrixXd &dSdc,const MatrixXd &dcdd, MatrixXd &dSdd){
	dSdd=dSdc*dcdd;
}

void Deriv2ndChainScalar(const MatrixXd &dSdc,const MatrixXd &dcdd,
	const MatrixXd &HSc,const MatrixXd &Hcd,MatrixXd &HSd){
	
	VecBy3DimArray(dSdc, Hcd, HSd);
	HSd+=(dcdd.transpose()*HSc*dcdd);

}

/// INSPIRATION
/*
void SnakeSurfaceCentroid_fun(coordvec &coord,const surf &surfin, const mesh& meshin){ 
	int ii,n;
	vector<int> vertind;
	vector<vector<double> const *> veccoord;
	SurfCentroid tempCalc;
	ArrayVec<double> tempCoord,jac,hes;

	coord.assign(0,0,0);
	n=int(surfin.edgeind.size());

	veccoord.reserve(n);
	ConnVertFromConnEdge(meshin, surfin.edgeind,vertind);

	for(ii=0; ii<n; ++ii){
		veccoord.push_back(&(meshin.verts.isearch(vertind[ii])->coord));
	}

	tempCalc.assign(veccoord);
	tempCalc.Calc();

	tempCalc.ReturnDat(tempCoord,jac,hes);
	coord.assign(tempCoord[0][0],tempCoord[1][0],tempCoord[2][0]);
}

void HybridSurfaceCentroid_fun(coordvec &coord,const trianglesurf &surfin, const mesh& meshin,
	const mesh& snakeconn){ 
	int ii,n;
	vector<int> vertind;
	vector<vector<double> const *> veccoord;
	SurfCentroid tempCalc;
	ArrayVec<double> tempCoord,jac,hes;

	coord.assign(0,0,0);

	n=surfin.indvert.size();
	for(ii=0; ii<n; ++ii){
		if(surfin.typevert[ii]==1){
			veccoord.push_back(&(meshin.verts.isearch(surfin.indvert[ii])->coord));
		} else if (surfin.typevert[ii]==2){
			veccoord.push_back(&(snakeconn.verts.isearch(surfin.indvert[ii])->coord));
		}
	}

	tempCalc.assign(veccoord);
	tempCalc.Calc();

	tempCalc.ReturnDat(tempCoord,jac,hes);
	coord.assign(tempCoord[0][0],tempCoord[1][0],tempCoord[2][0]);
}
*/

