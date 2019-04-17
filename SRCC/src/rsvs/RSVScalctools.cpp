#include <vector>
#include <Eigen>

#include "arraystructures.hpp"
#include "triangulate.hpp"
#include "RSVScalctools.hpp"
#include "RSVScalc.hpp"
#include "RSVSmath.hpp"


using namespace Eigen;
using namespace rsvs3d::constants;

vector<vector<double> const *> TrianglePointerCoordinates(const triangle &triIn,
	const triangulation& triRSVS){

	vector<vector<double> const *> veccoord;
	veccoord.reserve(3);

	int ni=3;
	for(int ii=0; ii<ni; ++ii){
		if(triIn.pointtype[ii]==meshtypes::mesh){
			veccoord.push_back(&(triRSVS.meshDep->verts.isearch(triIn.pointind[ii])->coord));
		} else if (triIn.pointtype[ii]==meshtypes::snake){
			veccoord.push_back(&(triRSVS.snakeDep->snakeconn.verts.isearch(triIn.pointind[ii])->coord));
		} else if (triIn.pointtype[ii]==meshtypes::triangulation){
			veccoord.push_back((triRSVS.trivert.isearch(triIn.pointind[ii])->coord.retPtr()));
		}
	}

	return veccoord;
}

HashedVector<int,int> TriangleActiveDesignVariables(const triangle &triIn,
	const triangulation& triRSVS, const HashedVector<int, int> &objDvMap){
	
	HashedVector<int,int> dvListMap;
	vector<int> vertsSurf;
	vector<int> &dvList = dvListMap.vec;

	int ni=3;
	for(int ii=0; ii<ni; ++ii){
		vertsSurf.clear();
		if (triIn.pointtype[ii]==meshtypes::snake 
			&& objDvMap.find(triIn.pointind[ii])!=__notfound){
			dvList.push_back(triIn.pointind[ii]);
		} else if (triIn.pointtype[ii]==meshtypes::triangulation && true){
			switch (triRSVS.trivert.isearch(triIn.pointind[ii])->parentType){
				case meshtypes::triangulation:
				{
					// // If parent surf is a trisurf
					// auto surfCurr = triRSVS.trisurf.isearch(
					// 	triRSVS.trivert.isearch(triIn.pointind[ii])->parentsurf);
					// int nj=surfCurr->indvert.size();
					// for(int jj=0;jj<nj;++jj){
					// 	if (surfCurr->typevert[jj]==meshtypes::snake){
					// 		if (objDvMap.find(surfCurr->indvert[jj])
					// 			!=__notfound){
					// 			dvList.push_back(surfCurr->indvert[jj]);
					// 		}
					// 	}
					// }
				}
				break;
				case meshtypes::snake:
				{
					// If parent surf is a snakesurf
					auto surfCurr =  triRSVS.snakeDep->snakeconn.surfs.isearch(
						triRSVS.trivert.isearch(triIn.pointind[ii])->parentsurf);
					vertsSurf = surfCurr->OrderedVerts(&(triRSVS.snakeDep->snakeconn));
					int nj = vertsSurf.size();
					for(int jj=0;jj<nj;++jj){
						if (objDvMap.find(vertsSurf[jj])!=__notfound){
							dvList.push_back(vertsSurf[jj]);
						}
					}
				}
				break;
			} 
		}
	}

	sort(dvListMap.vec);
	unique(dvListMap.vec);
	dvListMap.GenerateHash();

	return(dvListMap);
}

/**
 * This relies on the RSVSmath objects to provide a the positional derivatives
 * as well as the known derivatives of the snaxel 
 * ($\frac{\partial coord_i}{\partial d} = \Delta g_i$ and the second 
 * derivatives) are 0.
 */
void TrianglePositionDerivatives(const triangle &triIn, 
	const triangulation &triRSVS, const HashedVector<int,int> &dvListMap,
	MatrixXd &dPos, MatrixXd &HPos){

	std::vector<int> vertsSurf;
	int nDvAct=dvListMap.vec.size();

	// Positional Derivatives

	// HERE -> function to calculate SurfCentroid (dc/dd)^T Hm (dc/dd)
	HPos.setZero(nDvAct,nDvAct*9);
	dPos.setZero(9,nDvAct);
	int ni=3;
	for(int ii=0; ii<ni; ++ii){
		vertsSurf.clear();
		if (triIn.pointtype[ii]==meshtypes::snake
			&& dvListMap.find(triIn.pointind[ii])!=__notfound){
			auto tempSnax = triRSVS.snakeDep->snaxs.isearch(triIn.pointind[ii]);
			auto fromVert = triRSVS.meshDep->verts.isearch(tempSnax->fromvert);
			auto toVert = triRSVS.meshDep->verts.isearch(tempSnax->tovert);
			int dvSub = dvListMap.find(triIn.pointind[ii]);
			for(int jj=0; jj<3; ++jj){
				dPos(ii*3+jj, dvSub) 
					+= (toVert->coord[jj] - fromVert->coord[jj]);
			}
		} else if (triIn.pointtype[ii]==meshtypes::triangulation && true){
			switch (triRSVS.trivert.isearch(triIn.pointind[ii])->parentType){
				case meshtypes::triangulation:
				{
					// auto surfCurr=triRSVS.trisurf.isearch(
					// 	triRSVS.trivert.isearch(triIn.pointind[ii])->parentsurf);
					// auto surfCentre =  SurfaceCentroid_TriangleSurf(*surfCurr, 
					// 	*(triRSVS.meshDep), triRSVS.snakeDep->snakeconn);
					// surfCentre.Calc();
					// auto jacCentre = surfCentre.jac_ptr();
					// auto hesCentre = surfCentre.hes_ptr();
					// int count = surfCurr->indvert.size();
					// for (int j = 0; j < count; ++j)
					// {
					// 	int dvSub = dvListMap.find(surfCurr->indvert[j]);
					// 	if (dvSub!=__notfound 
					// 		&& surfCurr->typevert[j]==meshtypes::snake)
					// 	{
					// 		auto currSnax = triRSVS.snakeDep->snaxs.isearch(
					// 			surfCurr->indvert[j]);
					// 		auto fromVert = triRSVS.meshDep->verts.isearch(
					// 			currSnax->fromvert);
					// 		auto toVert = triRSVS.meshDep->verts.isearch(
					// 			currSnax->tovert);
					// 		for (int k = 0; k < 3; ++k)
					// 		{
					// 			for (int l = 0; l < 3; ++l)
					// 			{
					// 				dPos(ii*3+k,dvSub) += jacCentre[k][j+count*l]
					// 					*(toVert->coord[l]-fromVert->coord[l]);
					// 			}
					// 		}
					// 	}
					// }
				}
				break;
				case meshtypes::snake:
				{
					auto surfCurr =  triRSVS.snakeDep->snakeconn.surfs.isearch(
						triRSVS.trivert.isearch(triIn.pointind[ii])->parentsurf);
					auto surfCentre =  SurfaceCentroid_SnakeSurf(
						*surfCurr, triRSVS.snakeDep->snakeconn);
					surfCentre.Calc();
					auto jacCentre = surfCentre.jac_ptr();
					auto hesCentre = surfCentre.hes_ptr();
					ConnVertFromConnEdge(triRSVS.snakeDep->snakeconn, 
						surfCurr->edgeind, vertsSurf);
					// Check that the triangle and the surface 
					int count = vertsSurf.size();
					for (int j = 0; j < count; ++j)
					{
						int dvSub = dvListMap.find(vertsSurf[j]);
						auto currSnax = triRSVS.snakeDep->snaxs.isearch(vertsSurf[j]);
						auto fromVert = triRSVS.meshDep->verts.isearch(
							currSnax->fromvert);
						auto toVert = triRSVS.meshDep->verts.isearch(
							currSnax->tovert);

						if (dvSub!=__notfound){
							for (int k = 0; k < 3; ++k)
							{
								for (int l = 0; l < 3; ++l)
								{
									dPos(ii*3+k,dvSub) += 
										(isnan(jacCentre[k][j+count*l])?
											0.0:jacCentre[k][j+count*l])
										*(toVert->coord[l]-fromVert->coord[l]);
								}
							}
							// HPos(dvSub,dvSub) += 
						}
					}
				}
				break;
			}
		}
	}
}
