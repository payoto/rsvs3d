
#include <iostream>
#include <cstdlib>
#include "arraystructures.hpp"
#include "snake.hpp" 
#include "snakevel.hpp"

using namespace std;	

void CalculateSnakeVel(snake &snakein){

	int ii=0;

	for (ii=0;ii<int(snakein.snaxs.size());++ii){
		if (snakein.snaxs(ii)->isfreeze==1){
			snakein.snaxs[ii].v=(0.5-snakein.snaxs[ii].d)*0.3;
			snakein.snaxs[ii].isfreeze=0;
		}
		//snakein.snaxs[ii].v=(double(rand()%1001)/1000.0+0.5)*snakein.snaxs[ii].v;
		//snakein.snaxs[ii].v=(0.4*(double(rand()%1001)/1000.0)+0.8)*snakein.snaxs[ii].v;
	}
}



void TriangulateMesh(const mesh& meshin, triangulation &triangleRSVS){

	TriangulateContainer(meshin,triangleRSVS , 1);

}
void TriangulateSnake(const snake& snakein, triangulation &triangleRSVS){

		TriangulateContainer(snakein.snakeconn, triangleRSVS , 2);
}


void TriangulateContainer(const mesh& meshin, triangulation &triangleRSVS , const int typeMesh){
	int ii,n,maxInd;
	triarray triangulation::*mp;

	if (typeMesh==1){
		mp=&triangulation::stattri;
	} else {
		mp=&triangulation::dynatri;
	}

	maxInd=triangleRSVS.trivert.GetMaxIndex();
	n=meshin.surfs.size();

	for (ii=0; ii<n; ++ii){
		TriangulateSurface((*meshin.surfs(ii)),meshin,triangleRSVS.*mp, 
			triangleRSVS.trivert, typeMesh, maxInd+ii+1);
	}
}

void TriangulateSurface(const surf &surfin,const mesh& meshin, 
	triarray &triangul, tripointarray& trivert, const int typeMesh, int trivertMaxInd){
	// Generates the triangulation parts 
	// typeMesh=1 is a static mesh, type 2 is a dynamic one (snake)

	int ii,n;
	coordvec edgeCentre;
	double surfLength,edgeLength;
	trianglepoint surfCentre;
	triangle triangleEdge;

	// Loop around edges calculating centre and length of each edge
	surfLength=0;
	edgeLength=0;
	n=int(surfin.edgeind.size());

	triangleEdge.SetPointType(typeMesh,typeMesh,3);
	triangleEdge.pointind[2]=trivertMaxInd+1;
	triangleEdge.parentsurf=surfin.index;

	for(ii=0; ii<n; ++ii){
		meshin.edges.isearch(surfin.edgeind[ii])->GeometricProperties(&meshin,edgeCentre,edgeLength);
		edgeCentre.mult(edgeLength);
		surfCentre.coord.add(edgeCentre.usedata());
		surfLength+=edgeLength;

		triangleEdge.pointind[0]=meshin.edges.isearch(surfin.edgeind[ii])->vertind[0];
		triangleEdge.pointind[1]=meshin.edges.isearch(surfin.edgeind[ii])->vertind[1];
		triangul.push_back(triangleEdge);
	}

	surfCentre.coord.div(surfLength);
	surfCentre.index=trivertMaxInd+1;
	surfCentre.parentsurf=surfin.index;
	surfCentre.parentType=typeMesh;
}