/**
 * Provides a triangulation for snake and meshes.
 * 
 * This links an active snake and mesh to their triangulated representation 
 * needed to compute areas and volumes. 
 * 
 *@file
 */

//===============================================
// Include Guards
#ifndef TRIANGULATE_H_INCLUDED
#define TRIANGULATE_H_INCLUDED

//===============================================
// Levels of debuging Guards
#ifdef DEBUGLVL2 // All Debugging calls
#define DEBUGLVL1
#define TEST_ALL
#endif

#ifdef DEBUGLVL1

#endif

//========================================
// Forward declared dependencies
class RSVScalc;
class SurfCentroid;

//=================================
// included dependencies
#include <vector>
#include "arraystructures.hpp"
#include "mesh.hpp"
#include "snake.hpp"

//==================================
// Code
// NOTE: function in a class definition are IMPLICITELY INLINED 
//       ie replaced by their code at compile time
using namespace std;

template <class T> class TriStruct;

class triangle;
class trianglepoint;
class trianglesurf;
class tecplotfile; 

// typedef TriStruct<triangle> triarray;
// typedef TriStruct<trianglepoint>  tripointarray;
// typedef TriStruct<trianglesurf>  trisurfarray;

// Template Class
namespace rsvs3d {
	namespace constants{
		namespace meshtypes {
			static const int mesh = 1;
			static const int snake = 2;
			static const int triangulation = 3;
		}
	}
}
template <class T>
class TriStruct : public SnakStruct<T>{
	friend class triangulation;
};

class triarray : public TriStruct<triangle>{};
class tripointarray : public TriStruct<trianglepoint>{};
class trisurfarray : public TriStruct<trianglesurf>{};

// Base classes
class triangulation  
{
public:
	vector<int> acttri; // 
	triarray stattri;
	triarray dynatri;
	triarray intertri;
	
	tripointarray trivert;
	trisurfarray trisurf;

	snake* snakeDep=NULL; // Pointer to the Snake referred to in triangles
	mesh* meshDep=NULL; // Pointer to the Mesh referred to in triangles

	void disp() const;
	void PrepareForUse();
	void CleanDynaTri();
	void CalcTriVertPosDyna(int ii);
	void CalcTriVertPosDyna();
	void CalcTriVertPos(int ii);
	void CalcTriVertPos();
	void SetActiveStaticTri();

	void SetConnectivity();
	void SetConnectivityStat(int ii);
	void SetConnectivityInter(int ii);
	void SetConnectivityDyna(int ii);

	triangulation(){
		meshDep=NULL;
		snakeDep=NULL;
	}
	explicit triangulation(mesh &meshin){
		meshDep=&meshin;
		snakeDep=NULL;
	}
};


class tri2mesh {
public:
	vector<int> celltarg; // cell need to be indexed in a support array which maps to constraints etc
	vector<double> constrinfluence; // +1 -1 to indicate 
};


class triangle : public meshpart , public snakpart {	
private:
	bool isTriangleReady=false;
public:

	vector<int> pointtype; // 1=mesh vertex 2=snaxel 3=trianglepoint
	vector<int> pointind;
	int parentsurf=0; // Surface in the snakemesh() needs to be converted to constr
	int parenttype=0;

	//double constrinfluence=0; //usually 1 or -1 to do with the ordering
	tri2mesh connec;
 
	// interface functions
	void disp() const override;
	void disptree(const mesh &meshin, int n) const override;
	int Key() const override {return (index);};
	int KeyParent() const override {return (parentsurf);};
	void ChangeIndices (int nVert,int nEdge,int nSurf,int nVolu) override;
	void PrepareForUse() override {};
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wunused-parameter"
	bool isready(bool isInMesh) const override{
		return(isTriangleReady);
	}
	void SwitchIndex(int typeInd, int oldInd, int newInd){}
	#pragma GCC diagnostic pop
	void read(FILE * fid) override;
	void write(FILE * fid) const override;
	void TightenConnectivity() override{}
	void SetPointType(int a, int b, int c){	pointtype[0]=a;pointtype[1]=b;pointtype[2]=c;}
	triangle(){ // Constructor
		index=0;

		pointtype.assign(3,0);
		pointind.assign(3,0);
		isTriangleReady=false;
	}
};


class trianglepoint : public meshpart , public snakpart {	
public:
	
	coordvec coord;
	int parentsurf=0;
	int parentType=0;
	int nInfluences=0;

	// interface functions
	void disp() const override;
	void disptree(const mesh &meshin, int n) const override;
	int Key() const override {return (index);};
	int KeyParent() const override {return (parentsurf);};
	void ChangeIndices (int nVert,int nEdge,int nSurf,int nVolu) override;
	void ChangeIndicesSnakeMesh(int nVert,int nEdge,int nSurf,int nVolu);
	void PrepareForUse() override {};
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wunused-parameter"
	bool isready(bool isInMesh) const override {return(true);}
	void SwitchIndex(int typeInd, int oldInd, int newInd){}
	#pragma GCC diagnostic pop
	void read(FILE * fid) override;
	void write(FILE * fid) const override;
	void TightenConnectivity() override {}

};

class trianglesurf : public meshpart , public snakpart {

public:
	vector<int> indvert;
	vector<int> typevert;
	vector<int> voluind;
	int parentsurfmesh=0;
	
	// interface functions
	void disp() const override;
	void disptree(const mesh &meshin, int n) const override;
	int Key() const override {return (index);};
	int KeyParent() const override {return (parentsurfmesh);};
	void ChangeIndices (int nVert,int nEdge,int nSurf,int nVolu) override;
	void ChangeIndicesSnakeMesh(int nVert,int nEdge,int nSurf,int nVolu);
	void PrepareForUse() override {};
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wunused-parameter"
	bool isready(bool isInMesh) const override {return(true);}
	void SwitchIndex(int typeInd, int oldInd, int newInd){}
	#pragma GCC diagnostic pop
	void read(FILE * fid) override;
	void write(FILE * fid) const override;
	void TightenConnectivity() override {}

};

void CalculateSnakeVel(snake &snakein);
void CalculateSnakeVelRand(snake &snakein);
void CalculateSnakeVelUnit(snake &snakein);
void CalculateSnakeVelFast(snake &snakein);
void CalculateNoNanSnakeVel(snake &snakein, double deltaStep=0.01);
void TriangulateSurface(const surf &surfin,const mesh& meshin, 
	triarray &triangul, tripointarray& trivert, const int typeMesh, int trivertMaxInd);
void TriangulateTriSurface(const trianglesurf &surfin,
	triarray &triangul, tripointarray& trivert, const int typeMesh, int trivertMaxInd);
void TriangulateContainer(const mesh& meshin, triangulation &triangleRSVS , const int typeMesh, const vector<int> &subList={});
void TriangulateSnake(snake& snakein, triangulation &triangleRSVS);
void TriangulateMesh(mesh& meshin, triangulation &triangleRSVS);
void MeshTriangulation(mesh &meshout,const mesh& meshin,triarray &triangul, tripointarray& trivert);
void MaintainTriangulateSnake(triangulation &triangleRSVS);
void SnakeSurfaceCentroid_fun(coordvec &coord,const surf &surfin, const mesh& meshin);
void HybridSurfaceCentroid_fun(coordvec &coord,const trianglesurf &surfin, const mesh& meshin, const mesh& snakeconn);

void Test_stepalgoRSVS(snake &testSnake,triangulation &RSVStri , vector<double> &dt,
	vector<int> &isImpact, RSVScalc &calcObj, tecplotfile &outSnake2, double totT);
void BuildTriSurfGridSnakeIntersect(triangulation &triangleRSVS);
int FollowVertexConnection(int actVert, int prevEdge, const HashedVector<int,int> &edgeSurfInd,	const HashedVector<int,int> &vertSurfInd, const snake &snakeRSVS, const mesh &meshRSVS, int &returnIndex,int &returnType, int &nextEdge);
int FollowSnaxelDirection(int actSnax,const snake &snakeRSVS, int &returnIndex, int &returnType, int &actEdge);
bool FollowSnaxEdgeConnection(int actSnax, int actSurf,int followSnaxEdge,  const snake &snakeRSVS, vector<bool> &isSnaxEdgeDone, int & returnIndex);
mesh TriarrayToMesh(const triangulation& triangul, const triarray& triin);
void FlattenBoundaryFaces(mesh &meshin);
SurfCentroid SurfaceCentroid_TriangleSurf(const trianglesurf &surfin, 
	const mesh& meshin,	const mesh& snakeconn);
SurfCentroid SurfaceCentroid_SnakeSurf(const surf &surfin,
	const mesh& meshin);

int Test_snakeRSVS();
int Test_surfcentre();
int Test_snakeRSVS_singlevol();
int Test_RSVSalgo_singlevol();
int Test_MeshOrient();
#endif // TRIANGULATE_H_INCLUDED

