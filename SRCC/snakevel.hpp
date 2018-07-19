
//===============================================
// Include Guards
#ifndef SNAKEVEL_H_INCLUDED
#define SNAKEVEL_H_INCLUDED

//===============================================
// Levels of debuging Guards
#ifdef DEBUGLVL2 // All Debugging calls
#define DEBUGLVL1
#define TEST_ALL
#endif

#ifdef DEBUGLVL1 // Debugging of new features.

#endif

//=================================
// forward declared dependencies
// 		class foo; //when you only need a pointer not the actual object
// 		and to avoid circular dependencies


//=================================
// included dependencies
#include <vector>
#include "arraystructures.hpp"
#include "snake.hpp" 
#include "mesh.hpp" 


//==================================
// Code
// NOTE: function in a class definition are IMPLICITELY INLINED 
//       ie replaced by their code at compile time
using namespace std;

class triangle;
class trianglepoint;
typedef SnakStruct<triangle> triarray;
typedef SnakStruct<trianglepoint>  tripointarray;

// Template Class


// Base classes
class triangulation  
{
public:
	vector<int> acttri; // 
	triarray stattri;
	triarray dynatri;
	tripointarray trivert;

	snake* snakeDep=NULL; // Pointer to the Snake referred to in triangles
	mesh* meshDep=NULL; // Pointer to the Mesh referred to in triangles

	void disp() const;
	void PrepareForUse();
	void CleanDynaTri();
};




class triangle : public meshpart , public snakpart {	
private:
	bool isTriangleReady=false;
public:
	
	vector<int> pointtype; // 1=mesh vertex 2=snaxel 3=trianglepoint
	vector<int> pointind;
	int parentsurf=0; // Element in the constraint vector
	double constrinfluence=0; //usually 1 or -1 to do with the ordering

	// interface functions
	void disp() const;
	void disptree(const mesh &meshin, int n) const;
	int Key() const {return (index);};
	int KeyParent() const {return (parentsurf);};
	void ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu);
	void PrepareForUse(){};
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wunused-parameter"
	bool isready(bool isInMesh) const {
		return(isTriangleReady);
	}
	void SwitchIndex(int typeInd, int oldInd, int newInd){}
	#pragma GCC diagnostic pop
	void read(FILE * fid);
	void write(FILE * fid) const;
	void TightenConnectivity(){}
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
	// interface functions
	void disp() const;
	void disptree(const mesh &meshin, int n) const;
	int Key() const {return (index);};
	int KeyParent() const {return (parentsurf);};
	void ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu);
	void ChangeIndicesSnakeMesh(int nVert,int nEdge,int nSurf,int nVolu);
	void PrepareForUse(){};
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wunused-parameter"
	bool isready(bool isInMesh) const {return(true);}
	void SwitchIndex(int typeInd, int oldInd, int newInd){}
	#pragma GCC diagnostic pop
	void read(FILE * fid);
	void write(FILE * fid) const;
	void TightenConnectivity() override {}

};

void CalculateSnakeVel(snake &snakein);
void TriangulateSurface(const surf &surfin,const mesh& meshin, 
	triarray &triangul, tripointarray& trivert, const int typeMesh, int trivertMaxInd);
void TriangulateContainer(const mesh& meshin, triangulation &triangleRSVS , const int typeMesh, const vector<int> &subList={});
void TriangulateSnake(snake& snakein, triangulation &triangleRSVS);
void TriangulateMesh(mesh& meshin, triangulation &triangleRSVS);
void MeshTriangulation(mesh &meshout,const mesh& meshin,triarray &triangul, tripointarray& trivert);
void MaintainTriangulateSnake(triangulation &triangleRSVS);
void SurfaceCentroid_fun(coordvec &coord,const surf &surfin, const mesh& meshin);

void Test_stepalgoRSVS(snake &testSnake,triangulation &RSVStri , vector<double> dt, vector<int> isImpact, tecplotfile &outSnake);
int Test_snakeRSVS();
int Test_surfcentre();
#endif // SNAKEVEL_H_INCLUDED

