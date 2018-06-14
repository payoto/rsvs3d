
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
#include "snakstruct.hpp" 


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
	vector<int> acttri;
	triarray stattri;
	triarray dynatri;
	tripointarray trivert;
};




class triangle : public meshpart , public snakpart {	
private:
	bool isTriangleReady=false;
public:
	
	vector<int> pointtype; // 1=mesh vertex 2=snaxel 3=trianglepoint
	vector<int> pointind;
	int constrind; // Element in the constraint vector
	double constrinfluence; //usually 1 or -1

	// interface functions
	void disp() const;
	int Key() const {return (index);};
	int KeyParent() const {return (constrind);};
	void ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu);
	void PrepareForUse(){};
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wunused-parameter"
	bool isready(bool isInMesh) const {
		return(isTriangleReady);
	}
	void SwitchIndex(int typeInd, int oldInd, int newInd){}
	void disptree(const mesh &meshin, int n) const {};
	#pragma GCC diagnostic pop
	void read(FILE * fid);
	void write(FILE * fid) const;
	void TightenConnectivity(){}
	
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
	int parentsurf;

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
	void TightenConnectivity(){}

};

void CalculateSnakeVel(snake &snakein);


#endif // SNAKEVEL_H_INCLUDED

