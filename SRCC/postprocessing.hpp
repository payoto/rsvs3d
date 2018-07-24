
//===============================================
// Include Guards
#ifndef POSTPROCESSING_H_INCLUDED
#define POSTPROCESSING_H_INCLUDED

//===============================================
// Levels of debuging Guards
#ifdef DEBUGLVL2 // All Debugging calls
#define DEBUGLVL1

#endif

#ifdef DEBUGLVL1 // Debugging of new features.
#define TEST_POSTPROCESSING
#endif

//=================================
// forward declared dependencies
// 		class foo; //when you only need a pointer not the actual object
// 		and to avoid circular dependencies

class mesh;

//=================================
// included dependencies
#include <iostream>
#include <stdarg.h>
#include "arraystructures.hpp"
#include "snakevel.hpp" 


//==================================
// Code
// NOTE: function in a class definition are IMPLICITELY INLINED 
//       ie replaced by their code at compile time
using namespace std;


// Base classes

class tecplotfile {
private:
	FILE *fid;
	int lengthLine;
	int nZones=0;
public:
	int OpenFile(const char *str);
	void CloseFile();

	// Mesh out
	int PrintMesh(const mesh& meshout,int strandID=0, double timeStep=0, int forceOutType=0);
	int VolDataBlock(const mesh& meshout,int nVert,int nVolu, int nVertDat);
	int SurfDataBlock(const mesh& meshout,int nVert,int nSurf, int nVertDat);
	int LineDataBlock(const mesh &meshout,int nVert,int nEdge, int nVertDat,int nCellDat);
	int VolFaceMap(const mesh& meshout,int nSurf);
	int SurfFaceMap(const mesh& meshout,int nEdge);
	int LineFaceMap(const mesh& meshout,int nEdge);

	// Triangulation out
	int VolDataBlock(const triangulation &triout, triarray triangulation::*mp,int nVert,int nVolu, int nVertDat);
	int SurfDataBlock(const triangulation &triout, triarray triangulation::*mp,int nVert,int nSurf, int nVertDat);
	int LineDataBlock(const triangulation &triout, triarray triangulation::*mp,int nVert,int nEdge, int nVertDat,int nCellDat);
	int SurfFaceMap(const triangulation &triout, triarray triangulation::*mp,int nEdge);
	int LineFaceMap(const triangulation &triout, triarray triangulation::*mp,int nEdge);
	int VolFaceMap(const triangulation &triout, triarray triangulation::*mp,int nSurf);
	int PrintMesh(const triangulation &triout, triarray triangulation::*mp,int strandID, double timeStep, int forceOutType);

	void ZoneHeaderPolyhedron(int nVert, int nVolu, int nSurf, int totNumFaceNode,int nVertDat, int nCellDat);

	void ZoneHeaderPolygon(int nVert,int nEdge, int nSurf,int nVertDat, int nCellDat);
	void ZoneHeaderFelineseg(int nVert,int nEdge,int nVertDat, int nCellDat);

	tecplotfile(){
		fid=NULL;
		lengthLine=0;
		cout << "tecplot file object created" << endl;
	}
	~tecplotfile(){
		if (fid!=NULL){
			this->CloseFile();
			cout << "Object deleted - File Closed" << endl;
		} else {
			cout << "Object deleted - No File Closed" << endl;
		}
	}

	void NewZone(){
		fprintf(fid, "ZONE\n");
		nZones++;
	}
	void StrandTime(int strandID, double timeStep){
		fprintf(fid, "STRANDID = %i \nSOLUTIONTIME = %.10lf \n",strandID,timeStep);
	}
	
	int Print(const char *format, ...){ // Mimics the printf function hides fid

		va_list arg;
		int i;

		va_start (arg, format);
		i = vfprintf (fid, format, arg);
		va_end (arg);
		lengthLine+=i;
		if (lengthLine>25000){
			fprintf(fid, " \n");lengthLine=0;
		}
		return(lengthLine);
	}

	void ResetLine(){lengthLine=0;}

};

// Derived Classes

// functions

int Test_tecplotfile();
int TestCompareReadWrite(const char* fileToOpen, mesh &blockGrid, tecplotfile &outmesh1);
// member function definition template <class T> 



#endif // POSTPROCESSING_H_INCLUDED

