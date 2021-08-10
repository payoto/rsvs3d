/**
 * Provide tecplot file formating for mesh and snake outputs.
 *
 *@file
 */

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
class snake;
class RSVScalc;
class triangulation;
class trisurfarray;
class triarray;

//=================================
// included dependencies
#include <stdarg.h>

#include <iostream>

#include "arraystructures.hpp"

//==================================
// Code
// NOTE: function in a class definition are IMPLICITELY INLINED
//       ie replaced by their code at compile time

namespace rsvs3d
{
namespace constants
{
namespace tecplot
{
static const int autoselect = 0;
static const int polyhedron = 1;
static const int polygon = 2;
static const int line = 3;
static const int point = 4;

static const int nosharedzone = -1;
static const auto __issharedzone = [](int z) -> bool { return z > 0; };
namespace snakedata
{
static const std::string snaxel = "snaxel";
static const std::string normal = "normal";
static const std::string laplacian = "laplacian";
static const std::string direction = "direction";
static const std::string __default = snaxel;
} // namespace snakedata
} // namespace tecplot
} // namespace constants
} // namespace rsvs3d

// Base classes

class tecplotfile
{
  private:
    HashedVector<int, int> strands;
    FILE *fid;
    int lengthLine;
    int nZones = 0;
    bool isloud;

  public:
    bool isStrand(int strandID)
    {
        return rsvs3d::logicals::__isfound(this->strands.find(strandID));
    }
    int OpenFile(const char *str, const char *mode = "w");
    void CloseFile();
    int ZoneNum() const
    {
        return (nZones);
    }
    // Mesh out
    int PrintMesh(const mesh &meshout, int strandID = 0, double timeStep = 0,
                  int forceOutType = rsvs3d::constants::tecplot::autoselect, const std::vector<int> &vertList = {});
    int PrintSnakeInternalPts(const snake &snakein, int strandID = 0, double timeStep = 0);
    int VolDataBlock(const mesh &meshout, int nVert, int nVolu, int nVertDat, const std::vector<int> &voluList = {},
                     const std::vector<int> &vertList = {});
    int SurfDataBlock(const mesh &meshout, int nVert, int nSurf, int nVertDat);
    int LineDataBlock(const mesh &meshout, int nVert, int nEdge, int nVertDat, int nCellDat);
    int VertDataBlock(const mesh &meshout, int nVert, int nVertDat, int nCellDat,
                      const std::vector<int> &vertList = {});
    int VolFaceMap(const mesh &meshout, int nSurf);
    int VolFaceMap(const mesh &meshout, const std::vector<int> &surfList, const std::vector<int> &voluList,
                   const std::vector<int> &vertList);
    int SurfFaceMap(const mesh &meshout, int nEdge);
    int LineFaceMap(const mesh &meshout, int nEdge);
    // Partial Mesh Out with variable sharing
    int PrintVolumeDat(const mesh &meshout, int shareZone, int strandID, double timeStep);
    int DefShareZoneVolume(int shareZone, int nVertDat);

    // Triangulation out
    int VolDataBlock(const triangulation &triout, triarray triangulation::*mp, int nVert, int nVolu, int nVertDat);
    int SurfDataBlock(const triangulation &triout, triarray triangulation::*mp, int nVert, int nSurf, int nVertDat);
    int LineDataBlock(const triangulation &triout, triarray triangulation::*mp, int nVert, int nEdge, int nVertDat,
                      int nCellDat);
    int LineDataBlock(const triangulation &triout, triarray triangulation::*mp, int nVert, int nEdge, int nVertDat,
                      int nCellDat, const std::vector<int> &triList);
    int SurfFaceMap(const triangulation &triout, triarray triangulation::*mp);
    int LineFaceMap(const triangulation &triout, triarray triangulation::*mp);
    int LineFaceMap(const std::vector<int> &triList);
    int VolFaceMap(const triangulation &triout, triarray triangulation::*mp, int nSurf);
    int PrintTriangulation(const triangulation &triout, triarray triangulation::*mp, int strandID = 0,
                           double timeStep = 0, int forceOutType = rsvs3d::constants::tecplot::autoselect,
                           const std::vector<int> &triList = {});

    // Triangulation surface array out
    int VolDataBlock(const triangulation &triout, trisurfarray triangulation::*mp, int nVert, int nVolu, int nVertDat);
    int SurfDataBlock(const triangulation &triout, trisurfarray triangulation::*mp, int nVert, int nSurf, int nVertDat);
    int LineDataBlock(const triangulation &triout, trisurfarray triangulation::*mp, int nVert, int nEdge, int nVertDat,
                      int nCellDat);
    int SurfFaceMap(const triangulation &triout, trisurfarray triangulation::*mp);
    int LineFaceMap(const triangulation &triout, trisurfarray triangulation::*mp);
    int VolFaceMap(const triangulation &triout, trisurfarray triangulation::*mp, int nSurf);
    int PrintTriangulation(const triangulation &triout, trisurfarray triangulation::*mp, int strandID = 0,
                           double timeStep = 0, int forceOutType = rsvs3d::constants::tecplot::autoselect);

    // Snake specific functions
    int SnakeDataBlock(const snake &snakeout, int nVert, int nVertDat,
                       std::string = rsvs3d::constants::tecplot::snakedata::__default, bool printCoord = true);
    int PrintSnake(const snake &snakeout, int strandID = 0, double timeStep = 0,
                   int forceOutType = rsvs3d::constants::tecplot::autoselect, const std::vector<int> &vertList = {});
    int PrintSnake(std::string snakeData, const snake &snakeout, int strandID = 0, double timeStep = 0,
                   int forceOutType = rsvs3d::constants::tecplot::autoselect,
                   int coordConnShareZone = rsvs3d::constants::tecplot::nosharedzone,
                   const std::vector<int> &vertList = {});
    // Zone Headers
    void ZoneHeaderPolyhedron(int nVert, int nVolu, int nSurf, int totNumFaceNode, int nVertDat, int nCellDat);
    void ZoneHeaderPolygon(int nVert, int nEdge, int nSurf, int nVertDat, int nCellDat);
    void ZoneHeaderFelineseg(int nVert, int nEdge, int nVertDat, int nCellDat);
    void ZoneHeaderOrdered(int nVert, int nVertDat, int nCellDat, int nSensDat = 0);

    void ZoneHeaderPolyhedronSnake(int nVert, int nVolu, int nSurf, int totNumFaceNode, int nVertDat, int nCellDat,
                                   int nSensDat = 0);
    void ZoneHeaderPolygonSnake(int nVert, int nEdge, int nSurf, int nVertDat, int nCellDat, int nSensDat = 0);
    void ZoneHeaderFelinesegSnake(int nVert, int nEdge, int nVertDat, int nCellDat, int nSensDat = 0);

    // Snake sensitivity
    int PrintSnakeSensitivity(const triangulation &triRSVS, const RSVScalc &calcObj, int strandID = 0,
                              double timeStep = 0, int forceOutType = rsvs3d::constants::tecplot::autoselect,
                              const std::vector<int> &vertList = {});
    int RSVScalcDataBlock(const triangulation &triRSVS, const RSVScalc &calcObj, int nVert, int nSensDat,
                          int sensStart = 0, int methodProcess = 1);
    int RSVScalcVectorDataBlock(const triangulation &triRSVS, const RSVScalc &calcObj, int nVert, int numConstrPlot = 0,
                                int methodProcess = 1);
    int PrintSnakeSensitivityTime(const triangulation &triRSVS, const RSVScalc &calcObj, int strandID = 0,
                                  double timeStep = 0, int forceOutType = rsvs3d::constants::tecplot::autoselect,
                                  const std::vector<int> &vertList = {});
    int PrintSnakeSensitivityVector(const triangulation &triRSVS, const RSVScalc &calcObj, int strandID = 0,
                                    double timeStep = 0, int forceOutType = rsvs3d::constants::tecplot::autoselect,
                                    const std::vector<int> &vertList = {});
    int PrintSnakeGradients(const triangulation &triRSVS, const RSVScalc &calcObj, int strandID = 0,
                            double timeStep = 0, int forceOutType = rsvs3d::constants::tecplot::autoselect,
                            const std::vector<int> &vertList = {});

    tecplotfile()
    {
        this->fid = NULL;
        this->lengthLine = 0;
        this->isloud = false;
    }
    tecplotfile(bool isloudIn)
    {
        fid = NULL;
        lengthLine = 0;
        this->isloud = isloudIn;
        if (this->isloud)
        {
            std::cout << "tecplot file object created" << std::endl;
        }
    }
    ~tecplotfile()
    {
        if (fid != NULL)
        {
            this->CloseFile();
            if (this->isloud)
            {
                std::cout << "Object deleted - File Closed" << std::endl;
            }
        }
        else
        {
            if (this->isloud)
            {
                std::cout << "Object deleted - No File Closed" << std::endl;
            }
        }
    }

    void NewZone()
    {
        fprintf(fid, "ZONE\n");
        nZones++;
    }
    void StrandTime(int strandID, double timeStep)
    {
        fprintf(fid, "STRANDID = %i \nSOLUTIONTIME = %.10lf \n", strandID, timeStep);
        if (!rsvs3d::logicals::__isfound(this->strands.find(strandID)))
        {
            this->strands.push_back(strandID);
        }
    }

    int Print(const char *format, ...)
    { // Mimics the printf function hides fid

        va_list arg;
        int i;

        va_start(arg, format);
        i = vfprintf(fid, format, arg);
        va_end(arg);
        lengthLine += i;
        if (lengthLine > 25000)
        {
            fprintf(fid, " \n");
            lengthLine = 0;
        }
        return (lengthLine);
    }
    void ResetLine()
    {
        lengthLine = 0;
    }
    void NewLine()
    {
        fprintf(this->fid, "\n");
        this->ResetLine();
    };
};

// Derived Classes

// functions

namespace dataoutput
{
void Coord(tecplotfile &tecout, const mesh &meshout, int nVert, int nVertDat);
void Snaxel(tecplotfile &tecout, const snake &snakeout, int nVert);
void VertexNormal(tecplotfile &tecout, const mesh &meshin, int nVert);
void VertexLaplacian(tecplotfile &tecout, const mesh &meshin, int nVert);
} // namespace dataoutput

int Test_tecplotfile();
int TestCompareReadWrite(const char *fileToOpen, mesh &blockGrid, tecplotfile &outmesh1);
// member function definition template <class T>

#endif // POSTPROCESSING_H_INCLUDED
