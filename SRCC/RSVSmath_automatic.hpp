#ifndef RSVSMATH_AUTOMATIC_H_INCLUDED 
#define RSVSMATH_AUTOMATIC_H_INCLUDED 

#include <vector> 
#include <cmath> 
#include "vectorarray.hpp" 
using namespace std; 

void Volume_f(const vector<double>& p0 , const vector<double>& p1 , const vector<double>& p2 , double &   t0 ); 
void Volume_df(const vector<double>& p0 , const vector<double>& p1 , const vector<double>& p2 , ArrayVec<double> &   A0 ); 
void Volume_ddf(const vector<double>& p0 , const vector<double>& p1 , const vector<double>& p2 , ArrayVec<double> &   A0 ); 
void Area_f(const vector<double>& p0 , const vector<double>& p1 , const vector<double>& p2 , double &   t0 ); 
void Area_df(const vector<double>& p0 , const vector<double>& p1 , const vector<double>& p2 , ArrayVec<double> &   A0 ); 
void Area_ddf(const vector<double>& p0 , const vector<double>& p1 , const vector<double>& p2 , ArrayVec<double> &   A0 ); 
void LengthEdge_f(const vector<double>& p0 , const vector<double>& p1 , double &   t0 ); 
void LengthEdge_df(const vector<double>& p0 , const vector<double>& p1 , ArrayVec<double> &   A0 ); 
void LengthEdge_ddf(const vector<double>& p0 , const vector<double>& p1 , ArrayVec<double> &   A0 ); 
void SurfIntersect_f(const vector<double>& p1 , const vector<double>& p2 , const vector<double>& v0 , const vector<double>& v01 , const vector<double>& v02 , const vector<double>& v11 , const vector<double>& v12 , const vector<double>& c , double &   t0 ); 
void SurfIntersect_df(const vector<double>& p1 , const vector<double>& p2 , const vector<double>& v0 , const vector<double>& v01 , const vector<double>& v02 , const vector<double>& v11 , const vector<double>& v12 , const vector<double>& c , ArrayVec<double> &   A0 ); 
void SurfIntersect_ddf(const vector<double>& p1 , const vector<double>& p2 , const vector<double>& v0 , const vector<double>& v01 , const vector<double>& v02 , const vector<double>& v11 , const vector<double>& v12 , const vector<double>& c , ArrayVec<double> &   A0 ); 
void SurfCentroid4_f(const vector<double>& x , const vector<double>& y , const vector<double>& z , ArrayVec<double> &   A0 ); 
void SurfCentroid4_df(const vector<double>& x , const vector<double>& y , const vector<double>& z , ArrayVec<double> &   A0 ); 
void SurfCentroidDiv_f(const vector<double>& x , const vector<double>& y , const vector<double>& z , double &   t0 ); 
void SurfCentroidDiv_df(const vector<double>& x , const vector<double>& y , const vector<double>& z , ArrayVec<double> &   A0 ); 
void SurfCentroid5_f(const vector<double>& x , const vector<double>& y , const vector<double>& z , ArrayVec<double> &   A0 ); 
void SurfCentroid5_df(const vector<double>& x , const vector<double>& y , const vector<double>& z , ArrayVec<double> &   A0 ); 
void SurfCentroid6_f(const vector<double>& x , const vector<double>& y , const vector<double>& z , ArrayVec<double> &   A0 ); 
void SurfCentroid6_df(const vector<double>& x , const vector<double>& y , const vector<double>& z , ArrayVec<double> &   A0 ); 
void SurfCentroidConnec_f(const vector<double>& x , const vector<double>& y , const vector<double>& z , ArrayVec<double> &   A0 ); 
void SurfCentroidConnec_df(const vector<double>& x , const vector<double>& y , const vector<double>& z , ArrayVec<double> &   A0 ); 
void SurfCentroidNoConnec_f(const vector<double>& x , const vector<double>& y , const vector<double>& z , ArrayVec<double> &   A0 ); 
void SurfCentroidNoConnec_df(const vector<double>& x , const vector<double>& y , const vector<double>& z , ArrayVec<double> &   A0 ); 
#endif 
