#ifndef RSVSMATH_AUTOMATIC_H_INCLUDED 
#define RSVSMATH_AUTOMATIC_H_INCLUDED 

#include <vector> 
#include "vectorarray.hpp" 

using namespace std; 

void Volume_f(vector<double> p0 , vector<double> p1 , vector<double> p2 , double &   t0 ); 
void Volume_df(vector<double> p0 , vector<double> p1 , vector<double> p2 , ArrayVec<double> &   A0 ); 
void Volume_ddf(vector<double> p0 , vector<double> p1 , vector<double> p2 , ArrayVec<double> &   A0 ); 
void Area_f(vector<double> p0 , vector<double> p1 , vector<double> p2 , double &   t0 ); 
void Area_df(vector<double> p0 , vector<double> p1 , vector<double> p2 , ArrayVec<double> &   A0 ); 
void Area_ddf(vector<double> p0 , vector<double> p1 , vector<double> p2 , ArrayVec<double> &   A0 ); 
void LengthEdge_f(vector<double> p0 , vector<double> p1 , double &   t0 ); 
void LengthEdge_df(vector<double> p0 , vector<double> p1 , ArrayVec<double> &   A0 ); 
void LengthEdge_ddf(vector<double> p0 , vector<double> p1 , ArrayVec<double> &   A0 ); 
void SurfIntersect_f(vector<double> p1 , vector<double> p2 , vector<double> v0 , vector<double> v01 , vector<double> v02 , vector<double> v11 , vector<double> v12 , vector<double> c , double &   t0 ); 
void SurfIntersect_df(vector<double> p1 , vector<double> p2 , vector<double> v0 , vector<double> v01 , vector<double> v02 , vector<double> v11 , vector<double> v12 , vector<double> c , ArrayVec<double> &   A0 ); 
void SurfIntersect_ddf(vector<double> p1 , vector<double> p2 , vector<double> v0 , vector<double> v01 , vector<double> v02 , vector<double> v11 , vector<double> v12 , vector<double> c , ArrayVec<double> &   A0 ); 
void SurfCentroid_f(vector<double> x , vector<double> y , vector<double> z , ArrayVec<double> &   A0 ); 
void SurfCentroid_df(vector<double> x , vector<double> y , vector<double> z , ArrayVec<double> &   A0 ); 
void SurfCentroidDiv_f(vector<double> x , vector<double> y , vector<double> z , double &   t0 ); 
void SurfCentroidDiv_df(vector<double> x , vector<double> y , vector<double> z , ArrayVec<double> &   A0 ); 
void SurfCentroidDiv_f(vector<double> x , vector<double> y , vector<double> z , ArrayVec<double> &   A0 ); 
void SurfCentroidDiv_df(vector<double> x , vector<double> y , vector<double> z , ArrayVec<double> &   A0 ); 

#endif