#ifndef RSVSMATH_AUTOMATIC_H_INCLUDED 
#define RSVSMATH_AUTOMATIC_H_INCLUDED 

#include <vector> 
#include <cmath> 
#include "vectorarray.hpp" 
using namespace std; 

static const double rsvsmath_automatic_eps_edge = 1e-15;
static const double rsvsmath_automatic_eps_surf = 1e-15;
static const double rsvsmath_automatic_eps_volu = 1e-15;

void Volume_f(const vector<double>& p0 , const vector<double>& p1 , const vector<double>& p2 , double &   t0 ); 
void Volume_df(const vector<double>& p0 , const vector<double>& p1 , const vector<double>& p2 , ArrayVec<double> &   A0 ); 
void Volume_ddf(const vector<double>& p0 , const vector<double>& p1 , const vector<double>& p2 , ArrayVec<double> &   A0 ); 

void Volume2_f(double  d0 , double  d1 , double  d2 , const vector<double>& g0s , const vector<double>& g1s , const vector<double>& g2s , const vector<double>& g0e , const vector<double>& g1e , const vector<double>& g2e , double &   t0 ) ;
void Volume2_df(double  d0 , double  d1 , double  d2 , const vector<double>& g0s , const vector<double>& g1s , const vector<double>& g2s , const vector<double>& g0e , const vector<double>& g1e , const vector<double>& g2e , ArrayVec<double> &   A0 );
void Volume2_ddf(double  d0 , double  d1 , double  d2 , const vector<double>& g0s , const vector<double>& g1s , const vector<double>& g2s , const vector<double>& g0e , const vector<double>& g1e , const vector<double>& g2e , ArrayVec<double> &   A0 ) ;


void Area_f(const vector<double>& p0 , const vector<double>& p1 , const vector<double>& p2 , double &   t0 ); 
void Area_df(const vector<double>& p0 , const vector<double>& p1 , const vector<double>& p2 , ArrayVec<double> &   A0 ); 
void Area_ddf(const vector<double>& p0 , const vector<double>& p1 , const vector<double>& p2 , ArrayVec<double> &   A0 ); 
void LengthEdge_f(const vector<double>& p0 , const vector<double>& p1 , double &   t0 ); 
void LengthEdge_df(const vector<double>& p0 , const vector<double>& p1 , ArrayVec<double> &   A0 ); 
void LengthEdge_ddf(const vector<double>& p0 , const vector<double>& p1 , ArrayVec<double> &   A0 ); 
void SurfIntersect_f(const vector<double>& p1 , const vector<double>& p2 , const vector<double>& v0 ,
	const vector<double>& v01 , const vector<double>& v02 , const vector<double>& v11 , const vector<double>& v12 ,
	 const vector<double>& c , double &   t0 ); 
void SurfIntersect_df(const vector<double>& p1 , const vector<double>& p2 , const vector<double>& v0 ,
 const vector<double>& v01 , const vector<double>& v02 , const vector<double>& v11 , const vector<double>& v12 ,
  const vector<double>& c , ArrayVec<double> &   A0 ); 
void SurfIntersect_ddf(const vector<double>& p1 , const vector<double>& p2 , const vector<double>& v0 ,
 const vector<double>& v01 , const vector<double>& v02 , const vector<double>& v11 , const vector<double>& v12 ,
  const vector<double>& c , ArrayVec<double> &   A0 ); 
void SurfCentroid4_f(const vector<double>& x , const vector<double>& y , const vector<double>& z ,
 const double totD , const double X_dot_d , const double Y_dot_d , const double Z_dot_d , ArrayVec<double> &   A0 ); 
void SurfCentroid4_df(const vector<double>& x , const vector<double>& y , const vector<double>& z ,
 const double totD , const double X_dot_d , const double Y_dot_d , const double Z_dot_d , ArrayVec<double> &   A0 ); 
void SurfCentroid5_f(const vector<double>& x , const vector<double>& y , const vector<double>& z ,
 const double totD , const double X_dot_d , const double Y_dot_d , const double Z_dot_d , ArrayVec<double> &   A0 ); 
void SurfCentroid5_df(const vector<double>& x , const vector<double>& y , const vector<double>& z ,
 const double totD , const double X_dot_d , const double Y_dot_d , const double Z_dot_d , ArrayVec<double> &   A0 ); 
void SurfCentroid6_f(const vector<double>& x , const vector<double>& y , const vector<double>& z ,
 const double totD , const double X_dot_d , const double Y_dot_d , const double Z_dot_d , ArrayVec<double> &   A0 ); 
void SurfCentroid6_df(const vector<double>& x , const vector<double>& y , const vector<double>& z ,
 const double totD , const double X_dot_d , const double Y_dot_d , const double Z_dot_d , ArrayVec<double> &   A0 ); 
void SurfCentroidConnec_f(const vector<double>& x , const vector<double>& y , const vector<double>& z ,
 const double totD , const double X_dot_d , const double Y_dot_d , const double Z_dot_d , ArrayVec<double> &   A0 ); 
void SurfCentroidConnec_df(const vector<double>& x , const vector<double>& y , const vector<double>& z ,
 const double totD , const double X_dot_d , const double Y_dot_d , const double Z_dot_d ,
  ArrayVec<double> &   A0 ); 
void SurfCentroidConnec_ddf(const vector<double>& x , const vector<double>& y , const vector<double>& z ,
 const double totD , const double X_dot_d , const double Y_dot_d , const double Z_dot_d , ArrayVec<double> &   A0 ); 
void SurfCentroidNoConnec_f(const vector<double>& x , const vector<double>& y , const vector<double>& z ,
 const double totD , const double X_dot_d , const double Y_dot_d , const double Z_dot_d , ArrayVec<double> &   A0 ); 
void SurfCentroidNoConnec_df(const vector<double>& x , const vector<double>& y , const vector<double>& z ,
 const double totD , const double X_dot_d , const double Y_dot_d , const double Z_dot_d , ArrayVec<double> &   A0 ); 
void SurfCentroidNoConnec_ddf(const vector<double>& x , const vector<double>& y , const vector<double>& z ,
 const double totD , const double X_dot_d , const double Y_dot_d , const double Z_dot_d , ArrayVec<double> &   A0 ); 
void SurfCentroidSelf_f(const vector<double>& x , const vector<double>& y , const vector<double>& z ,
 const double totD , const double X_dot_d , const double Y_dot_d , const double Z_dot_d , ArrayVec<double> &   A0 ); 
void SurfCentroidSelf_df(const vector<double>& x , const vector<double>& y , const vector<double>& z ,
 const double totD , const double X_dot_d , const double Y_dot_d , const double Z_dot_d , ArrayVec<double> &   A0 ); 
void SurfCentroidSelf_ddf(const vector<double>& x , const vector<double>& y , const vector<double>& z ,
 const double totD , const double X_dot_d , const double Y_dot_d , const double Z_dot_d , ArrayVec<double> &   A0 ); 
#endif 
