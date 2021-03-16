#ifndef RSVSMATH_AUTOMATIC_H_INCLUDED
#define RSVSMATH_AUTOMATIC_H_INCLUDED

#include <cmath>
#include <vector>

#include "vectorarray.hpp"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
static double rsvsmath_automatic_eps_edge = 0.0;
static double rsvsmath_automatic_eps_surf = 1e-10;
static double rsvsmath_automatic_eps_volu = 0.0;
static double rsvsmath_automatic_eps_centre = 1e-10;
static double rsvsmath_automatic_eps_centre2 = 1e-10;
#pragma GCC diagnostic pop

void Volume_f(const std::vector<double> &p0, const std::vector<double> &p1, const std::vector<double> &p2, double &t0);
void Volume_df(const std::vector<double> &p0, const std::vector<double> &p1, const std::vector<double> &p2,
               ArrayVec<double> &A0);
void Volume_ddf(const std::vector<double> &p0, const std::vector<double> &p1, const std::vector<double> &p2,
                ArrayVec<double> &A0);

void Volume2_f(double d0, double d1, double d2, const std::vector<double> &g0s, const std::vector<double> &g1s,
               const std::vector<double> &g2s, const std::vector<double> &g0e, const std::vector<double> &g1e,
               const std::vector<double> &g2e, double &t0);
void Volume2_df(double d0, double d1, double d2, const std::vector<double> &g0s, const std::vector<double> &g1s,
                const std::vector<double> &g2s, const std::vector<double> &g0e, const std::vector<double> &g1e,
                const std::vector<double> &g2e, ArrayVec<double> &A0);
void Volume2_ddf(double d0, double d1, double d2, const std::vector<double> &g0s, const std::vector<double> &g1s,
                 const std::vector<double> &g2s, const std::vector<double> &g0e, const std::vector<double> &g1e,
                 const std::vector<double> &g2e, ArrayVec<double> &A0);

void Area_f(const std::vector<double> &p0, const std::vector<double> &p1, const std::vector<double> &p2, double &t0,
            double eps = 0.0);
void Area_df(const std::vector<double> &p0, const std::vector<double> &p1, const std::vector<double> &p2,
             ArrayVec<double> &A0, double eps = rsvsmath_automatic_eps_surf);
void Area_ddf(const std::vector<double> &p0, const std::vector<double> &p1, const std::vector<double> &p2,
              ArrayVec<double> &A0, double eps = rsvsmath_automatic_eps_surf);
void LengthEdge_f(const std::vector<double> &p0, const std::vector<double> &p1, double &t0);
void LengthEdge_df(const std::vector<double> &p0, const std::vector<double> &p1, ArrayVec<double> &A0);
void LengthEdge_ddf(const std::vector<double> &p0, const std::vector<double> &p1, ArrayVec<double> &A0);
void SurfIntersect_f(const std::vector<double> &p1, const std::vector<double> &p2, const std::vector<double> &v0,
                     const std::vector<double> &v01, const std::vector<double> &v02, const std::vector<double> &v11,
                     const std::vector<double> &v12, const std::vector<double> &c, double &t0);
void SurfIntersect_df(const std::vector<double> &p1, const std::vector<double> &p2, const std::vector<double> &v0,
                      const std::vector<double> &v01, const std::vector<double> &v02, const std::vector<double> &v11,
                      const std::vector<double> &v12, const std::vector<double> &c, ArrayVec<double> &A0);
void SurfIntersect_ddf(const std::vector<double> &p1, const std::vector<double> &p2, const std::vector<double> &v0,
                       const std::vector<double> &v01, const std::vector<double> &v02, const std::vector<double> &v11,
                       const std::vector<double> &v12, const std::vector<double> &c, ArrayVec<double> &A0);
void SurfCentroid4_f(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z,
                     const double totD, const double X_dot_d, const double Y_dot_d, const double Z_dot_d,
                     ArrayVec<double> &A0);
void SurfCentroid4_df(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z,
                      const double totD, const double X_dot_d, const double Y_dot_d, const double Z_dot_d,
                      ArrayVec<double> &A0);
void SurfCentroid5_f(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z,
                     const double totD, const double X_dot_d, const double Y_dot_d, const double Z_dot_d,
                     ArrayVec<double> &A0);
void SurfCentroid5_df(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z,
                      const double totD, const double X_dot_d, const double Y_dot_d, const double Z_dot_d,
                      ArrayVec<double> &A0);
void SurfCentroid6_f(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z,
                     const double totD, const double X_dot_d, const double Y_dot_d, const double Z_dot_d,
                     ArrayVec<double> &A0);
void SurfCentroid6_df(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z,
                      const double totD, const double X_dot_d, const double Y_dot_d, const double Z_dot_d,
                      ArrayVec<double> &A0);
void SurfCentroidConnec_f(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z,
                          const double totD, const double X_dot_d, const double Y_dot_d, const double Z_dot_d,
                          ArrayVec<double> &A0);
void SurfCentroidConnec_df(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z,
                           const double totD, const double X_dot_d, const double Y_dot_d, const double Z_dot_d,
                           ArrayVec<double> &A0);
void SurfCentroidConnec_ddf(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z,
                            const double totD, const double X_dot_d, const double Y_dot_d, const double Z_dot_d,
                            ArrayVec<double> &A0);
void SurfCentroidNoConnec_f(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z,
                            const double totD, const double X_dot_d, const double Y_dot_d, const double Z_dot_d,
                            ArrayVec<double> &A0);
void SurfCentroidNoConnec_df(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z,
                             const double totD, const double X_dot_d, const double Y_dot_d, const double Z_dot_d,
                             ArrayVec<double> &A0);
void SurfCentroidNoConnec_ddf(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z,
                              const double totD, const double X_dot_d, const double Y_dot_d, const double Z_dot_d,
                              ArrayVec<double> &A0);
void SurfCentroidSelf_f(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z,
                        const double totD, const double X_dot_d, const double Y_dot_d, const double Z_dot_d,
                        ArrayVec<double> &A0);
void SurfCentroidSelf_df(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z,
                         const double totD, const double X_dot_d, const double Y_dot_d, const double Z_dot_d,
                         ArrayVec<double> &A0);
void SurfCentroidSelf_ddf(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z,
                          const double totD, const double X_dot_d, const double Y_dot_d, const double Z_dot_d,
                          ArrayVec<double> &A0);
#endif
