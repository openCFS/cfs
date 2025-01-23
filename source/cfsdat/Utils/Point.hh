// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     Point.hh
 *       \brief    <Description>
 *
 *       \date     Sep 10, 2016
 *       \author   kroppert
 */
//================================================================================================
#pragma once

//Definition of a point for kd-search-tree CGAL of FLANN

#include "cfsdat/Utils/Defines.hh"
#include <def_use_cgal.hh>
#include <def_use_flann.hh>
#ifdef USE_CGAL
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits.h>
#endif
#include <list>
#include <cmath>

namespace CGAL
{

struct Point {
  double vec[3];
  double vel[3];
  CoupledField::Complex velZ[3];

  Point() {
    vec[0]= vec[1] = vec[2] = 0;
    vel[0]= vel[1] = vel[2] = 0;
  }
  Point (double x, double y, double z,
      double vx, double vy, double vz) {
    vec[0]=x; vec[1]=y; vec[2]=z;
    vel[0]=vx; vel[1]=vy; vel[2]=vz;
  }

  Point (double x, double y, double z,
      CoupledField::Complex vx, CoupledField::Complex vy, CoupledField::Complex vz) {
    vec[0]=x; vec[1]=y; vec[2]=z;
    velZ[0]=vx; velZ[1]=vy; velZ[2]=vz;
  }

  double x() const { return vec[ 0 ]; }
  double y() const { return vec[ 1 ]; }
  double z() const { return vec[ 2 ]; }

  void vx(double& ret) const { ret = vel[ 0 ]; }
  void vy(double& ret) const { ret = vel[ 1 ]; }
  void vz(double& ret) const { ret = vel[ 2 ]; }

  void vx(CoupledField::Complex& ret) const { ret = velZ[ 0 ]; }
  void vy(CoupledField::Complex& ret) const { ret = velZ[ 1 ]; }
  void vz(CoupledField::Complex& ret) const { ret = velZ[ 2 ]; }

  double& x() { return vec[ 0 ]; }
  double& y() { return vec[ 1 ]; }
  double& z() { return vec[ 2 ]; }

  bool operator==(const Point& p) const
    {
    return (x() == p.x()) && (y() == p.y()) && (z() == p.z())  ;
    }

  bool  operator!=(const Point& p) const { return ! (*this == p); }
}; //end of class

#ifdef USE_CGAL
template <>
struct Kernel_traits<CGAL::Point> {
  struct Kernel {
    typedef double FT;
    typedef double RT;
  };
};
#endif
} // end of namespace CGAL

#ifdef USE_CGAL
struct Construct_coord_iterator {
  typedef  const double* result_type;
  const double* operator()(const CGAL::Point& p) const
  { return static_cast<const double*>(p.vec); }

  const double* operator()(const CGAL::Point& p, int)  const
  { return static_cast<const double*>(p.vec+3); }
};
#endif

#ifdef USE_CGAL
struct Distance {
  typedef CGAL::Point Query_item;
  typedef CGAL::Point Point_d;
  typedef double FT;

  double transformed_distance(const CGAL::Point& p1, const CGAL::Point& p2) const {
    double distx= p1.x()-p2.x();
    double disty= p1.y()-p2.y();
    double distz= p1.z()-p2.z();
    return distx*distx+disty*disty+distz*distz;
  }

  template <class TreeTraits>
  double min_distance_to_rectangle(const CGAL::Point& p,
                                   const CGAL::Kd_tree_rectangle<TreeTraits>& b) const {
    double distance(0.0), h = p.x();
    if (h < b.min_coord(0)) distance += (b.min_coord(0)-h)*(b.min_coord(0)-h);
    if (h > b.max_coord(0)) distance += (h-b.max_coord(0))*(h-b.max_coord(0));
    h=p.y();
    if (h < b.min_coord(1)) distance += (b.min_coord(1)-h)*(b.min_coord(1)-h);
    if (h > b.max_coord(1)) distance += (h-b.max_coord(1))*(h-b.min_coord(1));
    h=p.z();
    if (h < b.min_coord(2)) distance += (b.min_coord(2)-h)*(b.min_coord(2)-h);
    if (h > b.max_coord(2)) distance += (h-b.max_coord(2))*(h-b.max_coord(2));
    return distance;
  }

  template <class TreeTraits>
  double max_distance_to_rectangle(const CGAL::Point& p,
                                   const CGAL::Kd_tree_rectangle<TreeTraits>& b) const {
    double h = p.x();

    double d0 = (h >= (b.min_coord(0)+b.max_coord(0))/2.0) ?
                (h-b.min_coord(0))*(h-b.min_coord(0)) : (b.max_coord(0)-h)*(b.max_coord(0)-h);

    h=p.y();
    double d1 = (h >= (b.min_coord(1)+b.max_coord(1))/2.0) ?
                (h-b.min_coord(1))*(h-b.min_coord(1)) : (b.max_coord(1)-h)*(b.max_coord(1)-h);
    h=p.z();
    double d2 = (h >= (b.min_coord(2)+b.max_coord(2))/2.0) ?
                (h-b.min_coord(2))*(h-b.min_coord(2)) : (b.max_coord(2)-h)*(b.max_coord(2)-h);
    return d0 + d1 + d2;
  }

  double new_distance(double& dist, double old_off, double new_off,
                      int )  const {
    return dist + new_off*new_off - old_off*old_off;
  }

  double transformed_distance(double d) const { return d*d; }

  double inverse_of_transformed_distance(double d) { return std::sqrt(d); }

}; // end of struct Distance
#endif

#ifdef USE_CGAL
typedef CGAL::Search_traits<double, CGAL::Point, const double*, Construct_coord_iterator> Traits;
typedef CGAL::Orthogonal_k_neighbor_search<Traits, Distance> K_neighbor_search;
typedef K_neighbor_search::Tree Tree;
#endif


