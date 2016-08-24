// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     NearesNeighbourInterpolator.hh
 *       \brief    <Description>
 *
 *       \date     Aug 8, 2016
 *       \author   kroppert
 */
//================================================================================================

#pragma once


#include "MeshBasedInterpolator.hh"
#include "DataInOut/SimInput.hh"
#include <boost/tr1/type_traits.hpp>
#include <def_use_cgal.hh>
#include <def_use_flann.hh>



#ifdef USE_CGAL
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits.h>
#include <list>
#include <cmath>

namespace CGAL
{

/** define this Point within the CGAL namespace as there is a Point.hh in CFS */
struct Point {
  double vec[3];
  double vel[3];
  Complex velZ[3];

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
      Complex vx, Complex vy, Complex vz) {
    vec[0]=x; vec[1]=y; vec[2]=z;
    velZ[0]=vx; velZ[1]=vy; velZ[2]=vz;
  }

  double x() const { return vec[ 0 ]; }
  double y() const { return vec[ 1 ]; }
  double z() const { return vec[ 2 ]; }

  void vx(double& ret) const { ret = vel[ 0 ]; }
  void vy(double& ret) const { ret = vel[ 1 ]; }
  void vz(double& ret) const { ret = vel[ 2 ]; }

  void vx(Complex& ret) const { ret = velZ[ 0 ]; }
  void vy(Complex& ret) const { ret = velZ[ 1 ]; }
  void vz(Complex& ret) const { ret = velZ[ 2 ]; }

  double& x() { return vec[ 0 ]; }
  double& y() { return vec[ 1 ]; }
  double& z() { return vec[ 2 ]; }

  bool operator==(const Point& p) const
    {
    return (x() == p.x()) && (y() == p.y()) && (z() == p.z())  ;
    }

  bool  operator!=(const Point& p) const { return ! (*this == p); }
}; //end of class

template <>
struct Kernel_traits<CGAL::Point> {
  struct Kernel {
    typedef double FT;
    typedef double RT;
  };
};
} // end of namespace CGAL


struct Construct_coord_iterator {
  typedef  const double* result_type;
  const double* operator()(const CGAL::Point& p) const
  { return static_cast<const double*>(p.vec); }

  const double* operator()(const CGAL::Point& p, int)  const
  { return static_cast<const double*>(p.vec+3); }
};

struct Distance {
  typedef CGAL::Point Query_item;
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
                      int /* cutting_dimension */)  const {
    return dist + new_off*new_off - old_off*old_off;
  }

  double transformed_distance(double d) const { return d*d; }

  double inverse_of_transformed_distance(double d) { return std::sqrt(d); }

}; // end of struct Distance


typedef CGAL::Search_traits<double, CGAL::Point, const double*, Construct_coord_iterator> Traits;
typedef CGAL::Orthogonal_k_neighbor_search<Traits, Distance> K_neighbor_search;
typedef K_neighbor_search::Tree Tree;

#endif // USE_CGAL


#ifdef USE_FLANN
#include <flann/flann.hpp>
#endif // USE_FLANN


namespace CFSDat{

//! Class for calculating a nearest neighbour interpolation using CGAL or FLANN
//! for neighbour search


class NearestNeighbourInterpolator : public MeshBasedInterpolator{

  struct InpolationStruct{
    CF::Vector<Double> localCoords;
    UInt tENum;
    UInt srcEqn;
    Double volume;

    InpolationStruct() : tENum(0),srcEqn(0),volume(.0){
      localCoords.Resize(3);
    }

    bool operator < (const InpolationStruct& str) const
    {
        return (srcEqn < str.srcEqn);
    }
  };

public:

  NearestNeighbourInterpolator(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan);

  virtual ~NearestNeighbourInterpolator();

  virtual bool Run();

  enum InterpolationAlgorithm
  {
    SHEPARD, NEAREST_NEIGHBOR
  };

/*  enum KNNLibary
  {
    CGAL, FLANN
  };
*/
protected:

  virtual void PrepareInterpolation();

  virtual ResultIdList SetUpstreamResults();

  virtual void AdaptFilterResults();

  // Read scattered data
  void ReadScatteredData(CF::StdVector< CF::Vector<Double> > elemCentroids, CF::StdVector< CF::Vector<Double> > scatteredData);

  // Coordinates of input data
  CF::StdVector< CF::Vector<double> > sourceCoords_;

  // Coordinates of target data
  CF::StdVector< CF::Vector<double> > targetCoords_;

  //! Dimension of input values (0=scalar, 1=two-dim vector, 2=three-dim vector).
  UInt inDim_;

  // Search radius for values.
  Double searchRadius_;


  //! Number of neighbor points to include in interpolation.
  UInt numNeighbors_;

  //! Exponent for calculation of interpolation weight function.
  Double p_;

  //! Library used to find the k nearest neighbors of a point.
  //KNNLibary knnLib_;
  //0 for CGAL, 1 for FLANN
  UInt knnLib_;

#ifdef USE_CGAL
    boost::shared_ptr<Tree> searchTree_;

    void KNNSearch_CGAL(const Vector<Double> globPoint,
                        StdVector< Vector<Double> >& neighbors,
                        StdVector< Double >& l2Distances,
                        StdVector< Vector<Double> >& vectors);
#endif

#ifdef USE_FLANN
    boost::shared_ptr< flann::Index<flann::L2<Double> > > index_;
    boost::shared_ptr< flann::Matrix<Double> > dataset_;

    void KNNSearch_FLANN(const Vector<Double> globPoint,
                         StdVector< Vector<Double> >& neighbors,
                         StdVector< Double >& l2Distances,
                         StdVector< Vector<Double> >& vectors,
                         StdVector< Vector<Double> > scatteredData);
#endif


private:

  std::vector<InpolationStruct> interpolData_;
  StdVector<UInt> nodeNeighbours_;

};

}
