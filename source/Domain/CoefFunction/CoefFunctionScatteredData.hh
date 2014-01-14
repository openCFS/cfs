//==============================================================================
/*!
 *       \file     CoefFunctionScatteredData.hh
 *       \brief    Coefficient function for (interpolation) of scattered data.
 *
 *       \date     09/12/2013
 *       \author   Simon Triebenbacher
 */
//==============================================================================

#ifndef FILE_COEFFUNCTION_SCATTEREDDATA_HH
#define FILE_COEFFUNCTION_SCATTEREDDATA_HH

#include <boost/tr1/type_traits.hpp>

#include <def_use_cgal.hh>

#ifdef USE_CGAL
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits.h>
#include <list>
#include <cmath>
#endif

#include "CoefFunction.hh"

#ifdef USE_CGAL
struct Point {
  double vec[3];
  double vel[3];

  Point() {
    vec[0]= vec[1] = vec[2] = 0;
    vel[0]= vel[1] = vel[2] = 0;
  }
  Point (double x, double y, double z,
         double vx, double vy, double vz) {
    vec[0]=x; vec[1]=y; vec[2]=z;
    vel[0]=vx; vel[1]=vy; vel[2]=vz;
  }

  double x() const { return vec[ 0 ]; }
  double y() const { return vec[ 1 ]; }
  double z() const { return vec[ 2 ]; }

  double vx() const { return vel[ 0 ]; }
  double vy() const { return vel[ 1 ]; }
  double vz() const { return vel[ 2 ]; }

  double& x() { return vec[ 0 ]; }
  double& y() { return vec[ 1 ]; }
  double& z() { return vec[ 2 ]; }

  bool operator==(const Point& p) const
  {
    return (x() == p.x()) && (y() == p.y()) && (z() == p.z())  ;
  }

  bool  operator!=(const Point& p) const { return ! (*this == p); }
}; //end of class



namespace CGAL {

  template <>
  struct Kernel_traits<Point> {
    struct Kernel {
      typedef double FT;
      typedef double RT;
    };
  };
}


struct Construct_coord_iterator {
  typedef  const double* result_type;
  const double* operator()(const Point& p) const
  { return static_cast<const double*>(p.vec); }

  const double* operator()(const Point& p, int)  const
  { return static_cast<const double*>(p.vec+3); }
};

struct Distance {
  typedef Point Query_item;
  typedef double FT;

  double transformed_distance(const Point& p1, const Point& p2) const {
    double distx= p1.x()-p2.x();
    double disty= p1.y()-p2.y();
    double distz= p1.z()-p2.z();
    return distx*distx+disty*disty+distz*distz;
  }

  template <class TreeTraits>
  double min_distance_to_rectangle(const Point& p,
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
  double max_distance_to_rectangle(const Point& p,
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


typedef CGAL::Search_traits<double, Point, const double*, Construct_coord_iterator> Traits;
typedef CGAL::Orthogonal_k_neighbor_search<Traits, Distance> K_neighbor_search;
typedef K_neighbor_search::Tree Tree;

#endif // USE_CGAL

namespace CoupledField {

  //! Interpolation of scattered data
  template<typename T, UInt DOFS = 2>
  class CoefFunctionScatteredData : public CoefFunction
  {
  public:
    
    //! Constructor
    CoefFunctionScatteredData(PtrParamNode& scatteredDataNode);
    //! Destructor
    virtual ~CoefFunctionScatteredData(){;}
    
    //! Return vector value at integration point
    virtual void GetVector( Vector<T>& vec, 
                            const LocPointMapped& lpm );
    

    //! Return size of vector in case coefficient function is a vector
    virtual UInt GetVecSize() const { return DOFS; }

    //! Dump coefficient function to string 
    virtual std::string ToString() const;
    
  protected:

    void Read();
    
    std::vector< std::vector<double> > scatteredData_;

    std::string fileName_;
          
    std::map<UInt, UInt> dof2CoordColumn_;
    std::map<UInt, UInt> dof2ValueColumn_;

    Double factor_;

    std::string format_;
    
#ifdef USE_CGAL
    boost::shared_ptr<Tree> searchTree_;
#endif
  };
}

#endif
