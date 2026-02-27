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


#include <def_use_cgal.hh>
#include <def_use_nanoflann.hh>

#include "CoefFunction.hh"

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
  typedef CGAL::Point Point_d;
  typedef double FT;

  double transformed_distance(const CGAL::Point& p1, const CGAL::Point& p2) const {
    double distx= p1.x()-p2.x();
    double disty= p1.y()-p2.y();
    double distz= p1.z()-p2.z();
    return distx*distx+disty*disty+distz*distz;
  }

  template <class TreeTraits>
  double min_distance_to_rectangle(const CGAL::Point &p,
                                   const CGAL::Kd_tree_rectangle<TreeTraits> &b, std::vector<FT> &dists)
  {
    double h = p.x();
    if (h < b.min_coord(0))
      dists[0] = (b.min_coord(0) - h) * (b.min_coord(0) - h);
    if (h > b.max_coord(0))
      dists[0] = (h - b.max_coord(0)) * (h - b.max_coord(0));
    h = p.y();
    if (h < b.min_coord(1))
      dists[1] = (b.min_coord(1) - h) * (b.min_coord(1) - h);
    if (h > b.max_coord(1))
      dists[1] = (h - b.max_coord(1)) * (h - b.min_coord(1));
    h = p.z();
    if (h < b.min_coord(2))
      dists[2] = (b.min_coord(2) - h) * (b.min_coord(2) - h);
    if (h > b.max_coord(2))
      dists[2] = (h - b.max_coord(2)) * (h - b.max_coord(2));
    return dists[0] + dists[1] + dists[2];
  }

  template <class TreeTraits>
  double max_distance_to_rectangle(const CGAL::Point &p,
                                   const CGAL::Kd_tree_rectangle<TreeTraits> &b,
                                   std::vector<FT> &dists)
  {
    double h = p.x();

    dists[0] = (h >= (b.min_coord(0) + b.max_coord(0)) / 2.0) ? (h - b.min_coord(0)) * (h - b.min_coord(0)) : (b.max_coord(0) - h) * (b.max_coord(0) - h);

    h = p.y();
    dists[1] = (h >= (b.min_coord(1) + b.max_coord(1)) / 2.0) ? (h - b.min_coord(1)) * (h - b.min_coord(1)) : (b.max_coord(1) - h) * (b.max_coord(1) - h);
    h = p.z();
    dists[2] = (h >= (b.min_coord(2) + b.max_coord(2)) / 2.0) ? (h - b.min_coord(2)) * (h - b.min_coord(2)) : (b.max_coord(2) - h) * (b.max_coord(2) - h);

    return dists[0] + dists[1] + dists[2];
  }

  template <class TreeTraits>
  double min_distance_to_rectangle(const CGAL::Point &p,
                                   const CGAL::Kd_tree_rectangle<TreeTraits> &b) const
  {
    std::vector<FT> dists(3, 0);
    return min_distance_to_rectangle(p, b, dists);
  }

  template <class TreeTraits>
  double max_distance_to_rectangle(const CGAL::Point &p,
                                   const CGAL::Kd_tree_rectangle<TreeTraits> &b) const
  {
    std::vector<FT> dists(3, 0);
    return max_distance_to_rectangle(p, b, dists);
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

#ifdef USE_NANOFLANN
#include <nanoflann.hpp>
#endif // USE_NANOFLANN


namespace CoupledField {

  //! Interpolation of scattered data
  template<typename T, UInt DOFS = 2>
  class CoefFunctionScatteredData : public CoefFunction
  {
  public:
    enum InterpolationAlgorithm
    {
      SHEPARD, NEAREST_NEIGHBOR
    };


    enum KNNLibary
    {
      CGAL, NANOFLANN
    };
    
    //! Constructor
    CoefFunctionScatteredData(PtrParamNode& scatteredDataNode);
    //! Destructor
    virtual ~CoefFunctionScatteredData(){;}
    
    virtual string GetName() const { return "CoefFunctionScatteredData"; }

    //! Return vector value at integration point
    virtual void GetVector( Vector<T>& vec, 
                            const LocPointMapped& lpm );
    
    //! Return scalar value at Integration point
    virtual void GetScalar( T & value,
                               const LocPointMapped& lpm );

    //! Return size of vector in case coefficient function is a vector
    virtual UInt GetVecSize() const { return DOFS; }

    //! Dump coefficient function to string 
    virtual std::string ToString() const;
    
    //! \copydoc CoefFunction::SetDerivativeOperation
    virtual void SetDerivativeOperation(CoefDerivativeType type){
      this->derivType_ = type;

      //make some checks here!
      switch(dimType_){
      case SCALAR:
        //only NONE is valid right now
        //if extended to gradient, this would be fine too
        if(type==VECTOR_DIVERGENCE){
          EXCEPTION("CoefFunctionScatteredData: VECTOR_DIVERGENCE is not a valid operator for scalar coefFunction");
        }
        break;
      case VECTOR:
        //this is fine in all cases right now
        if(type==VECTOR_DIVERGENCE){
          //change dim type to scalar
          this->dimType_ = SCALAR;
          WARN("Scattered Data Interpolation is currently very sensitive when it comes to divergences! Use with special care!")
        }
        break;
      case TENSOR:
        if(type==VECTOR_DIVERGENCE){
          EXCEPTION("CoefFunctionScatteredData: VECTOR_DIVERGENCE is not a valid operator for tensor coefFunction");
        }
        break;
      default:
        break;
      }
      return;
    }

  protected:

    void InterpolateVector(Vector<Double> globPoint, Vector<T> & vec);

    void Read(bool updateMode);
    void GetQuantityData(bool updateMode);
    void DumpData();
    
    std::vector< std::vector<double> > coordinates_;
    std::vector< std::vector<T> > scatteredData_;  // CHANGED

    std::string qid_;
          
    // Scale factor for values.
    Double factor_;

    // Search radius for values.
    Double searchRadius_;

    //! Type of interpolation algorithm.
    InterpolationAlgorithm interpolAlgo_;

    //! Number of neighbor points to include in interpolation.
    UInt numNeighbors_;

    //! Exponent for calculation of interpolation weight function.
    Double p_;

    //! Library used to find the k nearest neighbors of a point.
    KNNLibary knnLib_;

    PtrParamNode quantityNode_;


#ifdef USE_CGAL
    boost::shared_ptr<Tree> searchTree_;

    void KNNSearch_CGAL(const Vector<Double> globPoint,
                        StdVector< Vector<Double> >& neighbors,
                        StdVector< Double >& l2Distances,
                        StdVector< Vector<T> >& vectors);
#endif

#ifdef USE_NANOFLANN
    // the adaptor provides the interface between our own coordinates_ and nanoflann -> no copying! 
    struct PointCloudAdaptor {
      const std::vector<std::vector<double>>& pts;
      PointCloudAdaptor(const std::vector<std::vector<double>>& pts_) : pts(pts_) {}
      size_t kdtree_get_point_count() const { return pts.size(); }
      double kdtree_get_pt(const size_t idx, const size_t dim) const { return pts[idx][dim]; }
      template <class BBOX>
      bool kdtree_get_bbox(BBOX&) const { return false; } // let nanoflann compute it
    };

    using NanoFlannIndex = nanoflann::KDTreeSingleIndexAdaptor<
        nanoflann::L2_Simple_Adaptor<Double, PointCloudAdaptor>,
        PointCloudAdaptor, 3>;

    shared_ptr<PointCloudAdaptor> adaptor_;
    shared_ptr<NanoFlannIndex>    index_;

    void KNNSearch_FLANN(const Vector<Double> globPoint,
                         StdVector< Vector<Double> >& neighbors,
                         StdVector< Double >& l2Distances,
                         StdVector< Vector<T> >& vectors);
#endif

  };
}

#endif
