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
#include <def_use_flann.hh>

#include "CoefFunction.hh"

#ifdef USE_CGAL
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits.h>
#include <list>
#include <cmath>

#ifdef USE_ADOLC
#include "General/Exception.hh"
#endif

namespace CGAL
{

  /** define this Point within the CGAL namespace as there is a Point.hh in CFS */
  struct Point {
    Double vec[3];
    Double vel[3];
    Complex velZ[3];

    Point() {
      vec[0]= vec[1] = vec[2] = 0;
      vel[0]= vel[1] = vel[2] = 0;
    }
    Point (Double x, Double y, Double z,
        Double vx, Double vy, Double vz) {
      vec[0]=x; vec[1]=y; vec[2]=z;
      vel[0]=vx; vel[1]=vy; vel[2]=vz;
    }

    Point (Double x, Double y, Double z,
        Complex vx, Complex vy, Complex vz) {
      vec[0]=x; vec[1]=y; vec[2]=z;
      velZ[0]=vx; velZ[1]=vy; velZ[2]=vz;
    }

    Double x() const { return vec[ 0 ]; }
    Double y() const { return vec[ 1 ]; }
    Double z() const { return vec[ 2 ]; }

    void vx(Double& ret) const { ret = vel[ 0 ]; }
    void vy(Double& ret) const { ret = vel[ 1 ]; }
    void vz(Double& ret) const { ret = vel[ 2 ]; }

    void vx(Complex& ret) const { ret = velZ[ 0 ]; }
    void vy(Complex& ret) const { ret = velZ[ 1 ]; }
    void vz(Complex& ret) const { ret = velZ[ 2 ]; }

    Double& x() { return vec[ 0 ]; }
    Double& y() { return vec[ 1 ]; }
    Double& z() { return vec[ 2 ]; }

    bool operator==(const Point& p) const
      {
      return (x() == p.x()) && (y() == p.y()) && (z() == p.z())  ;
      }

    bool  operator!=(const Point& p) const { return ! (*this == p); }
  }; //end of class

  template <>
  struct Kernel_traits<CGAL::Point> {
    struct Kernel {
      typedef Double FT;
      typedef Double RT;
    };
  };
} // end of namespace CGAL


struct Construct_coord_iterator {
  typedef  const Double* result_type;
  const Double* operator()(const CGAL::Point& p) const
  { return static_cast<const Double*>(p.vec); }

  const Double* operator()(const CGAL::Point& p, int)  const
  { return static_cast<const Double*>(p.vec+3); }
};

struct Distance {
  typedef CGAL::Point Query_item;
  typedef Double FT;

  Double transformed_distance(const CGAL::Point& p1, const CGAL::Point& p2) const {
    Double distx= p1.x()-p2.x();
    Double disty= p1.y()-p2.y();
    Double distz= p1.z()-p2.z();
    return distx*distx+disty*disty+distz*distz;
  }

  template <class TreeTraits>
  Double min_distance_to_rectangle(const CGAL::Point& p,
                                   const CGAL::Kd_tree_rectangle<TreeTraits>& b) const {
    Double distance(0.0), h = p.x();
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
  Double max_distance_to_rectangle(const CGAL::Point& p,
                                   const CGAL::Kd_tree_rectangle<TreeTraits>& b) const {
    Double h = p.x();

    Double d0 = (h >= (b.min_coord(0)+b.max_coord(0))/2.0) ?
                (h-b.min_coord(0))*(h-b.min_coord(0)) : (b.max_coord(0)-h)*(b.max_coord(0)-h);

    h=p.y();
    Double d1 = (h >= (b.min_coord(1)+b.max_coord(1))/2.0) ?
                (h-b.min_coord(1))*(h-b.min_coord(1)) : (b.max_coord(1)-h)*(b.max_coord(1)-h);
    h=p.z();
    Double d2 = (h >= (b.min_coord(2)+b.max_coord(2))/2.0) ?
                (h-b.min_coord(2))*(h-b.min_coord(2)) : (b.max_coord(2)-h)*(b.max_coord(2)-h);
    return d0 + d1 + d2;
  }

  Double new_distance(Double& dist, Double old_off, Double new_off,
                      int /* cutting_dimension */)  const {
    return dist + new_off*new_off - old_off*old_off;
  }

  Double transformed_distance(Double d) const { return d*d; }

  Double inverse_of_transformed_distance(Double d) { return sqrt(d); }

}; // end of struct Distance


typedef CGAL::Search_traits<Double, CGAL::Point, const Double*, Construct_coord_iterator> Traits;
typedef CGAL::Orthogonal_k_neighbor_search<Traits, Distance> K_neighbor_search;
typedef K_neighbor_search::Tree Tree;

#endif // USE_CGAL

#ifndef USE_ADOLC
#ifdef USE_FLANN
#include <flann/flann.hpp>
#endif // USE_FLANN
#endif


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
      CGAL, FLANN
    };
    
    //! Constructor
    CoefFunctionScatteredData(PtrParamNode& scatteredDataNode);
    //! Destructor
    virtual ~CoefFunctionScatteredData(){;}
    
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
    
    std::vector< std::vector<Double> > coordinates_;
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

#ifdef USE_FLANN
#ifndef USE_ADOLC
    boost::shared_ptr< flann::Index<flann::L2<Double> > > index_;
    boost::shared_ptr< flann::Matrix<Double> > dataset_;

    void KNNSearch_FLANN(const Vector<Double> globPoint,
                         StdVector< Vector<Double> >& neighbors,
                         StdVector< Double >& l2Distances,
                         StdVector< Vector<T> >& vectors);
#endif
#endif
  };
}

#endif
