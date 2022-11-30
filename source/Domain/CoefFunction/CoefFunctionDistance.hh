/*
 * CoefFunctionDistance.hh
 *
 *  Created on: 18 Nov 2022
 *      Author: Dominik Mayrhofer
 */

#ifndef SOURCE_DOMAIN_COEFFUNCTION_COEFFUNCTIONDISTANCE_HH_
#define SOURCE_DOMAIN_COEFFUNCTION_COEFFUNCTIONDISTANCE_HH_

#include "CoefFunction.hh"
#include "FeBasis/FeFunctions.hh"


namespace CoupledField {

// ============================================================================
//  Coef Function Distance
// ============================================================================
//! Coefficient function for calculating distances between surface regions
//!
//! The coefFct gets passed a surface pair, the adjacent volume, the corresponding 
//! feFunction as well as a bool if nodal or surface midpoint evaluation should be 
//! used. Depending on the mode, we add these points to our target points.
//! By calling GetScalar we re-calculate the target points and evaluate the
//! distance to a certain query point. Since we perform this calculation for each
//! call and use the current location of the nodes, this function is applicable
//! for moving mesh applications.
//! 
//! \note This class only works for real-valued scalar data.

  class CoefFunctionDistance : public CoefFunction {

  public:

    //! Constructor
    CoefFunctionDistance(shared_ptr<BaseFeFunction> feFct,
                                std::string surf1,
                                std::string surf2,
                                std::string vol,
                                bool useSurfaceMidpoints);

    //! Destructor
    virtual ~CoefFunctionDistance();

    virtual string GetName() const { return "CoefFunctionDistance"; }

    //! Return real-valued vector at integration point
    void GetScalar(Double& scal, const LocPointMapped& surflpm );

    //! Add target point of the target surface
    void AddPoint(Vector<Double> targetPoint);

    //! Query point and return distance
    void QueryPoint(Double& dist, Vector<Double> queryPoint);

    //! Reset the point list
    void ResetPoints();


  protected:
    //! Source and target region
    RegionIdType srcRegion_; // here we comput the force
    RegionIdType targetRegion_; // used for the distance evaluation

    //! Set feFunction and xml input data
    shared_ptr<BaseFeFunction> feFct_;
    std::string surf1_;
    std::string surf2_;
    std::string vol_;
    bool useSurfaceMidpoints_;

    //! Storage for all regionIdTypes
    RegionIdType volRegion_;
    RegionIdType surfRegion1_;
    RegionIdType surfRegion2_;

    //! Vector of target points
    StdVector<Vector<Double>> targetPoints_;
  };

}


#endif /* SOURCE_DOMAIN_COEFFUNCTION_COEFFUNCTIONDISTANCE_HH_ */
