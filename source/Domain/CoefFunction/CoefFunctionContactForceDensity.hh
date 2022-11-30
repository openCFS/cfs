/*
 * CoefFunctionContactForceDensity.hh
 *
 *  Created on: 18 Nov 2022
 *      Author: Dominik Mayrhofer
 */

#ifndef SOURCE_DOMAIN_COEFFUNCTION_COEFFUNCTIONCONTACTFORCEDENSITY_HH_
#define SOURCE_DOMAIN_COEFFUNCTION_COEFFUNCTIONCONTACTFORCEDENSITY_HH_

#include "CoefFunction.hh"
#include "FeBasis/FeFunctions.hh"

namespace CoupledField {

// ============================================================================
//  Coef Function Contact Force Density
// ============================================================================
//! Coefficient function for calculating distance based force densities via
//! analytic contact laws
//!
//! The coefFct gets passed a vector of surface pairs, the adjacent volumes, the 
//! corresponding feFunctions, the contact laws as well as the bools stating if 
//! nodal or surface midpoint evaluation should be used. It has to be noted that 
//! each volume is allowed to have only one unique surface pair, otherwise, an 
//! exception is thrown.
//! For each volume and surface pair, a coefFunction for the distance evaluation
//! is created. Additionally, the contact law can be given analytically and is
//! evaluated based on the distance of the two surfaces. The GetVector function
//! returns the force density vector based on this contact law accordingly.
//! Since we perform this calculation for each call and use the current location 
//! of the nodes, this function is applicable for moving mesh applications.
//! 
//! \note This class only works for real-valued data.

  class CoefFunctionContactForceDensity : public CoefFunction {

  public:

    //! Constructor
    CoefFunctionContactForceDensity(shared_ptr<BaseFeFunction> feFct,
                                StdVector<std::string> surfList1,
                                StdVector<std::string> surfList2,
                                StdVector<std::string> volumeList,
                                StdVector<std::string> contactLawList,
                                StdVector<bool> useSurfaceMidpointsList);

    //! Destructor
    virtual ~CoefFunctionContactForceDensity();

    virtual string GetName() const { return "CoefFunctionContactForceDensity"; }

    //! Return real-valued vector at integration point
    void GetVector(Vector<Double>& vec, const LocPointMapped& lpm );

  protected:
    //! Source and target region
    RegionIdType srcRegion_; // here we comput the force
    RegionIdType targetRegion_; // used for the distance evaluation

    //! Set feFunction and xml input data
    shared_ptr<BaseFeFunction> feFct_;
    StdVector<std::string> surfList1_;
    StdVector<std::string> surfList2_;
    StdVector<std::string> volumeList_;
    StdVector<std::string> contactLawList_;
    StdVector<bool> useSurfaceMidpointsList_;
    StdVector<PtrCoefFct> coefs_;

    //! Set with all regionIdTypes
    std::set<RegionIdType> volRegions_;
    std::set<RegionIdType> surfRegions1_;
    std::set<RegionIdType> surfRegions2_;

    //! Pointer to MathParser object
    MathParser * mParser_ = nullptr;

    Vector<Double> rVals_;

    //! vector with related handles for math parser
    StdVector<unsigned int> rHandles_;
  };

}


#endif /* SOURCE_DOMAIN_COEFFUNCTION_COEFFUNCTIONCONTACTFORCEDENSITY_HH_ */
