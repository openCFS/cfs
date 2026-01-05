/*
 * CoefFunctionGeom.hh
 *
 *  Created on: 05 Jan 2026
 *      Author: Dominik Mayrhofer
 */

#ifndef SOURCE_DOMAIN_COEFFUNCTION_COEFFUNCTIONGEOM_HH_
#define SOURCE_DOMAIN_COEFFUNCTION_COEFFUNCTIONGEOM_HH_

#include "CoefFunction.hh"


namespace CoupledField {

// ============================================================================
//  Coef Function Geom
// ============================================================================
//! Coefficient function for calculating geometry related results of regions
//!
//! This coefFct gives access to fundamental mesh features as for example the
//! Jacobian and provides an interface to pass it as results during computations
//! as used in for example the smoothPDE.

  class CoefFunctionGeom : public CoefFunction {

  public:

    //! Constructor
    CoefFunctionGeom(std::string result);

    //! Destructor
    virtual ~CoefFunctionGeom();

    virtual string GetName() const { return "CoefFunctionGeom"; }

    //! Return real-valued vector at integration point
    void GetScalar(Double& scal, const LocPointMapped& surflpm );

  protected:
    //! Result to evaluate
    std::string resultType_;
  };

}


#endif /* SOURCE_DOMAIN_COEFFUNCTION_COEFFUNCTIONGEOM_HH_ */
