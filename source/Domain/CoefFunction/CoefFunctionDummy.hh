/*
 * CoefFunctionDummy.hh
 *
 *  Created on: 20 Jun 2024
 *      Author: Dominik Mayrhofer
 */

#ifndef SOURCE_DOMAIN_COEFFUNCTION_COEFFUNCTIONDUMMY_HH_
#define SOURCE_DOMAIN_COEFFUNCTION_COEFFUNCTIONDUMMY_HH_

#include "CoefFunction.hh"
#include "FeBasis/FeFunctions.hh"
#include "PDE/SinglePDE.hh"


namespace CoupledField {

// ============================================================================
//  Coef Function Dummy
// ============================================================================
//!
//! Dummy coefFunction that can be set to have a coefFunction set, where the
//! actual coefFunction used for the evaluation will be set at a later point
//! 
//! \note This class only works for real-valued scalar data.

  class CoefFunctionDummy : public CoefFunction {

  public:

    //! Constructor
    CoefFunctionDummy(SolutionType type, shared_ptr<EntityList> list, std::string pdeName);

    //! Destructor
    virtual ~CoefFunctionDummy();

    virtual string GetName() const { return "CoefFunctionDummy"; }

    //! Return real-valued scalar at integration point
    void GetScalar(Double& scal, const LocPointMapped& surflpm );

    //! Return real-valued vector at integration point
    void GetVector(Vector<Double>& coefVec, const LocPointMapped& lpm );

    //! Return scalar values at global coordinate locations
    void GetScalarValuesAtCoords( const StdVector<Vector<Double> >& globCoord,
                                  StdVector< Double >& values, 
                                  Grid* ptGrid,
                                  const StdVector<shared_ptr<EntityList> >& srcEntities =
                                  StdVector<shared_ptr<EntityList> >() );

    //! Return scalar values at global coordinate locations
    void GetScalarValuesAtCoords( const StdVector<Vector<Double> >& globCoord,
                                  StdVector< Complex >& values, 
                                  Grid* ptGrid,
                                  const StdVector<shared_ptr<EntityList> >& srcEntities =
                                  StdVector<shared_ptr<EntityList> >() );

    //! Return vectorial values at global coordinate locations
    void GetVectorValuesAtCoords( const StdVector<Vector<Double> >& globCoord,
                                  StdVector< Vector<Double> >& values,
                                  Grid* ptGrid,
                                  const StdVector<shared_ptr<EntityList> >& srcEntities =
                                  StdVector<shared_ptr<EntityList> >(),
                                  bool updatedGeo = false);

    //! Set the coefFunction used to evaluate the dummy coefFunction
    void SetCoef(PtrCoefFct coef) {
      coef_ = coef;
      dimType_ = coef->GetDimType();
    }

    //! Get all information necessary to get the coefFunction afterwards
    void GetCoefInfo(SolutionType &type, shared_ptr<EntityList> &list, std::string &pdeName) {
      type = type_;
      list = list_;
      pdeName = pdeName_;
    }

  protected:
    // actual coefFunction that gets passed through
    PtrCoefFct coef_;

    SolutionType type_;
    shared_ptr<EntityList> list_;
    std::string pdeName_;
  };

}


#endif /* SOURCE_DOMAIN_COEFFUNCTION_COEFFUNCTIONDUMMY_HH_ */
