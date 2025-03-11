/*
 * CoefFunctionDummy.cc
 *
 *  Created on: 20 Jun 2024
 *      Author: Dominik Mayrhofer
 */

#include "CoefFunctionDummy.hh"

#include "DataInOut/Logging/LogConfigurator.hh"

DEFINE_LOG(coeffctdummy, "coeffctdummy")

namespace CoupledField {
  CoefFunctionDummy::CoefFunctionDummy(SolutionType type,
                                        shared_ptr<EntityList> list,
                                        std::string pdeName) : CoefFunction() {

    // not sure about that one
    dependType_ = GENERAL;
    isAnalytic_ = false;
    isComplex_  = false;
    dimType_ = NO_DIM; // this has to be set later on

    type_ = type;
    list_ = list;
    pdeName_ = pdeName;

    //coef_ = feFct->GetPDE()->GetCoefFct(feFct->GetResultInfo()->resultType);
  }

  CoefFunctionDummy::~CoefFunctionDummy() {
    ;
  }

  void CoefFunctionDummy::GetScalar( Double& scal, const LocPointMapped& surfLpm ) {
    LOG_DBG(coeffctdummy) << "+++++ CoefFunctionDummy::GetScalar ++++++";
    coef_->GetScalar( scal, surfLpm );
  }

  void CoefFunctionDummy::GetVector( Vector<Double>& vec, const LocPointMapped& surfLpm ) {
    LOG_DBG(coeffctdummy) << "+++++ CoefFunctionDummy::GetVector ++++++";
    coef_->GetVector( vec, surfLpm );
  }

  void CoefFunctionDummy::GetScalarValuesAtCoords( const StdVector<Vector<Double> >& globCoord,
                                              StdVector< Double >& values,
                                              Grid* ptGrid,
                                              const StdVector<shared_ptr<EntityList> >& srcEntities) {
    LOG_DBG(coeffctdummy) << "+++++ CoefFunctionDummy::GetScalarValuesAtCoords ++++++";
    coef_->GetScalarValuesAtCoords( globCoord, values, ptGrid, srcEntities );
  }

  void CoefFunctionDummy::GetScalarValuesAtCoords( const StdVector<Vector<Double> >& globCoord,
                                              StdVector< Complex >& values,
                                              Grid* ptGrid,
                                              const StdVector<shared_ptr<EntityList> >& srcEntities) {
    LOG_DBG(coeffctdummy) << "+++++ CoefFunctionDummy::GetScalarValuesAtCoords ++++++";
    EXCEPTION("CoefFunctionDummy::GetScalarValuesAtCoords not implemented for complex type variables!");
  }
  
  void CoefFunctionDummy::GetVectorValuesAtCoords( const StdVector<Vector<Double> >& globCoord,
                                                   StdVector< Vector<Double> >& values, 
                                                   Grid* ptGrid,
                                                   const StdVector<shared_ptr<EntityList> >& srcEntities,
                                                   bool updatedGeo) {
    LOG_DBG(coeffctdummy) << "+++++ CoefFunctionDummy::GetVectorValuesAtCoords ++++++";
    coef_->GetVectorValuesAtCoords( globCoord, values, ptGrid, srcEntities, updatedGeo );
  }
}
