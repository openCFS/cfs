/*
 * CoefFunctionGeom.cc
 *
 *  Created on: 05 Jan 2026
 *      Author: Dominik Mayrhofer
 */

#include "CoefFunctionGeom.hh"

#include "DataInOut/Logging/LogConfigurator.hh"

DEFINE_LOG(coeffctgeom, "coeffctgeom")

namespace CoupledField {
  CoefFunctionGeom::CoefFunctionGeom(std::string resultType) : CoefFunction() {
    dimType_ = SCALAR; // should be determined from input
    dependType_ = GENERAL; // should be determined from input
    isAnalytic_ = false;
    isComplex_  = false;

    // Store volume list used for evaluation
    resultType_ = resultType;
  }

  CoefFunctionGeom::~CoefFunctionGeom() {
    ;
  }

  void CoefFunctionGeom::GetScalar(Double& scal, const LocPointMapped& lpm ) {
    LOG_DBG(coeffctgeom) << "+++++ CoefFunctionGeom::GetScalar ++++++";
    
    if( resultType_=="Jacobian" ){
      LOG_DBG(coeffctgeom) << "Calculating Jacobian";
      shared_ptr<ElemShapeMap> shapeMap = lpm.GetShapeMap();
      // not needed, but we have to define the return value
      Matrix<Double> mat;
      scal = shapeMap->CalcJDet(mat, lpm.lp, true);
    }
  }
}
