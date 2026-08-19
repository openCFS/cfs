/*
 * CoefFunctionSUPGDeriv.cc
  *  Created on: 20.08.2026
 *      Author: Lorenz Klimon
 */

#include "CoefFunctionSUPGDeriv.hh"
#include "CoefFunctionSUPG.hh"
#include <cmath>

#include "FeBasis/FeSpace.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

namespace CoupledField {
DEFINE_LOG(coeffunctionSUPGDeriv, "coeffunctionSUPGDeriv")

CoefFunctionSUPGDeriv::CoefFunctionSUPGDeriv( PtrCoefFct velocityField,  PtrCoefFct materialCoeff, shared_ptr<BaseFeFunction> FeFnc)
  : CoefFunction(),
    velField_(velocityField),
    matCoeff_(materialCoeff), 
    feFct_(FeFnc)       
{
  dimType_ = SCALAR;
  isAnalytic_ = false;
  isComplex_ = false;
  dependType_ = SPACE;
}

//! \see CoefFunction::GetScalar
  void CoefFunctionSUPGDeriv::GetScalar(Double& scal, const LocPointMapped& lpm) {
    LOG_DBG(coeffunctionSUPGDeriv) << "+++++ coeffunctionSUPGDeriv::GetScalar ++++++"; 


    //velocity 
    Vector<double> v;
    // material parameter
    double m = 0.0;

    velField_->GetVector(v,lpm);
    double velNorm = v.NormL2();
    int nDim = velField_->GetVecSize(); // Dimentions of space 2D or 3D

    // if velocity is 0, the derivative of the stabilization factor is also zero
    double epsvelNorm = 1e-13; 
    if(abs(velNorm) <= epsvelNorm){
      LOG_DBG(coeffunctionSUPGDeriv) << "velNorm=0, returning";
      scal = 0.0;
      return;
    }

    // length of the element
    double lElem = CoefFunctionSUPG::CalcElementLength(feFct_, lpm, v, nDim);    

    // To make it work robustly both for TENSOR and SCALAR material parameters
    switch (matCoeff_->GetDimType())
    {
      case TENSOR:
      {
        EXCEPTION("The material property is TENSOR. It is not implemented for the Derivative of SUPG");
        break;
      }
      case SCALAR:
      {
        double matrixM;
        matCoeff_->GetScalar(matrixM,lpm);
        m = matrixM;
        break;
      }
      case VECTOR:
      {
        EXCEPTION("The material property is VECTOR. It is not implemented for the Derivative of SUPG");
        break;
      }
      default:
      {
        EXCEPTION("The material property is an unknown type. It is not implemented for the Derivative of SUPG");
        break;
      }
    }
    
    //compute Peclet-Number
    double peclet = (velNorm * lElem) / (2 * m);

    //compute the stabilization factor
    scal = lElem / (2 * velNorm * m) * (peclet / std::pow(std::sinh(peclet), 2) - (1/peclet));

    // For debugging purposes
    LOG_DBG(coeffunctionSUPGDeriv) << "Calculated element size " << lElem << std::endl;
    LOG_DBG(coeffunctionSUPGDeriv) << "Calculated Peclet Number " << peclet << std::endl;
    LOG_DBG(coeffunctionSUPGDeriv) << "Calculated derivative of stabilization factor " << scal << std::endl;

  }
}