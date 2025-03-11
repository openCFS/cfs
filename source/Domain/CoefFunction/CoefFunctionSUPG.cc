/*
 * CoefFunctionSUPG.cc
 *
 *  Created on: 25.04.2024
 *      Author: Vladimir Filkin and Lukas Schafferhofer
 */

#include "CoefFunctionSUPG.hh"

#include "FeBasis/H1/H1Elems.hh"
#include "FeBasis/FeFunctions.hh"
#include "Forms/Operators/BaseBOperator.hh"
#include "FeBasis/FeSpace.hh"
#include "DataInOut/Logging/LogConfigurator.hh"



namespace CoupledField {
DEFINE_LOG(coeffunctionSUPG, "coeffunctionSUPG")

    static double coth(double x){;
    return 1/std::tanh(x);
    }   

  CoefFunctionSUPG::CoefFunctionSUPG( PtrCoefFct velocityField,  PtrCoefFct materialCoeff,
   shared_ptr<BaseFeFunction> FeFnc)
   : CoefFunction(),
    feFct_(FeFnc),
    matCoeff_(materialCoeff),
    velField_(velocityField)                
  {
    dimType_ = SCALAR; // should be determined from input
    isAnalytic_ = false;
    isComplex_  = false;
    dependType_ = SPACE;
  }

  //! \see CoefFunction::GetScalar
  void CoefFunctionSUPG::GetScalar(Double& scal, const LocPointMapped& lpm) {
    LOG_DBG(coeffunctionSUPG) << "+++++ coeffunctionSUPG::GetScalar ++++++"; 


    //velocity (for HeatPDE |u|/k)
    Vector<double> v;
    // material parameter
    double m = 0.0;
    // length of the element
    double lElem;

    velField_->GetVector(v,lpm);
    double velNorm = v.NormL2();
    int nDim = velField_->GetVecSize(); // Dimentions of space 2D or 3D

    // Calculate the dimention of an element in velocity direction
    // lElem = 2 / Sum((v / velNorm) \dot (\nabla shape_function) )
    // The idea is that we use projection of the unit length velocity field from global configuration into reference configuration.
    // So the length of the global unit vector in reference configuration is a reversed length 
    // in the global configuration of the dimention of the element in the velocity direction. 
    BaseFE* ptFe = feFct_->GetFeSpace()->GetFe(lpm.ptEl->elemNum);
    FeH1 *fe = ( static_cast<FeH1*>(ptFe) );
    Matrix<Double> xiDx;
    fe->GetGlobDerivShFnc(xiDx, lpm, lpm.shapeMap->GetElem(),1);
    double denom = 0;
    for (UInt i=0; i < xiDx.GetNumRows(); ++i){
      double scalar_product = 0;
      for (int j=0; j< nDim; ++j){
        std::cout << std::fixed;
        scalar_product += v[j]*xiDx[i][j];
      }
      denom += abs(scalar_product)/velNorm;
    }
    lElem = 2/(denom);

    // To make it work robustly both for TENSOR and SCALAR material parameters
    switch (matCoeff_->GetDimType())
    {
      case TENSOR:
      {
        Matrix<Double> matrixM;
        matCoeff_->GetTensor(matrixM,lpm);

        //material properties only in the direction of the velocity m = v * matrixM * v / |v|^2

        Vector<Double> resVec = Vector<Double>(nDim);
        matrixM.Mult(v, resVec);
        m = (v*resVec)/(velNorm*velNorm);
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
        EXCEPTION("The material property is VECTOR. It is not implemented for SUPG and Artificial diffusion");
        break;
      }
      default:
      {
        EXCEPTION("The material property is an unknown type. It is not implemented for SUPG and Artificial diffusion");
        break;
      }
    }
    
    //compute Peclet-Number
    double peclet = (velNorm * lElem) / (2 * m);
    
    // check if peclet == 0.0 with the tolerance epsilon
    double epsilon = 1e-13; 
    if (abs(peclet) <= epsilon)
      scal = 0.0;
    //compute the stabilization factor
    else 
      scal = lElem * (coth(peclet) - (1/peclet)) / (2 * velNorm);

    // For debugging purposes
    LOG_DBG(coeffunctionSUPG) << "Calculated element size " << lElem << std::endl;
    LOG_DBG(coeffunctionSUPG) << "Calculated Peclet Number " << peclet << std::endl;
    LOG_DBG(coeffunctionSUPG) << "Calculated stabilization factor " << scal << std::endl;

  }

  void CoefFunctionSUPG::GetVector(Vector<Double>& vec, const LocPointMapped& lpm ) {
  //Vector<Complex> DependentVec;
  //depCoef_->GetVector(DependentVec, lpm);

  Vector<Double> RealDependentVec;
  RealDependentVec.Resize(3);
  RealDependentVec.Init(0);

    for(UInt i = 0 ; i < 3; ++i){
      vec[i] = std::real(RealDependentVec[i]); 
    }
  }
}
