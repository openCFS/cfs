#include "Domain/CoefFunction/CoefFunctionOpt.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/Logging/log.hpp"

#ifdef _OPENMP
  #include <omp.h>
#endif


namespace CoupledField
{

DECLARE_LOG(coef)
DEFINE_LOG(coef, "coef")

CoefFunctionOpt::CoefFunctionOpt(DesignSpace* design, PtrCoefFct orgMat, SinglePDE* pde) : CoefFunction()
{
  // this type of coefficient is always constant
  dependType_ = GENERAL;
  isAnalytic_ = false;
  isComplex_ = false;
  supportDerivative_ = false;
  dimType_ = orgMat->GetDimType();

  this->design  = design;
  this->orgMat  = orgMat;
  this->form    = NULL; // can only be set later
  this->state = OPT;
  this-> subTensor = pde->GetSubTensorType();
}

void CoefFunctionOpt::SetToOptimization()
{
  state = OPT;
  shadowMat.reset(); // equivalent to = nullptr (C++11)
  direction = DesignElement::NO_DERIVATIVE;
}

void CoefFunctionOpt::SetToOrgMaterial()
{
  state = ORG;
  shadowMat.reset(); // equivalent to = nullptr (C++11)
  direction = DesignElement::NO_DERIVATIVE;
}

void CoefFunctionOpt::SetToShadow(PtrCoefFct shadow)
{
  state = SHADOW;
  shadowMat = shadow;
  direction = DesignElement::NO_DERIVATIVE;
}

void CoefFunctionOpt::SetToMaterialDerivative(DesignElement::Type dir)
{
  state = DIRECTION;
  direction = dir;
  assert(dir != DesignElement::NO_DERIVATIVE && dir != DesignElement::NO_MULTIMATERIAL && dir != DesignElement::NO_TYPE);
}



template <class T>
void CoefFunctionOpt::GetTensor(Matrix<T>& coefMat, const LocPointMapped& lpm)
{
 assert(this->dimType_ == TENSOR);
 // if no coordinate system is set, just
  // use internal vector
 if(coordSys_ != NULL)
   EXCEPTION("the rotation is not fully finished ':-(\n");

 switch(state)
 {
 case DIRECTION:
 case OPT:
   // DesignSpace::ApplyPhysicalDesign(Tensor) is not thread save yet :(
   #pragma omp critical (SAVE_OPT_TENSOR)
   {
   // the element does not necessarily lay in the design space!
   // if ApplyPhysicalDesign() returns true, coefMat is already set
   if(!design->ApplyPhysicalDesign<T>(shared_from_this(), coefMat, &lpm))
     orgMat->GetTensor(coefMat, lpm);
   }
   break;
 case ORG:
   orgMat->GetTensor(coefMat, lpm);
   break;
 case SHADOW:
   shadowMat->GetTensor(coefMat, lpm);
   break;
 }

  LOG_DBG3(coef) << "CFO:GT el=" << lpm.ptEl->elemNum  << " state=" << state << " shadow=" << (shadowMat ? "set" : "not set") << " -> " << coefMat.ToString(0, false);
}

template <class T>
void CoefFunctionOpt::GetScalar(T& scal, const LocPointMapped& lpm)
{
  assert(this->dimType_ == SCALAR);
  // if no coordinate system is set, just
   // use internal vector
  if(coordSys_ != NULL)
    EXCEPTION("the rotation is not fully finished ':-(\n");

  switch(state)
  {
  case DIRECTION:
  case OPT:
    // the element does not necessarily lay in the design space!
    // if ApplyPhysicalDesign() returns true, coefMat is already set
    if(!design->ApplyPhysicalDesign<T>(shared_from_this(), scal, &lpm))
      orgMat->GetScalar(scal, lpm);
    break;
  case ORG:
    orgMat->GetScalar(scal, lpm);
    break;
  case SHADOW:
    shadowMat->GetScalar(scal, lpm);
    break;
  }

  LOG_DBG3(coef) << "CFO:GT el=" << lpm.ptEl->elemNum  << " state=" << state << " shadow=" << (shadowMat ? "set" : "not set") << " -> " << scal;

}

}
