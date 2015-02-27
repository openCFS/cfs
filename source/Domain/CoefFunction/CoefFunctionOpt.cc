#include "Domain/CoefFunction/CoefFunctionOpt.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/Logging/log.hpp"

namespace CoupledField
{

DECLARE_LOG(coef)
DEFINE_LOG(coef, "coef")

CoefFunctionOpt::CoefFunctionOpt(DesignSpace* design, PtrCoefFct orgMat) : CoefFunction()
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
}

void CoefFunctionOpt::SetToOptimization()
{
  state = OPT;
  shadowMat.reset(); // equivalent to = nullptr (C++11)
}

void CoefFunctionOpt::SetToOrgMaterial()
{
  state = ORG;
  shadowMat.reset(); // equivalent to = nullptr (C++11)
}

void CoefFunctionOpt::SetToShadow(PtrCoefFct shadow)
{
  state = SHADOW;
  shadowMat = shadow;
}


void CoefFunctionOpt::GetTensor(Matrix<double>& coefMat, const LocPointMapped& lpm)
{
 assert(this->dimType_ == TENSOR);
 // if no coordinate system is set, just
  // use internal vector
 if(coordSys_ != NULL)
   EXCEPTION("the rotation is not fully finished ':-(\n");

 PtrCoefFct f = shared_from_this();
 if(f->IsComplex())
 {
   assert(false); // what does this case/ code mean??
   return;
 }

 switch(state)
 {
 case OPT:
   // the element does not necessarily lay in the design space!
   // if ApplyPhysicalDesign() returns true, coefMat is already set
   if(!design->ApplyPhysicalDesign(shared_from_this(), coefMat, &lpm))
     orgMat->GetTensor(coefMat, lpm);
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


void CoefFunctionOpt::GetScalar(double& scal, const LocPointMapped& lpm)
{
  assert(this->dimType_ == SCALAR);
  // if no coordinate system is set, just
   // use internal vector
  if(coordSys_ != NULL)
    EXCEPTION("the rotation is not fully finished ':-(\n");

  PtrCoefFct f = shared_from_this();
  if(f->IsComplex())
  {
    assert(false); // what does this case/ code mean??
    return;
  }


  switch(state)
  {
  case OPT:
    // the element does not necessarily lay in the design space!
    // if ApplyPhysicalDesign() returns true, coefMat is already set
    if(!design->ApplyPhysicalDesign(shared_from_this(), scal, &lpm))
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
