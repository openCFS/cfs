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
  this->enabled = true;
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


 // if ApplyPhysicalDesign() returns true, coefMat is already set
 if(!enabled || !design->ApplyPhysicalDesign(shared_from_this(), coefMat, &lpm))
   orgMat->GetTensor(coefMat, lpm);

  LOG_DBG3(coef) << "CFO:GT el=" << lpm.ptEl->elemNum  << " en=" << enabled << " -> " << coefMat.ToString(0, false);
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


  // if ApplyPhysicalDesign() returns true, coefMat is already set
  if(!enabled || !design->ApplyPhysicalDesign(shared_from_this(), scal, &lpm))
    orgMat->GetScalar(scal, lpm);

   LOG_DBG3(coef) << "CFO:GT el=" << lpm.ptEl->elemNum  << " en=" << enabled << " -> " << scal;

}

}
