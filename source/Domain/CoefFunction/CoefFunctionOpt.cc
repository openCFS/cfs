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

void CoefFunctionOpt::GetTensor(Matrix<double>& coefMat, const LocPointMapped& lpm )
{
 assert(this->dimType_ == TENSOR);
 // if no coordinate system is set, just
  // use internal vector
 if(coordSys_ != NULL)
   EXCEPTION("the rotation is not fully finished ':-(\n");

 PtrCoefFct f = shared_from_this();
 if(f->IsComplex())
   return;


 if(!enabled || !design->TryApplyPhysicalDesign(shared_from_this(), coefMat, &lpm))
   orgMat->GetTensor(coefMat, lpm);

  LOG_DBG3(coef) << "CFO:GT el=" << lpm.ptEl->elemNum  << " en=" << enabled << " -> " << coefMat.ToString(0, false);
}

/*
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template class CoefFunctionOpt<double>;
  template class CoefFunctionOpt<Complex>;
#endif
*/

}
