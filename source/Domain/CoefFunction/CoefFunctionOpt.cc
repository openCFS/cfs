#include "Domain/CoefFunction/CoefFunctionOpt.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

#ifdef _OPENMP
  #include <omp.h>
#endif


namespace CoupledField
{
DEFINE_LOG(coef, "coef")

CoefFunctionOpt::CoefFunctionOpt(DesignSpace* design, PtrCoefFct orgMat, MaterialType mt, SinglePDE* pde) : CoefFunction()
{
  isAnalytic_ = false;
  isComplex_ = false;
  supportDerivative_ = false;
  dimType_ = orgMat->GetDimType();
  dependType_ = SPACE;
  this->direction = DesignElement::NO_DERIVATIVE;
  this->design  = design;
  this->orgMat  = orgMat;
  this->materialType = mt;
  this->form    = NULL; // can only be set later
  this->formL   = NULL;
  this->state = OPT;
  this->pde = pde;
  this->subTensor = pde->GetSubTensorType();
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
  Matrix<T> locMatrix;

  switch(state)
  {
  case DIRECTION:
  case OPT:
    // the element does not necessarily lay in the design space!
    // if ApplyPhysicalDesign() returns true, coefMat is already set
    if(!design->ApplyPhysicalDesign<T>(this, locMatrix, &lpm))
      orgMat->GetTensor(locMatrix, lpm);
    //if (coefMat.GetNumCols() > 0) {
     //assert(design->TestTensorPosDef<T>(coefMat, &lpm , shared_from_this()->GetMaterialDerivative()));
    //}
    break;
  case ORG:
    orgMat->GetTensor(locMatrix, lpm);
    break;
  case SHADOW:
    shadowMat->GetTensor(locMatrix, lpm);
    break;
  }

  TransformTensorByCoordSys(coefMat, locMatrix, lpm);

  LOG_DBG3(coef) << "CFO:GT el=" << lpm.ptEl->elemNum  << " state=" << state << " shadow=" << (shadowMat ? "set" : "not set")
                 << " lM=" << locMatrix.ToString() << " -> " << coefMat.ToString();
}

template <class T>
void CoefFunctionOpt::GetScalar(T& scal, const LocPointMapped& lpm)
{
  assert(this->dimType_ == SCALAR);

  switch(state)
  {
  case DIRECTION:
  case OPT:
  {
    // the element does not necessarily lay in the design space!
    // if ApplyPhysicalDesign() returns true, coefMat is already set
    bool designset = design->ApplyPhysicalDesign<T>(this, scal, &lpm);
    if(!designset)
      orgMat->GetScalar(scal, lpm);
    break;
  }
  case ORG:
    orgMat->GetScalar(scal, lpm);
    break;
  case SHADOW:
    shadowMat->GetScalar(scal, lpm);
    break;
  }

  LOG_DBG3(coef) << "CFO:GS el=" << lpm.ptEl->elemNum  << " state=" << state << " shadow=" << (shadowMat ? "set" : "not set") << " -> " << scal;

}

template <class T>
void CoefFunctionOpt::GetVector(Vector<T>& vec, const LocPointMapped& lpm)
{
  assert(this->dimType_ == VECTOR);
  Vector<T> locVec;

  switch(state)
  {
  case DIRECTION:
  case OPT:
    // the element does not necessarily lay in the design space!
    // if ApplyPhysicalDesign() returns true, coefMat is already set
    if(!design->ApplyPhysicalDesign<T>(this, locVec, &lpm))
      orgMat->GetVector(locVec, lpm);
    break;
  case ORG:
    orgMat->GetVector(locVec, lpm);
    break;
  case SHADOW:
    shadowMat->GetVector(locVec, lpm);
    break;
  }

  TransformVectorByCoordSys(vec, locVec, lpm);

  LOG_DBG3(coef) << "CFO:GV node=" << lpm.lp.number  << " state=" << state << " shadow=" << (shadowMat ? "set" : "not set") << " -> " << vec;

}

template void CoefFunctionOpt::GetTensor<double>(Matrix<double>&, const LocPointMapped&);
template void CoefFunctionOpt::GetTensor<Complex>(Matrix<Complex>&, const LocPointMapped&);
template void CoefFunctionOpt::GetScalar<double>(Double&, const LocPointMapped&);
template void CoefFunctionOpt::GetScalar<Complex>(Complex&, const LocPointMapped&);
template void CoefFunctionOpt::GetVector<double>(Vector<Double>&, const LocPointMapped&);

}
