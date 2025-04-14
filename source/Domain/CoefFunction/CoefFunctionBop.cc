#include "Domain/CoefFunction/CoefFunctionBop.hh"

namespace CoupledField
{

CoefFunctionBop::CoefFunctionBop() : CoefFunction()
{
  isAnalytic_ = false;
  isComplex_ = false;
  dependType_ = GENERAL;
  this->form_ = NULL; // can only be set later
}


void CoefFunctionBop::CalcElementMatrixLpm(Matrix<Double>& elemMat, const LocPointMapped& lpm, BaseFE* ptFe){
  Matrix<Double> curElemMat;
  this->form_->CalcElementMatrixLpm(curElemMat, ptFe, lpm);
  elemMat = curElemMat;
}

}
