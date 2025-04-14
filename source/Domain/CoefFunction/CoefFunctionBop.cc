#include "Domain/CoefFunction/CoefFunctionBop.hh"

namespace CoupledField
{

CoefFunctionBop::CoefFunctionBop(shared_ptr<BaseFeFunction> feFct) : CoefFunction()
{
  isAnalytic_ = false;
  isComplex_ = false;
  dependType_ = GENERAL;
  this->form_ = NULL; // can only be set later
  this->feFct_ = feFct;
}


template <class T>
void CoefFunctionBop::CalcElementMatrixLpm(Matrix<T>& elemMat, const LocPointMapped& lpm){
  Matrix<T> curElemMat;
  this->form_->CalcElementMatrixLpm(curElemMat, this->feFct_->GetFeSpace()->GetFe( lpm.ptEl->elemNum ), lpm);
  elemMat = curElemMat;
}


template void CoefFunctionBop::CalcElementMatrixLpm<double>(Matrix<double>&, const LocPointMapped&);
template void CoefFunctionBop::CalcElementMatrixLpm<Complex>(Matrix<Complex>&, const LocPointMapped&);
}
