#include "CoefFunctionDiagTensorFromScalar.hh"

#include <limits>

namespace CoupledField {

  CoefFunctionDiagTensorFromScalar::
  CoefFunctionDiagTensorFromScalar(const StdVector<PtrCoefFct>& scalVals)
 : CoefFunction() {

   // this type of coefficient is never analytic
   isAnalytic_ = false;

   dimType_ = TENSOR;

   // Loop over all entries and check if all entries have the same
   assert(scalVals.GetSize() > 0 );
   size_ = scalVals.GetSize();
   dependType_ = scalVals[0]->GetDependency();
   isComplex_ = scalVals[0]->IsComplex();
   for( UInt i = 1; i < size_; ++i ) {
     if( dependType_ != scalVals[i]->GetDependency() ) {
       EXCEPTION("All diagonal functions must have the same dependency type");
     }
     if( isComplex_ != scalVals[i]->IsComplex() ) {
       EXCEPTION("All diagonal functions must have the same complex type");
     }
   }
   scalVals_ = scalVals;
 }

//! \copydoc CoefFunction::IsZero
bool CoefFunctionDiagTensorFromScalar::IsZero() const {
  bool ret = true;
  for( UInt i = 0; i < size_; ++i ) {
    ret |= scalVals_[i]->IsZero();
  }
  return ret;
}

//! \copydoc CoefFunction::ToString
std::string CoefFunctionDiagTensorFromScalar::ToString() const {
  std::string ret;
  ret = "DiagTensor(";
  for( UInt i = 0; i < size_-1; ++i ) {
    ret += scalVals_[i]->ToString();
    ret += ", ";
  }
  ret += scalVals_[size_-1]->ToString();
  ret += ")";
  return ret;
}
} // end of namespace
