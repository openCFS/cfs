#include "CoefFunctionDiagTensorFromScalar.hh"

#include <limits>

namespace CoupledField {

  CoefFunctionDiagTensorFromScalar::
  CoefFunctionDiagTensorFromScalar(const StdVector<PtrCoefFct>& scalVals, std::string subType)
 : CoefFunction() {

   // this type of coefficient is never analytic
   isAnalytic_ = false;
   subType_ = subType;

   if( subType_ == "" ) {
     dimType_ = TENSOR;
   } else {
     // if we specify a subtype we assume that we want to have the output in voigt-notation
     dimType_ = VECTOR;
   }

   // Loop over all entries and check if all entries have the same
   assert(scalVals.GetSize() > 0 );
   // For the vector-case we define the correct size depending on the sub-type
   if( subType_ == "3d" ) {
      //Components: "xx", "yy", "zz", "yz", "xz", "xy"
      size_ = 2*scalVals.GetSize();
    } else if( subType_ == "plane" ) {
      //Components: "xx", "yy", "xy"
      size_ = scalVals.GetSize()+1;
    } else if( subType_ == "axi" ) {
      //Components = "rr", "zz", "rz", "phiphi"
      size_ = 2*scalVals.GetSize();
    } else if( subType_ == ""){
      size_ = scalVals.GetSize();
    } else {
      EXCEPTION("Unkown subtype, can't convert to voigt notation")
    }
   dependType_ = scalVals[0]->GetDependency();
   isComplex_ = scalVals[0]->IsComplex();
   for( UInt i = 1; i < scalVals.GetSize(); ++i ) {
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
