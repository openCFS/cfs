#include "SingleEntryInt.hh"
#include "Domain/Domain.hh"


namespace CoupledField {


 SingleEntryInt::SingleEntryInt( PtrCoefFct& val )
    : LinearForm() {

    name_ = "SingleEntryInt";
    
    // check, if we have a constant expression coefficient function
    if( val->GetDependency() == CoefFunction::GENERAL) {
      EXCEPTION("SingleEntryInt only works with constant coefficients");
    }
    val_ = val;

  }

 SingleEntryInt::SingleEntryInt(const SingleEntryInt& right )
   : LinearForm(right){
   this->val_ = right.val_;
 }

  SingleEntryInt::~SingleEntryInt() {
  }
  
  void SingleEntryInt::CalcElemVector( Vector<Double>& elemVec,
                                       EntityIterator& ent1) {
    
    // we use just a dummy local point, as we assume constant
    // expression coefficient function
    LocPointMapped lpm; 
    if( val_->GetDimType() == CoefFunction::SCALAR) {
      elemVec.Resize(1);
      val_->GetScalar(elemVec[0], lpm);
    } else  if( val_->GetDimType() == CoefFunction::VECTOR) {
      val_->GetVector(elemVec, lpm);
    } else {
      EXCEPTION( "SingleEntryInt only works for SCALAR and VECTOR" );
    }
  }
  
  void SingleEntryInt::CalcElemVector( Vector<Complex>& elemVec,
                                       EntityIterator& ent1) {
    
    // we use just a dummy local point, as we assume constant
    // expression coefficient function
    LocPointMapped lpm; 
    if( val_->GetDimType() == CoefFunction::SCALAR) {
      elemVec.Resize(1);
      val_->GetScalar(elemVec[0], lpm);
    } else  if( val_->GetDimType() == CoefFunction::VECTOR) {
      val_->GetVector(elemVec, lpm);
    } else {
      EXCEPTION( "SingleEntryInt only works for SCALAR and VECTOR" );
    }
  }

  
}
