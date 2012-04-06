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


  SingleEntryInt::~SingleEntryInt() {
  }
  
  void SingleEntryInt::CalcElemVector( Vector<Double>& elemVec,
                                       EntityIterator& ent1) {
    
    // we use just a dummy local point, as we assume constant
    // expression coefficient function
    LocPointMapped lpm; 
    val_->GetVector(elemVec, lpm);
  }
  
  void SingleEntryInt::CalcElemVector( Vector<Complex>& elemVec,
                                       EntityIterator& ent1) {
    
    // we use just a dummy local point, as we assume constant
    // expression coefficient function
    LocPointMapped lpm; 
    val_->GetVector(elemVec, lpm);
  }

  
}
