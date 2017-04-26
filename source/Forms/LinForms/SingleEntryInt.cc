#include "SingleEntryInt.hh"
#include "Domain/Domain.hh"


namespace CoupledField {


 SingleEntryInt::SingleEntryInt( PtrCoefFct& val )
    : LinearForm() {

    name_ = "SingleEntryInt";
    // check, if we have a constant expression coefficient function
    if((val->GetDependency() != CoefFunction::CONSTANT) && (val->GetDependency() != CoefFunction::TIMEFREQ)) {
      EXCEPTION("SingleEntryInt only works with space independent coefficients");
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

      if (val_->GetDependency() == CoefFunction::GENERAL) { // if we have design dependent loads
          assert(ent1.GetType() == EntityList::NODE_LIST);
          elemVec.Resize(ent1.GetSize()); // resize to number of nodes

          LocPoint lp;
          lp.number = ent1.GetNode();
          Vector<double> coords(3);
          lp.coord = coords;
          shared_ptr<ElemShapeMap> esm = ent1.GetGrid()->GetElemShapeMap( ent1.GetGrid()->GetElem(1), this->coordUpdate_ );
          lpm.Set(lp, esm,1.0);
        }

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
    // std::cout << "val_->GetDimType(): " << CoefFunction::CoefDimType_.ToString(val_->GetDimType()) << std::endl;
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
