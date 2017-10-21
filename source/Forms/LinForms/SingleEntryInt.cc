#include "SingleEntryInt.hh"
#include "Domain/Domain.hh"
#include "Domain/CoefFunction/CoefFunctionOpt.hh"

namespace CoupledField {


 SingleEntryInt::SingleEntryInt( PtrCoefFct& val )
    : LinearForm() {

    name_ = "SingleEntryInt";

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
    if(val_->GetDimType() == CoefFunction::SCALAR)
    {
      elemVec.Resize(1);
      val_->GetScalar(elemVec[0], lpm);
      return;
    }
    if(val_->GetDimType() == CoefFunction::VECTOR)
    {
      // if we have design dependent loads for optimization
      if(typeid(*val_) == typeid(CoefFunctionOpt))
      {
        assert(ent1.GetType() == EntityList::NODE_LIST);
        elemVec.Resize(ent1.GetSize()); // resize to number of nodes

        LocPoint lp;
        lp.number = ent1.GetNode();
        lp.coord.Resize(3);
        lpm.lp = lp;
        // we don't set an ElemShapeMap as we are here purely in the regime of nodes an use lpm only to transport
        // the node to DesignSpace::ApplyPhysicalDesign(... Vector ...)
      }

      val_->GetVector(elemVec, lpm);
      return;
    }
    assert(false); // SingleEntryInt only works for SCALAR and VECTOR
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
