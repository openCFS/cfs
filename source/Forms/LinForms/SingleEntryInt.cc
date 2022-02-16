#include "SingleEntryInt.hh"
#include "Domain/Domain.hh"
#include "Domain/CoefFunction/CoefFunctionOpt.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

DEFINE_LOG(entryint, "entryint")

namespace CoupledField
{


SingleEntryInt::SingleEntryInt(PtrCoefFct& val) : LinearForm()
 {
    name_ = "SingleEntryInt";
    val_ = val;
  }

 SingleEntryInt::SingleEntryInt(const SingleEntryInt& right) : LinearForm(right)
 {
   this->val_ = right.val_;
 }

  SingleEntryInt::~SingleEntryInt() {
  }
  
  const std::string SingleEntryInt::ToString()
  {
    return GetName() + " coef=" + val_->GetName() + " val=" + val_->ToString();

  }

  void SingleEntryInt::CalcElemVector(Vector<Double>& elemVec, EntityIterator& ent1)
  {
    LOG_DBG2(entryint) << "SEI:CEV: coef=" << val_->GetName() << " val=" << val_->ToString() << " ent1=" << ent1.ToString()
        << " size=" << ent1.GetSize() << " dep=" << CoefFunction::coefDependType.ToString(val_->GetDependency());
    
    // not that elemVec gets here usually (e.g. mech rhs force) the nodal (vectorial) value
    // but when we have a NODE_LIST the LPM::shapeMap cannot be set.
    // usually we have coef const, the lpm is not used there anyway and just required for the function :(
    // special cases are optimization with design dependend loads or coef expression (force is sin(x))
    // we then set the coordinate and hope the coef::GetVector() can handle it
    // (implemented for CoefFunctionExp::GetVector())

    LocPointMapped lpm; 
    if(val_->GetDimType() == CoefFunction::SCALAR)
    {
      elemVec.Resize(1);
      val_->GetScalar(elemVec[0], lpm);
      return;
    }
    if(val_->GetDimType() == CoefFunction::VECTOR)
    {
      // if we have design dependent loads for optimization we need to set the lp number
      // for spatial displacement we set the coordinate
      if(typeid(*val_) == typeid(CoefFunctionOpt) || val_->IsSpacialDependent())
      {
        assert(ent1.GetType() == EntityList::NODE_LIST);
        lpm.lp.number = ent1.GetNode();
        assert(!lpm.shapeMap); // having no shapeMap makes usage of lpm.lp
      }
      val_->GetVector(elemVec, lpm); // gets nodal coordinate from lp.number via lpm::GetCoord() with fallbac

      LOG_DBG2(entryint) << "SEI:CEV: res=" << elemVec.ToString() << " coord=" << lpm.lp.coord.ToString();
      return;
    }
    assert(false); // SingleEntryInt only works for SCALAR and VECTOR
  }
  
  void SingleEntryInt::CalcElemVector( Vector<Complex>& elemVec,
                                       EntityIterator& ent1) {

    // we use just a dummy local point, as we assume constant
    // expression coefficient function
    LocPointMapped lpm;
    // std::cout << "val_->GetDimType(): " << CoefFunction::coefDimType.ToString(val_->GetDimType()) << std::endl;
    if( val_->GetDimType() == CoefFunction::SCALAR) {
      elemVec.Resize(1);
      val_->GetScalar(elemVec[0], lpm);
    } else  if( val_->GetDimType() == CoefFunction::VECTOR) {
      if((typeid(*val_) == typeid(CoefFunctionOpt) || val_->IsSpacialDependent()) && ent1.GetType() == EntityList::NODE_LIST)
        lpm.lp.number = ent1.GetNode();
      val_->GetVector(elemVec, lpm);
    } else {
      EXCEPTION( "SingleEntryInt only works for SCALAR and VECTOR" );
    }
  }

  
}
