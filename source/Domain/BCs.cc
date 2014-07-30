#include "BCs.hh"
#include "General/Environment.hh"
#include "Domain/ElemMapping/EntityLists.hh"
#include "Domain/CoefFunction/CoefFunction.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField {


  HomDirichletBc::HomDirichletBc() {
  }

  void HomDirichletBc::ToInfo(PtrParamNode in) const
  {
    in->Get("name")->SetValue(entities->GetName() );
    std::set<UInt>::const_iterator it = dofs.begin();
    StdVector<UInt> dofVec;
    for( ; it != dofs.end(); ++it ) {
      dofVec.Push_back(*it);
    }
    in->Get("dofs")->SetValue(dofVec.Serialize());
  }



  HomDirichletBc::~HomDirichletBc()
  {
     // we don't do anything in this (virtual) destructor, it's just there because
     // of the virtual Dump()
  }


  InhomDirichletBc::InhomDirichletBc() {
  }
  
  InhomDirichletBc::~InhomDirichletBc() {
    value.reset();
  }

  void InhomDirichletBc::ToInfo(PtrParamNode in) const
  {
    HomDirichletBc::ToInfo(in);
    in->Get("value")->SetValue(value->ToString());
  }


  Constraint::Constraint() :
    masterDof(1),
    slaveDof(1),
    periodic(false)
  {
  }
}
