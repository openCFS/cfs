#include "BCs.hh"
#include "General/Environment.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField {


  HomDirichletBc::HomDirichletBc() {
    dof = 1;
  }

  void HomDirichletBc::ToInfo(PtrParamNode in) const
  {
    in->Get("dof")->SetValue(dof);
  }

  std::string HomDirichletBc::ToString()
  {
    std::stringstream ss;
    ss << "dof=" << dof;
    return ss.str();
  }


  HomDirichletBc::~HomDirichletBc()
  {
     // we don't do anything in this (virtual) destructor, it's just there because
     // of the virtual Dump()
  }


  InhomDirichletBc::InhomDirichletBc() {
  }

  void InhomDirichletBc::ToInfo(PtrParamNode in) const
  {
    HomDirichletBc::ToInfo(in);
    in->Get("value")->SetValue(value);
    in->Get("phase")->SetValue(phase);
  }


  std::string InhomDirichletBc::ToString()
  {
    std::stringstream ss;
    ss << HomDirichletBc::ToString() << ", value=" << value << ", phase=" << phase;
    return ss.str();
  }


  Constraint::Constraint() :
    masterDof(1),
    slaveDof(1),
    periodic(false)
  {
  }
}
