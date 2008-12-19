// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "bcs.hh"
#include "General/environment.hh"
#include "DataInOut/ParamHandling/InfoNode.hh"

namespace CoupledField {


  HomDirichletBc::HomDirichletBc() {
    dof = 1;
  }

  void HomDirichletBc::ToInfo(InfoNode* in) const
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

  InhomNeumannBc::InhomNeumannBc() {
  }

  void InhomDirichletBc::ToInfo(InfoNode* in) const
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

  LoadBc::LoadBc() {
  }

  void LoadBc::ToInfo(InfoNode* in) const
  {
    HomDirichletBc::ToInfo(in);
    in->Get("value")->SetValue(value);
    in->Get("phase")->SetValue(phase);
    in->Get("weight")->SetValue(weight);
  }


  std::string LoadBc::ToString()
  {
    std::stringstream ss;
    ss << HomDirichletBc::ToString() << ", value=" << value << ", phase=" << phase
       << ", weight=" << weight;
    return ss.str();
  }


  Constraint::Constraint() {
    masterDof = 1;
    slaveDof = 1;
  }
}
