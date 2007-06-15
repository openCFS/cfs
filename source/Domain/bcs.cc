// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "bcs.hh"
#include "General/environment.hh"

namespace CoupledField {


  HomDirichletBc::HomDirichletBc() {
    ENTER_IFCN( "HomDirichletBc::HomDirichletBc" );

    dof = 1;
  }

  HomDirichletBc::~HomDirichletBc() 
  {
     // we don't do anything in this (virtual) destructor, it's just there because
     // of the virtual Dump()
  }


  InhomDirichletBc::InhomDirichletBc() {
    ENTER_IFCN( "InhomDirichletBc::InhomDirichletBc");
  }

  InhomNeumannBc::InhomNeumannBc() {
    ENTER_IFCN( "InhomNeumannBc::InhomNeumannBc" );
  }
  
  LoadBc::LoadBc() {
    ENTER_IFCN( "LoadBc::LoadBc" );
  }
  
  Constraint::Constraint() {
    ENTER_IFCN( "Constraint::Constraint" );

    masterDof = 1;
    slaveDof = 1;
  }
}
