// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "basePDE.hh"
#include "Driver/baseSolveStep.hh"


namespace CoupledField {


  // ***********************
  //   Default Constructor
  // ***********************
  BasePDE::BasePDE( ParamNode* paramNode ) {
    
    ENTER_FCN( "BasePDE::BasePDE" );
    
    myParam_ = paramNode;
    sequenceStep_ = 0;
  }

  // **************
  //   Destructor
  // **************
  BasePDE::~BasePDE() {
    ENTER_FCN( "BasePDE::~BasePDE" );
  }

}
