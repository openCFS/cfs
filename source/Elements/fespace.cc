// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "fespace.hh"

namespace CoupledField {

  //! Constructor
  FeSpace::FeSpace(){
    //build up the pointerMap
    refElems_[Elem::LINE2] = ptL1;
    refElems_[Elem::QUAD4] = ptQ1;
    refElems_[Elem::HEXA8] = ptHexa1;
    isFinalized_ = false;
    numEqns_ = 0;
    numUnknowns_ = 0;
    numFreeEquations_ = 0;

    bcCounter_[NOBC] = 0;
    bcCounter_[HDBC] = 0;
    bcCounter_[IDBC] = 0;
    bcCounter_[CONSTRAINT] = 0;
  }

  FeSpace::~FeSpace(){
  }

}
