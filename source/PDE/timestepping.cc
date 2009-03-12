// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <string>

#include "timestepping.hh"

#include "OLAS/algsys/basesystem.hh"

namespace CoupledField {

  TimeStepping::TimeStepping(BaseSystem * algebraicsystem )
  {

    algsys_  = algebraicsystem;
    rhsSize_ = 0 ;
    isDeriv1Set_ = false;
    isDeriv2Set_ = false;
    isSolTN1Set_ = false;
  }

  const std::map<FEMatrixType,Double>  &
  TimeStepping::GetEffSysMatFactors( ) const {
    return matrix_factors_;
  }
                                              

  TimeStepping::~TimeStepping()
  {
  }

  bool TimeStepping::FeMatrixPresent( FEMatrixType type) {
    std::set<FEMatrixType> matTypes;    
    algsys_->GetFEMatrixTypes(matTypes);
    
    if ( matTypes.find(type) != matTypes.end() ) {
      return true;
    } else {
      return false;
    }
  }
     


} // end of namespace
