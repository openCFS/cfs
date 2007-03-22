// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <string>
#include <iostream>

#include "matvec/stdvector.hh"

namespace OLAS {

  // =======================
  //   Default Constructor
  // =======================
  StdVector::StdVector() {
    ENTER_FCN( "StdVector::StdVector" );
    size_ = 0;
  }
  
  // ======================
  //   Default Destructor
  // ======================
  StdVector::~StdVector() {
    ENTER_FCN( "StdVector::~StdVector" );
  }


}
