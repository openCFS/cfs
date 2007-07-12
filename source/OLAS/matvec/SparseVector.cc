// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <string>
#include <iostream>

#include "matvec/SparseVector.hh"

namespace OLAS {

  // =======================
  //   Default Constructor
  // =======================
  SparseVector::SparseVector() {
    ENTER_FCN( "SparseVector::SparseVector" );
    size_ = 0;
  }
  
  // ======================
  //   Default Destructor
  // ======================
  SparseVector::~SparseVector() {
    ENTER_FCN( "SparseVector::~SparseVector" );
  }


}
