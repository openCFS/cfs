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
