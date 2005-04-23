#include "basePDE.hh"

#include "Driver/baseSolveStep.hh"

namespace CoupledField {
  
  BasePDE::BasePDE() {
    
    ENTER_FCN( "BasePDE::BasePDE" );
    
    bcSequenceIndex_ = 0;
    bcSequenceTag_   = "anyTag";
  }

  // **********************
  //   Default Destructor
  // **********************
  BasePDE::~BasePDE() {
    ENTER_FCN( "BasePDE::~BasePDE" );
  }
  

} // end of namespace

