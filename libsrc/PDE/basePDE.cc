#include "basePDE.hh"
#include "Driver/baseSolveStep.hh"


namespace CoupledField {


  // ***********************
  //   Default Constructor
  // ***********************
  BasePDE::BasePDE() {
    
    ENTER_FCN( "BasePDE::BasePDE" );
    
    bcSequenceIndex_ = 0;
    bcSequenceTag_   = "anyTag";
  }

  // **************
  //   Destructor
  // **************
  BasePDE::~BasePDE() {
    ENTER_FCN( "BasePDE::~BasePDE" );
  }

}
