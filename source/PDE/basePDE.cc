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
