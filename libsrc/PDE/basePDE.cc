#include "basePDE.hh"

#include "Driver/baseSolveStep.hh"

namespace CoupledField {

  BasePDE::BasePDE() {

    ENTER_FCN( "BasePDE::BasePDE" );

    solveStep_       = NULL;
    bcSequenceIndex_ = 0;
    bcSequenceTag_   = "anyTag";
  }
 

  // **********************
  //   Default Destructor
  // **********************
  BasePDE::~BasePDE() {
    ENTER_FCN( "BasePDE::~BasePDE" );
    delete solveStep_;
  }

} // end of namespace

