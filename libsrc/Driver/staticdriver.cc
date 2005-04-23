#include <fstream>
#include <iostream>
#include <string>

#include "staticdriver.hh"
#include "stdSolveStep.hh"

#include "PDE/StdPDE.hh"
#include "Domain/domain.hh"

namespace CoupledField {

  // ***************
  //   Constructor
  // ***************
  StaticDriver::StaticDriver(Domain * adomain, 
			     Integer stepOffset,
			     Double timeOffset,
			     std::string driverTag,
			     Boolean isPartOfSequence) 
    : SingleDriver(adomain, stepOffset, timeOffset, 
		   driverTag, isPartOfSequence) {
    ENTER_FCN( "StaticDriver::StaticDriver" );

  }


  // **********************
  //   Default destructor
  // **********************
  StaticDriver::~StaticDriver() {
    ENTER_FCN( "StaticDriver::~StaticDriver" );
  }


  // *****************
  //   Solve problem
  // *****************
  void StaticDriver::SolveProblem() {
    ENTER_FCN( "StaticDriver::SolveProblem" );

    Integer level=0;

    // if driver is not part of multiSequence Driver, get list
    // of pdes which have to be solved and intialize them
    if (isPartOfSequence_ == FALSE){     
      GetMyPDEs();
      Info->StartProgress ("Starting to solve problem", FALSE);
    }

    // Initialize 'TimeStepping'
    const Integer nstep = 1;
    Double  steptime = 0.0;
    Boolean reset = FALSE;
    
    ptPDE_->WriteGeneralPDEdefines();
    ptPDE_->GetSolveStep()->PreStepStatic(nstep, steptime, level, reset);
    ptPDE_->GetSolveStep()->SolveStepStatic(nstep, steptime, level,reset);
    ptPDE_->GetSolveStep()->PostStepStatic(nstep, steptime, level);
    
    ptPDE_->PostProcess(level);
    
    // if multiSequence is performed, the ms-driver
    // writes out the grid one time
    if (! isPartOfSequence_)
	ptdomain_->PrintGrid(level);
      
    ptPDE_->WriteResultsInFile(nstep, steptime, stepOffset_, timeOffset_);
  }


} // end of namespace
