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
  StaticDriver::StaticDriver( UInt stepOffset,
                              Double timeOffset,
                              std::string driverTag,
                              bool isPartOfSequence ) 
    : SingleDriver( driverTag, isPartOfSequence ) {
    
    ENTER_FCN( "StaticDriver::StaticDriver" );

    analysis_ = STATIC;

    stepOffset_ = stepOffset;
    timeOffset_ = timeOffset;
    
  }

  void StaticDriver::Init() {
    ENTER_FCN( "StaticDriver::Init" );

    InitializePDEs();
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

    // Initialize 'TimeStepping'
    const UInt nstep = 1;
    Double  steptime = 0.0;

    ptPDE_->GetSolveStep()->SetActTime(steptime);
    ptPDE_->GetSolveStep()->SetActStep(nstep);
    ptPDE_->WriteGeneralPDEdefines();
    ptPDE_->GetSolveStep()->PreStepStatic();
    ptPDE_->GetSolveStep()->SolveStepStatic();
    ptPDE_->GetSolveStep()->PostStepStatic();

    ptPDE_->PostProcess();

    
    // if multiSequence is performed, the ms-driver
    // writes out the grid one time
    if (! isPartOfSequence_)
      domain->PrintGrid();
      
    ptPDE_->WriteResultsInFile(nstep, steptime, stepOffset_, timeOffset_);
    ptPDE_->WriteHistoryInFile(nstep, steptime, stepOffset_, timeOffset_);

    SETPROFILE("After Static Step");
  }


} // end of namespace
