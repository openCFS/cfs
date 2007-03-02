#include <fstream>
#include <iostream>
#include <string>

#include "staticdriver.hh"
#include "stdSolveStep.hh"

#include "PDE/StdPDE.hh"
#include "Domain/domain.hh"
#include "DataInOut/resultHandler.hh"

namespace CoupledField {

  // ***************
  //   Constructor
  // ***************
  StaticDriver::StaticDriver( UInt stepOffset,
                              Double timeOffset,
                              UInt sequenceStep, 
                              bool isPartOfSequence ) 
    : SingleDriver( sequenceStep, isPartOfSequence ) {
    
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

    ResultHandler * resHandler = domain->GetResultHandler();

    // Initialize 'TimeStepping'
    const UInt nstep = 1;
    Double  steptime = 0.0;

    // notify resultHandler about beginning of new sequence step 
    if( !isPartOfSequence_ )
      resHandler->BeginMultiSequenceStep( 1, analysis_ );
    
    

    ptPDE_->GetSolveStep()->SetActTime(steptime);
    ptPDE_->GetSolveStep()->SetActStep(nstep);
    ptPDE_->WriteGeneralPDEdefines();
    ptPDE_->GetSolveStep()->PreStepStatic();
    ptPDE_->GetSolveStep()->SolveStepStatic();
    ptPDE_->GetSolveStep()->PostStepStatic();

    resHandler->BeginStep( nstep+stepOffset_, timeOffset_ + steptime );
    ptPDE_->WriteResultsInFile(nstep, steptime, stepOffset_, timeOffset_);
    resHandler->FinishStep();
    
    ptPDE_->Finalize();

    // notify resultHandler about finishing of current sequence step
    if( !isPartOfSequence_ )
      resHandler->FinishMultiSequenceStep();
    SETPROFILE("After Static Step");
  }


} // end of namespace
