// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

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
    

    analysis_ = STATIC;

    stepOffset_ = stepOffset;
    timeOffset_ = timeOffset;
    consecutiveRun_ = false; 
  }

  void StaticDriver::Init() {

    InitializePDEs();
  }


  // **********************
  //   Default destructor
  // **********************
  StaticDriver::~StaticDriver() {
  }


  // *****************
  //   Solve problem
  // *****************
  void StaticDriver::SolveProblem(int optimizationIteration) {
 
    lastOptimizationIteration_ = optimizationIteration;
    ResultHandler * resHandler = domain->GetResultHandler();

    // notify resultHandler about beginning of new sequence step 
    resHandler->BeginMultiSequenceStep( sequenceStep_, analysis_, 1 );
    
    // 'TimeStepping' is here the optimization iteration
    ptPDE_->GetSolveStep()->SetActTime(optimizationIteration);
    ptPDE_->GetSolveStep()->SetActStep(optimizationIteration);
    ptPDE_->WriteGeneralPDEdefines();
    ptPDE_->GetSolveStep()->PreStepStatic();
    ptPDE_->GetSolveStep()->SolveStepStatic();
    ptPDE_->GetSolveStep()->PostStepStatic();
  }

  void StaticDriver::StoreResults(double step_val)
  {
    // post process ??
    int loi = lastOptimizationIteration_;

    ResultHandler * resHandler = domain->GetResultHandler();

    // resHandler->BeginStep( nstep+stepOffset_, timeOffset_ + steptime );
    resHandler->BeginStep( loi + stepOffset_, step_val > 0 ? step_val : timeOffset_ + (loi-1) );
    // actually the parameters seem to be not used :(
    ptPDE_->WriteResultsInFile(loi , loi-1, stepOffset_, timeOffset_);
    resHandler->FinishStep();

    ptPDE_->Finalize();
    
    // notify resultHandler about finishing of current sequence step
    if(!optimization_ && !isPartOfSequence_ ) {
      resHandler->FinishMultiSequenceStep();
      resHandler->Finalize();
    }
    SETPROFILE("After Static Step");
    
  }

} // end of namespace
