// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <assert.h>
#include <iostream>
#include <string>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/programOptions.hh"
#include "DataInOut/resultHandler.hh"
#include "Domain/domain.hh"
#include "Driver/baseSolveStep.hh"
#include "Driver/singleDriver.hh"
#include "General/environment.hh"
#include "PDE/basePDE.hh"
#include "PDE/LatticeBoltzmannPDE.hh"
#include "staticdriver.hh"

namespace CoupledField {

  // ***************
  //   Constructor
  // ***************
class AdjointParameters;

  StaticDriver::StaticDriver( UInt sequenceStep,
                              bool isPartOfSequence ) 
    : SingleDriver( sequenceStep, isPartOfSequence ) {

    analysis_ = BasePDE::STATIC;
    restartIncr_= 0;

    // get parameter node
    PtrParamNode myNode = 
      param->GetByVal("sequenceStep", std::string("index"), sequenceStep)
      ->Get("analysis")->Get("static");
    
    // Get save increment for restart file (optional)
    myNode->GetValue( "writeRestartInc", restartIncr_, ParamNode::PASS );

    driverNode = driverNode->Get("static");
    driverNode->Get("sequenceStep")->SetValue(sequenceStep);
    lbm_ = false;
  }

  void StaticDriver::Init() {

    // Initialize first multisequence step, as the method "CheckStoreResults" 
    // relies on the result handler to know already about the current
    // sequence step. However, in case of optimization, the sequence step
    // gets initialized in Optimization::SolveProblem()
    if( !domain->GetOptimization()) {
//      if (lbm_)
        handler_->BeginMultiSequenceStep( sequenceStep_, BasePDE::TRANSIENT, 9999);
//      else
//        handler_->BeginMultiSequenceStep( sequenceStep_, analysis_, 1);
    }
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
  void StaticDriver::SolveProblem(bool write_results, PtrParamNode given_analysis_id, AdjointParameters* adjointParams)
  {
    if(ptPDE_->GetName() == "LatticeBoltzmann")
      lbm_ = true;
    // in the optimization case the step is given, otherwise it is created
    // store such that special steps can add non-lin stuff and optimization adjoints
    if(given_analysis_id == NULL)
    {
      // do we really want to create a new entry? Might blast up the output
      ParamNode::ActionType at = progOpts->DoDetailedInfo() ? ParamNode::APPEND : ParamNode::DEFAULT;
      analysis_id_ = driverNode->Get(ParamNode::PROCESS)->Get("step", at);
      analysis_id_->Get("analysis_id")->SetValue("0");
    }
    else
    {
      analysis_id_ = given_analysis_id;
      assert(analysis_id_->Has("analysis_id"));
    }
    
    //check, if we are just interested in getting the data for the approximation 
    //of nonlinear curves
    if ( !progOpts->DoApproxNLmatData() ) {
    
      // 'TimeStepping' is here the optimization iteration
      ptPDE_->GetSolveStep()->SetActTime(0.0);
      ptPDE_->GetSolveStep()->SetActStep(1);
      ptPDE_->GetSolveStep()->PreStepStatic();
      ptPDE_->GetSolveStep()->SolveStepStatic(analysis_id_, adjointParams);
      ptPDE_->GetSolveStep()->PostStepStatic();
    }

    // for LBM case we overwrite this data because we might write intermediate LBM iterations
    int step = 1;
    double value = 0.0; // step value
    if(lbm_) {
      dynamic_cast<LatticeBoltzmannPDE*>(ptPDE_)->Solve(); // might call many StoreResults() for intermediate steps
      step = dynamic_cast<LatticeBoltzmannPDE*>(ptPDE_)->GetNumWriteResults();
      value = step;
    }

    // in optimization we write the results via StoreResults() because
    // we don't write every forward step. 
    if(write_results)
      {
        if (lbm_)
          StoreResults(step+1,(double)value+1.0);
        else
          StoreResults(1,0.0);
        handler_->FinishMultiSequenceStep();

        if(!isPartOfSequence_)
          handler_->Finalize(); // to be called only once in a HDF5 lifetime!
      }
    // writing current PDE-state into the restart-file
    if (restartIncr_ >= 1){
      std::cout << std::endl << " *** Write a restart file *** " << std::endl;      
      ptPDE_->WriteRestart( );
    }
  }

  void StaticDriver::StoreResults(UInt stepNum, double step_val)
  {
    assert(analysis_ == BasePDE::STATIC);

    handler_->BeginStep(stepNum, step_val);
    ptPDE_->WriteResultsInFile(stepNum, step_val);
    ptPDE_->WriteGeneralPDEdefines();
    handler_->FinishStep();
  }

} // end of namespace
