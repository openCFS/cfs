#include <fstream>
#include <iostream>
#include <string>

#include "StaticDriver.hh"
#include "Driver/SolveSteps/StdSolveStep.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ProgramOptions.hh"
#include "PDE/StdPDE.hh"
#include "Domain/Domain.hh"
#include "DataInOut/ResultHandler.hh"
#include "DataInOut/SimState.hh"

namespace CoupledField {

  // ***************
  //   Constructor
  // ***************
  StaticDriver::StaticDriver( UInt sequenceStep,
                              bool isPartOfSequence,
                              shared_ptr<SimState> state, Domain* domain,
                              PtrParamNode paramNode, PtrParamNode infoNode ) 
    : SingleDriver( sequenceStep, isPartOfSequence, state, domain, paramNode, infoNode )
  {
    analysis_ = BasePDE::STATIC;
    param_ = param_->Get("static");
    info_ = info_->Get("static");
   
    // read flag if all results should get written to database file section
    // to allow e.g. for general postprocessing or result extraction
    param_->GetValue("allowPostProc", writeAllSteps_, ParamNode::PASS );
  }

  void StaticDriver::Init(bool restart)
  {
    InitializePDEs();

    // Initialize first multisequence step, as the method "CheckStoreResults"
    // relies on the result handler to know already about the current
    // sequencestep. However, in case of optimization, the sequence step
    // gets initialized in Optimization::SolveProblem()
    if(!domain->GetOptimization())
      handler_->BeginMultiSequenceStep( sequenceStep_, analysis_, 1);
  }


  // **********************
  //   Default destructor
  // **********************
  StaticDriver::~StaticDriver() {
  }


  // *****************
  //   Solve problem
  // *****************
  void StaticDriver::SolveProblem(bool write_results, PtrParamNode given_analysis_id)
  {

    // in the optimization case the step is given, otherwise it is created
    // store such that special steps can add non-lin stuff and optimization adjoints
    if(given_analysis_id == NULL)
    {
      // do we really want to create a new entry? Might blast up the output
      ParamNode::ActionType at = progOpts->DoDetailedInfo() ? ParamNode::APPEND : ParamNode::DEFAULT;
      analysis_id_ = info_->Get(ParamNode::PROCESS)->Get("step", at);
      analysis_id_->Get("analysis_id")->SetValue("0");
    }
    else
    {
      analysis_id_ = given_analysis_id;
      assert(analysis_id_->Has("analysis_id"));
    }

    // In case we allow general postprocessing or this analysis is part of 
    // a multisequence (in which case the subsequent run could need this
    // simulation as restart information)
    if(writeAllSteps_ || isPartOfSequence_)
      simState_->BeginMultiSequenceStep( sequenceStep_, analysis_ );

    // set dummy time to zero
    mathParser_->SetValue( MathParser::GLOB_HANDLER, "t", 0.0 );
    
    // 'TimeStepping' is here the optimization iteration
    ptPDE_->GetSolveStep()->SetActTime(0.0);
    ptPDE_->GetSolveStep()->SetActStep(1);
    ptPDE_->GetSolveStep()->PreStepStatic();
    ptPDE_->GetSolveStep()->SolveStepStatic(analysis_id_);
    ptPDE_->GetSolveStep()->PostStepStatic();

    // in optimization we write the results via StoreResults() because
    // we don't necessarily write every forward step.
    if(write_results)
    {
      StoreResults(1,0.0);
      handler_->FinishMultiSequenceStep();

      if(!isPartOfSequence_)
        handler_->Finalize(); // to be called only once in a HDF5 lifetime!
    }
    
    if(writeAllSteps_ || isPartOfSequence_ ) { 
      simState_->FinishMultiSequenceStep(true);
    }
  }
  
  void StaticDriver::SetToStepValue(UInt stepNum, Double stepVal ) {
    // ensure that this method is only called if simState has input
    if( ! simState_->HasInput()) {
      EXCEPTION( "Can only set external time step, if simulation state "
              << "is read from external file" );
    }
    
    if( stepNum != 1 || stepVal > EPS ) {
      EXCEPTION( "A static driver only has 1 step and its time value is 0.0!")
    }
    
    
    // For a static driver only stepNum == 1 and stepVal == 0.0 make sense
    mathParser_->SetValue( MathParser::GLOB_HANDLER, "t", 0.0 );
    mathParser_->SetValue( MathParser::GLOB_HANDLER, "step", 1 );    
  }

  void StaticDriver::StoreResults(UInt stepNum, double step_val)
  {
    assert(analysis_ == BasePDE::STATIC);

    handler_->BeginStep(stepNum, step_val);
    ptPDE_->WriteResultsInFile(stepNum, step_val);
    ptPDE_->WriteGeneralPDEdefines();
    handler_->FinishStep();
    if( writeAllSteps_ || isPartOfSequence_ )
      simState_->WriteStep(stepNum, step_val );
  }

} // end of namespace
