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

namespace CoupledField {

  // ***************
  //   Constructor
  // ***************
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
  }

  void StaticDriver::Init() {

    // Initialize first multisequence step, as the method "CheckStoreResults" 
    // relies on the result handler to know already about the current
    // sequencestep. However, in case of optimization, the sequence step
    // gets initialized in Optimization::SolveProblem()
    handler_->BeginMultiSequenceStep( sequenceStep_, analysis_, 1);

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
    // Set current value of time step and time step size in the mathParser
    domain->GetMathParser()->SetValue( MathParser::GLOB_HANDLER,
                                         "t", 0.0 );
    domain->GetMathParser()->SetValue( MathParser::GLOB_HANDLER,
                                         "dt", 0.0 );    
    domain->GetMathParser()->SetValue( MathParser::GLOB_HANDLER,
                                         "step", 0 );        

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
    
    // 'TimeStepping' is here the optimization iteration
    ptPDE_->GetSolveStep()->SetActTime(0.0);
    ptPDE_->GetSolveStep()->SetActStep(1);
    ptPDE_->GetSolveStep()->PreStepStatic();
    ptPDE_->GetSolveStep()->SolveStepStatic(analysis_id_, adjointParams);
    ptPDE_->GetSolveStep()->PostStepStatic();

    // in optimization we write the results via StoreResults() because
    // we don't write every forward step. 
    if(write_results)
    {
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
