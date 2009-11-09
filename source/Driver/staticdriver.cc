// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <string>

#include "staticdriver.hh"
#include "stdSolveStep.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "PDE/StdPDE.hh"
//#include "PDE/pdememento.hh"
#include "Domain/domain.hh"
#include "DataInOut/resultHandler.hh"

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
    ParamNode * myNode = 
      param->Get("sequenceStep", std::string("index"), sequenceStep)
      ->Get("analysis")->Get("static");
    
    // Get save increment for restart file (optional)
    myNode->Get( "writeRestartInc", restartIncr_, false );

    driverNode = driverNode->Get("static");
    driverNode->Get("sequenceStep")->SetValue(sequenceStep);
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
  void StaticDriver::SolveProblem(bool write_results, InfoNode* given_analysis_id, const bool reAssembleMatrices)
  {
    // Set curent value of timestep and time step size in the mathParser
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
      analysis_id_ = driverNode->Get(InfoNode::PROCESS)->Get("step", InfoNode::APPEND);
      analysis_id_->Get("analysis_id")->SetValue(0);
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
    ptPDE_->GetSolveStep()->SolveStepStatic(analysis_id_, reAssembleMatrices);
    ptPDE_->GetSolveStep()->PostStepStatic();

    // in optimization we write the results via StoreResults() because
    // we don't write every forward step. 
    if(write_results)
    {
      handler_->BeginMultiSequenceStep( sequenceStep_, analysis_, 1);      
      StoreResults(1);
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

  void StaticDriver::StoreResults(double step_val)
  {
    assert(analysis_ == BasePDE::STATIC);

    handler_->BeginStep((unsigned int) step_val, step_val);
    ptPDE_->WriteResultsInFile((unsigned int) step_val, step_val);
    ptPDE_->WriteGeneralPDEdefines();
    handler_->FinishStep();
  }

} // end of namespace
