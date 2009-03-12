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
      param->Get("sequenceStep", "index", GenStr( sequenceStep ) )
      ->Get("analysis")->Get("static");
    
    // Get save increment for restart file (optional)
    myNode->Get( "writeRestartInc", restartIncr_, false );

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
  void StaticDriver::SolveProblem(bool write_results, const std::string& comment)
  {
    // Set curent value of timestep and time step size in the mathParser
    domain->GetMathParser()->SetValue( MathParser::GLOB_HANDLER,
                                         "t", 0.0 );
    domain->GetMathParser()->SetValue( MathParser::GLOB_HANDLER,
                                         "dt", 0.0 );    
    domain->GetMathParser()->SetValue( MathParser::GLOB_HANDLER,
                                         "step", 0 );        

    // 'TimeStepping' is here the optimization iteration
    ptPDE_->GetSolveStep()->SetActTime(0.0);
    ptPDE_->GetSolveStep()->SetActStep(1);
    ptPDE_->GetSolveStep()->PreStepStatic();
    ptPDE_->GetSolveStep()->SolveStepStatic(comment);
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
    
    SETPROFILE("After Static Step");    
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
