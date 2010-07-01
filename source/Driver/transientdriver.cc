// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <string>
#include <stdio.h>
#include <list>
#include <math.h>


#include <boost/filesystem.hpp>

#include "transientdriver.hh"
#include "stdSolveStep.hh"

#include "DataInOut/programOptions.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "PDE/StdPDE.hh"
#include "PDE/pdememento.hh"
#include "PDE/SinglePDE.hh"
#include "Domain/domain.hh"
#include "DataInOut/resultHandler.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "Optimization/Optimization.hh"

using std::cout;
using std::endl;
namespace fs = boost::filesystem;

namespace CoupledField {

  DECLARE_LOG(trans_driver)
  DEFINE_LOG(trans_driver, "transient_driver")

  // ===============
  //   Constructor
  // ===============
  TransientDriver::TransientDriver( UInt sequenceStep,
                                    bool isPartOfSequence) 
    : SingleDriver( sequenceStep, isPartOfSequence ) {
    
    LOG_TRACE(trans_driver) << "TransientDriver():  "
                            << " sequenceStep: " << sequenceStep
                            << " isPartOfSequence: " << isPartOfSequence;
    // Set analysistype
    analysis_ = BasePDE::TRANSIENT;

    actTimeStep_ = 0;
    firstdt_ = 0.0;
    restartIncr_ = 0;
    restartStep_ = 0;

    // get parameter node
    PtrParamNode myNode = 
      param->GetByVal("sequenceStep", std::string("index"), sequenceStep)
      ->Get("analysis")->Get("transient");

    driverNode = driverNode->Get("transient");
    driverNode->Get("sequenceStep")->SetValue(sequenceStep);
    driverNode->Get(ParamNode::HEADER)->Get("unit")->SetValue("s");
    
    // Get time stepping information from parameter object
    myNode->GetValue( "numSteps", numstep_ );
    myNode->GetValue( "deltaT", firstdt_ );

    // Get save increment for restart file (optional)
    myNode->GetValue( "writeRestartInc", restartIncr_, ParamNode::PASS );
  
    // remove HALTCFS File at the beginning
    if(fs::exists("./HALTCFS")) 
      fs::remove("./HALTCFS");
  }

  // ==============
  //   Destructor
  // ==============
  TransientDriver::~TransientDriver()
  {
  }

  // ==================
  //   Initialization
  // ==================
  void TransientDriver::Init() 
  {
    InitializePDEs();
  }
    
  void TransientDriver::SolveProblem(bool write_results, PtrParamNode given_analysis_id, AdjointParameters* adjointParams, const bool reAssembleMatrices) 
  {
    // notify resultHandler about beginning of new sequence step 
    ResultHandler * resHandler = domain->GetResultHandler();

    int direction = adjointParams ? -1 : 1; // for adjoint the time is backwards
    UInt startStep = restartStep_ + 1;
    UInt endStep = restartStep_ + numstep_;
    Double  steptime  = firstdt_ * (adjointParams ? endStep : startStep);
    UInt nextRestart  = restartStep_ + restartIncr_ * direction;
    Double  dt = firstdt_;
    bool haltFlag=false;
    bool isFirstStep = true;
    
    Optimization* optimization = domain->GetOptimization();

    Double timeStepPercent = (double)numstep_/10;
    Double percentCounter = timeStepPercent;
    if(direction < 0){
      percentCounter = (double)numstep_ - timeStepPercent;
    }
  
    if(write_results){
      resHandler->BeginMultiSequenceStep( sequenceStep_,
                                          analysis_,
                                          numstep_+restartStep_-startStep+1 );
      if(optimization != NULL){ // we have to save everytime to a new multisequencestep
        sequenceStep_++;
      }
    }
  
    ptPDE_->WriteGeneralPDEdefines();

    ptPDE_->GetSolveStep()->SetStartStep( startStep );
    // Solve problem
    //ptPDE_->GetSolveStep()->SetTimeStep(dt);
    
    
    ptPDE_->GetSolveStep()->SetNumTimeSteps(restartStep_+numstep_);
    //---------------------------------------------------------------------------
    // to save the initial state
    // commented out, since all time-dependend examples in testsuite consider the first
    // not the zeroth time step
    //     resHandler->BeginStep( 0, 0 );
    //     ptPDE_->WriteResultsInFile(stepOffset_, timeOffset_);
    //     resHandler->FinishStep( );


    //---------------------------------------------------------------------------
    // Outer loop over all timesteps
    
    for (actTimeStep_ = adjointParams ? endStep : startStep; actTimeStep_ <= endStep && actTimeStep_ >= startStep; actTimeStep_ += direction) {
      LOG_DBG(trans_driver) << "loop over timestep " << actTimeStep_;
      // check for a HALTCFS File
      // if there exist a file with name HALTCFS in the executing directory
      // than CFS++ will end after the current time step
      std::ifstream readHALTCFS("HALTCFS", std::ios_base::in );
      if (readHALTCFS) {
        readHALTCFS.close();
        haltFlag = true;
        numstep_=actTimeStep_;
        ptPDE_->GetSolveStep()->SetNumTimeSteps(numstep_);
      }
      
      // Set curent value of timestep and time step size in the mathParser
      domain->GetMathParser()->SetValue( MathParser::GLOB_HANDLER,
                                         "t", steptime );
      domain->GetMathParser()->SetValue( MathParser::GLOB_HANDLER,
                                         "dt", dt );    
      domain->GetMathParser()->SetValue( MathParser::GLOB_HANDLER,
                                         "step", actTimeStep_ );    

      // Determine when to write logging information on terminal
      bool log = false;
      if(optimization != NULL){
        cout << ".";
        cout.flush();
      }else{
        if(numstep_ <= 50)
          log = true;
        if(numstep_ > 50 && numstep_ <= 500 && ! (actTimeStep_ % 10)) // every tenth step, not all but the tenth  
          log = true;
        if(numstep_ > 500 && ( (direction > 0 &&  (double) actTimeStep_ >= percentCounter) || (direction < 0 && (double) actTimeStep_ <= percentCounter) ) ){
          log = true;
          percentCounter += timeStepPercent * direction;
        }

        if(log)
        {
          if(progOpts->IsQuiet())
            cout << ptPDE_->GetName() << ": Time step " << actTimeStep_ << " time " << steptime << endl; 
          else
            // write std::out info    
            cout << endl << ptPDE_->GetName() << ": Time step " 
            << actTimeStep_ <<" ======================= " << endl;
        }
      }
      
      if(given_analysis_id == NULL)
      {
        // do we really want to create a new entry? Might blast up the output
        ParamNode::ActionType at = progOpts->DoDetailedInfo() ? ParamNode::APPEND : ParamNode::DEFAULT;
        analysis_id_ = driverNode->Get(ParamNode::PROCESS)->Get("step", at);
        analysis_id_->Get("analysis_id")->SetValue(actTimeStep_);
      }
      else
      {
        assert(given_analysis_id->Has("analysis_id"));
        analysis_id_ = CreateAnalysisIdChild(given_analysis_id, "step", actTimeStep_);
      }
      analysis_id_->Get("step")->SetValue(actTimeStep_);
      analysis_id_->Get("value")->SetValue(steptime);
      
      // Perform actions
      ptPDE_->GetSolveStep()->SetActTime(steptime);
      ptPDE_->GetSolveStep()->SetActStep(actTimeStep_);
      ptPDE_->GetSolveStep()->PreStepTrans();
      ptPDE_->GetSolveStep()->SolveStepTrans(analysis_id_, adjointParams, isFirstStep);
      ptPDE_->GetSolveStep()->PostStepTrans();
      
      if(optimization != NULL){
        optimization->TimeStepCalculated(actTimeStep_, adjointParams);
      }
      
      if(write_results){
        // writing results in output-file(s)
        resHandler->BeginStep( actTimeStep_, steptime );
        ptPDE_->WriteResultsInFile(actTimeStep_, steptime );
        resHandler->FinishStep( );
      }

      // writing current PDE-state into the restart-file
      if (restartIncr_ >= 1){
        if (  actTimeStep_ == nextRestart  ) { 
          std::cout << std::endl << "Write a restart file after step " 
                    << actTimeStep_ <<" *********** " << std::endl;      
          
          ptPDE_->WriteRestart( );
          nextRestart+=restartIncr_ * direction;
        }
      }
      if (haltFlag) {
        std::cout << std::endl << "Write a restart file after interuppting a simulation "
                  << "run with a HALTCFS-file at step:  " 
                  << actTimeStep_ <<" *********** " << std::endl;      

        ptPDE_->WriteRestart( );
      }    

      steptime+=dt*direction;
      isFirstStep = false;
    }
    if(optimization){
      cout << endl;
    }

    // notify resultHandler about finishing of current sequence step
    if(!isPartOfSequence_ && write_results) {
      resHandler->FinishMultiSequenceStep();
      if(optimization == NULL){ // in transient-optimization, finalization may only occur after the last iteration
        resHandler->Finalize();
      }
    }
    
  }
  
  void TransientDriver::ReadRestart() {
    
    if ( progOpts->GetRestart() ){
      ptPDE_->ReadRestart( restartStep_ );
      //startStep = restartStep_ + 1;
      //stepsave=startStep;
      //numstep_+=lastStepToRestartFrom;
      //steptime=(startStep)*dt;
      //restartStep  = startStep + restartIncr_ - 1;
      std::cout << std::endl << "Reading a restart file from step " 
                << restartStep_ <<" ********** " << std::endl;      
    }
    
  }

} // end of namespace
