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
#include "Domain/domain.hh"
#include "DataInOut/resultHandler.hh"
#include "DataInOut/Logging/cfslog.hh"

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
    ParamNode * myNode = 
      param->Get("sequenceStep", "index", GenStr( sequenceStep ) )
      ->Get("analysis")->Get("transient");
    
    // Get time stepping information from parameter object
    myNode->Get( "numSteps", numstep_ );
    myNode->Get( "deltaT", firstdt_ );

    // Get save increment for restart file (optional)
    myNode->Get( "writeRestartInc", restartIncr_, false );
  
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
    
  void TransientDriver::SolveProblem(bool write_results, InfoNode* given_analysis_id) 
  {
    // options not implemented
    assert(write_results == true);

    // notify resultHandler about beginning of new sequence step 
    ResultHandler * resHandler = domain->GetResultHandler();

    UInt startStep = restartStep_ + 1;
    Double  steptime  = startStep * firstdt_;
    UInt nextRestart  = restartStep_ + restartIncr_ ;
    Double  dt = firstdt_;
    bool haltFlag=false;

    Double timeStepPercent = (double)numstep_/10;
    Double percentCounter = timeStepPercent;
  
    resHandler->BeginMultiSequenceStep( sequenceStep_,
                                        analysis_,
                                        numstep_+restartStep_-startStep+1 );
  
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
    for (actTimeStep_ = startStep; actTimeStep_ <= numstep_+restartStep_; actTimeStep_++) {
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
      if(numstep_ <= 50)
        log = true;
      if(numstep_ > 50 && numstep_ <= 500 && actTimeStep_ % 10) 
        log = true;
      if(numstep_ > 500 && (double) actTimeStep_ >= percentCounter) {
        log = true;
        percentCounter += timeStepPercent;
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
      
      if(given_analysis_id == NULL)
      {
        analysis_id_ = info->Get("analysis")->Get(InfoNode::PROCESS)->Get("step", InfoNode::APPEND); 
        analysis_id_->Get("analysis_id")->SetValue(actTimeStep_);
      }
      else
      {
        analysis_id_ = given_analysis_id;
        assert(analysis_id_->Has("analysis_id"));
      }
      analysis_id_->Get("step")->SetValue(actTimeStep_);
      analysis_id_->Get("value")->SetValue(steptime);
      
      // Perform actions
      ptPDE_->GetSolveStep()->SetActTime(steptime);
      ptPDE_->GetSolveStep()->SetActStep(actTimeStep_);
      ptPDE_->GetSolveStep()->PreStepTrans();
      ptPDE_->GetSolveStep()->SolveStepTrans(analysis_id_);
      ptPDE_->GetSolveStep()->PostStepTrans();
      // writing results in output-file(s)
      resHandler->BeginStep( actTimeStep_, steptime );
      ptPDE_->WriteResultsInFile(actTimeStep_, steptime );
      resHandler->FinishStep( );

      // writing current PDE-state into the restart-file
      if (restartIncr_ >= 1){
        if (  actTimeStep_ == nextRestart  ) { 
          std::cout << std::endl << "Write a restart file after step " 
                    << actTimeStep_ <<" *********** " << std::endl;      
          
          ptPDE_->WriteRestart( );
          nextRestart+=restartIncr_;
        }
      }
      if (haltFlag) {
        std::cout << std::endl << "Write a restart file after interuppting a simulation "
                  << "run with a HALTCFS-file at step:  " 
                  << actTimeStep_ <<" *********** " << std::endl;      

        ptPDE_->WriteRestart( );
      }    

      steptime+=dt;        

      SETPROFILE("After Transient Step");
    }

    // notify resultHandler about finishing of current sequence step
    if(!isPartOfSequence_ ) {
      resHandler->FinishMultiSequenceStep();
      resHandler->Finalize();
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
