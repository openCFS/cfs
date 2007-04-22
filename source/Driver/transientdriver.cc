// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <string>
#include <stdio.h>
#include <list>
#include <math.h>


// include generic file system handling routines
// #include "boost/filesystem/operations.hpp"
// #include "boost/filesystem/path.hpp"
// using boost::filesystem;

#include "transientdriver.hh"
#include "stdSolveStep.hh"

#include "DataInOut/CommandLine/BaseCommandLineHandler.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "PDE/StdPDE.hh"
#include "PDE/pdememento.hh"
#include "Domain/domain.hh"
#include "DataInOut/resultHandler.hh"


namespace CoupledField {


  // ===============
  //   Constructor
  // ===============
  TransientDriver::TransientDriver( UInt stepOffset,
                                    Double timeOffset, 
                                    UInt sequenceStep,
                                    bool isPartOfSequence) 
    : SingleDriver( sequenceStep, isPartOfSequence ) {
    ENTER_FCN( "TransientDriver::TransientDriver" );
    
    // Set analysistype
    analysis_ = TRANSIENT;

    stepOffset_ = stepOffset;
    timeOffset_ = timeOffset;
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
    std::string command = "rm -f ./HALTCFS";
    std::system( command.c_str() );
    
//     boost::filesystem::path haltFile ( "./HALTCFS" );
//     if ( boost::filesystem::exists( haltFile ) ) {
//       boost::filesystem::remove( haltFile );
//     }
  }

  // ==============
  //   Destructor
  // ==============
  TransientDriver::~TransientDriver()
  {
    ENTER_FCN( "TransientDriver::~TransientDriver" );
  }

  // ==================
  //   Initialization
  // ==================
  void TransientDriver::Init() {
    ENTER_FCN( "TransientDriver::Init" );


    InitializePDEs();

  }
    

  // =================
  //   Solve Problem
  // =================
  void TransientDriver::SolveProblem()
  {
    ENTER_FCN( "TransientDriver::SolveProblem" );

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

    // Outer loop over all timesteps
    for (actTimeStep_ = startStep; actTimeStep_ <= numstep_+restartStep_; actTimeStep_++) {

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

      // Determine when to write logging information on terminal
      if ( numstep_ <= 50 )
        Info->WriteTimeStep(ptPDE_->GetName(), actTimeStep_+stepOffset_, 
                            steptime+timeOffset_);
      else if ( (numstep_ > 50) && (numstep_ <= 500) ) {
        if ( (actTimeStep_%10) == 0 )            
          Info->WriteTimeStep(ptPDE_->GetName(), actTimeStep_+stepOffset_, 
                              steptime+timeOffset_);
      }
      else if ( numstep_ > 500 ) {
        // Output in steps of ten percent 
        if ((double)actTimeStep_ >= percentCounter  ){           
          Info->WriteTimeStep(ptPDE_->GetName(), 
                              actTimeStep_+stepOffset_, 
                              steptime+timeOffset_);
          percentCounter += timeStepPercent;
        }
      }
      
      // Perform actions
      ptPDE_->GetSolveStep()->SetActTime(steptime);
      ptPDE_->GetSolveStep()->SetActStep(actTimeStep_);
      ptPDE_->GetSolveStep()->PreStepTrans();
      ptPDE_->GetSolveStep()->SolveStepTrans();
      ptPDE_->GetSolveStep()->PostStepTrans();

      
      // writing results in output-file(s)
      resHandler->BeginStep( stepOffset_ + actTimeStep_,
                              timeOffset_ + steptime );
      ptPDE_->WriteResultsInFile(actTimeStep_, steptime, 
                                 stepOffset_, timeOffset_);
      resHandler->FinishStep( );

      // writing current PDE-state into the restart-file
      if (restartIncr_ >= 1){
        if (  actTimeStep_ == nextRestart  ) { 
          std::cout << myEndl << "Write a restart file after step " 
                    << actTimeStep_ <<" *********** " << std::endl;      
          
          ptPDE_->WriteRestart( );
          nextRestart+=restartIncr_;
        }
      }
      if (haltFlag) {
        std::cout << myEndl << "Write a restart file after interuppting a simulation "
                  << "run with a HALTCFS-file at step:  " 
                  << actTimeStep_ <<" *********** " << std::endl;      

        ptPDE_->WriteRestart( );
      }    

      steptime+=dt;        

      SETPROFILE("After Transient Step");
    }

    ptPDE_->Finalize();

    // notify resultHandler about finishing of current sequence step
    if( !isPartOfSequence_ )
      resHandler->FinishMultiSequenceStep();
    
  }
  
  void TransientDriver::ReadRestart() {
    ENTER_FCN( "TransientDriver::ReadRestart" );
    
    if ( commandLine->GetRestart() ){
      ptPDE_->ReadRestart( restartStep_ );
      //startStep = restartStep_ + 1;
      //stepsave=startStep;
      //numstep_+=lastStepToRestartFrom;
      //steptime=(startStep)*dt;
      //restartStep  = startStep + restartIncr_ - 1;
      std::cout << myEndl << "Reading a restart file from step " 
                << restartStep_ <<" ********** " << std::endl;      
    }
    
  }

} // end of namespace
