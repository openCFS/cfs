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
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "PDE/StdPDE.hh"
#include "PDE/pdememento.hh"
#include "Domain/domain.hh"



namespace CoupledField {


  // ===============
  //   Constructor
  // ===============
  TransientDriver::TransientDriver( UInt stepOffset,
                                    Double timeOffset, 
                                    std::string driverTag,
                                    bool isPartOfSequence) 
    : SingleDriver( driverTag, isPartOfSequence ) {
    ENTER_FCN( "TransientDriver::TransientDriver" );
    
    // Set analysistype
    analysis_ = TRANSIENT;

    stepOffset_ = stepOffset;
    timeOffset_ = timeOffset;
    actTimeStep_ = 0;
    firstdt_ = 0.0;
    isavebegin_ = 0;
    isaveincr_ = 0;
    isaveend_ = 0;
    restartIncr_ = 0;
    restartStep_ = 0;

    // vecotrs for accessing parameters
    StdVector<std::string> keyVec, attrVec, valVec;

    attrVec = "tag";
    valVec = driverTag_;
   
    // Get time stepping information from parameter object
    keyVec = "transient", "numSteps";
    params->Get(keyVec, attrVec, valVec, numstep_);
  
    keyVec = "transient", "firstDt";
    params->Get(keyVec, attrVec, valVec, firstdt_);
  
    keyVec = "transient", "stepSaveBeg";
    params->Get(keyVec, attrVec, valVec, isavebegin_);

    keyVec = "transient", "stepSaveEnd";
    params->Get(keyVec, attrVec, valVec, isaveend_);

    keyVec = "transient", "stepSaveInc";
    params->Get(keyVec, attrVec, valVec, isaveincr_);

    keyVec = "transient", "writeRestartInc";
    params->Get(keyVec, attrVec, valVec, restartIncr_);

    // Make consistency check. In fact in the XML case the Schema should catch
    // this error. But one can never be sure.
    if(isavebegin_ <= 0)  {
      Error( "Value of stepsavebegin must be positive and nonzero!",
             __FILE__, __LINE__ );
    }
  
    ptMeshes_ = NULL;
  

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
  
    UInt startStep = restartStep_ + 1;

    Double  steptime  = startStep * firstdt_;
    UInt stepsave     = startStep > isavebegin_ ? startStep : isavebegin_;
    isaveend_ += restartStep_;
    UInt nextRestart  = isavebegin_ + restartStep_ + restartIncr_ - 1;
    Double  dt = firstdt_;
    bool haltFlag=false;
 
    Double timeStepPercent = (double)numstep_/10;
    Double percentCounter = timeStepPercent;
  
  
    // if multiSequence is performed, the ms-driver
    // writes out the grid one time
    if (! isPartOfSequence_)
      domain->PrintGrid();
  
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
      
      // Set curent value of timestep in the mathParser
      domain->GetMathParser()->SetValue( MathParser::GLOB_HANDLER,
                                         "t", steptime );
    
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

      ptPDE_->PostProcess();  
      
      //write history data
      ptPDE_->WriteHistoryInFile(actTimeStep_, steptime, 
                                 stepOffset_, timeOffset_);
    
      // writing results in output-file

      if (actTimeStep_ == stepsave && (actTimeStep_ <= isaveend_)) { 
        ptPDE_->WriteResultsInFile(actTimeStep_, steptime, 
                                   stepOffset_, timeOffset_);
        stepsave+=isaveincr_;
      }

      // writing current PDE-state into the restart-file
      if (restartIncr_ >= 1){
        if ( ( actTimeStep_ == nextRestart && (actTimeStep_ <= isaveend_) ) 
             || (actTimeStep_ == isaveend_) ) { 
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
