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
  TransientDriver::TransientDriver(Domain * adomain, 
                                   UInt stepOffset,
                                   Double timeOffset, 
                                   std::string driverTag,
                                   bool isPartOfSequence) 
    : SingleDriver(adomain, stepOffset, timeOffset, 
                   driverTag, isPartOfSequence)
  {
    ENTER_FCN( "TransientDriver::TransientDriver" );
  
    // Set analysistype
    analysis_ = TRANSIENT;

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


  // =================
  //   Solve Problem
  // =================
  void TransientDriver::SolveProblem()
  {
    ENTER_FCN( "TransientDriver::SolveProblem" );
  
    Double  steptime  = firstdt_;
    UInt stepsave     = isavebegin_;
    UInt restartStep  = isavebegin_ + restartIncr_ - 1;
  
    Double  dt = firstdt_;
    bool haltFlag=false;
  
    // if driver is not part of multiSequence Driver, get list
    // of pdes which have to be solved and intialize them
    if (isPartOfSequence_ == false) {     
      GetMyPDEs();
      Info->StartProgress ("Starting to solve problem", false);
    }

    Double timeStepPercent = (double)numstep_/10;
    Double percentCounter = timeStepPercent;
  
  
    // if multiSequence is performed, the ms-driver
    // writes out the grid one time
    if (! isPartOfSequence_)
      ptdomain_->PrintGrid();
  
    ptPDE_->WriteGeneralPDEdefines();

    UInt startStep=1, lastStepToRestartFrom=0;
  
    if ( commandLine->GetRestart() ){
      ptPDE_->ReadRestart(lastStepToRestartFrom);
      startStep = lastStepToRestartFrom + 1;
      stepsave=startStep;
      numstep_+=lastStepToRestartFrom;
      ptPDE_->GetSolveStep()->SetStartStep(startStep);
      steptime=(startStep)*dt;
      restartStep  = startStep + restartIncr_ - 1;
      std::cout << myEndl << "Reading a restart file from step " 
                << lastStepToRestartFrom <<" ********** " << std::endl;      
      ptPDE_->GetSolveStep()->SetStartStep(startStep);
    }

    else {
      ptPDE_->GetSolveStep()->SetStartStep(1);
    }
    // Solve problem
    ptPDE_->GetSolveStep()->SetTimeStep(dt);

    ptPDE_->GetSolveStep()->SetNumTimeSteps(numstep_);

    UInt nstep;

    // Out loop over all timesteps
    for (nstep = startStep; nstep <= numstep_; nstep++) {

      // check for a HALTCFS File
      // if there exist a file with name HALTCFS in the executing directory
      // than CFS++ will end after the current time step
      std::ifstream readHALTCFS("HALTCFS", std::ios_base::in );
      if (readHALTCFS) {
        readHALTCFS.close();
        haltFlag = true;
        numstep_=nstep;
        ptPDE_->GetSolveStep()->SetNumTimeSteps(numstep_);
      }
      
      // Set curent value of timestep in the mathParser
      domain->GetMathParser()->SetValue( MathParser::GLOB_HANDLER,
                                         "t", steptime );
    
      // Determine when to write logging information on terminal
      if ( numstep_ <= 50 )
        Info->WriteTimeStep(ptPDE_->GetName(), nstep+stepOffset_, 
                            steptime+timeOffset_);
      else if ( (numstep_ > 50) && (numstep_ <= 500) ) {
        if ( (nstep%10) == 0 )            
          Info->WriteTimeStep(ptPDE_->GetName(), nstep+stepOffset_, 
                              steptime+timeOffset_);
      }
      else if ( numstep_ > 500 ) {
        // Output in steps of ten percent 
        if ((double)nstep >= percentCounter  ){           
          Info->WriteTimeStep(ptPDE_->GetName(), nstep+stepOffset_, 
                              steptime+timeOffset_);
          percentCounter += timeStepPercent;
        }
      }
      
      // Perform actions
      ptPDE_->GetSolveStep()->SetActTime(steptime);
      ptPDE_->GetSolveStep()->SetActStep(nstep);
      ptPDE_->GetSolveStep()->PreStepTrans();
      ptPDE_->GetSolveStep()->SolveStepTrans();
      ptPDE_->GetSolveStep()->PostStepTrans();

      ptPDE_->PostProcess();  
      
      //write history data
      ptPDE_->WriteHistoryInFile(nstep, steptime, stepOffset_, timeOffset_);
    
      // writing results in output-file

      if (nstep == stepsave && (nstep <= isaveend_)) { 
        ptPDE_->WriteResultsInFile(nstep, steptime, stepOffset_, timeOffset_);
        stepsave+=isaveincr_;
      }

      // writing current PDE-state into the restart-file
      if (restartIncr_ >= 1){
        if ( ( nstep == restartStep && (nstep <= isaveend_) ) || (nstep == isaveend_) ) { 
          std::cout << myEndl << "Write a restart file after step " 
                    << nstep <<" *********** " << std::endl;      

          ptPDE_->WriteRestart(nstep);
          restartStep+=restartIncr_;
        }
      }
      if (haltFlag) {
        std::cout << myEndl << "Write a restart file after interuppting a simulation run with a HALTCFS-file at step:  " 
                  << nstep <<" *********** " << std::endl;      

        ptPDE_->WriteRestart(nstep);
      }    

      steptime+=dt;        

      SETPROFILE("After Transient Step");
    }
  
  }

} // end of namespace
