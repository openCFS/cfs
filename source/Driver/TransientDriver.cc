#include <fstream>
#include <iostream>
#include <string>
#include <stdio.h>
#include <list>
#include <cmath>

// signal handling for catching Ctrl-C
#include <signal.h>

#include "Utils/mathParser/mathParser.hh"
#include "TransientDriver.hh"
#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Utils/Timer.hh"
#include "DataInOut/SimState.hh"
#include "DataInOut/ProgramOptions.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "PDE/StdPDE.hh"
#include "PDE/SinglePDE.hh"
#include "Domain/Domain.hh"
#include "DataInOut/ResultHandler.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

using std::cout;
using std::endl;

namespace CoupledField {

  DEFINE_LOG(trans_driver, "transient_driver")

  // Define pointer to transient driver instance, needed for the signal handler
  // to communicate with
  TransientDriver * instance = NULL;
  
  // ===============
  //   Constructor
  // ===============
  TransientDriver::TransientDriver( UInt sequenceStep,
                                    bool isPartOfSequence,
                                    shared_ptr<SimState> state, Domain* domain,
                                    PtrParamNode paramNode, PtrParamNode infoNode) 
    : SingleDriver( sequenceStep, isPartOfSequence, state, domain, paramNode, infoNode ), 
      timer_(new Timer())
    {
    
    LOG_TRACE(trans_driver) << "TransientDriver():  "
                            << " sequenceStep: " << sequenceStep
                            << " isPartOfSequence: " << isPartOfSequence;
    // Set analysistype
    analysis_ = BasePDE::TRANSIENT;

    actTime_ = 0.0;
    actTimeStep_ = 0;
    initialTime_ = 0.0;
    firstdt_ = 0.0;
    restartStep_ = 0;
    endStep_ = 0;
    isRestarted_ = false;

    // get parameter node
    param_ = param_->Get("transient");

    info_ = info_->Get("transient");
    info_->Get(ParamNode::HEADER)->Get("unit")->SetValue("s");
    
    // for the evaluation of deltaT, we make use of math Parser to
    // allow variable definitions of time step size
    firstdt_ = param_->Get( "deltaT")->MathParse<Double>();
    
    // Get time stepping information from parameter object
    numstep_ = param_->Get( "numSteps")->MathParse<UInt>();
    
    // query if accumulated time should be used as initial time
    std::string initTimeString = param_->Get("initialTime")->As<std::string>();
    useAccumulatedTime_ = initTimeString == "accumulated" ? true : false;
    
    // Check for presence of restart flag.
    // Note; in case of a multisequence analysis we always have to 
    //       write a "restart", as it will be the basis for the
    //       next multisequence step
    writeRestart_ = true;
    PtrParamNode restartNode = param_->Get("writeRestart", ParamNode::PASS);
    if (restartNode && !isPartOfSequence)
      writeRestart_ = restartNode->As<bool>();
  
    // read flag if all results should get written to database file section
    // to allow e.g. for general postprocessing or result extraction
    param_->GetValue("allowPostProc", writeAllSteps_, ParamNode::PASS );
    
    // in the end, directly register the global transient variables
    mathParser_->SetValue( MathParser::GLOB_HANDLER, "step", 1 );
    
    // Accumulated time so far (including time of all previous MS steps)
    accTime_ = 0.0;

    // Estimated time per step
    timePerStep_ = 0.0;

    // register signal handler only, if it is a child driver
    if( !simState_->HasInput() ) {
      if( signal( SIGINT, TransientDriver::SignalHandler) == SIG_ERR ) {
        EXCEPTION( "Could not register Signal Handler");
      }

      // store pointer to global instance variable, if not yet set
      if( !instance ) {
        instance = this;
      }
    }
  }

  void TransientDriver::SetAccumulatedTime(Double accTime ) {
    accTime_ = accTime;
    
    if( !useAccumulatedTime_ ) 
      return;

    initialTime_ = accTime;
    actTime_ = accTime;
    mathParser_->SetValue( MathParser::GLOB_HANDLER, "t", actTimeStep_ );
    mathParser_->SetValue( MathParser::GLOB_HANDLER, "t0", initialTime_ );
    
  }
  // ==============
  //   Destructor
  // ==============
  TransientDriver::~TransientDriver()
  { 
    if( !simState_->HasInput() ) {
      // unregister signal handler and use default action
      // register signal handler
      if( signal( SIGINT, SIG_DFL) == SIG_ERR ) {
        std::cerr << "Could not assign default signal action" << std::endl; // to exceptions is destructors with gcc 6
        domain->GetInfoRoot()->ToFile();
        exit(-1);
      }

      // set global pointer to zero
      instance = NULL;
    }
  }

  // ==================
  //   Initialization
  // ==================
  void TransientDriver::Init( bool restart) 
  {
    isRestarted_ = restart;
    InitializePDEs();
  }
    
  void TransientDriver::SolveProblem()
  {
    // notify resultHandler about beginning of new sequence step 
    ResultHandler * resHandler = domain_->GetResultHandler();

    if(writeRestart_ || writeAllSteps_ || isPartOfSequence_ )
      simState_->BeginMultiSequenceStep(sequenceStep_, analysis_);
    
    // Check for restart
    ReadRestart();
    
    // correct numstep due to restart
    numstep_ = numstep_ - restartStep_;
    
    UInt startStep = restartStep_ + 1;
    endStep_ = numstep_ + restartStep_;
    actTime_  = firstdt_ * startStep + initialTime_;
    Double  dt = firstdt_;
    Double timeStepPercent = (double)numstep_/10;
    Double percentCounter = timeStepPercent;
  
 
   
    ptPDE_->WriteGeneralPDEdefines();
    ptPDE_->GetSolveStep()->SetStartStep( startStep );
    ptPDE_->GetSolveStep()->SetNumTimeSteps(restartStep_+numstep_);
    
    //just initialize some variables
    ptPDE_->GetSolveStep()->InitTimeStepping();
    
    // Ensure, that at least one step has to be computed, otherwise leave
    if( numstep_ == 0 ) {
     return;
    }
    
    resHandler->BeginMultiSequenceStep( sequenceStep_, analysis_, numstep_+restartStep_ );
    
    //only used if AMG is set
    ptPDE_->GetSolveStep()->SetAuxMat(false);

    // Outer loop over all timesteps
    UInt count = 0;
    for (actTimeStep_ = startStep; 
         actTimeStep_ <= endStep_; 
         actTimeStep_ += 1, count++) {

      LOG_DBG(trans_driver) << "loop over timestep " << actTimeStep_;
      
      // start timer only in 2nd loop, as in the 1st step normally the
      // factorization is incorporated,
      if( actTimeStep_ == startStep+1) {
        timer_->Start();
      }
      
      // Set current value of timestep and time step size in the mathParser
      mathParser_->SetValue( MathParser::GLOB_HANDLER, "t", actTime_ );
      mathParser_->SetValue( MathParser::GLOB_HANDLER, "dt", dt );    
      mathParser_->SetValue( MathParser::GLOB_HANDLER, "step", actTimeStep_ );    

      // Determine when to write logging information on terminal
      bool log = false;
      if(numstep_ <= 50)
        log = true;
      if(numstep_ > 50 && numstep_ <= 500 && ! (actTimeStep_ % 10)) // every tenth step, not all but the tenth
        log = true;
      if(numstep_ > 500 && (  (double) count >= percentCounter)   ){
        log = true;
        percentCounter += timeStepPercent;
      }

      /*
       * for debugging -> remove later one
       * TODO: remove!
       */
      if(true){
        log = true;
      }


      if(log)
      {
        if(progOpts->IsQuiet())
          cout << ptPDE_->GetName() << ": Time step " << actTimeStep_ 
               << " time " << actTime_ << endl;
        else
          // write std::out info
          cout << endl << ptPDE_->GetName() << ": Time step "
          << actTimeStep_ <<" ======================= " << endl;
      }
      

      analysis_id_.time = actTime_;
      analysis_id_.step = actTimeStep_;

      // Perform actions
      ptPDE_->GetSolveStep()->SetActTime(actTime_);
      ptPDE_->GetSolveStep()->SetActStep(actTimeStep_);
      ptPDE_->GetSolveStep()->PreStepTrans();
      ptPDE_->GetSolveStep()->SolveStepTrans();
      ptPDE_->GetSolveStep()->PostStepTrans();

      // writing results in output-file(s)
      resHandler->BeginStep( actTimeStep_, actTime_ );
      ptPDE_->WriteResultsInFile(actTimeStep_, actTime_ );
      resHandler->FinishStep( );
      
      // write out re-start only in the last step
      if( actTimeStep_ == endStep_ || abortSimulation_  || writeAllSteps_ ) {
        if( writeRestart_ || writeAllSteps_ || isPartOfSequence_)
         simState_->WriteStep( actTimeStep_, actTime_ );
      }
        
      // leave loop, if simulation should be aborted
      if ( abortSimulation_ ) {
        break;
      }

      // increase current time step
      actTime_+=dt;

      // perform runtime estimation (after 1st step)
      if( actTimeStep_ > 1 ) {
        Double totalTime = timer_->GetWallTime();
        timePerStep_ = totalTime / (Double) count;
        Double remainingTime = (endStep_ - actTimeStep_) * timePerStep_;
        auto time = std::chrono::system_clock::now();
        time += std::chrono::seconds(static_cast<long int>(remainingTime));
        PtrParamNode envNode = info_->GetRoot()->Get(ParamNode::HEADER)->Get("environment");
        envNode->Get("estimatedEnd")->SetValue(Timer::TimeStamp(time));
        envNode->Get("remainingTime")->SetValue(remainingTime);
        envNode->Get("timePerStep")->SetValue(timePerStep_);
      }
      
    }

    // notify resultHandler about finishing of current sequence step
    resHandler->FinishMultiSequenceStep();
    if(!isPartOfSequence_) {
      resHandler->Finalize();
    }
    
    if(writeRestart_ || writeAllSteps_ || isPartOfSequence_ ) {
      simState_->FinishMultiSequenceStep( !abortSimulation_, 
                                          accTime_+GetDuration());
    }

  }

  void TransientDriver::SetToStepValue(UInt stepNum, Double stepVal ) {
    // ensure that this method is only called if simState has input
    if( ! simState_->HasInput()) {
      EXCEPTION( "Can only set external time step, if simulation state "
              << "is read from external file" );
    }
    
    actTime_ = stepVal;
    actTimeStep_ = stepNum;
    mathParser_->SetValue( MathParser::GLOB_HANDLER, "t", actTime_ );
    mathParser_->SetValue( MathParser::GLOB_HANDLER, "step", actTimeStep_ );    
  }
  
  void TransientDriver::ReadRestart() {

    if ( isRestarted_ ) {

      // Ensure simState is present
      assert( simState_ );

      // Create input reader from current output reader
      bool hasInput = simState_->SetInputReaderToSameOutput();
      if( !hasInput)  {
        EXCEPTION( "Can not perform restarted simulation, as HDF5 file "
            << "contains no restart information.");
      }
      
      // Obtain last step
      UInt lastStepNum;
      Double lastStepVal;
      simState_->GetLastStepNum(sequenceStep_, lastStepNum, lastStepVal );
      restartStep_ = lastStepNum;
      
      // if lastStep is 0, no restart possibility
      if( lastStepNum == 0 ) {
        EXCEPTION( "Can not perform restarted simulation, as HDF5 file "
            << "contains no restart information.");
      }

      // Restore coefficients from last restart step
      simState_->UpdateToStep(sequenceStep_, lastStepNum);

      if (simState_->IsCompleted( sequenceStep_ ) && restartStep_ == numstep_) {
        std::cout << "\n\n";
        std::cout << "*******************************************************\n";
        std::cout << " No restart necessary, as the desired number of \n";
        std::cout << " time steps are already computed. \n";
        std::cout << "*******************************************************\n\n";
        return;
      }
      else {
        std::cout << "\n\n";
        std::cout << "*******************************************************\n";
        std::cout << " Continuing simulation from step " << restartStep_  << std::endl;
        std::cout << "*******************************************************\n";
      }
    }

  }
  
  void TransientDriver::SignalHandler( int sig ) {

    if( !instance->abortSimulation_) {
      
      // in addition check, if the current step is the last step, so we
      // do not have to print out any message
      if( instance->actTimeStep_ == instance->endStep_ )
        return;
      
      instance->abortSimulation_ = true;
      
      std::cout << "\n\n";
      std::cout << "*******************************************************\n";
      std::cout << " Simulation will be halted after the current time \n";
      std::cout << " step " << instance->actTimeStep_ << " / " << instance->actTime_
                << " s.\n\n";
      std::cout << " Estimated time before end of simulation run: " << 
          int(instance->timePerStep_) << " s" << std::endl;
      std::cout << "*******************************************************\n\n";

    }
  }

} // end of namespace
