#include <fstream>
#include <iostream>
#include <string>
#include <stdio.h>
#include <list>
#include <math.h>

// signal handling for post-poning Ctr-C
#include <signal.h>

#include <boost/filesystem.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

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
namespace fs = boost::filesystem;
namespace pt = boost::posix_time;

namespace CoupledField {

  DECLARE_LOG(trans_driver)
  DEFINE_LOG(trans_driver, "transient_driver")

  // Define pointer to transient driver instance, needed for the signal handler
  // to communicate with
  TransientDriver * instance = NULL;
  
  // ===============
  //   Constructor
  // ===============
  TransientDriver::TransientDriver( UInt sequenceStep,
                                    bool isPartOfSequence,
                                    shared_ptr<SimState> state, Domain* domain ) 
    : SingleDriver( sequenceStep, isPartOfSequence, state, domain ), 
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
    abortSimulation_ = false;

    // get parameter node
    PtrParamNode myNode = 
      domain_->GetParamRoot()->GetByVal("sequenceStep", std::string("index"), sequenceStep)
      ->Get("analysis")->Get("transient");

    driverNode = driverNode->Get("transient");
    driverNode->Get("sequenceStep")->SetValue(sequenceStep);
    driverNode->Get(ParamNode::PN_HEADER)->Get("unit")->SetValue("s");
    
    // for the evaluation of deltaT, we make use of math Parser to
    // allow variable definitions of time step size
    firstdt_ = myNode->Get( "deltaT")->MathParse<Double>();
    
    // Get time stepping information from parameter object
    numstep_ = myNode->Get( "numSteps")->MathParse<UInt>();
    
    // query if accumulated time should be used as initial time
    std::string initTimeString = myNode->Get("initialTime")->As<std::string>();
    useAccumulatedTime_ = initTimeString == "accumulated" ? true : false;
    
    // Check for presence of restart flag
    writeRestart_ = true;
    PtrParamNode restartNode = myNode->Get("writeRestart", ParamNode::PASS);
    if (restartNode)
      writeRestart_ = restartNode->As<bool>();
  
    // in the end, directly register the global transient variables
    domain_->GetMathParser()->SetValue( MathParser::GLOB_HANDLER,
                                       "t", 0 );
    domain_->GetMathParser()->SetValue( MathParser::GLOB_HANDLER,
                                        "t0", 0 );
    domain_->GetMathParser()->SetValue( MathParser::GLOB_HANDLER,
                                       "dt", 0.0 );    
    domain_->GetMathParser()->SetValue( MathParser::GLOB_HANDLER,
                                       "step", 1 );
    
    // register signal handler
    if( signal( SIGINT, TransientDriver::SignalHandler) == SIG_ERR ) {
      EXCEPTION( "Could not register Signal Handler");
    }
    
    // store pointer to global instance variable, if not yet set
    if( !instance ) {
      instance = this;
    }
  }

  void TransientDriver::SetAccumulatedTime(Double accTime ) {
    if( !useAccumulatedTime_ ) 
      return;

    initialTime_ = accTime;
    actTime_ = accTime;
    domain_->GetMathParser()->SetValue( MathParser::GLOB_HANDLER,
                                        "t", actTimeStep_  );
    domain_->GetMathParser()->SetValue( MathParser::GLOB_HANDLER,
                                        "t0", initialTime_ );
    
  }
  // ==============
  //   Destructor
  // ==============
  TransientDriver::~TransientDriver()
  { 
    // unregister signal handler and use default action
    // register signal handler
    if( signal( SIGINT, SIG_DFL) == SIG_ERR ) {
      EXCEPTION( "Could not assign default signal action");
    }
    
    // set global pointer to zero
    instance = NULL;
    
  }

  // ==================
  //   Initialization
  // ==================
  void TransientDriver::Init() 
  {
    InitializePDEs();
  }
    
  void TransientDriver::SolveProblem(bool write_results, PtrParamNode given_analysis_id) 
  {
    // notify resultHandler about beginning of new sequence step 
    ResultHandler * resHandler = domain_->GetResultHandler();

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
    
    if(write_results){
       resHandler->BeginMultiSequenceStep( sequenceStep_, analysis_, numstep_+restartStep_ );
     }
    
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
      domain_->GetMathParser()->SetValue( MathParser::GLOB_HANDLER,
                                         "t", actTime_ );
      domain_->GetMathParser()->SetValue( MathParser::GLOB_HANDLER,
                                         "dt", dt );    
      domain_->GetMathParser()->SetValue( MathParser::GLOB_HANDLER,
                                         "step", actTimeStep_ );    

      // Determine when to write logging information on terminal
      bool log = false;
      if(numstep_ <= 50)
        log = true;
      if(numstep_ > 50 && numstep_ <= 500 && ! (actTimeStep_ % 10)) // every tenth step, not all but the tenth
        log = true;
      if(numstep_ > 500 && (  (double) actTimeStep_ >= percentCounter)   ){
        log = true;
        percentCounter += timeStepPercent;
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
      

      if(given_analysis_id == NULL)
      {
        // do we really want to create a new entry? Might blast up the output
        ParamNode::ActionType at = progOpts->DoDetailedInfo() ? ParamNode::APPEND : ParamNode::DEFAULT;
        analysis_id_ = driverNode->Get(ParamNode::PN_PROCESS)->Get("step", at);
        analysis_id_->Get("analysis_id")->SetValue(actTimeStep_);
      }
      else
      {
        assert(given_analysis_id->Has("analysis_id"));
        analysis_id_ = CreateAnalysisIdChild(given_analysis_id, given_analysis_id->Get("analysis_id")->As<std::string>(), actTimeStep_, "step", actTimeStep_);
      }
      analysis_id_->Get("step")->SetValue(actTimeStep_);
      analysis_id_->Get("value")->SetValue(actTime_);
      
      // Perform actions
      ptPDE_->GetSolveStep()->SetActTime(actTime_);
      ptPDE_->GetSolveStep()->SetActStep(actTimeStep_);
      ptPDE_->GetSolveStep()->PreStepTrans();
      ptPDE_->GetSolveStep()->SolveStepTrans(analysis_id_);
      ptPDE_->GetSolveStep()->PostStepTrans();

      if(write_results){
        // writing results in output-file(s)
        resHandler->BeginStep( actTimeStep_, actTime_ );
        ptPDE_->WriteResultsInFile(actTimeStep_, actTime_ );
        resHandler->FinishStep( );
      }
      
      // write out re-start only in the last step
      if( actTimeStep_ == endStep_ || abortSimulation_ ) {
        if( writeRestart_)
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
        pt::ptime now = pt::second_clock::local_time();
        now += pt::seconds(static_cast<long int>(remainingTime));
        
        analysis_id_->Get("timePerStep")->SetValue( timePerStep_ );
        PtrParamNode envNode = driverNode->GetRoot()->
            Get(ParamNode::PN_HEADER)->Get("environment");
        envNode->Get("estimatedEnd")->SetValue(pt::to_simple_string( now ));
        envNode->Get("remainingTime")->SetValue(remainingTime);
      }
      
    }

    // notify resultHandler about finishing of current sequence step
    if(!isPartOfSequence_ && write_results) {
      resHandler->FinishMultiSequenceStep();
      // notify resultHandler about finishing of current sequence step
      if(!isPartOfSequence_) resHandler->Finalize();
    }
    simState_->FinishMultiSequenceStep();
    
  }

  void TransientDriver::ReadRestart() {

    if ( progOpts->GetRestart() ){

      // Ensure simState is present
      assert( simState_ );

      // Create input reader from current output reader 
      simState_->SetInputReaderToSameInput();

      // Obtain last step
      UInt lastStep = simState_->GetLastStepNum();
      restartStep_ = lastStep;
      
      // if lastStep is 0, no restart possibility
      if( lastStep == 0 ) {
        EXCEPTION( "Can not perform restarted simulation, as HDF5 file "
            << "contains no restart information.");
      }

      // Store restart step
      simState_->UpdateToStep(lastStep);

      if( lastStep == numstep_) {

        std::cerr << "\n\n";
        std::cerr << "*******************************************************\n";
        std::cerr << " No restart necessary, as the desired number of \n";
        std::cerr << " time steps are already computed. \n";
        std::cerr << "*******************************************************\n\n";
        return;
      } else{

        std::cerr << "\n\n";
        std::cerr << "*******************************************************\n";
        std::cerr << " Continuing simulation from step " << restartStep_  << std::endl;
        std::cerr << "*******************************************************\n";
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
      
      std::cerr << "\n\n";
      std::cerr << "*******************************************************\n";
      std::cerr << " Simulation will be halted after the current time \n";
      std::cerr << " step " << instance->actTimeStep_ << " / " << instance->actTime_
                << " s.\n\n";
      std::cerr << " Estimated time before end of simulation run: " << 
          int(instance->timePerStep_) << " s" << std::endl;
      std::cerr << "*******************************************************\n\n";

    }
  }

} // end of namespace
