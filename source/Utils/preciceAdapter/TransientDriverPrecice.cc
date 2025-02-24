#include <fstream>
#include <iostream>
#include <string>
#include <stdio.h>
#include <list>
#include <cmath>

// signal handling for catching Ctrl-C
#include <signal.h>

#include <boost/date_time/posix_time/posix_time.hpp>
#include "Utils/mathParser/mathParser.hh"
#include "TransientDriverPrecice.hh"
#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Utils/Timer.hh"
#include "IPreciceAdapter.hh"
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
namespace pt = boost::posix_time;

namespace CoupledField {

  DEFINE_LOG(trans_driverprecice, "transient_driverprecice")


  
  // ===============
  //   Constructor
  // ===============
  TransientDriverPrecice::TransientDriverPrecice( UInt sequenceStep,
                                    bool isPartOfSequence,
                                    shared_ptr<SimState> state, Domain* domain,
                                    PtrParamNode paramNode, PtrParamNode infoNode) 
    : SingleDriver(sequenceStep, isPartOfSequence, state, domain, paramNode, infoNode), // initialize virtual base
      TransientDriver(sequenceStep, isPartOfSequence, state, domain, paramNode, infoNode)
    {
    TransientDriver::instance = this;  

    LOG_TRACE(trans_driverprecice) << "TransientDriverPrecice():  "
                            << " sequenceStep: " << sequenceStep
                            << " isPartOfSequence: " << isPartOfSequence;

    preciceAdapter_ = domain_->GetPreciceAdapter();
  
  }

  // ==============
  //   Destructor
  // ==============
  TransientDriverPrecice::~TransientDriverPrecice()
  { 
  }


  void TransientDriverPrecice::CheckpointingToTimestep(UInt stepNum) {
      // Ensure simState is present
      assert( simState_ );

      // Create input reader from current output reader
      bool hasInput = simState_->SetInputReaderToSameOutput();
      if( !hasInput)  {
        EXCEPTION( "Can not perform restarted simulation, as HDF5 file "
            << "contains no restart information.");
      }
      simState_->UpdateToStep(sequenceStep_, stepNum);

  }
    
  void TransientDriverPrecice::SolveProblem()
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
  
    // init precice coupling adapter
    //preciceAdapter_->initialize(domain_);
    // when using precice, we need direct access to the solvestep
    preciceAdapter_->RegisterSolveStep(ptPDE_->GetSolveStep());

   
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
    


    // Outer loop over all timesteps
    UInt count = 0;
    for (actTimeStep_ = startStep; 
         actTimeStep_ <= endStep_; 
         actTimeStep_ += 1) {

      LOG_DBG(trans_driverprecice) << "loop over timestep " << actTimeStep_;
      
      // start timer only in 2nd loop, as in the 1st step normally the
      // factorization is incorporated,
      if( actTimeStep_ == startStep+1) {
        timer_->Start();
      }
      
      if((actTimeStep_ == 8) && (count == 0)){
        count=1;
        UInt checkstep = 2;
        this->CheckpointingToTimestep(checkstep);
        //actTimeStep_ = checkstep;
      
        // UInt startStep = checkstep + 1;
        // endStep_ = numstep_;
        // actTime_  = firstdt_ * startStep + initialTime_;
        // Double  dt = firstdt_;
        // Double timeStepPercent = (double)numstep_/10;
        // Double percentCounter = timeStepPercent;
      
        // init precice coupling adapter
        //preciceAdapter_->initialize(domain_);
        // when using precice, we need direct access to the solvestep
        preciceAdapter_->RegisterSolveStep(ptPDE_->GetSolveStep());

        //ptPDE_->WriteGeneralPDEdefines();
        ptPDE_->GetSolveStep()->SetStartStep( startStep );
        ptPDE_->GetSolveStep()->SetNumTimeSteps(numstep_);
        
        //just initialize some variables
        ptPDE_->GetSolveStep()->InitTimeStepping();
   

      }

      // Set current value of timestep and time step size in the mathParser
      mathParser_->SetValue( MathParser::GLOB_HANDLER, "t", actTime_ );
      mathParser_->SetValue( MathParser::GLOB_HANDLER, "dt", dt );    
      mathParser_->SetValue( MathParser::GLOB_HANDLER, "step", actTimeStep_ );    
      

      analysis_id_.time = actTime_;
      analysis_id_.step = actTimeStep_;

      // only if precice is activated and used
      preciceAdapter_->RegisterTimeStepReadData();

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
      
      // only if precice is activated and used
      preciceAdapter_->RegisterTimeStepWriteData();

      // write out re-start only in the last step
      simState_->WriteStep( actTimeStep_, actTime_ );
        
      // leave loop, if simulation should be aborted
      if ( abortSimulation_ ) {
        break;
      }

      // increase current time step
      actTime_+=dt;
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


} // end of namespace
