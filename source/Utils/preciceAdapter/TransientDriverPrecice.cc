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

    std::cout << "Using TransientDriverPrecice" << std::endl;
    if(!preciceAdapter_)
      EXCEPTION("No precice adapter set in domain!")
      
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

    simState_->BeginMultiSequenceStep(sequenceStep_, analysis_);
    
    actTimeStep_ = 1;
    UInt startStep = actTimeStep_;
    endStep_ = numstep_ ;
    actTime_  = firstdt_ * startStep + initialTime_;
    Double  cfs_dt = firstdt_;

    std::cout << "actTimeStep_: " << actTimeStep_ << ", actTime_: " << actTime_ << ", firstdt_: " << firstdt_ << ", initialTime_: " << initialTime_ << std::endl;
    std::cout << "endStep_: " << endStep_ << std::endl;
  
    UInt checkpointStep = 0;
    Double checkpointTime = 0.0;

    // init precice coupling adapter...already done in SinglePDE.cc domain_->InitPreciceAdapter(this)
    // when using precice, we need direct access to the solvestep
    preciceAdapter_->RegisterSolveStep(ptPDE_->GetSolveStep());

   
    ptPDE_->WriteGeneralPDEdefines();
    ptPDE_->GetSolveStep()->SetStartStep( startStep );
    ptPDE_->GetSolveStep()->SetNumTimeSteps(restartStep_+numstep_);
    
    //just initialize some variables
    ptPDE_->GetSolveStep()->InitTimeStepping();
    resHandler->BeginMultiSequenceStep( sequenceStep_, analysis_, numstep_+restartStep_ );

    // following the description of https://precice.org/couple-your-code-implicit-coupling.html
    // but we do not explicitely save the state because we store all steps and overwrite
    // actTimeStep_ is zero at this point (we do not consider restart when using precice ... currently. But it 
    // is not a big deal, just handling the step values correctly, the functionality is there)
    while(preciceAdapter_->IsCouplingOngoing()){

      if(preciceAdapter_->RequiresWritingCheckpoint()){
        LOG_DBG(trans_driverprecice) << "checkpointing timestep " << actTimeStep_;
        checkpointStep = actTimeStep_;
        checkpointTime = actTime_;
      }


      LOG_DBG(trans_driverprecice) << "loop over timestep " << actTimeStep_;
      
      // start timer only in 2nd loop, as in the 1st step normally the
      // factorization is incorporated,
      if( actTimeStep_ == startStep+1) {
        timer_->Start();
      }
      // timer_->Start();

      cfs_dt = this->GetDeltaT();
      double precice_dt =  preciceAdapter_->GetMaxTimeStepSize();
      double dt = (cfs_dt <= precice_dt) ? cfs_dt : precice_dt;

      std::cout << "cfs_dt: " << cfs_dt << ", precice_dt: " << precice_dt << ", using dt: " << dt << std::endl;



    // Set current value of timestep and time step size in the mathParser
      mathParser_->SetValue( MathParser::GLOB_HANDLER, "t", actTime_ );
      mathParser_->SetValue( MathParser::GLOB_HANDLER, "dt", cfs_dt );    
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
      // Mark PreCICE read results as updated so ResultHandler writes them
      // even though they have no functor (data was filled by RegisterTimeStepReadData)
      preciceAdapter_->MarkReadResultsUpdated();
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

      preciceAdapter_->Advance(dt);


      if(preciceAdapter_->RequiresReadingCheckpoint()){ // iteration did not converge
        // ==========================================
        // back rollback to checkpoint state
        this->CheckpointingToTimestep(checkpointStep);
        actTimeStep_ = checkpointStep;
        actTime_ = checkpointTime;

        // when using precice, we need direct access to the solvestep
        preciceAdapter_->RegisterSolveStep(ptPDE_->GetSolveStep());

        //ptPDE_->WriteGeneralPDEdefines();
        ptPDE_->GetSolveStep()->SetStartStep( checkpointStep );
        ptPDE_->GetSolveStep()->SetNumTimeSteps(numstep_);
        
        //just initialize some variables
        ptPDE_->GetSolveStep()->InitTimeStepping();
        // ==========================================
      }
      else{ // iteration converged
        actTimeStep_ += 1;
        actTime_ += cfs_dt;
      }


    }

    // precice finalize() is handled in top CFS.cc

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
