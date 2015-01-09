#include "MultiSequenceDriver.hh"

#include <set>

// include all kinds of single drivers
#include "StaticDriver.hh"
#include "TransientDriver.hh"
#include "HarmonicDriver.hh"
#include "EigenFrequencyDriver.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "PDE/SinglePDE.hh"
#include "Domain/Domain.hh"
#include "DataInOut/ResultHandler.hh"

#include "DataInOut/Logging/LogConfigurator.hh"

namespace CoupledField {

DECLARE_LOG(msDriver)
DEFINE_LOG(msDriver, "msDriver")

  // ***************
  //   Constructor
  // ***************
  MultiSequenceDriver::MultiSequenceDriver( shared_ptr<SimState> state,
                                            Domain* domain,
                                            PtrParamNode paramNode, PtrParamNode infoNode
                                            ) 
    : BaseDriver(state, domain, paramNode, infoNode) {

    analysis_ = BasePDE::MULTI_SEQUENCE;
    
    sequenceStep_ = 1;
    numSteps_ = 0;
    accumulatedTime_ = 0.0;
    actDriver_ = NULL;
    isRestarted_ = false;
  }


  // **********************
  //   Default destructor
  // **********************
  MultiSequenceDriver::~MultiSequenceDriver() {
    
    // It can happen, that the driver was not delete (e.g. is SolveProblem
    // was not called)
    if( actDriver_ ) {
      delete actDriver_;
      actDriver_ = NULL;
    }
  }


  void MultiSequenceDriver::SolveProblem(bool write_results) {

    std::cout << "++ Starting to solve problem" << std::endl;

    // get resultHandler
    //ResultHandler * resHandler = domain_->GetResultHandler();

    // initiate restart step, which basically sets the sequenceStep_
    // to the step currently to be solved.
    ReadRestart();
    
    // outer loop over all single sequences
    for (UInt iStep = sequenceStep_; iStep <= numSteps_; iStep++ ) {

      // log info about new step to info-class
      WriteMultiSequenceStep(iStep, analysisPerStep_[iStep]);
      
      // Setup step (create driver, create PDEs etc.)
      SetupStep( iStep );
      
      // remember current sequence step
      sequenceStep_ = iStep;
      
      // Solve Problem
      actDriver_->SolveProblem();

      // add up accumulated time
      if (analysisPerStep_[iStep] == BasePDE::TRANSIENT) {
        TransientDriver * tD = dynamic_cast<TransientDriver*>(actDriver_);
        accumulatedTime_ += tD->GetDuration();
      }
      
//      // finish sequence step with resulthandler
//      resHandler->FinishMultiSequenceStep();
//      simState_->FinishMultiSequenceStep( !abortSimulation_, accumulatedTime_ );
      
      // Get solution for next step and delete all PDEs
      if (iStep < numSteps_) {
        // delete PDEs
        domain_->ResetPDEs();
      }

      // delete analysis types
      delete actDriver_;
      actDriver_ = NULL;
      
      // after the first run, we do not have a "restarted" simulation anymore
      isRestarted_ = false;
      
      // Leave, if a supordinate driver has set the abortSimulation_ flag
      if( abortSimulation_ )
        break;
    } // iStep

    // trigger finalization
    domain_->GetResultHandler()->Finalize();
  }


  UInt MultiSequenceDriver::GetActStep( const std::string& pdename ) {
    assert( actDriver_);
    return actDriver_->GetActStep( pdename );
  }

  // *****************
  //   Initializer
  // *****************
  void MultiSequenceDriver::Init(bool restart) 
  {
    LOG_DBG(msDriver) << "Initializing MultiSequenceDriver (restart: " 
        << (restart ? "true" : "false");
    
    isRestarted_ = restart;

    // get nodes for all sequencesteps
    ParamNodeList seqNodes = domain_->GetParamRoot()->GetList("sequenceStep");

    // 1.) Fill vector with step indices and ensure that all occur
    std::set<UInt> stepIndices;
    for( UInt i = 0; i < seqNodes.GetSize(); i++ ) {
      UInt actStepIndex = seqNodes[i]->Get("index")->As<UInt>();
      if( stepIndices.find( actStepIndex ) != stepIndices.end() ) {
        EXCEPTION( "Multisequence step with index " << actStepIndex
                   << " occurs more than one time!") ;
      } else {
        stepIndices.insert( actStepIndex );
      }
    }
    for( UInt i = 1; i <= seqNodes.GetSize(); i++ ) {
      if( stepIndices.find( i ) == stepIndices.end() ) {
        EXCEPTION( "Multisequence step  with index "
            << i << " is not defined!" );
      }
    }

    numSteps_ = stepIndices.size();


    // 3.) Read in all pdes and analysis types

    // iterate over all sequence Steps
    for( UInt iStep = 1; iStep <= numSteps_; iStep++) {

      // get current step node
      PtrParamNode actStepNode = 
          domain_->GetParamRoot()->GetByVal("sequenceStep", std::string("index"), iStep);

      // get current analysistype
      std::string analysisString;
      analysisString = 
          actStepNode->Get("analysis")->GetChild()->GetName();
      analysisPerStep_[iStep] = BasePDE::analysisType.Parse(analysisString);

      // get all pde-nodes in current sequence step
      ParamNodeList pdeNodes;
      pdeNodes = actStepNode->Get("pdeList")->GetChildren();

      // iterate over all pdes
      pdesPerStep_[iStep].Resize( pdeNodes.GetSize() );
      for( UInt iPde = 0; iPde < pdeNodes.GetSize(); iPde++ ) {

        // get pdeName
        std::string pdeName = pdeNodes[iPde]->GetName();
        pdesPerStep_[iStep][iPde] = pdeName;

      } // loop over pdes
    } // loop over sequence steps

    // 4.) Perform final consistency checks
    // Not much to do here yet ...
    
  }
  
  void MultiSequenceDriver::SetSequenceStep(UInt sequenceStep ) {
    LOG_TRACE(msDriver) << "Calling SetSequenceStep with step " << sequenceStep;
    if( sequenceStep > numSteps_) {
      EXCEPTION( "Can not update to mulstistep " << sequenceStep 
                 << " as only " << numSteps_ << " are defined" );
    }
    SetupStep(sequenceStep);
    sequenceStep_ = sequenceStep;
    
  }
  
  void MultiSequenceDriver::SetupStep(UInt sequenceStep) {

    LOG_TRACE(msDriver) << "SetupStep for MS " << sequenceStep;
    // Only setup step once
    if( sequenceStep_ == sequenceStep && 
        actDriver_ !=  NULL) {
      return;
    }
    
    if( actDriver_ ) {
      delete actDriver_;
      actDriver_ = NULL;
    }
    sequenceStep_ = sequenceStep;
    
    // obtain paramNode and infoNode
    PtrParamNode seqNode = param_->GetByVal("sequenceStep","index",sequenceStep)
        ->Get("analysis");
    PtrParamNode info = info_->Get("sequenceStep",ParamNode::APPEND);
    info->SetComment("=== SEQUENCE STEP ===");
    info->Get("index")->SetValue(sequenceStep);
    
    
    // Since per time step only one type of 
    // analysis is allowed, we simple access
    // the first entry fo analysisPerStep_
    if (analysisPerStep_[sequenceStep] == BasePDE::STATIC) {
      actDriver_ = new StaticDriver( sequenceStep, true, 
                                     simState_, domain_, seqNode, info );
    }
    else if (analysisPerStep_[sequenceStep] == BasePDE::TRANSIENT) {      
      TransientDriver * tD = new TransientDriver( sequenceStep,true, 
                                                 simState_, domain_, seqNode, info  );
      // pass accumulated time
      tD->SetAccumulatedTime( accumulatedTime_ );
      actDriver_ = tD;
    }
    else if (analysisPerStep_[sequenceStep] == BasePDE::HARMONIC) {
      actDriver_ = new HarmonicDriver( sequenceStep, true, 
                                       simState_, domain_, seqNode, info  );
    }
    else if( analysisPerStep_[sequenceStep] == BasePDE::EIGENFREQUENCY ) {
      actDriver_ = new EigenFrequencyDriver( sequenceStep, true, 
                                             simState_, domain_, seqNode, info  );
    }
    LOG_DBG(msDriver) << "\tCreated driver of type " 
        << BasePDE::analysisType.ToString(actDriver_->GetAnalysisType())
        << "(adress: " << actDriver_ << ")";

    // Initialize all PDEs
    LOG_DBG(msDriver) << "\tCreating PDEs"; 
    domain_->CreatePDEs( sequenceStep, info );

    domain_->SetDriver( actDriver_ );
    actDriver_->SetPDE(domain_->GetBasePDE());

    // Initialize driver objects
    LOG_DBG(msDriver) << "\tInitializing single driver (isRestarted = " 
                        << (isRestarted_ ? "true" :"false") << ")";
    actDriver_->Init(isRestarted_);

    //! Initialize Pdes, after having set the memento object
    LOG_DBG(msDriver) << "\tInitializing PDEs";
    domain_->InitPDEs( sequenceStep);
    
    LOG_TRACE(msDriver) << "Finished SetupStep";

  }
  
  void MultiSequenceDriver::ReadRestart() {
    
    // check if we are in restart mode
    if( !isRestarted_)
      return;
    LOG_TRACE(msDriver) << "Reading Restart information...";
    
    // Delete old driver if present
    if( actDriver_ ) {
        delete actDriver_;
        actDriver_ = NULL;
      }
    
    // Create input reader from current output reader
    bool hasInput = simState_->SetInputReaderToSameOutput();
    if( !hasInput)  {
      EXCEPTION( "Can not perform restarted simulation, as HDF5 file "
          << "contains no restart information.");
    }
    
    // Get all analysis types
    std::map<UInt, BasePDE::AnalysisType> analysis;
    std::map<UInt, Double> accTime;
    std::map<UInt, bool> isFinished;
    simState_->GetSequenceSteps( analysis, accTime, isFinished);
    std::map<UInt, BasePDE::AnalysisType>::const_iterator it;
    it = analysis.begin();
    UInt lastStepFinished = 0;
    for( ; it != analysis.end(); ++it ) {
      UInt actMsStep = it->first;
      // set sequenceStep
      if( isFinished[actMsStep]) {
        accumulatedTime_ = accTime[actMsStep];
        lastStepFinished = actMsStep;
      }
      // we always set the last multisequence step to continue
      sequenceStep_ = actMsStep;
    }
    
    // additional check: In case the last multisequence "finished" in the 
    // HDF5 file is the last one of the total simulation, it could still be
    // extended (e.g. by adding more steps by hand in the xml file). 
    // In this case we set the accumulated time to the previous step and
    // let the last SingeleDriver determine, if mroe steps have to be simulated.
    if( lastStepFinished == numSteps_ ) {
      accumulatedTime_ = accTime[lastStepFinished-1];
    }
    
    LOG_DBG(msDriver) << "\tSequence Step to continue: " << sequenceStep_;
    LOG_DBG(msDriver) << "\tAccumulated time so far: " << accumulatedTime_;
    // Note: the following piece of code is not necessary, as always
    //       the SingleDriver is repsonsible to determine, if a 
    //       a step has to be continued.
    
    // Check in the end, if all multisequence steps are already finished
    //if( lastFinishedStep == numSteps_ ) {
    //  std::cout << "\n\n";
    //  std::cout<< "*******************************************************\n";
    //  std::cout << " No restart necessary, as the desired number of \n";
    //  std::cout << " multi sequence steps are already computed. \n";
    //  std::cout << "*******************************************************\n\n";
      
    //  sequenceStep_ = numSteps_+1;
    //}
    
  }
  
  void MultiSequenceDriver::WriteMultiSequenceStep(const UInt sequenceStep, 
                                  const BasePDE::AnalysisType analysis) {
    std::string analysisString = BasePDE::analysisType.ToString(analysis);
    // write std::out info 
    std::cout << std::endl 
        << " ***************************** " << std::endl
        << " MultiSequenceStep: " << sequenceStep << std::endl
        << " AnalysisType:      " << analysisString << std::endl
        << " ***************************** " << std::endl << std::endl;
  }

} // end of namespace
