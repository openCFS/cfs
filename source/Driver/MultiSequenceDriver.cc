#include "MultiSequenceDriver.hh"

#include <set>

// include all kinds of single drivers
#include "EigenValueDriver.hh"
#include "StaticDriver.hh"
#include "TransientDriver.hh"
#include "HarmonicDriver.hh"
#include "MultiHarmonicDriver.hh"
#include "EigenFrequencyDriver.hh"
#include "InverseSourceDriver.hh"
#include "Utils/preciceAdapter/IPreciceAdapter.hh"
#include "Utils/preciceAdapter/TransientDriverPrecice.hh"
#include "BucklingDriver.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "PDE/SinglePDE.hh"
#include "Domain/Domain.hh"
#include "DataInOut/ResultHandler.hh"

#include "DataInOut/Logging/LogConfigurator.hh"

namespace CoupledField {

DEFINE_LOG(msDriver, "msDriver")

  // ***************
  //   Constructor
  // ***************
  MultiSequenceDriver::MultiSequenceDriver(shared_ptr<SimState> state, Domain* domain,
                                            PtrParamNode paramNode, PtrParamNode infoNode, bool keep)
    : BaseDriver(state, domain, paramNode, infoNode)
  {
    analysis_ = BasePDE::MULTI_SEQUENCE;
    
    sequenceStep_ = 1;
    numSteps_ = 0;
    accumulatedTime_ = 0.0;
    actDriver_ = NULL;
    isRestarted_ = false;
    keep_ = domain->GetParamRoot()->Has("optimization");
  }

  // **********************
  //   Default destructor
  // **********************
  MultiSequenceDriver::~MultiSequenceDriver()
  {
    // the keep_ results in set kept*_ arrays. Handle only the non-actual stuff
    for(unsigned int i = 0; i < keptDrivers_.GetSize(); i++)
      if(i != sequenceStep_-1) // the actDriver_ does not need to be keptDrivers_[0]
        delete keptDrivers_[i]; // includes the driver in actDriver_

    for(unsigned int i = 0; i < keptPDEs_.GetSize(); i++)
      if(i != sequenceStep_-1) // the currently set pdes are deleteted in ~Domain()
        for(unsigned int p = 0; p < keptPDEs_[i].GetSize(); p++)
          delete keptPDEs_[i][p];

    delete actDriver_; // It can happen, that the driver was not deleted (e.g. is SolveProblem was not called)
  }


  void MultiSequenceDriver::SolveProblem()
  {
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
      
      // only update in the multisequence case
      domain_->GetPreciceAdapter()->UpdateDomain(domain_);
      
      // Solve Problem
      actDriver_->SolveProblem();

      // add up accumulated time
      if (analysisPerStep_[iStep] == BasePDE::TRANSIENT) {
        TransientDriver * tD = dynamic_cast<TransientDriver*>(actDriver_);
        accumulatedTime_ += tD->GetDuration();
      }
      
//      // finish sequence step with result handler
//      resHandler->FinishMultiSequenceStep();
//      simState_->FinishMultiSequenceStep( !abortSimulation_, accumulatedTime_ );
      
      // Get solution for next step and delete all PDEs
      if (iStep < numSteps_) {
        // delete PDEs
        domain_->ResetPDEs(keep_);
      }

      // delete analysis types
      if(!keep_)
        delete actDriver_;
      actDriver_ = NULL; // in the keep_ case we have a copy in keptDrivers_
      
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
    assert(actDriver_);
    return actDriver_->GetActStep( pdename );
  }

  // *****************
  //   Initializer
  // *****************
  void MultiSequenceDriver::Init(bool restart) 
  {
    LOG_DBG(msDriver) << "Init: restart= " << restart;
    assert(actDriver_ == NULL && keptDrivers_.IsEmpty()); // make sure call this only once!
    
    isRestarted_ = restart;

    // get nodes for all sequence steps
    ParamNodeList seqNodes = domain_->GetParamRoot()->GetList("sequenceStep");

    if(keep_) {
      keptDrivers_.Resize(seqNodes.GetSize(), NULL); // important to initialize with zeros!
      keptPDEs_.Resize(seqNodes.GetSize());
    }

    // Fill vector with step indices and ensure that all occur.
    // One based, compact but arbitrary order
    std::set<UInt> stepIndices;
    for(UInt i = 0; i < seqNodes.GetSize(); i++ )
    {
      UInt actStepIndex = seqNodes[i]->Get("index")->As<UInt>();
      paramPerStep[actStepIndex] = seqNodes[i];
      if(stepIndices.find(actStepIndex) != stepIndices.end())
        EXCEPTION("Multisequence step with index " << actStepIndex << " occurs more than once!");
      stepIndices.insert(actStepIndex);
    }
    // verify the 1-based, compatness requirement
    for(UInt i = 1; i <= seqNodes.GetSize(); i++ )
      if(stepIndices.find(i) == stepIndices.end())
        EXCEPTION("Multisequence step  with index " << i << " is not defined!" );

    numSteps_ = stepIndices.size();

    // 3.) Read in all pdes and analysis types

    // iterate over all sequence Steps
    for( UInt iStep = 1; iStep <= numSteps_; iStep++) {

      // get current step node
      PtrParamNode actStepNode = paramPerStep[iStep];

      // get current analysistype
      std::string analysisString = actStepNode->Get("analysis")->GetChild()->GetName();
      analysisPerStep_[iStep] = BasePDE::analysisType.Parse(analysisString);

      // get all pde-nodes in current sequence step
      ParamNodeList pdeNodes;
      pdeNodes = actStepNode->Get("pdeList")->GetChildren();

      // iterate over all pdes
      pdesPerStep_[iStep].Resize( pdeNodes.GetSize() );
      for( UInt iPde = 0; iPde < pdeNodes.GetSize(); iPde++ )
      {
        // get pdeName
        std::string pdeName = pdeNodes[iPde]->GetName();
        pdesPerStep_[iStep][iPde] = pdeName;
      } // loop over pdes
    } // loop over sequence steps

    // 4.) Perform final consistency checks
    // Not much to do here yet ...
    
  }
  
  void MultiSequenceDriver::SetSequenceStep(UInt sequenceStep )
  {
    LOG_TRACE(msDriver) << "Calling SetSequenceStep with step " << sequenceStep;
    if(sequenceStep > numSteps_)
      EXCEPTION( "Can not update to multistep " << sequenceStep << " as only " << numSteps_ << " are defined" );

    SetupStep(sequenceStep);
    sequenceStep_ = sequenceStep;
    
  }
  
  void MultiSequenceDriver::SetupStep(UInt sequenceStep) {

    LOG_TRACE(msDriver) << "SS: SetupStep for MS=" << sequenceStep << " current=" << sequenceStep_ << " driver=" << (actDriver_ != NULL ? actDriver_->GetAnalysisType() : -1);
    // Only setup step once
    if(sequenceStep_ == sequenceStep && actDriver_ !=  NULL)
      return;
    
    if(!keep_)
      delete actDriver_;
    actDriver_ = NULL;

    sequenceStep_ = sequenceStep;

    // in the keep_ case we create only once in the non keep_ case we create always
    bool create_new = !keep_ || keptPDEs_[sequenceStep-1].IsEmpty();
    
    // obtain paramNode and infoNode
    PtrParamNode seqNode = paramPerStep[sequenceStep]->Get("analysis");
    PtrParamNode info = info_->GetByVal("sequenceStep","sequence", sequenceStep);
    
    if(create_new)
    {
      // Since per sequence step only one type of analysis is allowed, we simple access the first entry of analysisPerStep_
      if (analysisPerStep_[sequenceStep] == BasePDE::STATIC) {
        actDriver_ = new StaticDriver(sequenceStep, true, simState_, domain_, seqNode, info);
      }
      else if (analysisPerStep_[sequenceStep] == BasePDE::TRANSIENT) {
        if((!domain_->GetPreciceAdapter()->IsPreciceDummy()) && (domain_->GetPreciceAdapter()->GetSequenceStep() == sequenceStep)){
          TransientDriverPrecice * tD = new TransientDriverPrecice( sequenceStep, true, simState_, domain_, seqNode, info );
          actDriver_ = tD;
        }else{
          TransientDriver * tD = new TransientDriver( sequenceStep,true, simState_, domain_, seqNode, info  );
          // pass accumulated time
          tD->SetAccumulatedTime( accumulatedTime_ );
          actDriver_ = tD;
        }
        
      }
      else if (analysisPerStep_[sequenceStep] == BasePDE::HARMONIC) {
        actDriver_ = new HarmonicDriver( sequenceStep, true, simState_, domain_, seqNode, info  );
      }
      else if( analysisPerStep_[sequenceStep] == BasePDE::EIGENFREQUENCY ) {
        actDriver_ = new EigenFrequencyDriver( sequenceStep, true, simState_, domain_, seqNode, info  );
      }
      else if( analysisPerStep_[sequenceStep] == BasePDE::BUCKLING ) {
        actDriver_ = new BucklingDriver( sequenceStep, true, simState_, domain_, seqNode, info  );
      }
      else if( analysisPerStep_[sequenceStep] == BasePDE::EIGENVALUE ) {
        actDriver_ = new EigenValueDriver( sequenceStep, true, simState_, domain_, seqNode, info  );
      }
      else if( analysisPerStep_[sequenceStep] == BasePDE::MULTIHARMONIC ) {
        actDriver_ = new MultiHarmonicDriver( sequenceStep, true, simState_, domain_, seqNode, info  );
      }else{
        EXCEPTION("MultiSequenceDriver::SetupStep Could not find the correct driver");
      }

      LOG_DBG(msDriver) << "SS: created driver of type "  << BasePDE::analysisType.ToString(actDriver_->GetAnalysisType()) << "(adress: " << actDriver_ << ")";
      if(keep_)
        keptDrivers_[sequenceStep-1] = actDriver_;
    }
    else
    {
      actDriver_ = keptDrivers_[sequenceStep-1]; // reuse!
      LOG_DBG(msDriver) << "SS: reuse existing driver "  << BasePDE::analysisType.ToString(actDriver_->GetAnalysisType()) << " (address: " << actDriver_ << ")";
    }

    domain_->SetDriver(actDriver_);

    // Initialize all PDEs
    if(create_new) {
      LOG_DBG(msDriver) << "SS: creating PDEs.";
      domain_->CreatePDEs( sequenceStep, info );
      if(keep_)
        keptPDEs_[sequenceStep-1] = domain_->GetSinglePDEs();
    }
    else {
      LOG_DBG(msDriver) << "SS: restoring PDEs.";
      domain_->RestorePDEs(keptPDEs_[sequenceStep-1]);
    }

    actDriver_->SetPDE(domain_->GetBasePDE());

    if(create_new)
    {
      // Initialize driver objects
      LOG_DBG(msDriver) << "\tInitializing single driver (isRestarted = " << (isRestarted_ ? "true" :"false") << ")";
      actDriver_->Init(isRestarted_);

      //! Initialize Pdes, after having set the memento object

      LOG_DBG(msDriver) << "\tInitializing PDEs";
      domain_->InitPDEs( sequenceStep);
    }
    
    LOG_TRACE(msDriver) << "Finished SetupStep";
  }
  
  void MultiSequenceDriver::ReadRestart() {
    
    // check if we are in restart mode
    if( !isRestarted_)
      return;
    LOG_TRACE(msDriver) << "Reading Restart information...";
    
    // Delete old driver if present
    if(!keep_)
      delete actDriver_;

    actDriver_ = NULL;
    
    // Create input reader from current output reader
    if(!simState_->SetInputReaderToSameOutput())
      EXCEPTION("Can not perform restarted simulation, as HDF5 file contains no restart information.");

    
    // Get all analysis types
    std::map<UInt, BasePDE::AnalysisType> analysis;
    std::map<UInt, Double> accTime;
    std::map<UInt, bool> isFinished;
    simState_->GetSequenceSteps( analysis, accTime, isFinished);
    std::map<UInt, BasePDE::AnalysisType>::const_iterator it;
    it = analysis.begin();
    UInt lastStepFinished = 0;
    for( ; it != analysis.end(); ++it )
    {
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
