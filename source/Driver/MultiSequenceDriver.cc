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

namespace CoupledField {

  // ***************
  //   Constructor
  // ***************
  MultiSequenceDriver::MultiSequenceDriver( shared_ptr<SimState> state,
                                            Domain* domain) 
    : BaseDriver(state, domain) {

    analysis_ = BasePDE::MULTI_SEQUENCE;
    
    numSteps_ = 0;
    accumulatedTime_ = 0.0;
    actDriver_ = 0;
  }


  // **********************
  //   Default destructor
  // **********************
  MultiSequenceDriver::~MultiSequenceDriver() {
  }


  void MultiSequenceDriver::SolveProblem(bool write_results, PtrParamNode given_analysis_id) {
    // options not implemented
    assert(write_results == true);


    // vector containing pointer to current set of PDEs
    StdVector<SinglePDE *> ptPDEs;

    std::cout << "++ Starting to solve problem" << std::endl;

    // get resultHandler
    ResultHandler * resHandler = domain->GetResultHandler();

    // outer loop over all single sequences
    for (UInt iStep = 0; iStep < numSteps_; iStep++ ) {

      // remember current sequence step
      curSequenceStep_ = iStep;
      
      // log info about new step to info-class
      WriteMultiSequenceStep(iStep+1, analysisPerStep_[iStep]);

      // Since per time step only one type of 
      // analysis is allowed, we simple access
      // the first entry fo analysisPerStep_
      if (analysisPerStep_[iStep] == BasePDE::STATIC) {
        actDriver_ = new StaticDriver( iStep+1, true, 
                                       simState_, domain_ );
      }
      else if (analysisPerStep_[iStep] == BasePDE::TRANSIENT) {      
        TransientDriver * tD = new TransientDriver( iStep+1,true, 
                                                   simState_, domain_  );
        // pass accumulated time
        tD->SetAccumulatedTime( accumulatedTime_ );
        actDriver_ = tD;
      }
      else if (analysisPerStep_[iStep] == BasePDE::HARMONIC) {
        actDriver_ = new HarmonicDriver( iStep+1, true, 
                                         simState_, domain_  );
      }
      else if( analysisPerStep_[iStep] == BasePDE::EIGENFREQUENCY ) {
        actDriver_ = new EigenFrequencyDriver( iStep+1, true, 
                                               simState_, domain_  );
      }

      // Initialize all PDEs
      domain->CreatePDEs( iStep+1 );

      domain->SetDriver( actDriver_ );
      actDriver_->SetPDE(domain->GetBasePDE());

      // Give the according pdes to the driver
      ptPDEs.Clear();
      ptPDEs.Resize(pdesPerStep_[iStep].GetSize());
      for (UInt iPde=0; iPde < pdesPerStep_[iStep].GetSize(); iPde++)
        ptPDEs[iPde]=domain->GetSinglePDE( pdesPerStep_[iStep][iPde] );

      // Initialize driver objects
      actDriver_->Init();


      //! Initialize Pdes, after having set the memento object
      domain->InitPDEs( iStep+1 );

      // Solve Problem
      actDriver_->SolveProblem(write_results, given_analysis_id);

      // finish sequence step with resulthandler
      resHandler->FinishMultiSequenceStep();
      
      // Get solution for next step and delete all PDEs
      if (iStep < numSteps_-1) {
        // delete PDEs
        domain->ResetPDEs();
      }

      // add up accumulated time
      if (analysisPerStep_[iStep] == BasePDE::TRANSIENT) {
        TransientDriver * tD = dynamic_cast<TransientDriver*>(actDriver_);
        accumulatedTime_ += tD->GetDuration();
      }
      
      // delete analysis types
      delete actDriver_;
      actDriver_ = NULL;
    } // iStep

    // trigger finalization
    domain->GetResultHandler()->Finalize();
  }


  UInt MultiSequenceDriver::GetActStep( const std::string& pdename ) {
    assert( actDriver_);
    return actDriver_->GetActStep( pdename );
  }

  // *****************
  //   Initializer
  // *****************
  void MultiSequenceDriver::Init() 
  {

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

    // 2.) Resize 'outer' vectors
    pdesPerStep_.Resize(numSteps_);
    analysisPerStep_.Resize(numSteps_);


    // 3.) Read in all pdes and analysis types

    // iterate over all sequence Steps
    for( UInt iStep = 0; iStep < numSteps_; iStep++) {

      // get current step node
      PtrParamNode actStepNode = 
          domain_->GetParamRoot()->GetByVal("sequenceStep", std::string("index"), iStep+1);

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
