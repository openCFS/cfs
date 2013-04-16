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
    actStep_ = 0;
    actTime_ = 0;
    actDriver_ = 0;
  }


  // **********************
  //   Default destructor
  // **********************
  MultiSequenceDriver::~MultiSequenceDriver() {
  }


  void MultiSequenceDriver::SolveProblem(bool write_results, PtrParamNode given_analysis_id) {
    REFACTOR;
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
    usageDirichletPerStep_.Resize( numSteps_ );
    

    // 3.) Read in all pdes and analysis types

    // iterate over all sequence Steps
    for( UInt iStep = 0; iStep < numSteps_; iStep++) {
      
      // get current step node
      PtrParamNode actStepNode = 
          domain_->GetParamRoot()->GetByVal("sequenceStep", std::string("index"), iStep+1);

      // get current usage type 
      std::string usageString;
      actStepNode->GetValue( "usage", usageString );
      if( usageString == "startValue") {
        usageDirichletPerStep_[iStep] = false;
      } else {
        usageDirichletPerStep_[iStep] = true;
      }
      
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
  
} // end of namespace
