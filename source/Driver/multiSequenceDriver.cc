// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <assert.h>
#include <iostream>
#include <set>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/resultHandler.hh"
#include "Domain/domain.hh"
#include "Driver/basedriver.hh"
#include "Driver/singleDriver.hh"
#include "General/Enum.hh"
#include "General/environment.hh"
#include "General/exception.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/pdememento.hh"
#include "eigenFrequencyDriver.hh"
#include "harmonicDriver.hh"
#include "multiSequenceDriver.hh"
// include all kinds of single drivers
#include "staticdriver.hh"
#include "transientdriver.hh"

namespace CoupledField {

  // ***************
  //   Constructor
  // ***************
class AdjointParameters;

  MultiSequenceDriver::MultiSequenceDriver( ) 
    : BaseDriver() {

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


  void MultiSequenceDriver::SolveProblem(bool write_results, PtrParamNode given_analysis_id, AdjointParameters* adjointParams) {
    // options not implemented
    assert(write_results == true);


    // vector containing pointer to current set of PDEs
    StdVector<SinglePDE *> ptPDEs;
  
    // Vector of memento objects, which save the internal state
    // of a PDE
    StdVector<shared_ptr<PDEMemento> > memento;

    std::cout << "++ Starting to solve problem" << std::endl;

    // get resultHandler
    ResultHandler * resHandler = domain->GetResultHandler();

    // outer loop over all single sequences
    for (UInt iStep = 0; iStep < numSteps_; iStep++ ) {
    
      // remeber current sequence step
      curSequenceStep_ = iStep;

      // log info about new step to info-class
      Info->WriteMultiSequenceStep(iStep+1, 
                                   analysisPerStep_[iStep]);
      
      // Since per time step only one type of 
      // analysis is allowed, we simple access
      // the first entry fo analysisPerStep_
      if (analysisPerStep_[iStep] == BasePDE::STATIC) {
        actDriver_ = new StaticDriver( iStep+1, true );
      }
      else if (analysisPerStep_[iStep] == BasePDE::TRANSIENT) {      
        actDriver_ = new TransientDriver( iStep+1,true );
      }
      else if (analysisPerStep_[iStep] == BasePDE::HARMONIC) {
        actDriver_ = new HarmonicDriver( iStep+1, true );
      }
      else if( analysisPerStep_[iStep] == BasePDE::EIGENFREQUENCY ) {
        actDriver_ = new EigenFrequencyDriver( iStep+1, true );
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

      // After the first run, initialize this PDE
      // with the solution of the previous run
      if (iStep > 0) {
        for (UInt iPde = 0; iPde < pdesPerStep_[iStep].GetSize(); iPde++ ) {
          if( memento[iPde] != NULL) {
            if (memento[iPde]->IsSet()) {
              ptPDEs[iPde]->SetMemento( memento[iPde], 
                                        usageDirichletPerStep_[iStep] );
            }
          }
        }
      }

      //! Initialize Pdes, after having set the memento object
      domain->InitPDEs( iStep+1 );
    
      // Solve Problem
      actDriver_->SolveProblem(write_results, given_analysis_id);

      // Get solution for next step and delete
      // all PDEs
      if (iStep < numSteps_-1) {
        memento.Resize(pdesPerStep_[iStep+1].GetSize());

        // Iterate over all PDEs in the next step
        for (UInt iPde = 0; iPde < pdesPerStep_[iStep+1].GetSize(); iPde++) {
          // Iterate over all PDEs in the current step
          for(UInt kPde = 0; kPde < pdesPerStep_[iStep].GetSize(); kPde++) {
            // If both match, then save the result of this step
            // for the next step
            if (pdesPerStep_[iStep+1][iPde] == pdesPerStep_[iStep][kPde])
              //dynamic_cast<const NodeStoreSol<Double>& >
              ptPDEs[kPde]->GetMemento(memento[iPde]);
          }
        }
      
        // delete PDEs
        domain->ResetPDEs();
      }


      // finish sequence step with resulthandler
      resHandler->FinishMultiSequenceStep();

      // delete analysistypes
      delete actDriver_;
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
    ParamNodeList seqNodes = param->GetList("sequenceStep");


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
        param->GetByVal("sequenceStep", std::string("index"), iStep+1);

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
