// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "eigenFrequencyDriver.hh"
#include "stdSolveStep.hh"

#include <iostream>
#include <iomanip>

#include "PDE/StdPDE.hh"
#include "Domain/domain.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/InfoNode.hh"
#include "DataInOut/resultHandler.hh"

namespace CoupledField {

  // ***************
  //   Constructor
  // ***************
  EigenFrequencyDriver::EigenFrequencyDriver(UInt sequenceStep,
                                             bool isPartOfSequence) 
    : SingleDriver( sequenceStep, isPartOfSequence ) {

    // set correct analysistype
    analysis_ = BasePDE::EIGENFREQUENCY;

    numFreq_ = 0;
    freqShift_ = 0.0;
    writeModes_ = true;
    isQuadratic_ = false;
    // replace with a concrete element
    driverNode = driverNode->Get("eigenFrequency");
  }


  void EigenFrequencyDriver::Init() {
    
    // get parameter node
    ParamNode * myNode = 
      param->Get("sequenceStep", std::string("index"), sequenceStep_)
        ->Get("analysis")->Get("eigenFrequency");

    // read required parameters from parameter node
    myNode->Get( "numModes", numFreq_ );
    myNode->Get( "freqShift", freqShift_ );
    myNode->Get( "writeModes", writeModes_, true );
    myNode->Get( "isQuadratic", isQuadratic_, true );

    InitializePDEs();
  }


  // **********************
  //   Default destructor
  // **********************
  EigenFrequencyDriver::~EigenFrequencyDriver() {
  }


  // *****************
  //   Solve problem
  // *****************
  void EigenFrequencyDriver::SolveProblem(bool write_results, InfoNode* given_analysis_id, const bool reAssembleMatrices) {
    // options not implemented
    assert(write_results == true);
    assert(given_analysis_id == NULL); // not implemented yet
    
    ResultHandler* resHandler = domain->GetResultHandler();

    // ------------------------------
    // Phase 1: calculate eigenvalues( generalized problem)
    // ------------------------------
    
    // the eigenfrequencies are complex in the quadratic case
    SingleVector* eigenFreqs = NULL;
    if(isQuadratic_) eigenFreqs = new Vector<Complex>(numFreq_); 
                else eigenFreqs = new Vector<Double>(numFreq_);
    Vector<Double> errBounds( numFreq_ );

    // Trigger calculation 
    UInt numConverged = 0;
    ptPDE_->WriteGeneralPDEdefines();
    BaseSolveStep* step = ptPDE_->GetSolveStep();
    if(isQuadratic_) {
      numConverged = step->CalcEigenFrequencies( dynamic_cast<Vector<Complex>& >(*eigenFreqs),
                                                 errBounds,numFreq_, freqShift_ );
      PrintResult<Double>(eigenFreqs, errBounds, resHandler, numConverged);
    }
    else  {
      numConverged = step->CalcEigenFrequencies( dynamic_cast<Vector<Double>& >(*eigenFreqs),
                                                 errBounds,numFreq_, freqShift_ );
      PrintResult<Double>(eigenFreqs, errBounds, resHandler, numConverged);
    }
    
    // notify resultHandler about finishing of current sequence step
    if( !isPartOfSequence_ ) {
      resHandler->FinishMultiSequenceStep();
      resHandler->Finalize();
    }
  }

  
  template <class T>
  void EigenFrequencyDriver::PrintResult(SingleVector* freq_ptr, Vector<Double> errBounds, 
                                         ResultHandler* resHandler, UInt numConverged)
  {
    Vector<T>& eigenFreqs = dynamic_cast<Vector<T>&>(*freq_ptr);
    
    
    // notify resultHandler about beginning of new sequence step
    resHandler->BeginMultiSequenceStep( sequenceStep_,
        analysis_,
        numConverged );

    // Print out eigenfrequencies
    std::cout << std::endl << std::endl;

    std::cout << std::setw(20) << "Frequency [Hz]" << " | ";
    std::cout << std::setw(20) << "Errorbound";
    std::cout << "\n";
    std::cout << std::setfill('-') << std::setw(40) << "" << std::setfill(' ');
    std::cout << "\n";

    for( UInt i=0; i<numConverged; i++ ) {
      std::cout << std::setw(20) << eigenFreqs[i];
      std::cout <<" | ";
      std::cout << std::setw(20) << errBounds[i] <<  "\n";

      // also log via info node
      InfoNode* result = driverNode->Get("result",InfoNode::APPEND);
      result->Get("frequency")->SetValue(eigenFreqs[i]);
      result->Get("errorbound")->SetValue(errBounds[i]);
    }

    
    // ------------------------------
    // Phase 2: calculate eigenmodes
    // ------------------------------
    if ( writeModes_ == true ) {

      for ( UInt i = 0 ; i < numConverged; i++ ) {
        ptPDE_->GetSolveStep()->SetActStep(i);
        ptPDE_->GetSolveStep()->SetActFreq(std::abs(eigenFreqs[i]));
        ptPDE_->GetSolveStep()->CalcEigenMode( i );
        resHandler->BeginStep( i+1, std::abs(eigenFreqs[i]) );
        ptPDE_->WriteResultsInFile(i+1, std::abs(eigenFreqs[i]) );
        resHandler->FinishStep( );
      }
    }
  }
  
  
} // end of namespace
