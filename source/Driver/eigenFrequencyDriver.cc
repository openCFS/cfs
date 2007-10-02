// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "eigenFrequencyDriver.hh"
#include "stdSolveStep.hh"

#include "PDE/StdPDE.hh"
#include "Domain/domain.hh"

#include "DataInOut/CommandLine/BaseCommandLineHandler.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/resultHandler.hh"

namespace CoupledField {

  // ***************
  //   Constructor
  // ***************
  EigenFrequencyDriver::EigenFrequencyDriver(UInt sequenceStep,
                                             bool isPartOfSequence) 
    : SingleDriver( sequenceStep, isPartOfSequence ) {

    // set correct analysistype
    analysis_ = EIGENFREQUENCY;

    numFreq_ = 0;
    freqShift_ = 0.0;
    writeModes_ = true;
    isQuadratic_ = false;
  }


  void EigenFrequencyDriver::Init() {
    
    // get parameter node
    ParamNode * myNode = 
      param->Get("sequenceStep", "index", GenStr( sequenceStep_ ) )
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
  void EigenFrequencyDriver::SolveProblem() {

    ResultHandler * resHandler = domain->GetResultHandler();

    // ------------------------------
    // Phase 1: calculate eigenvalues( generalized problem)
    // ------------------------------
    
    if ( isQuadratic_ ) {

      // === QUADRATIC CASE ===
      // Create vector for eigenvalues and error bounds
      Vector<Complex> eigenFreqs( numFreq_ );
      Vector<Double> errBounds( numFreq_ );
      
      // Trigger calculation 
      UInt numConverged = 0;
      ptPDE_->WriteGeneralPDEdefines();
      numConverged = ptPDE_->GetSolveStep()->
        CalcEigenFrequencies( eigenFreqs, errBounds, 
                              numFreq_, freqShift_ );
      
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
          ptPDE_->Finalize();
        }
      }
    } else {

      // === GENERALIZED CASE ===
      // Create vector for eigenvalues and error bounds
      Vector<Double> eigenFreqs( numFreq_ );
      Vector<Double> errBounds( numFreq_ );
      
      // Trigger calculation 
      UInt numConverged = 0;
      ptPDE_->WriteGeneralPDEdefines();
      numConverged = ptPDE_->GetSolveStep()->
        CalcEigenFrequencies( eigenFreqs, errBounds, 
                              numFreq_, freqShift_);
      
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
        
      }
      
      // ------------------------------
      // Phase 2: calculate eigenmodes
      // ------------------------------
      if ( writeModes_ == true ) {
        
        for ( UInt i = 0 ; i < numConverged; i++ ) {
          ptPDE_->GetSolveStep()->SetActStep(i);
          ptPDE_->GetSolveStep()->SetActTime(eigenFreqs[i]);
          ptPDE_->GetSolveStep()->CalcEigenMode( i );
          resHandler->BeginStep( i+1, eigenFreqs[i] );
          ptPDE_->WriteResultsInFile(i+1, eigenFreqs[i] );
          resHandler->FinishStep( );
          
        }
        ptPDE_->Finalize();
      }
    }

    // notify resultHandler about finishing of current sequence step
    if( !isPartOfSequence_ ) {
      resHandler->FinishMultiSequenceStep();
      resHandler->Finalize();
    }
  }


} // end of namespace
