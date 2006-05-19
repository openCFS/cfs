#include "eigenFrequencyDriver.hh"
#include "stdSolveStep.hh"

#include "PDE/StdPDE.hh"
#include "Domain/domain.hh"

#include "DataInOut/CommandLine/BaseCommandLineHandler.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"

namespace CoupledField {

  // ***************
  //   Constructor
  // ***************
  EigenFrequencyDriver::EigenFrequencyDriver(Domain * adomain, 
                             UInt stepOffset,
                             Double timeOffset,
                             std::string driverTag,
                             Boolean isPartOfSequence) 
    : SingleDriver(adomain, stepOffset, timeOffset, 
                   driverTag, isPartOfSequence) {
    ENTER_FCN( "EigenFrequencyDriver::EigenFrequencyDriver" );

    // vecotrs for accessing parameters
    StdVector<std::string> keyVec, attrVec, valVec;

    attrVec = "tag";
    valVec = driverTag_;

    // Get number of eigenmodes
    keyVec = "eigenFrequency", "numModes";
    params->Get( keyVec, attrVec, valVec, numFreq_ );

    // Get frequency shift
    keyVec = "eigenFrequency", "freqShift";
    params->Get( keyVec, attrVec, valVec, freqShift_ );

    // Get flag for writing out the modes
    std::string temp;
    keyVec = "eigenFrequency", "writeModes";
    params->Get( keyVec, attrVec, valVec, temp );

    if ( temp == "yes" ) {
      writeModes_ = TRUE;
    } else {
      writeModes_ = FALSE;
    }
    
    // Get flag for writing out the modes
    keyVec = "eigenFrequency", "isQuadratic";
    params->Get( keyVec, attrVec, valVec, temp );
    
    if ( temp == "yes" ) {
      isQuadratic_ = TRUE;
    } else {
      isQuadratic_ = FALSE;
    }
    
    
  }


  // **********************
  //   Default destructor
  // **********************
  EigenFrequencyDriver::~EigenFrequencyDriver() {
    ENTER_FCN( "EigenFrequencyDriver::~EigenFrequencyDriver" );
  }


  // *****************
  //   Solve problem
  // *****************
  void EigenFrequencyDriver::SolveProblem() {
    ENTER_FCN( "EigenFrequencyDriver::SolveProblem" );

    // if driver is not part of multiSequence Driver, get list
    // of pdes which have to be solved and intialize them
    if (isPartOfSequence_ == FALSE){     
             GetMyPDEs();
             Info->StartProgress ("Starting to solve problem", FALSE);
    }
    
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
      if ( writeModes_ == TRUE ) {
        
        // Iterate over all frequencies an calculate according mode
        if (! isPartOfSequence_)
          ptdomain_->PrintGrid();
        
        for ( UInt i = 0 ; i < numConverged; i++ ) {
          ptPDE_->GetSolveStep()->SetActStep(i);
          ptPDE_->GetSolveStep()->SetActFreq(std::abs(eigenFreqs[i]));
          ptPDE_->GetSolveStep()->CalcEigenMode( i );
          ptPDE_->WriteResultsInFile(i+1, std::abs(eigenFreqs[i]), 
                                     stepOffset_, timeOffset_);
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
      if ( writeModes_ == TRUE ) {
        
        // Iterate over all frequencies an calculate according mode
        if (! isPartOfSequence_)
          ptdomain_->PrintGrid();
        
        for ( UInt i = 0 ; i < numConverged; i++ ) {
          ptPDE_->GetSolveStep()->SetActStep(i);
          ptPDE_->GetSolveStep()->SetActTime(eigenFreqs[i]);
          ptPDE_->GetSolveStep()->CalcEigenMode( i );
          ptPDE_->WriteResultsInFile(i+1, eigenFreqs[i], 
                                     stepOffset_, timeOffset_);
        }
      }
    }
    
  }


} // end of namespace
