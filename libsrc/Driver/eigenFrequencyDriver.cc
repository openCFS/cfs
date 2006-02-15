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
    
    // Get flag for shift-mode
    keyVec = "eigenFrequency", "shiftMode";
    params->Get( keyVec, attrVec, valVec, temp );

    if ( temp == "yes" ) {
      shiftMode_ = TRUE;
    } else {
      shiftMode_ = FALSE;
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
    // Phase 1: calculate eigenvalues
    // ------------------------------
    
    // Create vector for eigenvalues
    Vector<Double> eigenFreqs( numFreq_ );

    // Trigger calculation 
    UInt numConverged = 0;
    ptPDE_->WriteGeneralPDEdefines();
    numConverged = ptPDE_->GetSolveStep()->
      CalcEigenFrequencies( eigenFreqs, numFreq_, freqShift_, shiftMode_ );

    // Temporary
    std::cout << std::endl << std::endl;
    std::cout << "Eigenfrequencies are: \n" << eigenFreqs << std::endl;

    // ------------------------------
    // Phase 2: calculate eigenmodes
    // ------------------------------
    if ( writeModes_ == TRUE ) {

      // Iterate over all frequencies an calculate according mode
      if (! isPartOfSequence_)
        ptdomain_->PrintGrid();
      
      for ( UInt i = 0 ; i < numConverged; i++ ) {
        ptPDE_->GetSolveStep()->CalcEigenMode( i );
        ptPDE_->WriteResultsInFile(i, eigenFreqs[i], 
                                   stepOffset_, timeOffset_);
      }
    }
    
  }


} // end of namespace
