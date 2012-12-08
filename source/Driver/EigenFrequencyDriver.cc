#include "EigenFrequencyDriver.hh"

#include <iostream>
#include <iomanip>

#include "Driver/SolveSteps/StdSolveStep.hh"

#include "Domain/Domain.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ResultHandler.hh"

#include "PDE/StdPDE.hh"

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
    PtrParamNode myNode = 
      param->GetByVal("sequenceStep", std::string("index"), sequenceStep_)
        ->Get("analysis")->Get("eigenFrequency");

    // read required parameters from parameter node
    myNode->GetValue( "numModes", numFreq_ );
    myNode->GetValue( "freqShift", freqShift_ );
    myNode->GetValue( "writeModes", writeModes_, ParamNode::PASS );
    myNode->GetValue( "isQuadratic", isQuadratic_, ParamNode::PASS );

    InitializePDEs();
  }


  // **********************
  //   Default destructor
  // **********************
  EigenFrequencyDriver::~EigenFrequencyDriver() {
  }


  template <>
  void EigenFrequencyDriver::PrintResult<Double>(SingleVector* freq_ptr, Vector<Double> errBounds, 
                                                 ResultHandler* resHandler, UInt numConverged)
                                                 {
    Vector<Double>& eigenFreqs = dynamic_cast<Vector<Double>&>(*freq_ptr);

    // If no frequency at all converged, just leave
    if( numConverged == 0) {
      WARN( "No eigenfrequency converged, so no output will be written to "
          "the result files." );
      return;
    }

    // Issue warning, if number of converged eigenvalues differs from
    // number of requested ones
    if( numConverged != numFreq_ ) {
      WARN( "Only " << numConverged << " eigenfrequencies of " 
            << numFreq_ << " converged. To improve convergence, either "
            << "reduce the number of eigenfrequencies or the tolerance." );
    }

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
      PtrParamNode result = driverNode->Get("result",ParamNode::APPEND);
      result->Get("frequency")->SetValue(eigenFreqs[i]);
      result->Get("errorbound")->SetValue(errBounds[i]);
    }


    // ------------------------------
    // Phase 2: calculate eigenmodes
    // ------------------------------
    if ( writeModes_ == true ) {

      for ( UInt i = 0 ; i < numConverged; i++ ) {
        // Set current frequency value in the mathParser
        domain->GetMathParser()->SetValue( MathParser::GLOB_HANDLER, "f", abs(eigenFreqs[i]) );
        domain->GetMathParser()->SetValue( MathParser::GLOB_HANDLER, "step", i+1 );
        
        ptPDE_->GetSolveStep()->SetActStep(i);
        ptPDE_->GetSolveStep()->SetActFreq(std::abs(eigenFreqs[i]));
        ptPDE_->GetSolveStep()->CalcEigenMode( i );
        resHandler->BeginStep( i+1, std::abs(eigenFreqs[i]) );
        ptPDE_->WriteResultsInFile(i+1, std::abs(eigenFreqs[i]) );
        resHandler->FinishStep( );
      }
    }
                                                 }

  template <>
  void EigenFrequencyDriver::PrintResult<Complex>(SingleVector* freq_ptr, Vector<Double> errBounds, 
                                                  ResultHandler* resHandler, UInt numConverged)
                                                  {
    Vector<Complex>& eigenFreqs = dynamic_cast<Vector<Complex>&>(*freq_ptr);

    // If no frequency at all converged, just leave
    if( numConverged == 0) {
      WARN( "No eigenfrequency converged, so no output will be written to "
          "the result files." );
      return;
    }

    // Issue warning, if number of converged eigenvalues differs from
    // number of requested ones
    if( numConverged != numFreq_ ) {
      WARN( "Only " << numConverged << " eigenfrequencies of " 
            << numFreq_ << " converged. To improve convergence, either "
            << "reduce the number of eigenfrequencies or the tolerance." );
    }

    // notify resultHandler about beginning of new sequence step
    resHandler->BeginMultiSequenceStep( sequenceStep_,
                                        analysis_,
                                        numConverged );

    // Print out eigenfrequencies
    std::cout << std::endl << std::endl;

    std::cout << std::setw(20) << "Frequency [Hz]" << " | ";
    std::cout << std::setw(20) << "Damping       " << " | ";
    std::cout << std::setw(20) << "Errorbound";
    std::cout << "\n";
    std::cout << std::setfill('-') << std::setw(70) << "" << std::setfill(' ');
    std::cout << "\n";

    for( UInt i=0; i<numConverged; i++ ) {
      Double freq = eigenFreqs[i].imag()/(8.0*atan(1.0));
      Double damp = eigenFreqs[i].real();

      std::cout << std::setw(20) <<  freq;
      std::cout <<" | ";
      std::cout << std::setw(20) << damp;
      std::cout <<" | ";
      std::cout << std::setw(20) << errBounds[i] <<  "\n";

      // also log via info node
      PtrParamNode result = driverNode->Get("result",ParamNode::APPEND);
      result->Get("frequency")->SetValue(freq);
      result->Get("damping")->SetValue(damp);
      result->Get("errorbound")->SetValue(errBounds[i]);
    }

    // ------------------------------
    // Phase 2: calculate eigenmodes
    // ------------------------------
    if ( writeModes_ == true ) {

      for ( UInt i = 0 ; i < numConverged; i++ ) {
        Double actFreq = eigenFreqs[i].imag()/(8.0*atan(1.0));
        
        // Set current frequency value in the mathParser
        domain->GetMathParser()->SetValue( MathParser::GLOB_HANDLER, "f", actFreq );
        domain->GetMathParser()->SetValue( MathParser::GLOB_HANDLER, "step", i+1 );
        
        ptPDE_->GetSolveStep()->SetActStep(i);
        ptPDE_->GetSolveStep()->SetActFreq(actFreq);
        ptPDE_->GetSolveStep()->CalcEigenMode( i );
        resHandler->BeginStep( i+1, actFreq );
        ptPDE_->WriteResultsInFile(i+1, actFreq );
        resHandler->FinishStep( );
      }
    }
                                                  }
  
  // *****************
  //   Solve problem
  // *****************
  void EigenFrequencyDriver::SolveProblem(bool write_results, PtrParamNode given_analysis_id, AdjointParameters* adjointParams) {
    // options not implemented
    assert(write_results == true);
    
    if(given_analysis_id == NULL)
      analysis_id_ = driverNode->Get(ParamNode::PROCESS);
    else
      analysis_id_ = given_analysis_id;

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
    try {
      if(isQuadratic_) {
        numConverged = step->CalcEigenFrequencies( dynamic_cast<Vector<Complex>& >(*eigenFreqs),
                                                   errBounds,numFreq_, freqShift_ );
        PrintResult<Complex>(eigenFreqs, errBounds, resHandler, numConverged);
      }
      else  {
        numConverged = step->CalcEigenFrequencies( dynamic_cast<Vector<Double>& >(*eigenFreqs),
                                                   errBounds,numFreq_, freqShift_ );
        PrintResult<Double>(eigenFreqs, errBounds, resHandler, numConverged);
      } 
    } catch (Exception &ex ) {
      RETHROW_EXCEPTION(ex, "Could not calculate eigenfrequencies of setup");
    }
    
    // notify resultHandler about finishing of current sequence step
    if( !isPartOfSequence_ ) {
      resHandler->FinishMultiSequenceStep();
      resHandler->Finalize();
    }
  }

  
 
  
  
} // end of namespace
