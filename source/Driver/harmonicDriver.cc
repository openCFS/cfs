// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <string>

#include "harmonicDriver.hh"
#include "stdSolveStep.hh"
#include "assemble.hh"

#include "DataInOut/CommandLine/BaseCommandLineHandler.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/WriteInfo.hh"
#include "PDE/StdPDE.hh"
#include "Domain/domain.hh"
#include "DataInOut/resultHandler.hh"

namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
  HarmonicDriver::HarmonicDriver( UInt sequenceStep,
                                  bool isPartOfSequence )

    : SingleDriver( sequenceStep, isPartOfSequence ) {
    
    
    // Set correct analysistype
    analysis_ = HARMONIC;

    startFreq_ = 0.0;
    stopFreq_ = 0.0;
    numFreq_ = 0;
    actFreqStep_ = 0;
  }

  void HarmonicDriver::Init() {

    // get parameter node
    ParamNode * myNode = 
      param->Get("sequenceStep", "index", GenStr( sequenceStep_ ) )
      ->Get("analysis")->Get("harmonic");
    
    // get start/stop/num frequencies
    myNode->Get( "startFreq", startFreq_ );
    myNode->Get( "stopFreq", stopFreq_ );
    myNode->Get( "numFreq", numFreq_ );

    // read sampling type (optional)
    std::string sampling = "linear";
    myNode->Get( "sampling", sampling, false );
    String2Enum( sampling, samplingType_ );

    // Be verbose
    Info->PrintF( "", "\n HarmonicDriver\n" ); 
    Info->PrintF( "", " -> Sampling strategy ....... '%s'\n",
		  sampling.c_str() );
    Info->PrintF( "", " -> starting frequency ...... '%f'\n", startFreq_ );
    Info->PrintF( "", " -> stopping frequency ...... '%f'\n", stopFreq_  );
    Info->PrintF( "", " -> number of frequencies ... '%d'\n", numFreq_   );
    Info->PrintF( "", "\n" );


    // ---------------------------------
    //  Perform some consistency checks
    // ---------------------------------

    // We cannot have negative frequencies (Schema avoids this, if
    // validated parsing is used)
    if ( startFreq_ < 0.0 || stopFreq_ < 0.0 ) {
      (*error) << "Found negative frequency (startFreq = "
               << startFreq_ << ", stopFreq = " << stopFreq_
               << ") in xml-file";
      Error( __FILE__, __LINE__ );
    }

    // We cannot have negative or zero frequency number (Schema avoids this,
    // if validated parsing is used)
    if ( numFreq_ < 1 ) {
      (*error) << "Found numFreq = " << numFreq_ << " in xml-File! "
               << "This is probably not what you want!";
      Error( __FILE__, __LINE__ );
    }

    // If only one step, check that start and stop are equal
    if ( numFreq_ == 1 ) {
      if ( startFreq_ != stopFreq_ ) {
        (*error) << "Found numFreq = " << numFreq_ << " in xml-File, "
                 << "but startFreq_ = " << startFreq_ << " != "
                 << "stopFreq_ = " << stopFreq_;
      Error( __FILE__, __LINE__ );
      }
    }

    // We do not allow smaller stopping than starting frequencies
    if ( startFreq_ > stopFreq_ ) {
      (*error) << "startFreq = " << startFreq_ << " > stopFreq = "
               << stopFreq_ << " in xml-file!";
      Error( __FILE__, __LINE__ );
    }

    // Check that starting frequency is positive for non-linear sampling
    if ( startFreq_ == 0.0 && samplingType_ != LINEAR_SAMPLING ) {
      (*error) << "You specified sampling = '" << sampling
               << "' which conflicts with startFreq = " << startFreq_;
      Error( __FILE__, __LINE__ );
    }

    // Check for single frequency computation
    if ( startFreq_ == stopFreq_ && numFreq_ != 1 ) {
      if ( numFreq_ > 0 ) {
        (*warning) << "Re-setting numFreq to 1, since startFreq = stopFreq"
                   << " = " << startFreq_;
        Warning( __FILE__, __LINE__ );
        numFreq_ = 1;
      }
      if ( samplingType_ != LINEAR_SAMPLING ) {
        (*warning) << "Re-setting sampling type to 'linear', since "
                   << "startFreq = stopFreq" << " = " << startFreq_;
        Warning( __FILE__, __LINE__ );
        samplingType_ = LINEAR_SAMPLING;
      }
    }

    InitializePDEs();
  }


  // **************
  //   Destructor
  // **************
  HarmonicDriver::~HarmonicDriver() {
  }


  // ****************
  //   SolveProblem
  // ****************
  void HarmonicDriver::SolveProblem() {


    // notify resultHandler about beginning of new sequence step 
    ResultHandler * resHandler = domain->GetResultHandler();
    resHandler->BeginMultiSequenceStep( sequenceStep_,
                                        analysis_,
                                        numFreq_ );

    ptPDE_->WriteGeneralPDEdefines();

    // Perform one simulation for each desired frequency
    for ( actFreqStep_ = 1; actFreqStep_ <= numFreq_; actFreqStep_++ ) {

      // Determine next frequency value
      actFreq_ = ComputeNextFrequency( actFreqStep_ );

      // Set curent frequency value in the mathParser
      domain->GetMathParser()->SetValue( MathParser::GLOB_HANDLER,
                                         "f", actFreq_ );

      // Log info for this frequency
      Info->WriteHarmonicStep( ptPDE_->GetName(), actFreqStep_, actFreq_ );

      // Perform steps for the solution
      ptPDE_->GetSolveStep()->SetActFreq( actFreq_ );
      ptPDE_->GetSolveStep()->SetActStep( actFreqStep_ );
      ptPDE_->GetSolveStep()->PreStepHarmonic();
      ptPDE_->GetSolveStep()->SolveStepHarmonic();
      ptPDE_->GetSolveStep()->PostStepHarmonic();

      // Write results into output-file(s)
      resHandler->BeginStep( actFreqStep_, actFreq_ );
      ptPDE_->WriteResultsInFile( actFreqStep_, actFreq_ );
      resHandler->FinishStep( );

    }

    ptPDE_->Finalize();

    // notify resultHandler about finishing of current sequence step
    if( !isPartOfSequence_ ) {
      resHandler->FinishMultiSequenceStep();
      resHandler->Finalize();
    }
  }


  // ************************
  //   ComputeNextFrequency
  // ************************
  Double HarmonicDriver::ComputeNextFrequency( UInt freqIndex ) {


    Double retFreq = 0.0;

    // Check for single step
    if ( numFreq_ == 1 ) {
      retFreq = startFreq_;
    }
    else {

      switch( samplingType_ ) {

        // Linear sampling
      case LINEAR_SAMPLING:
        retFreq = (freqIndex - 1) * (stopFreq_ - startFreq_) /
          (Double)( numFreq_ - 1 ) + startFreq_;
        break;

        // Logarithmic sampling
      case LOG_SAMPLING:
        {
          Double fac = stopFreq_ / startFreq_;
          fac = std::pow( fac, (Double)(freqIndex - 1) / (numFreq_ - 1) );
          retFreq = startFreq_ * fac;
        }
        break;

        // Reverse logarithmic sampling
      case REVERSE_LOG_SAMPLING:
        {
          Double fac = stopFreq_ / startFreq_;
          fac = std::pow( fac, (Double)(numFreq_ - freqIndex)
                          / (numFreq_ - 1));
          retFreq = stopFreq_ + startFreq_ * ( 1.0 - fac );
        }
        break;

        // Something's wrong
      default:
        std::string damp;
        Enum2String( samplingType_, damp );
        (*error) << "HarmonicDriver::ComputeNextFrequency: '"
                 << damp
                 << "' is not supported as sampling type";
        Error( __FILE__, __LINE__ );
      }
    }

    return retFreq;
  }

}
