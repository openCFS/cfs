#include <fstream>
#include <iostream>
#include <string>

#include "harmonicDriver.hh"
#include "stdSolveStep.hh"
#include "assemble.hh"

#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/CommandLine/BaseCommandLineHandler.hh"
#include "DataInOut/WriteInfo.hh"
#include "PDE/StdPDE.hh"
#include "Domain/domain.hh"


namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
  HarmonicDriver::HarmonicDriver( Domain * adomain, 
                                  UInt stepOffset,
                                  Double timeOffset,
                                  std::string driverTag,
                                  Boolean isPartOfSequence )

    : SingleDriver( adomain, stepOffset, timeOffset, driverTag,
                    isPartOfSequence ) {

    ENTER_FCN( " HarmonicDriver::HarmonicDriver" );

    // vectors for accessing parameters
    StdVector<std::string> keyVec, attrVec, valVec;

    attrVec = "tag";
    valVec = driverTag_;

    // --------------------------------------------------------
    // Get frequency sampling information from parameter object
    // --------------------------------------------------------

    // In order to avoid problems with non-matching tag attributes,
    // we do a GetList() first, since this does not check for defaults
    StdVector<Double> auxVec;
    keyVec = "harmonic", "startFreq";
    params->GetList( keyVec, attrVec, valVec, auxVec );
    if ( auxVec.GetSize() == 1 ) {
      startFreq_ = auxVec[0];
    }
    else {
      (*error) << "HarmonicDriver::HarmonicDriver:\n "
               << "Could not find a section <harmonic> with tag = '"
               << driverTag_ << "' in " << commandLine->GetParamFile();
      Error( __FILE__, __LINE__ );
    }

    keyVec = "harmonic", "stopFreq";
    params->Get( keyVec, attrVec, valVec, stopFreq_ );

    keyVec = "harmonic", "numFreq";
    params->Get( keyVec, attrVec, valVec, numFreq_ );

    std::string damp;
    keyVec = "harmonic", "adjustDamping";
    params->Get( keyVec, attrVec, valVec, damp );
    adjustDamping_ = damp == "yes" ? TRUE : FALSE;

    std::string sampling;
    keyVec = "harmonic", "sampling";
    params->Get( keyVec, attrVec, valVec, sampling );
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
  }


  // **************
  //   Destructor
  // **************
  HarmonicDriver::~HarmonicDriver() {
    ENTER_FCN( " HarmonicDriver::~HarmonicDriver" );
  }


  // ****************
  //   SolveProblem
  // ****************
  void HarmonicDriver::SolveProblem() {

    ENTER_FCN( " HarmonicDriver::SolveProblem" );

    Boolean reset = TRUE;

    // There are some special things to do in the case that we are not
    // part of a multi sequence analysis
    if ( isPartOfSequence_ == FALSE ) {

      ptdomain_->PrintGrid();

      // Get list of PDEs which have to be solved and intialize them
      GetMyPDEs();
      Info->StartProgress ( "Starting to solve problem", FALSE );
    }

    ptPDE_->WriteGeneralPDEdefines();

    // Perform one simulation for each desired frequency
    for ( UInt fstep = 1; fstep <= numFreq_; fstep++ ) {

      // Determine next frequency value
      actFreq_ = ComputeNextFrequency( fstep );

      // Set curent frequency value in the mathParser
      domain->GetMathParser()->SetValue( MathParser::GLOB_HANDLER,
                                         "f", actFreq_ );

      // Log info for this frequency
      Info->WriteHarmonicStep( ptPDE_->GetName(), fstep, actFreq_ );

      // Perform steps for the solution
      ptPDE_->GetSolveStep()->SetActFreq( actFreq_ );
      ptPDE_->GetSolveStep()->SetActStep( fstep );
      ptPDE_->GetSolveStep()->PreStepHarmonic( reset );
      ptPDE_->GetSolveStep()->SolveStepHarmonic( reset );
      ptPDE_->GetSolveStep()->PostStepHarmonic( reset );

      // Write results into output-file
      ptPDE_->PostProcess();
      ptPDE_->WriteResultsInFile(fstep, actFreq_);

      //write history data
      ptPDE_->WriteHistoryInFile(fstep, actFreq_);
    }
  }


  // ************************
  //   ComputeNextFrequency
  // ************************
  Double HarmonicDriver::ComputeNextFrequency( UInt freqIndex ) {

    ENTER_FCN( "HarmonicDriver::ComputeNextFrequency" );

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
