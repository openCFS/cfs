// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <assert.h>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <string>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/programOptions.hh"
#include "DataInOut/resultHandler.hh"
#include "Domain/domain.hh"
#include "Driver/baseSolveStep.hh"
#include "Driver/harmonicDriver.hh"
#include "Driver/singleDriver.hh"
#include "General/exception.hh"
#include "PDE/basePDE.hh"
#include "Utils/Timer.hh"
#include "Utils/mathParser/mathParser.hh"
#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/lexical_cast.hpp"

namespace CoupledField {
class AdjointParameters;
}  // namespace CoupledField


using std::cout;
using std::endl;
namespace pt = boost::posix_time;

namespace CoupledField
{
  // ***************
  //   Constructor
  // ***************
  HarmonicDriver::HarmonicDriver( UInt sequenceStep, bool isPartOfSequence )
    : SingleDriver( sequenceStep, isPartOfSequence ), timer_(new Timer())
  {
    // Set correct analysistype
    analysis_ = BasePDE::HARMONIC;

    // replace our info node by a more detailed level
    driverNode = driverNode->Get("harmonic");
    driverNode->Get("sequenceStep")->SetValue(sequenceStep);
    driverNode->Get(ParamNode::HEADER)->Get("unit")->SetValue("Hz");

    pn_ = param->GetByVal("sequenceStep", "index", boost::lexical_cast<std::string>(sequenceStep_))
         ->Get("analysis")->Get("harmonic");

    startFreq_ = 0.0;
    stopFreq_ = 0.0;
    numFreq_ = 0;
    actFreq_ = 0.0;
    actFreqStep_ = 0;
    samplingType_ = NO_SAMPLING_TYPE;
  }

  HarmonicDriver::~HarmonicDriver()
  {  }


  void HarmonicDriver::Init()
  {
    // We do not know yet if frequencies are given with start, stop, ... or as list
    bool params = ReadParametrizedFrequencies();
    bool list   = ReadFrequencyList();

    if(params && list)
      EXCEPTION("'analysis/harmonic' contains 'numFreq/startFreq/stopFreq' and 'frequencyList' concurrently");

    if(!params && !list)
      EXCEPTION("'analysis/harmonic' contains neither 'numFreq/startFreq/stopFreq' nor 'frequencyList'");

    PtrParamNode in = driverNode->Get(ParamNode::HEADER);
    in->Get("start")->SetValue(startFreq_);
    in->Get("end")->SetValue(stopFreq_);
    in->Get("numFreq")->SetValue(numFreq_);

    // Initialize first multisequence step, as the method "CheckStoreResults" 
    // relies on the result handler to know already about the current
    // sequencestep. However, in case of optimization, the sequence step
    // gets initialized in Optimization::SolveProblem()
    if( !domain->GetOptimization()) {
      handler_->BeginMultiSequenceStep( sequenceStep_, analysis_, numFreq_ );
    }
    InitializePDEs();
  }

  bool HarmonicDriver::ReadFrequencyList()
  {
    if(!pn_->Has("frequencyList")) return false;

    ParamNodeList& list = pn_->Get("frequencyList")->GetChildren();
    freqs.Resize(list.GetSize());
    numFreq_ = freqs.GetSize();
    if(freqs.GetSize() == 0)
      EXCEPTION("cannot have empty frequency list");

    driverNode->Get(ParamNode::HEADER)->Get("sampling")->SetValue("frequency list given");

    for(int fi = 0; fi < (int) list.GetSize(); fi++)
    {
      PtrParamNode pn = list[fi];
      assert(pn->GetName() == "freq");
      Frequency& f = freqs[fi];
      f.step = fi+1;
      f.freq = pn->Get("value")->MathParse<Double>();
      f.weight = pn->Get("weight")->MathParse<Double>();

      // set bounds (we keep unsorted)
      startFreq_ = std::min(startFreq_, f.freq);
      stopFreq_  = std::max(stopFreq_, f.freq);

      // plausibility checking
      if(f.freq < 0.0)   EXCEPTION("The " << (fi+1) << ". frequencyList entry is negative: " << f.freq);
      if(f.weight < 0.0) EXCEPTION("The " << (fi+1) << ". frequencyList entry weight is negative: " << f.weight);
      for(int o = 0; o < fi-1; o++)
        if(freqs[o].freq == f.freq)
          EXCEPTION("Multiple occurence of in frequencyList: f=" << f.freq << " at position " << (o+1) << " and " << (fi+1));
    }
    return true;
  }

  bool HarmonicDriver::ReadParametrizedFrequencies()
  {
    // check for existence
    if(!pn_->Has("startFreq") || !pn_->Has("stopFreq") || !pn_->Has("numFreq"))
      return false;

    // get start/stop/num frequencies
    startFreq_ = pn_->Get( "startFreq" )->MathParse<Double>();
    stopFreq_ = pn_->Get( "stopFreq" )->MathParse<Double>();
    numFreq_ = pn_->Get( "numFreq" )->MathParse<UInt>();

    // read sampling type (optional)
    std::string sampling = "linear";
    pn_->GetValue( "sampling", sampling, ParamNode::PASS );
    String2Enum( sampling, samplingType_ );

    // store only the sampling strategy
    driverNode->Get(ParamNode::HEADER)->Get("sampling")->SetValue(sampling);

    // ---------------------------------
    //  Perform some consistency checks
    // ---------------------------------

    // We cannot have negative frequencies (Schema avoids this, if
    // validated parsing is used)
    if ( startFreq_ < 0.0 || stopFreq_ < 0.0 ) {
      EXCEPTION( "Found negative frequency (startFreq = "
               << startFreq_ << ", stopFreq = " << stopFreq_
               << ") in xml-file" );
    }

    // We cannot have negative or zero frequency number (Schema avoids this,
    // if validated parsing is used)
    if ( numFreq_ < 1 ) {
      EXCEPTION( "Found numFreq = " << numFreq_ << " in xml-File! "
               << "This is probably not what you want!" );
    }

    // If only one step, check that start and stop are equal
    if ( numFreq_ == 1 && startFreq_ != stopFreq_) {
      if ( startFreq_ != stopFreq_ ) {
        EXCEPTION( "Found numFreq = " << numFreq_ << " in xml-File, "
                 << "but startFreq_ = " << startFreq_ << " != "
                 << "stopFreq_ = " << stopFreq_ );
      }
    }

    // We do not allow smaller stopping than starting frequencies
    if ( startFreq_ > stopFreq_ ) {
      EXCEPTION( "startFreq = " << startFreq_ << " > stopFreq = "
               << stopFreq_ << " in xml-file!" );
    }

    // Check that starting frequency is positive for non-linear sampling
    if ( startFreq_ == 0.0 && samplingType_ != LINEAR_SAMPLING ) {
      EXCEPTION( "You specified sampling = '" << sampling
               << "' which conflicts with startFreq = " << startFreq_ );
    }

    // Check for single frequency computation
    if(startFreq_ == stopFreq_ && numFreq_ > 1)
    {
      driverNode->Get(ParamNode::HEADER)->Get(ParamNode::WARNING)->SetValue("Re-setting numFreq to 1, since startFreq = stopFreq");
      numFreq_ = 1;

      if(samplingType_ != LINEAR_SAMPLING)
      {
        driverNode->Get(ParamNode::HEADER)->Get(ParamNode::WARNING)->SetValue("Re-setting sampling type to 'linear', since startFreq = stopFreq");
        samplingType_ = LINEAR_SAMPLING;
      }
    }

    // pre calculate the list of frequencies
    freqs.Resize(numFreq_);

    for (unsigned int i = 1; i <= numFreq_; i++)
    {
      // Determine next frequency value
      freqs[i-1].step = i; // one based :(
      freqs[i-1].freq = ComputeNextFrequency(i);
      freqs[i-1].weight = 1.0;
    }

    return true; // valid values
  }




  // ****************
  //   SolveProblem
  // ****************
  void HarmonicDriver::SolveProblem(bool write_results, PtrParamNode analysis_id, AdjointParameters* adjointParams)
  {
    // in harmonics one cannot extraxt the result writing to StoreResults() as
    // we have multiple frequencies. (exceptions is optimization)

    // be 'silent' and don't do output only for optimization
    if(write_results) {
      // info stuff
      ptPDE_->WriteGeneralPDEdefines();
    }


    // Perform one simulation for each desired frequency
    for ( actFreqStep_ = 1; actFreqStep_ <= numFreq_; actFreqStep_++ )
    {
      // Determine next frequency value
      ComputeFrequencyStep(actFreqStep_, analysis_id, adjointParams);

      // Write results into output-file(s) if we don't do optimization
      if(write_results) {
        // Log info for this frequency - suppress in Optimization due to search steps
        if(progOpts->IsQuiet())
          cout << ptPDE_->GetName() << ": Harmonic step " << actFreqStep_ << " frequency " << actFreq_ << endl; 
        else
          cout << endl << ptPDE_->GetName() << ": Harmonic step " 
               << actFreqStep_ <<" ======================= " << endl;

        handler_->BeginStep( actFreqStep_, actFreq_ );
        ptPDE_->WriteResultsInFile( actFreqStep_, actFreq_ );
        handler_->FinishStep( );
      }
      
      // perform runtime estimation
      Double totalTime = timer_->GetWallTime();
      Double timePerStep = totalTime / (Double) actFreqStep_;
      Double remainingTime = (numFreq_ - actFreqStep_) * timePerStep;
      pt::ptime now = pt::second_clock::local_time();
      now += pt::seconds(static_cast<long int>(remainingTime));
      analysis_id_->Get("timePerStep")->SetValue( timePerStep );
      PtrParamNode envNode = info->Get(ParamNode::HEADER)->Get("environment");
      envNode->Get("estimatedEnd")->SetValue(pt::to_simple_string( now ));
      envNode->Get("remainingTime")->SetValue(remainingTime);
    }

    if(write_results) {
      handler_->FinishMultiSequenceStep();
      // notify resultHandler about finishing of current sequence step
      if(!isPartOfSequence_) handler_->Finalize();
    }
  }

  Double HarmonicDriver::ComputeFrequencyStep(UInt actFreqStep, PtrParamNode given_analysis_id, AdjointParameters* adjointParams)
  {
    assert(actFreqStep >= 1);
    assert(actFreqStep <= numFreq_);

    actFreqStep_ = actFreqStep;

    // Determine next frequency value from precalculated list
    actFreq_ = freqs[actFreqStep-1].freq; // 1 based!
    assert(freqs[actFreqStep-1].step == actFreqStep);

    if(given_analysis_id == NULL)
    {
      analysis_id_ = driverNode->Get(ParamNode::PROCESS)->Get("step", ParamNode::APPEND);
      analysis_id_->Get("analysis_id")->SetValue(actFreqStep);
    }
    else
    {
      analysis_id_ = given_analysis_id;
      assert(analysis_id_->Has("analysis_id"));
    }
    analysis_id_->Get("step")->SetValue(actFreqStep_);
    analysis_id_->Get("value")->SetValue(actFreq_);

    // Set curent frequency value in the mathParser
    mathParser_->SetValue( MathParser::GLOB_HANDLER, "f", actFreq_ );
    mathParser_->SetValue( MathParser::GLOB_HANDLER, "step", actFreqStep_ );

    // Perform steps for the solution
    ptPDE_->GetSolveStep()->SetActFreq( actFreq_ );
    ptPDE_->GetSolveStep()->SetActStep( actFreqStep_ );
    ptPDE_->GetSolveStep()->PreStepHarmonic();
    ptPDE_->GetSolveStep()->SolveStepHarmonic(analysis_id_, adjointParams);
    ptPDE_->GetSolveStep()->PostStepHarmonic();

    return actFreq_;
  }


  void HarmonicDriver::StoreResults(UInt stepNum, double step_val)
  {
    assert(analysis_ == BasePDE::HARMONIC);

    // Write results into output-file(s)
    handler_->BeginStep(stepNum, step_val);
    ptPDE_->WriteResultsInFile(stepNum, step_val);
    handler_->FinishStep( );
  }



  // ************************
  //   ComputeNextFrequency
  // ************************
  Double HarmonicDriver::ComputeNextFrequency( UInt freqIndex ) const
  {
    Double retFreq = -1.0;

    // Check for single step
    if ( numFreq_ == 1 )
    {
      retFreq = startFreq_;
    }
    else
    {
      switch( samplingType_ )
      {
      // Linear sampling
      case LINEAR_SAMPLING:
        retFreq = (freqIndex - 1) * (stopFreq_ - startFreq_) /  (Double)( numFreq_ - 1 ) + startFreq_;
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
        fac = std::pow( fac, (Double)(numFreq_ - freqIndex) / (numFreq_ - 1));
        retFreq = stopFreq_ + startFreq_ * ( 1.0 - fac );
      }
      break;

      // Something's wrong
      default:
        std::string damp;
        Enum2String( samplingType_, damp );
        EXCEPTION( "HarmonicDriver::ComputeNextFrequency: '"
                 << damp
                 << "' is not supported as sampling type" );
        break;
      }
    }

    return retFreq;
  }

}
