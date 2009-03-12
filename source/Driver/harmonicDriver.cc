// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <string>

#include <boost/lexical_cast.hpp>

#include "Driver/harmonicDriver.hh"
#include "Driver/stdSolveStep.hh"
#include "Driver/assemble.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/InfoNode.hh"
#include "DataInOut/WriteInfo.hh"
#include "PDE/StdPDE.hh"
#include "Domain/domain.hh"
#include "DataInOut/resultHandler.hh"


namespace CoupledField
{
  // ***************
  //   Constructor
  // ***************
  HarmonicDriver::HarmonicDriver( UInt sequenceStep, bool isPartOfSequence )
    : SingleDriver( sequenceStep, isPartOfSequence )
  {
    // Set correct analysistype
    analysis_ = BasePDE::HARMONIC;

    // replace our info node by a more detailed level
    infoNode_ = infoNode_->Get("harmonic");

    pn_ = param->Get("sequenceStep", "index", boost::lexical_cast<std::string>(sequenceStep_))
         ->Get("analysis")->Get("harmonic");

    startFreq_ = 0.0;
    stopFreq_ = 0.0;
    numFreq_ = 0;
    actFreqStep_ = 0;
  }

  HarmonicDriver::~HarmonicDriver()
  {
  }


  void HarmonicDriver::Init()
  {
    // We do not know yet if frequencies are given with start, stop, ... or as list
    bool params = ReadParametrizedFrequencies();
    bool list   = ReadFrequencyList();

    if(params && list)
      EXCEPTION("'analysis/harmonic' contains 'numFreq/startFreq/stopFreq' and 'frequencyList' concurrently");

    if(!params && !list)
      EXCEPTION("'analysis/harmonic' contains neither 'numFreq/startFreq/stopFreq' nor 'frequencyList' concurrently");

    InfoNode* in = infoNode_->Get(InfoNode::HEADER);
    in->Get("start", "starting frequency")->SetValue(startFreq_);
    in->Get("end", "stopping frequency")->SetValue(stopFreq_);
    in->Get("numFreq", "number of frequencies")->SetValue(numFreq_);

    InitializePDEs();
  }

  bool HarmonicDriver::ReadFrequencyList()
  {
    if(!pn_->Has("frequencyList")) return false;

    StdVector<ParamNode*>& list = pn_->Get("frequencyList")->GetChildren();
    freqs.Resize(list.GetSize());
    numFreq_ = freqs.GetSize();
    if(freqs.GetSize() == 0)
      EXCEPTION("cannot have empty frequeny list");

    infoNode_->Get(InfoNode::HEADER)->Get("sampling", "sampling strategy")->SetValue("frequency list given");

    for(int fi = 0; fi < (int) list.GetSize(); fi++)
    {
      ParamNode* pn = list[fi];
      assert(pn->GetName() == "freq");
      Frequency& f = freqs[fi];
      f.step = fi+1;
      f.freq = pn->Get("value")->AsDouble();
      f.weight = pn->Get("weight")->AsDouble();

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
    pn_->Get( "startFreq", startFreq_ );
    pn_->Get( "stopFreq", stopFreq_ );
    pn_->Get( "numFreq", numFreq_ );

    // read sampling type (optional)
    std::string sampling = "linear";
    pn_->Get( "sampling", sampling, false );
    String2Enum( sampling, samplingType_ );

    // store only the sampling strategy
    infoNode_->Get(InfoNode::HEADER)->Get("sampling", "sampling strategy")->SetValue(sampling);

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
      infoNode_->Get(InfoNode::HEADER)->Get(InfoNode::WARNING)->SetValue("Re-setting numFreq to 1, since startFreq = stopFreq");
      numFreq_ = 1;

      if(samplingType_ != LINEAR_SAMPLING)
      {
        infoNode_->Get(InfoNode::HEADER)->Get(InfoNode::WARNING)->SetValue("Re-setting sampling type to 'linear', since startFreq = stopFreq");
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
  void HarmonicDriver::SolveProblem(bool write_results, const std::string& comment)
  {
    // in harmonics one cannot extraxt the result writing to StoreResults() as
    // we have multiple frequencies. (exceptions is optimization)

    // be 'silent' and don't do output only for optimization
    if(write_results) {
      handler_->BeginMultiSequenceStep( sequenceStep_, analysis_, numFreq_ );
      // info stuff
      ptPDE_->WriteGeneralPDEdefines();
    }


    // Perform one simulation for each desired frequency
    for ( actFreqStep_ = 1; actFreqStep_ <= numFreq_; actFreqStep_++ )
    {
      // Determine next frequency value
      ComputeFrequencyStep(actFreqStep_);

      // Write results into output-file(s) if we don't do optimization
      if(write_results) {
        // Log info for this frequency - suppress in Optimization due to search steps
        Info->WriteHarmonicStep( ptPDE_->GetName(), actFreqStep_, actFreq_ );
        handler_->BeginStep( actFreqStep_, actFreq_ );
        ptPDE_->WriteResultsInFile( actFreqStep_, actFreq_ );
        handler_->FinishStep( );
      }
    }

    if(write_results) {
      handler_->FinishMultiSequenceStep();
      // notify resultHandler about finishing of current sequence step
      if(!isPartOfSequence_) handler_->Finalize();
    }
  }

  Double HarmonicDriver::ComputeFrequencyStep(UInt actFreqStep, const std::string& comment)
  {
    assert(actFreqStep >= 1);
    assert(actFreqStep <= numFreq_);

    actFreqStep_ = actFreqStep;

    // Determine next frequency value from precalculated list
    actFreq_ = freqs[actFreqStep-1].freq; // 1 based!
    assert(freqs[actFreqStep-1].step == actFreqStep);

    // log this
    // in this current step node also the solvers can write their stuff
    // do it here as ComputeFrequencyStep can be called directly from Optimization!
    currentStepInfoNode_ = infoNode_->Get(InfoNode::PROCESS)->Get("solveStep", InfoNode::APPEND);
    currentStepInfoNode_->Get("currentFrequency")->SetValue(actFreq_);
    currentStepInfoNode_->Get("frequencyStep")->SetValue(actFreqStep_);
    currentStepInfoNode_->Get("systemName")->SetValue(comment);

    // Set curent frequency value in the mathParser
    domain->GetMathParser()->SetValue( MathParser::GLOB_HANDLER, "f", actFreq_ );
    domain->GetMathParser()->SetValue( MathParser::GLOB_HANDLER, "step", actFreqStep_ );

    // Perform steps for the solution
    ptPDE_->GetSolveStep()->SetActFreq( actFreq_ );
    ptPDE_->GetSolveStep()->SetActStep( actFreqStep_ );
    ptPDE_->GetSolveStep()->PreStepHarmonic();
    ptPDE_->GetSolveStep()->SolveStepHarmonic(comment);
    ptPDE_->GetSolveStep()->PostStepHarmonic();

    return actFreq_;
  }


  void HarmonicDriver::StoreResults(double step_val)
  {
    assert(analysis_ == BasePDE::HARMONIC);

    // Write results into output-file(s)
    handler_->BeginStep((unsigned int) step_val, step_val);
    ptPDE_->WriteResultsInFile((unsigned int) step_val, step_val);
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
      }
    }

    return retFreq;
  }

}
