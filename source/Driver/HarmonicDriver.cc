#include <fstream>
#include <iostream>
#include <string>

#include <boost/lexical_cast.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>


// signal handling for catching Ctr-C
#include <signal.h>

#include "Driver/HarmonicDriver.hh"
#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Driver/Assemble.hh"
#include "DataInOut/SimState.hh"
#include "Utils/Timer.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ResultHandler.hh"
#include "DataInOut/ProgramOptions.hh"

#include "PDE/StdPDE.hh"

#include "Domain/Domain.hh"


using std::cout;
using std::endl;
namespace pt = boost::posix_time;


// Define pointer to transient driver instance, needed for the signal handler
// to communicate with
HarmonicDriver * instance = NULL;


namespace CoupledField
{
  // ***************
  //   Constructor
  // ***************
  HarmonicDriver::HarmonicDriver( UInt sequenceStep, bool isPartOfSequence,
                                  shared_ptr<SimState> state, Domain* domain,
                                  PtrParamNode paramNode, PtrParamNode infoNode)
    : SingleDriver( sequenceStep, isPartOfSequence, state, domain,
                    paramNode, infoNode ), 
      timer_(new Timer())
  {
    // Set correct analysistype
    analysis_ = BasePDE::HARMONIC;

    param_ = param_->Get("harmonic");

        // replace our info node by a more detailed level
    info_ = info_->Get("harmonic");
    info_->Get(ParamNode::PN_HEADER)->Get("unit")->SetValue("Hz");

    startFreq_ = 0.0;
    stopFreq_ = 0.0;
    numFreq_ = 0;
    actFreqStep_ = 0;
    actFreq_ = 0.0;
    restartStep_ = 0;
    
    isRestarted_ = false;
    
    // Check for presence of restart flag.
    writeRestart_ = true;
    PtrParamNode restartNode = param_->Get("writeRestart", ParamNode::PASS);
    if (restartNode)
      writeRestart_ = restartNode->As<bool>();
    
    // read flag if all results should get written to database file section
    // to allow e.g. for general postprocessing or result extraction
    param_->GetValue("allowPostProc", writeAllSteps_, ParamNode::PASS );
    
  }

  HarmonicDriver::~HarmonicDriver()
  { 
    if( !simState_->HasInput() ) {
      // unregister signal handler and use default action
      // register signal handler
      if( signal( SIGINT, SIG_DFL) == SIG_ERR ) {
        EXCEPTION( "Could not assign default signal action");
      }

      // set global pointer to zero
      instance = NULL;
    }
  }


  void HarmonicDriver::Init(bool restart)
  {
    isRestarted_ = restart;
    
    // We do not know yet if frequencies are given with start, stop, ... or as list
    bool params = ReadParametrizedFrequencies();
    bool list   = ReadFrequencyList();
    
    if(params && list)
      EXCEPTION("'analysis/harmonic' contains 'numFreq/startFreq/stopFreq' and 'frequencyList' concurrently");

    if(!params && !list)
      EXCEPTION("'analysis/harmonic' contains neither 'numFreq/startFreq/stopFreq' nor 'frequencyList' concurrently");

    PtrParamNode in = info_->Get(ParamNode::PN_HEADER);
    in->Get("start")->SetValue(startFreq_);
    in->Get("end")->SetValue(stopFreq_);
    in->Get("numFreq")->SetValue(numFreq_);

    // Initialize first multisequence step, as the method "CheckStoreResults" 
    // relies on the result handler to know already about the current
    // sequencestep. However, in case of optimization, the sequence step
    // gets initialized in Optimization::SolveProblem()
    InitializePDEs();
  }

  bool HarmonicDriver::ReadFrequencyList()
  {
    if(!param_->Has("frequencyList")) return false;

    ParamNodeList& list = param_->Get("frequencyList")->GetChildren();
    freqs.Resize(list.GetSize());
    numFreq_ = freqs.GetSize();
    if(freqs.GetSize() == 0)
      EXCEPTION("cannot have empty frequeny list");

    info_->Get(ParamNode::PN_HEADER)->Get("sampling")->SetValue("frequency list given");

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
    stopFreqStep_ = numFreq_;
    return true;
  }

  bool HarmonicDriver::ReadParametrizedFrequencies()
  {
    // check for existence
    if(!param_->Has("startFreq") || !param_->Has("stopFreq") || !param_->Has("numFreq"))
      return false;

    // get start/stop/num frequencies
    startFreq_ = param_->Get( "startFreq" )->MathParse<Double>();
    stopFreq_ = param_->Get( "stopFreq" )->MathParse<Double>();
    numFreq_ = param_->Get( "numFreq" )->MathParse<UInt>();

    // read sampling type (optional)
    std::string sampling = "linear";
    param_->GetValue( "sampling", sampling, ParamNode::PASS );
    String2Enum( sampling, samplingType_ );

    // store only the sampling strategy
    info_->Get(ParamNode::PN_HEADER)->Get("sampling")->SetValue(sampling);

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
      info_->Get(ParamNode::PN_HEADER)->Get(ParamNode::PN_WARNING)->SetValue("Re-setting numFreq to 1, since startFreq = stopFreq");
      numFreq_ = 1;

      if(samplingType_ != LINEAR_SAMPLING)
      {
        info_->Get(ParamNode::PN_HEADER)->Get(ParamNode::PN_WARNING)->SetValue("Re-setting sampling type to 'linear', since startFreq = stopFreq");
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
    stopFreqStep_ = numFreq_;

    return true; // valid values
  }




  // ****************
  //   SolveProblem
  // ****************
  void HarmonicDriver::SolveProblem()
  {
    // in harmonics one cannot extraxt the result writing to StoreResults() as
    // we have multiple frequencies. (exceptions is optimization)

    ptPDE_->WriteGeneralPDEdefines();
    handler_->BeginMultiSequenceStep( sequenceStep_, analysis_, numFreq_ );
    
    if(writeRestart_ || writeAllSteps_ )
      simState_->BeginMultiSequenceStep( sequenceStep_, analysis_ );
    
    // Read restart information
    ReadRestart();
    numFreq_ = numFreq_ - restartStep_;
    stopFreqStep_ = numFreq_ + restartStep_;
    
    // Perform one simulation for each desired frequency
    for ( actFreqStep_ = restartStep_+1; actFreqStep_ <= numFreq_+restartStep_; actFreqStep_++ )
    {
      // register signal handler only, if it is a child driver
      if( actFreqStep_ > restartStep_+1 && !simState_->HasInput() ) {
        if( signal( SIGINT, HarmonicDriver::SignalHandler) == SIG_ERR ) {
          EXCEPTION( "Could not register Signal Handler");
        }

        // store pointer to global instance variable, if not yet set
        if( !instance ) {
          instance = this;
        }
      }
      
      // Determine next frequency value
      ComputeFrequencyStep(actFreqStep_);

      // Log info for this frequency - suppress in Optimization due to search steps
      if(progOpts->IsQuiet())
        cout << ptPDE_->GetName() << ": Harmonic step " << actFreqStep_ << " frequency " << actFreq_ << endl; 
      else
        cout << endl << ptPDE_->GetName() << ": Harmonic step " 
        << actFreqStep_ <<" ======================= " << endl;

      handler_->BeginStep( actFreqStep_, actFreq_ );
      ptPDE_->WriteResultsInFile( actFreqStep_, actFreq_ );
      handler_->FinishStep( );

      // write out re-start in case of aborted simulation or if all steps should be written
      if(  abortSimulation_  || writeAllSteps_ ) {
        if( writeRestart_ || writeAllSteps_ )
         simState_->WriteStep( actFreqStep_, actFreq_);
      }
      
      // leave loop, if simulation should be aborted
      if ( abortSimulation_ ) {
        break;
      }
        
      // perform runtime estimation
      Double totalTime = timer_->GetWallTime();
      timePerStep_ = totalTime / (Double) actFreqStep_;
      Double remainingTime = (numFreq_ - actFreqStep_) * timePerStep_;
      pt::ptime now = pt::second_clock::local_time();
      now += pt::seconds(static_cast<long int>(remainingTime));
      analysis_id_->Get("timePerStep")->SetValue( timePerStep_ );
      PtrParamNode envNode = info_->GetRoot()->Get(ParamNode::PN_HEADER)->Get("environment");
      envNode->Get("estimatedEnd")->SetValue(pt::to_simple_string( now ));
      envNode->Get("remainingTime")->SetValue(remainingTime);
    } // loop: frequencies

    handler_->FinishMultiSequenceStep();
    if(writeRestart_ || writeAllSteps_ )
      simState_->FinishMultiSequenceStep( !abortSimulation_ );

    // Perform finalization only if not part of sequence
    if(!isPartOfSequence_) 
      handler_->Finalize();
  }

  Double HarmonicDriver::ComputeFrequencyStep(UInt actFreqStep)
  {
    assert(actFreqStep >= 1);
    assert(actFreqStep <= numFreq_+restartStep_);

    actFreqStep_ = actFreqStep;

    // Determine next frequency value from precalculated list
    actFreq_ = freqs[actFreqStep-1].freq; // 1 based!
    assert(freqs[actFreqStep-1].step == actFreqStep);

    analysis_id_ = info_->Get(ParamNode::PN_PROCESS)->Get("step", ParamNode::APPEND);
    analysis_id_->Get("analysis_id")->SetValue(actFreqStep);
    analysis_id_->Get("step")->SetValue(actFreqStep_);
    analysis_id_->Get("value")->SetValue(actFreq_);

    // Set current frequency value in the mathParser
    mathParser_->SetValue( MathParser::GLOB_HANDLER, "f", actFreq_ );
    mathParser_->SetValue( MathParser::GLOB_HANDLER, "step", actFreqStep_ );

    // Perform steps for the solution
    ptPDE_->GetSolveStep()->SetActFreq( actFreq_ );
    ptPDE_->GetSolveStep()->SetActStep( actFreqStep_ );
    ptPDE_->GetSolveStep()->PreStepHarmonic();
    ptPDE_->GetSolveStep()->SolveStepHarmonic(analysis_id_);
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


  void HarmonicDriver::ReadRestart() {

    if ( isRestarted_ ) {

      // Ensure simState is present
      assert( simState_ );

      // Create input reader from current output reader
      bool hasInput = simState_->SetInputReaderToSameOutput();
      if( !hasInput)  {
        EXCEPTION( "Can not perform restarted simulation, as HDF5 file "
            << "contains no restart information.");
      }

      if( simState_->IsCompleted( sequenceStep_ )) {
        std::cout << "\n\n";
        std::cout << "*******************************************************\n";
        std::cout << " No restart necessary, as the desired number of \n";
        std::cout << " frequency steps are already computed. \n";
        std::cout << "*******************************************************\n\n";
        restartStep_ = stopFreqStep_ +1; 
        return;

      } else{

        // Obtain last step
        UInt lastStepNum;
        Double lastStepVal;
        simState_->GetLastStepNum(sequenceStep_, lastStepNum, lastStepVal );
        restartStep_ = lastStepNum;

        // if lastStep is 0, no restart possibility
        if( lastStepNum == 0 ) {
          EXCEPTION( "Can not perform restarted simulation, as HDF5 file "
              << "contains no restart information.");
        }
        std::cout << "\n\n";
        std::cout << "*******************************************************\n";
        std::cout << " Continuing simulation from step " << restartStep_  << std::endl;
        std::cout << "*******************************************************\n";
      }
    }

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

  void HarmonicDriver::SetToStepValue(UInt stepNum, Double stepVal )  {
    // ensure that this method is only called if simState has input
    if( ! simState_->HasInput()) {
      EXCEPTION( "Can only set external time step, if simulation state "
              << "is read from external file" );
    }
    
    actFreqStep_ = stepNum;;
    actFreq_ = stepVal;

    // Set current frequency value in the mathParser
    domain_->GetMathParser()->SetValue( MathParser::GLOB_HANDLER, "f", actFreq_ );
    domain_->GetMathParser()->SetValue( MathParser::GLOB_HANDLER, "step", actFreqStep_ );
  }
  
  
  void HarmonicDriver::SignalHandler( int sig ) {

    if( !instance->abortSimulation_) {
      
      // in addition check, if the current step is the last step, so we
      // do not have to print out any message
      if( instance->actFreqStep_ == instance->stopFreq_ )
        return;
      
      instance->abortSimulation_ = true;
      
      std::cout << "\n\n";
      std::cout << "*******************************************************\n";
      std::cout << " Simulation will be halted after the current time \n";
      std::cout << " step " << instance->actFreqStep_ << " / " << instance->actFreq_
                << " Hz.\n\n";
      std::cout << " Estimated time before end of simulation run: " << 
          int(instance->timePerStep_) << " s" << std::endl;
      std::cout << "*******************************************************\n\n";

    }
  }
}
