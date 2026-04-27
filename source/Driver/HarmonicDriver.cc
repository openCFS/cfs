#include <fstream>
#include <iostream>
#include <string>

#include <boost/lexical_cast.hpp>


// signal handling for catching Ctr-C
#include <signal.h>

#include "Driver/HarmonicDriver.hh"
#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Driver/Assemble.hh"
#include "DataInOut/SimState.hh"
#include "Utils/Timer.hh"
#include "Utils/mathParser/mathParser.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ResultHandler.hh"
#include "DataInOut/ProgramOptions.hh"

#include "DataInOut/Logging/LogConfigurator.hh"

#include "PDE/StdPDE.hh"

#include "Domain/Domain.hh"


using std::cout;
using std::endl;


// Define pointer to transient driver instance, needed for the signal handler
// to communicate with
HarmonicDriver * instance = NULL;


namespace CoupledField
{
  DEFINE_LOG(harmonicDriver, "harmonicDriver")

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
    info_->Get(ParamNode::HEADER)->Get("unit")->SetValue("Hz");

    startFreq_ = 0.0;
    stopFreq_ = 0.0;
    stopFreqStep_ = 0;
    numFreq_ = 0;
    actFreqStep_ = 0;
    actFreq_ = 0.0;
    restartStep_ = 0;
    timePerStep_ = 0.0;
    samplingType_ = NO_SAMPLING_TYPE;
    
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
        std::cerr << "Could not assign default signal action" << std::endl; // no exceptions in destructors!
        domain->GetInfoRoot()->ToFile();
        exit(-1);
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

    PtrParamNode in = info_->Get(ParamNode::HEADER);
    in->Get("start")->SetValue(startFreq_);
    in->Get("end")->SetValue(stopFreq_);
    in->Get("numFreq")->SetValue(numFreq_);

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

    info_->Get(ParamNode::HEADER)->Get("sampling")->SetValue("frequency list given");

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

  unsigned int HarmonicDriver::GetNumFreq(PtrParamNode node)
  {
    if(node->Has("numFreq"))
      return node->Get("numFreq")->As<unsigned int>();
    else
      return node->Get("frequencyList")->GetChildren().GetSize();
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
    treeStepping_ = false;
    if(param_->Has("sampling")){
      PtrParamNode samplingNode = param_->Get( "sampling" );
      param_->GetValue( "sampling", sampling, ParamNode::PASS );

      treeStepping_ = samplingNode->Get("treeStepping")->As<bool>();
      maxTreeLevel_ = samplingNode->Get("maxTreeLevel")->As<UInt>();
    }
    if (treeStepping_){
      numSteps_ = std::pow(2,maxTreeLevel_) + 1;
    }
    else{
      numSteps_ = numFreq_;
    }
    String2Enum( sampling, samplingType_ );

    // store only the sampling strategy
    info_->Get(ParamNode::HEADER)->Get("sampling")->SetValue(sampling);

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

    // When tree stepping is used, tree level must be high enough to contain all steps at least
    // If not high enough, set such that 1 further refinement level is possible
    if ( treeStepping_ && numFreq_ > numSteps_){
        maxTreeLevel_ = (std::ceil(std::log2(numFreq_ - 1)) + 1);
        numSteps_ = std::pow(2,maxTreeLevel_) + 1;
        WARN( "Number of Frequencies (" << numFreq_ 
           << ") exceed number of steps for defined tree level. Setting max tree level to " << maxTreeLevel_);
    }

    // We do not allow smaller stopping than starting frequencies
    if ( startFreq_ > stopFreq_ ) {
      EXCEPTION( "startFreq = " << startFreq_ << " > stopFreq = "
               << stopFreq_ << " in xml-file!" );
    }

    // Check that starting frequency is positive for non-linear sampling
    if ( startFreq_ == 0.0 && (samplingType_ != LINEAR_SAMPLING) ) {
      EXCEPTION( "You specified sampling = '" << sampling
               << "' which conflicts with startFreq = " << startFreq_ );
    }

    // Check for single frequency computation
    if(startFreq_ == stopFreq_ && numFreq_ > 1)
    {
      info_->Get(ParamNode::HEADER)->SetWarning("Re-setting numFreq to 1, since startFreq = stopFreq");
      numFreq_ = 1;
      numSteps_ = 1;

      if(samplingType_ != LINEAR_SAMPLING)
      {
        info_->Get(ParamNode::HEADER)->SetWarning("Re-setting sampling type to 'linear', since startFreq = stopFreq");
        samplingType_ = LINEAR_SAMPLING;
      }
    }

    // pre calculate the list of frequencies
    freqs.Resize(numSteps_);

    for (unsigned int i = 1; i <= numSteps_; i++)
    {
      // Determine next frequency value
      freqs[i-1].step = i; // one based :(
      freqs[i-1].freq = ComputeNextFrequency(i);
      freqs[i-1].weight = 1.0;
    }
    // Base List
    LOG_DBG(harmonicDriver) << "Base Frequency List:" << std::endl;
    for(auto const & v : freqs)
      LOG_DBG(harmonicDriver) << v << std::endl;
    // Sort list of frequencies to obtain list of steps in case of tree sampling
    if (treeStepping_)
      {
        // Copy list of frequency steps for sorting
        StdVector<Frequency> freqs_tmp(freqs);
        // Sort list regarding frequency
        std::sort(freqs_tmp.begin(),freqs_tmp.end(),[](Frequency const & f0, Frequency const & f1){return (f0.freq < f1.freq);});
        // Sorted list (frequency)
        LOG_DBG(harmonicDriver) << "Sorted List:" << std::endl;
        for(auto const & v : freqs_tmp)
          LOG_DBG(harmonicDriver) << v << std::endl;
        // Set stepping such that increasing step numbers correspond to increasing frequencies
        for (UInt i=0;i<freqs.size();i++)
        {
          auto f_tmp = std::find_if(freqs_tmp.begin(), freqs_tmp.end(), [&](Frequency const & f_tmp) { return f_tmp.freq == freqs[i].freq; });
          UInt idx = std::distance(freqs_tmp.begin(), f_tmp);   
          freqs[i].step = idx+1;
        }
        // Reordered step numbering
        LOG_DBG(harmonicDriver) << "Renumbered Frequency List:" << std::endl;
        for(auto const & v : freqs)
          LOG_DBG(harmonicDriver) << v << std::endl;
      }
    
    stopFreqStep_ = freqs[numFreq_-1].step;

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
    UInt stpIndexStart = 0;
    // Find stpIndex to start in case of restart
    if (isRestarted_)
    {
      auto lastFreq = std::find_if(freqs.begin(), freqs.end(), [&](Frequency const & f) { return f.step == restartStep_; });
      stpIndexStart = std::distance(freqs.begin(), lastFreq) +1;
      numFreq_ = numFreq_ - stpIndexStart;
      LOG_DBG(harmonicDriver) << "Resuming from Frequency Step: " << freqs[stpIndexStart-1] << std::endl;
    }   
    //only used if AMG is set
    ptPDE_->GetSolveStep()->SetAuxMat(false);
    // Perform one simulation for each desired frequency
    for ( UInt stpIndex = stpIndexStart; stpIndex < numFreq_+stpIndexStart; stpIndex++ )
    {
      actFreqStep_ = freqs[stpIndex].step;
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
      ComputeFrequencyStep(freqs[stpIndex]);

      // Log info for this frequency - suppress in Optimization due to search steps
      if(progOpts->IsQuiet())
        cout << ptPDE_->GetName() << ": Harmonic step " << actFreqStep_ << " (" << stpIndex+1-stpIndexStart << "/" << numFreq_ << ")" << " frequency " << actFreq_ << endl; 
      else
        cout << endl << ptPDE_->GetName() << ": Harmonic step " << actFreqStep_ << " (" << stpIndex+1-stpIndexStart << "/" << numFreq_ << ")" <<" ======================= " << endl;

      analysis_id_.step = actFreqStep_;
      analysis_id_.freq = actFreq_;
      // analysis_id_->Get("timePerStep")->SetValue( timePerStep_ );

      handler_->BeginStep( actFreqStep_, actFreq_ );
      ptPDE_->WriteResultsInFile( actFreqStep_, actFreq_ );
      handler_->FinishStep( );

      // write out re-start in case of aborted simulation or if all steps should be written
      if(  actFreqStep_ == stopFreqStep_ || abortSimulation_  || writeAllSteps_ ) {
        if( writeRestart_ || writeAllSteps_ || isPartOfSequence_)
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
      auto time = std::chrono::system_clock::now();
      time += std::chrono::seconds(static_cast<long int>(remainingTime));

      PtrParamNode envNode = info_->GetRoot()->Get(ParamNode::HEADER)->Get("environment");
      envNode->Get("estimatedEnd")->SetValue(Timer::TimeStamp(time));
      envNode->Get("remainingTime")->SetValue(remainingTime);
      envNode->Get("timePerStep")->SetValue(timePerStep_);
    } // loop: frequencies

    handler_->FinishMultiSequenceStep();
    if(writeRestart_ || writeAllSteps_ )
      simState_->FinishMultiSequenceStep( !abortSimulation_ );

    // Perform finalization only if not part of sequence
    if(!isPartOfSequence_) 
      handler_->Finalize();
  }

  Double HarmonicDriver::ComputeFrequencyStep(Frequency const& freqStp)
  {
    assert(freqStp.step >= 1);
    if (!treeStepping_)
      assert(freqStp.step <= stopFreqStep_);

    actFreqStep_ = freqStp.step;
    actFreq_ = freqStp.freq;

    LOG_DBG(harmonicDriver) << "Step: " << actFreqStep_ << ", Frequency: " << actFreq_ << std::endl;

    this->analysis_id_.step = actFreqStep_;
    this->analysis_id_.time = actFreq_;

    // analysis_id_ = info_->Get(ParamNode::PROCESS)->Get("step", ParamNode::APPEND);
    // analysis_id_->Get("analysis_id")->SetValue(actFreqStep);
    // analysis_id_->Get("step")->SetValue(actFreqStep_);
    // analysis_id_->Get("value")->SetValue(actFreq_);

    SetMathParserFreq();

    // Perform steps for the solution
    ptPDE_->GetSolveStep()->SetActFreq( actFreq_ );
    ptPDE_->GetSolveStep()->SetActStep( actFreqStep_ );
    ptPDE_->GetSolveStep()->PreStepHarmonic();
    ptPDE_->GetSolveStep()->SolveStepHarmonic();
    ptPDE_->GetSolveStep()->PostStepHarmonic();

    return actFreq_;
  }

  void HarmonicDriver::SetMathParserFreq(Double freq, UInt freq_step) {
    mathParser_->SetValue( MathParser::GLOB_HANDLER, "f", freq);
    mathParser_->SetValue( MathParser::GLOB_HANDLER, "step", freq_step);
  }
  void HarmonicDriver::SetMathParserFreq() {
    SetMathParserFreq(actFreq_, actFreqStep_);
  }

  unsigned int HarmonicDriver::StoreResults(UInt stepNum, double step_val)
  {
    assert(analysis_ == BasePDE::HARMONIC);

    // Write results into output-file(s)
    handler_->BeginStep(stepNum, step_val);
    ptPDE_->WriteResultsInFile(stepNum, step_val);
    handler_->FinishStep( );

    return stepNum;
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
        restartStep_ = stopFreqStep_; 
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
      // Functions for tree sampling
      // calculates the first distance from freqStart
      auto calc_first_step = [](unsigned int const n)
      {
        return 1.0 / (std::pow(2, n));
      };
      // calculates the every other distance from freqStart
      auto calc_further_step = [&calc_first_step](unsigned int const n, unsigned int const idx)
      {
        return calc_first_step(n) * (1 + 2 * ((idx - 2) - std::pow(2, (n - 1))));
      };
      // calculates step (0.0 to 1.0) based on tree refining (interval splitting)
      auto calc_tree_step = [&calc_first_step,&calc_further_step](unsigned int const idx)
      {
        double step = 0.0;
        if (idx == 1)
          step = 0.0;
        else if (idx == 2)
          step = 1.0;
        else
        {
          int const n = std::floor(std::log2(idx - 2)) + 1;
          if (((int(idx) - 2) % int(std::pow(2, n))) == 0)
            step = calc_first_step(n);
          else
            step = calc_further_step(n, idx);
        }
        return step;
      };


      switch( samplingType_ )
      {
      // Linear sampling
      case LINEAR_SAMPLING:
      {
        if (treeStepping_){
          double const interval = (stopFreq_ - startFreq_);
          double const step = calc_tree_step(freqIndex);
          retFreq = startFreq_ + step*interval;
        }
        else{
          retFreq = (freqIndex - 1) * (stopFreq_ - startFreq_) /  (Double)( numFreq_ - 1 ) + startFreq_;
        }
      }
      break;

      // Logarithmic sampling
      case LOG_SAMPLING:
      {
        if (treeStepping_){
          double const interval = (stopFreq_ - startFreq_);
          double const step = calc_tree_step(freqIndex);
          double const log_start = std::log(startFreq_)/std::log(interval);
          double const log_interval = std::log(stopFreq_/startFreq_)/std::log(interval);
          retFreq = std::pow(interval,log_start + log_interval*step);
        }
        else{
          Double fac = stopFreq_ / startFreq_;
          fac = std::pow( fac, (Double)(freqIndex - 1) / (numFreq_ - 1) );
          retFreq = startFreq_ * fac;
        }
      }
      break;

      // Reverse logarithmic sampling
      case REVERSE_LOG_SAMPLING:
      {
        if (treeStepping_){
          EXCEPTION("Tree stepping for reverse log sampling not implemented yet.");
        }
        else{
          Double fac = stopFreq_ / startFreq_;
          fac = std::pow( fac, (Double)(numFreq_ - freqIndex) / (numFreq_ - 1));
          retFreq = stopFreq_ + startFreq_ * ( 1.0 - fac );
        }
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
      std::cout << " Simulation will be halted after the current frequency\n";
      std::cout << " step " << instance->actFreqStep_ << " / " << instance->actFreq_
                << " Hz.\n\n";
      std::cout << " Estimated time before end of simulation run: " << 
          int(instance->timePerStep_) << " s" << std::endl;
      std::cout << "*******************************************************\n\n";

    }
  }
}
