#include <fstream>
#include <iostream>
#include <string>

#include <boost/lexical_cast.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>


// signal handling for catching Ctr-C
#include <signal.h>

#include "Driver/InverseSourceDriver.hh"
#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Driver/Assemble.hh"
#include "DataInOut/SimState.hh"
#include "Utils/Timer.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ResultHandler.hh"
#include "DataInOut/ProgramOptions.hh"

#include "PDE/StdPDE.hh"
#include "FeBasis/BaseFE.hh"
#include "Domain/Domain.hh"


using std::cout;
using std::endl;
namespace pt = boost::posix_time;


// Define pointer to driver instance, needed for the signal handler
// to communicate with
//InverseSourceDriver * instance = NULL;


namespace CoupledField
{
  // ***************
  //   Constructor
  // ***************
  InverseSourceDriver::InverseSourceDriver( UInt sequenceStep, bool isPartOfSequence,
                                  shared_ptr<SimState> state, Domain* domain,
                                  PtrParamNode paramNode, PtrParamNode infoNode)
    : SingleDriver( sequenceStep, isPartOfSequence, state, domain,
                    paramNode, infoNode ), 
      timer_(new Timer())
  {
    // Set correct analysistype
    analysis_ = BasePDE::INVERSESOURCE;

    param_ = param_->Get("inverseSource");

    // replace our info node by a more detailed level
    info_ = info_->Get("inverseSource");
    info_->Get(ParamNode::PN_HEADER)->Get("unit")->SetValue("Hz");
    freq_ = param_->Get( "freq")->MathParse<Double>();

    numFreq_ = 1;
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

  InverseSourceDriver::~InverseSourceDriver()
  { 
    if( !simState_->HasInput() ) {
      // unregister signal handler and use default action
      // register signal handler
      if( signal( SIGINT, SIG_DFL) == SIG_ERR ) {
        EXCEPTION( "Could not assign default signal action");
      }

      // set global pointer to zero
      //instance = NULL;
    }
  }


  void InverseSourceDriver::Init(bool restart)
  {
    isRestarted_ = restart;
    
//    PtrParamNode in = info_->Get(ParamNode::PN_HEADER);
//    in->Get("freq")->SetValue(freq_);
//    std::cout << "Read freq: " << freq_ << std::endl;

    actFreq_ = freq_;

    // Initialize first multisequence step, as the method "CheckStoreResults" 
    // relies on the result handler to know already about the current
    // sequencestep. However, in case of optimization, the sequence step
    // gets initialized in Optimization::SolveProblem()
    InitializePDEs();
  }


  // ****************
  //   SolveProblem
  // ****************
  void InverseSourceDriver::SolveProblem()
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
        if( signal( SIGINT, InverseSourceDriver::SignalHandler) == SIG_ERR ) {
          EXCEPTION( "Could not register Signal Handler");
        }

        // store pointer to global instance variable, if not yet set
//        if( !instance ) {
//          instance = this;
//        }
      }
      
      // Determine next frequency value
      ComputeFrequencyStep(actFreqStep_);

      // Log info for this frequency - suppress in Optimization due to search steps
      if(progOpts->IsQuiet())
        cout << ptPDE_->GetName() << ": InverseSource step " << actFreqStep_ << " frequency " << actFreq_ << endl; 
      else
        cout << endl << ptPDE_->GetName() << ": InverseSource step " 
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

  Double InverseSourceDriver::ComputeFrequencyStep(UInt actFreqStep)
  {
    assert(actFreqStep >= 1);
    assert(actFreqStep <= numFreq_+restartStep_);

    actFreqStep_ = actFreqStep;

    // Determine next frequency value from precalculated list
    actFreq_ = freq_;
    std::cout << "Freq: " <<  actFreq_  << std::endl;

    //get pointers to CoefFncs
    bool isRHSsource = false;
    bool isRHSmeas = false;
    UInt num = 0;
    std::map<SolutionType, shared_ptr<BaseFeFunction> > rhsFeFunctions;
    rhsFeFunctions = ptPDE_->GetRhsFeFunctions();
    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator rFncIt= rhsFeFunctions.begin();
    while ( rFncIt != rhsFeFunctions.end() ) {
     	LoadCoefList loadCoefs = rFncIt->second->GetLoadCoefFunctions();
     	LoadCoefList::iterator it = loadCoefs.begin();
     	for ( ; it != loadCoefs.end(); ++it  ) {
     		PtrCoefFct ptCoef = it->first;
     		CoefFunction::CoefInverseType type = ptCoef->GetInverseType();
     		if ( type == CoefFunction::INVSOURCE ) {
     			rhsSource_ = ptCoef;
     			isRHSsource = true;
     			num++;
     		}
     		else if ( type == CoefFunction::INVMEASURE ) {
     			rhsMeas_ = ptCoef;
     			isRHSmeas = true;
     			num++;
     		}
     	}
     	rFncIt++;
    }
    // check
    if ( num != 2 || isRHSsource == false || isRHSmeas == false ) {
    	std::string message = "You have to define two 'rhsValues', one for the searched";
    	message += " for sources and one for the measurements";
    	EXCEPTION(message);

    }

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

    Double optAmp, optPhase;
    UInt iter = 0;
    bool notConverged = true;

    while ( notConverged && iter < 5 ) {
    	//compute with RHS being the conjugate difference of computed acoustic pressure
    	//and measured pressure in the microphone points
    	rhsSource_->SetActive(false);
    	rhsMeas_->SetActive(true);
    	ptPDE_->GetSolveStep()->PreStepHarmonic();
    	ptPDE_->GetSolveStep()->SolveStepHarmonic(analysis_id_);
    	ptPDE_->GetSolveStep()->PostStepHarmonic();
    	rhsSource_->ComputeOptCondition(optAmp, optPhase);

    	if ( optAmp < 0.1 && optPhase < 0.1 )
    		notConverged = false;

    	std::cout << "\n OptCond, Amp: " << optAmp << "   Phase: " << optPhase << std::endl;

    	//compute with RHS as the current identified sources
    	rhsMeas_->SetActive(false);
      	rhsSource_->SetActive(true);
      	ptPDE_->GetSolveStep()->PreStepHarmonic();
      	ptPDE_->GetSolveStep()->SolveStepHarmonic(analysis_id_);
      	ptPDE_->GetSolveStep()->PostStepHarmonic();
      	rhsSource_->SetActive(false);

    	iter++;
    }


    return actFreq_;
  }


  void InverseSourceDriver::StoreResults(UInt stepNum, double step_val)
  {
    assert(analysis_ == BasePDE::INVERSESOURCE);

    // Write results into output-file(s)
    handler_->BeginStep(stepNum, step_val);
    ptPDE_->WriteResultsInFile(stepNum, step_val);
    handler_->FinishStep( );

  }


  void InverseSourceDriver::ReadRestart() {

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


  void InverseSourceDriver::SetToStepValue(UInt stepNum, Double stepVal )  {
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
  
  
  void InverseSourceDriver::SignalHandler( int sig ) {

//    if( !instance->abortSimulation_) {
//
//      instance->abortSimulation_ = true;
//
//      std::cout << "\n\n";
//      std::cout << "*******************************************************\n";
//      std::cout << " Simulation will be halted after the current frequency\n";
//      std::cout << " step " << instance->actFreqStep_ << " / " << instance->actFreq_
//                << " Hz.\n\n";
//      std::cout << " Estimated time before end of simulation run: " <<
//          int(instance->timePerStep_) << " s" << std::endl;
//      std::cout << "*******************************************************\n\n";
//
//    }

  }
}
