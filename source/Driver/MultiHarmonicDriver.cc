// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     MultiHarmonicDriver.cc
 *       \brief    Driver for nonlinear frequency domain analysis
 *
 *       \date     Jan 2, 2017
 *       \author   kroppert
 */
//================================================================================================


#include <fstream>
#include <iostream>
#include <string>

#include <boost/lexical_cast.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>


// signal handling for catching Ctr-C
#include <signal.h>

#include "Driver/MultiHarmonicDriver.hh"
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



namespace CoupledField
{
  // ***************
  //   Constructor
  // ***************
  MultiHarmonicDriver::MultiHarmonicDriver( UInt sequenceStep, bool isPartOfSequence,
                                  shared_ptr<SimState> state, Domain* domain,
                                  PtrParamNode paramNode, PtrParamNode infoNode)
    : SingleDriver( sequenceStep, isPartOfSequence, state, domain,
                    paramNode, infoNode ), 
      timer_(new Timer())
  {
    // Set correct analysistype
    analysis_ = BasePDE::MULTIHARMONIC;

    param_ = param_->Get("multiharmonic");

    // replace our info node by a more detailed level
    info_ = info_->Get("multiharmonic");
    info_->Get(ParamNode::HEADER)->Get("unit")->SetValue("Hz");

    baseFreq_ = param_->Get("baseFreq")->MathParse<Double>();

    numHarmonics_N_ = param_->Get("numHarmonics_N")->MathParse<Double>();
    
    // default for numHarmonics_M is numHarmonics_N
    numHarmonics_M_ = ( param_->Has("numHarmonics_M") ) ?  param_->Get("numHarmonics_M")->MathParse<Double>() : numHarmonics_N_;
    
    // read flag if all results should get written to database file section
    // to allow e.g. for general postprocessing or result extraction
    param_->GetValue("allowPostProc", writeAllSteps_, ParamNode::PASS );
    
  }

  MultiHarmonicDriver::~MultiHarmonicDriver()
  { 
    if( !simState_->HasInput() ) {
      // unregister signal handler and use default action
      // register signal handler
      if( signal( SIGINT, SIG_DFL) == SIG_ERR ) {
        std::cerr << "Could not assign default signal action" << std::endl; // no exceptions in destructors!
        domain->GetInfoRoot()->ToFile();
        exit(-1);
      }

    }
  }


  void MultiHarmonicDriver::Init(bool restart)
  {
	  InitializePDEs();
  }


  // ****************
  //   SolveProblem
  // ****************
  void MultiHarmonicDriver::SolveProblem()
  {
/*

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
    
    //only used if AMG is set
    ptPDE_->GetSolveStep()->SetAuxMat(false);
    // Perform one simulation for each desired frequency
    for ( actFreqStep_ = restartStep_+1; actFreqStep_ <= numFreq_+restartStep_; actFreqStep_++ )
    {
      // register signal handler only, if it is a child driver
      if( actFreqStep_ > restartStep_+1 && !simState_->HasInput() ) {
        if( signal( SIGINT, MultiHarmonicDriver::SignalHandler) == SIG_ERR ) {
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
        cout << endl << ptPDE_->GetName() << ": Harmonic step " << actFreqStep_ <<" ======================= " << endl;

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
      pt::ptime now = pt::second_clock::local_time();
      now += pt::seconds(static_cast<long int>(remainingTime));

      PtrParamNode envNode = info_->GetRoot()->Get(ParamNode::HEADER)->Get("environment");
      envNode->Get("estimatedEnd")->SetValue(pt::to_simple_string( now ));
      envNode->Get("remainingTime")->SetValue(remainingTime);
      envNode->Get("timePerStep")->SetValue(timePerStep_);
    } // loop: frequencies

    handler_->FinishMultiSequenceStep();
    if(writeRestart_ || writeAllSteps_ )
      simState_->FinishMultiSequenceStep( !abortSimulation_ );

    // Perform finalization only if not part of sequence
    if(!isPartOfSequence_) 
      handler_->Finalize();
*/
  }
/*

  Double MultiHarmonicDriver::ComputeFrequencyStep(UInt actFreqStep)
  {

	  assert(actFreqStep >= 1);
    assert(actFreqStep <= numFreq_+restartStep_);

    actFreqStep_ = actFreqStep;

    // Determine next frequency value from precalculated list
    actFreq_ = freqs[actFreqStep-1].freq; // 1 based!
    assert(freqs[actFreqStep-1].step == actFreqStep);

    this->analysis_id_.step = actFreqStep;
    this->analysis_id_.time = actFreq_;

    // analysis_id_ = info_->Get(ParamNode::PROCESS)->Get("step", ParamNode::APPEND);
    // analysis_id_->Get("analysis_id")->SetValue(actFreqStep);
    // analysis_id_->Get("step")->SetValue(actFreqStep_);
    // analysis_id_->Get("value")->SetValue(actFreq_);

    // Set current frequency value in the mathParser
    mathParser_->SetValue( MathParser::GLOB_HANDLER, "f", actFreq_ );
    mathParser_->SetValue( MathParser::GLOB_HANDLER, "step", actFreqStep_ );

    // Perform steps for the solution
    ptPDE_->GetSolveStep()->SetActFreq( actFreq_ );
    ptPDE_->GetSolveStep()->SetActStep( actFreqStep_ );
    ptPDE_->GetSolveStep()->PreStepHarmonic();
    ptPDE_->GetSolveStep()->SolveStepHarmonic();
    ptPDE_->GetSolveStep()->PostStepHarmonic();

    return actFreq_;
  }



  void MultiHarmonicDriver::StoreResults(UInt stepNum, double step_val)
  {
    assert(analysis_ == BasePDE::HARMONIC);

    // Write results into output-file(s)
    handler_->BeginStep(stepNum, step_val);
    ptPDE_->WriteResultsInFile(stepNum, step_val);
    handler_->FinishStep( );

  }
*/

  
}
