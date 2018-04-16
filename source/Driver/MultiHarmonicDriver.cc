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


#include <boost/smart_ptr/shared_ptr.hpp>
#include <DataInOut/ParamHandling/ParamNode.hh>
#include <DataInOut/ResultHandler.hh>
#include <DataInOut/SimState.hh>
#include <Domain/Domain.hh>
#include <Driver/MultiHarmonicDriver.hh>
#include <Driver/SolveSteps/BaseSolveStep.hh>
#include <General/defs.hh>
#include <General/Environment.hh>
#include <General/Exception.hh>
#include <PDE/BasePDE.hh>
#include <signal.h>
#include <Utils/StdVector.hh>
#include <Utils/Timer.hh>
#include <cstdlib>
#include <iostream>


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
    
    numFFT_ = param_->Get("numFFTPoints")->MathParse<Double>();

    // read flag if all results should get written to database file section
    // to allow e.g. for general postprocessing or result extraction
    param_->GetValue("allowPostProc", writeAllSteps_, ParamNode::PASS );

    actHarm_ = 0;

    harmFreq_.Resize(2 * numHarmonics_N_ + 1);
    // static harmonic (zero frequency) at position numHarmonics_N_
    harmFreq_[numHarmonics_N_] = 0.0;
    Double dummy = -1;
    for(UInt i = 0; i < numHarmonics_N_; ++i  ){
      // negative frequencies at the beginning
      harmFreq_[numHarmonics_N_ - 1 - i] = dummy * (i+1) * baseFreq_;
      // positive frequencies starting at numHarmonics_N_ + 1
      harmFreq_[numHarmonics_N_ + 1 + i] = (i + 1) * baseFreq_;
    }

    // Set the flag for the harmonic callback mechanism
    mathParser_->SetValue( MathParser::GLOB_HANDLER, "mhFlag", 0 );
    // Set the flag for the solution-cache callback mechanism
    mathParser_->SetValue( MathParser::GLOB_HANDLER, "updateTrigger", false );

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

    ptPDE_->WriteGeneralPDEdefines();
    
    UInt numFreq = 2 * numHarmonics_N_ + 1;
    handler_->BeginMultiSequenceStep( sequenceStep_, analysis_, numFreq );

    if( writeAllSteps_ )
      simState_->BeginMultiSequenceStep( sequenceStep_, analysis_ );
    

    // In multiharmonic analysis we speak in terms of multiples of base-harmonics.
    // The system matrices of the single harmonics get inserted into the global matrix
    // and then the system is solved.

    // Therefore we don't consider one frequency isolated from the others
    // because they are coupled



    //===========================================================================================
    //        WHAT ARE WE DOING WITH THE MATHPARSER EXPRESSIONS???
    //===========================================================================================
    // Set current frequency value in the mathParser
    //mathParser_->SetValue( MathParser::GLOB_HANDLER, "f", actFreq_ );
    //mathParser_->SetValue( MathParser::GLOB_HANDLER, "step", actFreqStep_ );

    // Perform steps for the solution
    ptPDE_->GetSolveStep()->SetMultHarmonicFreq( harmFreq_ );

    //===========================================================================================
    //        CONTINUE HERE
    //===========================================================================================
    // next step is to set a frequency loop, probably pass harmFreq_to
    // SolveStepHarmonic(), since it also has a poiter to mathparser

    ptPDE_->GetSolveStep()->PreStepHarmonic();

    ptPDE_->GetSolveStep()->SolveStepHarmonic();
    ptPDE_->GetSolveStep()->PostStepHarmonic();



    for(UInt i = 0; i < harmFreq_.GetSize(); ++i  ){
      // current frequency
      Double actFreq = harmFreq_[i];
      // which harmonic are we considering
      Integer h = -numHarmonics_N_ + i;

      // the harmonics number h can be negative but the frequency step can't
      // therefore use the index i
      // i = [  0     1     2  ...   N    N+1   N+2 ...  2N ]
      // h = [ -N   -N+1  -N+2 ...   0     1     2  ...   N ]
      UInt actFreqStep = i;
      analysis_id_.step = actFreqStep;
      analysis_id_.freq = actFreq;

      // Now we need to set the actual solution and rhs vector (depending on
      // the current harmonic). This basically stores the solution and rhs back to
      // the PDE itself
      // Usually this is done in the StdSolveStep when algsys_->GetSolutionVal(solVec_)
      // is called but in the multiharmonic case this must be done here.
      //TODO do we even need this ?!

      ptPDE_->GetSolveStep()->GetRHSValMultHarm(i);
      ptPDE_->GetSolveStep()->GetSolutionValMultHarm(i);

      handler_->BeginStep( actFreqStep, actFreq );
      ptPDE_->WriteResultsInFile( actFreqStep, actFreq );
      handler_->FinishStep( );
    }



/*
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
*/


    handler_->FinishMultiSequenceStep();
    if( writeAllSteps_ )
      simState_->FinishMultiSequenceStep( !abortSimulation_ );

    // Perform finalization only if not part of sequence
    if(!isPartOfSequence_) 
      handler_->Finalize();

  }



  /*
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
