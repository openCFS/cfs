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
#include "Utils/mathParser/mathParser.hh"
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
    if(numHarmonics_N_ %2 == 0){
      EXCEPTION("Please provide an odd number of harmonics N")
    }
    
    // default for numHarmonics_M is numHarmonics_N
    numHarmonics_M_ = ( param_->Has("numHarmonics_M") ) ?  param_->Get("numHarmonics_M")->MathParse<Double>() : numHarmonics_N_;
    if(numHarmonics_M_ > numHarmonics_N_){
      EXCEPTION("Number of harmonics M is larger than number of harmonics N !!")
    }

    numFFT_ = param_->Get("numFFTPoints")->MathParse<Double>();

    // read flag if all results should get written to database file section
    // to allow e.g. for general postprocessing or result extraction
    param_->GetValue("allowPostProc", writeAllSteps_, ParamNode::PASS );

    actHarm_ = 0;

    // Do we need to incorporate the zero harmonic?
    param_->GetValue("fullSystem", fullSystem_, ParamNode::PASS );


    if(fullSystem_){
      harmFreq_.Resize(2*numHarmonics_N_ + 1 );
    }else{
      harmFreq_.Resize(numHarmonics_N_ + 1 );
    }
    harmFreq_.Init(0.0);
    for(UInt i = 0; i < harmFreq_.GetSize(); ++i  ){
      harmFreq_[i] = HarmonicOfIndex(i) * baseFreq_;
    }

    numFreq_ = harmFreq_.GetSize();

    // Set the flag for the harmonic callback mechanism
    mathParser_->SetValue( MathParser::GLOB_HANDLER, "finishCash", 0 );
    mathParser_->SetValue( MathParser::GLOB_HANDLER, "harmonicHandle", -1 );

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



  UInt MultiHarmonicDriver::IndexOfHarmonic(const Integer& harmonic){
    UInt ind;
    if(fullSystem_){
      ind = harmonic + numHarmonics_N_;
    }else{
      ind = (numHarmonics_N_ + harmonic)/2;
    }
    return ind;
  }

  Integer MultiHarmonicDriver::HarmonicOfIndex(const UInt& ind){
    Integer harmonic;
    if(fullSystem_){
      harmonic = -numHarmonics_N_ + ind;
    }else{
      harmonic = -numHarmonics_N_ + 2*ind;
    }
    return harmonic;
  }


  void MultiHarmonicDriver::SetToStepValue(UInt stepNum, Double stepVal ){
    // ensure that this method is only called if simState has input
    if( ! simState_->HasInput()) {
      EXCEPTION( "Can only set external time step, if simulation state "
              << "is read from external file" );
    }

    // current frequency
    Double actFreq = harmFreq_[stepNum];

    UInt actFreqStep = stepNum;
    analysis_id_.step = actFreqStep;
    analysis_id_.freq = actFreq;

    Integer h = HarmonicOfIndex(stepNum);
    // We need to activate the correct harmonic results in CoefFunctionHarmBalance
    mathParser_->SetValue(MathParser::GLOB_HANDLER, "f", actFreq);
    mathParser_->SetValue(MathParser::GLOB_HANDLER, "harmonicHandle", h);
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
    
    handler_->BeginMultiSequenceStep( sequenceStep_, analysis_, numFreq_ );

    //if( writeAllSteps_ )
      simState_->BeginMultiSequenceStep( sequenceStep_, analysis_ );
    
    // In multiharmonic analysis we speak in terms of multiples of base-harmonics.
    // The system matrices of the single harmonics get inserted into the global matrix
    // and then the system is solved.
    // Therefore we don't consider one frequency isolated from the others
    // because they are coupled

    // Perform steps for the solution
    ptPDE_->GetSolveStep()->SetMultHarmonicFreq( harmFreq_ );

    ptPDE_->GetSolveStep()->PreStepHarmonic();
    ptPDE_->GetSolveStep()->SolveStepHarmonic();
    ptPDE_->GetSolveStep()->PostStepHarmonic();

    for(UInt i = 0; i < harmFreq_.GetSize(); ++i  ){
      // current frequency
      Double actFreq = harmFreq_[i];

      // the harmonics number h can be negative but the frequency step can't
      // therefore use the index i
      // i = [  0     1     2  ...   N    N+1   N+2 ...  2N ]
      // h = [ -N   -N+1  -N+2 ...   0     1     2  ...   N ]
      // NOTE: If the new (performance optimized) version is used,
      // only odd harmonics are considered and therefore the above mapping looks like
      // index:     [  0     1     2  ... (N-1)/2   (N+1)/2   (N-1)/2+2   ...  N+1 ]
      // harmonic:  [ -N   -N+2  -N+4 ...    -1        1           3      ...   N ]

      UInt actFreqStep = i;
      analysis_id_.step = actFreqStep;
      analysis_id_.freq = actFreq;

      Integer h = this->HarmonicOfIndex(i);
      // We need to activate the correct harmonic results in CoefFunctionHarmBalance
      mathParser_->SetValue(MathParser::GLOB_HANDLER, "f", actFreq);
      mathParser_->SetValue(MathParser::GLOB_HANDLER, "harmonicHandle", h);

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

      simState_->WriteStep( actFreqStep, actFreq);
    }

    handler_->FinishMultiSequenceStep();
    simState_->FinishMultiSequenceStep( !abortSimulation_ );

    // Perform finalization only if not part of sequence
    if(!isPartOfSequence_)
      handler_->Finalize();
  }
  
}
