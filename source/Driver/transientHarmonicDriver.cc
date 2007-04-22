// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "transientHarmonicDriver.hh"
#include "transientdriver.hh"
#include "stdSolveStep.hh"
#include "harmonicDriver.hh"

#include "DataInOut/CommandLine/BaseCommandLineHandler.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "PDE/StdPDE.hh"
#include "PDE/pdememento.hh"
#include "Domain/domain.hh"
#include "DataInOut/resultHandler.hh"

namespace CoupledField {


  // ===============
  //   Constructor
  // ===============
  TransientHarmonicDriver::TransientHarmonicDriver( UInt stepOffset,
                                                    Double timeOffset, 
                                                    UInt sequenceStep,
                                                    bool isPartOfSequence)
    : SingleDriver( sequenceStep, isPartOfSequence ),
      TransientDriver( stepOffset, timeOffset, sequenceStep, isPartOfSequence),
      HarmonicDriver( stepOffset, 0.0, sequenceStep, isPartOfSequence )
  {
    ENTER_FCN( "TransientHarmonicDriver::TransientHarmonicDriver" );
    

  }

  void TransientHarmonicDriver::Init() {
    ENTER_FCN( "TransientHarmonicDriver::Init" );

    // get parameter node
    ParamNode * myNode = 
      param->Get("sequenceStep", "index", GenStr( sequenceStep_ ) )
      ->Get("analysis")->Get("transientHarmonic");

    // Be verbose!
    Info->PrintF( "", "\n TransientHarmonicDriver:\n" ); 
    
    StdVector<ParamNode*> pdeNodes;
    pdeNodes =  param->Get("sequenceStep", "index", GenStr(sequenceStep_) )
      ->Get("pdeList")->GetChildren();    
    
    std::string pdeName, analysis;
    for ( UInt k=0; k < pdeNodes.GetSize(); k++) {

      // get name and analysistype of current pde
      pdeName = pdeNodes[k]->GetName();
      analysis = myNode->Get("analysisPart", "pdeName", pdeName)
        ->Get("type")->AsString();

      if ( analysis == "transient" ) {
        Info->PrintF( "", "\nTransient simulation for '%s'", pdeName.c_str() );
      }
      else if ( analysis == "harmonic" ) {
        Info->PrintF( "", "\nHarmonic simulation for '%s'", pdeName.c_str() );
        Info->PrintF( "", "\n -> frequency:  %f Hz\n", startFreq_ );
      }

    } 
    Info->PrintF( "", "\n" );

    Info->StartProgress ("Starting to solve problem", false);
  }

  // ==============
  //   Destructor
  // ==============
  TransientHarmonicDriver::~TransientHarmonicDriver()
  {
    ENTER_FCN( "TransientHarmonicDriver::~TransientHarmonicDriver" );
  }

  
  AnalysisType TransientHarmonicDriver::
  GetAnalysisType( const std::string& pdeName ) {
    ENTER_FCN( "TransientHarmonicDriver::GetAnalysisType" );

    // get parameter node
    ParamNode * myNode = 
      param->Get("sequenceStep", "index", GenStr( sequenceStep_ ) )
      ->Get("analysis")->Get("transientHarmonic");

    std::string actTypeString;
    actTypeString = myNode->Get("analysisPart", "pdeName", pdeName)->AsString();
    AnalysisType actType;
    String2Enum(actTypeString, actType);

    return actType;
  }

  UInt TransientHarmonicDriver::GetActStep( const std::string& pdename ) {
    ENTER_FCN( "TransientHarmonicDriver::GetActStep" );

    if( GetAnalysisType( pdename) == HARMONIC ) {
      return 1;
    } else {
      return actTimeStep_;
    }
  }

  // =================
  //   Solve Problem
  // =================
  void TransientHarmonicDriver::SolveProblem()
  {
    ENTER_FCN( "TransientHarmonicDriver::SolveProblem" );
  
    // notify resultHandler about beginning of new sequence step 
    ResultHandler * resHandler = domain->GetResultHandler();
    if( !isPartOfSequence_ )
      resHandler->BeginMultiSequenceStep( sequenceStep_,
                                          analysis_,
                                          numstep_);

    // transient settings
    Double  steptime  = firstdt_;
    Double  dt = firstdt_;
    UInt startStep=1;

    // harmonic settings
    actFreq_ = startFreq_; // Works for only one frequency!!

    // Output time steps in steps of ten percent 
    Double timeStepPercent = (double)numstep_/10;
    Double percentCounter = timeStepPercent;

    ptPDE_->WriteGeneralPDEdefines();

    // Solve problem
    ptPDE_->GetSolveStep()->SetNumTimeSteps(numstep_);
    ptPDE_->GetSolveStep()->SetStartStep(startStep);



    for (UInt actTimeStep_ = 1; actTimeStep_ <= numstep_; actTimeStep_++) {

      if ((double)actTimeStep_ >= percentCounter  ){           
        Info->WriteTimeStep(ptPDE_->GetName(), 
                            actTimeStep_+TransientDriver::stepOffset_, 
                            steptime+timeOffset_);
        percentCounter += timeStepPercent;
      }

      // Set curent frequency and time value in the mathParser
      domain->GetMathParser()->SetValue( MathParser::GLOB_HANDLER,
                                         "f", actFreq_ );
      domain->GetMathParser()->SetValue( MathParser::GLOB_HANDLER,
                                         "t", steptime );

      // transient settings
      ptPDE_->GetSolveStep()->SetActTime(steptime);

      // harmonic settings
      ptPDE_->GetSolveStep()->SetActFreq( actFreq_ );

      ptPDE_->GetSolveStep()->SetActStep( actTimeStep_ );

      ptPDE_->GetSolveStep()->PreStepTransHarmonic();
      ptPDE_->GetSolveStep()->SolveStepTransHarmonic();
      ptPDE_->GetSolveStep()->PostStepTransHarmonic();

      // writing results in output-file
      resHandler->BeginStep( TransientDriver::stepOffset_ + actTimeStep_,
                             timeOffset_ + steptime );
      ptPDE_->WriteResultsInFile(actTimeStep_, steptime, 
                                 TransientDriver::stepOffset_, timeOffset_);
      resHandler->FinishStep();

      steptime+=dt;        

      SETPROFILE("After Transient Step");
    }

    ptPDE_->Finalize();
  
    // notify resultHandler about finishing of current sequence step
    if( !isPartOfSequence_ )
      resHandler->FinishMultiSequenceStep();
    

  }
} // end of namespace
