#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "transientHarmonicDriver.hh"
#include "transientdriver.hh"
#include "stdSolveStep.hh"
#include "harmonicDriver.hh"

#include "DataInOut/CommandLine/BaseCommandLineHandler.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "PDE/StdPDE.hh"
#include "PDE/pdememento.hh"
#include "Domain/domain.hh"

namespace CoupledField {


  // ===============
  //   Constructor
  // ===============
  TransientHarmonicDriver::TransientHarmonicDriver( UInt stepOffset,
                                                    Double timeOffset, 
                                                    std::string driverTag,
                                                    bool isPartOfSequence)
    : TransientDriver( stepOffset, timeOffset, driverTag, isPartOfSequence),
      HarmonicDriver( stepOffset, 0, driverTag, isPartOfSequence )
      
  {
    ENTER_FCN( "TransientHarmonicDriver::TransientHarmonicDriver" );
    

  }

  void TransientHarmonicDriver::Init() {
    ENTER_FCN( "TransientHarmonicDriver::Init" );

    // vectors for accessing parameters
    StdVector<std::string> keyVec, attrVec, valVec;
    attrVec = "tag";
    valVec  = driverTag_;

    // Be verbose!
    Info->PrintF( "", "\n TransientHarmonicDriver:\n" ); 

    StdVector<std::string> pdeNames, analysisType;
    params->GetPDEList( pdeNames );
    for ( UInt k=0; k < pdeNames.GetSize(); k++) {

      keyVec = "transientHarmonic", "analysisPart", "type";
      attrVec = "", "pdeName";
      valVec = "", pdeNames[k];
      params->GetList(keyVec, attrVec, valVec, analysisType);

      if ( (analysisType.GetSize() > 1) || analysisType.IsEmpty() ) {
        (*error) << "Please specify for PDE " << pdeNames[k]
                 << " either transient or harmonic calculation";
        Error( __FILE__, __LINE__ );
      }
      else if ( analysisType[0] == "transient" ) {
        Info->PrintF( "", "\nTransient simulation for '%s'", pdeNames[k].c_str() );
      }
      else if ( analysisType[0] == "harmonic" ) {
        Info->PrintF( "", "\nHarmonic simulation for '%s'", pdeNames[k].c_str() );
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

    StdVector<std::string> keyVec, attrVec, valVec;
    keyVec = "transientHarmonic", "analysisPart", "type";
    attrVec = "", "pdeName";
    valVec = "", pdeName;
    std::string actTypeString;
    AnalysisType actType;
    
    params->Get(keyVec, attrVec, valVec, actTypeString );
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
  
    // transient settings
    Double  steptime  = firstdt_;
    UInt stepsave     = isavebegin_;
    Double  dt = firstdt_;
    UInt startStep=1;

    // harmonic settings
    actFreq_ = startFreq_; // Works for only one frequency!!

    // Output time steps in steps of ten percent 
    Double timeStepPercent = (double)numstep_/10;
    Double percentCounter = timeStepPercent;

  
    // if driver is not part of multiSequence Driver, get list
    // of pdes which have to be solved and intialize them
    if (isPartOfSequence_ == false) {

      // writes out the grid one time
      domain->PrintGrid();
    }

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

      ptPDE_->PostProcess();   
      // write history data
      ptPDE_->WriteHistoryInFile(actTimeStep_, steptime, 
                                 TransientDriver::stepOffset_, timeOffset_);

    
      // writing results in output-file
      if (actTimeStep_ == stepsave && (actTimeStep_ <= isaveend_)) { 
        ptPDE_->WriteResultsInFile(actTimeStep_, steptime, 
                                   TransientDriver::stepOffset_, timeOffset_);

        stepsave+=isaveincr_;
      }

      steptime+=dt;        

      SETPROFILE("After Transient Step");
    }
  
  }
} // end of namespace
