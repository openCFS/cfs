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
  TransientHarmonicDriver::TransientHarmonicDriver(Domain * adomain, 
                                                   UInt stepOffset,
                                                   Double timeOffset, 
                                                   std::string driverTag,
                                                   bool isPartOfSequence)
    : SingleDriver(adomain, stepOffset, timeOffset, 
                   driverTag, isPartOfSequence)
  {
    ENTER_FCN( "TransientHarmonicDriver::TransientHarmonicDriver" );
    
    // set correct analysistype
    analysis_ = TRANSIENTHARMONIC;

    // vectors for accessing parameters
    StdVector<std::string> keyVec, attrVec, valVec;
    attrVec = "tag";
    valVec  = driverTag_;



    // Get harmonic information from parameter object
    // ==============================================
    StdVector<Double> auxVec;
    keyVec = "harmonic", "startFreq";
    params->GetList( keyVec, attrVec, valVec, auxVec );
    if ( auxVec.GetSize() == 1 ) {
      startFreq_ = auxVec[0];
    }
    else {
      (*error) << "Could not find a section <harmonic> with tag = '"
               << driverTag_ << "' in " << commandLine->GetParamFile();
      Error( __FILE__, __LINE__ );
    }

    keyVec = "harmonic", "stopFreq";
    params->Get( keyVec, attrVec, valVec, stopFreq_ );

    keyVec = "harmonic", "numFreq";
    params->Get( keyVec, attrVec, valVec, numFreq_ );

    // Dirty Hack!!!
    // We allow only one frequency, so check that start and stop are equal
    if ( ( numFreq_ != 1 ) && ( startFreq_ != stopFreq_ ) ) {
      (*error) << "Found numFreq = " << numFreq_ << " in xml-File, "
               << "and startFreq_ = " << startFreq_ << " != "
               << "stopFreq_ = " << stopFreq_;
      Error( __FILE__, __LINE__ );

    }


    // Get time stepping information from parameter object
    // ===================================================
    keyVec = "transient", "firstDt";
    params->GetList( keyVec, attrVec, valVec, auxVec );
    if ( auxVec.GetSize() == 1 ) {
      firstdt_ = auxVec[0];
    }
    else {
      (*error) << "Could not find a section <transient> with tag = '"
               << driverTag_ << "' in " << commandLine->GetParamFile();
      Error( __FILE__, __LINE__ );
    }

    keyVec  = "transient", "numSteps";
    params->Get(keyVec, attrVec, valVec, numstep_);

    keyVec = "transient", "stepSaveBeg";
    params->Get(keyVec, attrVec, valVec, isavebegin_);

    keyVec = "transient", "stepSaveEnd";
    params->Get(keyVec, attrVec, valVec, isaveend_);

    keyVec = "transient", "stepSaveInc";
    params->Get(keyVec, attrVec, valVec, isaveincr_);

    keyVec = "transient", "writeRestartInc";
    params->Get(keyVec, attrVec, valVec, restartIncr_);

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
      ptdomain_->PrintGrid();

      GetMyPDEs();
      Info->StartProgress ("Starting to solve problem", false);
    }

    ptPDE_->WriteGeneralPDEdefines();

    // Solve problem
    ptPDE_->GetSolveStep()->SetTimeStep(dt);
    ptPDE_->GetSolveStep()->SetNumTimeSteps(numstep_);
    ptPDE_->GetSolveStep()->SetStartStep(startStep);



    for (UInt nstep = 1; nstep <= numstep_; nstep++) {

      if ((double)nstep >= percentCounter  ){           
        Info->WriteTimeStep(ptPDE_->GetName(), nstep+stepOffset_, 
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

      ptPDE_->GetSolveStep()->SetActStep( nstep );

      ptPDE_->GetSolveStep()->PreStepTransHarmonic();
      ptPDE_->GetSolveStep()->SolveStepTransHarmonic();
      ptPDE_->GetSolveStep()->PostStepTransHarmonic();

      ptPDE_->PostProcess();   
      // write history data
      ptPDE_->WriteHistoryInFile(nstep, steptime, stepOffset_, timeOffset_);

    
      // writing results in output-file
      if (nstep == stepsave && (nstep <= isaveend_)) { 
        ptPDE_->WriteResultsInFile(nstep, steptime, stepOffset_, timeOffset_);

        stepsave+=isaveincr_;
      }

      steptime+=dt;        

      SETPROFILE("After Transient Step");
    }
  
  }
} // end of namespace
