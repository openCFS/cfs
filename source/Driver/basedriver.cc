// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <string>

#include "basedriver.hh"
#include "Utils/tools.hh"
#include "Utils/StdVector.hh"
#include "General/exception.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/resultHandler.hh"
#include "Driver/driver_header.hh"
#include "Domain/domain.hh"

using namespace CoupledField;

BaseDriver::BaseDriver( )
{
  actSequenceStep_ = 1;
  nummeshes_=0;
  handler_ = domain->GetResultHandler();
}

BaseDriver::~BaseDriver()
{
  //delete ptdomain_;
}


UInt BaseDriver::GetActSequenceStep() {
  return actSequenceStep_;
}

// for computation with adaptivity
bool BaseDriver::printMeshesOrNot() {


  EXCEPTION( "Currently not working, need change to XML-Standard" );
  bool meshesInfo=false;

  return meshesInfo;
}

void BaseDriver::PrintSeqMeshes()
{
  Warning( "Not implemented anymore", __FILE__, __LINE__ );
}
// static stuff
BaseDriver* BaseDriver::CreateInstance()
{
#ifdef ADAPTGRID
  flags = new Flags();
#endif

  // Generate log message
  Info->StartProgress( "Generating driver" );

  BaseDriver   *ptdriver = NULL;
  StdVector<std::string>  analysisTypes;

  // Get number of occurences of step-variable
  BasePDE::AnalysisType type;
  UInt numSteps = param->Count( "sequenceStep" );

  // a) if count is bigger than one -> create multiSequence
  if( numSteps == 1 ) {
    std::string analysisString;
    UInt seqStep = 1;
    std::string name = "sequenceStep";
    std::string idx = "index";
    std::string one = "1";

    analysisString = 
      param->Get(name, idx, one)->Get("analysis")->GetChild()->GetName();
    type = BasePDE::analysisType.Parse(analysisString);

    // Generate driver
    switch( type ) {
    case BasePDE::STATIC:

      ptdriver = new StaticDriver( seqStep );
      break;

    case BasePDE::TRANSIENT:
      ptdriver = new TransientDriver( seqStep );
      break;

    case BasePDE::HARMONIC:
      // calls Driver for parameter identification, using harmonic analysis
      if ( analysisString == "paramIdent" ) {
        ptdriver = new piezoParamIdent(0, 0.0 );
      }
      else
        ptdriver = new HarmonicDriver( seqStep );
      break;

    case BasePDE::EIGENFREQUENCY:
      ptdriver = new EigenFrequencyDriver( seqStep);
      break;

    default:
      EXCEPTION( "Could not create driver" );
    }

    // b) create multiSequence driver    
  } else if( numSteps > 1 ) {
    ptdriver = new MultiSequenceDriver( );
  } else {
    EXCEPTION( "At least one sequenceStep has to be provided" );
  }
  // Pass driver to domain object
  return ptdriver;
}
