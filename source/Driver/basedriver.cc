// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <sstream>
#include <string>

#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/programOptions.hh"
#include "Domain/domain.hh"
#include "Driver/basedriver.hh"
#include "Driver/eigenFrequencyDriver.hh"
#include "Driver/harmonicDriver.hh"
#include "Driver/multiSequenceDriver.hh"
#include "Driver/staticdriver.hh"
#include "Driver/transientdriver.hh"
#include "General/Enum.hh"
#include "General/environment.hh"
#include "General/exception.hh"
#include "ParamIdent/piezoParamIdent.hh"
#include "Utils/StdVector.hh"
#include "basedriver.hh"

using namespace CoupledField;

DECLARE_LOG(driver)
DEFINE_LOG(driver, "driver")

BaseDriver::BaseDriver( )
{
  actSequenceStep_ = 1;
  nummeshes_=0;
  handler_ = domain->GetResultHandler();
  driverNode = info->Get("analysis"); // analysis step set in singleDriver
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
  WARN( "Not implemented anymore" );
}

PtrParamNode BaseDriver::CreateAnalysisId(const std::string& child_name, int child_id, 
                                       const std::string& child_2_name, int child_2_id)
{
  // do we really want to create a new entry? Might blast up the output
  ParamNode::ActionType at = progOpts->DoDetailedInfo() ? ParamNode::APPEND : ParamNode::DEFAULT;
  PtrParamNode child = driverNode->Get(ParamNode::PROCESS)->Get("step", at);
  child->Get("analysis_id")->SetValue(ConcatAnalysisId(child, child_name, child_id, child_2_name, child_2_id));
  return child;
}

PtrParamNode BaseDriver::CreateAnalysisIdChild(PtrParamNode base, const std::string& child_name, int child_id, 
    const std::string& child_2_name, int child_2_id)
{
  // create a child
  PtrParamNode child = base->Get(child_name);
  std::string val = domain->GetDriver()->ConcatAnalysisId(base, child_name, child_id, child_2_name, child_2_id);
  child->Get("analysis_id")->SetValue(val);
  return child;
}


std::string BaseDriver::ConcatAnalysisId(PtrParamNode analysis_id, const std::string& child_name, int child_id, 
                                    const std::string& child_2_name, int child_2_id)
{
  assert(!(child_name == "" && child_id != -1));
  assert(!(child_name == "" && child_2_name != ""));
  assert(!(child_2_name != "" && child_id == -1));
  
  std::string old = (analysis_id->Has("analysis_id") ? analysis_id->Get("analysis_id")->As<std::string>() : "");

  std::stringstream ss;

  // if we simply concatenate, the string becomes very long for optimization and other iterative stuff.
  // therefore create a new string if child_name already is in the string.
  if(old.find(child_name) == std::string::npos) // when child_name is not given we use the old string
    ss << old;

  if(child_name != "")
    ss << (ss.str().length() > 0 ? ":" : "") << child_name;

  if(child_id != -1) 
    ss << ":" << child_id;
  
  if(child_2_name != "")
    ss << ":" << child_2_name; 
  
  if(child_2_id != -1) 
      ss << ":" << child_2_id;

  LOG_DBG3(driver) << "BD:CAI -> " << ss.str();

  return ss.str();
}

// static stuff
BaseDriver* BaseDriver::CreateInstance()
{
#ifdef ADAPTGRID
  flags = new Flags();
#endif

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
      param->GetByVal(name, idx, one)->Get("analysis")->GetChild()->GetName();
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
