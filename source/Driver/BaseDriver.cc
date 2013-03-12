#include "BaseDriver.hh"
#include "driver_header.hh"

#include <fstream>
#include <iostream>
#include <string>

#include "Utils/tools.hh"
#include "Utils/StdVector.hh"
#include "General/Exception.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ProgramOptions.hh"
#include "DataInOut/ResultHandler.hh"

#include "Domain/Domain.hh"

using namespace CoupledField;

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
  PtrParamNode child = driverNode->Get(ParamNode::PN_PROCESS)->Get("step", at);
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
  
  std::stringstream ss;
  ss << (analysis_id->Has("analysis_id") ? analysis_id->Get("analysis_id")->As<std::string>() : ""); 

  if(child_name != "")
  ss << ":" << child_name;

  if(child_id != -1) 
    ss << ":" << child_id;
  
  if(child_2_name != "")
    ss << ":" << child_2_name; 
  
  if(child_2_id != -1) 
      ss << ":" << child_2_id;

  return ss.str();
}

// static stuff
BaseDriver* BaseDriver::CreateInstance()
{

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
        REFACTOR;
//        ptdriver = new piezoParamIdent(0, 0.0 );
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
