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

BaseDriver::BaseDriver( shared_ptr<SimState> simState, Domain* myDom,
                        PtrParamNode paramNode, PtrParamNode infoNode)
{
  actSequenceStep_ = 1;
  domain_ = myDom;
  handler_ = domain_->GetResultHandler();
  simState_ = simState;

  param_ = paramNode;
  info_ = infoNode;
}

BaseDriver::~BaseDriver()
{
  //delete ptdomain_;
}


UInt BaseDriver::GetActSequenceStep() {
  return actSequenceStep_;
}


PtrParamNode BaseDriver::CreateAnalysisId(const std::string& child_name, int child_id, 
                                       const std::string& child_2_name, int child_2_id)
{
  // do we really want to create a new entry? Might blast up the output
  ParamNode::ActionType at = progOpts->DoDetailedInfo() ? ParamNode::APPEND : ParamNode::DEFAULT;
  PtrParamNode child = info_->Get(ParamNode::PN_PROCESS)->Get("step", at);
  child->Get("analysis_id")->SetValue(ConcatAnalysisId(child, child_name, child_id, child_2_name, child_2_id));
  return child;
}

PtrParamNode BaseDriver::CreateAnalysisIdChild(PtrParamNode base, const std::string& child_name, int child_id, 
    const std::string& child_2_name, int child_2_id)
{
  // create a child
  PtrParamNode child = base->Get(child_name);
  std::string val = domain_->GetDriver()->ConcatAnalysisId(base, child_name, child_id, child_2_name, child_2_id);
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
BaseDriver* BaseDriver::CreateInstance(shared_ptr<SimState> state, Domain* myDom,
                                       PtrParamNode paramNode, PtrParamNode infoNode)
{

  BaseDriver *ptdriver = NULL;
  StdVector<std::string>  analysisTypes;

  // Get number of occurences of step-variable
  BasePDE::AnalysisType type;
  UInt numSteps = paramNode->Count( "sequenceStep" );

  // a) if count is bigger than one -> create multiSequence
  if( numSteps == 1 ) {
    std::string analysisString;
    UInt seqStep = 1;
    std::string name = "sequenceStep";
    std::string idx = "index";
    std::string one = "1";

    PtrParamNode seqNode = paramNode->Get("sequenceStep")->Get("analysis");
    PtrParamNode info = infoNode->Get("sequenceStep",ParamNode::APPEND);
    infoNode->Get("index")->SetValue(1);
    analysisString =
        myDom->GetParamRoot()->GetByVal(name, idx, one)->Get("analysis")->GetChild()->GetName();
    type = BasePDE::analysisType.Parse(analysisString);

    // Generate driver
    switch( type ) {
      case BasePDE::STATIC:

        ptdriver = new StaticDriver( seqStep, false, state, myDom, 
                                     seqNode, info );
        break;

      case BasePDE::TRANSIENT:
        ptdriver = new TransientDriver( seqStep, false, state, myDom, 
                                        seqNode, info );
        break;

      case BasePDE::HARMONIC:
        ptdriver = new HarmonicDriver( seqStep, false, state, myDom, 
                                       seqNode, info  );
        break;

      case BasePDE::EIGENFREQUENCY:
        ptdriver = new EigenFrequencyDriver( seqStep, false, state, myDom, 
                                             seqNode, info );
        break;

      default:
        EXCEPTION( "Could not create driver" );
    }

    // b) create multiSequence driver
  } else if( numSteps > 1 ) {
    ptdriver = new MultiSequenceDriver(state, myDom, paramNode, infoNode);
  } else {
    EXCEPTION( "At least one sequenceStep has to be provided" );
  }
  // Pass driver to domain object
  return ptdriver;
}
