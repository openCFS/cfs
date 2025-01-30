#include "BaseDriver.hh"
#include "driver_header.hh"

#include <fstream>
#include <iostream>
#include <string>

#include "Utils/ToolsFull.hh"
#include "Utils/preciceAdapter/IPreciceAdapter.hh"
#include "Utils/preciceAdapter/TransientDriverPrecice.hh"
#include "Utils/StdVector.hh"
#include "General/Exception.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ProgramOptions.hh"
#include "DataInOut/ResultHandler.hh"
#include "Domain/Domain.hh"

using namespace CoupledField;

// initialize static HALT flag
bool BaseDriver::abortSimulation_ = false;

BaseDriver::BaseDriver( shared_ptr<SimState> simState, Domain* myDom,
                        PtrParamNode paramNode, PtrParamNode infoNode)
{
  sequenceStep_ = 0;
  domain_ = myDom;
  handler_ = domain_->GetResultHandler();
  simState_ = simState;

  param_ = paramNode;
  info_ = infoNode;
  analysis_ = BasePDE::NO_ANALYSIS;
}

BaseDriver::~BaseDriver()
{
}


UInt BaseDriver::GetActSequenceStep() {
  return sequenceStep_;
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
    PtrParamNode info = infoNode->GetByVal("sequenceStep","sequence", 1);

    analysisString = myDom->GetParamRoot()->GetByVal(name, idx, one)->Get("analysis")->GetChild()->GetName();
    type = BasePDE::analysisType.Parse(analysisString);

    // Generate driver
    switch( type ) {
      case BasePDE::STATIC:
        ptdriver = new StaticDriver( seqStep, false, state, myDom, seqNode, info );
        break;

      case BasePDE::TRANSIENT:
        // A registered, non-dummy preCICE adapter -> use the coupling driver. Secondary/source
        // domains (e.g. a previous sequence step loaded via SimState::GetDomain during a field
        // read) never get an adapter registered, so preciceAdapter_ is null there: guard it and
        // fall back to the plain transient driver.
        if(myDom->GetPreciceAdapter() && !myDom->GetPreciceAdapter()->IsPreciceDummy())
          ptdriver = new TransientDriverPrecice( seqStep, false, state, myDom, seqNode, info );
        else
          ptdriver = new TransientDriver( seqStep, false, state, myDom, seqNode, info );
        break;

      case BasePDE::HARMONIC:
        ptdriver = new HarmonicDriver( seqStep, false, state, myDom, seqNode, info  );
        break;

      case BasePDE::MULTIHARMONIC:
        ptdriver = new MultiHarmonicDriver( seqStep, false, state, myDom, seqNode, info  );
        break;

      case BasePDE::INVERSESOURCE:
        ptdriver = new InverseSourceDriver( seqStep, false, state, myDom, seqNode, info  );
        break;

      case BasePDE::EIGENFREQUENCY:
        ptdriver = new EigenFrequencyDriver( seqStep, false, state, myDom, seqNode, info );
        break;

      case BasePDE::BUCKLING:
        ptdriver = new BucklingDriver( seqStep, false, state, myDom, seqNode, info );
        break;

      case BasePDE::EIGENVALUE:
        ptdriver = new EigenValueDriver( seqStep, false, state, myDom, seqNode, info);
        break;

      case BasePDE::HARMONIC25D:
        ptdriver = new Harmonic25DDriver( seqStep, false, state, myDom, seqNode, info  );
        break;

      default:
        EXCEPTION( "Could not create driver" );
    }

  } else if( numSteps > 1 ) {
    bool keep = domain->GetParamRoot()->Has("optimization");
    // serve the entries for the info.xml to have optimization las
    for(unsigned int i = 0; i < numSteps; i++)
      infoNode->GetByVal("sequenceStep","sequence", i+1); // 1-based

    ptdriver = new MultiSequenceDriver(state, myDom, paramNode, infoNode, keep);
  } else {
    EXCEPTION( "At least one sequenceStep has to be provided" );
  }
  
  // Pass driver to domain object
  return ptdriver;
}
