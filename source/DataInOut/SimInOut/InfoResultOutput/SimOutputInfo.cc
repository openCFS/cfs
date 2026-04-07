#include "SimOutputInfo.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Utils/tools.hh"
#include "Domain/Domain.hh"

using namespace CoupledField;
using std::string;

SimOutputInfo::SimOutputInfo(PtrParamNode outputNode, PtrParamNode infoNode,
                             bool isRestart ) : 
    SimOutput("", outputNode, infoNode, isRestart )
{
  
  // The restart case is currently not implemented, i.e. resuls from a 
  // partial simulation get overwritten.
  if( isRestart_ ) {
    WARN( "The INFO-Writer is currently not adapted to write restarted "
          "results correctly, thus the results of the previous run get"
          " overwritten." );
  }
  
  formatName_ = "info";
  capabilities_.insert(HISTORY);
  dirName_ = ".";
  
  
  info_root = myInfo_->Get("calculation")->Get(ParamNode::PROCESS);
}

SimOutputInfo::~SimOutputInfo() 
{
}

void SimOutputInfo::Init(Grid * ptGrid, bool printGridOnly) 
{
}



//! Begin multisequence step
void SimOutputInfo::BeginMultiSequenceStep(UInt step, BasePDE::AnalysisType type, UInt numSteps)
{
  actMSStep_ = step;
}

//! Register result (within one multisequence step)
void SimOutputInfo::RegisterResult(shared_ptr<BaseResult> br, UInt saveBegin, UInt saveInc, UInt saveEnd, bool isHistory)
{
   shared_ptr<ResultInfo> ri = br->GetResultInfo();
   PtrParamNode in = info_root->GetByVal("sequence", "step", boost::lexical_cast<string>(actMSStep_));
   // create new
   in = in->Get("result", ParamNode::APPEND);
   br->SetInfoNode(in);
   in->Get("data")->SetValue(ToValidXML(ri->resultName));
   in->Get("location")->SetValue(br->GetEntityList()->GetName());
   string loc_type;
   ResultInfo::Enum2String(ri->definedOn, loc_type);
   in->Get("definedOn")->SetValue(loc_type);
}

void SimOutputInfo::BeginStep( UInt stepNum, Double stepVal) 
{
  actStep_ = stepNum;
  actStepVal_ = stepVal;
}


void SimOutputInfo::AddResult(shared_ptr<BaseResult> br) 
{

  shared_ptr<ResultInfo> ri = br->GetResultInfo();
  // we'll get either a real or complex result set
  Vector<Double>*  d_vec = NULL;
  Vector<Complex>* c_vec = NULL;
  if(br->GetEntryType() ==  BaseMatrix::DOUBLE) {
    d_vec = dynamic_cast<Vector<Double>*>(br->GetSingleVector()); }
  else {
    c_vec = dynamic_cast<Vector<Complex>*>(br->GetSingleVector()); }

  // loop over individual regions
  shared_ptr<EntityList> list = br->GetEntityList();
  EntityIterator it = list->GetIterator();
  // we also need the index for strange code below
  for(it.Begin(); !it.IsEnd(); it++) 
  {
    PtrParamNode in = br->GetInfoNode()->Get("item", ParamNode::APPEND);
    in->Get("step")->SetValue(actStep_);
    in->Get("step_val")->SetValue(actStepVal_);

    // print value(s)
    StdVector<string>& dofs = ri->dofNames; 
    for(unsigned int i = 0; i < dofs.GetSize(); i++) 
    {
      PtrParamNode in_ = in->Get(dofs[i] == "" ? "value" : dofs[i], ParamNode::APPEND); 
      if(br->GetEntryType() ==  BaseMatrix::DOUBLE)
        in_->SetValue((*d_vec)[it.GetPos() * dofs.GetSize() + i]);
      else
      {
        in_->SetValue((*c_vec)[it.GetPos() * dofs.GetSize() + i]);
      }
    }
    if(ri->unit != "") in->Get("unit")->SetValue(ri->unit);

    if(ri->definedOn != ResultInfo::REGION && ri->definedOn != ResultInfo::SURF_REGION) 
      in->Get("id")->SetValue(it.GetIdString()); // element/node/... number
  }
}

