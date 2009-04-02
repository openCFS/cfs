#include "SimOutputInfo.hh"

#include "DataInOut/ParamHandling/InfoNode.hh"

using namespace CoupledField;
using std::string;

SimOutputInfo::SimOutputInfo(ParamNode * outputNode ) : SimOutput("", outputNode )
{
  formatName_ = "info";
  capabilities_.insert(HISTORY);
  dirName_ = ".";
  
  output = outputNode != NULL ?outputNode->Get("output")->AsBool() : true;
  
  info_root = info->Get("calculation")->Get(InfoNode::PROCESS);
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
   InfoNode* in = info_root->Get("sequence", "step", boost::lexical_cast<string>(actMSStep_));
   // create new
   in = in->Get("result", InfoNode::APPEND);
   br->SetInfoNode(in);
   in->Get("data")->SetValue(ri->resultName);
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
  if(!output) return; // xml trigger

  shared_ptr<ResultInfo> ri = br->GetResultInfo();

  // we'll get either a real or complex result set
  Vector<Double>*  d_vec = NULL;
  Vector<Complex>* c_vec = NULL;
  if(br->GetEntryType() ==  BaseMatrix::DOUBLE) {
    d_vec = dynamic_cast<Vector<Double>*>(br->GetSingleVector()); }
  else {
    c_vec = dynamic_cast<Vector<Complex>*>(br->GetSingleVector()); }

  // loop over indiviual regions
  shared_ptr<EntityList> list = br->GetEntityList();
  EntityIterator it = list->GetIterator();
  // we also need the index for strange code below
  for(it.Begin(); !it.IsEnd(); it++) 
  {
    InfoNode* in = br->GetInfoNode()->Get("item", InfoNode::APPEND);

    in->Get("step_nr")->SetValue(actStep_);
    in->Get("step_val")->SetValue(actStepVal_);

    // print value(s)
    StdVector<string>& dofs = ri->dofNames; 
    for(unsigned int i = 0; i < dofs.GetSize(); i++) 
    {
      InfoNode* in_ = in->Get(dofs[i] == "" ? "value" : dofs[i], InfoNode::APPEND); 
      if(br->GetEntryType() ==  BaseMatrix::DOUBLE)
        in_->SetValue((*d_vec)[it.GetPos() * dofs.GetSize() + i]);
      else
        in_->SetValue((*c_vec)[it.GetPos() * dofs.GetSize() + i]);
    }
    if(ri->unit != "") in->Get("unit")->SetValue(ri->unit);

    if(ri->definedOn != ResultInfo::REGION && ri->definedOn != ResultInfo::SURF_REGION) 
      in->Get("id")->SetValue(it.GetIdString()); // element/node/... number
  }
}

