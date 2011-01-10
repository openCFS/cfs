#include "Optimization/Excitation.hh"
#include "Optimization/Design/DesignStructure.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Objective.hh"
#include "Optimization/Condition.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "General/exception.hh"
#include "Domain/domain.hh"

using namespace CoupledField;

Enum<ObjectiveContainer::StoppingRule::Type>  ObjectiveContainer::StoppingRule::type;

Objective::Objective(PtrParamNode pn, PtrParamNode pn_type, unsigned int idx)
 : Function(pn_type)
{
  // multiple excitation is handled in Optimization itself!

  // the current value -> check <Get/Set>Value() when altering the presets!
  this->index_       = idx;

  this->penalty_ = pn_type->Has("penalty") ? pn_type->Get("penalty")->As<Double>() : 1.0;

  get<0>(coord) = -1;
  get<1>(coord) = -1;
  get<2>(coord) = 1.0;
  if(pn_type->Has("coord"))
  {
    if(pn_type->Get("coord")->As<std::string>() == "all" && type_ == HOMOGENIZATION_TENSOR)
      EXCEPTION("homogenization tensor as objective does not support coord='all'!");
    ParseCoord(pn_type, coord);
  }

  this->pn = pn;
}

std::string Objective::GetName() const
{
  if(get<0>(coord) == -1)
    return type.ToString(type_);
  else
    return type.ToString(type_) + "E" + lexical_cast<std::string>(get<0>(coord)) + lexical_cast<std::string>(get<1>(coord));
}


void Objective::ToInfo(PtrParamNode info)
{
  Function::ToInfo(info);
  if(tensor_.GetNumRows() > 1)
    info->Get("tensor")->SetValue(new Matrix<double>(tensor_));

}

ObjectiveContainer::StoppingRule::StoppingRule()
{
  // sync with XSL defaults!
  value = 1e-3;
  queue = 5;
  type_ = DESIGN_CHANGE;
}

void ObjectiveContainer::StoppingRule::Init(PtrParamNode pn)
{
  if(pn == NULL) return;

  type_ = type.Parse(pn->Get("type")->As<std::string>());
  value = pn->Get("value")->As<Double>();
  queue = pn->Get("queue")->As<Integer>();
}

ObjectiveContainer::ObjectiveContainer()
{
  last_design_ = -121354;
}


ObjectiveContainer::~ObjectiveContainer()
{
  for(unsigned int i = 0; i < data.GetSize(); i++) {
    delete data[i];
    data[i] = NULL;
  }
}

void ObjectiveContainer::Read(PtrParamNode obj_node)
{
  this->minimize_ = obj_node->Get("task")->As<std::string>() == "minimize";

  // depending on the costFunction attribute type we read the multiObjective list.
  bool mo = Objective::type.Parse(obj_node->Get("type")->As<std::string>()) == Objective::MULTI_OBJECTIVE;

  // set to default if it is not set
  stop.Init(obj_node->Get("stopping", ParamNode::PASS));

  if(!mo)
  {
    //Objective* tmp = new Objective(obj_node, obj_node, 0);
    data.Resize(1, NULL);
    data[0] = new Objective(obj_node, obj_node, 0);
  }
  else
  {
    if(!obj_node->Has("multiObjective"))
      throw Exception("For costFunction type 'multiObjective' element with 'objective' childs is needed");

    ParamNodeList list = obj_node->Get("multiObjective")->GetList("objective");
    if(list.GetSize() == 0)
      throw Exception("For costFunction type 'multiObjective' element with 'objective' childs is needed");

    data.Resize(list.GetSize(), NULL);

    for(unsigned int i = 0; i < list.GetSize(); i++)
      data[i] = new Objective(obj_node, list[i], i);

    if(Has(Objective::MULTI_OBJECTIVE))
      throw Exception("special objective type 'multiObjective' not allowed in 'multiObjective' list");
  }
}

void ObjectiveContainer::PostProc(DesignSpace* space, DesignStructure* structure, MultipleExcitation* me)
{
  for(unsigned int i = 0; i < data.GetSize(); i++)
  {
    data[i]->PostProc(space, structure);
    data[i]->SetExcitation(me);
  }
}

void ObjectiveContainer::ToInfo(PtrParamNode in)
{
  if(data.GetSize() > 1)
  {
    PtrParamNode m = in->Get("multiObjective");
    for(unsigned int i = 0; i < data.GetSize(); i++)
    {
      PtrParamNode o = m->Get("objective", ParamNode::APPEND);
      Objective* f = data[i];
      f->ToInfo(o);
      o->Get("penalty")->SetValue(f->penalty_); // always for multiobjective
    }
  }
  else
  {
    data[0]->ToInfo(in);
    if(data[0]->GetPenalty() != 1.0) // only when it is set
      in->Get("penalty")->SetValue(data[0]->GetPenalty());
  }

  in->Get("task")->SetValue(minimize_ ? "minimize" : "maximize");
}


bool ObjectiveContainer::Has(Objective::Type type) const
{
  for(unsigned int i = 0, os = data.GetSize(); i < os; i++)
    if(data[i]->GetType() == type) return true;

  return false;
}


Objective* ObjectiveContainer::Get(Objective::Type type, bool throw_exception)
{
  for(unsigned int i = 0, os = data.GetSize(); i < os; i++)
    if(data[i]->GetType() == type) return data[i];

  if(!throw_exception) return NULL;
                  else EXCEPTION("No objective of type " << Objective::type.ToString(type) << " stored");
}


double ObjectiveContainer::GetHistoryValue(bool penalty, int index)
{
  double result = 0.0;

  for(unsigned int i = 0; i < data.GetSize(); i++)
  {
    double val = index == -1 ? data[i]->history_.Last() : data[i]->history_[index];
    result += (penalty ? data[i]->penalty_ : 1.0) * val;
  }

  return result;
}

unsigned int ObjectiveContainer::GetHistorySize()
{
  return data[0]->history_.GetSize();
}

void ObjectiveContainer::PushBackHistory()
{
  for(unsigned int i = 0; i < data.GetSize(); i++)
    data[i]->history_.Push_back(data[i]->value_);
}

void ObjectiveContainer::PushBackDesign(const DesignSpace* space)
{
  // don't push back if the design is the same -> e.g. last commit after convergence
  if(space->GetCurrentDesignId() == last_design_) return;

  // save this iteration - we need a temporary copy to calculate the distance
  Vector<double> curr_design(space->GetNumberOfVariables());
  space->WriteDesignToExtern(curr_design.GetPointer());

  // first iteration?
  double change = last_iteration_.GetSize() == 0 ? curr_design.NormMax() : curr_design.NormMax(last_iteration_);
  last_iteration_ = curr_design;

  design_change.Push_back(change);

  last_design_ = space->GetCurrentDesignId();
}
