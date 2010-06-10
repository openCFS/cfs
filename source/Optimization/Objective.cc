#include "Optimization/Objective.hh"
#include "Optimization/Condition.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "General/exception.hh"
#include "Domain/domain.hh"

using namespace CoupledField;

Objective::Objective(PtrParamNode pn, PtrParamNode pn_type, unsigned int idx)
 : Function(pn_type)
{
  // multiple excitation is handled in Optimization itself!

  // the current value -> check <Get/Set>Value() when altering the presets!
  this->index_       = idx;

  this->penalty_ = pn_type->Has("penalty") ? pn_type->Get("penalty")->As<Double>() : 1.0;

  // currently we have only relative implemented up to now!
  stop.value = pn->Has("stopping") ? pn->Get("stopping/value")->As<Double>() : 1e-5;
  stop.queue = pn->Has("stopping") ? pn->Get("stopping/queue")->As<Integer>() : 5; // is >= 1!

  coord.first = -1;
  coord.second = -1;
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
  if(coord.first == -1)
    return type.ToString(type_);
  else
    return type.ToString(type_) + "E" + lexical_cast<std::string>(coord.first) + lexical_cast<std::string>(coord.first);
}


ObjectiveContainer::ObjectiveContainer()
{
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


void ObjectiveContainer::ToInfo(PtrParamNode in) const
{
  bool mo = data.GetSize() > 1;
  if(mo)
  {
    PtrParamNode m = in->Get("multiObjective");
    for(unsigned int i = 0; i < data.GetSize(); i++)
    {
      PtrParamNode o = m->Get("objective", ParamNode::APPEND);
      Objective* f = data[i];
      o->Get("type")->SetValue(f->GetName());
      o->Get("penalty")->SetValue(f->penalty_);
      o->Get("evaluate")->SetValue(f->evaluateOnce_ ? "once" : "always");
    }
  }
  else
  {
    in->Get("type")->SetValue(data[0]->GetName());
    in->Get("evaluate")->SetValue(data[0]->evaluateOnce_ ? "once" : "always");
  }

  in->Get("task")->SetValue(minimize_ ? "minimize" : "maximize");
  if(data[0]->harmonic_)
    in->Get("factor")->Get("omega_omega")->SetValue(data[0]->omega_omega_);
  if(data[0]->tensor_.GetNumRows() > 1)
    in->Get("tensor")->SetValue(new Matrix<double>(data[0]->tensor_));
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

