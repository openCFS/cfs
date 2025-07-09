#include <stddef.h>
#include <ostream>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "General/defs.hh"
#include "General/Environment.hh"
#include "General/Exception.hh"
#include "MatVec/Matrix.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Objective.hh"
#include "Optimization/Condition.hh"

namespace CoupledField {
class DesignStructure;
}  // namespace CoupledField

using namespace CoupledField;

DEFINE_LOG(obj, "objective")

Enum<StoppingRule::Type>       StoppingRule::type;
Enum<StoppingRule::Condition>  StoppingRule::condition;
Enum<Objective::Term>          Objective::term;


Objective::Objective(PtrParamNode pn, PtrParamNode pn_type, unsigned int idx)
 : Function(pn_type)
{
  // multiple excitation is handled in Optimization itself!

  // the current value -> check <Get/Set>Value() when altering the presets!
  index_       = idx;

  // a tune scale increases the scale parameter automatically, depending on the tune configuration
  if(pn_type->Has("scale"))
  {
    if(pn_type->Get("scale")->As<string>() == "tune")
    {
      tune.Init(pn_type->Get("tune", ParamNode::PASS), Tune::FUNC_SCALE);
      scale = tune.start;
      assert(tune.IsSet());
    }
    else
      scale = pn_type->Get("scale")->As<double>();
  }
  else
    scale = 1.0;

  term_ = pn_type->Has("term") ? term.Parse(pn_type->Get("term")->As<string>()) : LINEAR;

  if(term_ == PENALTY || term_ == POWER)
  {
    string msg = term_ == PENALTY ? "penalty (scale * max(0, func-parameter)^2)" : "power (scale * func^parameter)";
    if(!pn_type->Has("parameter")) // read by Function.cc
      throw Exception("objective attribute 'term' " + msg + " requires 'parameter'");
  }

  get<0>(coord) = -1;
  get<1>(coord) = -1;
  get<2>(coord) = 1.0;
  if(pn_type->Has("coord"))
  {
    if(pn_type->Get("coord")->As<std::string>() == "all" && type_ == HOM_TENSOR)
      EXCEPTION("homogenization tensor as objective does not support coord='all'!");
    ParseCoord(pn_type, coord);
  }

  this->pn = pn;
}

Objective::Objective(Type type, double parameter, Access acc)
{
  Init();
  this->type_ = type;
  this->parameter_ = parameter;
  this->access_ = acc;
  this->scale = 1.0;
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
  info->Get("term")->SetValue(term.ToString(term_));

  if(tune.IsSet())
    tune.ToInfo(info->Get("scale"));
  else
    info->Get("scale")->SetValue(scale);

  if(tensor_.GetNumRows() > 1)
    info->Get("tensor")->SetValue(tensor_);
}

void StoppingRule::Init(PtrParamNode pn)
{
  if(pn == NULL) return;

  type_ = type.Parse(pn->Get("type")->As<std::string>());
  condition_ = condition.Parse(pn->Get("condition")->As<std::string>());
  value = pn->Get("value")->As<double>();
  queue = pn->Get("queue")->As<int>();
  if(pn->Has("function"))
  {
    if(type_ != ABOVE_FUNCTION && type_ != BELOW_FUNCTION)
      throw Exception("'stopping' attribute 'function' only for 'above/belowFunction'");
    function = pn->Get("function")->As<string>();

  }
  if((type_ == ABOVE_FUNCTION || type_ == BELOW_FUNCTION) && function == "")
   throw Exception("attribute 'function' mandatory for 'stopping' type '" + type.ToString(type_) + "'");
}

string StoppingRule::DoStop(ObjectiveContainer* oc, ConditionContainer* cc, StdVector<double>& time)
{
  int hs = oc->GetHistorySize();
  std::stringstream ss;
  switch(type_)
  {
  case MAX_HOURS:
  {
    if(time.GetSize() < 3)
      return "";

    auto cfs_timer = domain->GetInfoRoot()->Get(ParamNode::SUMMARY)->Get("timer")->AsTimer();

    // avg. of the last three in seconds
    double avg = (time.Last() - time[time.GetSize() - 3])/3.0;
    // remaining time to max_hours in seconds
    double remaining = value * 3600 - cfs_timer->GetWallTime();

    //LOG_DBG(opt) << "DSO: wt=" << cfs_timer->GetWallTime() << " cmp=" << time[time.GetSize() - 3] << " avg=" << avg << " mh=" << value << " re=" << remaining << " ti=" << time.ToString();

    if(remaining < 1.5 * avg)
    {
      ss << "Not enough time left to finish within " << value << "h with avg iteration duration " << avg << "s and " << remaining << "s left.";
      return ss.str();
    }
    else
      return "";
  }

  case OSCILLATIONS:
  {
    int cnt = 0;
    for(unsigned int i = 3; i < oc->GetHistorySize(); i++)
    {
      // is different from Function::CalcOscillations() as this combines multi-objectice
      double pp = oc->GetHistoryValue(true, i-2);
      double p = oc->GetHistoryValue(true, i-1);
      double c = oc->GetHistoryValue(true, i);
      if((pp-p)*(p-c) < 0)
      {
        cnt++;
        if(cnt >= value)
        {
          ss << "detected " << cnt << " oscillations in objective function";
          return ss.str();
        }
      }
    }
    return "";
  }

  case REL_COST_CHANGE:
    if(hs <= queue)
      return "";

    for(int i = hs-1; i >= (hs - queue); i--)
    {
      double delta = oc->GetHistoryValue(true, i) - oc->GetHistoryValue(true, i-1);
      double rel = abs(delta / oc->GetHistoryValue(true, i));
      if(rel > value)
        return "";
    }
    ss << "relative change in objective function below " << value << " for " << queue << " iterations";
    return ss.str();

  case DESIGN_CHANGE:
    if(hs <= queue)
      return "";

    for(int n = (int) oc->design_change.GetSize() - 1, i = n; i >= std::max(0, n - queue); i--)
      if(oc->design_change[i] > value)
        return "";

    ss << "design change below " << value << " for " << queue << " iterations";
    return ss.str();

  case ABOVE_FUNCTION:
  case BELOW_FUNCTION:
    {
      Function* f = oc->Get(function, false); // no exception
      if(f == nullptr)
        f = (Function*) cc->Get(function, true); // now exception if we don't have it
      double val = f->GetValue();
      if(type_ == ABOVE_FUNCTION && val <= value)
        return "";
      if(type_ == BELOW_FUNCTION && val > value)
        return "";
      ss << "function " << f->ToString() << " value " << val;
      ss << (type_ == ABOVE_FUNCTION ? " above " : " below ");
      ss << "stopping criteria " << value;
      return ss.str();
    }
  default:
    throw Exception("unhandled stopping criteria");
  }
}

ObjectiveContainer::ObjectiveContainer()
{
  last_design_ = -121354;
  minimize_ = false; // set later
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
  auto sl = obj_node->GetList("stopping");
  stop.Resize(sl.GetSize());
  for(unsigned int i = 0; i < sl.GetSize(); i++)
    stop[i].Init(sl[i]);

  if(!mo)
  {
    //Objective* tmp = new Objective(obj_node, obj_node, 0);
    data.Resize(1, NULL);
    data[0] = new Objective(obj_node, obj_node, 0);
  }
  else
  {
    if(!obj_node->Has("multiObjective"))
      throw Exception("For costFunction type 'multiObjective' element with 'objective' child is needed");

    ParamNodeList list = obj_node->Get("multiObjective")->GetList("objective");
    if(list.GetSize() == 0)
      throw Exception("For costFunction type 'multiObjective' element with 'objective' child is needed");

    data.Resize(list.GetSize(), NULL);

    for(unsigned int i = 0; i < list.GetSize(); i++)
      data[i] = new Objective(obj_node, list[i], i);

    if(Has(Objective::MULTI_OBJECTIVE))
      throw Exception("special objective type 'multiObjective' not allowed in 'multiObjective' list");
  }
}

void ObjectiveContainer::PostProc(DesignSpace* space, DesignStructure* structure, MultipleExcitation* me)
{
  for(Objective* o : data)
  {
    assert(o->HasDenseJacobian());
    if(o->tune.IsSet())
      o->tune.Register(&(o->scale), domain->GetOptimization());
    o->SetDenseSparsityPattern(space);
    o->SetElements(space, o->region); // before Function::PostProc() !
    o->PostProc(space, structure);
    o->SetExcitation(me);
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
      f->ToInfo(o); // has scale with sensitivity to tune
    }
  }
  else
  {
    data[0]->ToInfo(in);
    if(data[0]->scale != 1.0) // only when it is set
      in->Get("scale")->SetValue(data[0]->scale);
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
  for(Objective* o : data)
    if(o->GetType() == type)
      return o;

  if(throw_exception)
    Exception("No objective of type '" + Objective::type.ToString(type) + "' stored");
  return nullptr;
}

Objective* ObjectiveContainer::Get(const std::string& name, bool throw_exception)
{
  for(Objective* o : data)
    if(o->ToString() == name)
      return o;

  if(throw_exception)
    Exception("No objective '" + name + "' stored");

  return nullptr;
}

double ObjectiveContainer::GetHistoryValue(bool penalty, int index)
{
  double result = 0.0;

  Vector<double> vals(data.GetSize());

  int idx = index >= 0 ? index : (int) GetHistorySize() + index; // -1 is last, ...
  LOG_DBG(obj) << "OC:GHV p=" << penalty << " index=" << index << " size=" << data.GetSize() << " idx=" << idx;
  assert(idx >= 0);
  for(unsigned int i = 0; i < data.GetSize(); i++)
  {
    double val = data[i]->history[idx];
    vals[i] = (penalty ? data[i]->scale : 1.0) * val;
    result += vals[i];
  }

  double beta;
  if(domain->GetOptimization()->GetMOType(beta) == Function::SMOOTH_MIN)
    result = SmoothMin(vals, beta);
  else if(domain->GetOptimization()->GetMOType(beta) == Function::SMOOTH_MAX)
    result = SmoothMax(vals, beta);

  return result;
}

unsigned int ObjectiveContainer::GetHistorySize()
{
  return data[0]->history.GetSize();
}

void ObjectiveContainer::PushBackHistory()
{
  for(Objective* f : data)
    f->history.Push_back(f->value_);
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
