#include "Optimization/Objective.hh"
#include "Optimization/Condition.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "General/exception.hh"
#include "Domain/domain.hh"

using namespace CoupledField;

// instantiation of the static elements
Enum<Objective::Type> Objective::type;

Objective::Objective(PtrParamNode pn, PtrParamNode pn_type, unsigned int idx)
{
  // multiple excitation is handled in Optimization itself!

  // the current value -> check <Get/Set>Value() when altering the presets!
  this->index_       = idx;
  this->value_       = -1.0;
  this->pn           = pn;
  this->harmonic_    = BasePDE::IsComplex(domain->GetDriver()->GetAnalysisType());
  this->omega_omega_ = pn->Has("factor") ? pn->Get("factor")->Get("omega_omega")->As<bool>() : false;
  this->volumePenaltyExponent = 1.0;
  if(!harmonic_ && omega_omega_)
    throw Exception("It makes no sense to set costFunction/factor/omega_omega in static optimization");

  this->type_ = type.Parse(pn_type->Get("type")->As<std::string>());
  this->penalty_ = pn_type->Has("penalty") ? pn_type->Get("penalty")->As<Double>() : 1.0;

  if(type_ == HOMOGENIZATION_TRACKING)
    if(!Condition::ReadTensor(pn, tensor))
      EXCEPTION("A 'tensor' element is mandatory  for 'homogenizationTracking'");

  // currently we have only relative implemented up to now!
  stop.value = pn->Has("stopping") ? pn->Get("stopping/value")->As<Double>() : 1e-5;
  stop.queue = pn->Has("stopping") ? pn->Get("stopping/queue")->As<Integer>() : 5; // is >= 1!

  // set how often to evaluate the objective for multiple excitations
  switch(type_)
  {
    case VOLUME:
      if(pn->Has("volumePenaltyExponent"))
        volumePenaltyExponent = pn->Get("volumePenaltyExponent")->Get("value")->As<Double>();
      // we do not break here!
      // because we also need the evaluateOnce_ value to be set!
    case HOMOGENIZATION_TENSOR:
    case HOMOGENIZATION_TRACKING:
    case HOMOGENIZATION_E11:
    case POISSONS_RATIO:
    case YOUNGS_MODULUS:
    case TYCHONOFF:
      this->evaluateOnce_ = true;
      break;

    case MULTI_OBJECTIVE: // only to make the switch complete
    case COMPLIANCE:
    case OUTPUT:
    case DYNAMIC_OUTPUT:
    case TRACKING:
    case ABS_DYN_OUTPUT_SQUARED:
    case CONJUGATE_COMPLIANCE:
    case GLOBAL_DYNAMIC_COMPLIANCE:
    case ELEC_ENERGY:
    case TEMPERATURE:
      this->evaluateOnce_ = false; // standard case
  }
}

double Objective::GetValue() const
{
  return value_;
}


void Objective::SetValue(double val)
{
  value_ = val;
}

bool Objective::IsHomogenization() const
{
  switch(type_)
  {
    case HOMOGENIZATION_TENSOR:
    case HOMOGENIZATION_TRACKING:
    case HOMOGENIZATION_E11:
    case POISSONS_RATIO:
    case YOUNGS_MODULUS:
      return true;

    default:
      return false;
  }
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
    data.Resize(1);
    data[0] = new Objective(obj_node, obj_node, 0);
  }
  else
  {
    if(!obj_node->Has("multiObjective"))
      throw Exception("For costFunction type 'multiObjective' element with 'objective' childs is needed");

    ParamNodeList list = obj_node->Get("multiObjective")->GetList("objective");
    if(list.GetSize() == 0)
      throw Exception("For costFunction type 'multiObjective' element with 'objective' childs is needed");

    data.Resize(list.GetSize());

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
      o->Get("type")->SetValue(Objective::type.ToString(data[i]->type_));
      o->Get("penalty")->SetValue(data[i]->penalty_);
      o->Get("evaluate")->SetValue(data[i]->evaluateOnce_ ? "once" : "always");
    }
  }
  else
  {
    in->Get("type")->SetValue(Objective::type.ToString(data[0]->type_));
    in->Get("evaluate")->SetValue(data[0]->evaluateOnce_ ? "once" : "always");
  }

  in->Get("task")->SetValue(minimize_ ? "minimize" : "maximize");
  if(data[0]->harmonic_)
    in->Get("factor")->Get("omega_omega")->SetValue(data[0]->omega_omega_);
  if(data[0]->tensor.GetNumRows() > 1)
    in->Get("tensor")->SetValue(new Matrix<double>(data[0]->tensor));
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

