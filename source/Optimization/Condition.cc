#include "Optimization/Condition.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Design/DesignStructure.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "General/exception.hh"
#include "General/environment.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"
#include "MatVec/denseMatrix.hh"
#include <sstream>

using namespace CoupledField;


DECLARE_LOG(conditions)

// instantiation of the static elements
Enum<Condition::Bound> Condition::bound;

Condition::Condition(PtrParamNode pn) : Function(pn)
{
  volume_fraction = 0.0;
  special_result_idx = -1;
  blown_up_ = false;
  index_ = -1; // to be set by ConditionContainer::Read()
  virtual_base_index_ = -1;
  bound_ = bound.Parse(pn->Get("bound")->As<std::string>());
  design = !pn->Has("design") ? DesignElement::DEFAULT :
           DesignElement::type.Parse(pn->Get("design")->As<std::string>());

  // the bound value is called value in the problem file!
  // there must not  be a value when a homogenization tensor is given
  this->boundValue_ = pn->Has("value") ? pn->Get("value")->As<double>() : -1.0;

  // special handling of scaling
  objective_scaling_ = pn->Get("scaling")->As<std::string>() == "objective";
  manual_scaling_value = objective_scaling_ ? -1.0 : pn->Get("scaling")->As<double>();

  observation_ = pn->Get("mode")->As<std::string>() == "observation";

  delta_logging_ignored_ = false;
  delta_logging = pn->Get("log_delta")->As<bool>();

  if(pn->Get("log_delta")->As<bool>() && (!pn->Has("value") && !pn->Has("tensor") && !pn->Has("isotropic")))
    delta_logging_ignored_ = true;

  if(pn->Has("coord"))
    ReadCoord(pn);

  penalty = pn->Get("penalty")->As<double>();
  region = ALL_REGIONS;

  if(pn->Has("region") && pn->Get("region")->As<std::string>() != "all")
    region = domain->GetGrid()->GetRegion().Parse(pn->Get("region")->As<std::string>());

  // value is not mandatory for all almost all constraints. Check for homogenization later
  if(!observation_)
  {
    switch(type_)
    {
    case HOMOGENIZATION_TENSOR:
    case HOMOGENIZATION_TRACKING:
      if(!pn->Has("tensor") && !pn->Has("value") && !pn->Has("isotropic"))
        throw Exception("Neither value nor tensor is given for constraint '" + type.ToString(type_) + "'");
      break;

    case ISOTROPY:
      if(pn->Has("value"))
        throw Exception("No value allowed for constraint '" + type.ToString(type_) + "'");
      break; // ok without value

    default:
      if(!pn->Has("value"))
        throw Exception("No value given for constraint '" + type.ToString(type_) + "'");
    }
  }
  
  linear_ = type_ == VOLUME || type_ == SLOPE ? true : false;
  //  snopt only makes a difference between linear and nonlinear constraints!
  if(pn->Has("linear"))
    linear_ = pn->Get("linear")->As<bool>();

  // check bound
  switch(type_)
  {
  case SLOPE:
  case GLOBAL_SLOPE:
  case CHECKERBOARD:
  case OSCILLATION:
    if(bound_ != UPPER_BOUND && IsActive())
      bound_ = UPPER_BOUND; // detected in PostProc() and warning is given
    break;

  default:
    break;
  }
}

void Condition::PostProc(DesignSpace* space, DesignStructure* structure)
{
  // note, meanwhile we have info_ set! but not yet in the constructor
  Function::PostProc(space, structure);

  // did we change the bound?
  switch(type_)
  {
  case SLOPE:
  case GLOBAL_SLOPE:
  case CHECKERBOARD:
  case OSCILLATION:
    if(pn->Get("bound")->As<std::string>() != bound.ToString(bound_) && IsActive())
      info_->Get(ParamNode::WARNING)->SetValue("changed bound for '" + type.ToString(type_) + "' to 'upperBound'");
    break;

  default:
    break;
  }
}

bool Condition::ReadCoord(PtrParamNode pn)
{
  std::string val = pn->Get("coord")->As<std::string>();
  if(val == "all") return false;

  assert(val.size() == 2);
  coords.Resize(1);

  ParseCoord(pn, coords[0]);

  return true;
}



void Condition::AddCondition(PtrParamNode pn, StdVector<Condition*>& list)
{
  Type t = type.Parse(pn->Get("type")->As<std::string>());
  list.Push_back(t == SLOPE || t == MOLE || t == CHECKERBOARD || t == OSCILLATION ? new LocalCondition(pn) : new Condition(pn));

  // note that the pointer becomes invalid by AddSubCondition()
  Condition* g = list.Last();
  g->index_ = -1; // to be defined by the virtual "all" list in the ConditionContainer::Read() method

  // homogenization are special constraints which constraints which "blow up" to more constraints
  if(g->type_ == HOMOGENIZATION_TENSOR)
  {
    // there has been a extensive test in the constructor
    // do we need to blow-up?
    if(!g->ReadCoord(pn)) // this is done if coord="all" is given in xml!
    {
      // the first entry is this constraint
      // for the conversion of the indices see Thesis from Ole Sigmund, p 30 and book of Manfred
      // p 42 (first edition)
      g->coords.Resize(1);
      g->coords[0].first = 1;      // 1. ijkl = 1111
      g->coords[0].second = 1;
      g->boundValue_ = g->tensor_[1-1][1-1]; // one-based!
      g = g->AppendSubCondition(list, 2,2); // 2. 2222
      if(domain->GetGrid()->GetDim() == 2)
      {
        g = g->AppendSubCondition(list, 1,2); // 3. 1122

        //g = g->AppendSubCondition(list, 3,1); // no covered by Sigmund
        //g = g->AppendSubCondition(list, 3,2); // no covered by Sigmund

        g = g->AppendSubCondition(list, 3,3); // 4. 1212 // would be 6,6 according to book
        
        g = g->AppendSubCondition(list, 1,3); // 4. 1122
        g = g->AppendSubCondition(list, 2,3); // 5. 2233
      }
      else
      {
        g = g->AppendSubCondition(list, 3,3); // 3. 3333
        g = g->AppendSubCondition(list, 1,2); // 4. 1122
        g = g->AppendSubCondition(list, 1,3); // 5. 1133
        g = g->AppendSubCondition(list, 2,3); // 6. 2233
        g = g->AppendSubCondition(list, 6,6); // 7. 1212
        g = g->AppendSubCondition(list, 5,5); // 8. 1313
        g = g->AppendSubCondition(list, 4,4); // 9. 2323
        
        g = g->AppendSubCondition(list, 1,4);
        g = g->AppendSubCondition(list, 1,5);
        g = g->AppendSubCondition(list, 1,6);
        g = g->AppendSubCondition(list, 2,4);
        g = g->AppendSubCondition(list, 2,5);
        g = g->AppendSubCondition(list, 2,6);
        g = g->AppendSubCondition(list, 3,4);
        g = g->AppendSubCondition(list, 3,5);
        g = g->AppendSubCondition(list, 3,6);
        
        g = g->AppendSubCondition(list, 4,5);
        g = g->AppendSubCondition(list, 4,6);
        g = g->AppendSubCondition(list, 5,6);
      }
    }
  }
  // isotropy is a special constraint which blows up special tensor entry constraints
  if(g->type_ == ISOTROPY)
  {
    if(pn->Has("coord"))
      throw Exception("don't use attribute 'coord' for constraint 'isotropy'");

    if(g->bound_ != EQUAL)
      throw Exception("the 'isotropy' constraint requires equality constraint type");

    // become an HOMOGENIZATION_TENSOR constraint!
    g->type_ = HOMOGENIZATION_TENSOR;

    if(domain->GetGrid()->GetDim() == 2)
    {
      // E11 - E22 = 0
      assert(g->coords.GetSize() == 0);
      g->boundValue_ = 0;
      g->coords.Push_back(std::make_pair(1,1));
      g->coords.Push_back(std::make_pair(2,2));

      // E11 - E12 - 2E33 = E11 - E12 - E33 - E33 = 0
      g = g->AppendSubCondition(list);
      g->coords.Clear(); // the copy constructor above copies old stuff
      g->coords.Push_back(std::make_pair(1,1));
      g->coords.Push_back(std::make_pair(1,2));
      g->coords.Push_back(std::make_pair(3,3));
      g->coords.Push_back(std::make_pair(3,3));

      // E13 = 0
      g = g->AppendSubCondition(list);
      g->coords.Clear();
      g->coords.Push_back(std::make_pair(1,3));

      // E23 = 0
      g = g->AppendSubCondition(list);
      g->coords.Clear();
      g->coords.Push_back(std::make_pair(2,3));
    }
    else
    {
      // non-shear diagonal is constant
      // E11 = E22 = E33 -> E11 - E22 = 0, E22 - E33 = 0
      assert(g->coords.GetSize() == 0);
      g->boundValue_ = 0;
      g->coords.Push_back(std::make_pair(1,1));
      g->coords.Push_back(std::make_pair(2,2));

      g = g->AppendSubCondition(list);
      g->coords.Clear();
      g->coords.Push_back(std::make_pair(2,2));
      g->coords.Push_back(std::make_pair(3,3));

      // upper non-shear triangle is constant
      // E12 = E13 = E23 -> E12 - E13 = 0, E13 - E23 = 0
      g = g->AppendSubCondition(list);
      g->coords.Clear();
      g->coords.Push_back(std::make_pair(1,2));
      g->coords.Push_back(std::make_pair(1,3));

      g = g->AppendSubCondition(list);
      g->coords.Clear();
      g->coords.Push_back(std::make_pair(1,3));
      g->coords.Push_back(std::make_pair(2,3));

      // shear diagonal is constant
      // E44 = E55 = E66 -> E44 - E55 = 0, E55 - E66 = 0
      g = g->AppendSubCondition(list);
      g->coords.Clear();
      g->coords.Push_back(std::make_pair(4,4));
      g->coords.Push_back(std::make_pair(5,5));

      g = g->AppendSubCondition(list);
      g->coords.Clear();
      g->coords.Push_back(std::make_pair(5,5));
      g->coords.Push_back(std::make_pair(6,6));

      // relationship of the three unique values
      // E11 - E12 - 2E66 = E11 - E12 - E66 - E66 = 0
      g = g->AppendSubCondition(list);
      g->coords.Clear(); // the copy constructor above copies old stuff
      g->coords.Push_back(std::make_pair(1,1));
      g->coords.Push_back(std::make_pair(1,2));
      g->coords.Push_back(std::make_pair(6,6));
      g->coords.Push_back(std::make_pair(6,6));

      // the rest is zero
      // E14 = E15 = E16 = E24 = E25 = E26 = E34 = E35 = E36 = E45 = E46 = E56 = 0
      g = g->AppendSubCondition(list);
      g->coords.Clear();
      g->coords.Push_back(std::make_pair(1,4));

      g = g->AppendSubCondition(list);
      g->coords.Clear();
      g->coords.Push_back(std::make_pair(1,5));

      g = g->AppendSubCondition(list);
      g->coords.Clear();
      g->coords.Push_back(std::make_pair(1,6));

      g = g->AppendSubCondition(list);
      g->coords.Clear();
      g->coords.Push_back(std::make_pair(2,4));

      g = g->AppendSubCondition(list);
      g->coords.Clear();
      g->coords.Push_back(std::make_pair(2,5));

      g = g->AppendSubCondition(list);
      g->coords.Clear();
      g->coords.Push_back(std::make_pair(2,6));

      g = g->AppendSubCondition(list);
      g->coords.Clear();
      g->coords.Push_back(std::make_pair(3,4));

      g = g->AppendSubCondition(list);
      g->coords.Clear();
      g->coords.Push_back(std::make_pair(3,5));

      g = g->AppendSubCondition(list);
      g->coords.Clear();
      g->coords.Push_back(std::make_pair(3,6));

      g = g->AppendSubCondition(list);
      g->coords.Clear();
      g->coords.Push_back(std::make_pair(4,5));

      g = g->AppendSubCondition(list);
      g->coords.Clear();
      g->coords.Push_back(std::make_pair(4,6));

      g = g->AppendSubCondition(list);
      g->coords.Clear();
      g->coords.Push_back(std::make_pair(5,6));
    }
  }
}

Condition* Condition::AppendSubCondition(StdVector<Condition*>& list)
{
  list.Push_back(new Condition(*this)); // make a copy of this element by the (default) copy constructor
  Condition* sub = list.Last(); // copy this entry as reference
  sub->index_ = list.GetSize() - 1;
  sub->blown_up_ = true;
  return sub;
}


Condition* Condition::AppendSubCondition(StdVector<Condition*>& list, int pos_x, int pos_y)
{
  Condition* sub = AppendSubCondition(list);
  sub->coords.Resize(1);
  sub->coords[0].first = pos_x;
  sub->coords[0].second = pos_y;
  sub->boundValue_ = sub->tensor_[pos_x-1][pos_y-1];
  sub->blown_up_ = true;
  return sub;
}

void Condition::SetDenseSparsityPattern(DesignSpace* space)
{
  unsigned int size = space->GetNumberOfVariables();
  sparsity_.Resize(size);
  
  for(unsigned int i = 0; i < size; i++)
    sparsity_[i] = i;
}

StdVector<unsigned int>& Condition::GetSparsityPattern()
{
  assert(!sparsity_.IsEmpty()); // SetSparsityPattern() needs to be called before
  
  return sparsity_;
}

std::string Condition::ToString() const
{
  std::ostringstream os;
  
  if(delta_logging) os << "delta_";
  os << Function::ToString(); // includes physical
  if(design != DesignElement::DEFAULT)
    os << "_(" << DesignElement::type.ToString(design) << ")";
  if(type_ == HOMOGENIZATION_TENSOR)
    for(unsigned int i = 0; i < coords.GetSize(); i++)
      os << "_" << coords[i].first << coords[i].second;
  return os.str();  
}


void Condition::ToInfo(PtrParamNode in)
{
  Function::ToInfo(in);

  in->Get("mode")->SetValue(observation_ ? "observation" : "constraint");
  in->Get("design")->SetValue(DesignElement::type.ToString(design));
  if(IsActive())
  {
    in->Get("bound")->SetValue(bound.ToString(bound_));
    if(type_ != HOMOGENIZATION_TRACKING)
      in->Get("bound_value")->SetValue(boundValue_);
  }
  if(type_ == HOMOGENIZATION_TENSOR)
  {
    std::ostringstream os;
    for(unsigned int i = 0; i < coords.GetSize(); i++)
      os << (i > 0 ? "_" : "") << coords[i].first << coords[i].second;
    in->Get("tensor_entry")->SetValue(os.str());
  }

  if(delta_logging_ignored_)
    in->Get("delta_logging")->Get(ParamNode::WARNING)->SetValue("no value given");
  else
    in->Get("delta_logging")->SetValue(delta_logging);

  // TODO somehow scaling does not work ??
  // if(IsHomogenization() && !objective_scaling_ && !blown_up_) // warn only the first time!
  //  in->Get(ParamNode::WARNING)->SetValue("Doing homogenization without 'objective' scaling constraint '" + type.ToString(type_) + "'");

  if(!observation_)
    in->Get("linear")->SetValue(linear_);
}

bool Condition::IsForRegion(RegionIdType regionId)
{
  return(region == ALL_REGIONS || region == regionId);
}

LocalCondition::LocalCondition(PtrParamNode pn) : Condition(pn)
{
  current_view_index_ = -1;
}


Function::Local::Identifier& LocalCondition::GetCurrentVirtualContext()
{
  assert(IsLocal());

  unsigned int idx = current_view_index_ - virtual_base_index_;
  return local->virtual_elem_map[idx];
}

StdVector<unsigned int>& LocalCondition::GetSparsityPattern()
{
  /* for debug purposes you can enable the dense pattern by commenting out
   * the two lines */
  /*
  SetDenseSparsityPattern(space_);
  return sparsity_;
  */
  
  assert(IsLocal());

  Function::Local::Identifier& id = GetCurrentVirtualContext();

    // we shall sort the indices
  std::list<unsigned int> indices;
  for(int i = -1 ; i < (int) id.neighbor.GetSize(); i++)
  {
    DesignElement* de = id.GetElement(i);
    int other_idx = local->space->Find(de); // needs to be fast!
    indices.push_back(other_idx);
  }

  // sort and copy
  indices.sort();
  sparsity_.Resize(0); // keeps capacity, hence Push_back is cheap
  for(std::list<unsigned int>::const_iterator it = indices.begin(); it != indices.end(); ++it)
    sparsity_.Push_back(*it);

  LOG_DBG3(conditions) << "LC:GSP: current_view_index_=" << current_view_index_ << " -> " << sparsity_.ToString();
  return sparsity_;
}

double LocalCondition::CalcMeanValue() const
{
  double sum = 0.0;
  for(unsigned int i = 0, n = local->virtual_elem_map.GetSize(); i < n; i++)
    sum += std::abs(local->virtual_elem_map[i].EvalFunction(local));
  return sum / local->virtual_elem_map.GetSize();
}

double LocalCondition::CalcMaxValue() const
{
  double max = 0.0;
  for(unsigned int i = 0, n = local->virtual_elem_map.GetSize(); i < n; i++)
  {
    double v = std::abs(local->virtual_elem_map[i].EvalFunction(local));
    max = std::max(max, v);
  }

  return max;
}


double LocalCondition::GetValue() const
{
  if(IsLocal())
    return local->values[current_view_index_ - virtual_base_index_];
  else
    assert(false);
    return value_;
}

void LocalCondition::SetValue(double val)
{
  if(IsLocal())
  {
    local->values[current_view_index_ - virtual_base_index_] = val;
    value_ = -1.0; // invalidated
  }
  else
    value_ = val; // only for Done() allowed!!
}


std::string LocalCondition::ToString() const
{
  std::stringstream ss;

  ss << Condition::ToString();

  if(IsLocal())
    ss << " cvi=" << current_view_index_;

  return ss.str();
}


ConditionContainer::ConditionContainer()
{
  view = NULL;
  space_ = NULL;
  // set in Read()
}

ConditionContainer::~ConditionContainer()
{
  delete(view);
}

void ConditionContainer::Read(ParamNodeList pn_list)
{
  // call only once
  assert(all.IsEmpty());

  // slope constraints need to be post processed in ErsatzMaterial
  for(unsigned int i = 0; i < pn_list.GetSize(); i++)
  {
    PtrParamNode pn = pn_list[i];
    bool act = pn->Get("mode")->As<std::string>() == "constraint";
    Condition::AddCondition(pn, act ? active : observe);
  }

  // process the virtual containers
  all.Reserve(active.GetSize() + observe.GetSize());

  for(unsigned int i = 0; i < active.GetSize(); i++)
    all.Push_back(active[i]);

  for(unsigned int i = 0; i < observe.GetSize(); i++)
    all.Push_back(observe[i]);

  // set index
  for(unsigned int i = 0; i < all.GetSize(); i++)
    all[i]->index_ = i;

  view = new VirtualView(this);
}


void ConditionContainer::PostProc(DesignSpace* space, DesignStructure* structure)
{
  this->space_ = space;

  // the conditions have no space 
  for(unsigned int i = 0; i < all.GetSize(); i++)
    if(all[i]->HasDenseJacobian())
      all[i]->SetDenseSparsityPattern(space);
  
  // check for special result index if there was a <result> for value=constraintGradient
  for(unsigned int i = 0; i < all.GetSize(); i++)
  {
    Condition* g = all[i];

    // the strings for Function::Type are (partially) repeated as DesignElement::Detail
    std::string constr_str = g->type.ToString(g->type_); // our type as string
    if(DesignElement::detail.IsValid(constr_str)) // is it defined for output?
    {
      DesignElement::Detail detail = DesignElement::detail.Parse(constr_str);
      g->special_result_idx = space->GetSpecialResultIndex(g->design, DesignElement::CONSTRAINT_GRADIENT, detail);
    } // -1 for else by default in constructor
  }

  for(unsigned int i = 0; i < all.GetSize(); i++)
  {
    all[i]->PostProc(space, structure);
  }

  view->Refresh(); // inform about the news if the slopes created a lot of virtual objectives!
}

void ConditionContainer::ToInfo(PtrParamNode in)
{
  for(unsigned int i = 0; i < all.GetSize(); i++)
    all[i]->ToInfo(in->Get("constraint", ParamNode::APPEND));
}

StdVector<Condition*> ConditionContainer::GetList(Condition::Type type, DesignElement::Type design, bool only_active)
{
  StdVector<Condition*> result;

  for(unsigned int i = 0; i < active.GetSize(); i++)
    if(active[i]->GetType() == type && (design != DesignElement::NO_TYPE ? active[i]->design == design : true))
      result.Push_back(active[i]);

  for(unsigned int i = 0; !only_active && i < observe.GetSize(); i++)
    if(observe[i]->GetType() == type && (design != DesignElement::NO_TYPE ? observe[i]->design == design : true))
      result.Push_back(observe[i]);

  return result;
}

bool ConditionContainer::Has(Condition::Type type, DesignElement::Type design)
{
  // be save and check for uniqueness!
  StdVector<Condition*> list = GetList(type, design, true);

  return !list.IsEmpty();
}

Condition* ConditionContainer::Get(Condition::Type type, DesignElement::Type design, bool only_active, bool throw_exception)
{
  // be save and check for uniqueness!
  StdVector<Condition*> list = GetList(type, design, only_active);

  if(list.GetSize() == 0)
  {
    if(throw_exception)
      throw Exception("no active constraint '" + Condition::type.ToString(type) + "' found");
    else
      return NULL;
  }

  if(list.GetSize() > 1 && throw_exception)
    throw Exception("constraint " + Condition::type.ToString(type) + "is not unique");

  return list[0];
}

ConditionContainer::VirtualView::VirtualView(ConditionContainer* constraints)
{
  container_ = constraints;
  Refresh();
}

void ConditionContainer::VirtualView::Refresh()
{
  // find the indices of LocalConditions, they need to be sorted.
  std::list<unsigned int> tmp;

  // search also for observe conditions!
  Condition* c = container_->Get(Condition::SLOPE, DesignElement::NO_TYPE, false, false);
  if(c != NULL) tmp.push_back(c->GetIndex());
  c = container_->Get(Condition::CHECKERBOARD, DesignElement::NO_TYPE, false, false);
  if(c != NULL) tmp.push_back(c->GetIndex());
  c = container_->Get(Condition::OSCILLATION, DesignElement::NO_TYPE, false, false);
  if(c != NULL) tmp.push_back(c->GetIndex());
  c = container_->Get(Condition::MOLE, DesignElement::NO_TYPE, false, false);
  if(c != NULL) tmp.push_back(c->GetIndex());

  tmp.sort();

  // copy sorted
  local_cond_index_.Resize(0);
  for(std::list<unsigned int>::iterator it = tmp.begin(); it != tmp.end(); ++it)
    local_cond_index_.Push_back(*it);

  // determine the virtual sizes for the container
  virtual_active_size_ = container_->active.GetSize();
  virtual_total_size_ = container_->all.GetSize();

  for(unsigned int i = 0; i < local_cond_index_.GetSize(); i++)
  {
    LocalCondition* lc = dynamic_cast<LocalCondition*>(container_->all[local_cond_index_[i]]);

    // replace the global slope by many local slopes -> if it is initialized!
    if(lc != NULL && lc->GetLocal() != NULL && lc->IsActive() && lc->GetConstraintSize() > 0)
      virtual_active_size_ += lc->GetConstraintSize() -1; // replace means remove ourselves

    virtual_total_size_ = container_->all.GetSize();
    if(lc != NULL && lc->GetLocal() != NULL && lc->GetConstraintSize() > 0)
      virtual_total_size_ += lc->GetConstraintSize() -1;
  }

  // set the virtual base indices
  int curr = 0;
  for(unsigned int i = 0; i < container_->all.GetSize(); i++)
  {
    Condition* g = container_->all[i];
    g->virtual_base_index_ = curr;
    if(g->IsLocalCondition() && g->GetLocal() != NULL) // does not need to be initialized yet
      curr += std::max((int) dynamic_cast<LocalCondition*>(g)->GetConstraintSize(), 1);
    else
      curr++;
  }
  assert(curr == virtual_total_size_);
}

Condition* ConditionContainer::VirtualView::Get(int view_index)
{
  StdVector<Condition*>& all = container_->all;;
  assert(all.GetSize() > 0);

  // traverse the conditions, if we are above virtual_base_index we gone one too far
  Condition* g = NULL;
  for(unsigned int i = 0; g == NULL && i < all.GetSize()-1; i++)
    if(all[i+1]->virtual_base_index_ > view_index) // we are right if the next is too far
      g = all[i];
  if(g == NULL) g = all.Last();
  assert(g->virtual_base_index_ <= view_index);

  if(g->IsLocalCondition())
  {
    dynamic_cast<LocalCondition*>(g)->SetCurrentViewIndex(view_index);
  }

  return g;
}

void ConditionContainer::VirtualView::Done()
{
  for(unsigned int g = 0; g < local_cond_index_.GetSize(); g++) // no local conditions, nothing to do
  {
    LocalCondition* lc = dynamic_cast<LocalCondition*>(container_->all[local_cond_index_[g]]);
    lc->SetCurrentViewIndex(-1); // reset to global mode

    // shall we give the values as special result?
    DesignElement::ValueSpecifier vs = DesignElement::MAX_CHECKERBOARD; // overwrite if necessary
    if(lc->GetType() == Function::OSCILLATION)  vs = DesignElement::MAX_OSCILLATION;
    if(lc->GetType() == Function::SLOPE)        vs = DesignElement::MAX_SLOPE;
    if(lc->GetType() == Function::MOLE)         vs = DesignElement::MAX_MOLE;

    int idx = container_->space_->GetSpecialResultIndex(DesignElement::DEFAULT, vs);
    if(idx >= 0)
    {
      StdVector<DesignElement>& des_data = container_->space_->data;

      // we add up the max value and not elements have a slope constraint, therefore reset
      for(unsigned int e = 0; e < des_data.GetSize(); e++)
        des_data[e].specialResult[idx] = 0.0; // initialize

      StdVector<Function::Local::Identifier>& vem = lc->GetLocal()->virtual_elem_map;

      for(unsigned int i = 0; i < vem.GetSize(); i++)
      {
        Function::Local::Identifier& id = vem[i];
        DesignElement* de =  id.element;
        double sv = id.EvalFunction(lc->local);

        // in checkerboard we must not use abs
        double corr = lc->GetType() == Function::CHECKERBOARD  || lc->GetType() == Function::OSCILLATION ? sv : std::abs(sv);

        de->specialResult[idx] = std::max(de->specialResult[idx], corr);
      }
    }
  }
}
