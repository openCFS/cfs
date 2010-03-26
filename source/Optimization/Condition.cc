#include "Optimization/Condition.hh"
#include "Optimization/DesignSpace.hh"
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
  linear_ = false;
  index_ = -1; // to be set by ConditionContainer::Read()
  bound_ = bound.Parse(pn->Get("bound")->As<std::string>());
  design = !pn->Has("design") ? DesignElement::DEFAULT :
           DesignElement::type.Parse(pn->Get("design")->As<std::string>());

  // the bound value is called value in the problem file!
  // there must not  be a value when a homogenization tensor is give->As<std::string>()n
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
      if(!pn->Has("tensor") && !pn->Has("value"))
        throw Exception("Neither value nor tensor is given for constraint '" + type.ToString(type_) + "'");
      break;

    case ISOTROPY:
      if(pn->Has("value"))
        throw Exception("No value allowed for constraint '" + type.ToString(type_) + "'");
      break; // ok without value

    case CHECKERBOARD:
      throw Exception("'checkerboard' may only be used in 'observe' mode");

    default:
      if(!pn->Has("value"))
        throw Exception("No value given for constraint '" + type.ToString(type_) + "'");
    }
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
  list.Push_back(t == SLOPE ? new SlopeCondition(pn) : new Condition(pn));

  // note that the pointer becomes invalid by AddSubCondition()
  Condition* g = list.Last();
  g->index_ = -1; // to be defined by the virtual "all" list in the ConditionContainer::Read() method

  // homogenization are special constraints which constraints which "blow up" to more constraints
  if(g->type_ == HOMOGENIZATION_TENSOR)
  {
    // there has been a extensive test in the constructor
    // do we need to blow-up?
    if(!g->ReadCoord(pn))
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


void Condition::ToInfo(PtrParamNode in) const
{
  in->Get("type")->SetValue(type.ToString(type_));
  in->Get("design")->SetValue(DesignElement::type.ToString(design));
  if(IsActive())
  {
    in->Get("bound")->SetValue(bound.ToString(bound_));
    if(type_ != HOMOGENIZATION_TRACKING && type_ != SLOPE)
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

  if(IsHomogenization() && !objective_scaling_ && !blown_up_) // warn only the first time!
    in->Get(ParamNode::WARNING)->SetValue("Doing homogenization without 'objective' scaling constraint '" + type.ToString(type_) + "'");
}

bool Condition::IsForRegion(RegionIdType regionId)
{
  return(region == ALL_REGIONS || region == regionId);
}

SlopeCondition::SlopeCondition(PtrParamNode pn) : Condition(pn)
{
  current_view_index_ = -1;
  space_ = NULL;
  dim_ = domain->GetGrid()->GetDim();
  linear_ = true;

  if(bound_ != UPPER_BOUND && IsActive())
    throw Exception("Slope constraint needs to be 'upper' bound");

  if(pn->Has("parameter")) // set in base class
    throw Exception("Parameter for slope constraint specified, but not required");
}

void SlopeCondition::PostInit(DesignSpace* space)
{
  space_ = space;

  virtual_elem_map_.Reserve(2 * space->GetNumberOfElements() * dim_);

  // traverse all elements and check for full neighborhood
  int base  = space->FindDesign(design);
  int elems = space->GetNumberOfElements();
  for(int e = base * elems; e < (base + 1) * elems; e++)
  {
    DesignElement& de = space->data[e];

    // do we have a full neighborhood? All or none as in the original paper
    bool full = true;
    if(de.vicinity->design[VicinityElement::X_P] == NULL) full = false;
    if(de.vicinity->design[VicinityElement::Y_P] == NULL) full = false;
    if(dim_ == 3 && de.vicinity->design[VicinityElement::Z_P] == NULL) full = false;

    LOG_DBG2(conditions) << "SC:PI e_num=" << de.elem->elemNum << " vicinity=" << de.vicinity->ToString() << " full=" << full;

    if(full)
    {
      virtual_elem_map_.Push_back(Identifier(e, VicinityElement::X_P, 1));
      virtual_elem_map_.Push_back(Identifier(e, VicinityElement::X_P, -1));
      virtual_elem_map_.Push_back(Identifier(e, VicinityElement::Y_P, 1));
      virtual_elem_map_.Push_back(Identifier(e, VicinityElement::Y_P, -1));
      if(dim_ == 3) {
        virtual_elem_map_.Push_back(Identifier(e, VicinityElement::Z_P, 1));
        virtual_elem_map_.Push_back(Identifier(e, VicinityElement::Z_P, -1));
      }
    }
  }

  // needs to be set prior CalcSlopeConstraint() as the optimizers need the size
  values_.Resize(virtual_elem_map_.GetSize(), 0.0);
}

unsigned int SlopeCondition::GetCurrentVirtualElement() const
{
  assert(IsLocal());

  unsigned int idx = current_view_index_ - index_;

  LOG_DBG3(conditions) << "SC:curVirEle: cvi=" << current_view_index_ << " -> idx=" << virtual_elem_map_[idx].element_idx;

  return virtual_elem_map_[idx].element_idx;
}

int SlopeCondition::GetCurrentVirtualNeighbor() const
{
  assert(IsLocal());

  unsigned int idx = current_view_index_ - index_;

  LOG_DBG3(conditions) << "SC:curVirNeigh: cvi =" << current_view_index_ << " -> neigh=" << virtual_elem_map_[idx].neighbor;

  return virtual_elem_map_[idx].neighbor;
}

int SlopeCondition::GetCurrentVirtualSign() const
{
  assert(IsLocal());

  unsigned int idx = current_view_index_ - index_;

  return virtual_elem_map_[idx].sign;
}

StdVector<unsigned int>& SlopeCondition::GetSparsityPattern()
{
  /* for debug purposes you can enable the dense pattern by commenting out
   * the two lines */
  /*
  SetDenseSparsityPattern(space_);
  return sparsity_;
  */
  
  assert(IsLocal());

  sparsity_.Resize(2);

  // we have two entries. sort them
  int own_idx = GetCurrentVirtualElement();
  DesignElement* other_elem = space_->data[own_idx].vicinity->design[GetCurrentVirtualNeighbor()];
  int other_idx = space_->Find(other_elem);

  // as we have X_P, Y_P and Z_P neighbors the sorting could be right but better invest 1e-9 sec
  sparsity_[0] = own_idx < other_idx ? own_idx : other_idx;
  sparsity_[1] = own_idx > other_idx ? own_idx : other_idx;
  
  LOG_DBG3(conditions) << "SC:GetSparPat: current_view_index_=" << current_view_index_ 
                       << " g[0]=" << (sparsity_.IsEmpty() ? "-1" : lexical_cast<std::string>(sparsity_[0]))
                       << " g[1]=" << (sparsity_.IsEmpty() ? "-1" : lexical_cast<std::string>(sparsity_[1]));
  return sparsity_;
}



double SlopeCondition::GetValue() const
{
  if(IsLocal())
    return values_[current_view_index_ - index_];
  else
    return value_;
}

void SlopeCondition::SetValue(double val)
{
  if(IsLocal())
  {
    values_[current_view_index_ - index_] = val;
    value_ = -1.0; // invalidated
  }
  else
    value_ = val; // only for Done() allowed!!
}


double SlopeCondition::GetMaxElementSlope(unsigned int element) const
{
  int size = 2 * dim_;
  int base = element * size;

  double max = 0.0;

  for(int i = 1; i < size; i++)
    max = std::max(max, abs(values_[base + i]));

  return max;
}

std::string SlopeCondition::ToString() const
{
  std::stringstream ss;

  ss << Condition::ToString();

  if(IsLocal())
    ss << "_ve=" << GetCurrentVirtualElement()
       << "_vn=" << GetCurrentVirtualNeighbor()
       << "_vs=" << GetCurrentVirtualSign();

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
  
  // make a warning for penalized material when using penalized volumes
  // TODO this check works not for penalizedVolume as objective!
  StdVector<Condition*> list = GetList(Condition::PENALIZED_VOLUME, DesignElement::DEFAULT, false);
  if(!list.IsEmpty())
  {
    for(unsigned int i = 0; i < space->transfer.GetSize(); i++)
      if(space->transfer[i].IsPenalized())
      {
        PtrParamNode in = info_->Get("constraint")->Get("type")->Get(Function::type.ToString(list[0]->GetType()));
        in->Get(ParamNode::WARNING)->SetValue("transfer function '" + space->transfer[i].ToString() + " seems also to penalize");
      }
  }

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

  // check if we have a slope constraint -
  list = GetList(Condition::SLOPE, DesignElement::DEFAULT, false);

  if(!list.IsEmpty())
  {

    InitSlopeConstraint(list[0], space, structure);
    view->Refresh(); // inform about the news
  }

  // check if we have a checkerboard - up to now a simple meassure
  list = GetList(Condition::CHECKERBOARD, DesignElement::DEFAULT, false); // must not be active

  if(!list.IsEmpty())
  {
    assert(list[0]->IsObservation());
    VicinityElement::Init(space, structure);
    assert(space->IsRegular()); // VicinityElements work only on a regular grid
  }
}

void ConditionContainer::ToInfo(PtrParamNode in)
{
  for(unsigned int i = 0; i < all.GetSize(); i++)
    all[i]->ToInfo(in->Get("constraint", ParamNode::APPEND));

  // save for later use in InitSlopeConstraint()
  info_ = in;
}

void ConditionContainer::InitSlopeConstraint(Condition* g_in, DesignSpace* space, DesignStructure* structure)
{
  SlopeCondition* g = dynamic_cast<SlopeCondition*>(g_in);

  PtrParamNode in = info_->Get("constraint")->Get("type")->Get(g->type.ToString(g->GetType()));

  // the design elements require the vicinity element to be set which holds the direct
  // neighbors. Is save to call several times
  VicinityElement::Init(space, structure);

  assert(space->IsRegular()); // VicinityElements work only on a regular grid

  g->PostInit(space);

  in->Get("bound_value")->SetValue(g->GetBoundValue());
  in->Get("active_size")->SetValue(g->GetData().GetSize());

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

Condition* ConditionContainer::Get(Condition::Type type, DesignElement::Type design)
{
  // be save and check for uniqueness!
  StdVector<Condition*> list = GetList(type, design, true);

  if(list.GetSize() == 0)
    throw Exception("no active constraint '" + Condition::type.ToString(type) + "' found");

  if(list.GetSize() > 1)
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
  // set slope_index
  StdVector<Condition*> list = container_->GetList(Condition::SLOPE, DesignElement::NO_TYPE, false); // also observe

  SlopeCondition* slope = list.IsEmpty() ? NULL :dynamic_cast<SlopeCondition*>(list[0]);

  slope_index_ = slope == NULL ? -1 : slope->GetIndex();

  virtual_active_size_ = container_->active.GetSize();
  // replace the global slope by many local slopes -> if it is initialized!
  if(slope != NULL && slope->IsActive() && slope->GetConstraintSize() > 0)
    virtual_active_size_ += slope->GetConstraintSize() -1;

  virtual_total_size_ = container_->all.GetSize();
  if(slope != NULL && slope->GetConstraintSize() > 0)
    virtual_total_size_ += slope->GetConstraintSize() -1;
}

Condition* ConditionContainer::VirtualView::Get(int view_index)
{
  // simple case: no slope or before slope
  if(slope_index_ == -1 || view_index < slope_index_)
    return container_->all[view_index];

  SlopeCondition* slope = dynamic_cast<SlopeCondition*>(container_->all[slope_index_]);

  // after_slope index is the first constraint after slope
  int next = slope_index_ + slope->GetConstraintSize();
  // are we within
  if(view_index < next)
  {
    assert(slope->GetIndex() == slope_index_);
    // tell the slope what is actually requested
    slope->SetCurrentViewIndex(view_index);
    return slope;
  }

  // no - we want a constraint after the slope. Add 1 as the slop consists in the container of one entry
  int corrected = view_index - slope->GetConstraintSize() + 1;
  return container_->all[corrected];
}

void ConditionContainer::VirtualView::Done()
{
  if(slope_index_ == -1)
    return; // nothing to do

  int dim = domain->GetGrid()->GetDim();
  SlopeCondition* slope = dynamic_cast<SlopeCondition*>(container_->all[slope_index_]);

  // set global result
  slope->SetValue(slope->GetData().NormMax());

  // check for special result.
  int idx = container_->space_->GetSpecialResultIndex(DesignElement::DEFAULT, DesignElement::MAX_SLOPE);
  if(idx >= 0)
  {
    StdVector<DesignElement>& data = container_->space_->data;

    // elements with no full neighborhood don't exist as constraint, therefore loop twice
    for(unsigned int e = 0; e < data.GetSize(); e++)
      data[e].specialResult[idx] = 0.0; // initialize

    // constraints -> note we have full neighborhood!
    for(int i = 0, ni = slope->GetData().GetSize(); i < ni; i += dim * 2)
    {
      slope->SetCurrentViewIndex(i + slope->GetIndex());
      int elem_idx = slope->GetCurrentVirtualElement();
      data[elem_idx].specialResult[idx] = slope->GetMaxElementSlope(i / (dim * 2));
    }
  }

  // reset local, set to global
  slope->SetCurrentViewIndex(-1);

}
