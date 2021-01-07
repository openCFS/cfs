#include "DataInOut/Logging/LogConfigurator.hh"
#include "Optimization/Design/FeaturedDesign.hh"
#include "Optimization/Excitation.hh"

namespace CoupledField {

DEFINE_LOG(FD, "featuredDesign")

Enum<FeaturedDesign::IntStrategy> FeaturedDesign::intStrategy_;
unsigned int FeaturedDesign::dim_ = 99;

FeaturedDesign::FeaturedDesign(StdVector<RegionIdType>& regionIds, PtrParamNode pn, ErsatzMaterial::Method method)
: AuxDesign(regionIds, pn, method)
{
  SetEnums();

  this->dim_ = domain->GetGrid()->GetDim();
  this->export_fe_design_ = false; // we use the original design but don't communicate it via ReadDesignFromExtern(), ...
  this->tailing_aux_design_ = true; // we want our own design (e.g. SBD: param_, SMD: shape_param_ or better their opt_* versions) to take the role of DesignSpace::data

  this->mapping_timer_  = info_->Get("features/mapping/timer")->AsTimer();
  this->mapping_timer_->SetLabel("features_map");
  this->mapping_timer_->SetSub(); // already in eval_*
  this->gradient_timer_ = info_->Get("features/gradient/timer")->AsTimer();
  this->gradient_timer_->SetLabel("features_grad");
  this->gradient_timer_->SetSub(); // already in eval_*
}

void FeaturedDesign::SetEnums()
{
  intStrategy_.SetName("FeaturedDesign::IntStrategy");
  intStrategy_.Add(CONSTANT_FULL, "constant_full");
  intStrategy_.Add(FULL_OR_NOTHING, "full_or_nothing");

  combine.SetName("FeaturedDesign::Combine");
  combine.Add(MAX, "max");
  combine.Add(TANH_SUM, "tanh_sum");

  boundary.SetName("FeaturedDesign::Boundary");
  boundary.Add(TANH, "tanh");
  boundary.Add(LINEAR, "linear");
  boundary.Add(POLY, "poly");
}


void FeaturedDesign::PostInit(int objectives, int constraints)
{
  AuxDesign::PostInit(objectives, constraints);

  if(domain->GetOptimization() != NULL)
    opt_ = domain->GetOptimization();

  if(domain->GetOptimization() != NULL)
    CheckPlausibility();

  assert(objectives > 0);

  // holds also for opt_shape_param_
  for(unsigned int i = 0; i < shape_param_.GetSize(); ++i) {
    shape_param_[i]->PostInit(objectives, constraints);
  }
};

int FeaturedDesign::ReadDesignFromExtern(const double* space_in)
{
  assert(!std::isnan(scaling_));
  int old_design = design_id;

  // write aux design variables (slack and alpha if any) last
  assert(export_fe_design_ == false); // we do shape map
  assert(DesignSpace::GetNumberOfVariables() > 0); // we need this variables but they are hidden!

  bool new_design = false;

  for(unsigned int i = 0; i < opt_shape_param_.GetSize(); ++i) {
    double v = space_in[i] * scaling_;
    assert(!std::isnan(v));
    if(!new_design && v != opt_shape_param_[i]->GetPlainDesignValue()) {
      new_design = true;
    }

    opt_shape_param_[i]->SetDesign(v);

    LOG_DBG3(FD) << "RDFE: i=" << i << ", " << opt_shape_param_[i]->ToString() << " -> " << v;
  }

  // append aux design, might also change design_id
  AuxDesign::ReadDesignFromExtern(space_in);

  if(new_design && design_id <= old_design)
    design_id++;

  return design_id;
}

bool FeaturedDesign::CompareDesign(const double* space_in)
{
  for(unsigned int i=0; i < opt_shape_param_.GetSize(); i++)
  {
    double v = space_in[i] * scaling_;
    if(v != opt_shape_param_[i]->GetPlainDesignValue())
      return false;
  }

  // no change here, let aux_design decide
  return AuxDesign::CompareDesign(space_in);
}

int FeaturedDesign::WriteDesignToExtern(double* space_out, bool scaling) const
{
  double rscaling = scaling ? 1.0 / scaling_ : 1.0;

  for(unsigned int i=0; i < opt_shape_param_.GetSize(); ++i) {
    space_out[i] = opt_shape_param_[i]->GetPlainDesignValue() * rscaling;
    LOG_DBG3(FD) << "WDTE: out[" << i << "]=" << space_out[i];
  }

  AuxDesign::WriteDesignToExtern(space_out, scaling);

  LOG_DBG(FD) << "WDTE: di -> " << design_id;
  return design_id;
}

void FeaturedDesign::WriteGradientToExtern(StdVector<double>& out, DesignElement::ValueSpecifier vs, DesignElement::Access access, Function* f, bool scaling)
{
  LOG_DBG2(FD) << "WGTE: f=" << f->ToString() << " auxd=" << aux_design_.GetSize() << " d=" << shape_param_.GetSize()
      << " out=" << out.GetSize() << " outwindowstart=" << out.window.GetStart() << " outwindowsz=" << out.window.GetSize();
  assert(out.window.GetStart() + out.window.GetSize() <= out.GetSize());

  // MapFeatureGradient would be good to perform it for all functions concurrently, however this is not possible as it is not the case that first all
  // simp function gradients are called and then all exported. This would need rewriting some stuff in cfs!
  if(f->IsObjective() || !Function::IsLocal(f->GetType())) // don't map local functions
    MapFeatureGradient(f); // see comment above for what is necessary to cache the stuff

  assert(f != NULL);
  assert(opt_->objectives.data.GetSize() == 1); // implement multi objective and be careful!
  assert(!(vs == DesignElement::COST_GRADIENT && !f->IsObjective()));
  assert(vs == DesignElement::COST_GRADIENT || vs == DesignElement::CONSTRAINT_GRADIENT);
  assert(!opt_->GetMultipleExcitation()->DoMetaExcitation(f->ctxt)); // robustness and transformation don't make sense for feature map. In DesignSpace we do f->GetExcitation()->Apply()
  assert(design.GetSize() == 1); // only pseudo density, nothing else implemented yet
  assert(regions.GetSize() == 1); // needs to be as design shall be 1

  assert(out.window.Initialized());
  unsigned int base = out.window.GetStart();
  // out contains the Jacobian for possibly many functions. Snopt combines cost function (first) and then the constraints derivatives
  // Here we assume that out has a window set where to write to for the given function.
  if(f->HasDenseJacobian())
  {
    unsigned int end_opt = opt_shape_param_.GetSize();
    LOG_DBG(FD) << "WGTE: end_opt=" << end_opt << " ad=" << aux_design_.GetSize() << " w=" << out.window.GetSize();
    assert(end_opt + aux_design_.GetSize() == out.window.GetSize());
    for(unsigned int s = 0; s < end_opt; ++s)
    {
      assert(out.InWindow(base + s));

      double opt = opt_shape_param_[s]->GetPlainGradient(f);

      LOG_DBG3(FD) << "WGTE de=" << opt_shape_param_[s]->ToString();
      assert(!std::isnan(opt));

      out[base + s] = opt * scaling;

      LOG_DBG3(FD) << "WGTE f=" << f->ToString() << " ws=" << out.window.GetStart() << " s=" << s << " opt=" << opt << " -> " << out[base + s];
    }
    // add slack stuff. No need to cheat window size
    AuxDesign::WriteGradientToExtern(out, vs, access, f, scaling);
  }
  else
  {
    // uses opt_index_!
    StdVector<unsigned int>& sparsity = f->GetSparsityPattern();
    LOG_DBG2(FD) << "WGTE f=" << f->ToString() << " sparsity=" << sparsity.ToString();
    assert(out.window.GetSize() == sparsity.GetSize());
    for(unsigned int i = 0; i < sparsity.GetSize(); i++)
    {
      unsigned int s = sparsity[i];
      assert(s < opt_shape_param_.GetSize());
      LOG_DBG3(FD) << "WGTE i=" << i << " s=" << s << " base=" << base << " b+s=" << (base+s);

      assert(out.InWindow(base + i));
      double scale = scaling ? scaling_ : 1.0;
      assert(vs == BaseDesignElement::CONSTRAINT_GRADIENT);

      double opt = opt_shape_param_[s]->GetPlainGradient(f);

      out[base + i] = opt * scale;
    }
  }
}

void FeaturedDesign::Reset(DesignElement::ValueSpecifier vs, DesignElement::Type design)
{
  assert(design == BaseDesignElement::DEFAULT || design == BaseDesignElement::CP);
  for(unsigned int i=0; i < opt_shape_param_.GetSize(); i++)
  {
    opt_shape_param_[i]->Reset(vs);
//    assert(!opt_shape_param_[i].costGradient.IsEmpty());
  }

  AuxDesign::Reset(vs, design);
}

void FeaturedDesign::WriteBoundsToExtern(double* x_l, double* x_u) const
{
  for(unsigned int i=0; i < opt_shape_param_.GetSize(); i++)
  {
    x_l[i] = opt_shape_param_[i]->GetLowerBound() / scaling_;
    x_u[i] = opt_shape_param_[i]->GetUpperBound() / scaling_;
    LOG_DBG3(FD) << "WBTE: l[" << i << "]=" << x_l[i] << " u[" << i << "]=" << x_u[i];
  }

  AuxDesign::WriteBoundsToExtern(x_l, x_u);
}

void FeaturedDesign::SetupMeshStructure()
{
  Grid* grid = domain->GetGrid();
  assert(GetRegionIds().GetSize() == 1);

  StdVector<unsigned int> t = grid->GetBoundaries(GetRegionIds().First());
  assert(t.GetSize() == 3);

  n_.Fill(t.GetPointer(), t.GetSize());
  nx_ = n_[0];
  ny_ = n_[1];
  nz_ = n_[2];

  assert(!(dim_ == 2 && nz_ != 1));
}

void FeaturedDesign::SetupMapping()
{
  //set physical design, i.e. density field
  map_.Resize(data.GetSize());
  StdVector<Elem*> designElems;
  assert(GetRegionIds().GetSize() == 1);

  domain->GetGrid()->GetElems(designElems, GetRegionIds().First()); // FIXME assumes elements in designElems are ordered!
  assert(map_.GetSize() == designElems.GetSize());



  for(unsigned int i = 0, n = map_.GetSize(); i < n; i++)
  {
    Item& item = map_[i];
    // no idea what is going on there, even if I once coded it ?! - Fabian
    item.rho = &(data[Find(designElems[i]->elemNum)]); // is very fast and gives a layer for arbitrary element ordering in the mesh
    item.min_corner_value.Resize(1);
    item.max_corner_value.Resize(1);

    LOG_DBG3(FD) << "SM i=" << i << " elem=" << item.rho->elem->elemNum << " de_elem=" << designElems[i]->elemNum
                 << " coord=" << domain->GetGrid()->GetElemNodesCoord(item.rho->elem).ToString(1,false);
  }
}


int FeaturedDesign::FindDesign(DesignElement::Type dt, bool throw_exception) const
{
  // check for DENSITY, ...
  int idx = DesignSpace::FindDesign(dt, false);
  if(idx >= 0)
    return idx;

  assert(dt != DesignElement::SPAGHETTI); // don't know yet how to handle

  assert(dt == DesignElement::CP || dt == DesignElement::NODE
      || dt == DesignElement::NODE || dt == DesignElement::PROFILE
      || dt == DesignElement::SPLINE_BOX);

  if(dt == DesignElement::SPLINE_BOX)
    dt = DesignElement::CP; // return the node index

  for(unsigned int i = 0; i < shape_param_.GetSize(); i++)
    if(shape_param_[i]->GetType() == dt)
      return i;

  if(throw_exception)
    EXCEPTION("Design " << DesignElement::type.ToString(dt) << " no FEM based and no feature mapping design.");
  return -1;
}



void FeaturedDesign::CheckPlausibility()
{
  assert(opt_ != NULL);
  if(  (opt_->constraints.Has(Function::VOLUME) && opt_->constraints.Get(Function::VOLUME)->IsLinear())
      || (opt_->objectives.Has(Function::VOLUME) && opt_->objectives.Get(Function::VOLUME)->IsLinear()))
    throw Exception("Set 'volume' function to non-linear for the parametrization.");
}


double FeaturedDesign::Item::MaxDiffCornerValue() const
{
  assert(min_corner_value.GetSize() == max_corner_value.GetSize());
  assert(min_corner_value.GetSize() > 0);

  double diff = 0;

  for(unsigned int i = 0; i < min_corner_value.GetSize(); i++)
    diff = std::max(max_corner_value[i] - min_corner_value[i], diff);

  return diff;
}


int FeaturedDesign::Item::GetOrder(Vector<int>& order, const FeaturedDesign::NumInt& ni) const
{
  assert(order.GetSize() == min_corner_value.GetSize());
  assert(order.GetSize() == max_corner_value.GetSize());
  int max = 0;
  for(unsigned int s = 0; s < order.GetSize(); s++)
  {
    double min_val = min_corner_value[s];
    double max_val = max_corner_value[s];

    assert(min_val >= 0);
    assert(max_val >= min_val);
    assert(max_val <= 1);

    // prevent unused variable warning
    (void)(min_val);

    switch(ni.strategy_)
    {
    case FeaturedDesign::CONSTANT_FULL:
      order[s] = ni.max_order_;
      break;
    case FeaturedDesign::FULL_OR_NOTHING:
      order[s] = max_val < ni.sensitivity_ ? 0 : ni.max_order_;
      break;
    case FeaturedDesign::TAILORED:
      throw Exception("Not implemented");
      break;
    }

    max = std::max(max, order[s]);
  }
  return max;
}


void FeaturedDesign::NumInt::Init(FeaturedDesign* sbd, PtrParamNode pn, PtrParamNode info)
{
  this->sensitivity_ = pn->Get("sensitivity")->As<double>();
  assert(sensitivity_ > 0);
  this->max_order_ = pn->Get("integration_order")->As<int>();

  if(max_order_ < 1)
    info->SetWarning("minimal value for 'max_order' is '1'");

  this->strategy_ = intStrategy_.Parse(pn->Get("integration_strategy")->As<string>());

  cells_ = domain->GetGrid()->GetNumElems();
}


void FeaturedDesign::NumInt::ToInfo(PtrParamNode info) const
{
  info->Get("max_order")->SetValue(max_order_);
  info->Get("sensitivity")->SetValue(sensitivity_);
  info->Get("integration")->SetValue(FeaturedDesign::intStrategy_.ToString(strategy_));
  assert(cells_ > 0);
  PtrParamNode cells = info->Get("cells");
  cells->Get("integrate_fraction")->SetValue(int_cells_cnt_ / (double) cells_);
  cells->Get("avg_order")->SetValue(int_cells_order_sum_ / (double) int_cells_cnt_);
  cells->Get("total_int")->SetValue(std::pow(int_cells_order_sum_, domain->GetGrid()->GetDim()));
}


} // end of namespace


