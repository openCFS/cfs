#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ProgramOptions.hh"
#include "Optimization/Design/FeaturedDesign.hh"
#include "Optimization/Excitation.hh"
#include "Utils/ToolsFull.hh"

namespace CoupledField {

DEFINE_LOG(fd, "featuredDesign")

Enum<FeaturedDesign::IntStrategy> FeaturedDesign::intStrategy_;
unsigned int FeaturedDesign::dim_ = 99;

using DE = DesignElement; // shortcut

FeaturedDesign::FeaturedDesign(StdVector<RegionIdType>& regionIds, PtrParamNode pn, ErsatzMaterial::Method method)
: AuxDesign(regionIds, pn, method)
{
  setup_timer_->Start();

  SetEnums();

  this->dim_ = domain->GetDim();
  this->export_fe_design_ = false; // we use the original design but don't communicate it via ReadDesignFromExtern(), ...
  this->tailing_aux_design_ = true; // we want our own design (e.g. SBD: param_, SMD: shape_param_ or better their opt_* versions) to take the role of DesignSpace::data

  this->mapping_timer_  = info_->Get("features/mapping/timer")->AsTimer();
  this->mapping_timer_->SetLabel("features_map");
  this->mapping_timer_->SetSub(); // already in eval_*
  this->gradient_timer_ = info_->Get("features/gradient/timer")->AsTimer();
  this->gradient_timer_->SetLabel("features_grad");
  this->gradient_timer_->SetSub(); // already in eval_*

  setup_timer_->Stop();
}

void FeaturedDesign::SetEnums()
{
  intStrategy_.SetName("FeaturedDesign::IntStrategy");
  intStrategy_.Add(CONSTANT_FULL, "constant_full");
  intStrategy_.Add(FULL_OR_NOTHING, "full_or_nothing");

  combine.SetName("FeaturedDesign::Combine");
  combine.Add(MAX, "max");
  combine.Add(TANH_SUM, "tanh_sum");
  combine.Add(KS, "KS");
  combine.Add(P_NORM, "p-norm");
  combine.Add(SOFTMAX, "softmax");

  boundary.SetName("FeaturedDesign::Boundary");
  boundary.Add(TANH, "tanh");
  boundary.Add(LINEAR, "linear");
  boundary.Add(POLY, "poly");
  boundary.Add(QUINTIC, "quintic");
  boundary.Add(BEZIER, "bezier");

  orientation.Add(ROUNDED, "rounded");
  orientation.Add(STRAIGHT, "straight");
}


void FeaturedDesign::PostInit(int objectives, int constraints)
{
  AuxDesign::PostInit(objectives, constraints);

  setup_timer_->Start();

  if(domain->GetOptimization() != NULL)
    opt = domain->GetOptimization();

  if(domain->GetOptimization() != NULL)
    CheckPlausibility();

  assert(objectives > 0);

  // holds also for opt_shape_param_
  for(unsigned int i = 0; i < shape_param_.GetSize(); ++i) {
    shape_param_[i]->PostInit(objectives, constraints);
  }

  setup_timer_->Stop();
};

/** only used by SpaghettiDesign and FeaturedDesign(), others have their own */
void FeaturedDesign::SetupDesign(PtrParamNode base)
{
  assert(features_.IsEmpty());
  SetupParsedFeatures(base);
  assert(!features_.IsEmpty());
  int total = 0;
  int opt   = 0;
  for(Feature* feature :  features_)
  {
    total += feature->GetTotalVariables();
    opt   += feature->GetOptVariables();
  }

  LOG_DBG(fd) << "SD: total=" << total << " opt=" << opt;

  shape_param_.Reserve(total);
  opt_shape_param_.Reserve(opt);
  opt_indices.Resize(opt);
}

void FeaturedDesign::AddVariable(FeatureVariable* var)
{
  // we do exact space "estimation"

  assert(shape_param_.HasSpace());
  var->SetIndex(shape_param_.GetSize());
  shape_param_.Push_back(var);

  if(var->IsVariable())
  {
    assert(opt_shape_param_.HasSpace());
    opt_indices[opt_shape_param_.GetSize()] = shape_param_.GetSize()-1;
    var->SetOptIndex(opt_shape_param_.GetSize());
    opt_shape_param_.Push_back(var);
  }
}

int FeaturedDesign::ReadDesignFromExtern(const double* space_in, bool setAndWriteCurrent)
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

    LOG_DBG3(fd) << "RDFE: i=" << i << ", " << opt_shape_param_[i]->ToString() << " -> " << v;
  }

  // append aux design, might also change design_id
  AuxDesign::ReadDesignFromExtern(space_in, setAndWriteCurrent);

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
    LOG_DBG3(fd) << "WDTE: out[" << i << "]=" << space_out[i];
  }

  AuxDesign::WriteDesignToExtern(space_out, scaling);

  LOG_DBG(fd) << "WDTE: di -> " << design_id;
  return design_id;
}

void FeaturedDesign::WriteGradientToExtern(StdVector<double>& out, DesignElement::ValueSpecifier vs, DesignElement::Access access, Function* f, bool scaling)
{
  LOG_DBG(fd) << "WGTE: f=" << f->ToString() << " auxd=" << aux_design_.GetSize() << " d=" << shape_param_.GetSize()
      << " out=" << out.GetSize() << " outwindowstart=" << out.window.GetStart() << " outwindowsz=" << out.window.GetSize();
  assert(out.window.GetStart() + out.window.GetSize() <= out.GetSize());

  // MapFeatureGradient would be good to perform it for all functions concurrently, however this is not possible as it is not the case that first all
  // simp function gradients are called and then all exported. This would need rewriting some stuff in cfs!
  if(f->IsObjective() || !Function::IsLocal(f->GetType())) // don't map local functions
    MapFeatureGradient(f); // see comment above for what is necessary to cache the stuff

  assert(f != NULL);
  assert(opt->objectives.data.GetSize() == 1); // implement multi objective and be careful!
  assert(!(vs == DesignElement::COST_GRADIENT && !f->IsObjective()));
  assert(vs == DesignElement::COST_GRADIENT || vs == DesignElement::CONSTRAINT_GRADIENT);
  assert(!opt->GetMultipleExcitation()->DoMetaExcitation(f->ctxt)); // robustness and transformation don't make sense for feature map. In DesignSpace we do f->GetExcitation()->Apply()
  assert(!(regions.GetSize() > 1 && GetMethod() != ErsatzMaterial::SPAGHETTI_PARAM_MAT)); // only param mat can have multiple designs

  assert(out.window.Initialized());
  unsigned int base = out.window.GetStart();
  // out contains the Jacobian for possibly many functions. Snopt combines cost function (first) and then the constraints derivatives
  // Here we assume that out has a window set where to write to for the given function.
  if(f->HasDenseJacobian())
  {
    unsigned int end_opt = opt_shape_param_.GetSize();
    LOG_DBG(fd) << "WGTE: end_opt=" << end_opt << " ad=" << aux_design_.GetSize() << " w=" << out.window.GetSize();
    assert(end_opt + aux_design_.GetSize() == out.window.GetSize());
    for(unsigned int s = 0; s < end_opt; ++s)
    {
      assert(out.InWindow(base + s));

      double opt = opt_shape_param_[s]->GetPlainGradient(f);

      LOG_DBG3(fd) << "WGTE de=" << opt_shape_param_[s]->ToString();
      assert(!std::isnan(opt));

      out[base + s] = opt * scaling;

      LOG_DBG3(fd) << "WGTE f=" << f->ToString() << " ws=" << out.window.GetStart() << " s=" << s << " opt=" << opt << " -> " << out[base + s];
    }
    // add slack stuff. No need to cheat window size
    AuxDesign::WriteGradientToExtern(out, vs, access, f, scaling);
  }
  else
  {
    // uses opt_index_!
    StdVector<unsigned int>& sparsity = f->GetSparsityPattern();
    LOG_DBG2(fd) << "WGTE f=" << f->ToString() << " sparsity=" << sparsity.ToString();
    assert(out.window.GetSize() == sparsity.GetSize());
    for(unsigned int i = 0; i < sparsity.GetSize(); i++)
    {
      unsigned int s = sparsity[i];
      assert(s < opt_shape_param_.GetSize());
      LOG_DBG3(fd) << "WGTE i=" << i << " s=" << s << " base=" << base << " b+s=" << (base+s);

      assert(out.InWindow(base + i));
      double scale = scaling ? scaling_ : 1.0;
      assert(vs == BaseDesignElement::CONSTRAINT_GRADIENT);

      double opt = opt_shape_param_[s]->GetPlainGradient(f);

      out[base + i] = opt * scale;
    }
  }
}

void FeaturedDesign::OpenGradPlot(PtrParamNode pn)
{
  if(pn->Get("gradplot")->As<bool>())
    gradplot_.open((progOpts->GetSimName() + ".grad.dat").c_str()); // the auto destructor does the job.
}

void FeaturedDesign::WriteGradientFile()
{
  if(!gradplot_.is_open())
    return; // obviously the option was not set

  gradplot_.precision(5);
  gradplot_.flags(std::ios::scientific);

  int iter = opt->GetCurrentIteration();

  int cnt = 0;
  if(iter == 0)
  {
    gradplot_ << "# " << ++cnt << ": iter \t";
    for(const ShapeParamElement* spe : opt_shape_param_)
      gradplot_ << ++cnt << ": " << spe->GetLabel() << " \t";

    for(const Function* f : opt->GetFunctions(false))
      if(!f->IsLocal())
        for(const ShapeParamElement* spe : opt_shape_param_)
          gradplot_ << ++cnt << ": d_" << f->ToString() << " / d_" << spe->GetLabel() << " \t";

    gradplot_ << std::endl;
  }
  gradplot_ << iter << " \t";
  for(const ShapeParamElement* spe : opt_shape_param_)
    gradplot_ << spe->GetPlainDesignValue() << " \t";

  for(const Function* f : opt->GetFunctions(false))
    if(!f->IsLocal())
      for(const ShapeParamElement* spe : opt_shape_param_)
        gradplot_ << spe->GetPlainGradient(f) << " \t";

  gradplot_ << std::endl;
  // trust on the destructor to close the file
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
    LOG_DBG3(fd) << "WBTE: l[" << i << "]=" << x_l[i] << " u[" << i << "]=" << x_u[i];
  }

  AuxDesign::WriteBoundsToExtern(x_l, x_u);
}

void FeaturedDesign::SetupMeshStructure()
{
  Grid* grid = domain->GetGrid();

  LOG_DBG(fd) << "SMS: regions=" << GetRegionIds().ToString();

  assert(GetRegionIds().GetSize() == 1);
  StdVector<unsigned int> t = grid->GetRegularDiscretization(GetRegionIds().First());
  assert(t.GetSize() == 3);

  // note that spaghetti.py assumes its domain to be rectangular and start from (0,0,0) with a dx stepping
  // that implies that the point coordinates are not necessarily the mesh coordinates.
  // also nx_*ny_*nz might be larger than the data size for density!
  // @see SetupMapping()
  n_.Fill(t.GetPointer(), t.GetSize());
  nx_ = n_[0];
  ny_ = n_[1];
  nz_ = n_[2];
  LOG_DBG(fd) << "SMS: n=" << n_.ToString() << " #=" << (nx_ * ny_ * nz_) << " #data=" << data.GetSize();

  // We need the spacing of an element
  assert(data[0].elem != nullptr);
  shared_ptr<ElemShapeMap> esm = grid->GetElemShapeMap(data[0].elem);
  double minEdge, maxEdge;
  esm->GetMaxMinEdgeLength(maxEdge,minEdge);
  LOG_DBG(fd) << "SMS: dx_=min=" << minEdge << " max=" << maxEdge;
  assert(close(minEdge, maxEdge));
  dx_ = minEdge;

  assert(!(dim_ == 2 && nz_ != 1));
}

void FeaturedDesign::SetupMapping()
{
  // n_ and nx_, ny_, nz_
  SetupMeshStructure();

  // set physical design, which is usually the density but for spaghetti also angles.
  map.Resize(data.GetSize());
  StdVector<Elem*> designElems;

  assert(GetRegionIds().GetSize() == 1);
  domain->GetGrid()->GetElems(designElems, GetRegionIds().First()); // FIXME assumes elements in designElems are ordered!

  LOG_DBG(fd) << "SM: #map=" << map.GetSize() << " #designElems=" << designElems.GetSize();

  for(unsigned int i = 0, n = map.GetSize(); i < n; i++)
  {
    Item& item = map[i];
    item.elemval = &data[i]; // this could be the location to add an arbitrary element ordering layer towards the mesh
    item.lexicographic_pos = i; // this means we have a 1:1 map=data to tmp=virtual rectangular mesh for spaghetti -> corrected a few lines below
    //item.elemval = &(data[Find(designElems[i]->elemNum)]); // is very fast and gives a layer for arbitrary element ordering in the mesh
    item.min_corner_value.Resize(1);
    item.max_corner_value.Resize(1);

    LOG_DBG3(fd) << "SM i=" << i << " elem=" << item.elemval->elem->elemNum << " de_elem=" << designElems[i]->elemNum
                 << " coord=" << domain->GetGrid()->GetElemNodesCoord(item.elemval->elem).ToString();
  }

  assert(!(n_.GetSize() == 0));
  // #designElems can be smaller nx_*ny_*ny_ -> see comment in SetupMeshStructure()
  assert(!(designElems.GetSize() > n_.Product()));
  assert(elements == designElems.GetSize());
  assert(design.GetSize() * designElems.GetSize() == map.GetSize()); // # slack is in aux design

  if(designElems.GetSize() < n_.Product())
  {
    // size of map and designElems differ for density + angles
    Matrix<double> mm;
    domain->GetGrid()->CalcBoundingBoxOfRegion(GetRegionIds().First(), mm);
    LOG_DBG(fd) << "SM: bb=" << mm.ToString();
    assert(mm.GetNumRows() == dim_);
    Vector<double> lower(dim_);
    for(unsigned int i = 0; i < dim_; i++)
      lower[i] = mm[i][0];
    LOG_DBG(fd) << "SM: lower=" << lower.ToString();
    Vector<double>  glob(dim_);
    LocPoint loc;
    unsigned int virt_size = n_.Product();
    for(unsigned int z = 0; z < nz_; z++) {
      for(unsigned int y = 0; y < ny_; y++) {
        for(unsigned int x = 0; x < nx_; x++) {
          glob[0] = lower[0] + x * dx_ + .5*dx_;
          glob[1] = lower[1] + y * dx_ + .5*dx_;
          if(dim_ == 3)
            glob[2] = lower[2] + z * dx_ + .5*dx_;
          // being sure whe are regular and ordered as expected, we could directly
          // address the elements. This is our saveguard. With proper meshes, the 
          // default implementation is cached and fast. For wrong meshes it is totally
          // slow. With CGAL or LIBFBI it should make sense to use GetElemsAtGlobalCoords()
          const Elem* elem = domain->GetGrid()->GetElemAtGlobalCoord(glob, loc);
          assert(elem != nullptr);
          int idx = Find(elem, false); // design index or -1
          int pos = LexicographicPos(x,y,z);
          LOG_DBG3(fd) << "(" << x << "," << y << "," << z << ") g=" << glob.ToString() << " l=" << loc.coord.ToString() << " e=" << elem->elemNum << " idx=" << idx << " pos=" << pos;
          assert(pos < (int) virt_size);
          if(idx > -1)
            for(unsigned int d = 0; d < design.GetSize(); d++)
              map[d * elements + idx].lexicographic_pos = d * virt_size + pos;
        }
      }
    }
  }


}


int FeaturedDesign::FindDesign(DesignElement::Type dt, bool throw_exception) const
{
  // check for DENSITY, ...
  int idx = DesignSpace::FindDesign(dt, false);
  if(idx >= 0)
    return idx;

//  assert(dt == DesignElement::CP || dt == DesignElement::SPLINE_BOX || dt == DesignElement::SPAGHETTI
//      || dt == DesignElement::NODE || dt == DesignElement::PROFILE || dt == DesignElement::NORMAL);

  if(dt == DesignElement::SPLINE_BOX)
    dt = DesignElement::CP; // return the node index
  if(dt == DesignElement::SPAGHETTI)
    dt = DesignElement::NODE; // extremely ugly!! but to make bending work easily

  for(unsigned int i = 0; i < shape_param_.GetSize(); i++)
    if(shape_param_[i]->GetType() == dt)
      return i;

  if(throw_exception)
    EXCEPTION("Design " << DesignElement::type.ToString(dt) << " no FEM based and no feature mapping design.");
  return -1;
}

void FeaturedDesign::SetupVirtualShapeElementMap(Function* f, StdVector<Function::Local::Identifier>& vem, Function::Local::Locality locality)
{
  assert(f != NULL);
  assert(f->IsLocal(f->GetType()));
  // we shall be called by Local::PostInit() therefore local shall exist
  assert(f->GetLocal() != NULL);
  assert(f->GetType() == Function::DISTANCE || f->GetType() == Function::BENDING);

  // we assume fixed only for profile if at all
  assert(locality == Function::Local::FUNCTION_SPECIFIC || locality == Function::Local::FUNCTION_SPECIFIC_TWO_SIGNS);

  assert(!features_.IsEmpty());
  vem.Reserve(features_.GetSize()); // assume nothing fixed

  assert(dim_ == 2);
  StdVector<BaseDesignElement*> nodes;

  bool two_signs = locality == Function::Local::FUNCTION_SPECIFIC_TWO_SIGNS;

  int sign_1 = two_signs ? 1 : Function::Local::Identifier::NO_SIGN;
  int sign_2 = -1;

  for(Feature* s : features_)
  {
    nodes.Clear();
    StdVector<FeatureVariable>& P = s->points.First();
    StdVector<FeatureVariable>& Q = s->points.Last();

    // fixed nodes are constants: LocalCondition drops them from the sparsity pattern and Hessian
    // (EffectiveOptIndex() == -1). Without any free node the length is constant -> no constraint.
    if(f->GetType() == Function::DISTANCE)
    {
      if(FeatureVariable::CountRealVariables(P) + FeatureVariable::CountRealVariables(Q) == 0)
        continue; // next noodle
    }
    else // Bending. Note IsFixed() is true for any fixed component
      if(FeatureVariable::IsFixed(P) && FeatureVariable::IsFixed(Q) && !s->IsExtended())
        continue; // won't add empty constraint if all points are fixed -> next noodle

    // px is element, then py, then qx then qy
    nodes.Push_back(&P[1]);
    nodes.Push_back(&Q[0]);
    nodes.Push_back(&Q[1]);

    if(f->GetType() == Function::DISTANCE)
    {
      if(f->GetDesignType() != BaseDesignElement::NODE)
        throw Exception("Configuration error: distance needs to be defined for node.");

      vem.Push_back(Function::Local::Identifier(&P[0], nodes));
      LOG_DBG(fd) << "SVSEM: f=" << f->ToString() << " s=" << s->idx << " id=" << vem.Last().ToString();
    }
    else // Bending
    {
      assert(f->GetType() == Function::BENDING);
      SetupVirtualShapeElementMapBending(s, vem, nodes, two_signs, sign_1, sign_2);
    } // end bending
  } // end noodle
}


void FeaturedDesign::ToInfo(ErsatzMaterial* em)
{
  AuxDesign::ToInfo(em);

  PtrParamNode in = info_->Get("features");

  in->Get("combine/type")->SetValue(combine.ToString(combine_));
  in->Get("boundary/type")->SetValue(boundary.ToString(boundary_));
  
  PtrParamNode msh = in->Get("mesh");
  msh->Get("n")->SetValue(n_.ToString());

  PtrParamNode base = info_->Get("designVariables");
  for(Feature* f : features_)
    f->ToInfo(base->Get(f->GetName(), ParamNode::APPEND));
}

void FeaturedDesign::CheckPlausibility()
{
  assert(opt != NULL);
  if(  (opt->constraints.Has(Function::VOLUME) && opt->constraints.Get(Function::VOLUME)->IsLinear())
      || (opt->objectives.Has(Function::VOLUME) && opt->objectives.Get(Function::VOLUME)->IsLinear()))
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
    assert(max_val <= 1.01);

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
  cells->Get("total_int")->SetValue(std::pow(int_cells_order_sum_, domain->GetDim()));
}

void FeaturedDesign::Feature::Parse(PtrParamNode pn, int idx)
{
  this->idx = idx;
  opt_variables_ = 0;

  points.Resize(2); // start and end point

  // add tip
  StdVector<FeatureVariable>& start = points.First();
  start.Resize(dim_);
  start[0].Parse(pn->GetByVal("node", "dof", "x", "tip", "start"),idx, DE::FEATURE_MAPPING_PX);
  start[1].Parse(pn->GetByVal("node", "dof", "y", "tip", "start"),idx, DE::FEATURE_MAPPING_PY);
  if(dim_ == 3)
    start[2].Parse(pn->GetByVal("node", "dof", "z", "tip", "start"),idx, DE::FEATURE_MAPPING_PZ);
  opt_variables_ += FeatureVariable::CountRealVariables(start);

  StdVector<FeatureVariable>& end = points.Last();
  end.Resize(dim_);
  end[0].Parse(pn->GetByVal("node", "dof", "x", "tip", "end"),idx, DE::FEATURE_MAPPING_QX);
  end[1].Parse(pn->GetByVal("node", "dof", "y", "tip", "end"),idx, DE::FEATURE_MAPPING_QY);
  if(dim_ == 3)
    end[2].Parse(pn->GetByVal("node", "dof", "z", "tip", "end"),idx, DE::FEATURE_MAPPING_QZ);
  opt_variables_ += FeatureVariable::CountRealVariables(end);

  p.Parse(pn->Get("profile"),idx, DE::FEATURE_MAPPING_P); // exactly one
  opt_variables_ += p.IsVariable() ? 1 : 0;

  LOG_DBG(fd) << "F:P ov=" << opt_variables_ << " " << ToString();
}

int FeaturedDesign::Feature::GetTotalVariables() const
{
  int vars = 1; // profile p
  vars += dim_ * points.GetSize();
  return vars;
}


string FeaturedDesign::Feature::ToString() const
{
  std::stringstream ss;

  ss << "idx=" << idx;
  if(points.GetSize() > 1)
  {
    ss << " start=" << FeatureVariable::ToString(points.First());
    ss << " end=" << FeatureVariable::ToString(points.Last());
    for(unsigned int i = 1; i < points.GetSize()-1; i++) // 3d spaghetti
      ss << " num_" << i << "=" << FeatureVariable::ToString(points[i]);
  }
  ss << " p=" << (p.fixed ? "*" : "") << p.GetPlainDesignValue();
  return ss.str();
}

void FeaturedDesign::Feature::ToInfo(PtrParamNode in)
{
  in->Get("id")->SetValue(idx);

  // one could combine P and Q as coordinates.
  auto& P = points.First();
  auto& Q = points.Last();

  for(unsigned int d = 0; d < dim_; d++)
    P[d].ToInfo(in->Get("p" + Dof(d), ParamNode::APPEND));

  for(unsigned int d = 0; d < dim_; d++)
    Q[d].ToInfo(in->Get("p" + Dof(d), ParamNode::APPEND));

  p.ToInfo(in->Get("p", ParamNode::APPEND));
}


void FeaturedDesign::Feature::ToInfo(PtrParamNode in, FeatureVariable::Type type, const std::string& label, int dim) const
{
  std::string fixed;
  std::string lower;
  std::string upper;
  const FeatureVariable* v0 = nullptr;
  int size = -1;

  assert(type == FeatureVariable::NODE);
  assert(dim >= 0 && dim <= 2);
  if(points.GetSize() <= 2)
    return;

  v0 = &(points[1][dim]);
  assert(v0->tip == FeatureVariable::INNER);
  size = points.GetSize() - 2; // minus P and Q

  for(unsigned int i = 1; i < points.GetSize()-1; i++)
    FeatureVariable::CompareToInfoHelper(v0, &points[i][dim], fixed, lower, upper);
 
  v0->ToInfo(in->Get(label, ParamNode::APPEND));
  in->Get(label + "/size")->SetValue(size);
  if(v0->fixed)
    in->Get(label + "/fixed")->SetValue(fixed);
  else {
    in->Get(label + "/lower")->SetValue(lower);
    in->Get(label + "/upper")->SetValue(upper);
  }
}

void FeaturedDesign::Feature::GetAllVariables(StdVector<FeatureVariable*>& out) const
{
  unsigned int total = points.GetSize() * domain->GetDim() + 1; // start + end + p
  out.Clear();
  out.Reserve(total); 

  for(const StdVector<FeatureVariable>& vec : points) 
    for(const FeatureVariable& var : vec)
      out.Push_back(const_cast<FeatureVariable*>(&var));
  out.Push_back(const_cast<FeatureVariable*>(&p));

  assert(out.GetSize() == total); // we predicted properly
}

StdVector<FeatureVariable*> FeaturedDesign::Feature::GetAllVariables() const
{
  StdVector<FeatureVariable*> vec;
  GetAllVariables(vec);
  return vec;
}




} // end of namespace


