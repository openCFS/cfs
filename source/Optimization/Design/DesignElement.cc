#include "General/exception.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Design/DesignStructure.hh"
#include "Optimization/Condition.hh"
#include "Optimization/Objective.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"
#include "Elements/basefe.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "boost/lexical_cast.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"
#include "Optimization/LevelSet.hh"
#include "Optimization/Optimizer/ShapeOptimizer.hh"
#include "Utils/Timer.hh"

using namespace std;
using namespace CoupledField;

using boost::posix_time::ptime;
using boost::posix_time::second_clock;
using boost::posix_time::microsec_clock;

DECLARE_LOG(del)
DEFINE_LOG(del, "designElement")

// the static enum
Enum<Filter::Type>                  Filter::type;
Enum<Filter::Density>               Filter::density;
Enum<Filter::Sensitivity>           Filter::sensitivity;
Enum<DesignElement::Type>           DesignElement::type;
Enum<DesignElement::ValueSpecifier> DesignElement::valueSpecifier;
Enum<DesignElement::Access>         DesignElement::access;
Enum<DesignElement::Detail>         DesignElement::detail;

// is a static attribute
DesignSpace* DesignElement::space_(NULL);

BaseDesignElement::BaseDesignElement() {
  design          = 0.0;
  upper_          = 0.0;
  lower_          = 0.0;
}


void BaseDesignElement::PostInit(int objectives, int constraints)
{
  // resize and init with 0.0 so constraint, which only act on one design variable, do not have to set all others explicitly to zero
  costGradient.Resize(objectives, 0.0);
  constraintGradient.Resize(constraints, 0.0);
}


/** Get the gradient values for either objective or constraint */
double BaseDesignElement::GetPlainGradient(const Objective* c, const Condition* g) const
{
  assert(c == NULL || g == NULL);

  if(g != NULL) return constraintGradient[g->GetIndex()];
  if(c != NULL) return costGradient[c->GetIndex()];

  return SumObjectiveGradient();
}

/** Get the gradient values for either objective or constraint */
double BaseDesignElement::GetPlainGradient(const Function* f) const
{
  assert(!f->IsObjective() || (f->IsObjective() && dynamic_cast<const Objective*>(f) != NULL));
  assert( f->IsObjective() || (!f->IsObjective() && dynamic_cast<const Condition*>(f) != NULL));

  return GetPlainGradient(f->IsObjective() ? static_cast<const Objective*>(f) : NULL,
                           f->IsObjective() ? NULL : static_cast<const Condition*>(f));
}


/** Sum app the old value (get and set together) */
void BaseDesignElement::AddGradient(const Objective* f, const Condition* g, double value)
{
  assert(f == NULL || g == NULL);
  LOG_DBG3(del) << "AddGradient: f=" << (f == NULL ? "null" : f->type.ToString(f->GetType()))
                << " g=" << (g == NULL ? "null" : g->ToString()) << " val=" << value
                << " penalty=" << (f != NULL ? boost::lexical_cast<std::string>(f->GetPenalty()) : "-")
                << " old= " <<  (f != NULL ? costGradient[f->GetIndex()] : constraintGradient[g->GetIndex()])
                << " add=" << (f != NULL ? value * f->GetPenalty() : value)
                << " -> " << (f != NULL ? costGradient[f->GetIndex()] + value * f->GetPenalty() : constraintGradient[g->GetIndex()] + value);
  if(f != NULL) costGradient[f->GetIndex()] += value * f->GetPenalty();
           else constraintGradient[g->GetIndex()] += value;
}

void BaseDesignElement::AddGradient(const Function* f, double value)
{
  assert(!f->IsObjective() || (f->IsObjective() && dynamic_cast<const Objective*>(f) != NULL));
  assert( f->IsObjective() || (!f->IsObjective() && dynamic_cast<const Condition*>(f) != NULL));

  AddGradient(f->IsObjective() ? static_cast<const Objective*>(f) : NULL,
              f->IsObjective() ? NULL : static_cast<const Condition*>(f), value);
}

void BaseDesignElement::Reset(ValueSpecifier vs, Function*  f)
{
  switch(vs)
  {
  case COST_GRADIENT:
    for(unsigned int i = 0; i < costGradient.GetSize(); i++)
      costGradient[i] = 0.0;
    break;
  case CONSTRAINT_GRADIENT:
    if(f != NULL)
      constraintGradient[f->GetIndex()] = 0.0;
    else
      for(unsigned int i = 0; i < constraintGradient.GetSize(); i++)
        constraintGradient[i] = 0.0;
    break;
  default:
    EXCEPTION("invalid value specifier " << vs);
  }
}

double BaseDesignElement::SumObjectiveGradient() const
{
  double result = 0.0;
  for(unsigned int i = 0; i < costGradient.GetSize(); i++)
    result += costGradient[i];

  return result;
}

/** The default constructor for StdVector and ghost elements*/
DesignElement::DesignElement() : BaseDesignElement()
{
  Init();
}

DesignElement::DesignElement(PtrParamNode pn, Elem* elem) : BaseDesignElement()
{
  Init();
  this->elem = elem;
  this->specialResult.Resize(9, 0.0);

  // it is a little slow to perform this code for every DesignElement but the
  // implementations are rater fast and it should be not measurable in the end
  type_ = type.Parse(pn->Get("name")->As<std::string>());

  upper_ = 1.0;
  // eventually overwrite
  pn->GetValue("upper", upper_, ParamNode::INSERT);

  lower_ = type_ == POLARIZATION ? -1.0 : 0.001;
  pn->GetValue("lower", lower_, ParamNode::INSERT);
}


DesignElement::~DesignElement()
{
  delete simp; simp = NULL;
  delete vicinity; vicinity = NULL;
  delete location_; location_ = NULL;
}

void DesignElement::Init()
{
  simp            = NULL;
  vicinity       = NULL;
  lse_            = NULL;
  tge             = NULL;
  location_       = NULL;
  elem            = NULL;
  type_           = NO_TYPE;
}

Point* DesignElement::GetLocation()
{
  if(location_ != NULL) return location_;
  
  location_ = new Point();
  // calc location
  
  // check for ghost elements
  if(elem == NULL) EXCEPTION("location_ not set but elem is NULL");
    
  Matrix<double>  coords; // we ignore the n times constructs

  StdVector<unsigned int>& connect = elem->connect;
  // do not use updated coordinates up to now!!
  domain->GetGrid()->GetElemNodesCoord(coords, connect, false );
  // a barycenter is simply the average of all coordinates -> location_ gets the Point out
  BaseFE::CalcBarycenter(coords, *location_);

  LOG_DBG3(del) << "DesignElement::GetLocation() find " << location_->ToString() << " for " << ToString();
  return location_;
}

void DesignElement::GetValue(ResultDescription& rd, StdVector<double>& out, unsigned int dofs) const
{
  // check for special result
  if(    rd.value == OBJECTIVE
      || (rd.value == COST_GRADIENT && rd.detail != NONE)
      || rd.value == CONSTRAINT_GRADIENT
      || rd.value == MAX_SLOPE
      || rd.value == MAX_OSCILLATION
      || rd.value == MAX_MOLE
      || rd.value == MAX_JUMP
      || rd.value == PENALIZED_STRESS)
  {
    if(dofs != 1) throw Exception("special results is only defined for scalar values");
    // note, that on EACH_FORWARD/ADJOINT we need excitation based results
    switch(rd.solutionType)
    {
    case OPT_RESULT_1:
      out[0] = specialResult[0];
      break;
    case OPT_RESULT_2:
      out[0] = specialResult[1];
      break;
    case OPT_RESULT_3:
      out[0] = specialResult[2];
      break;
    case OPT_RESULT_4:
      out[0] = specialResult[3];
      break;
    case OPT_RESULT_5:
      out[0] = specialResult[4];
      break;
    case OPT_RESULT_6:
      out[0] = specialResult[5];
      break;
    case OPT_RESULT_7:
      out[0] = specialResult[6];
      break;
    case OPT_RESULT_8:
      out[0] = specialResult[7];
      break;
    case OPT_RESULT_9:
      out[0] = specialResult[8];
      break;

    default: throw Exception("solution type not handled");
    }
  }
  else
  {
    if(dofs == 1)
    {
      if(rd.solutionType == PHYSICAL_PSEUDO_DENSITY)
        out[0] = GetPhysicalDesign();
      else
        out[0] = GetValue(rd.value, rd.access);
    }
    else
    {
      throw Exception("we only have one dof, only scalar values");
    }
  }
  return;
}

double DesignElement::GetValue(ValueSpecifier vs, Access access, Condition* g) const
{
  // pessimistic assumption :)
  bool sens_filter = false;
  bool design_filter = false;
  bool design_filter_grad = false;

  if(access == SMART && simp != NULL && simp->filter.type_ != Filter::NO_FILTERING)
  {
    Filter& f = simp->filter;
    if(f.type_ == Filter::DENSITY)
    {
      if(vs == DesignElement::DESIGN)
        design_filter = true;
      if(vs == DesignElement::COST_GRADIENT) // TODO also check for separate costs! not only sum!
        design_filter_grad = true;
      if(vs == DesignElement::CONSTRAINT_GRADIENT)
        design_filter_grad = g->ForDensityFiltering();
    }
    else
    {
      if(vs == DesignElement::COST_GRADIENT)// TODO also check for separate costs! not only sum!
        sens_filter = true;
      if(vs == DesignElement::CONSTRAINT_GRADIENT)
        sens_filter = g->ForSensitivityFiltering();
    }
  }

  // we are silent if one wants filtering and it is not possible
  double val = 0.0;
  if(sens_filter)
    val = simp->GetSensitivityFilteredValue(vs, g);
  if(design_filter)
    val = simp->GetDensityFilteredValue(vs, simp->filter.density_);
  if(design_filter_grad)
    val = simp->GetDensityFilteredGradient(vs, g);
  if(!sens_filter && !design_filter && !design_filter_grad)
    val = GetPlainValue(vs, g);

  LOG_DBG3(del) << "GV: " << elem->elemNum << " (" << valueSpecifier.ToString(vs) << ", "
                << (access == PLAIN ? "plain" : "smart") << ", " << (g == NULL ? "null" : g->ToString())
                << ") sf=" << sens_filter << " df=" << design_filter << " dfg=" << design_filter_grad << " -> " << val;
  return val;
}


double DesignElement::GetPlainValue(ValueSpecifier sp, Condition* g) const
{
  // validate first:
  switch(sp)
  {
  case DESIGN:                return design;

  case COST_GRADIENT:         return SumObjectiveGradient();
  case WEIGHT:
    if(simp == NULL) throw Exception("'" + valueSpecifier.ToString(sp) + "' is specific to SIMP");
    return simp->weight;

  case NUM_NEIGHBOURS:       
    if(simp == NULL) throw Exception("'" + valueSpecifier.ToString(sp) + "' is specific to SIMP");
    return simp->neighborhood.GetSize();

  case CONSTRAINT_GRADIENT:
    assert(g != NULL);
    return constraintGradient[g->GetIndex()];

  case MAX_SLOPE:
  case MAX_MOLE:
  case MAX_OSCILLATION:
  case MAX_JUMP:
  case PENALIZED_STRESS:
    assert(false); // should be covered before by special result index

  case TOPGRAD_VALUE:
    if(tge == 0) throw Exception("'" + valueSpecifier.ToString(sp) + "' is specific to TopGrad");
    return tge->value;

  case SHAPEGRAD_VALUE:
    if(lse_ == 0) throw Exception("'" + valueSpecifier.ToString(sp) + "' is specific to Levelset");
    return lse_->shapeGradValue;
  default: throw Exception(valueSpecifier.ToString(sp) + " is no scalar value");
  }
}


double DesignElement::GetDesign(Access access) const
{
  return GetValue(DESIGN, access);
}

double DesignElement::GetDesign() const
{
  EXCEPTION("use DesignElement::GetDesign(Access)");
}

double DesignElement::GetPhysicalDesign() const
{
  TransferFunction* tf = space_->GetTransferFunction(type_, ErsatzMaterial::ToApp(type_), true);

  return tf->Transform(this, SMART);

  //const TransferFunction* tf = const_cast<const TransferFunction*>(space_->GetTransferFunction(type_, ErsatzMaterial::MECH, true));
 // return tf->Transform(this);

  // we need the transder function
}

bool DesignElement::HasPhysicalDesign() const
{
  return(type_ == DENSITY || type_ == POLARIZATION);
}


void DesignElement::ToInfo(PtrParamNode in) const
{
  in->Get("type")->SetValue(type.ToString(type_));
  in->Get("upperBound")->SetValue(upper_);
  in->Get("lowerBound")->SetValue(lower_);
}

std::string DesignElement::ToString(const DesignElement* de)
{
  if(de == NULL) return "null";

  std::stringstream ss;
  if(de->elem == NULL)
  {
    ss << "ghost";
    if(de->vicinity != NULL) ss << de->vicinity->ToString();
  }
  else ss << boost::lexical_cast<std::string>(de->elem->elemNum);
  
  return ss.str();
}

std::string DesignElement::ToString(const StdVector<DesignElement*>& vec)
{
  std::stringstream ss;
  ss << "[";
  for(unsigned int i = 0; i < vec.GetSize(); i++)
    ss << ToString(vec[i]) << (i < vec.GetSize() - 1 ? "," : "");
  ss << "]";

  return ss.str();
}

void DesignElement::SetEnums()
{
  Filter::type.SetName("Filter::Type");
  Filter::type.Add(Filter::NO_FILTERING, "no_filtering");
  Filter::type.Add(Filter::SENSITIVITY, "sensitivity");
  Filter::type.Add(Filter::DENSITY, "density");

  Filter::sensitivity.SetName("Filter::Sensitivity");
  Filter::sensitivity.Add(Filter::PLAIN, "plain");
  Filter::sensitivity.Add(Filter::SHARP_PLAIN, "sharp_plain");
  Filter::sensitivity.Add(Filter::SIGMUND, "sigmund");
  Filter::sensitivity.Add(Filter::SHARP_SIGMUND, "sharp_sigmund");
  Filter::sensitivity.Add(Filter::BORRVALL, "borrvall");

  Filter::density.SetName("Filter::Density");
  Filter::density.Add(Filter::STANDARD, "standard");
  Filter::density.Add(Filter::HEAVISIDE, "heaviside");
  Filter::density.Add(Filter::MOD_HEAVISIDE, "mod_heaviside");

  type.SetName("DesignElement::Type");
  type.Add(TENSOR_TRACE, "tensor_trace");
  type.Add(DEFAULT, "default");
  type.Add(DENSITY, "density");
  type.Add(POLARIZATION, "polarization");
  type.Add(EMODUL, "emodul");
  type.Add(POISSON, "poisson");
  type.Add(LAMELAMBDA, "lamelambda");
  type.Add(LAMEMU, "lamemu");
  type.Add(EMODULISO, "emodul-iso");
  type.Add(POISSONISO, "poisson-iso");
  type.Add(GMODUL, "gmodul");
  type.Add(MASS, "mass");
  type.Add(DAMPINGALPHA, "damping-alpha");
  type.Add(DAMPINGBETA, "damping-beta");
  type.Add(UNITY, "unity");

  access.SetName("DesignElement::Access");
  access.Add(PLAIN, "plain");
  access.Add(SMART, "smart");

  valueSpecifier.SetName("DesignElement::ValueSpecifier");
  valueSpecifier.Add(DESIGN, "design");
  valueSpecifier.Add(COST_GRADIENT, "costGradient");
  valueSpecifier.Add(CONSTRAINT_GRADIENT, "constraintGradient");
  valueSpecifier.Add(MAX_SLOPE, "maxSlope");
  valueSpecifier.Add(MAX_OSCILLATION, "maxOscillation");
  valueSpecifier.Add(MAX_MOLE, "maxMole");
  valueSpecifier.Add(MAX_JUMP, "maxJump");
  valueSpecifier.Add(PENALIZED_STRESS, "penalizedStress");
  valueSpecifier.Add(WEIGHT, "weight");
  valueSpecifier.Add(OBJECTIVE, "objective");
  valueSpecifier.Add(NUM_NEIGHBOURS, "neighbours");
  valueSpecifier.Add(LEVEL_SET_VALUE, "levelSetValue");
  valueSpecifier.Add(LEVEL_SET_STATE, "levelSetState");
  valueSpecifier.Add(TOPGRAD_VALUE, "topGradValue");
  valueSpecifier.Add(SHAPEGRAD_VALUE, "shapeGradValue");
  valueSpecifier.Add(SHAPEGRAD_NODE_VALUE, "shapeGradNodeValue");
  valueSpecifier.Add(LEVEL_SET_GRAD_XP, "levelSetGradXP");
  valueSpecifier.Add(LEVEL_SET_GRAD_XN, "levelSetGradXN");
  valueSpecifier.Add(LEVEL_SET_GRAD_YP, "levelSetGradYP");
  valueSpecifier.Add(LEVEL_SET_GRAD_YN, "levelSetGradYN");
  valueSpecifier.Add(LEVEL_SET_GRAD_ZP, "levelSetGradZP");
  valueSpecifier.Add(LEVEL_SET_GRAD_ZN, "levelSetGradZN");

  detail.SetName("DesignElement::Detail");
  detail.Add(NONE, "none");
  detail.Add(SYMMETRY, "symmetry");
  detail.Add(FINITE_DIFF_COST_GRADIENT, "finiteDiffCostGrad");
  detail.Add(ERROR_COST_GRADIENT, "finiteDiffCostGradRelError");
  detail.Add(MECH_MECH, "mech_mech");
  detail.Add(ELEC_ELEC, "elec_elec");
  detail.Add(ELEC_ELEC_QUAD, "elec_elec_quad");
  detail.Add(ELEC_MECH, "elec_mech");
  detail.Add(MECH_ELEC, "mech_elec");
  // this ia a selection of constraints for constraintGradient
  detail.Add(COMPLIANCE, "compliance");
  detail.Add(VOLUME, "volume");
  detail.Add(PENALIZED_VOLUME, "penalizedVolume");
  detail.Add(GAP, "gap");
  detail.Add(REALVOLUME, "realvolume");
  detail.Add(TRACKING, "tracking");
  detail.Add(HOMOGENIZATION_TRACKING, "homTracking");
  detail.Add(POISSONS_RATIO, "poissonsRatio");
  detail.Add(YOUNGS_MODULUS, "youngsModulus");
  detail.Add(TYCHONOFF, "tychonoff");
  detail.Add(GREYNESS, "greyness");
  detail.Add(GLOBAL_SLOPE, "globalSlope");
  detail.Add(GLOBAL_CHECKERBOARD, "globalCheckerboard");
  detail.Add(STRESS, "stress");

}


SIMPElement::SIMPElement(DesignElement* base)
{
  this->de_ = base;
  this->filter = Filter();
  this->weight  = 1.0;
}

double SIMPElement::CalcWeightSum(bool include_this) const
{
  double res = 0.0;

  for(unsigned int i = 0, n = neighborhood.GetSize(); i < n; i++)
    res += neighborhood[i].weight;

  if(include_this)
    res += this->weight;

  return res;
}


double SIMPElement::GetSensitivityFilteredValue(DesignElement::ValueSpecifier sp, Condition* g) const
{
  // We filter over this element and the neighbors.
  assert(de_->simp != NULL);

  Filter* f = &de_->simp->filter;
  assert(f->type_ == Filter::SENSITIVITY);

  // See Filter::Sensitivity: w = weight, p is density, f' is cost gradient
  // Sigmund  = sum_i w(x_i) * p_i * f_i' / p_e * sum_i w(x_i)
  // Sharp Sigmund  = sum_i (i=e ? 1:0 : w(x_i)) * p_i * f_i' / p_e * sum_i w(x_i) + "bug" in normalized weighting
  // Borrvall = sum_i w(x_i) * p_i * f_i' / sum_i p_i * w(x_i)
  // plain    = sum_i w(x_i) * f_i' / sum_i w(x_i)
  // sharp plain = plain, just the filter is setup with normalized weights with a "bug"

  // factor design in numerator (SIGMUND and BORRVALL) for densitiy filtering the value is any in
  bool numerator_design = f->sensitivity_ != Filter::PLAIN && f->sensitivity_ != Filter::SHARP_PLAIN;
  // factor design in denominator sum (BORRVALL)
  bool denominator_design = f->sensitivity_ == Filter::BORRVALL;
  // weight the denominator by this design (SIGMUND)
  bool sigmund_denominator = f->sensitivity_ == Filter::SIGMUND || f->sensitivity_ == Filter::SHARP_SIGMUND;
  // short-cut for fake in Sharp Sigmund: i=e ? 1:0 : w(x_i)
  bool cheat_this_weight = f->sensitivity_ == Filter::SHARP_SIGMUND || f->sensitivity_ == Filter::SHARP_PLAIN;

  LOG_DBG3(del) << "FGV: el=" << de_->elem->elemNum
                << " sp=" << DesignElement::valueSpecifier.ToString(sp)
                << " dens=" << Filter::density.ToString(f->density_)
                << " numerator_design=" << numerator_design
                << " denominator_design=" << denominator_design
                << " sigmund_denominator=" << sigmund_denominator
                << " cheat_this_weight=" << cheat_this_weight;


  double numerator = 0.0;
  double denominator = 0.0;

  // mathematically the neighborhood includes this element, but this is not in the structure
  for(int i = -1, ni = (int) neighborhood.GetSize(); i < ni; i++)
  {
    const NeighbourElement* ne = i == -1 ? NULL : &neighborhood[i];
    const DesignElement* de = i == -1 ? this->de_ : ne->neighbour;

    double w = i == -1 ?  this->weight : ne->weight;
    double nw = cheat_this_weight && i == -1 ? 1.0 : w;
    double v = de->GetPlainValue(sp, g);
    double x = de->GetPlainValue(DesignElement::DESIGN); // cheap if not used

    double numerator_summand   = nw * v * (numerator_design ? x : 1.0);
    double denominator_summand = w * (denominator_design ? x : 1.0);

    numerator   += numerator_summand;
    denominator += denominator_summand;

    LOG_DBG3(del) << "FGV: el=" << de_->elem->elemNum << ": curr=" << de->elem->elemNum
                  << " w= " << w  << " nw=" << nw << " x=" << x << " v=" << v << " ns=" << numerator_summand
                  << " ds=" << denominator_summand << " num=" << numerator << " den=" << denominator;
  }

  if(sigmund_denominator)
    denominator *= de_->GetDesign(DesignElement::PLAIN);

  LOG_DBG3(del) << "FGV: el=" << de_->elem->elemNum << ": result=" << numerator
                << "/" << denominator << " -> " << (numerator / denominator);

  return numerator / denominator;
}


double SIMPElement::GetDensityFilteredValue(DesignElement::ValueSpecifier sp, Filter::Density fd) const
{
  // We filter over this element and the neighbors.
  assert(de_->simp != NULL);

  assert(sp == DesignElement::DESIGN);

  // All equations from Sigmund; Morphology based black and white filters for topology optimization; 2007
  // p = rho. P is filtered rho (rho tilde)
  // P = sum_(i in N_e) w(x_i) p_i / sum_(i in N_e) w(x_i)

  double numerator = 0.0;
  double denominator = 0.0;

  // mathematically the neighborhood includes this element, but this is not in the structure
  for(int i = -1, ni = (int) neighborhood.GetSize(); i < ni; i++)
  {
    const NeighbourElement* ne = i == -1 ? NULL : &neighborhood[i];
    const DesignElement* de = i == -1 ? this->de_ : ne->neighbour;

    double w = i == -1 ?  this->weight : ne->weight;
    double x = de->GetPlainValue(DesignElement::DESIGN); // cheap if not used

    numerator   += w * x;
    denominator += w;

    LOG_DBG3(del) << "GDFV: el=" << de_->elem->elemNum << ": curr=" << de->elem->elemNum
                  << " w= " << w  << " p=" << x << " num=" << numerator << " den=" << denominator;
  }

  double p_filt = numerator / denominator;

  LOG_DBG3(del) << "GDFV: el=" << de_->elem->elemNum << " filtered_density=" << p_filt;

  assert(fd == Filter::STANDARD || fd == Filter::HEAVISIDE || fd == Filter::MOD_HEAVISIDE);

  // do we need to post proc?
  if(fd != Filter::STANDARD)
  {
    p_filt = CalcHeaviside(p_filt);

    assert(p_filt <= this->de_->GetUpperBound());
    assert(p_filt >= 0.99 * this->de_->GetLowerBound()); // relax the assert a little, cause of heaviside correction
  }

  LOG_DBG3(del) << "GDFV: el=" << de_->elem->elemNum << " design=" << Filter::density.ToString(de_->simp->filter.density_)
                << ": plain=" << this->de_->GetPlainValue(DesignElement::DESIGN) << " -> "<< p_filt;

  return p_filt;
}

double SIMPElement::CalcHeaviside(double input_value) const
{
  Filter* f = &de_->simp->filter;
  assert(f->type_ == Filter::DENSITY);
  assert(f->density_ == Filter::HEAVISIDE || f->density_ == Filter::MOD_HEAVISIDE);

  double b = f->GetBeta();
  assert(b >= 0.0 && b < 2000);

  if(f->density_ == Filter::HEAVISIDE)
  {
    // we apply the correction factor in a way that H(rho_min) = rho_min and H(1) = 1
    double corr = (1.0 - (1.0 - input_value) * f->heaviside_corr) * input_value;
    double result = 1.0 - std::exp(-1.0 * b * corr) + corr * std::exp(-1.0 * b);

    LOG_DBG3(del) << "CH: de=" << de_->elem->elemNum << " f=" << f->density.ToString(f->density_)
                  << " hc=" << f->heaviside_corr << " corr=" << corr << " iv=" << input_value << " -> " << result;
    return result;
  }
  else // if(f->density_ == Filter::MOD_HEAVISIDE)
  {
    // make sure we are within the bounds
    double ub = this->de_->GetUpperBound();
    double lb = this->de_->GetLowerBound();

    double first    = std::exp(-1.0 * b * (1.0 - input_value));
    double second   = -1.0 * (1.0 - input_value) * std::exp(-1.0 * b);

    return (ub-lb) * (first + second) + lb;

    //LOG_DBG3(del) << "GDFV: el=" << de_->elem->elemNum << " b=" << b << " lb=" << lb << " (ub-lb)=" << (ub-lb)
    //              << " c=" << constant << " +1st=" << first << " +2nd=" << second << " =" << (constant + first + second) << " -> " << p_filt;
  }
}

double SIMPElement::GetDensityFilteredGradient(DesignElement::ValueSpecifier sp, Condition* g) const
{
  // We filter over this element and the neighbors.
  assert(de_->simp != NULL);
  Filter& f = de_->simp->filter;
  assert(f.type_ == Filter::DENSITY);
  assert(sp == DesignElement::COST_GRADIENT || sp == DesignElement::CONSTRAINT_GRADIENT);
  assert((sp == DesignElement::COST_GRADIENT && g == NULL) || (sp == DesignElement::CONSTRAINT_GRADIENT && g != NULL));
  assert(g == NULL || g->ForDensityFiltering());

  // Density filtering for gradient is (Sigmund; Morphology-based black and white filters for topology optimization; 2007; eqn (35). (36)
  // p is rho and P is rho filtered! d f/d p_e = sum_i(in N_e) d f/d P_i * d P_i/d p_e with d P_i/d p_e = w(x_e)/ sum_j(in N_i) w(x_j)
  // note, that the stored value is already v = d f/d P_i

  LOG_DBG3(del) << "GDFG: el=" << de_->elem->elemNum
                << " sp=" << DesignElement::valueSpecifier.ToString(sp)
                << " g=" << (g != NULL ? Condition::type.ToString(g->GetType()) : "null");

  double sum = 0.0;

  // mathematically the neighborhood includes this element, but this is not in the structure
  for(int i = -1, ni = (int) neighborhood.GetSize(); i < ni; i++)
  {
    const NeighbourElement* ne = i == -1 ? NULL : &neighborhood[i];
    const DesignElement* de = i == -1 ? this->de_ : ne->neighbour;

    double v = de->GetPlainValue(sp, g); // d f/d P_i
    double h = 1.0; // for not-standard filters this is the additional derivative
    double x_n = 0.0;

    if(f.density_ != Filter::STANDARD)
    {
      double b = f.GetBeta();

      // we need the filtered density -> but the real filtered value!!
      x_n = de->simp->GetDensityFilteredValue(DesignElement::DESIGN, Filter::STANDARD);

      if(f.density_ == Filter::HEAVISIDE)
      {
        // corr = (1 - (1 - x) * hc) * x;
        //      = x^2 - hc*x  + x
        // H = 1 - exp(-b * corr) + corr * exp(-b)

        // let the compiler optimize!
        double corr  = f.heaviside_corr * x_n * f.heaviside_corr * x_n  -  f.heaviside_corr * x_n + x_n;
        double dcorr = 2.0 * f.heaviside_corr * x_n - f.heaviside_corr + 1;

        h *= b * dcorr * exp(-b * corr) + dcorr * exp(-1.0*b);
      }
      if(f.density_ == Filter::MOD_HEAVISIDE)
      {
        // general scaling
        h = de->GetUpperBound() - de->GetLowerBound();
        h *= b * exp(b*(x_n-1.0)) + exp(-1.0*b);
      }
    }

    double w = i == -1 ? this->weight : ne->weight;
    double w_sum = de->simp->CalcWeightSum(true);

    double summand = v * h * w / w_sum;
    sum += summand;

    LOG_DBG3(del) << "GDFG: el=" << de_->elem->elemNum << ": curr=" << de->elem->elemNum
                  << " v= " << v  << " h=" << h << " w=" << w << " x_n=" << x_n << " w_sum=" << w_sum
                  << " summand=" << summand << " sum=" << sum;
  }

  return sum;
}

std::string SIMPElement::ToString() const
{
  std::stringstream ss;
  ss << "el=" << de_->elem->elemNum << " #n=" << neighborhood.GetSize() << "(";
  for(unsigned int i = 0; i < neighborhood.GetSize(); i++)
  {
    ss << " " << neighborhood[i].neighbour->elem->elemNum;
    //ss << " n_" << neighborhood[i].neighbour->elem->elemNum << "_w=" << neighborhood[i].weight;
    //ss << " n_" << neighborhood[i].neighbour->elem->elemNum << "_d=" << neighborhood[i].distance;
  }
  ss << ")";
  return ss.str();
}

void SIMPElement::Dump()
{
  double weight_sum = weight;
  double distance_avg = 0.0;
  for(unsigned int i = 0; i < neighborhood.GetSize(); i++)
  {
    weight_sum   += neighborhood[i].weight;
    distance_avg += neighborhood[i].distance;
  }
  distance_avg /= neighborhood.GetSize();

  std::cout << "\nelement: " << de_->elem->elemNum << " location " << de_->GetLocation()->ToString() 
            << " weight sum " << weight_sum << " this weight " << weight <<" distance avg " << distance_avg << std::endl;
  for(unsigned int i = 0; i < neighborhood.GetSize(); i++)
    std::cout << "  n[" << i << "]: elem " << neighborhood[i].neighbour->elem->elemNum << " location "
              << neighborhood[i].neighbour->GetLocation()->ToString()
              << " dist=" << neighborhood[i].distance << " w=" << neighborhood[i].weight << std::endl;
}


VicinityElement::VicinityElement()
{
  design.Resize(domain->GetGrid()->GetDim() == 2 ? 4 : 6);
  design.Init(NULL);
}


void VicinityElement::Init(DesignSpace* space, DesignStructure* structure)
{
  // do it only once
  if(space->data[0].vicinity != NULL)
    return;

  Grid* grid = domain->GetGrid();
  
  if(!space->IsRegular())
    throw Exception("A regular design domain is required to use VicinityElements");

  // eventually the barycenters are already calculated, we need them to identify the orientation
  // we will need the barycenters in FindNeibhborhood()
  for(unsigned int i = 0; i < space->regions.GetSize(); i++)
    grid->SetElementBarycenters(space->regions[i].regionId, false); // no updated coordinates

  // let CFS find the neighborhood of *all* elements. With some luck this was done
  // anyway already and we get it for free
  // This does not include periodic b.c. neighbors, this is done by structure
  grid->FindElementNeighorhood();

  // todo: Only LevelSet does not call this method with structure but always NULL
  bool periodic = structure != NULL ? structure->IsPeriodic() : false;

  // we are either linear or quadratic quadrilaterals or hexahedrons. This is the linear
  // limit of common nodes in the elem->neighborhood pair
  int common = grid->GetDim() == 2 ? 2 : 4;

  // We need the spacing of a element to detect periodic elements later
  StdVector<double> spacing; // the output
  Matrix<double> coords; // temporary
  Elem* elem = space->data[0].elem;
  domain->GetGrid()->GetElemNodesCoord(coords, elem->connect);
  elem->ptElem->GetEdgeLength(coords, spacing);

  // in the periodic case the neighborhood is enlarged
  StdVector<std::pair<Elem*, int> > enlarged_data;

  // traverse over all elements
  // double cfs elements for multiple design variables are not handled special
  for(unsigned int e = 0, n = space->data.GetSize(); e < n; e++)
  {
    DesignElement& de = space->data[e];
    de.vicinity = new VicinityElement();
    // here we store the neighbors in a sorted way
    StdVector<DesignElement*>& ve_data = de.vicinity->design; // has proper size of NULLs

    // reference case
    StdVector<std::pair<Elem*, int> >& neighbors = *(de.elem->neighborhood);
    // reuse the enlarged_data element for the periodic case only
    if(periodic)
    {
      if(structure->ExtendPeriodicNeighborhood(de.elem, common, enlarged_data))
        neighbors = enlarged_data; // in the non-periodic case there is no need to copy the element data
    }

    LOG_DBG(del) << "VE:Init elem=" << de.elem->elemNum << " neighbors=" << neighbors.ToString();

    for(unsigned int n = 0; n < neighbors.GetSize(); n++)
    {
      // we consider only direct (edge/face) neighbors
      if(neighbors[n].second < common) continue;
      // now we have to find the relative position of candidate
      Elem* candidate = neighbors[n].first;

      LOG_DBG3(del) << "VE:Init elem=" << de.elem->elemNum << " e.bc=" << de.elem->barycenter.ToString() << " e.dim="
                    << de.elem->ptElem->GetDim() << " e.r=" << de.elem->regionId << " n=" << n << " o.el=" << candidate->elemNum
                    << " o.bc=" << candidate->barycenter.ToString() << " o.dim=" << candidate->ptElem->GetDim() << " o.r=" << candidate->regionId;

      // if the neighbor is a surface element we don't want to play with it
      if(de.elem->ptElem->GetDim() != candidate->ptElem->GetDim())
        continue;
      // same if the region does not match. E.g. if there is fixed region not subject to optimization, e.g. bruggi_two_bar
      if((de.elem->regionId != candidate->regionId) && (space->regions.GetSize() == 1 || !space->Contains(candidate->regionId)))
        continue;
      assert(space->regions.GetSize() == 1); // extend for multi-region

      // the spacing allows to identify periodic elements
      int idx = FindRelativeNeighborLocation(de.elem->barycenter, candidate->barycenter, spacing);
      ve_data[idx] = space->Find(candidate->elemNum, de.GetType(), false);
      LOG_DBG2(del) << "VE:Init elem=" << de.elem->elemNum << " idx=" << idx << " dat=" << ve_data[idx];
    }
  }
 }


VicinityElement::Neighbour VicinityElement::FindRelativeNeighborLocation(Point& ref, Point& other, StdVector<double> spacing)
{
  // -------------------------
  // |       |       |       |
  // |   2   |  1    |   8   |
  // |       |  +y   |       |
  // -------------------------
  // |       |       |       |
  // |   3   |  0    |   7   |
  // |  -x   |       |  +x   |
  // -------------------------
  // |       |       |       |
  // |   4   |   5   |  6    |
  // |       |  -y   |       |
  // -------------------------
  //

  // we do not return after the first result is found.
  // This allows the asserts and also it should be of the
  // similar performance as the code can be vectorized
  Neighbour res = NONE;
  Point diff = ref - other;

  if(!IsNoise(diff[0]))
  {
    // the nodes can only be larger than the regular spacing in the periodic case
    if(std::abs(diff[0]) > spacing[0] * 1.5) // 1.5 reference size makes it robust
     res = diff[0] < 0 ? X_N : X_P; // flip due to periodic b.c.
    else
     res = diff[0] < 0 ? X_P : X_N; // normal, non-periodic case

    assert(IsNoise(diff[1]));
    assert(IsNoise(diff[2]));
  }
  if(!IsNoise(diff[1]))
  {
    assert(IsNoise(diff[0]));

    if(std::abs(diff[1]) > spacing[1] * 1.5)
      res = diff[1] < 0 ? Y_N : Y_P;
    else
      res = diff[1] < 0 ? Y_P : Y_N;

    assert(IsNoise(diff[2]));
  }
  if(!IsNoise(diff[2]))
  {
    assert(IsNoise(diff[0]));
    assert(IsNoise(diff[1]));

    if(std::abs(diff[2]) > spacing[2] * 1.5)
      res = diff[2] < 0 ? Z_N : Z_P;
    else
      res = diff[2] < 0 ? Z_P : Z_N;
  }
  
  LOG_DBG(del) << "VE:FRNL ref =" << ref.ToString() << " other=" << other.ToString() << " -> " << res;

  if(res == NONE)
    EXCEPTION("cannot identify relative neighborhood of " << ref.ToString() << " and " << other.ToString());

  return res;
}

DesignElement* VicinityElement::GetNeighbour(DesignElement* base, Neighbour idx, int n, bool throw_exception)
{
  assert(n > 0);

  DesignElement* tmp  = base;
  for(int i = 0; i < n; i++)
  {
    tmp = tmp->vicinity->GetNeighbour(idx);
    if(tmp == NULL)
    {
      LOG_DBG3(del) << "VE:GN base=" << base->ToString() << " idx=" << idx << " n=" << n << " max neighbor=" << i;
      if(throw_exception)
        EXCEPTION("no neighbor in " << idx << " direction " << i << " elements ways for element " << tmp->ToString())
      else
        return NULL;
    }
  }
  LOG_DBG3(del) << "VE:GN base=" << base->ToString() << " idx=" << idx << " n=" << n << " -> " << (tmp != NULL ? tmp->ToString() : "null");
  return tmp;
}

bool VicinityElement::HasNeighbor(DesignElement* base, Neighbour idx, int n)
{
  assert(n > 0);

  for(int i = 0; i < n; i++)
  {
    base = base->vicinity->GetNeighbour(idx);
    if(base == NULL) return false;
  }
  return true;
}


int VicinityElement::GetNumberOfEntries() const
{
  int count(0);
  for(unsigned int i = 0; i < design.GetSize(); ++i)
    if(design[i] != NULL) ++count;

  return count;
}

string VicinityElement::ToString() const
{
  stringstream ss;
  ss << "(";
  unsigned int max = domain->GetGrid()->GetDim() == 2 ? 4 : 6;
  for(unsigned int i = 0; i < max; i++)
  {
    if(design[i] == NULL) ss << "null";
    else
    {
      if(design[i]->elem == NULL) ss << "ghost";
                             else ss << design[i]->elem->elemNum;
    }
    if(i < max-1) ss << ",";
  }
  ss << ")";
  return ss.str();
}


ResultDescription::ResultDescription()
{
  // does not set all!!
  access = DesignElement::PLAIN;
  value  = DesignElement::DESIGN;
  design = DesignElement::DEFAULT;
  excitation = "";
}

ResultDescription::ResultDescription(PtrParamNode pn)
{
  solutionType = SolutionTypeEnum.Parse(pn->Get("id")->As<std::string>());

  design = DesignElement::DEFAULT;
  if(pn->Has("design"))
    design = DesignElement::type.Parse(pn->Get("design")->As<std::string>());

  access = DesignElement::access.Parse(pn->Get("access")->As<std::string>());

  value = DesignElement::valueSpecifier.Parse(pn->Get("value")->As<std::string>());

  detail = DesignElement::detail.Parse(pn->Get("detail")->As<std::string>());

  excitation = pn->Get("excitation")->As<std::string>();
}
