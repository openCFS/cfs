#include <cmath>
#include <iostream>
#include <utility>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Domain/Domain.hh"
#include "Domain/ElemMapping/Elem.hh"
#include "Domain/Mesh/Grid.hh"
#include "FeBasis/BaseFE.hh"
#include "General/defs.hh"
#include "General/Exception.hh"
#include "MatVec/Matrix.hh"
#include "Optimization/Condition.hh"
#include "Optimization/Excitation.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Design/ShapeMapDesign.hh"
#ifdef USE_EMBEDDED_PYTHON 
  #include "Optimization/Design/SpaghettiDesign.hh"
#endif
#include "Optimization/Design/SplineBoxDesign.hh"
#include "Optimization/Design/DesignStructure.hh"
#include "Optimization/Function.hh"
#include "Optimization/LevelSet.hh"
#include "Optimization/Objective.hh"
#include "Optimization/Context.hh"
#include "Optimization/TransferFunction.hh"
#include "PDE/SinglePDE.hh"
#include "Utils/Point.hh"
#include "Utils/tools.hh"
#include "Utils/mathParser/mathParser.hh"
#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/math/special_functions/sign.hpp"
#include "boost/lexical_cast.hpp"


using namespace CoupledField;

using std::string;
using boost::math::sign;
using boost::posix_time::ptime;
using boost::posix_time::second_clock;
using boost::posix_time::microsec_clock;

DEFINE_LOG(desel, "designElement")

// the static enum
Enum<Filter::Type>                  Filter::type;
Enum<Filter::Density>               Filter::density;
Enum<Filter::Sensitivity>           Filter::sensitivity;
Enum<Filter::FilterSpace>           Filter::filterSpace;
Enum<BaseDesignElement::Type>       BaseDesignElement::type;
Enum<DesignElement::ValueSpecifier> DesignElement::valueSpecifier;
Enum<DesignElement::Access>         DesignElement::access;
Enum<DesignElement::Detail>         DesignElement::detail;
Enum<ShapeParamElement::Dof>        ShapeParamElement::dof;

Enum<FeatureVariable::Tip> FeatureVariable::tip_enum;

Enum<ShapeMapDesign::Type>          ShapeMapDesign::type;

// is a static attribute
DesignSpace* DesignElement::space_(NULL);

string BaseDesignElement::ToString() const
{
  std::stringstream ss;
  ss << "idx=" << index_ << " t=" << type.ToString(type_);
  return ss.str();
}


bool BaseDesignElement::IsCompatible(Type super, Type test)
{
  switch(super)
  {
  case MECH_TRACE:
  {
    switch(test)
    {
    case MECH_11:
    case MECH_22:
    case MECH_33:
    case MECH_44:
    case MECH_55:
    case MECH_66:

      // Tensor trace for param mat
    case STIFF1:
    case STIFF2:
    case STIFF3:
    case SHEAR1:
    //for mod_red
    case ROTANGLE:
    // Batian's stuff
    // FIXMI!!
    case POISSON:
    case POISSONISO:
    case EMODUL:
    case EMODULISO:
      return true;
    default:
      return false;
    }
    break;
  }

  case MECH_ALL:
  {
    switch(test)
    {
    case MECH_11:
    case MECH_12:
    case MECH_13:
    case MECH_14:
    case MECH_15:
    case MECH_16:
    case MECH_22:
    case MECH_23:
    case MECH_24:
    case MECH_25:
    case MECH_26:
    case MECH_33:
    case MECH_34:
    case MECH_35:
    case MECH_36:
    case MECH_44:
    case MECH_45:
    case MECH_46:
    case MECH_55:
    case MECH_56:
    case MECH_66:
      return true;
    default:
      return false;
    }
    break;
  }

  case DIELEC_TRACE:
  {
    switch(test)
    {
    case DIELEC_11:
    case DIELEC_22:
      return true;
    default:
      return false;
    }
    break;
  }

  case DIELEC_ALL:
  {
    switch(test)
    {
    case DIELEC_11:
    case DIELEC_12:
    case DIELEC_22:
      return true;
    default:
      return false;
    }
    break;
  }

  case PIEZO_ALL:
  {
    switch(test)
    {
    case PIEZO_11:
    case PIEZO_12:
    case PIEZO_13:
    case PIEZO_21:
    case PIEZO_22:
    case PIEZO_23:
      return true;
    default:
      return false;
    }
    break;
  }

  case DEFAULT:
  case ALL_DESIGNS:
    return true;
    break;

  default:
    if(super == test)
      return true;
    else
      return false;
  }

}

BaseDesignElement::Type BaseDesignElement::MapSolutionType(SolutionType soltype, bool throw_exception)
{
  switch(soltype)
  {
  case MECH_PSEUDO_DENSITY:
  case PSEUDO_DENSITY:
  case PHYSICAL_PSEUDO_DENSITY:
    return DENSITY;

  case RHS_PSEUDO_DENSITY:
  case PHYSICAL_RHS_PSEUDO_DENSITY:
    return RHS_DENSITY;

  default:
    break;
  }

  if(throw_exception)
    EXCEPTION("no DesignElement::Type for SolutionType " << soltype);

  return NO_TYPE;
}

bool BaseDesignElement::IsPhysical(SolutionType soltype)
{
  switch(soltype)
  {
  case PHYSICAL_PSEUDO_DENSITY:
  case PHYSICAL_RHS_PSEUDO_DENSITY:
  case ELEC_PHYSICAL_PSEUDO_DENSITY:
    return true;

  default:
    return false;
  }
}

BaseDesignElement::BaseDesignElement(Type t) {
  design          = 0.0;
  upper_          = 0.0;
  lower_          = 0.0;
  type_           = t;
  index_          = std::numeric_limits<unsigned int>::max();
}


/** Get the gradient values for either objective or constraint */
double BaseDesignElement::GetPlainGradient(const Function* f) const
{
  assert(f != NULL); // Call SumObjectiveGradient() if you don't want to specify!

  assert(!f->IsObjective() || (f->IsObjective() && dynamic_cast<const Objective*>(f) != NULL));
  assert( f->IsObjective() || (!f->IsObjective() && dynamic_cast<const Condition*>(f) != NULL));

  LOG_DBG3(desel) << "GPG idx=" << this->index_ << " v=" << design << " f=" << f->ToString() << " dcost=" << costGradient.ToString() << " dconst=" << constraintGradient.ToString();

  return f->IsObjective() ? costGradient[f->GetIndex()] : constraintGradient[f->GetIndex()];
}

double BaseDesignElement::GetPlainGradient(const Objective* c) const
{
  return costGradient[c->GetIndex()];
}

double BaseDesignElement::GetPlainGradient(const Condition* g) const
{
  return constraintGradient[g->GetIndex()];
}

/** Sum app the old value (get and set together) */
void BaseDesignElement::AddGradient(const Objective* f, const Condition* g, double value)
{
  assert(f == NULL || g == NULL);
  LOG_DBG3(desel) << "AddGradient: f=" << (f == NULL ? "null" : f->type.ToString(f->GetType()))
                << " g=" << (g == NULL ? "null" : g->ToString()) << " val=" << value
                << " penalty=" << (f != NULL ? std::to_string(f->scale) : "-")
                << " old= " <<  (f != NULL ? costGradient[f->GetIndex()] : constraintGradient[g->GetIndex()])
                << " add=" << (f != NULL ? value * f->scale : value)
                << " -> " << (f != NULL ? costGradient[f->GetIndex()] + value * f->scale : constraintGradient[g->GetIndex()] + value);
  if(f != NULL)
    costGradient[f->GetIndex()] += value * f->scale;
  else
    constraintGradient[g->GetIndex()] += value;
}

void BaseDesignElement::AddGradient(const Function* f, double d_value)
{
  assert(!std::isnan(d_value));
  assert(!std::isinf(d_value));
  assert(( f->IsObjective() && dynamic_cast<const Objective*>(f) != NULL)
      || (!f->IsObjective() && dynamic_cast<const Condition*>(f) != NULL) );

  if(f->IsObjective())
  {
    const Objective* c = static_cast<const Objective*>(f);
    // we need to differentiate between c->GetValue() and c->GetOrgValue() -> See Optimization::CalcObjective()
    // there we have c->SetValue(weight * tov), therefore we need for all term but LINEAR the c->GetOrgValue()
    // but in that case we also need to add the weight manually -> so check if you need to implement it!
    double res = 0.0;

    switch(c->GetTerm())
    {
     case Objective::LINEAR:
       res =  c->scale * d_value;
       break;
     case Objective::PENALTY:
       // scale * max(0, value - parameter)^2
       res = c->scale * 2.0 * std::max(0.0, c->GetOrgValue() - c->GetParameter()) * d_value;
       break;
     case Objective::LN1P:
       // f = scale * sign(v(x)) * ln(1 + |v(x)|)
       res = c->scale * (sign(c->GetOrgValue()) * d_value) / (1 + std::abs(c->GetOrgValue()));
       break;
     case Objective::POWER:
       // scale * value^parameter
       res = c->scale * c->GetParameter() * std::pow(c->GetOrgValue(), c->GetParameter()-1.0) * d_value;
       break;
    }
    LOG_DBG3(desel) << "AG: i=" << this->index_ << " o=" << c->type.ToString(f->GetType()) << " dv="  << d_value << " v=" << c->GetValue()
                    << " s=" << c->scale << " t=" << c->term.ToString(c->GetTerm()) << " -> " << res;
    costGradient[f->GetIndex()] += res;

  }
  else
    constraintGradient[f->GetIndex()] += d_value;
}

void BaseDesignElement::Reset(ValueSpecifier vs, Function*  f)
{
  if(vs == FUNCTION_GRADIENT)
  {
    assert(f != nullptr);
    vs = f->IsObjective() ? COST_GRADIENT : CONSTRAINT_GRADIENT;
  }

  switch(vs)
  {
  case COST_GRADIENT:
    for(unsigned int i = 0, s = costGradient.GetSize(); i < s; ++i)
      costGradient[i] = 0.0;
    break;
  case CONSTRAINT_GRADIENT:
    if(f != nullptr)
      constraintGradient[f->GetIndex()] = 0.0;
    else
      for(unsigned int i = 0, s = constraintGradient.GetSize(); i < s; ++i)
        constraintGradient[i] = 0.0;
    break;
  default:
    EXCEPTION("invalid value specifier " << vs);
  }
}

double BaseDesignElement::SumObjectiveGradient() const
{
  double result = 0.0;
  for(unsigned int i = 0, s = costGradient.GetSize(); i < s; ++i)
    result += costGradient[i];

  return result;
}

string BaseDesignElement::ToString(const StdVector<BaseDesignElement*>& vec, bool print_type)
{
  std::stringstream ss;
  ss << "[";
  for(unsigned int i = 0, s = vec.GetSize(); i < s; ++i) {
    if(typeid(vec[i]) == typeid(DesignElement*)){
      ss << DesignElement::ToString(dynamic_cast<DesignElement*>(vec[i]));
    }else{
      ss << " BaseDesignElement ";
    }
    if(print_type) ss << "=" << vec[i]->type_;
    if(i < s-1) ss << ",";
  }
  ss << "]";

  return ss.str();
}

/** Bastian's shape optimization */
ShapeDesignElement::ShapeDesignElement(unsigned int index) : BaseDesignElement() {
  index_ = index;
}

/* Parametric shape optimization */
ShapeParamElement::ShapeParamElement(Type type, unsigned int index) : BaseDesignElement(type)
{
  index_ = index;
}

string ShapeParamElement::ToString() const
{
  std::stringstream ss;
  ss << "(idx=" << index_ << " opt_idx=" << ((int) opt_index_) << " t=" << type.ToString(type_) << " d=" << dof.ToString(dof_) << " v=" << design << ")";
  return ss.str();
}

string ShapeParamElement::GetLabel() const
{
  return type.ToString(type_) +  "_" + dof.ToString(dof_) + "_" + std::to_string(opt_index_);
}

void FeatureVariable::CompareToInfoHelper(const FeatureVariable* v0, const FeatureVariable* test, std::string& fixed, std::string& lower, std::string& upper) 
{
  assert(v0 != nullptr && test != nullptr);
  assert(v0->fixed == test->fixed);

  if(v0->GetPlainDesignValue() == test->GetPlainDesignValue()) {
    if(fixed != "not constant")
      fixed = std::to_string(v0->GetPlainDesignValue());
  } else
    fixed = "not constant";

  if(v0->GetLowerBound() == test->GetLowerBound()) {
    if(lower != "not constant")
      lower = std::to_string(v0->GetLowerBound());
  } else
    lower = "not constant";

  if(v0->GetUpperBound() == test->GetUpperBound()) {
    if(upper != "not constant")
      upper = std::to_string(v0->GetUpperBound());
  } else
    upper = "not constant";
}


bool FeatureVariable::IsFixed(const StdVector<FeatureVariable>& point) 
{
  assert(point.GetSize() == domain->GetDim());
  for(const FeatureVariable& v : point)
    if(v.fixed)
      return true;
  return false;
}

int FeatureVariable::CountRealVariables(const StdVector<FeatureVariable>& point) 
{
  int sum = 0;
  assert(point.GetSize() == domain->GetDim());
  for(const FeatureVariable& v : point)
    sum += v.IsVariable() ? 1 : 0;
  return sum;
}

std::string FeatureVariable::ToString(const StdVector<FeatureVariable>& vec, bool show_fixed)
{
  std::stringstream ss;
  ss << "[";
  for(unsigned int i = 0; i < vec.GetSize(); i++)
  {
    const FeatureVariable& v = vec[i];
    if(show_fixed && v.fixed)
      ss << "*";
    ss << v.GetPlainDesignValue();
    if(i < vec.GetSize()-1)
      ss << ", ";
  }
  ss << "]";
  return ss.str();
}

Vector<double> FeatureVariable::AsVector(const StdVector<FeatureVariable>& vec)
{
  Vector<double> res(vec.GetSize());
  for(unsigned int i = 0; i < vec.GetSize(); i++)
    res[i] = vec[i].GetPlainDesignValue();
  return res;
}


void FeatureVariable::Parse(PtrParamNode pn, int feature_idx, Type var, double interpolate_value)
{
  this->feature = feature_idx;
  this->var = var;

  MathParser* mp = domain->GetMathParser();
  unsigned int handle = mp->GetNewHandle();

  // mapping via "key" and "map" needs to be validated and processed outside (e.g. FeatureMappingDesign)
  if(pn->Has("key")) 
    key = pn->Get("key")->As<string>();
  if(pn->Has("map")) 
    map = pn->Get("map")->As<string>();

  fixed = false;  
  string value; // subject to mathparser evaluation

  if(pn->Has("fixed")) 
  {
    value = pn->Get("fixed")->As<string>();

    if(pn->Has("key") || pn->Has("map")) // only pill has key/map, not spaghetti - so don't mix in messages
      throw Exception("When a '" + pn->GetName() + "' is 'fixed', don't use 'key' or 'map'");
    if(pn->Has("upper") || pn->Has("lower") || pn->Has("initial"))
      throw Exception("When a '" + pn->GetName() + "' is 'fixed', don't use 'initial', 'lower' or 'upper'");
    fixed = true;
  } 
  if(!fixed)
  {
    if(map == "")
    { 
      if(!pn->Has("upper") || !pn->Has("lower") || !pn->Has("initial"))
        throw Exception("'" + pn->GetName() + "' needs 'initial', 'lower' and 'upper' when not 'fixed'");
      value = pn->Get("initial")->As<string>();
      mp->SetExpr(handle, pn->Get("lower")->As<string>());
      SetLowerBound(mp->Eval(handle));
      mp->SetExpr(handle, pn->Get("upper")->As<string>());
      SetUpperBound(mp->Eval(handle));
    }
    else
      if(pn->Has("upper") || pn->Has("lower") || pn->Has("initial") || pn->Has("fixed"))
        throw Exception("When a '" + pn->GetName() + "' is mapped, don't use 'initial', 'lower', 'upper', or 'fixed'");
  }

  // this might contain a formula or simply a value
  if(value == "interpolated")
  {
    assert(interpolate_value != -12.34);
    SetDesign(interpolate_value);
  }
  else
  {
    assert(value != "" || map != ""); // either value or map must be set
    assert(!(value != "" && map != "")); // not both
    if(value != "")
    {
      mp->SetExpr(handle, value);
      SetDesign(mp->Eval(handle));
    }
  }

  mp->ReleaseHandle(handle);

  if(pn->Has("tip")) // node requires tip
    tip = FeatureVariable::tip_enum.Parse(pn->Get("tip")->As<string>());

  if(pn->Has("dof")) // node requires top
    dof_ = dof.Parse(pn->Get("dof")->As<string>());

  type_ = type.Parse(pn->GetName());

  LOG_DBG(desel) << "FV:P fi=" << feature_idx << " v=" << value << " f=" << fixed << " tip=" << tip << " dof=" << dof_ << " type=" << type_;

  if(type_ != NODE && (tip != NO_TIP || dof_ != NOT_SET))
    throw Exception("Within 'pill' use 'dof' and 'tip' only for 'node'");
}

void FeatureVariable::ToInfo(PtrParamNode in) const
{
  in->Get("type")->SetValue(type.ToString(type_));
  if(type_ == NODE)
    in->Get("dof")->SetValue(ShapeParamElement::dof.ToString(dof_));
  if(type_ == NODE && tip != NO_TIP) 
    in->Get("tip")->SetValue(FeatureVariable::tip_enum.ToString(tip));

  if(fixed)
   in->Get("fixed")->SetValue(GetPlainDesignValue());
  else 
  {
   // key and map are only used for pills, not for noodles - shape map has own legacy stuff
   if(map == "") 
   {
     in->Get("lower")->SetValue(lower_);
     in->Get("upper")->SetValue(upper_);
     if(key != "")
       in->Get("key")->SetValue(key);
   }
   else
    in->Get("map")->SetValue(map);
  }
}

std::string FeatureVariable::ToString() const
{
  std::stringstream ss;
  ss << ShapeParamElement::ToString() << " feature=" << feature << " type=" << type.ToString(type_) << " dof=" << dof.ToString(dof_) << " tip=" << tip;
  if(map != "" )
    ss << " map=" << map;
  if(key != "")
    ss << " key=" << key; 
  return ss.str();
}

void FeatureVariable::InitInnerNode(int feature_idx, Dof dof, double val, double lower, double upper)
{
  this->type_ = NODE;
  this->fixed = false;
  this->feature = feature_idx;
  this->dof_ = dof;
  this->tip = INNER;
  SetDesign(val);
  SetLowerBound(lower);
  SetUpperBound(upper);
  LOG_DBG(desel) << "V:IIN " << ToString();
}

std::string FeatureVariable::GetLabel() const
{
  std::stringstream ss;
  ss << "s" << feature << "_"; // 0-based
  switch(type_)
  {
  case NODE:
    if(tip == START)
      ss << "p";
    else if(tip == END)
      ss << "q";
    else // We increment after x/y/z. 
      assert(false);
    ss << dof.ToString(dof_);
    break;
  case PROFILE:
    ss << "p";
    break;
  case NORMAL:
  case RADIUS:
  {
#ifdef USE_EMBEDDED_PYTHON 
    // we are in a noodle
    assert(domain->GetOptimization() != NULL);
    SpaghettiDesign* sd = dynamic_cast<SpaghettiDesign*>(domain->GetOptimization()->GetDesign());
    const SpaghettiDesign::Noodle& n = sd->spaghetti[feature];
    ss << (type_ == NORMAL ? "a" : "r");
    assert(!(n.a.GetSize() > 0 && n.r.GetSize() > 0)); // we are in GetLable() and shall have it
    unsigned int pos = GetIndex() - (type_ == NORMAL ? n.a.First().GetIndex() : n.r.First().GetIndex());
    ss << (pos + 1); // 1-based
#endif    
    break;    
  }
  default:
    assert(false);
    break;
  }

  return ss.str();
}

/** The default constructor for StdVector and ghost elements*/
DesignElement::DesignElement() : BaseDesignElement()
{
  Init();
  this->interfaceDrivenLoadGrad_.Resize(4 * (domain->GetDim()-1),0.0);
}

/** The default constructor for StdVector and ghost elements*/
DesignElement::DesignElement(Elem* elem, Type type, unsigned int index, int pseudoElementIndex) : BaseDesignElement()
{
  Init();
  this->elem = elem;

  if(!elem->extended)
    this->elem->extended = new ExtendedElementInfo;

  this->type_ = type;
  this->index_ = index;
  this->pseudoElementIndex_ = pseudoElementIndex;
  this->upper_ = 1.0;
  this->lower_ = 1.0;
  this->specialResult.Resize(OPT_RESULT_BOUND - OPT_RESULT_1, 0.0);
  this->interfaceDrivenLoadGrad_.Resize(4 * (domain->GetDim()-1),0.0);
}


DesignElement::DesignElement(Type dt, double lower, double upper, Elem* elem, unsigned int index, MultiMaterial* mm) : BaseDesignElement()
{
  Init();
  this->elem = elem;

  if(!elem->extended)
    this->elem->extended = new ExtendedElementInfo;

  this->specialResult.Resize(OPT_RESULT_BOUND - OPT_RESULT_1, 0.0);
  this->index_ = index;
  this->multimaterial = mm;
  this->interfaceDrivenLoadGrad_.Resize(4 * (domain->GetDim()-1),0.0);

  type_ = dt;
  upper_ = upper;
  lower_ = lower;
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
  vicinity        = NULL;
  lse_            = NULL;
  tge             = NULL;
  location_       = NULL;
  elem            = NULL;
  type_           = NO_TYPE;
  pseudoElementIndex_ = -1;
  elemVol_        = -1.0;
  elemPorosity_   = -1.0;
}


DesignElement::Type DesignElement::Default(const Context* ctxt)
{
  switch(ctxt->ToApp()) // will fail for piezo ?!
  {
  case App::MECH:
  case App::MAG:
  case App::BUCKLING:
  case App::ACOUSTIC:
    return DENSITY;
  case App::ELEC:
    return POLARIZATION;
  default:
    throw Exception("invalid");
  }
}


const Point* DesignElement::GetLocation()
{
  if(location_ != NULL) return location_;

  // calc location
  location_ = new Point();
  // check for ghost elements
  if(elem == NULL) EXCEPTION("location_ not set but elem is NULL");

  domain->GetGrid()->GetElemShapeMap(elem, false)->CalcBarycenter(*location_);

  LOG_DBG3(desel) << "DesignElement::GetLocation() find " << location_->ToString() << " for " << ToString();
  return location_;
}

double DesignElement::CalcVolume()
{
  // precalculated?
  if(elemVol_ >= 0)
    return elemVol_;

  elemVol_ = domain->GetGrid()->GetElemShapeMap(elem, false)->CalcVolume();
  return elemVol_;
}

void DesignElement::SetElemPorosity(double vol)
{
  elemPorosity_ = vol;
}

double DesignElement::GetElemPorosity() {
  return elemPorosity_;
}

unsigned int DesignElement::GetElementSolutionIndex() const
{
  // easy case is pseudo design element
  if(pseudoElementIndex_ != -1)
    return pseudoElementIndex_;

  // if we have one design it is easily index_ - we use the modulus
  return index_ % space_->GetNumberOfElements();
}

int DesignElement::GetOptResultIndex(SolutionType st)
{
  if(st < OPT_RESULT_1 || st >= OPT_RESULT_BOUND)
    return -1;
  else
    return static_cast<int>(st)-static_cast<int>(OPT_RESULT_1); // OPT_RESULT_1 -> 0
}

void DesignElement::GetValue(ResultDescription& rd, StdVector<double>& out, unsigned int dofs) const
{
  // check for special result1
  if(    rd.value == OBJECTIVE
      || (rd.value == COST_GRADIENT && rd.detail != NONE)
      || rd.value == CONSTRAINT_GRADIENT
      || rd.value == MAX_SLOPE
      || rd.value == MAX_OSCILLATION
      || rd.value == MAX_MOLE
      || rd.value == MAX_JUMP
      || rd.value == QUADRATIC_VM_STRESS
      || rd.value == LOCAL_LOAD_FACTOR
      || rd.value == DESIGN_TRACKING
      || rd.value == PROJECTION
      || rd.value == SHAPE_MAP_GRAD
      || rd.value == SHAPE_MAP_ORDER
      || rd.value == SHAPE_MAP_CORNER
      || rd.value == MMA_ASYMPTOTE
      || rd.value == MMA_LOWER_VAL
      || rd.value == MMA_UPPER_VAL
      || rd.value == MMA_OBJ_GRADIANT
      || rd.value == MMA_CON_GRADIANT_1
      || rd.value == MMA_CON_GRADIANT_2
      || rd.value == SPLINE_BOX_GRAD_X
      || rd.value == SPLINE_BOX_GRAD_Y
      || rd.value == SPLINE_BOX_GRAD_Z
      || rd.value == SPLINE_BOX_INT_ORDER
      || rd.value == SPLINE_BOX_INT_CORNER
      || rd.value == GENERIC_ELEM
      || rd.value == FILTERED_DESIGN
      || rd.value == DIFF_FILTERED_DESIGN
      || rd.value == FEATURE_GRAD)
  {
    if(dofs != 1) throw Exception("special results is only defined for scalar values");
    // note, that on EACH_FORWARD/ADJOINT we need excitation based results
    int idx = GetOptResultIndex(rd.solutionType);

    if(idx >= 0)
      out[0] = specialResult[idx];
    else
     throw Exception("solution type not handled");
  }
  else
  {
    if(dofs == 1)
    {
      if(rd.solutionType == PHYSICAL_PSEUDO_DENSITY)
        out[0] = GetPhysicalDesign(Optimization::context);
      else if(rd.solutionType == ELEC_PHYSICAL_PSEUDO_DENSITY || rd.solutionType == LBM_PHYSICAL_PSEUDO_DENSITY)
        out[0] = GetPhysicalDesign(Optimization::context);
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


double DesignElement::GetValue(ValueSpecifier vs, Access access, Function* f) const
{
  // pessimistic assumption :)
  bool sens_filter = false;
  bool design_filter = false;
  bool design_filter_grad = false;

  unsigned int fix = simp != NULL ? simp->DetermineFilterIndex() : 0;

  if(access == SMART && simp != NULL && (simp->filter.GetSize() > 0 && simp->filter[fix].GetType() != Filter::NO_FILTERING))
  {
    if(simp->filter[fix].GetType() == Filter::DENSITY)
    {
      if(vs == DesignElement::DESIGN)
        design_filter = true;
      if(vs == DesignElement::COST_GRADIENT)
        design_filter_grad = f != NULL ? f->isFiltered() : true;
      if(vs == DesignElement::CONSTRAINT_GRADIENT)
        design_filter_grad = f != NULL ? f->isFiltered() : true;
    }
    else
    {
      if(vs == DesignElement::COST_GRADIENT)// TODO also check for separate costs! not only sum!
        sens_filter = true;
      if(vs == DesignElement::CONSTRAINT_GRADIENT)
        sens_filter = f != NULL ? f->isFiltered() : false;
    }
  }

  // we are silent if one wants filtering and it is not possible
  double val = 0.0;
  if(sens_filter)
    val = simp->GetSensitivityFilteredValue(vs, f);
  if(design_filter)
    val = simp->GetDensityFilteredValue(vs, simp->filter[fix].global->density);
  if(design_filter_grad)
    val = simp->GetDensityFilteredGradient(vs, f);
  if(!sens_filter && !design_filter && !design_filter_grad)
    val = GetPlainValue(vs, dynamic_cast<Condition*>(f));

  LOG_DBG3(desel) << "GV: " << elem->elemNum << " (" << valueSpecifier.ToString(vs) << ", "
                << (access == PLAIN ? "plain" : "smart") << ", " << (f == NULL ? "null" : f->ToString())
                << ") sf=" << sens_filter << " df=" << design_filter << " dfg=" << design_filter_grad << " -> " << val;
  return val;
}



#ifdef _WIN32
  inline double DesignElement::GetPlainValue(ValueSpecifier sp, Condition* g) const
#else
  __attribute__((always_inline)) inline double DesignElement::GetPlainValue(ValueSpecifier sp, Condition* g) const
#endif
{
  // validate first:
  switch(sp)
  {
  case DESIGN:
    return design;

  case COST_GRADIENT:
    return GetPlainCostGradient();

  case WEIGHT:
    if(simp == NULL)
      throw Exception("'" + valueSpecifier.ToString(sp) + "' is specific to SIMP");
    return simp->DetermineFilter().weight;

  case NUM_NEIGHBOURS:
    if(simp == NULL)
      throw Exception("'" + valueSpecifier.ToString(sp) + "' is specific to SIMP");
    return simp->DetermineFilter().neighborhood.GetSize();

  case CONSTRAINT_GRADIENT:
    assert(g != NULL);
    return constraintGradient[g->GetIndex()];

  case TOPGRAD_VALUE:
    if(tge == NULL)
      throw Exception("'" + valueSpecifier.ToString(sp) + "' is specific to TopGrad");
    return tge->value;

  case SHAPEGRAD_VALUE:
    if(lse_ == NULL)
      throw Exception("'" + valueSpecifier.ToString(sp) + "' is specific to Levelset");
    return lse_->shapeGradValue;

  default:
    throw Exception(valueSpecifier.ToString(sp) + " shall either be special result index or is no scalar value");
  }
  return -1.0; // cannot happen due to default in switch but to please the compiler :(
}

inline double DesignElement::GetPlainCostGradient() const {
  Optimization* opt = domain->GetOptimization();
  assert(opt != NULL);

  double beta;

  // this is also true for single objective
  if(opt->GetMOType(beta) == Function::WEIGHTED_SUM)
    return SumObjectiveGradient();

  // get values of all objectives
  // opt->CalcObjective() has to been had called before, such that the values are set!
  assert(opt->CalcObjectiveCalled());
  size_t n = opt->objectives.data.GetSize();
  Vector<double> ov(n);
  for(size_t i = 0; i < n; ++i)
    ov[i] = opt->objectives.data[i]->GetValue();
  double sum = 0.0;

  // GetMOType writes to beta
  if(opt->GetMOType(beta) == Function::SMOOTH_MIN)
    for(unsigned int i=0; i < costGradient.GetSize(); ++i)
      sum += DerivSmoothMin(ov, beta, i) * costGradient[i];

  else if(opt->GetMOType(beta) == Function::SMOOTH_MAX)
    for(unsigned int i=0; i < costGradient.GetSize(); ++i)
      sum += DerivSmoothMax(ov, beta, i) * costGradient[i];

  return sum;
}

double DesignElement::GetPlainMechTrace() const
{

  assert(space_ != NULL);

  double sum = 0.0;

  for(unsigned int d = 0; d < space_->design.GetSize(); d++)
  {
    if(IsCompatible(MECH_TRACE, space_->design[d].design) && IsCompatible(MECH_ALL, space_->design[d].design))
    {
      BaseDesignElement* de = space_->Find(elem->elemNum, space_->design[d].design, true);
      assert(de->GetType() == space_->design[d].design);
      sum += de->GetPlainDesignValue();
    }
  }

  return sum;
}



double DesignElement::GetDesign(Access access) const
{
  return GetValue(DESIGN, access);
}

double DesignElement::GetDesign(const Function* f, const TransferFunction* tf) const
{
  assert(f != NULL);
  assert(f->GetAccess() == Function::PHYSICAL || f->GetAccess() == Function::FILTERED || f->GetAccess() == Function::PLAIN);

  if(f->GetAccess() == Function::PHYSICAL && tf == NULL)
    return GetPhysicalDesign();

  // smart means that we filter when  filter is activated. Desired for PYHYSICAL and FILTERED
  DesignElement::Access dac = f->GetAccess() == Function::PLAIN ? DesignElement::PLAIN : DesignElement::SMART;

  double v = GetValue(DESIGN, dac);
  if(f->GetAccess() != Function::PHYSICAL)
    return v;
  assert(tf != NULL);
  return tf->Transform(v);
}

double DesignElement::GetDesign() const
{
  EXCEPTION("use DesignElement::GetDesign(Access)");
}

double DesignElement::GetPhysicalDesign(const Context* ctxt) const
{
  assert(space_ != NULL);
  TransferFunction* tf = space_->GetTransferFunction(type_, TransferFunction::Default(type_, ctxt), true);

  // when we have a transformation we want the physical value for the source design
  DesignElement* trans = space_->ApplyTransformations(this, NULL, NULL);

  return tf->Transform((trans != NULL ? trans : this), SMART, false); // if there is a transformation return the transformed stuff
}


bool DesignElement::HasPhysicalDesign() const
{
  return(type_ == DENSITY || type_ == RHS_DENSITY || type_ == POLARIZATION || (!simp->filter.IsEmpty() && simp->filter[0].GetType() == Filter::DENSITY));
}


void DesignElement::ToInfo(PtrParamNode in, TransferFunction* tf, ErsatzMaterial* em) const
{
  in->Get("type")->SetValue(type.ToString(type_));
  in->Get("upperBound")->SetValue(upper_);
  in->Get("lowerBound")->SetValue(lower_);
  if(tf != NULL)
    in->Get("physicalLowerBound")->SetValue(tf->Transform(lower_));
  if(multimaterial != NULL)
  {
    in->Get("material")->SetValue(multimaterial->name);
    in->Get("mm_index")->SetValue(multimaterial->index);
    multimaterial->ToInfo(in);
  }
}

string DesignElement::ToString(const DesignElement* de, bool barycenter)
{
  if(de == NULL) return "null";

  std::stringstream ss;
  if(de->elem == NULL)
  {
    ss << "ghost";
    if(de->vicinity != NULL) ss << de->vicinity->ToString();
  }
  else
  {
    ss << "e=" << boost::lexical_cast<string>(de->elem->elemNum);
    if(barycenter)
      ss << " bc=" << de->elem->extended->barycenter.ToString();
    else
      ss << " t=" << type.ToString(de->type_);
  }

  return ss.str();
}

string DesignElement::ToString(const StdVector<DesignElement*>& vec, bool print_type)
{
  std::stringstream ss;
  ss << "[";
  for(unsigned int i = 0, s = vec.GetSize(); i < s; ++i)
  {
    ss << ToString(vec[i]);
    if(print_type) ss << "=" << vec[i]->type_;
    if(i < s-1) ss << ",";
  }
  ss << "]";

  return ss.str();
}

string DesignElement::ToString(const StdVector<DesignElement>& vec, bool print_val, bool print_type)
{
  std::stringstream ss;
  ss << "[";
  for(unsigned int i = 0, s = vec.GetSize(); i < s; ++i)
  {
    ss << vec[i].elem->elemNum;
    if(print_val)  ss << ":" << vec[i].GetPlainDesignValue();
    if(print_type) ss << "=" << vec[i].type_;
    if(i < s-1) ss << ",";
  }
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
  Filter::sensitivity.Add(Filter::SIGMUND_TRACE, "sigmund_trace");
  Filter::sensitivity.Add(Filter::SHARP_SIGMUND, "sharp_sigmund");
  Filter::sensitivity.Add(Filter::BORRVALL, "borrvall");

  Filter::density.SetName("Filter::Density");
  Filter::density.Add(Filter::STANDARD, "standard");
  Filter::density.Add(Filter::SOLID_HEAVISIDE, "solid_heaviside");
  Filter::density.Add(Filter::VOID_HEAVISIDE, "void_heaviside");
  Filter::density.Add(Filter::TANH, "tanh");
  Filter::density.Add(Filter::MATERIAL, "material");
  Filter::density.Add(Filter::EXPRESSION, "expression");

  Filter::filterSpace.SetName("Filter::FilterSpace");
  Filter::filterSpace.Add(Filter::RADIUS, "radius");
  Filter::filterSpace.Add(Filter::VOLUME_RADIUS, "volumeRadius");
  Filter::filterSpace.Add(Filter::MAX_EDGE, "maxEdge");

  ShapeMapDesign::type.SetName("ShapeMapDesign::Type");
  ShapeMapDesign::type.Add(ShapeMapDesign::CENTER, "center");
  ShapeMapDesign::type.Add(ShapeMapDesign::NODE, "node");
  ShapeMapDesign::type.Add(ShapeMapDesign::PROFILE, "profile");

  ShapeParamElement::dof.SetName("ShapeParamElement::Dof");
  ShapeParamElement::dof.Add(ShapeParamElement::NOT_SET, "not_set");
  ShapeParamElement::dof.Add(ShapeParamElement::X, "x");
  ShapeParamElement::dof.Add(ShapeParamElement::Y, "y");
  ShapeParamElement::dof.Add(ShapeParamElement::Z, "z");

  FeatureVariable::tip_enum.SetName("FeatureMappingDesign::Tip");
  FeatureVariable::tip_enum.Add(FeatureVariable::START, "start");
  FeatureVariable::tip_enum.Add(FeatureVariable::END, "end"); 
  FeatureVariable::tip_enum.Add(FeatureVariable::INNER, "inner"); // only SpaghettiDesign

  type.SetName("BaseDesignElement::Type");
  type.Add(NO_TYPE, "no_type");
  type.Add(NO_MULTIMATERIAL, "no_multimaterial");
  type.Add(NO_DERIVATIVE, "no_derivative");
  type.Add(SHAPE_MAP, "shape_map");
  type.Add(SPAGHETTI, "spaghetti");
  type.Add(SPLINE_BOX, "spline_box");
  type.Add(MECH_TRACE, "mech_trace");
  type.Add(MECH_ALL, "mech_all");
  type.Add(DIELEC_TRACE, "dielec_trace");
  type.Add(DIELEC_ALL, "dielec_all");
  type.Add(PIEZO_ALL, "piezo_all");
  type.Add(DEFAULT, "default");
  type.Add(DENSITY, "density");
  type.Add(RHS_DENSITY, "rhsDensity");
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
  type.Add(MECH_11, "mech_11");
  type.Add(MECH_12, "mech_12");
  type.Add(MECH_13, "mech_13");
  type.Add(MECH_14, "mech_14");
  type.Add(MECH_15, "mech_15");
  type.Add(MECH_16, "mech_16");
  type.Add(MECH_22, "mech_22");
  type.Add(MECH_23, "mech_23");
  type.Add(MECH_24, "mech_24");
  type.Add(MECH_25, "mech_25");
  type.Add(MECH_26, "mech_26");
  type.Add(MECH_33, "mech_33");
  type.Add(MECH_34, "mech_34");
  type.Add(MECH_35, "mech_35");
  type.Add(MECH_36, "mech_36");
  type.Add(MECH_44, "mech_44");
  type.Add(MECH_45, "mech_45");
  type.Add(MECH_46, "mech_46");
  type.Add(MECH_55, "mech_55");
  type.Add(MECH_56, "mech_56");
  type.Add(MECH_66, "mech_66");
  type.Add(DIELEC_11, "dielec_11");
  type.Add(DIELEC_12, "dielec_12");
  type.Add(DIELEC_22, "dielec_22");
  type.Add(PIEZO_11, "piezo_11");
  type.Add(PIEZO_12, "piezo_12");
  type.Add(PIEZO_13, "piezo_13");
  type.Add(PIEZO_21, "piezo_21");
  type.Add(PIEZO_22, "piezo_22");
  type.Add(PIEZO_23, "piezo_23");
  type.Add(ROTANGLE, "rotAngle");
  type.Add(ROTANGLEFIRST, "rotAngleFirst");
  type.Add(ROTANGLESECOND, "rotAngleSecond");
  type.Add(ROTANGLETHIRD, "rotAngleThird");
  type.Add(STIFF1, "stiff1");
  type.Add(STIFF2, "stiff2");
  type.Add(STIFF3, "stiff3");
  type.Add(SHEAR1, "shear1");
  type.Add(SLACK, "slack");
  type.Add(ALPHA, "alpha");
  type.Add(LOWER_EIG_BOUND, "lowerEigenBound");
  type.Add(MULTIMATERIAL, "multimaterial");
  type.Add(INTERPOLATION, "interpolation");
  type.Add(NODE, "node");
  type.Add(PROFILE, "profile");
  type.Add(NORMAL, "normal");
  type.Add(RADIUS, "radius");
  type.Add(FEATURE_MAPPING_PX, "feature_var_Px");
  type.Add(FEATURE_MAPPING_PY, "feature_var_Py");
  type.Add(FEATURE_MAPPING_QX, "feature_var_Qx");
  type.Add(FEATURE_MAPPING_QY, "feature_var_Qy");
  type.Add(FEATURE_MAPPING_P,  "feature_var_P");        
  type.Add(CP, "controlpoint");
  type.Add(ALL_FEATURES, "allFeatures");   // e.g. featureDistance -> min of all distances
  type.Add(FEATURE, "feature");            // 0-based index specified in "generic" as number
  type.Add(ALL_DESIGNS, "allDesigns");     // ALL_DESIGNS needs last enum value   

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
  valueSpecifier.Add(QUADRATIC_VM_STRESS, "quadraticVMStress");
  valueSpecifier.Add(LOCAL_LOAD_FACTOR, "localBucklingLoadFactor");
  valueSpecifier.Add(DESIGN_TRACKING, "designTracking");
  valueSpecifier.Add(WEIGHT, "weight");
  valueSpecifier.Add(OBJECTIVE, "objective");
  valueSpecifier.Add(PROJECTION, "projection");
  valueSpecifier.Add(NUM_NEIGHBOURS, "neighbours");
  valueSpecifier.Add(LEVEL_SET_VALUE, "levelSetValue");
  valueSpecifier.Add(LEVEL_SET_STATE, "levelSetState");
  valueSpecifier.Add(TOPGRAD_VALUE, "topGradValue");
  valueSpecifier.Add(SHAPEGRAD_VALUE, "shapeGradValue");
  valueSpecifier.Add(SHAPEGRAD_NODE_VALUE, "shapeGradNodeValue");
  valueSpecifier.Add(SHAPE_MAP_GRAD, "shapeMapGrad");
  valueSpecifier.Add(SHAPE_MAP_ORDER, "shapeMapIntOrder");
  valueSpecifier.Add(SHAPE_MAP_CORNER, "shapeMapMinMaxCorner");
  valueSpecifier.Add(SPLINE_BOX_GRAD_X, "splineBoxGradX");
  valueSpecifier.Add(SPLINE_BOX_GRAD_Y, "splineBoxGradY");
  valueSpecifier.Add(SPLINE_BOX_GRAD_Z, "splineBoxGradZ");
  valueSpecifier.Add(SPLINE_BOX_INT_ORDER, "splineBoxIntOrder");
  valueSpecifier.Add(SPLINE_BOX_INT_CORNER, "splineBoxMinMaxCorner");
  valueSpecifier.Add(LEVEL_SET_GRAD_XP, "levelSetGradXP");
  valueSpecifier.Add(LEVEL_SET_GRAD_XN, "levelSetGradXN");
  valueSpecifier.Add(LEVEL_SET_GRAD_YP, "levelSetGradYP");
  valueSpecifier.Add(LEVEL_SET_GRAD_YN, "levelSetGradYN");
  valueSpecifier.Add(LEVEL_SET_GRAD_ZP, "levelSetGradZP");
  valueSpecifier.Add(LEVEL_SET_GRAD_ZN, "levelSetGradZN");
  valueSpecifier.Add(HEAT_NODAL_TRACK_VAL, "heatNodalTrackValue");
  valueSpecifier.Add(TEMP_AT_INTERFACE, "tempAtInterface");
  valueSpecifier.Add(MMA_ASYMPTOTE, "mmaAsymptote");
  valueSpecifier.Add(MMA_LOWER_VAL, "mmaLowerVal");
  valueSpecifier.Add(MMA_UPPER_VAL, "mmaUpperVal");
  valueSpecifier.Add(MMA_OBJ_GRADIANT, "mmaGradiant_0");
  valueSpecifier.Add(MMA_CON_GRADIANT_1, "mmaGradiant_1");
  valueSpecifier.Add(MMA_CON_GRADIANT_2, "mmaGradiant_2");
  valueSpecifier.Add(GENERIC_ELEM, "genericElem");
  valueSpecifier.Add(FILTERED_DESIGN, "filteredDesign");
  valueSpecifier.Add(DIFF_FILTERED_DESIGN, "diffFilteredDesign");
  valueSpecifier.Add(FEATURE_DISTANCE, "featureDistance");
  valueSpecifier.Add(FEATURE_PROJECTED, "featureProjected");
  valueSpecifier.Add(FEATURE_GRAD, "featureGrad");

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
  // this is a selection of constraints for constraintGradient
  detail.Add(DYNAMIC_OUTPUT, "dynamicOutput");
  detail.Add(REFLECTED_WAVE, "reflectedWave");
  detail.Add(COMPLIANCE, "compliance");
  detail.Add(VOLUME, "volume");
  detail.Add(PENALIZED_VOLUME, "penalizedVolume");
  detail.Add(GAP, "gap");
  detail.Add(REALVOLUME, "realvolume");
  detail.Add(TRACKING, "tracking");
  detail.Add(HOMOGENIZATION_TRACKING, "homTracking");
  detail.Add(POISSONS_RATIO, "poissonsRatio");
  detail.Add(YOUNGS_MODULUS, "youngsModulus");
  detail.Add(YOUNGS_MODULUS_E1, "youngsModulusE1");
  detail.Add(YOUNGS_MODULUS_E2, "youngsModulusE2");
  detail.Add(TYCHONOFF, "tychonoff");
  detail.Add(GREYNESS, "greyness");
  detail.Add(GLOBAL_SLOPE, "globalSlope");
  detail.Add(GLOBAL_CHECKERBOARD, "globalCheckerboard");
  detail.Add(GLOBAL_DESIGN, "globalDesign");
  detail.Add(STRESS, "stress");
  
  // secific data - e.g. for debuging   
  detail.Add(PROJECTION_FILTER, "projectionFilter");
  detail.Add(SM_NODE, "node");
  detail.Add(SM_NODE_A, "node_a");
  detail.Add(SM_NODE_B, "node_b");
  detail.Add(SM_PROFILE, "profile");
  detail.Add(SP_CP, "controlpoint");
  detail.Add(BUCKLING_LOAD_FACTOR, "bucklingLoadFactor");
  detail.Add(LOCAL_BUCKLING_LOAD_FACTOR, "localBucklingLoadFactor");
  detail.Add(GLOBAL_BUCKLING_LOAD_FACTOR, "globalBucklingLoadFactor");
  detail.Add(SIN, "sin");
  detail.Add(COS, "cos");
  detail.Add(GRAD_DISTANCE, "gradDistance"); // nodal derivative of featureDistance/Projection by 'design' feature_var_Px, ...
  detail.Add(FEATURE_PART, "featurePart");   // nodal 'featureDistance' feature part idx 0,1,2 for given feature
  detail.Add(FEATURE_RHO, "featureRho");     // for element 'featureGrad'
  detail.Add(COMBINE, "combine");            // for element 'featureGrad'
}


SIMPElement::SIMPElement(DesignElement* base)
{
  this->de_ = base;
}

double SIMPElement::GetSensitivityFilteredValue(DesignElement::ValueSpecifier sp, Function* g) const
{
  // We filter over this element and the neighbors.
  assert(de_->simp != NULL);
  assert(this == de_->simp);
  assert(de_->simp->filter.GetSize() == 1); // no robustness for sensitivity filtering worth implementing

  const Filter& f = filter[DetermineFilterIndex()];
  const GlobalFilter* gf = f.global;

  assert(f.GetType() == Filter::SENSITIVITY);

  // See Filter::Sensitivity: w = weight, p is density, f' is cost gradient
  // Sigmund  = sum_i w(x_i) * p_i * f_i' / p_e * sum_i w(x_i)
  // Sharp Sigmund  = sum_i (i=e ? 1:0 : w(x_i)) * p_i * f_i' / p_e * sum_i w(x_i) + "bug" in normalized weighting
  // Borrvall = sum_i w(x_i) * p_i * f_i' / sum_i p_i * w(x_i)
  // Safe Borrvall if |p_i| = epsilon for all i use plain filter, otherwise  Borrvall filter, usefull for filtering
  // of variables which are not bounded away from zero

  // plain    = sum_i w(x_i) * f_i' / sum_i w(x_i)
  // sharp plain = plain, just the filter is setup with normalized weights with a "bug"

  // factor design in numerator (SIGMUND and BORRVALL) for densitiy filtering the value is any in
  bool numerator_design = gf->sensitivity != Filter::PLAIN && gf->sensitivity != Filter::SHARP_PLAIN;
  // factor design in denominator sum (BORRVALL)
  bool denominator_design = (gf->sensitivity == Filter::BORRVALL);
  // weight the denominator by this design (SIGMUND)
  bool sigmund_denominator = gf->sensitivity == Filter::SIGMUND || gf->sensitivity == Filter::SHARP_SIGMUND;
  // short-cut for fake in Sharp Sigmund: i=e ? 1:0 : w(x_i)
  bool cheat_this_weight = gf->sensitivity == Filter::SHARP_SIGMUND || gf->sensitivity == Filter::SHARP_PLAIN;
  // weight the denominator by trace (SIGMUND_TRACE)
  bool sigmund_trace = gf->sensitivity == Filter::SIGMUND_TRACE;

  LOG_DBG3(desel) << "GSFV: el=" << de_->elem->elemNum
                << " t=" << de_->ToString()
                << " sp=" << DesignElement::valueSpecifier.ToString(sp)
                << " dens=" << Filter::density.ToString(gf->density)
                << " numerator_design=" << numerator_design
                << " denominator_design=" << denominator_design
                << " sigmund_denominator=" << sigmund_denominator
                << " cheat_this_weight=" << cheat_this_weight
                << " sigmund_trace=" << sigmund_trace;


  double numerator = 0.0;
  double denominator = 0.0;

  // mathematically the neighborhood includes this element, but this is not in the structure
  for(int i = -1, ni = (int) f.neighborhood.GetSize(); i < ni; i++)
  {
    const Filter::NeighbourElement* ne = i == -1 ? NULL : &f.neighborhood[i];
    const DesignElement* de = i == -1 ? this->de_ : ne->neighbour;

    double w = i == -1 ?  f.weight : ne->weight;
    double nw = cheat_this_weight && i == -1 ? 1.0 : w;
    double v = de->GetPlainValue(sp, dynamic_cast<Condition*>(g));
    double x = sigmund_trace ?  de->GetPlainMechTrace() : de->GetPlainValue(DesignElement::DESIGN); // cheap if not used

    double numerator_summand   = nw * v * (numerator_design ? x : 1.0);
    double denominator_summand = w * (denominator_design ? x : 1.0);

    numerator   += numerator_summand;
    denominator += denominator_summand;

    LOG_DBG3(desel) << "GSFV: el=" << de_->elem->elemNum << ": curr=" << de->elem->elemNum
                  << " w= " << w  << " nw=" << nw << " x=" << x << " v=" << v << " ns=" << numerator_summand
                  << " ds=" << denominator_summand << " num=" << numerator << " den=" << denominator;
  }

  if(sigmund_denominator)
    denominator *= de_->GetDesign(DesignElement::PLAIN);
  if(sigmund_trace)
    denominator *= de_->GetPlainMechTrace();

  LOG_DBG3(desel) << "GSFV: el=" << de_->elem->elemNum << ": result=" << numerator << "/" << denominator << " -> " << (numerator / denominator);

  return (fabs(denominator) < std::numeric_limits<float>::epsilon()) ? de_->GetPlainValue(sp, dynamic_cast<Condition*>(g)) : (numerator / denominator);
}

double SIMPElement::CalcNonLinFilter(double filtered_value, const GlobalFilter* global, double plain_value)
{
  switch(global->density)
  {
  case Filter::SOLID_HEAVISIDE:
  case Filter::VOID_HEAVISIDE:
    return CalcHeaviside(filtered_value, global);
  case Filter::TANH:
    return CalcTanh(filtered_value, global);
  case Filter::MATERIAL:
    return CalcMaterial(filtered_value, plain_value, global);
  case Filter::EXPRESSION:
    return global->EvalExpressionFunction(filtered_value, true);
  default:
    assert(false);
    return -1.0;
  }
}


double SIMPElement::GetDensityFilteredValue(DesignElement::ValueSpecifier sp, Filter::Density fd) const
{
  // We filter over this element and the neighbors.
  assert(sp == DesignElement::DESIGN);
  assert(!de_->simp->filter.IsEmpty());
  assert(de_->simp != NULL);

  unsigned int fix = DetermineFilterIndex();
  const Filter& f = DetermineFilter();
  const GlobalFilter* gf=f.global;
  DesignSpace * space = de_->GetDesignSpace();

  // All equations from Sigmund; Morphology based black and white filters for topology optimization; 2007
  // p = rho. P is filtered rho (rho tilde)
  // P = sum_(i in N_e) w(x_i) p_i / sum_(i in N_e) w(x_i)
  // mathematically the neighborhood includes this element, but this is not in the structure
  // we initialize numerator and denominator with the values obtained from this element

  double p_filt = 0.0;
  if (space->is_matrix_filt){
    // (idx within DesignSpace::data, dt, FE ElemNum):
    // we assume that in DesignSpace::data, an equally large number of design variables exists for each defined design type
    // here: 3 x density, 3 x angle
    // (0, density, 1),(1, density, 2), (2, density, 3), (3, angle, 1), (4, angle, 2), (5, angle, 3)
    LOG_DBG3(desel) << "de idx=" << de_->GetIndex() << " space->FindDesign(de_->GetType()):" << space->FindDesign(de_->GetType()) << "  space->GetNumberOfElements():" << space->GetNumberOfElements();
    unsigned int elem_idx = de_->GetIndex() - space->FindDesign(de_->GetType()) * space->GetNumberOfElements();
    p_filt =  space->density_filter[fix].filtered_vec[elem_idx];
    LOG_DBG3(desel) << "elemNum=" << de_->elem->elemNum << " filtered value=" << p_filt;
  }
  else
  {
    double rho_e = this->de_->GetPlainValue(DesignElement::DESIGN);
    if(fd == Filter::MATERIAL || fd == Filter::MATERIAL_PART)
      rho_e = gf->mat_filter.Transform(rho_e);
    double numerator = f.weight * rho_e;
    LOG_DBG3(desel) << "GDFV: el=" << de_->elem->elemNum << ": curr=" << de_->elem->elemNum
                      << " w= " << f.weight << " rho_e=" << this->de_->GetPlainValue(DesignElement::DESIGN)
                      << " num=" << numerator << " fix=" << fix;

    for(int i = 0, ni = (int) f.neighborhood.GetSize(); i < ni; i++)
    {
      const Filter::NeighbourElement* ne = &f.neighborhood[i];
      const DesignElement* de = ne->neighbour;

      double w = ne->weight;
      double rho_i = de->GetPlainDesignValue();
      if(fd == Filter::MATERIAL || fd == Filter::MATERIAL_PART)
        rho_i = gf->mat_filter.Transform(rho_i);
      numerator   += w * rho_i;
      // LOG_DBG3(desel) << "GDFV: el=" << de_->elem->elemNum << ": curr=" << de->elem->elemNum  << " w= " << w  << " rho_i=" << rho_i
      //                 << " num=" << numerator << " -> " << (numerator / f.weight_sum);
    }
    p_filt = numerator / f.weight_sum; // includes the element de
  }

  // all the Calc*() apply scaling scaling to correct nonlinear filters
  switch(fd)
  {
  case Filter::TANH:
    p_filt = CalcTanh(p_filt, gf);
    break;
  case Filter::SOLID_HEAVISIDE:
  case Filter::VOID_HEAVISIDE:
    p_filt = CalcHeaviside(p_filt, gf);
    break;
  case Filter::EXPRESSION:
    p_filt = gf->EvalExpressionFunction(p_filt, true);
    break;
  case Filter::MATERIAL:
  {
    double plain = de_->GetPlainDesignValue();
    LOG_DBG3(desel) << "GDFV: el=" << de_->elem->elemNum << " MATERIAL ms(" << p_filt << ")=" << gf->mat_scale.Transform(p_filt) << " *mp(" << plain<< ")=" << gf->mat_phase.Transform(plain)
                    << " offset=" << gf->non_lin_offset << " scale=" << gf->non_lin_scale;
    p_filt = CalcMaterial(p_filt, plain, gf);
    break;
  }
  case Filter::STANDARD:
  case Filter::MATERIAL_PART: // special case where we are almost like standard
    break; // nothing to do
  default:
    assert(false);
  }

  LOG_DBG3(desel) << "GDFV: el=" << de_->elem->elemNum << " filter_matrix=" << space->is_matrix_filt
                  << " design=" << Filter::density.ToString(gf->density)
                  << ": plain " << this->de_->GetPlainValue(DesignElement::DESIGN) << " offset=" << gf->non_lin_offset << " scale=" << gf->non_lin_scale << " -> "<< p_filt;
  assert(p_filt <= 1.000001 * this->de_->GetUpperBound());
  if(fd != Filter::STANDARD && fd != Filter::MATERIAL_PART)
  {
    if(p_filt < 0.7 * this->de_->GetLowerBound()) {
      LOG_DBG(desel) << "GDFV: el=" << de_->elem->elemNum << " design=" << Filter::density.ToString(gf->density)
                      << ": plain " << this->de_->GetPlainValue(DesignElement::DESIGN) << " p_filt=" << p_filt << " lb=" << this->de_->GetLowerBound() << " ub=" << this->de_->GetUpperBound();
    }
    assert(p_filt <= 1.000001 * this->de_->GetUpperBound());
    // standard is save and we might often have the case with too small loaded initial designs and enforce_bounds not set
    assert(p_filt >= 0.7 * this->de_->GetLowerBound()); // in the < 07.2021 times there was a check of filter[fix].GetLowerBound(this->de_))
  }
  return p_filt;

}


double SIMPElement::GetDensityFilteredGradient(DesignElement::ValueSpecifier sp, Function* func) const
{
  // We filter over this element and the neighbors.
  assert(de_->simp != NULL);
  assert(this == de_->simp);

  unsigned int fix = DetermineFilterIndex();
  const Filter& f = filter[fix];
  const GlobalFilter* gf = f.global;

  Condition* g = dynamic_cast<Condition*>(func);

  assert(f.GetType() == Filter::DENSITY);
  assert(sp == DesignElement::COST_GRADIENT || sp == DesignElement::CONSTRAINT_GRADIENT);
  assert((func == NULL || (func->IsObjective() && sp == DesignElement::COST_GRADIENT)) || (func == NULL || (!func->IsObjective() && sp == DesignElement::CONSTRAINT_GRADIENT)) || (func == NULL || (func->IsObjective())));
  // projection has density filtering only in the fake filter problem but not in the original problem (which should not be density filtered anyway)
  assert(func == NULL || func->isFiltered());


  // Density filtering for gradient is (Sigmund; Morphology-based black and white filters for topology optimization; 2007; eqn (35). (36)
  // p is rho and P is rho filtered!
  // d f/d p_e = sum_i(in N_e) d f/d P_i * d P_i/d p_e with d P_i/d p_e = w(x_e)/ sum_j(in N_i) w(x_j)
  // note, that the stored value is already v = d f/d P_i

  // LOG_DBG3(desel) << "GDFG: el=" << de_->elem->elemNum << " sp=" << DesignElement::valueSpecifier.ToString(sp) << " g=" << (g != NULL ? Condition::type.ToString(g->GetType()) : "null");

  double sum = 0.0;

  // used there in every loop
  double rho_e = gf->density == Filter::MATERIAL ? this->de_->GetPlainDesignValue() : -1.0; // plain rho_e

  for(int i = -1, ni = (int) f.neighborhood.GetSize(); i < ni; i++)
  {
    const Filter::NeighbourElement* ne = i == -1 ? NULL : &f.neighborhood[i];
    const DesignElement* de = i == -1 ? this->de_ : ne->neighbour;
    // note, that de might be another region than this!

    // d f/d P_i
    double df = de->GetPlainValue(sp, g);

    double w = i == -1 ? f.weight : ne->weight;
    assert(de->simp->filter[fix].weight_sum >= 0.0); // set in DesignStructure

    // d filtered_density_i/ d_rho_e: using not the weight sum of de_i but of de_e fails in gradient check and math, see Ole's morphology paper
    double dP = w / de->simp->filter[fix].weight_sum;

    // all non-standard filters have a scaling to match rho_min
    // when putting scale out off loop, make sure not to change it
    double nlf_scale = gf->non_lin_scale;  // for not-standard filters this a factor for the derivative

    // LOG_DBG3(desel) << "GDFG: el=" << de_->elem->elemNum << ": curr=" << de->elem->elemNum << " i=" << i
    //                      << " df=" << df << " w=" << w << " Ew=" << de->simp->filter[fix].weight_sum << " dP=" << dP
    //                      << " nlf_scale=" << nlf_scale;

    double summand = 0.0;

    switch(gf->density)
    {
    case Filter::STANDARD:
      assert(nlf_scale == 1.0);
      assert(gf->non_lin_offset == 0.0);
      summand = df * dP;
      break;
    case Filter::MATERIAL_PART:
      assert(false); // shall not be called
      break;
    case Filter::MATERIAL:
    {
      // the material filter is in principal rho_e * filtered_rho_e,
      // the derivative is product rule instead of chain rule
      // we have a slight complication by transfer functions mphase,mscale,mfilter
      //
      // F_t(rho)_e = mphase(rho_e)*mscale(P(mfilter(rho)))

      // F_t'= mphase'(rho_e) * mscale(P(mfilter(rho)))
      //     + mphase(rho_e)*mscale'*F(mfilter(rho))*F'(mfilter(rho))*mfilter'(rho)

      double mphase = gf->mat_phase.Transform(rho_e);
      double dmphase = i == -1 ? gf->mat_phase.Derivative(rho_e) : 0.0;

      // our filtered design is subject mat_filter, therefore we cannot use P_i and calling with MATERIAL would use phase and scale
      double mfilt_i = de->simp->GetDensityFilteredValue(DesignElement::DESIGN, Filter::MATERIAL_PART);
      double mscale =  i == -1 ? gf->mat_scale.Transform(mfilt_i) : 0.0; // not of interest when i != e
      double dmscale = gf->mat_scale.Derivative(mfilt_i);

      double rho_i = de->GetPlainDesignValue();
      double dmfilter = gf->mat_filter.Derivative(rho_i);

      // df = sum_i^N_e df/drho_i * (d_rho_e/d_rho_i  * rho_filt_e + rho_e * dP)
      summand = df * nlf_scale * (dmphase * mscale + mphase * dmscale * dP * dmfilter);

      // LOG_DBG3(desel) << "GDFG: el=" << de_->elem->elemNum << ": curr=" << de->elem->elemNum << " i=" << i << " MATERIAL "
      //                << " rho_e=" << rho_e << " rho_i=" << rho_i << " dmphase=" << dmphase << " mscale=" << mscale
      //                << " mphase=" << mphase << " dmscale=" << dmscale << " dmfilter=" << dmfilter << " -> " << summand;
      break;
    }
    case Filter::SOLID_HEAVISIDE:
    {
      // F = 1.0 - std::exp(-b * input_value) + input_value * std::exp(-b)
      // H = scale * F + offset
      double P_i = de->simp->GetDensityFilteredValue(DesignElement::DESIGN, Filter::STANDARD);
      double b = gf->beta;
      nlf_scale *= (b * exp(-b * P_i) + exp(-b));
      summand = df * nlf_scale* dP;
      break;
    }
    case Filter::VOID_HEAVISIDE:
    {
      // general scaling
      double P_i = de->simp->GetDensityFilteredValue(DesignElement::DESIGN, Filter::STANDARD);
      double b = gf->beta;
      nlf_scale *= (b * exp(b*(P_i-1.0)) + exp(-1.0*b));
      summand = df * nlf_scale* dP;
      break;
    }
    case Filter::TANH:
    {
      // f(x)  =  1 - 1/(exp(2*beta*(x-param)) + 1)
      // f'(x) =  (exp(2*beta*(x-param)+1)^-2 * 2 * beta * exp(2*beta*(x-param))
      double P_i = de->simp->GetDensityFilteredValue(DesignElement::DESIGN, Filter::STANDARD);
      double b = gf->beta;
      double eta = gf->eta;
      double e = exp(2.0 * b * (P_i - eta));
      nlf_scale *= (1.0/((e+1.0)*(e+1.0)) * 2.0 * b * e);
      summand = df * nlf_scale* dP;
      break;
    }
    case Filter::EXPRESSION:
    {
      double P_i = de->simp->GetDensityFilteredValue(DesignElement::DESIGN, Filter::STANDARD);
      nlf_scale *= gf->EvalExpressionDerivative(P_i, false);
      summand = df * nlf_scale * dP;
    }
    } // end of switch

    sum += summand;
    // LOG_DBG3(desel) << "GDFG: el=" << de_->elem->elemNum << ": curr=" << de->elem->elemNum << " t=" << f.type.ToString(f.GetType())
    //                << " f=" << (func ? func->ToString() : "-") << " df=" << df
    //                << " w=" << w << " Ew=" << de->simp->filter[fix].weight_sum << " dP=" << dP
    //                << " summand=" << summand << " sum=" << sum;

  } // end of i-loop

  return sum;
}


string SIMPElement::ToString(int level) const
{
  std::stringstream ss;
  const Filter& f = filter[DetermineFilterIndex()];
  ss << "el=" << de_->elem->elemNum << " #n=" << f.neighborhood.GetSize() << "(";
  for(unsigned int i = 0; !filter.IsEmpty() && i < f.neighborhood.GetSize(); i++)
  {
    if(level == 0)
      ss << " " << f.neighborhood[i].neighbour->elem->elemNum;
    else {
      ss << " n_" << f.neighborhood[i].neighbour->elem->elemNum << "_w=" << f.neighborhood[i].weight;
      ss << " n_" << f.neighborhood[i].neighbour->elem->elemNum << "_d=" << f.neighborhood[i].distance;
    }
  }
  ss << ")";
  return ss.str();
}

inline unsigned int SIMPElement::DetermineFilterIndex() const
{
  return Optimization::context->GetExcitation()->robust_filter_idx;
}

unsigned int SIMPElement::DetermineFilterIndexNonInlined() const
{
  return Optimization::context->GetExcitation()->robust_filter_idx;
}

inline const Filter& SIMPElement::DetermineFilter() const
{
  return filter[DetermineFilterIndex()];
}

void SIMPElement::Dump()
{
  std::cout << "\nelement: " << de_->elem->elemNum << " location " << de_->GetLocation()->ToString() << std::endl;
  for(unsigned int f = 0; f < filter.GetSize(); f++)
    filter[f].Dump();

}


VicinityElement::VicinityElement()
{
  design.Resize(domain->GetDim() == 2 ? 4 : 6);
  design.Init(NULL);
  periodic = false;
}


void VicinityElement::Init(DesignSpace* space, DesignStructure* structure)
{
  // do it only once
  if(space->data[0].vicinity != NULL)
    return;

  Grid* grid = domain->GetGrid();

  // we only hope the elements are aligned and rectangular, the size might vary
  // if(!space->IsRegular())
  //  throw Exception("A regular design domain is required to use VicinityElements");

  // eventually the barycenters are already calculated, we need them to identify the orientation
  // we will need the barycenters in FindNeibhborhood()
  for(unsigned int i = 0, s = space->regions[0].GetSize(); i < s; i++)
    grid->SetElementBarycenters(space->regions[0][i].regionId, false); // no updated coordinates
  // Handle also the off-design barycenters
  if(space->DoNonDesignVicinity())
    for(unsigned int i = 0; i < space->GetPseudoDesignRegions().GetSize(); i++)
      grid->SetElementBarycenters(space->GetPseudoDesignRegions()[i][0].elem->regionId, false);

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
  Elem* elem = space->data[0].elem;
  domain->GetGrid()->GetElemShapeMap(elem, false)->GetEdgeLength(spacing);

  // in the periodic case the neighborhood is enlarged
  StdVector<std::pair<Elem*, int> > enlarged_data;

  // traverse over all elements
  // double cfs elements for multiple design variables are not handled special
  for(unsigned int e = 0, n = space->GetTotalElements().GetSize(); e < n; e++)
  {
    DesignElement* de = space->GetTotalElements()[e];
    de->vicinity = new VicinityElement();
    // here we store the neighbors in a sorted way
    StdVector<DesignElement*>& ve_data = de->vicinity->design; // has proper size of NULLs

    StdVector<std::pair<Elem*, int> >& neighbors = *(de->elem->extended->neighborhood);
    // reuse the enlarged_data element for the periodic case only
    if(periodic)
    {
      if(structure->ExtendPeriodicNeighborhood(de->elem, common, enlarged_data))
      {
        neighbors = enlarged_data; // in the non-periodic case there is no need to copy the element data
        de->vicinity->periodic = true;
      }
    }

    for(unsigned int n = 0; n < neighbors.GetSize(); n++)
    {
      // we consider only direct (edge/face) neighbors
      if(neighbors[n].second < common) continue;
      // now we have to find the relative position of candidate
      Elem* candidate = neighbors[n].first;


      LOG_DBG3(desel) << "VE:Init elem=" << de->elem->elemNum << " e.bc=" << de->elem->extended->barycenter.ToString() << " e.dim="
                    << de->elem->GetShape().dim << " e.r=" << de->elem->regionId << " n=" << n << " o.el=" << candidate->elemNum
                    << " o.bc=" << candidate->extended->barycenter.ToString() << " o.dim=" << candidate->GetShape().dim << " o.r=" << candidate->regionId
                    << " e.c=" << de->elem->connect.ToString() << " o.c=" << candidate->connect.ToString();


      // if the neighbor is a surface element we don't want to play with it
      if(de->elem->GetShape().dim != candidate->GetShape().dim)
        continue;
      // we include pseudo design regions. This is for off-design optimization or non_design_vicinity=true in <ersatzMaterial>
      if(!space->Contains(candidate->regionId, space->DoNonDesignVicinity()))
        continue;

      // the spacing allows to identify periodic elements
      int idx = FindRelativeNeighborLocation(de->elem->extended->barycenter, candidate->extended->barycenter, spacing);
      ve_data[idx] = space->Find(candidate->elemNum, de->GetType(), false, space->DoNonDesignVicinity());
      LOG_DBG2(desel) << "VE:Init elem=" << de->elem->elemNum << " idx=" << idx << " val=" << ve_data[idx]->ToString() << " neigh=" << DesignElement::ToString(ve_data);
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

  LOG_DBG2(desel) << "VE:FRNL ref =" << ref.ToString() << " other=" << other.ToString() << " -> " << res;

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
      LOG_DBG3(desel) << "VE:GN base=" << base->ToString() << " idx=" << idx << " n=" << n << " max neighbor=" << i;
      if(throw_exception)
        EXCEPTION("no neighbor in " << idx << " direction " << i << " elements ways for element " << tmp->ToString())
      else
        return NULL;
    }
  }
  LOG_DBG3(desel) << "VE:GN base=" << base->ToString() << " idx=" << idx << " n=" << n << " -> " << (tmp != NULL ? tmp->ToString() : "null");
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
  std::stringstream ss;
  ss << "(";
  unsigned int max = domain->GetDim() == 2 ? 4 : 6;
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
  ss << " periodic=" << periodic;
  ss << ")";
  return ss.str();
}


ResultDescription::ResultDescription(PtrParamNode pn)
{
  solutionType = SolutionTypeEnum.Parse(pn->Get("id")->As<string>());

  design = pn->Has("design") ? DesignElement::type.Parse(pn->Get("design")->As<string>()) : DesignElement::DEFAULT;

  access = DesignElement::access.Parse(pn->Get("access")->As<string>());

  value = DesignElement::valueSpecifier.Parse(pn->Get("value")->As<string>());

  detail = DesignElement::detail.Parse(pn->Get("detail")->As<string>());

  excitation = pn->Has("excitation") ? pn->Get("excitation")->As<int>() : -1;

  if(pn->Has("generic"))
    generic = pn->Get("generic")->As<string>(); // otherwise default

  if(value == DesignElement::GENERIC_ELEM && generic.size() == 0)
    throw Exception("a result 'generic' needs the 'generic' attribute set");

  LOG_DBG(desel) << "RD:RD " << ToString();
}


string ResultDescription::ToString()
{
  std::stringstream ss;
  ss << "RD: design=" << DesignElement::type.ToString(design)
     << " access=" << DesignElement::access.ToString(access)
     << " value=" << DesignElement::valueSpecifier.ToString(value)
     //<< " detail=" << DesignElement::detail.ToString(detail)
     << " generic=" << generic
     << " ex=" << excitation;
  return ss.str();
}
