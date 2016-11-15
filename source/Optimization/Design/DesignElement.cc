#include <cmath>
#include <iostream>
#include <utility>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/Logging/log.hpp"
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
#include "Optimization/Design/DesignStructure.hh"
#include "Optimization/Function.hh"
#include "Optimization/LevelSet.hh"
#include "Optimization/Objective.hh"
#include "Optimization/Context.hh"
#include "Optimization/TransferFunction.hh"
#include "PDE/SinglePDE.hh"
#include "Utils/Point.hh"
#include "Utils/tools.hh"
#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/lexical_cast.hpp"

using namespace std;
using namespace CoupledField;

using boost::posix_time::ptime;
using boost::posix_time::second_clock;
using boost::posix_time::microsec_clock;

DECLARE_LOG(desel)
DEFINE_LOG(desel, "designElement")

// the static enum
Enum<Filter::Type>                  Filter::type;
Enum<Filter::Density>               Filter::density;
Enum<Filter::Sensitivity>           Filter::sensitivity;
Enum<BaseDesignElement::Type>       BaseDesignElement::type;
Enum<DesignElement::ValueSpecifier> DesignElement::valueSpecifier;
Enum<DesignElement::Access>         DesignElement::access;
Enum<DesignElement::Detail>         DesignElement::detail;
Enum<ShapeMapDesign::Type>          ShapeMapDesign::type;
Enum<ShapeMapDesign::Symmetry>      ShapeMapDesign::symmetry;

// is a static attribute
DesignSpace* DesignElement::space_(NULL);

std::string BaseDesignElement::ToString() const
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
    // Tensor trace for param mat
    case STIFF1:
    case STIFF2:
    case STIFF3:
    case SHEAR1:
    //for mod_red
    case SCALING1:
    case SCALING2:
    case ROTANGLE:
    case ROTANGLE2:
    case G11:
    case G12:
    case G21:
    case G22:
    case G_MAP_X:
    case G_MAP_Y:
    case GX_0:
    case GY_0:
    case GX_PX:
    case GY_PX:
    case GX_PY:
    case GY_PY:
    case GX_PXY:
    case GY_PXY:
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
    case MECH_22:
    case MECH_23:
    case MECH_33:
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

  case G_ALL:
  {
    switch(test)
    {
    case G11:
    case G12:
    case G21:
    case G22:
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

BaseDesignElement::BaseDesignElement(Type t) {
  design          = 0.0;
  upper_          = 0.0;
  lower_          = 0.0;
  type_           = t;
  index_          = numeric_limits<unsigned int>::max();
}


/** Get the gradient values for either objective or constraint */
double BaseDesignElement::GetPlainGradient(const Function* f) const
{
  assert(f != NULL); // Call SumObjectiveGradient() if you don't want to specify!

  assert(!f->IsObjective() || (f->IsObjective() && dynamic_cast<const Objective*>(f) != NULL));
  assert( f->IsObjective() || (!f->IsObjective() && dynamic_cast<const Condition*>(f) != NULL));

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
                << " penalty=" << (f != NULL ? boost::lexical_cast<std::string>(f->GetPenalty()) : "-")
                << " old= " <<  (f != NULL ? costGradient[f->GetIndex()] : constraintGradient[g->GetIndex()])
                << " add=" << (f != NULL ? value * f->GetPenalty() : value)
                << " -> " << (f != NULL ? costGradient[f->GetIndex()] + value * f->GetPenalty() : constraintGradient[g->GetIndex()] + value);
  if(f != NULL) costGradient[f->GetIndex()] += value * f->GetPenalty();
           else constraintGradient[g->GetIndex()] += value;
}

void BaseDesignElement::AddGradient(const Function* f, double value)
{
  assert(( f->IsObjective() && dynamic_cast<const Objective*>(f) != NULL)
      || (!f->IsObjective() && dynamic_cast<const Condition*>(f) != NULL) );

  if(f->IsObjective())
    costGradient[f->GetIndex()] += value * static_cast<const Objective*>(f)->GetPenalty();
  else
    constraintGradient[f->GetIndex()] += value;
}

void BaseDesignElement::Reset(ValueSpecifier vs, Function*  f)
{
  switch(vs)
  {
  case COST_GRADIENT:
    for(unsigned int i = 0, s = costGradient.GetSize(); i < s; ++i)
      costGradient[i] = 0.0;
    break;
  case CONSTRAINT_GRADIENT:
    if(f != NULL)
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

std::string BaseDesignElement::ToString(const StdVector<BaseDesignElement*>& vec, bool print_type)
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
  opt_index_ = std::numeric_limits<unsigned int>::max();
  dof = -1;
  coord.Resize(domain->GetGrid()->GetDim(), -1.0);
  idx.Resize(domain->GetGrid()->GetDim(), -1);
}

std::string ShapeParamElement::ToString() const
{
  std::stringstream ss;
  ss << " idx=" << index_ << " opt_idx=" << opt_index_ << " t=" << type.ToString(type_);
  return ss.str();
}


/** The default constructor for StdVector and ghost elements*/
DesignElement::DesignElement() : BaseDesignElement()
{
  Init();
  this->interfaceDrivenLoadGrad_.Resize(4 * (domain->GetGrid()->GetDim()-1),0.0);
}

/** The default constructor for StdVector and ghost elements*/
DesignElement::DesignElement(Elem* elem, Type type, unsigned int index, int pseudoElementIndex) : BaseDesignElement()
{
  Init();
  this->elem = elem;
  this->type_ = type;
  this->index_ = index;
  this->pseudoElementIndex_ = pseudoElementIndex;
  this->upper_ = 1.0;
  this->lower_ = 1.0;
  this->multimaterial = NULL;
  this->specialResult.Resize(9, 0.0);
  this->interfaceDrivenLoadGrad_.Resize(4 * (domain->GetGrid()->GetDim()-1),0.0);
}


DesignElement::DesignElement(Type dt, double lower, double upper, Elem* elem, unsigned int index, MultiMaterial* mm) : BaseDesignElement()
{
  Init();
  this->elem = elem;
  this->specialResult.Resize(9, 0.0);
  this->index_ = index;
  this->multimaterial = mm;
  this->interfaceDrivenLoadGrad_.Resize(4 * (domain->GetGrid()->GetDim()-1),0.0);

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
}


DesignElement::Type DesignElement::Default(const Context* ctxt)
{
  switch(ctxt->ToApp()) // will fail for piezo ?!
  {
  case App::MECH:
    return DENSITY;
  case App::ACOUSTIC:
    return ACOU_DENSITY;
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
  switch(st)
  {
  case OPT_RESULT_1:
    return 0;
  case OPT_RESULT_2:
    return 1;
  case OPT_RESULT_3:
    return 2;
  case OPT_RESULT_4:
    return 3;
  case OPT_RESULT_5:
    return 4;
  case OPT_RESULT_6:
    return 5;
  case OPT_RESULT_7:
    return 6;
  case OPT_RESULT_8:
    return 7;
  case OPT_RESULT_9:
    return 8;
  case OPT_RESULT_10:
    return 9;
  case OPT_RESULT_11:
    return 10;
  case OPT_RESULT_12:
    return 11;
  case OPT_RESULT_13:
    return 12;
  case OPT_RESULT_14:
    return 13;
  case OPT_RESULT_15:
    return 14;
  case OPT_RESULT_16:
    return 15;
  case OPT_RESULT_17:
    return 16;
  case OPT_RESULT_18:
    return 17;
  case OPT_RESULT_19:
    return 18;
  case OPT_RESULT_20:
    return 19;
  default:
    return -1;
  }
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
      || rd.value == PENALIZED_STRESS
      || rd.value == DESIGN_TRACKING
      || rd.value == PROJECTION
      || rd.value == TRANSFO_MATRIX
      || rd.value == SHAPE_MAP_GRAD
      || rd.value == SHAPE_MAP_RELEVANT)
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
        out[0] = GetPhysicalDesign(NULL);
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
        design_filter_grad = f != NULL ? f->ForDensityFiltering() : true;
      if(vs == DesignElement::CONSTRAINT_GRADIENT)
        design_filter_grad = f != NULL ? f->ForDensityFiltering() : true;
    }
    else
    {
      if(vs == DesignElement::COST_GRADIENT)// TODO also check for separate costs! not only sum!
        sens_filter = true;
      if(vs == DesignElement::CONSTRAINT_GRADIENT)
        sens_filter = f != NULL ? f->ForSensitivityFiltering() : false;
    }
  }

  // we are silent if one wants filtering and it is not possible
  double val = 0.0;
  if(sens_filter)
    val = simp->GetSensitivityFilteredValue(vs, f);
  if(design_filter)
    val = simp->GetDensityFilteredValue(vs, simp->filter[fix].density_);
  if(design_filter_grad)
    val = simp->GetDensityFilteredGradient(vs, f);
  if(!sens_filter && !design_filter && !design_filter_grad)
    val = GetPlainValue(vs, dynamic_cast<Condition*>(f));

  LOG_DBG3(desel) << "GV: " << elem->elemNum << " (" << valueSpecifier.ToString(vs) << ", "
                << (access == PLAIN ? "plain" : "smart") << ", " << (f == NULL ? "null" : f->ToString())
                << ") sf=" << sens_filter << " df=" << design_filter << " dfg=" << design_filter_grad << " -> " << val;
  return val;
}


__attribute__((always_inline)) inline double DesignElement::GetPlainValue(ValueSpecifier sp, Condition* g) const
{
  // validate first:
  switch(sp)
  {
  case DESIGN:                return design;

  case COST_GRADIENT:         return SumObjectiveGradient();
  case WEIGHT:
    if(simp == NULL) throw Exception("'" + valueSpecifier.ToString(sp) + "' is specific to SIMP");
    return simp->DetermineFilter().weight;

  case NUM_NEIGHBOURS:       
    if(simp == NULL) throw Exception("'" + valueSpecifier.ToString(sp) + "' is specific to SIMP");
    return simp->DetermineFilter().neighborhood.GetSize();

  case CONSTRAINT_GRADIENT:
    assert(g != NULL);
    return constraintGradient[g->GetIndex()];

  case MAX_SLOPE:
  case MAX_MOLE:
  case MAX_OSCILLATION:
  case MAX_JUMP:
  case PENALIZED_STRESS:
  case TRANSFO_MATRIX:
    assert(false); // should be covered before by special result index
    break; // only for the compiler

  case TOPGRAD_VALUE:
    if(tge == 0) throw Exception("'" + valueSpecifier.ToString(sp) + "' is specific to TopGrad");
    return tge->value;

  case SHAPEGRAD_VALUE:
    if(lse_ == 0) throw Exception("'" + valueSpecifier.ToString(sp) + "' is specific to Levelset");
    return lse_->shapeGradValue;
  default: throw Exception(valueSpecifier.ToString(sp) + " is no scalar value");
  }
  return -1.0; // cannot happen due to default in switch but to please the compiler :(
}


double DesignElement::GetDesign(Access access) const
{
  return GetValue(DESIGN, access);
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

  return tf->Transform(trans != NULL ? trans : this, SMART, false); // if there is a transformation return the transformed stuff
}


bool DesignElement::HasPhysicalDesign() const
{
  return(type_ == DENSITY || type_ == POLARIZATION || type_ == ACOU_DENSITY || (!simp->filter.IsEmpty() && simp->filter[0].GetType() == Filter::DENSITY));
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

std::string DesignElement::ToString(const DesignElement* de, bool barycenter)
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
    ss << "e=" << boost::lexical_cast<std::string>(de->elem->elemNum);
    if(barycenter)
      ss << " bc=" << de->elem->barycenter.ToString();
    else
      ss << " t=" << type.ToString(de->type_);
  }
  
  return ss.str();
}

std::string DesignElement::ToString(const StdVector<DesignElement*>& vec, bool print_type)
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

std::string DesignElement::ToString(const StdVector<DesignElement>& vec, bool print_val, bool print_type)
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
  Filter::sensitivity.Add(Filter::SHARP_SIGMUND, "sharp_sigmund");
  Filter::sensitivity.Add(Filter::BORRVALL, "borrvall");

  Filter::density.SetName("Filter::Density");
  Filter::density.Add(Filter::STANDARD, "standard");
  Filter::density.Add(Filter::SOLID_HEAVISIDE, "solid_heaviside");
  Filter::density.Add(Filter::VOID_HEAVISIDE, "void_heaviside");
  Filter::density.Add(Filter::TANH, "tanh");

  ShapeMapDesign::type.SetName("ShapeMapDesign::Type");
  ShapeMapDesign::type.Add(ShapeMapDesign::NODE, "node");
  ShapeMapDesign::type.Add(ShapeMapDesign::PROFILE, "profile");

  ShapeMapDesign::symmetry.SetName("ShapeMapDesign::Symmetry");
  ShapeMapDesign::symmetry.Add(ShapeMapDesign::NONE, "none");
  ShapeMapDesign::symmetry.Add(ShapeMapDesign::MIRROR, "mirror");


  type.SetName("BaseDesignElement::Type");
  type.Add(NO_TYPE, "no_type");
  type.Add(NO_MULTIMATERIAL, "no_multimaterial");
  type.Add(NO_DERIVATIVE, "no_derivative");
  type.Add(MECH_TRACE, "mech_trace");
  type.Add(MECH_ALL, "mech_all");
  type.Add(DIELEC_TRACE, "dielec_trace");
  type.Add(DIELEC_ALL, "dielec_all");
  type.Add(G_ALL, "G_all");
  type.Add(PIEZO_ALL, "piezo_all");
  type.Add(DEFAULT, "default");
  type.Add(DENSITY, "density");
  type.Add(ACOU_DENSITY, "acouDensity");
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
  type.Add(MECH_22, "mech_22");
  type.Add(MECH_33, "mech_33");
  type.Add(MECH_23, "mech_23");
  type.Add(MECH_13, "mech_13");
  type.Add(MECH_12, "mech_12");
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
  type.Add(ROTANGLE2, "rotAngle2");
  type.Add(SCALING1, "scaling1");
  type.Add(SCALING2, "scaling2");
  type.Add(G11, "G11");
  type.Add(G12, "G12");
  type.Add(G21, "G21");
  type.Add(G22, "G22");
  type.Add(G_MAP_X, "G_MAP_X");
  type.Add(G_MAP_Y, "G_MAP_Y");
  type.Add(GX_0, "GX_0");
  type.Add(GY_0, "GY_0");
  type.Add(GX_PX, "GX_PX");
  type.Add(GY_PX, "GY_PX");
  type.Add(GX_PY, "GX_PY");
  type.Add(GY_PY, "GY_PY");
  type.Add(GX_PXY, "GX_PXY");
  type.Add(GY_PXY, "GY_PXY");
  type.Add(ROTANGLEX, "rotAngleX");
  type.Add(ROTANGLEY, "rotAngleY");
  type.Add(ROTANGLEZ, "rotAngleZ");
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
  type.Add(ALL_DESIGNS, "allDesigns");

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
  valueSpecifier.Add(DESIGN_TRACKING, "designTracking");
  valueSpecifier.Add(WEIGHT, "weight");
  valueSpecifier.Add(OBJECTIVE, "objective");
  valueSpecifier.Add(PROJECTION, "projection");
  valueSpecifier.Add(TRANSFO_MATRIX, "transfoMatrix");
  valueSpecifier.Add(NUM_NEIGHBOURS, "neighbours");
  valueSpecifier.Add(LEVEL_SET_VALUE, "levelSetValue");
  valueSpecifier.Add(LEVEL_SET_STATE, "levelSetState");
  valueSpecifier.Add(TOPGRAD_VALUE, "topGradValue");
  valueSpecifier.Add(SHAPEGRAD_VALUE, "shapeGradValue");
  valueSpecifier.Add(SHAPEGRAD_NODE_VALUE, "shapeGradNodeValue");
  valueSpecifier.Add(SHAPE_MAP_GRAD, "shapeMapGrad");
  valueSpecifier.Add(SHAPE_MAP_RELEVANT, "shapeMapRelevant");
  valueSpecifier.Add(LEVEL_SET_GRAD_XP, "levelSetGradXP");
  valueSpecifier.Add(LEVEL_SET_GRAD_XN, "levelSetGradXN");
  valueSpecifier.Add(LEVEL_SET_GRAD_YP, "levelSetGradYP");
  valueSpecifier.Add(LEVEL_SET_GRAD_YN, "levelSetGradYN");
  valueSpecifier.Add(LEVEL_SET_GRAD_ZP, "levelSetGradZP");
  valueSpecifier.Add(LEVEL_SET_GRAD_ZN, "levelSetGradZN");
  valueSpecifier.Add(HEAT_NODAL_TRACK_VAL, "heatNodalTrackValue");
  valueSpecifier.Add(TEMP_AT_INTERFACE, "tempAtInterface");

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
  detail.Add(PROJECTION_FILTER, "projectionFilter");
  detail.Add(TRANSFO_MATRIX11, "transfoMatrix11");
  detail.Add(TRANSFO_MATRIX12, "transfoMatrix12");
  detail.Add(TRANSFO_MATRIX21, "transfoMatrix21");
  detail.Add(TRANSFO_MATRIX22, "transfoMatrix22");
  detail.Add(SM_NODE, "node");
  detail.Add(SM_PROFILE, "profile");
}


SIMPElement::SIMPElement(DesignElement* base)
{
  this->de_ = base;
}

double SIMPElement::GetSensitivityFilteredValue(DesignElement::ValueSpecifier sp, Function* g) const
{
  // We filter over this element and the neighbors.
  assert(de_->simp != NULL);
  assert(de_->simp->filter.GetSize() == 1); // no robustness for sensitivity filtering worth implementing

  const Filter& f = filter[DetermineFilterIndex()];
  assert(f.GetType() == Filter::SENSITIVITY);

  // See Filter::Sensitivity: w = weight, p is density, f' is cost gradient
  // Sigmund  = sum_i w(x_i) * p_i * f_i' / p_e * sum_i w(x_i)
  // Sharp Sigmund  = sum_i (i=e ? 1:0 : w(x_i)) * p_i * f_i' / p_e * sum_i w(x_i) + "bug" in normalized weighting
  // Borrvall = sum_i w(x_i) * p_i * f_i' / sum_i p_i * w(x_i)
  // plain    = sum_i w(x_i) * f_i' / sum_i w(x_i)
  // sharp plain = plain, just the filter is setup with normalized weights with a "bug"

  // factor design in numerator (SIGMUND and BORRVALL) for densitiy filtering the value is any in
  bool numerator_design = f.sensitivity_ != Filter::PLAIN && f.sensitivity_ != Filter::SHARP_PLAIN;
  // factor design in denominator sum (BORRVALL)
  bool denominator_design = f.sensitivity_ == Filter::BORRVALL;
  // weight the denominator by this design (SIGMUND)
  bool sigmund_denominator = f.sensitivity_ == Filter::SIGMUND || f.sensitivity_ == Filter::SHARP_SIGMUND;
  // short-cut for fake in Sharp Sigmund: i=e ? 1:0 : w(x_i)
  bool cheat_this_weight = f.sensitivity_ == Filter::SHARP_SIGMUND || f.sensitivity_ == Filter::SHARP_PLAIN;

  LOG_DBG3(desel) << "FGV: el=" << de_->elem->elemNum
                << " sp=" << DesignElement::valueSpecifier.ToString(sp)
                << " dens=" << Filter::density.ToString(f.density_)
                << " numerator_design=" << numerator_design
                << " denominator_design=" << denominator_design
                << " sigmund_denominator=" << sigmund_denominator
                << " cheat_this_weight=" << cheat_this_weight;


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
    double x = de->GetPlainValue(DesignElement::DESIGN); // cheap if not used

    double numerator_summand   = nw * v * (numerator_design ? x : 1.0);
    double denominator_summand = w * (denominator_design ? x : 1.0);

    numerator   += numerator_summand;
    denominator += denominator_summand;

    LOG_DBG3(desel) << "FGV: el=" << de_->elem->elemNum << ": curr=" << de->elem->elemNum
                  << " w= " << w  << " nw=" << nw << " x=" << x << " v=" << v << " ns=" << numerator_summand
                  << " ds=" << denominator_summand << " num=" << numerator << " den=" << denominator;
  }

  if(sigmund_denominator)
    denominator *= de_->GetDesign(DesignElement::PLAIN);

  LOG_DBG3(desel) << "FGV: el=" << de_->elem->elemNum << ": result=" << numerator
                << "/" << denominator << " -> " << (numerator / denominator);

  return numerator / denominator;
}




double SIMPElement::GetDensityFilteredValue(DesignElement::ValueSpecifier sp, Filter::Density fd) const
{
  // We filter over this element and the neighbors.
  assert(de_->simp != NULL);
  assert(sp == DesignElement::DESIGN);
  assert(!de_->simp->filter.IsEmpty());

  unsigned int fix = DetermineFilterIndex();
  const Filter& f = filter[fix];


  // All equations from Sigmund; Morphology based black and white filters for topology optimization; 2007
  // p = rho. P is filtered rho (rho tilde)
  // P = sum_(i in N_e) w(x_i) p_i / sum_(i in N_e) w(x_i)


  // mathematically the neighborhood includes this element, but this is not in the structure
  // we initialize numerator and denominator with the values obtained from this element
  double numerator = f.weight * this->de_->GetPlainValue(DesignElement::DESIGN);
  double denominator = f.weight;

  LOG_DBG3(desel) << "GDFV: el=" << de_->elem->elemNum << ": curr=" << de_->elem->elemNum
                   << " w= " << f.weight << " x=" << this->de_->GetPlainValue(DesignElement::DESIGN)
                   << " num=" << numerator << " den=" << denominator << " fix=" << fix;

  for(int i = 0, ni = (int) f.neighborhood.GetSize(); i < ni; i++)
  {
    const Filter::NeighbourElement* ne = &f.neighborhood[i];
    const DesignElement* de = ne->neighbour;

    double w = ne->weight;
    double x = de->GetPlainDesignValue();

    numerator   += w * x;
    denominator += w;

     LOG_DBG3(desel) << "GDFV: el=" << de_->elem->elemNum << ": curr=" << de->elem->elemNum  << " w= " << w  << " x=" << x << " num=" << numerator << " den=" << denominator;
  }

  double p_filt = numerator / denominator;

  LOG_DBG3(desel) << "GDFV: el=" << de_->elem->elemNum << " filtered_density=" << p_filt;

  assert(fd == Filter::STANDARD || fd == Filter::SOLID_HEAVISIDE || fd == Filter::VOID_HEAVISIDE || fd == Filter::TANH);

  // do we need to post proc?
  if(fd != Filter::STANDARD)
  {
    assert(fd == Filter::SOLID_HEAVISIDE || fd == Filter::VOID_HEAVISIDE  || fd == Filter::TANH);

    if(fd == Filter::TANH)
      p_filt = CalcTanh(p_filt, fix);
    else
      p_filt = CalcHeaviside(p_filt, fix);

    assert(p_filt <= 1.000001 * this->de_->GetUpperBound());
    assert(p_filt >= 0.7 * this->de_->simp->filter[fix].GetLowerBound(this->de_)); // relax the assert a little, cause of heaviside correction
  }

  LOG_DBG3(desel) << "GDFV: el=" << de_->elem->elemNum << " design=" << Filter::density.ToString(de_->simp->filter[fix].density_)
                   << ": plain=" << this->de_->GetPlainValue(DesignElement::DESIGN) << " -> "<< p_filt;

  return p_filt;
}

double SIMPElement::GetDensityFilteredGradient(DesignElement::ValueSpecifier sp, Function* func) const
{
  // We filter over this element and the neighbors.
  assert(de_->simp != NULL);

  unsigned int fix = DetermineFilterIndex();
  const Filter& f = filter[fix];

  assert(f.GetType() == Filter::DENSITY);
  assert(sp == DesignElement::COST_GRADIENT || sp == DesignElement::CONSTRAINT_GRADIENT);
  assert((func == NULL || (func->IsObjective() && sp == DesignElement::COST_GRADIENT)) || (func == NULL || (!func->IsObjective() && sp == DesignElement::CONSTRAINT_GRADIENT)) || (func == NULL || (func->IsObjective())));
  // projection has density filtering only in the fake filter problem but not in the original problem (which should not be density filtered anyway)
  assert(func == NULL || func->ForDensityFiltering());


  // Density filtering for gradient is (Sigmund; Morphology-based black and white filters for topology optimization; 2007; eqn (35). (36)
  // p is rho and P is rho filtered! d f/d p_e = sum_i(in N_e) d f/d P_i * d P_i/d p_e with d P_i/d p_e = w(x_e)/ sum_j(in N_i) w(x_j)
  // note, that the stored value is already v = d f/d P_i

  //LOG_DBG3(desel) << "GDFG: el=" << de_->elem->elemNum
                //<< " sp=" << DesignElement::valueSpecifier.ToString(sp)
                //<< " g=" << (g != NULL ? Condition::type.ToString(g->GetType()) : "null");

  double sum = 0.0;

  if(f.density_ == Filter::STANDARD)
  {
    for(int i = -1, ni = (int) f.neighborhood.GetSize(); i < ni; i++)
    {
      const Filter::NeighbourElement* ne = i == -1 ? NULL : &f.neighborhood[i];
      const DesignElement* de = i == -1 ? this->de_ : ne->neighbour;

      double v = de->GetPlainValue(sp, dynamic_cast<Condition*>(func)); // d f/d P_i

      double w = i == -1 ? f.weight : ne->weight;

      if (de->simp->filter[fix].weight_sum < 0.0)
        de->simp->filter[fix].weight_sum = de->simp->filter[fix].CalcWeightSum(true);

      double summand = v * w / de->simp->filter[fix].weight_sum;
      sum += summand;

      // LOG_DBG3(desel) << "GDFG: el=" << de_->elem->elemNum << ": curr=" << de->elem->elemNum
      //                << " v= " << v  << " h=" << h << " w=" << w << " x_n=" << x_n << " w_sum=" << w_sum
      //                << " summand=" << summand << " sum=" << sum;
    }
  }
  else // the non Filter::STANDARD case
  {
    // mathematically the neighborhood includes this element, but this is not in the structure
    for(int i = -1, ni = (int) f.neighborhood.GetSize(); i < ni; i++)
    {
      const Filter::NeighbourElement* ne = i == -1 ? NULL : &f.neighborhood[i];
      const DesignElement* de = i == -1 ? this->de_ : ne->neighbour;

      double v = de->GetPlainValue(sp, dynamic_cast<Condition*>(func)); // d f/d P_i
      double h = de->simp->filter[fix].non_lin_scale; // for not-standard filters this is the additional derivative
      double x_n = 0.0;

      double b = f.beta;

      // we need the filtered density -> but the real filtered value!!
      x_n = de->simp->GetDensityFilteredValue(DesignElement::DESIGN, Filter::STANDARD);

      if(f.density_ == Filter::SOLID_HEAVISIDE)
      {
       // F = 1.0 - std::exp(-b * input_value) + input_value * std::exp(-b)
       // H = scale * F + offset
        h *= (b * exp(-b * x_n) + exp(-b));
      }
      if(f.density_ == Filter::VOID_HEAVISIDE)
      {
        // general scaling
        h *= (b * exp(b*(x_n-1.0)) + exp(-1.0*b));
      }
      if(f.density_ == Filter::TANH)
      {
        // f(x)  =  1 - 1/(exp(2*beta*(x-param)) + 1)
        // f'(x) =  (exp(2*beta*(x-param)+1)^-2 * 2 * beta * exp(2*beta*(x-param))
        double eta = f.eta;
        double e = std::exp(2.0 * b * ( x_n - eta));

        h *= (1.0/((e+1.0)*(e+1.0)) * 2.0 * b * e);
      }

      double w = i == -1 ? f.weight : ne->weight;

      if (de->simp->filter[fix].weight_sum < 0.0)
        de->simp->filter[fix].weight_sum = de->simp->filter[fix].CalcWeightSum(true);

      double summand = v * h * w / de->simp->filter[fix].weight_sum;
      sum += summand;

      // LOG_DBG3(desel) << "GDFG: el=" << de_->elem->elemNum << ": curr=" << de->elem->elemNum
      //                << " v= " << v  << " h=" << h << " w=" << w << " x_n=" << x_n << " w_sum=" << w_sum
      //                << " summand=" << summand << " sum=" << sum;
    }
  }

  return sum;
}

std::string SIMPElement::ToString(int level) const
{
  std::stringstream ss;
  const Filter& f = filter[DetermineFilterIndex()];
  ss << "el=" << de_->elem->elemNum << " #n=" << f.neighborhood.GetSize() << "(";
  for(unsigned int i = 0; i < !filter.IsEmpty() && f.neighborhood.GetSize(); i++)
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

inline Filter& SIMPElement::DetermineFilter()
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
  design.Resize(domain->GetGrid()->GetDim() == 2 ? 4 : 6);
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

    // reference case
    StdVector<std::pair<Elem*, int> >& neighbors = *(de->elem->neighborhood);
    // reuse the enlarged_data element for the periodic case only
    if(periodic)
    {
      if(structure->ExtendPeriodicNeighborhood(de->elem, common, enlarged_data))
      {
        neighbors = enlarged_data; // in the non-periodic case there is no need to copy the element data
        de->vicinity->periodic = true;
      }
    }

    // FIXME LOG_DBG(desel) << "VE:Init elem=" << de->elem->elemNum << " neighbors=" << neighbors.ToString();

    for(unsigned int n = 0; n < neighbors.GetSize(); n++)
    {
      // we consider only direct (edge/face) neighbors
      if(neighbors[n].second < common) continue;
      // now we have to find the relative position of candidate
      Elem* candidate = neighbors[n].first;


      LOG_DBG3(desel) << "VE:Init elem=" << de->elem->elemNum << " e.bc=" << de->elem->barycenter.ToString() << " e.dim="
                    << de->elem->GetShape().dim << " e.r=" << de->elem->regionId << " n=" << n << " o.el=" << candidate->elemNum
                    << " o.bc=" << candidate->barycenter.ToString() << " o.dim=" << candidate->GetShape().dim << " o.r=" << candidate->regionId;

      // if the neighbor is a surface element we don't want to play with it
      if(de->elem->GetShape().dim != candidate->GetShape().dim)
        continue;
      // we include pseudo design regions. This is for off-design optimization or non_design_vicinity=true in <ersatzMaterial>
      if(!space->Contains(candidate->regionId, space->DoNonDesignVicinity()))
        continue;

      // the spacing allows to identify periodic elements
      int idx = FindRelativeNeighborLocation(de->elem->barycenter, candidate->barycenter, spacing);
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
  ss << " periodic=" << periodic;
  ss << ")";
  return ss.str();
}


ResultDescription::ResultDescription()
{
  // does not set all!!
  access = DesignElement::PLAIN;
  value  = DesignElement::DESIGN;
  design = DesignElement::DEFAULT;
  excitation = -1;
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

  excitation = pn->Has("excitation") ? pn->Get("excitation")->As<int>() : -1;

  LOG_DBG(desel) << "RD:RD " << ToString();
}


std::string ResultDescription::ToString()
{
  std::stringstream ss;
  ss << "RD: design=" << DesignElement::type.ToString(design)
     << " access=" << DesignElement::access.ToString(access)
     << " value=" << DesignElement::valueSpecifier.ToString(value)
     // << " detail=" << DesignElement::detail.ToString(detail)
     << " ex=" << excitation;
  return ss.str();
}
