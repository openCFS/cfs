#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <ostream>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/Logging/log.hpp"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Domain/ElemMapping/Elem.hh"
#include "General/Exception.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/TransferFunction.hh"
#include "PDE/SinglePDE.hh"


using namespace CoupledField;


DECLARE_LOG(trans)
DEFINE_LOG(trans, "transferFunction")


Enum<TransferFunction::Type> TransferFunction::type;

TransferFunction::TransferFunction()
{
  // do nothing - is for StdVector only
  // this constructor is also used for providing a identity transfer function during parametric material optimization
  if(type.map.size() == 0) SetEnums();
  type_ = IDENTITY;
  orgType_ = NO_TYPE;
  design_ = DesignElement::DEFAULT;
  application_ = Optimization::NO_APP;
  param_ = 0.0;
  beta_ = -1.0;
  scaling_ = 1.0;
  offset_ = 0.0;
}


TransferFunction::TransferFunction(Optimization::Application app, TransferFunction::Type tf_type, double param, DesignElement::Type design)
{
  assert(type.map.size() > 0);

  this->type_ = tf_type;
  this->orgType_ = NO_TYPE;
  this->application_ = app;
  this->design_ = design;
  this->param_ = param;
  this->beta_ = -1.0;
  this->scaling_ = 1.0;
  this->offset_ = 0.0;
}


TransferFunction::TransferFunction(PtrParamNode pn, DesignElement::Type default_design)
{
  // initialize the static Enum the first time
  if(type.map.size() == 0) SetEnums();
  this->type_ = type.Parse(pn->Get("type")->As<std::string>());
  this->orgType_ = NO_TYPE;
  this->application_ = Optimization::application.Parse(pn->Get("application")->As<std::string>());
  this->scaling_ = 1.0;
  this->offset_ = 0.0;
  if(pn->Has("design"))
    this->design_ = DesignElement::type.Parse(pn->Get("design")->As<std::string>());
  else
  {
    if(default_design == DesignElement::NO_TYPE) 
      throw Exception("Set the 'design' attribute for 'transferFunction' when using multiple design variables.");
    this->design_ = default_design;
  }
  this->param_ = pn->Has("param") ? pn->Get("param")->As<double>() : 1.0;
  this->beta_ = pn->Has("beta") ? pn->Get("beta")->As<double>() : -1.0;
  
  // validate param
  if(!pn->Has("param") && (type_ != IDENTITY && type_ != FULL))  
    throw Exception("transfer function '" + type.ToString(type_) + "' requires a parameter");
  
  if((type_ == HEAVISIDE || type_ == TANH) && (!pn->Has("param") || !pn->Has("beta")))
    throw Exception("transfer function '" + type.ToString(type_) + "' requires a 'param' and 'beta' to be set");

  if(type_ == TANH && (param_ < 0.0 || param_ > 1.0 )) // ignore exotic design variables!
    throw Exception("'param' for transfer function 'tanh' out of bound [0:1]");

  // validate design/application
  if(design_ == DesignElement::POLARIZATION && application_ != Optimization::PIEZO_COUPLING)
    throw Exception("transfer functions for 'polarization' can only be '" 
        + Optimization::application.ToString(Optimization::PIEZO_COUPLING) + "'");

  LOG_DBG(trans) << "TF::TF " << ToString();
}



Optimization::Application TransferFunction::Default(const SinglePDE* pde)
{
  if(pde->GetName() == "electrostatic")    return Optimization::ELEC;
  if(pde->GetName() == "mechanic")         return Optimization::MECH;
  if(pde->GetName() == "heatConduction")   return Optimization::LAPLACE;
  if(pde->GetName() == "acoustic")         return Optimization::LAPLACE;
  if(pde->GetName() == "LatticeBoltzmann") return Optimization::LBM;
  throw Exception("invalid");
}

/** see the other Default */
Optimization::Application TransferFunction::Default(DesignElement::Type type, const SinglePDE* pde)
{
  switch(type)
  {
  case DesignElement::DENSITY:
  {
    if(pde)
      return Default(pde);
  }
  case DesignElement::EMODUL:
  case DesignElement::EMODULISO:
  case DesignElement::GMODUL:
  case DesignElement::POISSON:
  case DesignElement::POISSONISO:
  case DesignElement::ROTANGLE:
  case DesignElement::ROTANGLEX:
  case DesignElement::ROTANGLEY:
  case DesignElement::ROTANGLEZ:
  case DesignElement::STIFF1:
  case DesignElement::STIFF2:
  case DesignElement::STIFF3:
  case DesignElement::MECH_11:
  case DesignElement::MECH_12:
  case DesignElement::MECH_22:
  case DesignElement::MECH_33:
  case DesignElement::SHEAR1:
  case DesignElement::ROTANGLE2:
  case DesignElement::SCALING1:
  case DesignElement::SCALING2:
  case DesignElement::G_MAP_X:
  case DesignElement::G_MAP_Y:
  case DesignElement::MULTIMATERIAL:
  case DesignElement::INTERPOLATION:
    return Optimization::MECH;
  case DesignElement::ACOU_DENSITY:
    return Optimization::LAPLACE;
  case DesignElement::POLARIZATION:
    return Optimization::PIEZO_COUPLING;
  default:
    throw Exception("invalid request for transfer function");
  }
}

void TransferFunction::ToInfo(PtrParamNode in) const
{
  in->Get("type")->SetValue(type.ToString(type_));
  in->Get("application")->SetValue(Optimization::application.ToString(application_));
  in->Get("design")->SetValue(DesignElement::type.ToString(design_));
  if(type_ != IDENTITY && type_ != FULL)
    in->Get("param")->SetValue(param_);
  if(type_ == HEAVISIDE || type_ == TANH)
    in->Get("beta")->SetValue(beta_);
  if(scaling_ != 1.0 || offset_ != 0.0)
  {
    assert(type_ == HEAVISIDE || type_ == TANH);
    in->Get("scaling")->SetValue(scaling_);
    in->Get("offset")->SetValue(offset_);
  }

  
}

void TransferFunction::Enable(bool enable)
{
  // try to handle to much toggling cases
  if(enable)
  {
    type_ = orgType_;
    assert(type_ != NO_TYPE);
  }
  else
  {
     // only disable if we are enabled
     assert(type_ != NO_TYPE);
    orgType_ = type_;
    type_ = NO_TYPE;
  }
}


bool TransferFunction::IsPenalized() const
{
  switch(type_)
  {
  case SIMP_TYPE:
    return param_ != 1.0;

  case RAMP:
    return param_ != 0.0;

  case FIXED:
  case FULL:
  case IDENTITY:
    return false;

  case HEAVISIDE:
  case TANH:
    return true;

  case NO_TYPE:
    EXCEPTION("wrong");
  }

  return false; // stupid C++
}

std::string TransferFunction::ToString()
{
  std::ostringstream os;
  os << "TransferFunction[type= " << type.ToString(type_) 
     << "; application=" << Optimization::application.ToString(application_) 
     << "; design=" << DesignElement::type.ToString(design_)
     << "; param=" << param_  
     << (orgType_ == NO_TYPE ? "; enabled" : "; DISABLED!") << "]";
  return os.str();   
}

double TransferFunction::Transform(const DesignElement* de, DesignElement::Access access, bool lower_bimat) const
{
  double value = de->GetValue(DesignElement::DESIGN, access);
  return this->Transform(value, lower_bimat, de);
}

double TransferFunction::Transform(double value, bool lower_bimat, const DesignElement* de) const
{
  if(lower_bimat)
    value = 1.0 - value;

  double result;
  switch(type_)
  {
  case NO_TYPE: // if disabled
    result = 1.0;
    break;
  
  case IDENTITY: 
    result = value;
    break;

  case SIMP_TYPE:
    assert(param_ >= 0);
    result = std::pow(value, param_);
    break;

  case RAMP:
    assert(param_ >= 0);
    result = value / (1.0 + param_ * (1.0 - value));
    break;

  case FIXED:
    result = param_;
    break;

  case FULL:
    assert(de != NULL);
    result = de->GetUpperBound();
    break;

  case HEAVISIDE:
    assert(!lower_bimat); // first check what we do!
    // some options and the derivatives
    // plot (1-exp(-20*x)), 20*x*exp(-20*x), 4*(1-exp(-10*x))**3 * 10*x*exp(-10*x), (1-exp(-10*x))**4, 1-exp(-20*x**6), 20*x**6*6*x**5*exp(-20*x**6)
    assert(beta_ >= 0.0);
    // we optionally scale the stuff when we have physical design
    result = scaling_ * (std::pow(1.0 - std::exp(-1.0 * beta_ * value), param_)) + offset_;
    break;

  case TANH:
    assert(!lower_bimat); // first check what we do!
    assert(beta_ >= 0.0);
    // tf(x) =  1 - 1/(exp(2*beta*(x-param)) + 1)
    // we optionally scale the stuff when we have physical design
    // tf(x) =  scaling * (1 - 1/(exp(2*beta*(x-param)) + 1)) + offset
    result = scaling_ * (1.0 - 1.0/(std::exp(2.0 * beta_ * ( value - param_)) + 1.0)) + offset_;
    break;

  default: throw Exception("type not implemented");
  }          
  
  //LOG_DBG3(trans) << "Transform de=" << (de != NULL ? (int) de->elem->elemNum : -1) << " value=" << value << " type=" << type.ToString(type_) << " param=" << param_ << " -> " << result;
  return result;
}     

double TransferFunction::Derivative(const DesignElement* de, DesignElement::Access access, bool lower_bimat) const
{
  double value = de->GetValue(DesignElement::DESIGN, access);
#ifdef CHECK_INDEX
  if(de->GetType() != design_ && (design_ == DesignElement::DEFAULT && de->GetDesignSpace() != NULL && de->GetDesignSpace()->design.GetSize() > 1))
    throw Exception("type mismatch for the transfer function");
#endif
  return this->Derivative(value, lower_bimat);
}

double TransferFunction::Derivative(double value, bool lower_bimat) const
{
    switch(type_)
    {
    case NO_TYPE:
    case IDENTITY:
      if(!lower_bimat)
        return 1.0;
      else
        return -1.0;

    case SIMP_TYPE:
      if (!lower_bimat)
        return param_ * std::pow(value, param_ - 1.0);
      else
        return - param_ * std::pow(1.0-value, param_ - 1.0);

    case RAMP: 
    {
      double p2 = param_ * param_;
      double p = param_;
      double x = value;
      if(!lower_bimat)
        return (1.0 + p) / (x*x*p2 - 2*x*(p2+p)+p2+2*p+1);
      else
        return -1.0/(p*x+1)-(p*(1-x))/pow(p*x+1,2);
    } 
      
    case HEAVISIDE:
    {
      assert(!lower_bimat);
      // f = (1-hs)^hp,
      assert(beta_ > 0.0);
      double hs = std::exp(-1.0 * beta_ * value);

      return scaling_ * (param_ * std::pow(1.0 - hs, param_ - 1.0) * beta_ * hs);
    }
    case TANH:
    {
      assert(!lower_bimat);
      // tf(x)  =  1 - 1/(exp(2*beta*(x-param)) + 1)
      // tf'(x) =  (exp(2*beta*(x-param)+1)^-2 * 2 * beta * exp(2*beta*(x-param))
      double e = std::exp(2.0 * beta_ * ( value - param_));
      return scaling_ * (1.0/((e+1.0)*(e+1.0)) * 2.0 * beta_ * e);
    }

    case FIXED:  
    case FULL:  
      return 0.0; // the derivative of a constant is 0 
      
    default: throw Exception("type not implemented");
  }            
}     
     
void TransferFunction::SetEnums()
{
  type.SetName("TransferFunction::Type");
  type.Add(NO_TYPE, "no_type");  
  type.Add(SIMP_TYPE, "simp");
  type.Add(IDENTITY, "identity");
  type.Add(RAMP, "ramp");
  type.Add(FIXED, "fixed");
  type.Add(FULL, "full");
  type.Add(HEAVISIDE, "heaviside");
  type.Add(TANH, "tanh");
}     
