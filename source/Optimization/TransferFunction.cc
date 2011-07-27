#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/TransferFunction.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "PDE/SinglePDE.hh"
#include "General/exception.hh"
#include "Utils/StdVector.hh"


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
}


TransferFunction::TransferFunction(PtrParamNode pn, DesignElement::Type default_design)
{
  // initialize the static Enum the first time
  if(type.map.size() == 0) SetEnums();
  this->type_ = type.Parse(pn->Get("type")->As<std::string>());
  this->orgType_ = NO_TYPE;
  this->application_ = Optimization::application.Parse(pn->Get("application")->As<std::string>());
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
}



Optimization::Application TransferFunction::Default(const SinglePDE* pde)
{
  if(pde->GetName() == "electrostatic") return Optimization::ELEC;
  if(pde->GetName() == "mechanic") return Optimization::MECH;
  if(pde->GetName() == "heatConduction") return Optimization::LAPLACE;
  if(pde->GetName() == "acoustic") return Optimization::LAPLACE;
  throw Exception("invalid");
}

/** see the other Default */
Optimization::Application TransferFunction::Default(DesignElement::Type type)
{
  switch(type)
  {
  case DesignElement::DENSITY:
    return Optimization::MECH;
  case DesignElement::ACOU_DENSITY:
    return Optimization::LAPLACE;
  case DesignElement::POLARIZATION:
    return Optimization::PIEZO_COUPLING;
  default:
    throw Exception("invalid");
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

double TransferFunction::Transform(const DesignElement* de, DesignElement::Access access, double external_value) const
{
  assert(!(external_value != -13.456 && access == DesignElement::SMART));
  double value = external_value == -13.456 ? de->GetValue(DesignElement::DESIGN, access) : external_value;
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
    result = de->GetUpperBound();
    break;
    
  case HEAVISIDE:
    // some options and the derivatives
    // plot (1-exp(-20*x)), 20*x*exp(-20*x), 4*(1-exp(-10*x))**3 * 10*x*exp(-10*x), (1-exp(-10*x))**4, 1-exp(-20*x**6), 20*x**6*6*x**5*exp(-20*x**6)
    assert(beta_ >= 0.0);
    result = std::pow(1.0 - std::exp(-1.0 * beta_ * value), param_);
    break;

  case TANH:
    assert(beta_ >= 0.0);
    // tf(x) =  1 - 1/(exp(2*beta*(x-param)) + 1)
    result = 1.0 - 1.0/(std::exp(2.0 * beta_ * ( value - param_)) + 1.0);
    break;

  default: throw Exception("type not implemented");
  }          
  
  LOG_DBG3(trans) << "Transform de=" << de->elem->elemNum << " value=" << value << " type=" << type.ToString(type_) << " param=" << param_ << " -> " << result;
  return result;
}     

double TransferFunction::Derivative(const DesignElement* de, DesignElement::Access access) const
{
  double value = de->GetValue(DesignElement::DESIGN, access);

  #ifdef CHECK_INDEX
    if(de->GetType() != design_) throw Exception("type missmatch");
  #endif
    switch(type_)
    {
    case NO_TYPE:
    case IDENTITY: 
      return 1.0;

    case SIMP_TYPE:
      return param_ * std::pow(value, param_ - 1.0);

    case RAMP: 
    {
      double p2 = param_ * param_;
      double p = param_;
      double x = value;
      return (1.0 + p) / (x*x*p2 - 2*x*(p2+p)+p2+2*p+1);
    } 
      
    case HEAVISIDE:
    {
      // f = (1-hs)^hp,
      assert(beta_ > 0.0);
      double hs = std::exp(-1.0 * beta_ * value);

      return param_ * std::pow(1.0 - hs, param_ - 1.0) * beta_ * hs;
    }
    case TANH:
    {
      // tf(x)  =  1 - 1/(exp(2*beta*(x-param)) + 1)
      // tf'(x) =  (exp(2*beta*(x-param)+1)^-2 * 2 * beta * exp(2*beta*(x-param))
      double e = std::exp(2.0 * beta_ * ( value - param_));
      return 1.0/((e+1.0)*(e+1.0)) * 2.0 * beta_ * e;
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
