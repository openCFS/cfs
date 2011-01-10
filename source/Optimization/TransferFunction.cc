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
}


TransferFunction::TransferFunction(Optimization::Application app, TransferFunction::Type tf_type, double param, DesignElement::Type design)
{
  assert(type.map.size() > 0);

  this->type_ = tf_type;
  this->orgType_ = NO_TYPE;
  this->application_ = app;
  this->design_ = design;
  this->param_ = param;
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
  this->param_ = pn->Has("param") ? pn->Get("param")->As<Double>() : 1.0;
  
  // validate param
  if(!pn->Has("param") && (type_ != IDENTITY && type_ != FULL))  
    throw Exception("transfer function '" + type.ToString(type_) + "' requires a parameter");
  
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
    return Optimization::ELEC;
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
  
}

bool TransferFunction::IsPenalized() const
{
  switch(type_)
  {
  case SIMP_TYPE:
    return param_ != 1.0;

  case RAMP:
    return param_ != 0.0;

  default:
    return false;
  }
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
}     
