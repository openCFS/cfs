#include "Optimization/DesignSpace.hh"
#include "Optimization/DesignElement.hh"
#include "Optimization/TransferFunction.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/InfoNode.hh"
#include "DataInOut/Logging/cfslog.hh"
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
}

TransferFunction::TransferFunction(ParamNode* pn, DesignElement::Type default_design)
{
  // initialize the static Enum the first time
  if(type.map.size() == 0) SetEnums();
  this->type_ = type.Parse(pn->Get("type"));
  this->orgType_ = NO_TYPE;
  this->application_ = Optimization::application.Parse(pn->Get("application"));
  if(pn->Has("design"))
    this->design_ = DesignElement::type.Parse(pn->Get("design"));
  else
  {
    if(default_design == DesignElement::NO_TYPE) 
      throw Exception("Set the 'design' attribute for 'transferFunction' when using multiple design variables.");
    this->design_ = default_design;
  }
  this->param_ = pn->Has("param") ? pn->Get("param")->AsDouble() : 1.0;
  
  // validate param
  if(type_ == IDENTITY && pn->Has("param"))
    throw Exception("it makes no sense to give a parameter for an identity transfer function");  
  if(type_ == FULL && pn->Has("param"))
    throw Exception("it makes no sense to give a parameter for the 'full' = 1 transfer function");  
  if(!pn->Has("param") && (type_ != IDENTITY && type_ != FULL))  
    throw Exception("transfer function '" + type.ToString(type_) + "' requires a parameter");
  
  // validate design/application
  if(design_ == DesignElement::POLARIZATION && application_ != Optimization::PIEZO_COUPLING)
    throw Exception("transfer functions for 'polarization' can only be '" 
        + Optimization::application.ToString(Optimization::PIEZO_COUPLING) + "'");
}

void TransferFunction::ToInfo(InfoNode* in) const
{
  in->Get("type")->SetValue(type.ToString(type_));
  in->Get("application")->SetValue(Optimization::application.ToString(application_));
  in->Get("design")->SetValue(DesignElement::type.ToString(design_));
  if(type_ != IDENTITY && type_ != FULL)
    in->Get("param")->SetValue(param_);
  
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

double TransferFunction::Transform(DesignElement* de)
{
  double value = de->GetDesign(DesignElement::PLAIN);
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
    assert(param_ == 0.0);
    result = de->GetUpperBound();
    break;
    
  default: throw Exception("type not implemented");
  }          
  
  LOG_DBG3(trans) << "Transform de=" << de->elem->elemNum << " value=" << value << " type=" << type.ToString(type_) << " param=" << param_ << " -> " << result;
  return result;
}     

double TransferFunction::Derivative(DesignElement* de)
{
  double value = de->GetDesign(DesignElement::PLAIN);

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
