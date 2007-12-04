#include "Optimization/DesignSpace.hh"
#include "Optimization/DesignElement.hh"
#include "Optimization/TransferFunction.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
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
}

TransferFunction::TransferFunction(ParamNode* pn)
{
  // initialize the static Enum the first time
  if(type.map.size() == 0) SetEnums();
  this->type_ = type.Parse(pn->Get("type"));
  this->orgType_ = NO_TYPE;
  this->application_ = Optimization::application.Parse(pn->Get("application"));
  this->design_ = DesignElement::type.Parse(pn->Get("design"));
  this->param_ = pn->Has("param") ? pn->Get("param")->AsDouble() : 1.0;
  if(type_ == IDENTITY && pn->Has("param"))
    throw Exception("it makes no sense to give a parameter for an identity transfer function");  
  if(type_ != IDENTITY && !pn->Has("param"))  
    throw Exception("Non identity transfer functions need a parameter");
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
    assert(value >= 0);
    assert(param_ >= 0);
    result = std::pow(value, param_);
    break;

  case RAMP:
    assert(param_ >= 0);
    result = value / (1.0 + param_ * (1.0 - value));
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
      return (1.0 / (1.0 + param_ * (1.0 - value))) - (value / param_);
      
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
}     
