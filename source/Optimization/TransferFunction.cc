#include "Optimization/DesignSpace.hh"
#include "Optimization/DesignElement.hh"
#include "Optimization/TransferFunction.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/exception.hh"
#include "Utils/StdVector.hh"

using namespace CoupledField;

Enum TransferFunction::type;

TransferFunction::TransferFunction()
{
  // do nothing - is for StdVector only
}

TransferFunction::TransferFunction(ParamNode* pn)
{
  // initialize the static Enum the first time
  if(type.map.size() == 0) SetEnums();
  this->type_ = (Type) type.Parse(pn->Get("type")->AsString());
  this->orgType_ = NO_TYPE;
  this->application_ = (Optimization::Application) Optimization::application.Parse(pn->Get("application")->AsString());
  this->design_ = (DesignElement::Type) DesignElement::type.Parse(pn->Get("design")->AsString());
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
  switch(type_)
  {
    case IDENTITY: 
         return value;
         
    case SIMP_TYPE:
         return std::pow(value, param_);
         
    default: throw Exception("type not implemented");
  }            
}     

double TransferFunction::Derivative(DesignElement* de)
{
  double value = de->GetDesign(DesignElement::PLAIN);

  #ifdef CHECK_INDEX
    if(de->GetType() != design_) throw Exception("type missmatch");
  #endif
  switch(type_)
  {
    case IDENTITY: 
         return 1.0;
         
    case SIMP_TYPE:
         return std::pow(value, param_ - 1.0);
         
    default: throw Exception("type not implemented");
  }            
}     
     
void TransferFunction::SetEnums()
{
  type.SetName("TransferFunction::Type");
  type.Add(SIMP_TYPE, "simp");
  type.Add(IDENTITY, "identity");
  type.Add(RAMP, "ramp");
}     
