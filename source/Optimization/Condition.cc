#include "Optimization/Condition.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/exception.hh"

using namespace CoupledField;

// instantiation of the static elements
Enum<Condition::Name> Condition::name;
Enum<Condition::Type> Condition::type;


Condition::Condition(ParamNode* pn, int index)
{
  index_ = index;
  name_ = name.Parse(pn->Get("name"));
  type_ = type.Parse(pn->Get("type"));
  design = !pn->Has("design") ? DesignElement::DEFAULT : 
           DesignElement::type.Parse(pn->Get("design"));
  value = pn->Get("value")->AsDouble();
  scaling = pn->Get("scaling")->AsDouble();
  parameter = pn->Has("parameter") ? pn->Get("parameter")->AsDouble() : 0.0;
  active = pn->Get("mode")->AsString() == "constraint";
  
  if(name_ == GAUSS_GREYNESS && parameter == 0.0)
    throw Exception("for gauss greyness an parameter h is mandatory (e.g. 0.2)"); 
}


std::string Condition::ToString() const
{
  std::ostringstream os;
  
  os << name.ToString(name_);
  if(design != DesignElement::DEFAULT)
    os << "_(" << DesignElement::type.ToString(design) << ")";
    
  return os.str();  
}

