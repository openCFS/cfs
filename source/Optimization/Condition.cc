#include "Optimization/Condition.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/exception.hh"

using namespace CoupledField;

// instantiation of the static elements
Enum Condition::name;
Enum Condition::type;


Condition::Condition(ParamNode* pn, int index)
{
  index_ = index;
  name_ = (Condition::Name) name.Parse(pn->Get("name")->AsString());
  type_ = (Condition::Type) type.Parse(pn->Get("type")->AsString());
  design = !pn->Has("design") ? DesignElement::DEFAULT : 
          (DesignElement::Type) DesignElement::type.Parse(pn->Get("design")->AsString());
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

