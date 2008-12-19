#include "Optimization/Condition.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/InfoNode.hh"
#include "General/exception.hh"
#include "General/environment.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"

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
  penalty = pn->Get("penalty")->AsDouble();
  region = ALL_REGIONS;
  if(pn->Has("region") && pn->Get("region")->AsString() != "all"){
    region = domain->GetGrid()->RegionNameToId(pn->Get("region")->AsString());
  }
  volume_fraction_ = 0.0;
  
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


void Condition::ToInfo(InfoNode* in) const
{
  in->Get("name")->SetValue(name.ToString(name_));
  in->Get("design")->SetValue(DesignElement::type.ToString(design));
  std::string mode = active ? "constraint" : "observation";
  in->Get("mode")->SetValue(mode);
  if(active)
  {
    in->Get("type")->SetValue(type.ToString(type_));
    in->Get("value")->SetValue(value);
  }
}

bool Condition::IsForRegion(RegionIdType regionId){
  return(region == ALL_REGIONS || region == regionId);
}

