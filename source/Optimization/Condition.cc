#include "Optimization/Condition.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/ParamHandling/InfoNode.hh"
#include "General/exception.hh"
#include "General/environment.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"
#include "MatVec/denseMatrix.hh"

using namespace CoupledField;
using boost::lexical_cast;

// instantiation of the static elements
Enum<Condition::Name> Condition::name;
Enum<Condition::Type> Condition::type;


Condition::Condition(ParamNode* pn)
{
  volume_fraction_ = 0.0;
  name_ = name.Parse(pn->Get("name"));
  type_ = type.Parse(pn->Get("type"));
  design = !pn->Has("design") ? DesignElement::DEFAULT :
           DesignElement::type.Parse(pn->Get("design"));
  scaling = pn->Get("scaling")->AsDouble();
  parameter = pn->Has("parameter") ? pn->Get("parameter")->AsDouble() : 0.0;
  active = pn->Get("mode")->AsString() == "constraint";

  // value is not mandatory for all constraints. Check for homogenization later
  if(!pn->Has("value") && active && name_ != HOMOGENIZATION_TENSOR && name_ != HOMOGENIZATION_TRACKING)
    EXCEPTION("constraint '" << name.ToString(name_) << "' is active but no value is given.");
  value = pn->Has("value") ? pn->Get("value")->AsDouble() : -1.0;

  penalty = pn->Get("penalty")->AsDouble();
  region = ALL_REGIONS;

  ReadTensor(pn, matrix_); // is save and sets default

  coord.first  = 0;
  coord.second = 0;

  if(pn->Has("region") && pn->Get("region")->AsString() != "all"){
    region = domain->GetGrid()->RegionNameToId(pn->Get("region")->AsString());
  }

  if(name_ == HOMOGENIZATION_TENSOR || name_ == HOMOGENIZATION_TRACKING)
  {
    // we must not give a value when there is a tensor
    if(name_ == HOMOGENIZATION_TENSOR && pn->Has("tensor") && pn->Has("value"))
      throw Exception("a value must not be given when a tensor is used in a homogenization constraint");

    if(name_ == HOMOGENIZATION_TRACKING && !pn->Has("tensor"))
      throw Exception("a 'tensor' is mandatory for homogenization tracking");

    if(active && !pn->Has("tensor") && !pn->Has("value"))
      throw Exception("a value is mandatory when no tensor is given in a homogenization constraint");
  }

  if(name_ == GAUSS_GREYNESS && parameter == 0.0)
    throw Exception("for gauss grayness an parameter h is mandatory (e.g. 0.2)");

}

bool Condition::ReadTensor(ParamNode* pn, Matrix<double>& matrix)
{
  matrix.Resize(1,1); // minimal size, as 0,0 is not defined.

  if(!pn->Has("tensor")) return false;

  ParamNode* tens = pn->Get("tensor");

  int dim = tens->Get("dim1")->AsInt();
  if(dim != 3 && dim != 6)
    EXCEPTION("The Voigt 'tensor' for homogenizations needs to be 3x3 or 6x6");
  if(tens->Has("dim2") && dim != tens->Get("dim2")->AsInt())
    EXCEPTION("The 'tensor' for homogenization needs to be symmetric");
  if((domain->GetGrid()->GetDim() == 2 && dim != 3) || (domain->GetGrid()->GetDim() == 3 && dim != 6))
    EXCEPTION("The 'tensor' for homogenization needs to be 3x3 for 2D and 6x6 for 3D");

  matrix.Resize(dim, dim);
  ParamTools::AsTensor<double>(tens->Get("real"), dim , dim, matrix);

  return true;
}

void Condition::AddCondition(ParamNode* pn, StdVector<Condition>& constraints, StdVector<Condition>& observation)
{
  Condition tmp = Condition(pn);
  // which list to append
  StdVector<Condition>& list = tmp.active ? constraints : observation;
  // the index is defined only for active constraints
  tmp.index_ = list.GetSize();
  list.Push_back(tmp);

  // note that the pointer becomes invalid by AddSubCondition()
  Condition* g = &list.Last();

  // homogenization are special constraints which constraints which "blow up" to more constraints
  if(g->name_ == HOMOGENIZATION_TENSOR)
  {
    // there has been a extensive test in the constructor
    // do we need to blow-up?
    if(!g->ReadCoord(pn))
    {
      // the first entry is this constraint
      // for the conversion of the indices see Thesis from Ole Sigmund, p 30 and book of Manfred
      // p 42 (first edition)
      g->coord.first = 1;      // 1. ijkl = 1111
      g->coord.second = 1;
      g->value = g->matrix_[1-1][1-1]; // one-based!
      g = g->AppendSubCondition(list, 2,2); // 2. 2222
      if(domain->GetGrid()->GetDim() == 2)
      {
        g = g->AppendSubCondition(list, 1,2); // 3. 1122
        g = g->AppendSubCondition(list, 3,3); // 4. 1212 // would be 6,6 according to book
      }
      else
      {
        g = g->AppendSubCondition(list, 3,3); // 3. 3333
        g = g->AppendSubCondition(list, 1,2); // 4. 1122
        g = g->AppendSubCondition(list, 1,3); // 5. 1133
        g = g->AppendSubCondition(list, 2,3); // 6. 2233
        g = g->AppendSubCondition(list, 6,6); // 7. 1212
        g = g->AppendSubCondition(list, 5,5); // 8. 1313
        g = g->AppendSubCondition(list, 4,4); // 9. 2323
      }
    }
  }
}

Condition* Condition::AppendSubCondition(StdVector<Condition>& list, unsigned int pos_x, unsigned int pos_y)
{
  list.Push_back(*this);
  Condition* sub = &list.Last(); // copy this entry as reference
  sub->coord.first = pos_x;
  sub->coord.second = pos_y;
  sub->value = sub->matrix_[pos_x-1][pos_y-1];
  sub->index_ = list.GetSize() - 1;
  return sub;
}

bool Condition::ReadCoord(ParamNode* pn)
{
  std::string val = pn->Get("coord")->AsString();
  if(val == "all") return false;

  assert(val.size() == 2);
  coord.first  = lexical_cast<unsigned int>(val.at(0));
  coord.second = lexical_cast<unsigned int>(val.at(1));

  return true;
}



std::string Condition::ToString() const
{
  std::ostringstream os;
  
  os << name.ToString(name_);
  if(design != DesignElement::DEFAULT)
    os << "_(" << DesignElement::type.ToString(design) << ")";
  if(name_ == HOMOGENIZATION_TENSOR)
    os << "_" << coord.first << coord.second;
  return os.str();  
}


void Condition::ToInfo(InfoNode* in) const
{
  in->Get("name")->SetValue(name.ToString(name_));
  in->Get("design")->SetValue(DesignElement::type.ToString(design));
  in->Get("mode")->SetValue(active ? "constraint" : "observation");
  if(active)
  {
    in->Get("type")->SetValue(type.ToString(type_));
    if(name_ != HOMOGENIZATION_TRACKING)
      in->Get("value")->SetValue(value);
  }
  if(name_ == HOMOGENIZATION_TENSOR)
    in->Get("tensor_entry")->SetValue(lexical_cast<std::string>(coord.first) + lexical_cast<std::string>(coord.second));
}

Matrix<double>& Condition::GetTensor()
{
  return matrix_;
}


bool Condition::IsForRegion(RegionIdType regionId)
{
  return(region == ALL_REGIONS || region == regionId);
}

