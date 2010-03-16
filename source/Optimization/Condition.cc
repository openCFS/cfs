#include "Optimization/Condition.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/exception.hh"
#include "General/environment.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"
#include "MatVec/denseMatrix.hh"
#include "Materials/mechanicMaterial.hh"
#include <sstream>

using namespace CoupledField;
using boost::lexical_cast;

// instantiation of the static elements
Enum<Condition::Name> Condition::name;
Enum<Condition::Type> Condition::type;


Condition::Condition(PtrParamNode pn)
{
  volume_fraction_ = 0.0;
  name_ = name.Parse(pn->Get("name")->As<std::string>());
  type_ = type.Parse(pn->Get("type")->As<std::string>());
  design = !pn->Has("design") ? DesignElement::DEFAULT :
           DesignElement::type.Parse(pn->Get("design")->As<std::string>());
  scaling = pn->Get("scaling")->As<Double>();
  parameter = pn->Has("parameter") ? pn->Get("parameter")->As<Double>() : 0.0;
  active = pn->Get("mode")->As<std::string>() == "constraint";

  delta_logging_ignored_ = false;
  delta_logging = pn->Get("log_delta")->As<bool>();

  if(pn->Get("log_delta")->As<bool>() && (!pn->Has("value") && !pn->Has("tensor") && !pn->Has("isotropic")))
    delta_logging_ignored_ = true;

  // value is not mandatory for all almost all constraints. Check for homogenization later
  if(!pn->Has("value") && active && name_ != HOMOGENIZATION_TENSOR && name_ != HOMOGENIZATION_TRACKING && name_ != ISOTROPY)
    EXCEPTION("constraint '" << name.ToString(name_) << "' is active but no value is given.");
  value = pn->Has("value") ? pn->Get("value")->As<Double>() : -1.0;

  penalty = pn->Get("penalty")->As<Double>();
  region = ALL_REGIONS;

  ReadTensor(pn, matrix_); // is save and sets default

  if(pn->Has("region") && pn->Get("region")->As<std::string>() != "all"){
    region = domain->GetGrid()->GetRegion().Parse(pn->Get("region")->As<std::string>());
  }

  if(name_ == ISOTROPY && pn->Has("value"))
    throw Exception("a value must not be given for a 'isotropy' constraint");

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

bool Condition::ReadTensor(PtrParamNode pn, Matrix<double>& matrix)
{
  matrix.Resize(1,1); // minimal size, as 0,0 is not defined.

  // sanity checks
  if(!pn->Has("tensor") && !pn->Has("isotropic")) return false;
  if(pn->Has("tensor") && pn->Has("isotropic"))
    EXCEPTION("please specify either <tensor> or <isotropic>, not both");

  // check for tensor element
  PtrParamNode tens = pn->Get("tensor", ParamNode::PASS);
  if(tens != NULL)
  {

    int dim = tens->Get("dim1")->As<Integer>();
    if(dim != 3 && dim != 6)
      EXCEPTION("The Voigt 'tensor' for homogenizations needs to be 3x3 or 6x6");
    if(tens->Has("dim2") && dim != tens->Get("dim2")->As<Integer>())
      EXCEPTION("The 'tensor' for homogenization needs to be symmetric");
    if((domain->GetGrid()->GetDim() == 2 && dim != 3) || (domain->GetGrid()->GetDim() == 3 && dim != 6))
      EXCEPTION("The 'tensor' for homogenization needs to be 3x3 for 2D and 6x6 for 3D");

    matrix.Resize(dim, dim);

    ParamTools::AsTensor<double>(tens->Get("real"), dim , dim, matrix);

    // check for a scaling factor
    const double factor(tens->Get("factor")->As<Double>());
    if(factor != 1.0) matrix *= factor;
  }
  
  tens = pn->Get("isotropic", ParamNode::PASS);
  if(tens != NULL)
  {
    double emod = tens->Get("real")->Get("elasticityModulus")->As<Double>();
    double poisson = tens->Get("real")->Get("poissonNumber")->As<Double>();
    
    MechanicMaterial::CalcIsotropicStiffnessTensorFromEAndPoisson(matrix, emod, poisson, domain->GetGrid()->GetDim());
  }
  
  std::cout << "Tensor = " << matrix.ToString() << std::endl;
  
  return true;
}

void Condition::AddCondition(PtrParamNode pn, StdVector<Condition>& constraints, StdVector<Condition>& observation)
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
      g->coord.Resize(1);
      g->coord[0].first = 1;      // 1. ijkl = 1111
      g->coord[0].second = 1;
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
  // isotropy is a special constraint which blows up special tensor entry constraints
  if(g->name_ == ISOTROPY)
  {
    if(pn->Has("coord"))
      throw Exception("don't use attribute 'coord' for constraint 'isotropy'");

    if(g->type_ != EQUAL)
      throw Exception("the 'isotropy' constraint requires equality constraint type");

    // become an HOMOGENIZATION_TENSOR constraint!
    g->name_ = HOMOGENIZATION_TENSOR;

    if(domain->GetGrid()->GetDim() == 2)
    {
      // E11 - E22 = 0
      assert(g->coord.GetSize() == 0);
      g->value = 0;
      g->coord.Push_back(std::make_pair(1,1));
      g->coord.Push_back(std::make_pair(2,2));

      // E11 - E12 - 2E33 = E11 - E12 - E33 - E33 = 0
      g = g->AppendSubCondition(list);
      g->coord.Clear(); // the copy constructor above copies old stuff
      g->coord.Push_back(std::make_pair(1,1));
      g->coord.Push_back(std::make_pair(1,2));
      g->coord.Push_back(std::make_pair(3,3));
      g->coord.Push_back(std::make_pair(3,3));

      // E13 = 0
      g = g->AppendSubCondition(list);
      g->coord.Clear();
      g->coord.Push_back(std::make_pair(1,3));

      // E23 = 0
      g = g->AppendSubCondition(list);
      g->coord.Clear();
      g->coord.Push_back(std::make_pair(2,3));
    }
    else
    {
      // non-shear diagonal is constant
      // E11 = E22 = E33 -> E11 - E22 = 0, E22 - E33 = 0
      assert(g->coord.GetSize() == 0);
      g->value = 0;
      g->coord.Push_back(std::make_pair(1,1));
      g->coord.Push_back(std::make_pair(2,2));

      g = g->AppendSubCondition(list);
      g->coord.Clear();
      g->coord.Push_back(std::make_pair(2,2));
      g->coord.Push_back(std::make_pair(3,3));

      // upper non-shear triangle is constant
      // E12 = E13 = E23 -> E12 - E13 = 0, E13 - E23 = 0
      g = g->AppendSubCondition(list);
      g->coord.Clear();
      g->coord.Push_back(std::make_pair(1,2));
      g->coord.Push_back(std::make_pair(1,3));

      g = g->AppendSubCondition(list);
      g->coord.Clear();
      g->coord.Push_back(std::make_pair(1,3));
      g->coord.Push_back(std::make_pair(2,3));

      // shear diagonal is constant
      // E44 = E55 = E66 -> E44 - E55 = 0, E55 - E66 = 0
      g = g->AppendSubCondition(list);
      g->coord.Clear();
      g->coord.Push_back(std::make_pair(4,4));
      g->coord.Push_back(std::make_pair(5,5));

      g = g->AppendSubCondition(list);
      g->coord.Clear();
      g->coord.Push_back(std::make_pair(5,5));
      g->coord.Push_back(std::make_pair(6,6));

      // relationship of the three unique values
      // E11 - E12 - 2E66 = E11 - E12 - E66 - E66 = 0
      g = g->AppendSubCondition(list);
      g->coord.Clear(); // the copy constructor above copies old stuff
      g->coord.Push_back(std::make_pair(1,1));
      g->coord.Push_back(std::make_pair(1,2));
      g->coord.Push_back(std::make_pair(6,6));
      g->coord.Push_back(std::make_pair(6,6));

      // the rest is zero
      // E14 = E15 = E16 = E24 = E25 = E26 = E34 = E35 = E36 = E45 = E46 = E56 = 0
      g = g->AppendSubCondition(list);
      g->coord.Clear();
      g->coord.Push_back(std::make_pair(1,4));

      g = g->AppendSubCondition(list);
      g->coord.Clear();
      g->coord.Push_back(std::make_pair(1,5));

      g = g->AppendSubCondition(list);
      g->coord.Clear();
      g->coord.Push_back(std::make_pair(1,6));

      g = g->AppendSubCondition(list);
      g->coord.Clear();
      g->coord.Push_back(std::make_pair(2,4));

      g = g->AppendSubCondition(list);
      g->coord.Clear();
      g->coord.Push_back(std::make_pair(2,5));

      g = g->AppendSubCondition(list);
      g->coord.Clear();
      g->coord.Push_back(std::make_pair(2,6));

      g = g->AppendSubCondition(list);
      g->coord.Clear();
      g->coord.Push_back(std::make_pair(3,4));

      g = g->AppendSubCondition(list);
      g->coord.Clear();
      g->coord.Push_back(std::make_pair(3,5));

      g = g->AppendSubCondition(list);
      g->coord.Clear();
      g->coord.Push_back(std::make_pair(3,6));

      g = g->AppendSubCondition(list);
      g->coord.Clear();
      g->coord.Push_back(std::make_pair(4,5));

      g = g->AppendSubCondition(list);
      g->coord.Clear();
      g->coord.Push_back(std::make_pair(4,6));

      g = g->AppendSubCondition(list);
      g->coord.Clear();
      g->coord.Push_back(std::make_pair(5,6));
    }
  }
}

Condition* Condition::AppendSubCondition(StdVector<Condition>& list)
{
  list.Push_back(*this);
  Condition* sub = &list.Last(); // copy this entry as reference
  sub->index_ = list.GetSize() - 1;
  return sub;
}


Condition* Condition::AppendSubCondition(StdVector<Condition>& list, unsigned int pos_x, unsigned int pos_y)
{
  Condition* sub = AppendSubCondition(list);
  sub->coord.Resize(1);
  sub->coord[0].first = pos_x;
  sub->coord[0].second = pos_y;
  sub->value = sub->matrix_[pos_x-1][pos_y-1];
  return sub;
}



bool Condition::ReadCoord(PtrParamNode pn)
{
  std::string val = pn->Get("coord")->As<std::string>();
  if(val == "all") return false;

  assert(val.size() == 2);
  coord.Resize(1);
  coord[0].first  = lexical_cast<unsigned int>(val.at(0));
  coord[0].second = lexical_cast<unsigned int>(val.at(1));

  return true;
}



bool Condition::IsHomogenization() const
{
  switch(name_)
  {
    case HOMOGENIZATION_TENSOR:
    case HOMOGENIZATION_TRACKING:
    case POISSONS_RATIO:
    case YOUNGS_MODULUS:
    case ISOTROPY:
      return true;

    default:
      return false;
  }
}

std::string Condition::ToString() const
{
  std::ostringstream os;
  
  os << name.ToString(name_);
  if(design != DesignElement::DEFAULT)
    os << "_(" << DesignElement::type.ToString(design) << ")";
  if(name_ == HOMOGENIZATION_TENSOR)
    for(unsigned int i = 0; i < coord.GetSize(); i++)
      os << "_" << coord[i].first << coord[i].second;
  return os.str();  
}


void Condition::ToInfo(PtrParamNode in) const
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
  {
    std::ostringstream os;
    for(unsigned int i = 0; i < coord.GetSize(); i++)
      os << (i > 0 ? "_" : "") << coord[i].first << coord[i].second;
    in->Get("tensor_entry")->SetValue(os.str());
  }

  if(delta_logging_ignored_)
    in->Get("deltaLogging")->Get(ParamNode::WARNING)->SetValue("no value given");
  else
    in->Get("deltaLogging")->SetValue(delta_logging);
}

Matrix<double>& Condition::GetTensor()
{
  return matrix_;
}


bool Condition::IsForRegion(RegionIdType regionId)
{
  return(region == ALL_REGIONS || region == regionId);
}

