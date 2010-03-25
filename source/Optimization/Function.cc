#include "Optimization/Function.hh"
#include "Optimization/Condition.hh"
#include "Optimization/Objective.hh"
#include "General/exception.hh"
#include "General/environment.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"
#include "Materials/mechanicMaterial.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"


// instantiation of the static elements
Enum<Function::Type> Function::type;
using boost::lexical_cast;

Function::Function(PtrParamNode pn)
{
  this->pn = pn;

  this->type_ = type.Parse(pn->Get("type")->As<std::string>());

  // function value to be evaluated
  this->value_ = -1.0;

  this->physical_ = pn->Has("physical") ? pn->Get("physical")->As<bool>() : false;

  this->parameter_ = pn->Has("parameter") ? pn->Get("parameter")->As<double>() : 0.0;

  bool tensor_ok = ReadTensor(pn, this->tensor_); // is save and sets default

  if(type_ == HOMOGENIZATION_TRACKING && !tensor_ok)
   EXCEPTION("A 'tensor' element is mandatory  for 'homTracking'");

  if(type_ == HOMOGENIZATION_TENSOR || type_ == HOMOGENIZATION_TRACKING)
  {
    // we must not give a value when there is a tensor
    if(type_ == HOMOGENIZATION_TENSOR && pn->Has("tensor") && pn->Has("value"))
      throw Exception("a value must not be given when a tensor is used in a homogenization constraint");

    if(type_ == HOMOGENIZATION_TRACKING && (!pn->Has("tensor") && !pn->Has("isotropic")))
      throw Exception("a 'tensor' is mandatory for homogenization tracking");
  }

  if(pn->Has("parameter") && type_ != PENALIZED_VOLUME && type_ != GAP)
   throw Exception("function '" + type.ToString(type_) + "' does not accept 'parameter' attribute");

  if(!pn->Has("parameter") && (type_ == PENALIZED_VOLUME || type_ == GAP))
   throw Exception("function '" + type.ToString(type_) + "' requires the 'parameter' attribute");

  if(physical_ && !(type_ == VOLUME || type_ == GREYNESS))
    throw Exception("'physical' not defined for '" + type.ToString(type_) + "'");

  // some functions need to be evaluated only once (last) for multiple excitations
  // multiple excitations are:
  // * static load cases
  // * different frequencies
  // * several homogenization test strains
  // * time steps
  switch(type_)
  {
    case VOLUME:
    case PENALIZED_VOLUME:
    case GAP:
    case REALVOLUME:
    case TYCHONOFF:
    case GREYNESS:
    case HOMOGENIZATION_TENSOR:
    case HOMOGENIZATION_TRACKING:
    case POISSONS_RATIO:
    case YOUNGS_MODULUS:
    case SLOPE:
    case ISOTROPY:
    case CHECKERBOARD:
      this->evaluateOnce_ = true;
      break;

    case MULTI_OBJECTIVE: // only to make the switch complete
    case COMPLIANCE:
    case OUTPUT:
    case DYNAMIC_OUTPUT:
    case ACOU_NEAR_FIELD:
    case TRACKING:
    case ABS_DYN_OUTPUT_SQUARED:
    case CONJUGATE_COMPLIANCE:
    case GLOBAL_DYNAMIC_COMPLIANCE:
    case ELEC_ENERGY:
    case TEMPERATURE:
      this->evaluateOnce_ = false; // standard case
  }

}

bool Function::ReadTensor(PtrParamNode pn, Matrix<double>& matrix)
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

    int dim = tens->Get("dim1")->As<int>();
    if(dim != 3 && dim != 6)
      EXCEPTION("The Voigt 'tensor' for homogenizations needs to be 3x3 or 6x6");
    if(tens->Has("dim2") && dim != tens->Get("dim2")->As<int>())
      EXCEPTION("The 'tensor' for homogenization needs to be symmetric");
    if((domain->GetGrid()->GetDim() == 2 && dim != 3) || (domain->GetGrid()->GetDim() == 3 && dim != 6))
      EXCEPTION("The 'tensor' for homogenization needs to be 3x3 for 2D and 6x6 for 3D");

    matrix.Resize(dim, dim);

    ParamTools::AsTensor<double>(tens->Get("real"), dim , dim, matrix);

    // check for a scaling factor
    const double factor(tens->Get("factor")->As<double>());
    if(factor != 1.0) matrix *= factor;
  }

  tens = pn->Get("isotropic", ParamNode::PASS);
  if(tens != NULL)
  {
    double emod = tens->Get("real")->Get("elasticityModulus")->As<double>();
    double poisson = tens->Get("real")->Get("poissonNumber")->As<double>();

    MechanicMaterial::CalcIsotropicStiffnessTensorFromEAndPoisson(matrix, emod, poisson, domain->GetGrid()->GetDim());
  }

  return true;
}

void Function::ParseCoord(PtrParamNode pn, std::pair<int, int>& coord)
{
  std::string val = pn->Get("coord")->As<std::string>();
  coord.first  = lexical_cast<unsigned int>(val.at(0));
  coord.second = lexical_cast<unsigned int>(val.at(1));
}

std::string Function::ToString() const
{
  return physical_ ? "physical_" + type.ToString(type_) : type.ToString(type_);
}


Function* Function::GetFunction(Objective* f, Condition* g)
{
  assert(!(f != NULL && g != NULL) || (f == NULL && g == NULL));
  return f != NULL ? dynamic_cast<Function*>(f) : dynamic_cast<Function*>(g);
}


bool Function::IsHomogenization() const
{
  switch(type_)
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
