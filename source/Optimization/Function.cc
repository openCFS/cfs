#include "Optimization/Function.hh"
#include "Optimization/Condition.hh"
#include "Optimization/Objective.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "General/exception.hh"
#include "General/environment.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"
#include "Materials/mechanicMaterial.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/Logging/cfslog.hh"

DECLARE_LOG(ofunc)
DEFINE_LOG(ofunc, "opt_func")

// instantiation of the static elements is in Optimization::SetEnums()
Enum<Function::Type> Function::type;
Enum<Function::Locality> Function::locality;

using boost::lexical_cast;

Function::Function(PtrParamNode pn)
{
  this->harmonic_    = BasePDE::IsComplex(domain->GetDriver()->GetAnalysisType());

  this->pn = pn;
  this->local = NULL;

  this->type_ = type.Parse(pn->Get("type")->As<std::string>());
  this->locality_ = locality.Parse(pn->Get("local")->As<std::string>());

  // function value to be evaluated
  this->value_ = -1.0;

  this->physical_ = pn->Has("physical") ? pn->Get("physical")->As<bool>() : false;

  this->parameter_ = pn->Has("parameter") ? pn->Get("parameter")->As<double>() : 0.0;

  this->omega_omega_ = pn->Has("factor") ? pn->Get("factor/omega_omega")->As<bool>() : false;
  if(!harmonic_ && omega_omega_)
    throw Exception("It makes no sense to set costFunction/factor/omega_omega in static optimization");

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

  if(pn->Has("parameter") && type_ != PENALIZED_VOLUME && type_ != GAP && type_ != CHECKERBOARD)
   throw Exception("function '" + type.ToString(type_) + "' does not accept 'parameter' attribute");

  if(!pn->Has("parameter") && (type_ == PENALIZED_VOLUME || type_ == GAP || type_ == CHECKERBOARD))
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
    case GLOBAL_SLOPE:
    case ISOTROPY:
    case CHECKERBOARD:
      this->evaluateOnce_ = true;
      break;

    case MULTI_OBJECTIVE: // only to make the switch complete
    case COMPLIANCE:
    case OUTPUT:
    case DYNAMIC_OUTPUT:
    case ENERGY_FLUX:
    case TRACKING:
    case ABS_DYN_OUTPUT_SQUARED:
    case CONJUGATE_COMPLIANCE:
    case GLOBAL_DYNAMIC_COMPLIANCE:
    case ELEC_ENERGY:
    case TEMPERATURE:
      this->evaluateOnce_ = false; // standard case
  }
}

Function::~Function()
{
  if(local != NULL) { delete local; local = NULL; }
}

bool Function::ReadTensor(PtrParamNode pn, Matrix<double>& matrix)
{
  matrix.Resize(1,1); // minimal size, as 0,0 is not defined.

  // sanity checks
  if(!pn->Has("tensor") && !pn->Has("isotropic")) 
    return false;

  if(pn->Has("tensor") && pn->Has("isotropic"))
    EXCEPTION("please specify either <tensor> or <isotropic>, not both");
  
  bool tensor_read(false);

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

    tensor_read = true;
  }

  tens = pn->Get("isotropic", ParamNode::PASS);
  if(tens != NULL)
  {
    double emod = tens->Get("real")->Get("elasticityModulus")->As<double>();
    double poisson = tens->Get("real")->Get("poissonNumber")->As<double>();

    MechanicMaterial::CalcIsotropicStiffnessTensorFromEAndPoisson(matrix, emod, poisson, domain->GetGrid()->GetDim());

    tensor_read = true;
  }
  
  if(tensor_read)
  {
    // output the target tensor to xml-file
    // to get a reference and be able to do quick checks
    PtrParamNode in_mat = info->Get("target_tensor");
    in_mat->SetType(ParamNode::ELEMENT);
    in_mat->Get("tensor")->SetValue(matrix);
  }
  
  return tensor_read;
}

void Function::ParseCoord(PtrParamNode pn, std::pair<int, int>& coord)
{
  std::string val = pn->Get("coord")->As<std::string>();
  coord.first  = lexical_cast<unsigned int>(val.at(0));
  coord.second = lexical_cast<unsigned int>(val.at(1));
}

void Function::ToInfo(PtrParamNode info)
{
  info_ = info;
  info->Get("type")->SetValue(type.ToString(type_));
  if(harmonic_)
    info->Get("omega_omega")->SetValue(omega_omega_);

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

bool Function::ForDensityFiltering() const
{
  switch(type_)
  {
  case MULTI_OBJECTIVE:
    EXCEPTION("Invalid query: " << type.ToString(type_));

  default:
    return true; // actually true for almost all!
  }
}

bool Function::ForSensitivityFiltering() const
{
  switch(type_)
  {
  // pure objective
  case OUTPUT:
  case DYNAMIC_OUTPUT:
  case ABS_DYN_OUTPUT_SQUARED:
  case CONJUGATE_COMPLIANCE:
  case GLOBAL_DYNAMIC_COMPLIANCE:
  case ELEC_ENERGY:
  case ENERGY_FLUX:
  // objective and constraint
  case COMPLIANCE:
  case TRACKING:
  case HOMOGENIZATION_TENSOR:
  case HOMOGENIZATION_TRACKING:
  case POISSONS_RATIO:
  case YOUNGS_MODULUS:
  case TEMPERATURE:
    return true;

  case VOLUME:
  case PENALIZED_VOLUME:
  case GAP:
  case TYCHONOFF:
  case GREYNESS:
  case REALVOLUME:
  case SLOPE:
  case GLOBAL_SLOPE:
  case CHECKERBOARD:
    return false;

  case ISOTROPY:
  case MULTI_OBJECTIVE:
    EXCEPTION("Invalid query: " << type.ToString(type_));
  }

  EXCEPTION("can never reach! Stupid C++");
}


void Function::PostProc(DesignSpace* space, DesignStructure* structure)
{
  // pre-init step
  switch(type_)
  {
  case CHECKERBOARD:
  case SLOPE:
  case GLOBAL_SLOPE:
    assert(space->IsRegular()); // VicinityElements work only on a regular grid
    // the design elements require the vicinity element to be set which holds the direct
    // neighbors. Is save to call several times
    VicinityElement::Init(space, structure);
    InitLocal(space);

    if(type_ == SLOPE)
      info_->Get("active_size")->SetValue(local->values.GetSize());

    break;

  case PENALIZED_VOLUME:
    for(unsigned int i = 0; i < space->transfer.GetSize(); i++)
      if(space->transfer[i].IsPenalized())
        info_->Get(ParamNode::WARNING)->SetValue("transfer function '" + space->transfer[i].ToString() + " seems also to penalize");

  default: // do nothing
    break;
  }
}

Function::Local* Function::InitLocal(DesignSpace* space)
{
  if(local == NULL) local = new Local(this, space);
  return local;
}

Function::Local::Local(Function* func, DesignSpace* space)
{
  this->space = space;
  this->func_ = func;

  assert(func->locality_ == Function::NEXT_BIDIR || func->locality_ == Function::NEXT);

  bool both = func->locality_ == Function::NEXT;
  unsigned int  dim  = domain->GetGrid()->GetDim();

  element_dimension_ = dim * (both ? 2.0 : 1.0);

  virtual_elem_map.Reserve(element_dimension_ * space->GetNumberOfElements());

  // traverse all elements and check for full neighborhood
  space->AssertOneDesignOnly();
  int elems = space->GetNumberOfElements();
  for(int e = 0, ss = elems; e < ss; ++e)
  {
    DesignElement& de = space->data[e];

    // do we have a full neighborhood? All or none as in the original paper
    bool full = true;
    if(de.vicinity->design[VicinityElement::X_P] == NULL) full = false;
    if(de.vicinity->design[VicinityElement::Y_P] == NULL) full = false;
    if(dim == 3 && de.vicinity->design[VicinityElement::Z_P] == NULL) full = false;

    LOG_DBG2(ofunc) << "Local::Local e_num=" << de.elem->elemNum << " vicinity=" << de.vicinity->ToString() << " full=" << full;

    if(full)
    {
      // for the slope constraint bounding box no sign is necessary!
      if(!both)
      {
        virtual_elem_map.Push_back(Identifier(e, VicinityElement::X_P));
        virtual_elem_map.Push_back(Identifier(e, VicinityElement::Y_P));
        if(dim == 3)
          virtual_elem_map.Push_back(Identifier(e, VicinityElement::Z_P));
      }

      if(both)
      {
        virtual_elem_map.Push_back(Identifier(e, VicinityElement::X_P, 1));
        virtual_elem_map.Push_back(Identifier(e, VicinityElement::X_P, -1));

        virtual_elem_map.Push_back(Identifier(e, VicinityElement::Y_P, 1));
        virtual_elem_map.Push_back(Identifier(e, VicinityElement::Y_P, -1));

        if(dim == 3)
        {
          virtual_elem_map.Push_back(Identifier(e, VicinityElement::Z_P, 1));
          virtual_elem_map.Push_back(Identifier(e, VicinityElement::Z_P, -1));
        }
      }
    }
  }

  // needs to be set prior CalcSlopeConstraint() as the optimizers need the size
  values.Resize(virtual_elem_map.GetSize(), -1.0);
}
