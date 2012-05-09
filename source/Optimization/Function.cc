#include <algorithm>
#include <cmath>
#include <ostream>

#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "Domain/domain.hh"
#include "Domain/elem.hh"
#include "Domain/grid.hh"
#include "Driver/basedriver.hh"
#include "Elements/basefe.hh"
#include "General/Enum.hh"
#include "General/defs.hh"
#include "General/environment.hh"
#include "General/exception.hh"
#include "MatVec/matrix.hh"
#include "MatVec/vector.hh"
#include "Materials/mechanicMaterial.hh"
//#include "Optimization/Function.hh"
#include "Optimization/Condition.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Design/DesignStructure.hh"
#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/Excitation.hh"
#include "Optimization/Function.hh"
#include "Optimization/Objective.hh"
#include "Optimization/TransferFunction.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/basePDE.hh"
#include "Utils/StdVector.hh"
#include "Utils/tools.hh"

DECLARE_LOG(func)
DEFINE_LOG(func, "opt_func")

// instantiation of the static elements is in Optimization::SetEnums()
Enum<Function::Type> Function::type;
Enum<Function::StressType>  Function::stressType;
Enum<Function::Local::Locality> Function::Local::locality;
Enum<Function::Local::Phase> Function::Local::phase;

// speed up by sharing
StdVector<double> Function::Local::Identifier::tmp1;
StdVector<double> Function::Local::Identifier::tmp2;

// sync the values with Local::Phase
const int Function::Local::Identifier::NO_SIGN       = -1000;
const int Function::Local::Identifier::VOID_SIGN     = -1;
const int Function::Local::Identifier::MATERIAL_SIGN = 1;

using boost::lexical_cast;

Function::Function(PtrParamNode pn)
{
  Init();

  this->preInfo_ = PtrParamNode(new ParamNode(ParamNode::INSERT, ParamNode::ELEMENT ));
  this->pn = pn;

  this->type_ = type.Parse(pn->Get("type")->As<std::string>());

  this->physical_ = pn->Has("physical") ? pn->Get("physical")->As<bool>() : false;

  this->parameter_ = pn->Has("parameter") ? pn->Get("parameter")->As<double>() : 0.0;

  this->omega_omega_ = pn->Has("factor") ? pn->Get("factor/omega_omega")->As<bool>() : false;
  if(!harmonic_ && omega_omega_)
    throw Exception("It makes no sense to set costFunction/factor/omega_omega in static optimization");

  bool tensor_ok = ReadTensor(pn, this->tensor_); // is save and sets default

  if((type_ == HOM_TRACKING || type_ == HOM_FROBENIUS_PRODUCT) && !tensor_ok)
    EXCEPTION("A 'tensor' element is mandatory  for 'homTracking'");

  if(type_ == HOM_TENSOR || type_ == HOM_TRACKING)
  {
    // we must not give a value when there is a tensor
    if(type_ == HOM_TENSOR && pn->Has("tensor") && pn->Has("value"))
      throw Exception("a value must not be given when a tensor is used in a homogenization constraint");
  }

  if(type_ == HOM_TRACKING && (!pn->Has("tensor") && !pn->Has("isotropic")))
    throw Exception("a 'tensor' is mandatory for homogenization tracking");

  HasSelectionTensor_ = false;
  HasSelectionTensor_ = ReadMaxwellTensor(pn, this->selectionTensor_, true);
  bool maxwell_tensor_ok = ReadMaxwellTensor(pn, this->maxwellTensor_); // is save and sets default
  if (HasSelectionTensor()){
    maxwellTensor_.SetPart(Global::REAL, maxwellTensor_.GetPart(Global::REAL).EntryMult(GetSelectionTensor().GetPart(Global::REAL)));
    maxwellTensor_.SetPart(Global::IMAG, maxwellTensor_.GetPart(Global::IMAG).EntryMult(GetSelectionTensor().GetPart(Global::IMAG)));
  }

  if(type_ == MAXWELL_HOM_TRACKING && !maxwell_tensor_ok)
    EXCEPTION("A 'maxwellTensor' element is mandatory  for 'maxwellHomTracking'");

  if(type_ == MAXWELL_HOM_TENSOR || type_ == MAXWELL_HOM_TRACKING)
  {
    // we must not give a value when there is a tensor
//    if(type_ == MAXWELL_HOM_TENSOR && pn->Has("maxwellTensor") && pn->Has("value"))
//      throw Exception("a value must not be given when a tensor is used in a homogenization constraint");

    if(type_ == MAXWELL_HOM_TRACKING && (!pn->Has("maxwellTensor") && !pn->Has("isotropic")))
      throw Exception("a 'maxwellTensor' is mandatory for homogenization tracking");
  }

  // check parameter
  switch(type_)
  {
  case PENALIZED_VOLUME:
  case GAP:
  case GLOBAL_SLOPE:
  case GLOBAL_OSCILLATION:
  case GLOBAL_MOLE:
  case STRESS:
  case STRESS_DENSITY:
    if(!pn->Has("parameter"))
      throw Exception("function '" + type.ToString(type_) + "' requires the 'parameter' attribute");
    break;

  case PROJECTION:
    if(!pn->Has("filter"))
      throw Exception("function '" + type.ToString(type_) + "' requires the 'filter' element");
    if(region != ALL_REGIONS)
      throw Exception("function '" + type.ToString(type_) + "' cannot be region restricted");
    break;


  default:
    break;
  }

  if(physical_ && !(type_ == VOLUME || type_ == GREYNESS))
    throw Exception("'physical' is no option for '" + type.ToString(type_) + "'");

}

Function::~Function()
{
  if(local != NULL) { delete local; local = NULL; }
  if(projectionDesign_ != NULL) { delete projectionDesign_; projectionDesign_ = NULL; }
}

void Function::Init()
{
  this->harmonic_    = BasePDE::IsComplex(domain->GetDriver()->GetAnalysisType());
  this->design = DesignElement::DEFAULT; // overwritten eventually in Condition
  this->region = ALL_REGIONS;  // overwritten eventually in Condition

  this->local = NULL;
  this->projectionDesign_ = NULL;

    // function value to be evaluated
  this->value_ = -1.0;

  // -2 is unset, -1 is all, >= 0 the excitation index
  this->excite_ = -1;

  this->stressType_ = MECH; // set in Condition

  this->omega_omega_ = false;
  this->index_ = -1;


}

Function* Function::Cast(Objective* c, Condition* g)
{
  assert((c != NULL && g == NULL) || (c == NULL && g != NULL));
  assert((c != NULL && dynamic_cast<Function*>(c) != NULL) || (g != NULL && dynamic_cast<Function*>(g) != NULL));
  return c != NULL ? static_cast<Function*>(c) : static_cast<Function*>(g);
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

    Matrix<double> tmp(6,6); // always 3D first
    MechanicMaterial::CalcIsotropicStiffnessTensorFromEAndPoisson(tmp, emod, poisson);
    MechanicMaterial::ComputeSubTensor(matrix, domain->GetSinglePDE("mechanic")->GetSubTensorType(), tmp);

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

bool Function::ReadMaxwellTensor(PtrParamNode pn, Matrix<Complex>& matrix, bool selectionTensor)
{
  matrix.Resize(1,1); // minimal size, as 0,0 is not defined.

  bool tensor_read(false);
  PtrParamNode tens;
  if (!selectionTensor){
    // sanity checks
    if(!pn->Has("maxwellTensor") && !pn->Has("isotropic"))
      return false;

    if(pn->Has("maxwellTensor") && pn->Has("isotropic"))
      EXCEPTION("please specify either <maxwellTensor> or <isotropic>, not both");

//    tens = pn->Get("isotropic", ParamNode::PASS);
//    if(tens != NULL)
//    {
//      UInt dim = domain->GetDim();
//      double real = tens->Get("permReal")->As<double>();
//      double imag = tens->Get("permImag")->As<double>();
//      matrix.Resize(dim, dim);
//      matrix.Init();
//      for (UInt i=0; i<dim; ++i){
//        matrix(i, i) = Complex(real, imag);
//      }
//      tensor_read = true;
//    }

  // check for tensor element
    tens = pn->Get("maxwellTensor", ParamNode::PASS);
  }
  else
    tens = pn->Get("selectionTensor", ParamNode::PASS);
  if(tens != NULL)
  {
    int dim = tens->Get("dim1")->As<int>();
    if(dim != 2 && dim != 3)
      EXCEPTION("The 'tensor' for Maxwell homogenization needs to be 2x2 or 3x3");
    if(tens->Has("dim2") && dim != tens->Get("dim2")->As<int>())
      EXCEPTION("The 'tensor' for homogenization needs to be symmetric");

    matrix.Resize(dim, dim);
    matrix.Init();
    Matrix<double> tmp_mat(dim, dim);
    ParamTools::AsTensor<double>(tens->Get("real"), dim , dim, tmp_mat);
    matrix.SetPart(Global::REAL, tmp_mat);
    if (tens->Has("imag")){
      ParamTools::AsTensor<double>(tens->Get("imag"), dim , dim, tmp_mat);
      matrix.SetPart(Global::IMAG, tmp_mat);
    }

    // check for a scaling factor
    const double factor(tens->Get("factor")->As<double>());
    if(factor != 1.0) matrix *= factor;

    tensor_read = true;
  }


  if(tensor_read)
  {
    if (!selectionTensor){
      // output the target tensor to xml-file
      // to get a reference and be able to do quick checks
      PtrParamNode in_mat = info->Get("target_tensor");
      in_mat->SetType(ParamNode::ELEMENT);
      in_mat->Get("maxwellTensor")->SetValue(matrix);
    }
    else {
      // output the selection tensor to xml-file
      // to get a reference and be able to do quick checks
      PtrParamNode in_mat = info->Get("selection_tensor");
      in_mat->SetType(ParamNode::ELEMENT);
      in_mat->Get("selectionTensor")->SetValue(matrix);
    }
  }

  return tensor_read;
}

void Function::ParseCoord(PtrParamNode pn, tuple<int, int, double>& coord)
{
  std::string val = pn->Get("coord")->As<std::string>();
  get<0>(coord) = lexical_cast<unsigned int>(val.at(0));
  get<1>(coord) = lexical_cast<unsigned int>(val.at(1));
  get<2>(coord) = 1.0; // default
}

void Function::ToInfo(PtrParamNode info)
{
  info_ = info;

  // there might be set something, i.g. in PostProc
  info_->SetValue(preInfo_, false); // don't do tricks with name

  info->Get("type")->SetValue(type.ToString(type_));
  if(harmonic_)
    info->Get("omega_omega")->SetValue(omega_omega_);
  // we check for valid ocurence of paramter in the constructor
  if(pn->Has("parameter"))
    info->Get("parameter")->SetValue(parameter_);

  // We might have non-standard stresses
  if(type_ == STRESS || type_ == STRESS_DENSITY)
    info->Get("stress")->SetValue(stressType.ToString(stressType_));

  if(local != NULL)
    local->ToInfo(info_);
}

std::string Function::ToString(MultipleExcitation* me) const
{
  // optional for oscillation
  if(local != NULL && local->GetPhase() != Local::BOTH)
    return Local::phase.ToString(local->GetPhase()) + "_" + type.ToString(type_);

  if(physical_)
    return "physical_" + type.ToString(type_);

  return type.ToString(type_);
}


Function* Function::GetFunction(Objective* f, Condition* g)
{
  assert(!(f != NULL && g != NULL) || (f == NULL && g == NULL));
  return f != NULL ? dynamic_cast<Function*>(f) : dynamic_cast<Function*>(g);
}


void Function::SetExcitation(MultipleExcitation* me, int excite_index)
{
  assert(me != NULL && me->excitations.GetSize() > 0);

  // some functions need to be evaluated only once (first) for multiple excitations
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
    case HOM_TENSOR:
    case MAXWELL_HOM_TENSOR:
    case HOM_TRACKING:
    case MAXWELL_HOM_TRACKING:
    case HOM_FROBENIUS_PRODUCT:
    case BITENSOR:
    case POISSONS_RATIO:
    case YOUNGS_MODULUS:
    case YOUNGS_MODULUS_E1:
    case YOUNGS_MODULUS_E2:
    case SLOPE:
    case GLOBAL_SLOPE:
    case ISOTROPY:
    case ISO_ORTHOTROPY:
    case ORTHOTROPY:
    case MOLE:
    case GLOBAL_MOLE:
    case OSCILLATION:
    case GLOBAL_OSCILLATION:
    case JUMP:
    case GLOBAL_JUMP:
    case MAXWELL_ISOTROPY:
    case BIISOTROPY:
    case BUMP:
    case DESIGN_TRACKING:
    case PROJECTION:
    case SUM_MODULI:
    case GLOBAL_SUM_MODULI:
    case PARAM_PS_POS_DEF:
      assert(excite_index < 0);
      excite_ = me->excitations.GetSize() - 1; // once only at the last excitation
      break;

    case MULTI_OBJECTIVE: // only to make the switch complete
    case COMPLIANCE:
    case OUTPUT:
    case DYNAMIC_OUTPUT:
    case ENERGY_FLUX:
    case TRACKING:
    case ABS_OUTPUT:
    case CONJUGATE_COMPLIANCE:
    case GLOBAL_DYNAMIC_COMPLIANCE:
    case ELEC_ENERGY:
    case TEMPERATURE:
      assert(excite_index < 0);
      excite_ = -1; // all excitations
      break;

    case STRESS:
    case STRESS_DENSITY:
      // there might be the optional excitation index set
      if(pn->Get("excitation")->As<std::string>() == "all")
      {
        excite_ = excite_index == -2 ? -1 : excite_index;
      }
      else
      {
        assert(excite_index == -2); // assert there is no conflict
        excite_ = me->GetExcitation(pn->Get("excitation")->As<std::string>())->index;
      }
      break;
  }
}

/** Shall/must we evaluate this objective at this excitation?
 * Stress constraints in homogenization are triggered for a single constraint only. */
bool Function::DoEvaluate(const Excitation* excite) const
{
  if(DoEvaluateAlways())
    return true;

  return excite->index == excite_;
}

bool Function::DoEvaluateAlways() const
{
  return excite_ == -1;
}

bool Function::IsExcitationSensitive() const
{
  if(type_ == STRESS || type_ == STRESS_DENSITY)
    return true;
  else
    return false;
}

bool Function::IsAdjointBased() const
{
  switch(type_)
  {
    case COMPLIANCE: // only in the transient case
    case TRACKING:
    case OUTPUT:
    case CONJUGATE_COMPLIANCE:
    case ABS_OUTPUT:
    case GLOBAL_DYNAMIC_COMPLIANCE:
    case DYNAMIC_OUTPUT:
    case ELEC_ENERGY:
    case ENERGY_FLUX:
    case STRESS:
    case STRESS_DENSITY:
    return true;

  default:
    return false;
  }
}

bool Function::NeedsSelectionVector() const
{
  switch(type_)
  {
    case OUTPUT:
//    case CONJUGATE_COMPLIANCE: ??
    case ABS_OUTPUT:
    case DYNAMIC_OUTPUT:
    case ENERGY_FLUX:
    return true;

  default:
    return false;
  }
}


bool Function::IsHomogenization() const
{
  switch(type_)
  {
    case HOM_TENSOR:
    case HOM_TRACKING:
    case HOM_FROBENIUS_PRODUCT:
    case POISSONS_RATIO:
    case YOUNGS_MODULUS:
    case YOUNGS_MODULUS_E1:
    case YOUNGS_MODULUS_E2:
    case ISOTROPY:
    case ISO_ORTHOTROPY:
    case ORTHOTROPY:
      return true;

    default:
      return false;
  }
}


bool Function::IsMaxwellHomogenization() const
{
  switch(type_)
  {
    case MAXWELL_HOM_TENSOR:
    case MAXWELL_HOM_TRACKING:
    case BITENSOR:
    case MAXWELL_ISOTROPY:
    case BIISOTROPY:
      return true;

    default:
      return false;
  }
}


bool Function::ForDensityFiltering() const
{
  switch(type_)
  {
  case PROJECTION:
    // for the projection case we have a density filter manually on Function::projectionDesign only
  case SUM_MODULI:
    return false;

  case MULTI_OBJECTIVE:
    EXCEPTION("Invalid query: " << type.ToString(type_));
    break;

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
  case CONJUGATE_COMPLIANCE:
  case GLOBAL_DYNAMIC_COMPLIANCE:
  case ELEC_ENERGY:
  case ENERGY_FLUX:
  // objective and constraint
  case COMPLIANCE:
  case TRACKING:
  case HOM_TENSOR:
  case MAXWELL_HOM_TENSOR:
  case HOM_TRACKING:
  case MAXWELL_HOM_TRACKING:
  case HOM_FROBENIUS_PRODUCT:
  case BITENSOR:
  case POISSONS_RATIO:
  case YOUNGS_MODULUS:
  case YOUNGS_MODULUS_E1:
  case YOUNGS_MODULUS_E2:
  case TEMPERATURE:
  case ABS_OUTPUT:
  case STRESS:
  case STRESS_DENSITY:
    return true;

  case VOLUME:
  case PENALIZED_VOLUME:
  case GAP:
  case TYCHONOFF:
  case GREYNESS:
  case REALVOLUME:
  case SLOPE:
  case GLOBAL_SLOPE:
  case MOLE:
  case GLOBAL_MOLE:
  case OSCILLATION:
  case GLOBAL_OSCILLATION:
  case JUMP:
  case GLOBAL_JUMP:
  case BUMP:
  case DESIGN_TRACKING:
  case PROJECTION:
  case SUM_MODULI:
  case GLOBAL_SUM_MODULI:
  case PARAM_PS_POS_DEF:
    return false;

  case ISOTROPY:
  case ISO_ORTHOTROPY:
  case ORTHOTROPY:
  case MULTI_OBJECTIVE:
  case MAXWELL_ISOTROPY:
  case BIISOTROPY:
    EXCEPTION("Invalid query: " << type.ToString(type_));
    break;
  }

  EXCEPTION("can never reach! Stupid C++");
}

StdVector<DesignElement>& Function::GetProjectionDesignClone()
{
  assert(type_ == PROJECTION);
  return projectionDesign_->data;
}

void Function::SetElements(DesignSpace* space, RegionIdType region)
{
  assert(elements.GetSize() == 0);
  Grid* grid = domain->GetGrid();

  // Bastian's multiple design test cases have situations where design is DEFAULT as it is not
  // set in the objective

  // if ALL_REGIONS for condition use what we define as design space which
  // this is still not good enough
  int nd = 1;
  if(design == DesignElement::TENSOR_TRACE || design == DesignElement::DEFAULT) {
    nd = space->design.GetSize();
    if(design == DesignElement::ALL_DESIGNS)
      nd++;
  }
  elements.Reserve(nd * (region == ALL_REGIONS ? space->GetNumberOfElements() : grid->GetNumElems(region)));

  if(region == ALL_REGIONS || space->Contains(region))
  {
    for(unsigned int i = 0; i < space->data.GetSize(); i++)
    {
      DesignElement* de = &(space->data[i]);
      if((design == DesignElement::DEFAULT || design == DesignElement::TENSOR_TRACE || design == de->GetType())
          && (region == ALL_REGIONS || de->elem->regionId == region))
        elements.Push_back(de);
    }
    if(design == DesignElement::ALL_DESIGNS)
    {
      for(unsigned int i = 0; i < space->GetNumberOfElements(); i++)
          {
            DesignElement* de = &(space->data[i]);
            elements.Push_back(de);
          }
    }
  }
  else
  {
    // this is a special case where the constraint does not act on the design space
    if(type_ != STRESS && type_ != STRESS_DENSITY)
    {
      std::string msg = "region " + grid->GetRegion().ToString(region) + " of condition " + type.ToString(type_) + " not within design domain";
      info_->Get(ParamNode::WARNING)->SetValue(msg);
    }

    assert(elements.GetSize() == 0);

    // this creates the pseudo design elements and both indices are hopefully properly set!
    space->RegisterPseudoDesignRegion(region, design, &elements);
  }

//  assert(elements.GetSize() == elements.Capacity());
}


void Function::PostProc(DesignSpace* space, DesignStructure* structure, ErsatzMaterial* em)
{
  // pre-init step
  switch(type_)
  {
  case SLOPE:
  case GLOBAL_SLOPE:
  case MOLE:
  case GLOBAL_MOLE:
  case OSCILLATION:
  case GLOBAL_OSCILLATION:
  case JUMP:
  case GLOBAL_JUMP:
  case BUMP:
  case STRESS:
  case STRESS_DENSITY:
  case SUM_MODULI:
  case GLOBAL_SUM_MODULI:
  case PARAM_PS_POS_DEF:
    // assert(space->IsRegular()); // VicinityElements work only on a regular grid
    // the design elements require the vicinity element to be set which holds the direct
    // neighbors. Is save to call several times
    VicinityElement::Init(space, structure);
    InitLocal(space);
    break;

  case PENALIZED_VOLUME:
    for(unsigned int i = 0; i < space->transfer.GetSize(); i++)
      if(space->transfer[i].IsPenalized())
        preInfo_->Get(ParamNode::WARNING)->SetValue("transfer function '" + space->transfer[i].ToString() + " seems also to penalize");
    break;

  case PROJECTION:
  {
    // We have to create a deep copy of the original design space. Also the neighborhood must not point to the original design.
    projectionDesign_ = space->Clone(); // the original regions if all

    DesignStructure ds(projectionDesign_, projectionDesign_->GetRegionIds());
    VicinityElement::Init(projectionDesign_, &ds);

    StdVector<DesignElement>& fake_data = GetProjectionDesignClone();

    PtrParamNode  reg = pn->Get("filter");
    ds.SetFilters(reg, preInfo_, &fake_data);

    assert(space->data.GetSize() == space->GetTotalElements().GetSize());
    assert(space->data.GetSize() == fake_data.GetSize());

    for(unsigned int i = 0, n = fake_data.GetSize(); i < n; i++)
    {
      // projectionDesign_[i].simp = new SIMPElement(&projectionDesign_[i]);
      // we need the gradient size for temporary storage in the gradient calculation
      fake_data[i].PostInit(em->objectives.data.GetSize(), em->constraints.active.GetSize());
    }
    break;
  }


  default: // do nothing
    break;
  }

  // don't define the elements here, it is specific for objective and conditions
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
  this->structure_ = NULL;
  this->infeasible = 0;

  // shortcuts
  Function::Type ftype = func->GetType();
  std::string    fname = Function::type.ToString(ftype);

  // read xml parameters -> might be null valued!
  PtrParamNode pn = func->pn->Get("local", ParamNode::PASS);

  this->beta_  = pn != NULL && pn->Has("beta") ? pn->Get("beta")->As<double>() : -3.14;
  this->eps_   = pn != NULL && pn->Has("eps") ? pn->Get("eps")->As<double>() : -3.14;
  this->power_ = pn != NULL && pn->Has("power") ? pn->Get("power")->As<double>() : 2.0;
  this->phase_ = pn != NULL && pn->Has("phase") ? phase.Parse(pn->Get("phase")->As<std::string>()) : BOTH; // only oscillation

  this->normalize_ = pn != NULL ? pn->Get("normalize")->As<bool>() : false;

  switch(ftype)
  {
  case GLOBAL_JUMP:
  case GLOBAL_MOLE:
  case GLOBAL_OSCILLATION:
  case GLOBAL_SLOPE:
  case STRESS:
  case STRESS_DENSITY:
  case GLOBAL_SUM_MODULI:
    this->globalized_ = true;
    break;

  default:
    this->globalized_ = false;
    break;
  }

  // check beta
  if(ftype == OSCILLATION || ftype == GLOBAL_OSCILLATION)
  {
    if(pn == NULL || !pn->Has("beta"))
      throw Exception("function '" + fname + "' requires the 'beta' attribute in a 'local' element");
    if((func->IsObjective() || dynamic_cast<Condition*>(func)->IsActive()) && beta_ < 0)
      throw Exception("'function '" + fname + "' allows beta=-1 only for condition in observe mode");
  }

  // check eps
  if((ftype == MOLE || ftype == GLOBAL_MOLE) && (pn == NULL || !pn->Has("eps")))
    throw Exception("function '" + fname + "' requires the 'eps' attribute in a 'local' element");

  // check phase
  if(phase_ != BOTH && ftype != OSCILLATION && ftype != GLOBAL_OSCILLATION)
    throw Exception("'phase' may only be set for (global) oscillation");

  // set locality
  this->locality_ = pn != NULL && pn->Has("locality") ?
      locality.Parse(pn->Get("locality")->As<std::string>()) : DEFAULT;
  Locality user = locality_; // default or set by user
  bool snopt = param->Get("optimization/optimizer/type")->As<std::string>() == "snopt";

  switch(ftype)
  {
  case SLOPE:
    if(user == DEFAULT && snopt)  locality_ = NEXT;
    if(user == DEFAULT && !snopt) locality_ = NEXT_AND_REVERSE;
    if(!snopt && locality_ != NEXT_AND_REVERSE)
      throw Exception("The optimizer has no bounds for constraints: your choice for 'local' is invalid");
    if(locality_ != NEXT && locality_ != NEXT_AND_REVERSE)
      throw Exception("Invalid locality '" + locality.ToString(locality_) + "' within '" + fname + "'");
    break;

  case GLOBAL_SLOPE:
    if(locality_ != NEXT && locality_ != DEFAULT)
      throw Exception("Invalid locality '" + locality.ToString(locality_) + "' within '" + fname + "'");
    locality_ = NEXT_AND_REVERSE;
    break;

  case OSCILLATION:
  case GLOBAL_OSCILLATION:
    if((phase_ == BOTH && locality_ != DEG_45_STAR_AND_REVERSE && locality_ !=  PREV_NEXT_AND_REVERSE && locality_ != DEFAULT)
      || (phase_ != BOTH && locality_ != DEG_45_STAR && locality_ !=  PREV_NEXT && locality_ != DEFAULT) )
      throw Exception("Invalid locality '" + locality.ToString(locality_) + "' within '" + fname + "' and phase '" + phase.ToString(phase_) + "'");
    if(locality_ == DEFAULT) locality_ = phase_ == BOTH ? DEG_45_STAR_AND_REVERSE : DEG_45_STAR;
    break;

  case MOLE:
  case GLOBAL_MOLE:
    if(locality_ != DEG_45_STAR && locality_ != DEFAULT)
      throw Exception("Invalid locality '" + locality.ToString(locality_) + "' within '" + fname + "'");
    locality_ = DEG_45_STAR;
    break;

  case JUMP:
  case GLOBAL_JUMP:
    if(locality_ != BOUNDARY && locality_ != DEFAULT)
      throw Exception("Invalid locality '" + locality.ToString(locality_) + "' within '" + fname + "'");
    locality_ = BOUNDARY;
    break;

  case STRESS:
  case STRESS_DENSITY:
    if(locality_ != ELEMENT && locality_ != DEFAULT)
      throw Exception("Invalid locality '" + locality.ToString(locality_) + "' within '" + fname + "'");
    locality_ = ELEMENT;
    break;

  case SUM_MODULI:
  case GLOBAL_SUM_MODULI:
  case PARAM_PS_POS_DEF:
    if(locality_ != MULT_DESIGNS_ELEMENT && locality_ != DEFAULT)
      throw Exception("Invalid locality '" + locality.ToString(locality_) + "' within '" + fname + "'");
    locality_ = MULT_DESIGNS_ELEMENT;
    break;

  case BUMP:
    if(locality_ != PREV_NEXT && locality_ != DEFAULT)
      throw Exception("Invalid locality '" + locality.ToString(locality_) + "' within '" + fname + "'");
    locality_ = PREV_NEXT;
    break;

  default: // no locality
    assert(false);
    break;
  }

  // this is actually pure constructor work, just extracted to handle function size
  switch(locality_)
  {
  case DEG_45_STAR:
  case DEG_45_STAR_AND_REVERSE:
  case BOUNDARY:
    if(!pn)
      throw Exception("sub element 'local' with neighborhood information mandatory for '" + fname + "'");
    structure_ = new NeighborhoodStructure(this, pn);
    if(locality_ == BOUNDARY)
      SetupBoundaryElementMap();
    else
      SetupStarLocalityElementMap(phase_);
    break;

  case ELEMENT:
    SetupSingularElementMap();
    break;

  case MULT_DESIGNS_ELEMENT:
    SetupMultDesignsElementMap();
    break;

  default:
    SetupVirtualElementMap(phase_);
    break;
  }

  if(virtual_elem_map.GetSize() == 0) throw Exception("mesh too small for locality of function '" + fname + "'");

  // needs to be set prior CalcSlopeConstraint() as the optimizers need the size
  values.Resize(virtual_elem_map.GetSize(), -1.0);
}

Function::Local::~Local()
{
  if(structure_ != NULL) { delete structure_; structure_ = NULL; }
}

void Function::Local::SetupVirtualElementMap(Phase ph)
{
  // we construct locality_ into reverse, prev and next
  // reverse means we have a REVERSE option which makes two constraints with different signs
  int  dim     = domain->GetGrid()->GetDim();
  bool prev    = locality_ == PREV_NEXT_AND_REVERSE || locality_ == PREV_NEXT;
  bool next    = true; // always
  bool two_signs = locality_ == NEXT_AND_REVERSE || locality_ == PREV_NEXT_AND_REVERSE;
  // assert((ph == BOTH && two_signs) || (!two_signs && ph != BOTH));
  // assume ph is set correctly and Phase is in sync with the signs
  int sign_1 = ph != BOTH ? (int) ph : two_signs ? 1 : Identifier::NO_SIGN;
  int sign_2 = ph != BOTH ? (int) ph : -1;

  element_dimension_ = dim * (two_signs ? 2 : 1);

  virtual_elem_map.Reserve(element_dimension_ * space->GetNumberOfElements());


  // traverse all elements and check for full neighborhood
  for(int e = 0, ss = space->data.GetSize(); e < ss; ++e)
  {
    DesignElement* de = &(space->data[e]);
    if(de->GetType() == ( (func_->design == DesignElement::DEFAULT) ? DesignElement::DENSITY : func_->design)){
      VicinityElement* ve = de->vicinity;

      // do we have a full neighborhood? All or none as in the original slope paper
      bool full = true;
      if(prev)
      {
        if(ve->design[VicinityElement::X_N] == NULL) full = false;
        if(ve->design[VicinityElement::Y_N] == NULL) full = false;
        if(dim == 3 && ve->design[VicinityElement::Z_N] == NULL) full = false;
      }
      if(next)
      {
        if(ve->design[VicinityElement::X_P] == NULL) full = false;
        if(ve->design[VicinityElement::Y_P] == NULL) full = false;
        if(dim == 3 && ve->design[VicinityElement::Z_P] == NULL) full = false;
      }

      LOG_DBG2(func) << "Local::Local e_num=" << de->elem->elemNum << " vicinity=" << ve->ToString() << " full=" << full;

      if(full)
      {
        for(int a = 0; a < dim; a++)
        {
          DesignElement* prev_de = prev ? ve->GetNeighbour(VicinityElement::ToNeighbour(a, -1)) : NULL;
          DesignElement* next_de = ve->GetNeighbour(VicinityElement::ToNeighbour(a, 1));

          virtual_elem_map.Push_back(Identifier(de, prev_de, next_de, sign_1));
          if(two_signs)
            virtual_elem_map.Push_back(Identifier(de, prev_de, next_de, sign_2));
        }
      }
    }
  }
}

void Function::Local::SetupStarLocalityElementMap(Phase ph)
{
  unsigned int dim  = domain->GetGrid()->GetDim();
  // oscillation has with BOTH DEG_45_STAR_AND_REVERSE,
  // mole has always BOTH and DEG_45_STAR.
  // oscillation w/o BOTH needs to be DEG_45_STAR
  Function::Type ft = func_->type_;
  assert(ft == OSCILLATION || ft == GLOBAL_OSCILLATION || ft == MOLE || ft == GLOBAL_MOLE);
  assert(locality_ == DEG_45_STAR || locality_ == DEG_45_STAR_AND_REVERSE);
  assert((ph != BOTH && locality_ == DEG_45_STAR) || ph == BOTH);
  assert(structure_ != NULL);
  NeighborhoodStructure* struc = structure_;

  // mole has NO_SIGN and no reverse
  // oscillation has either the given phase or when BOTH we have to add -1 and 1 as signs
  int sign_1 = (ft == MOLE || ft == GLOBAL_MOLE) ? Identifier::NO_SIGN : ph == BOTH ? Identifier::VOID_SIGN : ph;
  // sign_2 is relevant for DEG_45_STAR_AND_REVERSE only
  int sign_2 = Identifier::MATERIAL_SIGN; // the only possibility: oscillation with BOTH and REVERSE
  // the *and* reverse mode? and is to be read as plus
  bool two_signs = locality_ == DEG_45_STAR_AND_REVERSE;

  LOG_DBG(func) << "SSLEM: phase=" << phase.ToString(ph) << " ft=" << func_->ToString()
                << " locality=" << locality.ToString(locality_) << " s1=" << sign_1 << " s2=" << sign_2;

  element_dimension_ = (dim == 2 ? 4 : 13) * (two_signs ? 1 : 2); // see paper
  virtual_elem_map.Reserve(element_dimension_ * space->GetNumberOfElements());

  space->AssertOneDesignOnly(); // can be extended we use the design from the condition
  int elems = space->GetNumberOfElements();
  for(int e = 0, ss = elems; e < ss; ++e)
  {
    DesignElement* de = &(space->data[e]);

    // do we have a full neighborhood? All or none
    bool full = true;

    // we only need to check orthogonal
    assert(VicinityElement::X_P == 0);
    for(int dir = VicinityElement::X_P; dir <= (dim == 2 ? VicinityElement::Y_N : VicinityElement::Z_N); dir++)
    {
      int a = VicinityElement::ToMainAxis( (VicinityElement::Neighbour) dir );
      unsigned int n = struc->orthogonal[a];
      assert(n > 0);
      if(!VicinityElement::HasNeighbor(de, (VicinityElement::Neighbour) dir, n))
        full = false;
    }

    if(!full) continue;
    // orthogonal first
    StdVector<DesignElement*> buddies;

    for(int dir = VicinityElement::X_P; dir <= (dim == 2 ? VicinityElement::Y_N : VicinityElement::Z_N); dir += 2)
    {
      VicinityElement::Neighbour pos = (VicinityElement::Neighbour) dir;
      VicinityElement::Neighbour neg = (VicinityElement::Neighbour) (dir + 1);
      unsigned int a = VicinityElement::ToMainAxis(pos);
      unsigned int n = struc->orthogonal[a];

      buddies.Resize(0);
      for(unsigned int i = n; i > 0; i--)  buddies.Push_back(VicinityElement::GetNeighbour(de, neg, i));
      for(unsigned int i = 1; i <= n; i++) buddies.Push_back(VicinityElement::GetNeighbour(de, pos, i));

      LOG_DBG3(func) << "L:SSLEM: de=" << de->ToString() << " dir=" << dir << " pos=" << pos << " neg=" << neg << " a=" << a << " n=" << n << " buddies=" << DesignElement::ToString(buddies);

      virtual_elem_map.Push_back(Identifier(de, buddies, sign_1));
      if(two_signs)
        virtual_elem_map.Push_back(Identifier(de, buddies, sign_2));
    }

    // In 2D we have 2 diagonals in the xy plane. In 3D also in the xz and the yz plane which makes
    // 6 diagonals in the planes plus 4 diagonals between the 8 corners. -> 10 diagonals in 3D

    // we make the assumption, that our elements are almost quadratic, such our diagonals are
    // simply alternating the directions. The total number of diagonals (we want to describe a
    // ball and not a cube) is determined by the ball radius and the element diagonals.
    // start with plane diagonals
    for(unsigned int d = 0; d < (dim == 2 ? 1 : 3); d++)
    {
      unsigned int n = struc->diagonal[d];
      // within a plane we go zig-zag. First the first axis then the second.
      // Plus first in negative direction then positive
      // the planes are xy, xy and yz
      int axis_first  = d != 2 ? 0 : 1; // xy, xz, yz
      int axis_second = d == 0 ? 1 : 2; // xy, xz, yz

      // for the first  direction we have the X_N/Y_N elements and the X_P/Y_P elements.
      // for the second direction we have the X_N/Y_P elements and the X_P/Y_N elements.
      for(int dir = -1; dir <= 1; dir +=2)
      {
        buddies.Resize(0);
        for(int e = n; e > 0; e--)
        {
          DesignElement* tmp = VicinityElement::GetNeighbour(de, VicinityElement::ToNeighbour(axis_first, -1), e);
          buddies.Push_back(VicinityElement::GetNeighbour(tmp,VicinityElement::ToNeighbour(axis_second, dir), e));
        }
        for(unsigned int e = 1; e <= n; e++)
        {
          DesignElement* tmp = VicinityElement::GetNeighbour(de, VicinityElement::ToNeighbour(axis_first, 1), e);
          buddies.Push_back(VicinityElement::GetNeighbour(tmp, VicinityElement::ToNeighbour(axis_second, dir == 1 ? -1 : 1), e));
        }

        LOG_DBG3(func) << "L:SSLEM: diag de=" << de->ToString() << " dir=" << dir << " buddies=" << DesignElement::ToString(buddies);

        virtual_elem_map.Push_back(Identifier(de, buddies, sign_1));
        if(two_signs)
          virtual_elem_map.Push_back(Identifier(de, buddies, sign_2));
      }
    }
    if(dim == 3)
    {
      // now the four corners
      // -x -y -z -> +x +y +z
      // -x -y +z -> +x +y -z
      // -x +y -z -> +x -y +z
      // -x +y +z -> +x -y -z
      assert(struc->diagonal.GetSize() == 4);
      unsigned int n = struc->diagonal[3];
      for(int dir_y = -1; dir_y <= 1; dir_y += 2)
      {
        for(int dir_z = -1; dir_z <= 1; dir_z += 2)
        {
          buddies.Resize(0);
          for(int e = n; e > 0; e--)
          {
            DesignElement* tmp_x = VicinityElement::GetNeighbour(de, VicinityElement::X_N, e);
            DesignElement* tmp_y = VicinityElement::GetNeighbour(tmp_x, VicinityElement::ToNeighbour(1, dir_y), e);
            buddies.Push_back(VicinityElement::GetNeighbour(tmp_y, VicinityElement::ToNeighbour(2, dir_z), e));
          }
          for(unsigned int e = 1; e <= n; e++)
          {
            DesignElement* tmp_x = VicinityElement::GetNeighbour(de, VicinityElement::X_P, e);
            DesignElement* tmp_y = VicinityElement::GetNeighbour(tmp_x, VicinityElement::ToNeighbour(1, dir_y == 1 ? -1 : 1), e);
            buddies.Push_back(VicinityElement::GetNeighbour(tmp_y, VicinityElement::ToNeighbour(2, dir_z == 1 ? -1 : 1), e));
          }

          LOG_DBG3(func) << "L:SSLEM: corner de=" << de->ToString() << " dir_y=" << dir_y << " dir_z=" << dir_z << " buddies=" << DesignElement::ToString(buddies);

          virtual_elem_map.Push_back(Identifier(de, buddies, sign_1));
          if(two_signs)
            virtual_elem_map.Push_back(Identifier(de, buddies, sign_2));
        }
      }
    }
  }
}

void Function::Local::SetupBoundaryElementMap()
{
  unsigned int dim  = domain->GetGrid()->GetDim();
  // oscillation has with BOTH DEG_45_STAR_AND_REVERSE,
  // mole has always BOTH and DEG_45_STAR.
  // oscillation w/o BOTH needs to be DEG_45_STAR
  assert(structure_ != NULL);
  NeighborhoodStructure* struc = structure_;

  element_dimension_ = dim * 2; // two boundary "stones" per dimension
  virtual_elem_map.Reserve(element_dimension_ * space->GetNumberOfElements());

  space->AssertOneDesignOnly(); // can be extended we use the design from the condition
  int elems = space->GetNumberOfElements();
  for(int e = 0, ss = elems; e < ss; ++e)
  {
    DesignElement* de = &(space->data[e]);

    // do we have a full neighborhood? All or none
    bool full = true;

    for(int dir = VicinityElement::X_P; dir <= (dim == 2 ? VicinityElement::Y_N : VicinityElement::Z_N); dir++)
    {
      int a = VicinityElement::ToMainAxis( (VicinityElement::Neighbour) dir );
      unsigned int n = struc->orthogonal[a];
      assert(n > 0);
      if(!VicinityElement::HasNeighbor(de, (VicinityElement::Neighbour) dir, n))
        full = false;
    }

    if(!full) continue;

    for(int dir = VicinityElement::X_P; dir <= (dim == 2 ? VicinityElement::Y_N : VicinityElement::Z_N); dir += 2)
    {
      VicinityElement::Neighbour pos = (VicinityElement::Neighbour) dir;
      VicinityElement::Neighbour neg = (VicinityElement::Neighbour) (dir + 1);
      unsigned int a = VicinityElement::ToMainAxis(pos);
      unsigned int n = struc->orthogonal[a];

      DesignElement* prev = VicinityElement::GetNeighbour(de, neg, n);
      DesignElement* next = VicinityElement::GetNeighbour(de, pos, n);

      LOG_DBG3(func) << "L:SBEM: de=" << de->ToString() << " dir=" << dir << " pos=" << pos << " neg=" << neg << " a=" << a
                     << " n=" << n << " prev=" << prev << " next=" << next;

      virtual_elem_map.Push_back(Identifier(de, prev, next));
    }
  }
}


void Function::Local::SetupSingularElementMap()
{
  // only this element!
  element_dimension_ = 1; // two boundary "stones" per dimension
  virtual_elem_map.Reserve(element_dimension_ * func_->elements.GetSize());

  StdVector<DesignElement*> empty;

  for(int e = 0, en = func_->elements.GetSize(); e < en; e++)
  {
    DesignElement* de = func_->elements[e];
    virtual_elem_map.Push_back(Identifier(de, empty));
  }
}


void Function::Local::SetupMultDesignsElementMap()
{
  // only this element!
  element_dimension_ = 1; // two boundary "stones" per dimension
  UInt numFElements = space->GetNumberOfElements();
  virtual_elem_map.Reserve(element_dimension_ * numFElements);

  for(int e = 0, en = numFElements; e < en; e++)
  {
    DesignElement* de = func_->elements[e];
    int base = space->Find(de->elem, true);
    StdVector<DesignElement*> neighbours;
    for(unsigned int index = base + numFElements; index < space->data.GetSize(); index += numFElements){
      neighbours.Push_back(&(space->data[index]));
    }
    virtual_elem_map.Push_back(Identifier(de, neighbours));
  }
}


void Function::Local::ToInfo(PtrParamNode in)
{
  Function::Type ft = func_->type_;
  in->Get("locality")->SetValue(locality.ToString(locality_));
  in->Get("local_size")->SetValue(virtual_elem_map.GetSize());

  if(IsGlobalized())
  {
    in->Get("normalize")->SetValue(normalize_);
    in->Get("power")->SetValue(power_);
  }
  if(ft == OSCILLATION || ft == GLOBAL_OSCILLATION)
  {
    in->Get("beta")->SetValue(beta_);
    in->Get("phase")->SetValue(phase.ToString(phase_));
  }

  if(ft == MOLE || ft == GLOBAL_MOLE)
    in->Get("eps")->SetValue(eps_);

  if(structure_ != NULL)
    structure_->ToInfo(in->Get("neighborhood"));
}

Function::Local::NeighborhoodStructure::NeighborhoodStructure(Local* local, PtrParamNode pn)
{
  unsigned int dim = domain->GetGrid()->GetDim();
  // sample design element -> assume regular grid
  DesignElement& de = local->space->data[0];

  value = pn->Get("neighbor_value")->As<double>();
  fs = DesignStructure::filterSpace.Parse(pn->Get("neighbor_type")->As<std::string>());
  radius = DesignStructure::FindFilterRadius(fs, &de, value);

  // find the orthogonal dimensions based on radius
  Matrix<double>  coords;
  domain->GetGrid()->GetElemNodesCoord(coords, de.elem->connect, false );
  StdVector<double> edges;
  de.elem->ptElem->GetEdgeLength(coords, edges);

  assert(edges.GetSize() == dim);
  orthogonal.Resize(dim);
  for(unsigned int i = 0; i < edges.GetSize(); i++)
  {
    orthogonal[i] = (int) ((radius / edges[0]) + 0.5); // proper rounding
    LOG_DBG(func) << "L:NS:NS orthogonal[" << i << "]: radius=" << radius << " edge=" << edges[i] << " -> " << orthogonal[i];
  }

  // validate
  for(unsigned int i = 0; i < orthogonal.GetSize(); i++)
    if(orthogonal[i] == 0)
      throw Exception("your local neighbor radius for '" + Function::type.ToString(local->func_->GetType()) + "' is too small");

  // diagonals are xy, xz and yz plane, only xy in 2D
  for(int i = 0; i < (dim == 2 ? 1 : 3); i++)
  {
    int x = i != 2 ? 0 : 1; // xy, xz, yz
    int y = i == 0 ? 1 : 2; // xy, xz, yz
    double diag = std::sqrt(edges[x] * edges[x] + edges[y] * edges[y]);
    diagonal.Push_back((int) ((radius / diag) + 0.5));
    LOG_DBG(func) << "L:NS:NS diagonal[0]: x=" << x << " y=" << y << " radius="
                  << radius << " diag=" << diag << " -> " << diagonal[0];
  }
  // all 4 "total" diagonals have the same size
  if(dim == 3)
  {
    double diag = std::sqrt(edges[0] * edges[0] + edges[1] * edges[1] + edges[2] * edges[2]);
    diagonal.Push_back((int) ((radius / diag) + 0.5));
  }

  assert((dim == 2 && diagonal.GetSize() == 1) || (dim == 3 && diagonal.GetSize() == 4));
}

void Function::Local::NeighborhoodStructure::ToInfo(PtrParamNode in)
{
  in->Get("type")->SetValue(DesignStructure::filterSpace.ToString(fs));
  in->Get("value")->SetValue(value);
  in->Get("radius")->SetValue(radius);

  StdVector<double> tmp;
  tmp.Push_back(orthogonal[0] * 2 + 1);
  tmp.Push_back(orthogonal[1] * 2 + 1);
  if(orthogonal.GetSize() == 3)
    tmp.Push_back(orthogonal[2] * 2 + 1);

  in->Get("orthogonal")->SetValue(tmp.ToString());

  tmp.Resize(0);
  for(unsigned int i = 0; i < diagonal.GetSize(); i++)
    tmp.Push_back(diagonal[i] * 2 + 1);

  in->Get("diagonal")->SetValue(tmp.ToString());
}

Function::Local::Identifier::Identifier(DesignElement* elem, DesignElement* prev, DesignElement* next, int si)
{
  this->element = elem;

  assert(next != NULL);
  bool has_prev = prev != NULL;

  this->neighbor.Resize(has_prev ? 2 : 1);

  if(has_prev)
    this->neighbor[0] = prev;

  this->neighbor[has_prev ? 1 : 0] = next;

  this->sign = si;
}

Function::Local::Identifier::Identifier(DesignElement* elem, StdVector<DesignElement*> buddies, int si)
{
  this->element = elem;
  this->neighbor = buddies;
  assert(si == NO_SIGN || si == -1 || si == 1);
  this->sign = si;
}


double Function::Local::Identifier::EvalFunction(const Local* local, bool grad_glob, double von_mises_stress) const
{
  // function value
  double fv = 0.0;
  Function* f = local->func_;

  switch(f->type_)
  {
  case STRESS:
  case STRESS_DENSITY:
    assert(von_mises_stress >= 0);
    fv = von_mises_stress;
    break;

  case SLOPE:
  case GLOBAL_SLOPE:
    fv = CalcSlope();
    break;

  case OSCILLATION:
  case GLOBAL_OSCILLATION:
    fv = CalcOscillation(local->GetBeta());
    break;

  case MOLE:
  case GLOBAL_MOLE:
    fv = CalcMole(local->GetEps());
    break;

  case JUMP:
  case GLOBAL_JUMP:
    fv = CalcJump();
    break;

  case BUMP:
    fv = CalcBump();
    break;

  case SUM_MODULI:
  case GLOBAL_SUM_MODULI:
    fv = CalcSumModuli();
    break;

  case PARAM_PS_POS_DEF:
    fv = CalcParamPSPosDef(-1, false);
    break;

  default:
    assert(false);
    break;
  }

  LOG_DBG2(func) << "L:I:EF: f=" << f->type.ToString(f->type_)
                 << " de=" << element->elem->elemNum << " sign=" << sign << " fv=" << fv;

  // handle globalization
  switch(f->type_)
  {
  case GLOBAL_SLOPE:
  case GLOBAL_OSCILLATION:
  case GLOBAL_MOLE:
  case GLOBAL_JUMP:
  case STRESS:
  case STRESS_DENSITY:
  {
    // we normalize all values by the number of "constraints". Note that it is
    // sufficient for the function value, the gradient is then also right
    double factor = local->DoNormalizeGlobal() ? 1.0/local->virtual_elem_map.GetSize() : 1.0;

    double v = std::max(0.0, fv - f->GetParameter());

    double p = local->GetPower();

    double res = grad_glob ? p * std::pow(v, p-1.0) : std::pow(v, p);

    res *= factor;

    LOG_DBG2(func) << "L:I:EF: global! bound=" << f->GetParameter() << " factor=" << factor
                   << " grad_glob=" << grad_glob << " power=" << std::pow(v, local->GetPower()) << " -> " << res;

    return res;
  }
  case GLOBAL_SUM_MODULI:
  {
    double factor = local->DoNormalizeGlobal() ? 1.0/local->virtual_elem_map.GetSize() : 1.0;
    return fv*factor;
  }
  default:
    return fv; // check is done before
  }
}

void Function::Local::Identifier::EvalGradient(const Local* local)
{
  // TODO the dynamic_cast might be to slow, check! and do faster by IsObjective()
  // we need this pointers, note that C++ makes NULL for an invalid dynamic cast
  Function* funct = local->func_;
  Function::Type ft = funct->type_;
  Condition* g = dynamic_cast<Condition*>(funct);
  Objective* f = dynamic_cast<Objective*>(funct);
  assert((f == NULL && g != NULL) || (f != NULL && g == NULL));

  LOG_DBG2(func) << "L:I:EvalGrad: f=" << funct->type.ToString(funct->type_) << " de="
                     << element->elem->elemNum << " sign=" << sign;

  // are we global? then we don't do anything if the globalization function gives zero
  // this applies the gradient of the globalization function (max(0, fv)^2)
  double grad_glob_fv = local->IsGlobalized() ? EvalFunction(local, true) : 0.0;

  if(local->IsGlobalized() && grad_glob_fv == 0.0)
  {
    LOG_DBG2(func) << "L:I:EvalGrad: f=" << funct->type.ToString(funct->type_) << " de="
                   << element->elem->elemNum << " sign=" << sign << " fv=0.0 -> return immediately";
    return;
  }
  assert(local->IsGlobalized() || g != NULL); // only constraints are local

  for(int n = -1, nn = neighbor.GetSize(); n < nn; n++)
  {
    double gv = -5.0;

    switch(ft)
    {
    case SLOPE:
    case GLOBAL_SLOPE:
      gv = CalcSlopeGradient(n);
      break;

    case OSCILLATION:
    case GLOBAL_OSCILLATION:
      gv = CalcOscillationGradient(n, local->beta_);
      break;

    case MOLE:
    case GLOBAL_MOLE:
      gv = CalcMoleGradient(n, local->eps_);
      break;

    case JUMP:
    case GLOBAL_JUMP:
      gv = CalcJumpGradient(n);
      break;

    case BUMP:
      gv = CalcBumpGradient(n);
      break;

    case STRESS:
    case STRESS_DENSITY:
      assert(false); // in SIMP::CalcVonMisesStressGradient() only!
      break;

    case SUM_MODULI:
      CalcSumModuliGradient(n, f, g, 1.0);
      continue;

    case GLOBAL_SUM_MODULI:
      CalcSumModuliGradient(n, f, g, local->DoNormalizeGlobal() ? 1.0/local->virtual_elem_map.GetSize() : 1.0);
      continue;

    case PARAM_PS_POS_DEF:
      gv = CalcParamPSPosDef(n, true);
      break;

    default:
      assert(false);
      break;
    }

    LOG_DBG2(func) << "L:I:EvalGrad: f=" << funct->type.ToString(funct->type_) << " de="
                   << element->elem->elemNum << " sign=" << sign << " n=" << n
                   << " curr=" << GetElement(n)->elem->elemNum << " gv=" << gv;

    // post process the globalized functions
    if(local->IsGlobalized())
    {
//      double factor = local->DoNormalizeGlobal() ? 1.0/local->virtual_elem_map.GetSize() : 1.0;

      gv  *= grad_glob_fv;
      LOG_DBG2(func) << "L:I:EvalGrad: f=" << funct->type.ToString(funct->type_) << " de="
                     << element->elem->elemNum << " sign=" << sign << " n=" << n
                     << " curr=" << GetElement(n)->elem->elemNum
                     << " bound! grad_glob_gv=" << grad_glob_fv << " new gv=" << gv;
    }

    DesignElement* de = GetElement(n);

    if(!local->IsGlobalized())
    {
      // reset the constraint data. Note, as we are local, there are no side effects by elements
      de->Reset(DesignElement::CONSTRAINT_GRADIENT, g);
    }

    de->AddGradient(f, g, gv);
    LOG_DBG2(func) << "L:I:EvalGrad: f=" << funct->type.ToString(funct->type_) << " de="
                   << element->elem->elemNum << " sign=" << sign << " n=" << n
                   << " curr=" << GetElement(n)->elem->elemNum << " gv=" << gv
                   << " stored_gv=" << de->GetPlainGradient(f, g);
  }
}


double Function::Local::Identifier::CalcSlope() const
{
  double mine  = element->GetDesign(DesignElement::SMART);
  assert(this->neighbor.GetSize() == 1);
  double other = neighbor[0]->GetDesign(DesignElement::SMART);

  double s = this->sign == -1 ? -1.0 : 1.0;

  LOG_DBG3(func) << "L:I:CS de=" << element->elem->elemNum << " other=" << neighbor[0]->elem->elemNum
                 << " sign=" << sign << " slope -> " << (s * (mine - other));
  return s * (mine - other);
}

double Function::Local::Identifier::CalcSlopeGradient(int neigh_idx) const
{
  assert(neigh_idx == -1 || neigh_idx == 0);
  // we have the cases sign=1, sign=-1, NO_SIGN. NO_SIGN is handled as sign=-1
  if(neigh_idx == -1)
    return sign == -1 ? -1.0 : 1.0;
  else
    return sign == -1 ? 1.0 : -1.0;
}



double Function::Local::Identifier::CalcOscillation(double beta) const
{
  assert(sign == 1 || sign == -1);

  double own  = element->GetDesign(DesignElement::SMART);
  // we divide the neighbors in lower and upper
  unsigned int half = neighbor.GetSize() / 2;

  tmp1.Resize(half);
  for(unsigned int i = 0; i < half; i++)
    tmp1[i] = neighbor[i]->GetDesign(DesignElement::SMART);
  double prev = sign == -1 ? SmoothMax(tmp1, beta) : SmoothMin(tmp1, beta);

  tmp2.Resize(half);
  for(unsigned int i = 0; i < half; i++)
    tmp2[i] = neighbor[i + half]->GetDesign(DesignElement::SMART);
  double next = sign == -1 ? SmoothMax(tmp2, beta) : SmoothMin(tmp2, beta);

  double min_max = 0.0; // min or max
  double res = 0.0;

  if(sign == 1)
  {
    // "Heaviside"(rho_i - max( rho_i-1, rho_i+1)
    // double smaller = std::max(0.0, own - std::max(left, right));
    min_max = beta < 0 ? std::max(prev, next) : SmoothMax(prev, next, beta);
    res = own - min_max;
  }
  else
  {
    // "Heaviside"(min( rho_i-1, rho_i+1) - rho_i)
    // double larger = std::max(0.0, std::min(left, right) - own);
    min_max = beta < 0 ? std::min(prev, next) : SmoothMin(prev, next, beta);
    res = min_max - own;
  }

  LOG_DBG3(func) << "L:I:CO de=" << element->ToString() << " neigh=" << DesignElement::ToString(neighbor)
                 << " vals=" << tmp1.ToString() << "; " << tmp2.ToString() << " sign=" << sign << " own=" << own
                 << " prev=" << prev << " next=" << next << " smooth=" << min_max << " hard="
                 << (sign == 1 ? (own - std::max(prev, next)) : (std::min(prev, next) - own)) << " -> " << res;
  return res;
}


double Function::Local::Identifier::CalcOscillationGradient(int neigh_idx, double beta)
{
  assert(beta >= 0);

  // the own value is not within the min/max stuff
  if(neigh_idx == -1)
    return sign == 1 ? 1 : -1;

  // we divide the neighbors in lower and upper
  unsigned int half = neighbor.GetSize() / 2;

  tmp1.Resize(half);
  for(unsigned int i = 0; i < half; i++)
    tmp1[i] = neighbor[i]->GetDesign(DesignElement::SMART);
  double prev = sign == -1 ? SmoothMax(tmp1, beta) : SmoothMin(tmp1, beta);

  tmp2.Resize(half);
  for(unsigned int i = 0; i < half; i++)
    tmp2[i] = neighbor[i + half]->GetDesign(DesignElement::SMART);
  double next = sign == -1 ? SmoothMax(tmp2, beta) : SmoothMin(tmp2, beta);

  // sign = -1: min( max(x_0 .. x_half-1), max(x_half .. x_max) ) - own
  // sign =  1: own - max( min(x_0 .. x_half-1), min(x_half .. x_max) )

  // example derivative for sign = 1: - max'(min(x_0 .. x_half-1), min(x_half .. x_max) * min'(max(x_0 .. x_half-1) * 1

  // are we within prev (-1) or next (1)
  int side = neigh_idx < (int) half ? -1 : 1;

  // outer is simple
  double outer = sign == 1 ? DerivSmoothMax(prev, next, beta, side) : DerivSmoothMin(prev, next, beta, side);
  // we do not handle neighbor.GetSize() == 2 special as the Smooth tools are fast for this special case
  // inner depends on neigh_idx within the first or the next
  double inner = 0.0;
  if(side == -1)
    inner = sign == -1 ? DerivSmoothMax(tmp1, beta, neigh_idx) : DerivSmoothMin(tmp1, beta, neigh_idx);
  else
    inner = sign == -1 ? DerivSmoothMax(tmp2, beta, neigh_idx - half) : DerivSmoothMin(tmp2, beta, neigh_idx - half);
  assert(neighbor.GetSize() > 2 || close(inner, 1.0));

  return (sign == 1 ? -1.0 : 1.0) * outer * inner;
}



double Function::Local::Identifier::CalcMole(double eps) const
{
  // the neighborhood is even and stars with the most previous element
  // and ends with the most next element.
  // the own element in the center is not in neighbor but is separate
  // thats not ideal for the mole constraint but it is the way it is.
  //
  // based on Poulsen; A new scheme for imposing a minimum length scale
  // in topology optimization; 2003
  //
  // the form is M(x) = ( sum_(i=1 to n-1) | x_i+1 - x_i | ) - | x_n - x_1 |
  // the abs are regularized by A(x) = sqrt(x^2 + eps^2) - eps

  // make sequence of data
  tmp1.Resize(0); // is static and keeps capacity
  for(unsigned int i = 0; i < neighbor.GetSize() / 2; i++)
  {
    tmp1.Push_back(GetElement(i)->GetDesign(DesignElement::SMART));
    LOG_DBG3(func) << "L:I:CalcMole de=" << element->ToString() << " other=" << GetElement(i)->ToString() << " v=" << tmp1.Last();
  }

  tmp1.Push_back(element->GetDesign(DesignElement::SMART));
  LOG_DBG3(func) << "L:I:CalcMole de=" << element->ToString() << " other=" << element->ToString() << " v=" << tmp1.Last();

  for(unsigned int i = neighbor.GetSize() / 2; i < neighbor.GetSize(); i++)
  {
    tmp1.Push_back(GetElement(i)->GetDesign(DesignElement::SMART));
    LOG_DBG3(func) << "L:I:CalcMole de=" << element->ToString() << " other=" << GetElement(i)->ToString() << " v=" << tmp1.Last();
  }

  assert(tmp1.GetSize() == neighbor.GetSize() + 1);

  double sum = 0.0;
  for(unsigned int i = 0; i < tmp1.GetSize() - 1; i++)
  {
    sum += SmoothAbs(tmp1[i+1] - tmp1[i], eps);
    LOG_DBG3(func) << "L:I:CalcMole de=" << element->ToString() << " i=" << i << "|" << tmp1[i+1] << " - " << tmp1[i] << "|="
                   << (std::abs(tmp1[i+1] - tmp1[i])) << " smoothed=" << SmoothAbs(tmp1[i+1] - tmp1[i], eps) << " sum=" << sum;
  }
  double result = sum - SmoothAbs(tmp1.Last() - tmp1[0], eps);
  LOG_DBG3(func) << "L:I:CalcMole de=" << element->ToString() << " bound=|" << tmp1.Last() << " - " << tmp1[0] << "| smoothed="
                 << SmoothAbs(tmp1.Last() - tmp1[0], eps) << " -> " << result;
  return result;
}

double Function::Local::Identifier::CalcMoleGradient(int neigh_idx, double eps)
{
  // see comments in the forward function implementation
  // three cases for the sensitivity analysis: first, intermediate, last

  // sort all values in a sequence -> TODO! optimize and do not copy and pase
  tmp1.Resize(0); // is static and keeps capacity
  for(unsigned int i = 0; i < neighbor.GetSize() / 2; i++)
    tmp1.Push_back(GetElement(i)->GetDesign(DesignElement::SMART));

  tmp1.Push_back(element->GetDesign(DesignElement::SMART));

  for(unsigned int i = neighbor.GetSize() / 2; i < neighbor.GetSize(); i++)
      tmp1.Push_back(GetElement(i)->GetDesign(DesignElement::SMART));

  assert(tmp1.GetSize() == neighbor.GetSize() + 1);

  // correct the index
  unsigned int half = neighbor.GetSize() / 2;
  int idx = neigh_idx == -1 ? half : neigh_idx < (int) half ?  neigh_idx : neigh_idx + 1;
  assert(idx >= 0 && idx < (int) tmp1.GetSize());

  double res = 0.0;

  if(idx == 0)
  {
    res = DerivSmoothAbs(tmp1.Last() - tmp1[0], eps) - DerivSmoothAbs(tmp1[1] - tmp1[0], eps);
    LOG_DBG3(func) << "L:I:CalcMoleGrad de=" << element->ToString() << " neigh_idx=" << neigh_idx << " idx=" << idx << " tmp=" << tmp1.ToString()
                   << " DA(" << tmp1.Last() << "-" << tmp1[0] << ") - DA(" << tmp1[1] << "-" << tmp1[0] << ") -> " << res;

  }
  else if(idx == (int) tmp1.GetSize() - 1)
  {
    res = DerivSmoothAbs(tmp1.Last() - tmp1[idx-1], eps) - DerivSmoothAbs(tmp1.Last() - tmp1[0], eps);
    LOG_DBG3(func) << "L:I:CalcMoleGrad de=" << element->ToString() << " neigh_idx=" << neigh_idx << " idx=" << idx << " tmp=" << tmp1.ToString()
                   << " DA(" << tmp1.Last() << "-" << tmp1[idx-1] << ") - DA(" << tmp1.Last() << "-" << tmp1[0] << ") -> " << res;

  }
  else
  {
    res = DerivSmoothAbs(tmp1[idx] - tmp1[idx-1], eps) - DerivSmoothAbs(tmp1[idx+1] - tmp1[idx], eps);
    LOG_DBG3(func) << "L:I:CalcMoleGrad de=" << element->ToString() << " neigh_idx=" << neigh_idx << " idx=" << idx << " tmp=" << tmp1.ToString()
                   << " DA(" << tmp1[idx] << "-" << tmp1[idx-1] << ") - DA(" << tmp1[idx+1] << "-" << tmp1[idx] << ") -> " << res;
  }

  return res;
}

double Function::Local::Identifier::CalcJump() const
{
  assert(sign == NO_SIGN);
  assert(neighbor.GetSize() == 2);

  // sin(pi*(x_i-1 - x_+1))^2
  // no own value!
  double prev = neighbor[0]->GetDesign(DesignElement::SMART);
  double next = neighbor[1]->GetDesign(DesignElement::SMART);

  double sin = std::sin(PI*(prev-next));

  LOG_DBG3(func) << "L:I:CJ de=" << element->ToString() << " prev=" << neighbor[0]->ToString() << "/" << prev
                 << " next=" << neighbor[1]->ToString() << "/" << next << " slope=" << (prev-next)
                 << " -> sin*sin";

  return sin*sin;
}

double Function::Local::Identifier::CalcJumpGradient(int neigh_idx) const
{
  // g(x)=sin(pi*(x_i-1 - x_+1))^2
  // d g(x)/d x_i-1 = 2 * sin(pi*(x_i-1 - x_+1)) * cos (pi*(x_i-1 - x_+1)) * PI

  // no own value!
  if(neigh_idx == -1)
    return 0.0;

  double prev = neighbor[0]->GetDesign(DesignElement::SMART);
  double next = neighbor[1]->GetDesign(DesignElement::SMART);

  double slope = prev-next;

  assert(neigh_idx == 0 || neigh_idx == 1);
  double factor = neigh_idx == 0 ? 1.0 : -1.0;

  return 2.0 * std::sin(PI*slope) * std::cos(PI*slope) * PI * factor;
}

double Function::Local::Identifier::CalcBump() const
{
  assert(sign == NO_SIGN);
  assert(neighbor.GetSize() == 2);

  // (x_i-1 - x_i)(x_i - x_i+1)
  // no own value!
  double prev = neighbor[0]->GetDesign(DesignElement::SMART);
  double mine = element->GetDesign(DesignElement::SMART);
  double next = neighbor[1]->GetDesign(DesignElement::SMART);

  double val = (prev - mine) * (mine - next);

  LOG_DBG3(func) << "L:I:CB de=" << element->ToString()
                 << " prev=" << neighbor[0]->ToString() << "/" << prev
                 << " next=" << neighbor[1]->ToString() << "/" << next << " mine=" << mine
                 << " -> " << (prev - mine) << " * " << (mine - next) << " = " << val;

  return val;
}


double Function::Local::Identifier::CalcBumpGradient(int neigh_idx) const
{
  double prev = neighbor[0]->GetDesign(DesignElement::SMART);
  double mine = element->GetDesign(DesignElement::SMART);
  double next = neighbor[1]->GetDesign(DesignElement::SMART);


  double res = 0.0;

  switch(neigh_idx)
  {
  case 0:
    res = mine - next;
    break;

  case -1:
    res = prev - 2.0 * mine + next;
    break;

  case 1:
    res = mine - prev;
    break;

  default:
    assert(false);
    break;
  }

  return res;
}


double Function::Local::Identifier::CalcSumModuli() const
{
  double E1(0.0), E3(0.0), G(0.0);
  for(int i=-1; i < (int) neighbor.GetSize(); ++i)
  {
    switch(GetElement(i)->GetType())
    {
    case DesignElement::EMODULISO:
      E1 = GetElement(i)->GetDesign(DesignElement::PLAIN);
      break;
    case DesignElement::EMODUL:
      E3 = GetElement(i)->GetDesign(DesignElement::PLAIN);
      break;
    case DesignElement::GMODUL:
      G = GetElement(i)->GetDesign(DesignElement::PLAIN);
      break;
    case DesignElement::POISSON:
    case DesignElement::POISSONISO:
    case DesignElement::DENSITY:
      break;

    default:
      assert(false);
      break;
    }
  }
  return E1+E3+G;
}

void Function::Local::Identifier::CalcSumModuliGradient(int neigh_idx, const Objective* f, const Condition* g, double value)
{
  DesignElement::Type type = GetElement(neigh_idx)->GetType();
  if (type == DesignElement::EMODULISO || type == DesignElement::EMODUL || type == DesignElement::GMODUL)
    GetElement(neigh_idx)->AddGradient(f, g, value);
}

double Function::Local::Identifier::CalcParamPSPosDef(int neigh_idx, bool derivative) const
{
  double E1(0.0), E3(0.0), nu31(0.0);
  for(int i=-1; i < (int) neighbor.GetSize(); ++i)
  {
    switch(GetElement(i)->GetType())
    {
    case DesignElement::EMODULISO:
      E1 = GetElement(i)->GetDesign(DesignElement::PLAIN);
      break;
    case DesignElement::EMODUL:
      E3 = GetElement(i)->GetDesign(DesignElement::PLAIN);
      break;
    case DesignElement::POISSON:
      nu31 = GetElement(i)->GetDesign(DesignElement::PLAIN);
      break;
    case DesignElement::GMODUL:
    case DesignElement::POISSONISO:
    case DesignElement::DENSITY:
      break;

    default:
      assert(false);
      break;
    }
  }
  if(derivative)
    switch(GetElement(neigh_idx)->GetType())
    {
    case DesignElement::EMODULISO:
      return -nu31*nu31;
    case DesignElement::EMODUL:
      return 1.0;
    case DesignElement::POISSON:
      return -2.0*nu31*E1;
    default:
      return 0.0;
    }
  else
  {
    LOG_DBG3(func) << "Local::Local e_num=" << element->elem->elemNum << ", E3-E1*nu31^2=" << E3-E1*nu31*nu31;
    return E3-E1*nu31*nu31;
  }
}


