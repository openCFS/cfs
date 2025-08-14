#include <algorithm>
#include <cmath>
#include <ostream>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/ParamHandling/XmlReader.hh"
#include "Domain/Domain.hh"
#include "Domain/ElemMapping/Elem.hh"
#include "Domain/Mesh/Grid.hh"
#include "Driver/BaseDriver.hh"
#include "FeBasis/BaseFE.hh"
#include "General/Enum.hh"
#include "General/defs.hh"
#include "General/Environment.hh"
#include "General/Exception.hh"
#include "MatVec/Matrix.hh"
#include "MatVec/Vector.hh"
#include "Materials/MechanicMaterial.hh"
#include "Optimization/Condition.hh"
#include "Optimization/Design/AuxDesign.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Design/DesignStructure.hh"
#include "Optimization/Design/ShapeDesign.hh"
#include "Optimization/Design/ShapeMapDesign.hh"
#include "Optimization/Design/SplineBoxDesign.hh"
#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/Excitation.hh"
#include "Optimization/Function.hh"
#include "Optimization/Objective.hh"
#include "Optimization/TransferFunction.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/BasePDE.hh"
#include "Utils/CubicInterpolate.hh"
#include "Utils/BiCubicInterpolate.hh"
#include "Utils/TriCubicInterpolate.hh"
#include "Utils/StdVector.hh"
#include "Utils/tools.hh"

DEFINE_LOG(func, "opt_func")

// instantiation of the static elements is in Optimization::SetEnums()
Enum<Function::Type> Function::type;
Enum<Function::SlackFnct> Function::slackFnct;
Enum<Function::MultiObjType> Function::multiObjType;
Enum<Function::Access> Function::access;
Enum<Function::StressType> Function::stressType;
Enum<Function::Local::Locality> Function::Local::locality;
Enum<Function::Local::Phase> Function::Local::phase;

// speed up by sharing
Vector<double> Function::Local::Identifier::tmp1;
Vector<double> Function::Local::Identifier::tmp2;

// sync the values with Local::Phase
const int Function::Local::Identifier::NO_SIGN = -1000;
const int Function::Local::Identifier::VOID_SIGN = -1;
const int Function::Local::Identifier::MATERIAL_SIGN = 1;

using boost::lexical_cast;
using std::string;

Function::Function()
{
  Init();
}

Function::Function(PtrParamNode pn)
{
  Init();

  this->type_ = type.Parse(pn->Get("type")->As<string>());

  this->preInfo_ = PtrParamNode(new ParamNode(ParamNode::INSERT, ParamNode::ELEMENT));
  this->pn = pn;

  if (pn->Has("design")) // will sometime be in Function, now the default is set to DEFAULT
    this->design_ = BaseDesignElement::type.Parse(
        pn->Get("design")->As<string>());

  this->parameter_ = pn->Has("parameter") ? pn->Get("parameter")->As<double>() : 0.0;

  this->omega_omega_ = pn->Has("factor") ? pn->Get("factor/omega_omega")->As<bool>() : false;
  // FIXME
  //if (!complex_ && omega_omega_)
  //  throw Exception("It makes no sense to set costFunction/factor/omega_omega in static optimization");

  this->eigenvalue_id_ = pn->Has("ev") ? pn->Get("ev")->As<unsigned int>() : 0;

  this->type_ = type.Parse(pn->Get("type")->As<string>());

  // we know the filter type only very late in some post init spaghetti.
  // so guess now and confirm in PostProc() when we have the design space
  this->filter_ = DesignStructure::GuessFilterType();
  this->access_ = access.Parse(pn->Get("access")->As<string>()); // default in xml
  if(access_ == DEFAULT) {
    LOG_DBG(func) << "F:F t=" << type.ToString(type_) << " access -> " << access.ToString(DefaultAccess(type_));
    access_ = DefaultAccess(type_);
  }
  assert(access_ != DEFAULT);

  slackFnct_ = slackFnct.Parse(pn->Get("function")->As<string>());

  // default is set in Function, may this moves later to Function, too
  if(pn->Has("region") && pn->Get("region")->As<string>() != "all")
    region = domain->GetGrid()->GetRegion().Parse(pn->Get("region")->As<string>());

  if(type_ == SLACK_FNCT && slackFnct_ == NO_FUNCTION)
    EXCEPTION("a function 'slackFunction' requires the attribute 'function' to be set");

  if(type_ != SLACK_FNCT && slackFnct_ != NO_FUNCTION)
    preInfo_->SetWarning("providing 'function' " + slackFnct.ToString(slackFnct_) + " names only sense for type " + type.ToString(SLACK_FNCT));

  if(type_ == BANDGAP)
  {
    if(!pn->Has("bandgap"))
      throw Exception("function 'bandgap' required child element 'bandgap'");
    bandgap.lower_ev = pn->Get("bandgap/lower_ev")->As<int>();
    bandgap.upper_ev = pn->Get("bandgap/upper_ev")->As<int>();
    if(bandgap.lower_ev >= bandgap.upper_ev)
      throw Exception("within 'bandgap' 'lower_ev' needs to be smaller than 'upper_ev'");
    if(bandgap.upper_ev - bandgap.lower_ev > 1)
      preInfo_->SetWarning("'bandgap' defines a gap non-adjacent modes");
  }

  int sequence = pn->Get("sequence")->As<int>();
  if(sequence > (int) Optimization::manager.context.GetSize()) // note 1-based!
    EXCEPTION("too high sequence number " << sequence << " for function " << type.ToString(type_));
  this->ctxt = &(Optimization::manager.context[sequence - 1]);

  notation_ = VOIGT;
  if(pn->Has("notation"))
    notation_ = tensorNotation.Parse(pn->Get("notation")->As<string>());

  bool tensor_ok = ReadTensor(ctxt, pn, this->tensor_); // is save and sets default

  if ((type_ == HOM_TRACKING || type_ == HOM_FROBENIUS_PRODUCT) && !tensor_ok)
    EXCEPTION("A 'tensor' element is mandatory  for 'homTracking'");

  if(type_ == HOM_TENSOR || type_ == HOM_TRACKING) {
    // we must not give a value when there is a tensor
    if (type_ == HOM_TENSOR && pn->Has("tensor") && pn->Has("value"))
      throw Exception("a value must not be given when a tensor is used in a homogenization constraint");
  }

  if (type_ == HOM_TRACKING && (!pn->Has("tensor") && !pn->Has("isotropic")))
    throw Exception("a 'tensor' is mandatory for homogenization tracking");

  // check parameter
  switch (type_) {
  case PENALIZED_VOLUME:
  case GAP:
  case GLOBAL_SLOPE:
  case GLOBAL_OSCILLATION:
  case GLOBAL_MOLE:
  case GLOBAL_STRESS:
  case GLOBAL_BUCKLING_LOAD_FACTOR:
  case PERIMETER:
    if (!pn->Has("parameter"))
      throw Exception("function '" + type.ToString(type_) + "' requires the 'parameter' attribute");
    break;

  case EIGENFREQUENCY:
    if(!pn->Has("ev"))
      throw Exception("function '" + type.ToString(type_) + "' requires the 'ev' with value >= 1");
    break;

  case TENSOR_TRACE:
  case GLOBAL_TENSOR_TRACE:
    if(design_ != DesignElement::DEFAULT && design_ != DesignElement::MECH_TRACE && design_ != DesignElement::DIELEC_TRACE && design_ != DesignElement::ALL_DESIGNS)
      throw Exception("function '" + type.ToString(type_) + "' has invalid design type " + DesignElement::type.ToString(design_));
    break;

  case POS_DEF_DET_MINOR_1:
  case POS_DEF_DET_MINOR_2:
  case POS_DEF_DET_MINOR_3:
  case BENSON_VANDERBEI_1:
  case BENSON_VANDERBEI_2:
  case BENSON_VANDERBEI_3:
    if (design_ != DesignElement::MECH_ALL && design_ != DesignElement::DIELEC_ALL)
      throw Exception("mandatory 'design' for '" + type.ToString(type_) + "' is 'elast_all' and 'dielec_all'");

    break;

  case OVERHANG_VERT:
  case OVERHANG_HOR:
    if(!BaseDesignElement::IsShapeMapType(design_))
      throw Exception("'overhang' function requires design to be set to shape variables ('shape_map')");
    break;

  case SQR_MAG_FLUX_DENS_X:
  case SQR_MAG_FLUX_DENS_Y:
    if(domain->GetGrid()->IsAxi())
      throw Exception("not for axis symmetric setting: " + type.ToString(type_));
    break;

  case SQR_MAG_FLUX_DENS_RZ:
  case LOSS_MAG_FLUX_RZ:
  if(!domain->GetGrid()->IsAxi())
      throw Exception("only for axis symmetric setting: " + type.ToString(type_));
  break;

  case CONES:
    if(!BaseDesignElement::IsSplineBoxType(design_))
      throw Exception("'cones' function requires design to be set to spline box variables");
    break;

  default:
    break;
  }

  // set linear, to be overwritten by xml below
  switch (type_)
  {
  case VOLUME: // the volume is not linear on heaviside densities and shape mapping
  case SLOPE:
  case SLACK:
  case MULTIMATERIAL_SUM:
  case SHAPE_INF:
  case PERIODIC:
  case CURVATURE:
  case OVERHANG_VERT:
  case CONES: // might be also nonlinear in the future
    linear_ = true;
    break;
//  case TENSOR_TRACE:
//  case GLOBAL_TENSOR_TRACE:
//    if (design_ != DesignElement::ALL_DESIGNS) // This actually depends on the material parametrization, not design_ !
//      linear_ = false;
//    else
//      linear_ = true; // rarely true, has to be set in xml now
//    break;
  default:
    linear_ = false;
    break;
  }

  // snopt only makes a difference between linear and nonlinear constraints!
  if(pn->Has("linear"))
    linear_ = pn->Get("linear")->As<bool>();

}

Function::~Function()
{
  if (local != NULL) {
    delete local;
    local = NULL;
  }

  // this might lead to problems when they are active in Assemble and ~Assemble deletes them
  output_forms.Clear();
}

void Function::Init()
{
  this->design_ = DesignElement::DEFAULT; // overwritten eventually from xml
  this->region = ALL_REGIONS;  // overwritten eventually in Condition

  this->local = NULL;

  // function value to be evaluated
  this->value_ = -1.0;

  this->excite_ = UNSET_EX;
  this->sample_excitation_ = NULL;
  this->excite_sensitive_ = false;

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

bool Function::ReadTensor(Context* f_ctxt, PtrParamNode pn, Matrix<double>& matrix)
{
  matrix.Resize(1, 1); // minimal size, as 0,0 is not defined.

  // sanity checks
  if (!pn->Has("tensor") && !pn->Has("isotropic"))
    return false;

  if (pn->Has("tensor") && pn->Has("isotropic"))
    EXCEPTION("please specify either <tensor> or <isotropic>, not both");

  bool tensor_read(false);

  // check for tensor element
  PtrParamNode tens = pn->Get("tensor", ParamNode::PASS);
  if (tens != NULL) {
    int dim = tens->Get("dim1")->As<int>();
    if (dim != 3 && dim != 6)
      EXCEPTION("The Voigt 'tensor' for homogenizations needs to be 3x3 or 6x6");
    if (tens->Has("dim2") && dim != tens->Get("dim2")->As<int>())
      EXCEPTION("The 'tensor' for homogenization needs to be symmetric");
    if ((domain->GetDim() == 2 && dim != 3) || (domain->GetDim() == 3 && dim != 6))
      EXCEPTION("The 'tensor' for homogenization needs to be 3x3 for 2D and 6x6 for 3D");

    matrix.Resize(dim, dim);

    ParamTools::AsTensor<double>(tens->Get("real"), dim, dim, matrix);

    // check for a scaling factor
    const double factor(tens->Get("factor")->As<double>());
    if (factor != 1.0)
      matrix *= factor;

    tensor_read = true;
  }

  tens = pn->Get("isotropic", ParamNode::PASS);
  if (tens != NULL) {
    double emod = tens->Get("real")->Get("elasticityModulus")->As<double>();
    double poisson = tens->Get("real")->Get("poissonNumber")->As<double>();

    Matrix<double> tmp(6, 6); // always 3D first
    MechanicMaterial::CalcIsotropicStiffnessTensorFromEAndPoisson(tmp, emod, poisson);
    assert(f_ctxt->stt != NO_TENSOR);
    MechanicMaterial::ComputeSubTensor(matrix, f_ctxt->stt, tmp);

    tensor_read = true;
  }

  if (tensor_read) {
    // output the target tensor to xml-file
    // to get a reference and be able to do quick checks
    // FIXME
    PtrParamNode in_mat = domain->GetInfoRoot()->Get("target_tensor");
    in_mat->SetType(ParamNode::ELEMENT);
    in_mat->Get("tensor")->SetValue(matrix);
  }

  return tensor_read;
}


void Function::ParseCoord(PtrParamNode pn, boost::tuple<int, int, double>& coord) {
  string val = pn->Get("coord")->As<string>();
  boost::get<0>(coord) = lexical_cast<unsigned int>(val.at(0));
  boost::get<1>(coord) = lexical_cast<unsigned int>(val.at(1));
  boost::get<2>(coord) = 1.0; // default
}

void Function::ToInfo(PtrParamNode info) {
  info_ = info;

  // there might be set something, i.g. in PostProc
  info_->SetValue(preInfo_, false, false); // don't do tricks with name and don't repeat warning prints

  info->Get("type")->SetValue(type.ToString(type_));

  if(type_ == SLACK_FNCT)
    info->Get("function")->SetValue(slackFnct.ToString(slackFnct_));

  info->Get("name")->SetValue(ToString());

  if(Optimization::context->IsComplex() && omega_omega_) // reduce output
    info->Get("omega_omega")->SetValue(omega_omega_);
  // we check for valid occurrence of parameter in the constructor
  if(pn->Has("parameter") || IsLocal(type_))
    info->Get("parameter")->SetValue(parameter_);

  if(type_ == GREYNESS)
    info->Get("gray_scale")->SetValue(CalcGraynessScaling());

  // We might have non-standard stresses
  if(type_ == GLOBAL_STRESS  || type_ == LOCAL_STRESS)
    info->Get("stress")->SetValue(stressType.ToString(stressType_));

  if(type_ == EIGENFREQUENCY || type_ == BUCKLING_LOAD_FACTOR || type_ == LOCAL_BUCKLING_LOAD_FACTOR)
    info->Get("ev")->SetValue(eigenvalue_id_);

  if(IsObjective() || !(dynamic_cast<Condition*>(this)->IsObservation()))
    info->Get("linear")->SetValue(linear_);

  // only a small fraction of functions as volume and grayness really handles this information.
  // when set to default in xml (default setting) we take the setting from Function::DefaultAccess()
  // when not default, the function name <access>_<type>
  info->Get("access")->SetValue(access.ToString(GetAccess()));

  if(region != ALL_REGIONS)
    info->Get("region")->SetValue(domain->GetGrid()->GetRegion().ToString(region));

  if(local != NULL)
    local->ToInfo(info_);
}

string Function::ToString() const
{
  if(type_ == SLACK_FNCT)
  {
    switch(slackFnct_)
    {
    case NORM_BANDGAP:
      return "sf_2s_by_a";
    case REL_BANDGAP:
      return "sf_2s_by_a-s";
    case ALPHA_SLACK_QUOTIENT:
      return "sf_a_by_s";
    case ALPHA_MINUS_SLACK:
      return "sf_a-s";
    case NO_FUNCTION:
      assert(false);
      return "no_funct";
    }
  }

  string tn = (type_ == PYTHON_FUNCTION || type_ == LOCAL_PYTHON_FUNCTION) ? py_name_ : type.ToString(type_);

  // optional for oscillation
  if(local != NULL && local->GetPhase() != Local::BOTH)
    return Local::phase.ToString(local->GetPhase()) + "_" + tn;

  std::ostringstream os;

  if(IsDefaultAccess())
    os << tn;
  else
    os << access.ToString(access_) + "_" + tn;

  if(type_ == BUCKLING_LOAD_FACTOR || type_ == LOCAL_BUCKLING_LOAD_FACTOR)
    os << "_" << eigenvalue_id_;

  return os.str();
}

int Function::CountOscillations() const
{
   int cnt = 0;
   for(unsigned int i = 2; i < history.GetSize(); i++)
   {
     double pp = history[i-2];
     double p = history[i-1];
     double c = history[i];
     if((p-pp)*(c-p) < 0)
       cnt++;
   }
   return cnt;
}



void Function::DescribeProperties(StdVector<std::pair<string,string> >& map) const
{
  map.Push_back(std::make_pair("name", ToString()));
  map.Push_back(std::make_pair("type", type.ToString(type_)));
  map.Push_back(std::make_pair("access", access.ToString(access_)));
  map.Push_back(std::make_pair("value", std::to_string(value_)));
  if(region != ALL_REGIONS)
    map.Push_back(std::make_pair("region",domain->GetGrid()->GetRegion().ToString(region)));
}


Function::Access Function::DefaultAccess(Function::Type type) const
{
  // filter_ needs to be set!!!

  // we ignore that filter might be NO_FILTERING. This is rarely, the neighborhood is zero and overhead small
  switch(type)
  {
  // PLAIN for density and sensitivity
  case SLACK:
  case SLACK_FNCT:
  case GLOBAL_SLOPE:
  case GLOBAL_MOLE:
  case GLOBAL_OSCILLATION:
  case GLOBAL_JUMP:
  case GLOBAL_CURVATURE:
  case GLOBAL_DESIGN:
  case DESIGN:
  case SLOPE:
  case MOLE:
  case OSCILLATION:
  case JUMP:
  case BUMP:
  case CURVATURE:
  case PERIODIC:
  case OVERHANG_VERT:
  case OVERHANG_HOR:
  case DISTANCE:
  case BENDING:
  case CONES:
  case SUM_MODULI:
  case GLOBAL_SUM_MODULI:
  case ORTHOTROPIC_TENSOR_TRACE:
  case GLOBAL_ORTHOTROPIC_TENSOR_TRACE:
  case TENSOR_TRACE:
  case TENSOR_NORM:
  case GLOBAL_TENSOR_TRACE:
  case PARAM_PS_POS_DEF:
  case POS_DEF_DET_MINOR_1:
  case POS_DEF_DET_MINOR_2:
  case POS_DEF_DET_MINOR_3:
  case BENSON_VANDERBEI_1:
  case BENSON_VANDERBEI_2:
  case BENSON_VANDERBEI_3:
  case MULTIMATERIAL_SUM:
  case SHAPE_INF:
  case EXPRESSION:
  case PYTHON_FUNCTION:
  case LOCAL_PYTHON_FUNCTION:
  case ARC_OVERLAP:
    return PLAIN;

  // filtered stuff different for sensitivity filtering
  // we do sensitivity almost never, so extend if you really need it for other stuff
  // one can check DesignSpace::GetFilterType() for background
  case VOLUME:
  case GREYNESS:
    if(filter_ == Filter::SENSITIVITY)
      return PLAIN;
    else
      return FILTERED;

  // filtered stuff, meant for density filtering
  case PENALIZED_VOLUME:
  case GAP:
  case TYCHONOFF: // not sure about it! Fabian
  case FILTERING_GAP:
  case TWO_SCALE_VOL:
  case GLOBAL_TWO_SCALE_VOL:
  case PERIMETER:
  case REALVOLUME:
  case LOCAL_STRESS: // not plain filtered but not the physical transfer function but own stress transfer function
    return FILTERED;

  // pyhsical is penalized and filtered for density filtering.
  // In the sensitivity filter case, the design is never filtered, only the gradient
  case BANDGAP:
  case OUTPUT:
  case SQUARED_OUTPUT:
  case DYNAMIC_OUTPUT:
  case REFLECTED_WAVE:
  case ABS_OUTPUT:
  case CONJUGATE_COMPLIANCE:
  case GLOBAL_DYNAMIC_COMPLIANCE:
  case ELEC_ENERGY:
  case ENERGY_FLUX:
  case COMPLIANCE:
  case TRACKING:
  case HOM_TENSOR:
  case HOM_TRACKING:
  case HOM_FROBENIUS_PRODUCT:
  case POISSONS_RATIO:
  case YOUNGS_MODULUS:
  case YOUNGS_MODULUS_E1:
  case YOUNGS_MODULUS_E2:
  case TEMPERATURE:
  case HEAT_ENERGY:
  case SQR_MAG_FLUX_DENS_X:
  case SQR_MAG_FLUX_DENS_Y:
  case SQR_MAG_FLUX_DENS_RZ:
  case LOSS_MAG_FLUX_RZ:
  case MAG_COUPLING:
  case TEMP_TRACKING_AT_INTERFACE:
  case GLOBAL_STRESS:
  case EIGENFREQUENCY:
  case BUCKLING_LOAD_FACTOR:
  case LOCAL_BUCKLING_LOAD_FACTOR:
  case GLOBAL_BUCKLING_LOAD_FACTOR:
  case PRESSURE_DROP:
  case ISOTROPY:
  case ISO_ORTHOTROPY:
  case ORTHOTROPY:
  case DESIGN_TRACKING: // according to comment against physical design
    return PHYSICAL;

  case MULTI_OBJECTIVE:
  case NO_TYPE:
    assert(false);
    break;
  }
  return NO_ACCESS;
}

Function* Function::GetFunction(Objective* f, Condition* g) {
  assert(!(f != NULL && g != NULL) || (f == NULL && g == NULL));
  return f != NULL ? dynamic_cast<Function*>(f) : dynamic_cast<Function*>(g);
}

void Function::SetExcitation(MultipleExcitation* me, int excite_index)
{
  assert(me != NULL && me->excitations.GetSize() > 0);

  // some functions need to be evaluated only once (first) for multiple excitations
  // however for meta excitations (rotations) they need to be be evaluates at the last base
  //
  // multiple excitations are:
  // * static load cases
  // * different frequencies
  // * several homogenization test strains
  // * time steps
  // * bloch mode analysis wave vectors

  switch(type_)
  {
  // this stuff is really to be evaluated only once, even for meta excitations or multi sequence, but we stick
  // to the (default) sequence value
  case VOLUME:
  case PENALIZED_VOLUME:
  case GAP:
  case REALVOLUME:
  case TYCHONOFF:
  case GREYNESS:
  case SLOPE:
  case GLOBAL_SLOPE:
  case PERIMETER:
  case MOLE:
  case GLOBAL_MOLE:
  case OSCILLATION:
  case GLOBAL_OSCILLATION:
  case JUMP:
  case GLOBAL_JUMP:
  case BUMP:
  case CURVATURE:
  case GLOBAL_CURVATURE:
  case OVERHANG_VERT:
  case OVERHANG_HOR:
  case DISTANCE:
  case BENDING:
  case CONES:
  case DESIGN:
  case GLOBAL_DESIGN:
  case PERIODIC:
  case DESIGN_TRACKING:
  case SUM_MODULI:
  case GLOBAL_SUM_MODULI:
  case TWO_SCALE_VOL:
  case GLOBAL_TWO_SCALE_VOL:
  case PARAM_PS_POS_DEF:
  case POS_DEF_DET_MINOR_1:
  case POS_DEF_DET_MINOR_2:
  case POS_DEF_DET_MINOR_3:
  case BENSON_VANDERBEI_1:
  case BENSON_VANDERBEI_2:
  case BENSON_VANDERBEI_3:
  case ORTHOTROPIC_TENSOR_TRACE:
  case GLOBAL_ORTHOTROPIC_TENSOR_TRACE:
  case TENSOR_TRACE:
  case TENSOR_NORM:
  case GLOBAL_TENSOR_TRACE:
  case SHAPE_INF:
  case MULTIMATERIAL_SUM:
  case SLACK:
  case BANDGAP: // similar to bloch=extremal
  case SLACK_FNCT:
  case EXPRESSION:
  case FILTERING_GAP:
  case PYTHON_FUNCTION: // check in case!
  case LOCAL_PYTHON_FUNCTION:
  case ARC_OVERLAP:
    assert(excite_index < 0);
    excite_ = ctxt->excitations.Last()->index;
    break;

  // this stuff is to be evaluated at the last base for meta excitations
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
    assert(excite_index < 0);
    if(!me->DoMetaExcitation(ctxt))
      excite_ = ctxt->excitations.Last()->index; // with respect to our context
    else
    {
      if(!pn->Has("excitation"))
        throw Exception("doing homogenization with meta excitations the excitation parameters is mandatory for " + ToString());
      // assert(!ctxt->DoMultiSequence());
      excite_ = ctxt->GetExcitation(me->GetNumberHomogenization(ctxt->ToApp())-1, pn->Get("excitation")->As<string>())->index; // -1 to access the last
    }
    break;

  // this stuff is to be evaluated always
  case COMPLIANCE:
  case OUTPUT:
  case SQUARED_OUTPUT:
  case DYNAMIC_OUTPUT:
  case REFLECTED_WAVE:
  case ENERGY_FLUX:
  case TRACKING:
  case ABS_OUTPUT:
  case CONJUGATE_COMPLIANCE:
  case GLOBAL_DYNAMIC_COMPLIANCE:
  case ELEC_ENERGY:
  case TEMPERATURE:
  case PRESSURE_DROP:
  case HEAT_ENERGY:
  case SQR_MAG_FLUX_DENS_Y:
  case SQR_MAG_FLUX_DENS_X:
  case SQR_MAG_FLUX_DENS_RZ:
  case LOSS_MAG_FLUX_RZ:
  case TEMP_TRACKING_AT_INTERFACE:
    assert(excite_index < 0);
    if(!pn->Has("excitation") || pn->Get("excitation")->As<string>() == "all") {
      excite_ = ALL_EX; // all excitations within this sequence/ context
      // why is there no excite_sensitive_ = true; ??
    } else {
      const std::string& ex_label = pn->Get("excitation")->As<string>();
      const std::string& label = ctxt->DoMultiSequence() ? "s_" + lexical_cast<string>(ctxt->sequence) + "-" + ex_label : ex_label;
      excite_ = me->GetExcitation(label)->index;
      excite_sensitive_ = true;
    }
    break;

  case GLOBAL_STRESS:
  case LOCAL_STRESS:
  case EIGENFREQUENCY: // at least in the bloch mode case! Otherwise there is no multiple excitation for standard ev
  case BUCKLING_LOAD_FACTOR:
  case LOCAL_BUCKLING_LOAD_FACTOR:
  case GLOBAL_BUCKLING_LOAD_FACTOR: // blubber
    // there might be the optional excitation index set
    if (!pn->Has("excitation") || pn->Get("excitation")->As<string>() == "all") {
      excite_ = excite_index == UNSET_EX ? ALL_EX : excite_index;
    } else {
      assert(excite_index == UNSET_EX); // assert there is no conflict
      const std::string& ex_label = pn->Get("excitation")->As<string>();
      const std::string& label = ctxt->DoMultiSequence() ? "s_" + lexical_cast<string>(ctxt->sequence) + "-" + ex_label : ex_label;
      excite_ = me->GetExcitation(label)->index;
    }
    excite_sensitive_ = true;
    break;

  case MAG_COUPLING:
    // enforces the excitation to be manually set to "0_1" for the first two excitations
    assert(excite_index < 0);
    if(!pn->Has("excitation") || pn->Get("excitation")->As<string>() != "0_1")
       throw Exception("function " + type.ToString(MAG_COUPLING) + " requires excitation='0_1'");
    excite_ = COMBINED_0_1_EX;
    excite_sensitive_ = true;
    break;

  case MULTI_OBJECTIVE: // only to make the switch complete
  case NO_TYPE:
    assert(false);
    break;
  }

  sample_excitation_ = excite_ >= 0 ? &me->excitations[excite_] : &me->excitations[0];
  LOG_DBG(func) << "SE f=" << ToString() << " exite_=" << excite_ << " ex=" << sample_excitation_->GetFullLabel();
}

/** Shall/must we evaluate this function at this excitation?
 * Stress constraints in homogenization are triggered for a single constraint only. */
bool Function::DoEvaluate(const Excitation* excite) const {
  assert(excite != NULL);
  if(DoEvaluateAlways(excite->sequence))
    return true;

  if(excite_ == COMBINED_0_1_EX)
    return excite->index == 0 || excite->index == 1;

  return excite->index == excite_;
}

bool Function::DoEvaluateAlways(int context_sequence) const {
  if(excite_ != ALL_EX)
    return false;

  return ctxt->sequence == context_sequence; // excite_ == -1 is already assured
}


bool Function::IsSlackFunction() const
{
  switch(type_)
  {
  case SLACK:
  case SLACK_FNCT:
  case EXPRESSION:
    return true;
  default:
    return false;
  }
}

bool Function::IsStateDependent() const
{
  if(IsAdjointBased())
    return true;

  // only a few state dependent functions are not adjoint based

  // be careful: IsStateDependend() is not much used and might be incomplete!
  switch(type_)
  {
  case BANDGAP:
  case CONJUGATE_COMPLIANCE:
  case GLOBAL_DYNAMIC_COMPLIANCE:
  case ELEC_ENERGY:
  case COMPLIANCE:
  case HOM_TENSOR:
  case HOM_TRACKING:
  case HOM_FROBENIUS_PRODUCT:
  case POISSONS_RATIO:
  case YOUNGS_MODULUS:
  case YOUNGS_MODULUS_E1:
  case YOUNGS_MODULUS_E2:
  case TEMPERATURE:
  case HEAT_ENERGY:
  case EIGENFREQUENCY:
  case PRESSURE_DROP:
    return true;

  default:
    return false;
  }
}



bool Function::IsAdjointBased() const {
  switch (type_) {
  case TRACKING:
  case OUTPUT:
  case SQUARED_OUTPUT:
  case CONJUGATE_COMPLIANCE:
  case ABS_OUTPUT:
  case GLOBAL_DYNAMIC_COMPLIANCE:
  case DYNAMIC_OUTPUT:
  case REFLECTED_WAVE:
  case ELEC_ENERGY:
  case ENERGY_FLUX:
  case GLOBAL_STRESS:
  case LOCAL_STRESS:
  case TEMP_TRACKING_AT_INTERFACE:
  case SQR_MAG_FLUX_DENS_X:
  case SQR_MAG_FLUX_DENS_Y:
  case SQR_MAG_FLUX_DENS_RZ:
  case LOSS_MAG_FLUX_RZ:
  case MAG_COUPLING:
  case BUCKLING_LOAD_FACTOR:
  case LOCAL_BUCKLING_LOAD_FACTOR:
  case GLOBAL_BUCKLING_LOAD_FACTOR:
    return true;

  case COMPLIANCE: // only in the transient case
    return false; // FIXME check for transient case here!

  default:
    return false;
  }
}

bool Function::NeedsSelectionVector() const {
  switch (type_) {
  case OUTPUT:
  case SQUARED_OUTPUT:
  case CONJUGATE_COMPLIANCE:
  case ABS_OUTPUT:
  case DYNAMIC_OUTPUT:
  case REFLECTED_WAVE:
  case ENERGY_FLUX:
    return true;

  default:
    return false;
  }
}

bool Function::IsHomogenization() const {
  switch (type_) {
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

bool Function::CouldDoubleBounded(Type type)
{
  switch(type)
  {
  case SHAPE_INF:
    assert(false);
    // is it necessary to check for locality == Local::SHAPE?!
    return true;

  case SLOPE:
  case CURVATURE:
  case BENDING:
    return true;

  default:
    return false;
  }
  assert(false); // stupid compiler :(
  return false;
}

bool Function::IsDoubleBounded() const
{
  // could be extended to non-local functions!!
  return CouldDoubleBounded(type_) && local != NULL && !Local::IsReverse(local->GetLocality());
}

bool Function::IsLocal(Type t) {
  switch (t) {
  case SLOPE:
  case MOLE:
  case OSCILLATION:
  case JUMP:
  case BUMP:
  case CURVATURE:
  case OVERHANG_VERT:
  case OVERHANG_HOR:
  case DISTANCE:
  case BENDING:
  case CONES:
  case PERIODIC:
  case DESIGN:
  case SUM_MODULI:
  case TWO_SCALE_VOL:
  case ORTHOTROPIC_TENSOR_TRACE:
  case TENSOR_TRACE:
  case TENSOR_NORM:
  case PARAM_PS_POS_DEF:
  case POS_DEF_DET_MINOR_1:
  case POS_DEF_DET_MINOR_2:
  case POS_DEF_DET_MINOR_3:
  case BENSON_VANDERBEI_1:
  case BENSON_VANDERBEI_2:
  case BENSON_VANDERBEI_3:
  case MULTIMATERIAL_SUM:
  case SHAPE_INF:
  case LOCAL_STRESS:
  case LOCAL_BUCKLING_LOAD_FACTOR:
  case LOCAL_PYTHON_FUNCTION:
    return true;
  default:
    return false;
  }
}


void Function::SetElements(DesignSpace* space, RegionIdType region)
{
  assert(elements.GetSize() == 0);
  Grid* grid = domain->GetGrid();
  
  if(type_ == SHAPE_INF)
  {
    AuxDesign* aspace = dynamic_cast<AuxDesign*>(space);
    int n = space->GetNumberOfAuxParameters();
    elements.Reserve(n);
    for(int i = 0; i < n; i++)
      elements.Push_back(static_cast<DesignElement*>(aspace->GetAuxDesignElement(i)));
  }
  else
  {
    // Bastian's multiple design test cases have situations where design is DEFAULT as it is not
    // set in the objective
    // if ALL_REGIONS for condition use what we define as design space which
    // this is still not good enough
    int nd = 1;
    if(design_ == DesignElement::DEFAULT || design_ == DesignElement::ALL_DESIGNS)
      nd = space->design.GetSize();
    if(design_ == DesignElement::MECH_TRACE)
      nd = 6; // TODO why no 3?
    if(design_ == DesignElement::MECH_ALL)
      nd = 6;
    if(design_ == DesignElement::DIELEC_TRACE)
      nd = 2;
    if(design_ == DesignElement::DIELEC_ALL)
      nd = 2;
    if(design_ == DesignElement::PIEZO_ALL)
      nd = 6;
    //assert((int) space->design.GetSize() >= nd);

    elements.Reserve(nd * (region == ALL_REGIONS ? space->GetNumberOfElements() : grid->GetNumElems(region)));

    if (region == ALL_REGIONS || space->Contains(region)) {
      if (design_ == DesignElement::ALL_DESIGNS) {
        // FIXME - what the hell??? :(
        for (unsigned int i = 0; i < space->GetNumberOfElements(); i++) {
          DesignElement* de = &(space->data[i]);
          elements.Push_back(de);
        }
      } else
      {
        for(unsigned int i = 0; i < space->data.GetSize(); i++)
        {
          DesignElement* de = &(space->data[i]);
          if(DesignElement::IsCompatible(design_, de->GetType()) && (region == ALL_REGIONS || de->elem->regionId == region))
            elements.Push_back(de);
        }
      }
    } else {
      // this is a special case where the constraint does not act on the design space
      switch(type_)
      {
      case GLOBAL_STRESS:
      case LOCAL_STRESS:
      case LOCAL_BUCKLING_LOAD_FACTOR:
      case GLOBAL_BUCKLING_LOAD_FACTOR: // blubber
      case SQR_MAG_FLUX_DENS_X:
      case SQR_MAG_FLUX_DENS_Y:
      case SQR_MAG_FLUX_DENS_RZ:
      case MAG_COUPLING:
      case LOSS_MAG_FLUX_RZ:
        // do nothing
        break;
      default: {
          string a = grid->GetRegion().ToString(region);
          string msg = "region " + grid->GetRegion().ToString(region)
              + " of condition " + type.ToString(type_)
              + " not within design domain";
          preInfo_->SetWarning(msg);
          break;
        }
      }

      assert(elements.GetSize() == 0);

      // this creates the pseudo design elements and both indices are hopefully properly set!
      space->RegisterPseudoDesignRegion(region, design_, &elements);
    }
  }
  // empty elements for shape mapping!
  //  assert(elements.GetSize() == elements.GetCapacity());
}

void Function::SetDenseSparsityPattern(DesignSpace* space) {
  jac_sparsity_.Resize(space->GetNumberOfVariables()); // this might include aux variables

  for (unsigned int i = 0; i < jac_sparsity_.GetSize(); i++)
    jac_sparsity_[i] = i;
}

StdVector<unsigned int>& Function::GetSparsityPattern() {
  assert(!jac_sparsity_.IsEmpty()); // SetSparsityPattern() needs to be called before

  return jac_sparsity_;
}

Matrix<unsigned int>& Function::GetHessianSparsityPattern() {
  assert(hess_sparsity_.GetNumRows() == 0);

  return hess_sparsity_;
}

void Function::CalcHessian(StdVector<double>& out, double factor) {
  assert(out.GetSize() == 0);
  assert(GetHessianSparsityPattern().GetNumRows() == 0);
  // this function appears to have no Hessian, e.g. because it is linear :).
  assert(IsLinear());
  return;
}

void Function::PostProc(DesignSpace* space, DesignStructure* structure, ErsatzMaterial* em)
{
  // in the constructor we had to guess the filter type, now confirm
  LOG_DBG(func) << "PP: f_=" << Filter::type.ToString(filter_) << " sft=" << Filter::type.ToString(space->GetFilterType());
  if (domain->GetOptimization()->GetOptimizerType() != Optimization::SGP_SOLVER)
    // if SGP: cfs assembles the filters and passes them (without applying to the design) to external SGP solver
    // BUT: GuessFilter() returns density filter (which is ok) and GetFilterType() returns NO_FILTERING (need this to prevent cfs from applying the filter)
    // thus, in SGP case, this assert will fail
    assert(filter_ == space->GetFilterType());

  if(BaseDesignElement::IsShapeMapType(design_))
  {
    if(space->GetNumberOfFeatureMappingVariables() == 0)
      EXCEPTION("Function " << ToString() << " has shape mapping design type '" << BaseDesignElement::type.ToString(design_) << "' but 'ersatzMaterial@method' is not 'shapeMap'");
    if(!IsLocal(type_))
      EXCEPTION("Shape mapping design type " << BaseDesignElement::type.ToString(design_) << " for non-local function " << ToString());
  }

  // pre-init step
  switch (type_) {
  case SLOPE:
  case GLOBAL_SLOPE:
  case PERIMETER:
  case MOLE:
  case GLOBAL_MOLE:
  case OSCILLATION:
  case GLOBAL_OSCILLATION:
  case JUMP:
  case GLOBAL_JUMP:
  case BUMP:
  case CURVATURE:
  case GLOBAL_CURVATURE:
  case OVERHANG_VERT:
  case OVERHANG_HOR:
  case DISTANCE:
  case BENDING:
  case CONES:
  case PERIODIC:
    // assert(space->IsRegular()); // VicinityElements work only on a regular grid
    // the design elements require the vicinity element to be set which holds the direct
    // neighbors. Is save to call several times
    VicinityElement::Init(space, structure);
    InitLocal(space);
    break;

  case TWO_SCALE_VOL:
  case GLOBAL_TWO_SCALE_VOL:
  case ORTHOTROPIC_TENSOR_TRACE:
  case GLOBAL_ORTHOTROPIC_TENSOR_TRACE:
  case TENSOR_TRACE:
  case TENSOR_NORM:
  case GLOBAL_TENSOR_TRACE:
  case GLOBAL_SUM_MODULI:
  case PARAM_PS_POS_DEF:
  case DESIGN:
  case GLOBAL_DESIGN:
  case MULTIMATERIAL_SUM:
  case GLOBAL_STRESS:
  case LOCAL_STRESS:
  case LOCAL_BUCKLING_LOAD_FACTOR:
  case GLOBAL_BUCKLING_LOAD_FACTOR:
  case SUM_MODULI:
  case SHAPE_INF:
    // we need no neighbors.
    InitLocal(space);
    break;

  case PYTHON_FUNCTION:
  case LOCAL_PYTHON_FUNCTION:
    InitPythonFunction(pn, space);
    if(IsLocal())
      InitLocal(space);
    break;

  case PENALIZED_VOLUME:
    for(unsigned int i = 0; i < space->transfer.GetSize(); i++)
      if(space->transfer[i].IsPenalized())
        preInfo_->SetWarning("transfer function '" + space->transfer[i].ToString() + " seems also to penalize");
    break;

  case SLACK:
    if(!space->HasSlackVariable())
      throw Exception("'slack' as objective function requires 'slack' design");
    break;

  default: // do nothing
    break;
  }

  if(slackFnct_ != NO_FUNCTION && (!space->HasSlackVariable() || !space->HasAlphaVariable()))
    throw Exception("for slackFunction '" + slackFnct.ToString(slackFnct_) + "' designs 'slack' and 'alpha' are required");

  if(type_ == DISTANCE && design_ != DesignElement::NODE)
    throw Exception("for constraint 'distance' use design 'node'");

  if(type_ == BENDING && design_ != DesignElement::SPAGHETTI)
    throw Exception("for constraint 'bending' use design 'spaghetti'");

  // don't define the elements here, it is specific for objective and conditions

  if(design_ == DesignElement::DEFAULT && space->design.GetSize() > 1 && !IsStateDependent())
    if(!IsSlackFunction())
      preInfo_->Get(ParamNode::WARNING)->SetValue("consider setting 'design' for function '" + ToString() + "'");
}

double Function::CalcGraynessScaling() const
{
  DesignSpace* design = domain->GetDesign();
  assert(design != nullptr); // when called in Optimization constructor
  TransferFunction* tf = IsPhysical() ? design->GetTransferFunction(GetDesignType(), App::MECH) : nullptr;

  return tf ? (1.0/std::pow(tf->Transform(0.5),2)) : 4.0;
}

Function::Local* Function::InitLocal(DesignSpace* space)
{
  if (local == NULL)
  {
    local = new Local(this, space);
    local->PostInit();
  }
  return local;
}

Function::Local::Local(Function* func, DesignSpace* space)
{
  this->space = space;
  this->func_ = func;
  this->structure_ = NULL;
  this->element_dimension_ = -1;

  // shortcuts
  Function::Type ftype = func->GetType();
  string fname = Function::type.ToString(ftype);

  // read xml parameters -> might be null valued!
  PtrParamNode pn = func->pn->Get("local", ParamNode::PASS);

  this->beta_ = pn != NULL && pn->Has("beta") ? pn->Get("beta")->As<double>() : -3.14;
  this->eps_ = pn != NULL && pn->Has("eps") ? pn->Get("eps")->As<double>() : -3.14;
  this->power_ = pn != NULL && pn->Has("power") ? pn->Get("power")->As<double>() : 2.0;
  this->phase_ = pn != NULL && pn->Has("phase") ? phase.Parse(pn->Get("phase")->As<string>()) : BOTH; // only oscillation

  this->normalize_ = pn != NULL ? pn->Get("normalize")->As<bool>() : false;

  bool enable = pn != NULL ? pn->Get("periodic")->As<bool>() : true; // enable/disable is handled in As<bool>() and true is default
  this->periodic = enable & domain->HasPerdiodicBC();
  int dim = domain->GetDim();
  if (dim == 3 && (domain->GetParamRoot()->Has("optimization/ersatzMaterial/paramMat/designMaterials/designMaterial/homRectC1") || domain->GetParamRoot()->Has("optimization/ersatzMaterial/paramMat/designMaterials/designMaterial/homIsoC1") || domain->GetParamRoot()->Has("optimization/ersatzMaterial/paramMat/designMaterials/designMaterial/heat"))) {
    //read interpolation data for volume calculation in 3D
    PtrParamNode dm_node = domain->GetParamRoot()->Get("optimization/ersatzMaterial/paramMat/designMaterials")->GetByVal("designMaterial", "sequence", Optimization::context->sequence);
    assert(dm_node != NULL);
    std::string material_type = "";
    DesignMaterial::Type dtype = space->GetDesignMaterialType();
    if (dtype == DesignMaterial::HOM_RECT_C1)
      material_type = "homRectC1";
    if (dtype == DesignMaterial::HOM_ISO_C1)
      material_type = "homIsoC1";
    if (dtype == DesignMaterial::HEAT)
      material_type = "heat";

    std::string file = dm_node->Get(material_type)->Get("file")->As<std::string>();
    PtrParamNode root = XmlReader::ParseFile(file);
    ParamTools::AsVector<double>(root->Get("param1/matrix/real"), this->vol_a_);
    if (root->Has("param2")) {
      ParamTools::AsVector<double>(root->Get("param2/matrix/real"), this->vol_b_);
    } else {
      vol_b_.Resize(0);
      vol_b_.Init();
    }
    if (root->Has("param3")) {
      ParamTools::AsVector<double>(root->Get("param3/matrix/real"), this->vol_c_);
    } else {
      vol_c_.Resize(0);
      vol_c_.Init();
    }
    if (root->Has("coeffvol")) {
      int dim1 = root->Get("coeffvol/matrix/dim1")->As<int>();
      int dim2 = root->Get("coeffvol/matrix/dim2")->As<int>();
      ParamTools::AsTensor<double>(root->Get("coeffvol/matrix/real"), dim1, dim2, this->vol_coeff_);

      // copy Vector to StdVector
      StdVector<double> a(this->vol_a_.GetSize()), b(this->vol_b_.GetSize()), c(this->vol_c_.GetSize());
      std::copy(this->vol_a_.GetPointer(), this->vol_a_.GetPointer() + this->vol_a_.GetSize(), a.Begin());
      std::copy(this->vol_b_.GetPointer(), this->vol_b_.GetPointer() + this->vol_b_.GetSize(), b.Begin());
      std::copy(this->vol_c_.GetPointer(), this->vol_c_.GetPointer() + this->vol_c_.GetSize(), c.Begin());

      // create interpolator for volume
      StdVector<double> dummy_data;
      if (c.GetSize() > 0) {
        assert(a.GetSize() > 0 && b.GetSize() > 0);
        dummy_data.Resize(a.GetSize()*b.GetSize()*c.GetSize());
        volumeInterpolator_ = new TriCubicInterpolate(dummy_data, a, b, c, this->vol_coeff_);
      } else if (b.GetSize() > 0) {
        assert(a.GetSize() > 0);
        dummy_data.Resize(a.GetSize()*b.GetSize());
        volumeInterpolator_ = new BiCubicInterpolate(dummy_data, a, b, this->vol_coeff_);
      } else {
        assert(a.GetSize() > 0);
        dummy_data.Resize(a.GetSize());
        volumeInterpolator_ = new CubicInterpolate(dummy_data, a, this->vol_coeff_);
      }
    }
  }

  //total volume in the non-regular case is needed for the volume calculations
  this->total_vol_ = 0.0;

  if(!space->IsCubic())
    for (unsigned int i = 0, n = this->func_->elements.GetSize(); i < n;i++)
     this->total_vol_ += this->func_->elements[i]->CalcVolume();
  else
   this->total_vol_ = 1.0;

  switch (ftype)
  {
  // the PERIMETER is not globalized in the sense of the sum max(f,0)^p function.
  case GLOBAL_JUMP:
  case GLOBAL_MOLE:
  case GLOBAL_OSCILLATION:
  case GLOBAL_SLOPE:
  case GLOBAL_DESIGN:
  case GLOBAL_STRESS:
  case GLOBAL_BUCKLING_LOAD_FACTOR:
    this->globalized_ = true;
    break;

  case GLOBAL_SUM_MODULI:
  case GLOBAL_TWO_SCALE_VOL:
  case GLOBAL_ORTHOTROPIC_TENSOR_TRACE:
  case GLOBAL_TENSOR_TRACE:
    if (power_ != 1.0)
      domain->GetInfoRoot()->Get("optimization/header")->SetWarning("function '" + fname + "' has local/power " + lexical_cast<string>(power_) + ", for sum one needs power=1");
    if (!normalize_)
      domain->GetInfoRoot()->Get("optimization/header")->SetWarning("function '" + fname + "' should be normalized");
    if (!this->func_->IsObjective()) {
      double parameter = dynamic_cast<Condition*>(this->func_)->GetParameter();
      if (parameter > 0)
        domain->GetInfoRoot()->Get("optimization/header")->SetWarning("'" + fname + "' computes as sum(max(f_i-" + lexical_cast<string>(parameter) + ",0). Are you sure about parameter=" + lexical_cast<string>(parameter) + "?");
    }
    this->globalized_ = true;
    break;

  default:
    this->globalized_ = false;
    break;
  }

  // check beta
  if(RequiresBeta(ftype) && (pn == NULL || !pn->Has("beta")))
    throw Exception("function '" + fname + "' requires the 'beta' attribute in a 'local' element");
  if(RequiresBeta(ftype) && (func->IsObjective() || dynamic_cast<Condition*>(func)->IsActive()) && beta_ < 0)
      throw Exception("'function '" + fname + "' allows beta=-1 only for condition in observe mode");

  // check eps
  if(RequiresEps(ftype) && (pn == NULL || !pn->Has("eps")))
    throw Exception("function '" + fname + "' requires the 'eps' attribute in a 'local' element");

  // check phase
  if (phase_ != BOTH && ftype != OSCILLATION && ftype != GLOBAL_OSCILLATION)
    throw Exception("'phase' may only be set for (global) oscillation");

  // set locality
  this->locality_ = pn != NULL && pn->Has("locality") ? locality.Parse(pn->Get("locality")->As<string>()) : DEFAULT;
  Locality user = locality_; // default or set by user

  // this is our only double bounded optimizer
  bool db_opt = domain->GetParamRoot()->Get("optimization/optimizer/type")->As<string>() == "snopt";

  switch (ftype) {

  case SLOPE:
    if(user == DEFAULT)
      locality_ = db_opt ? NEXT : NEXT_AND_REVERSE;
    if(locality_ != NEXT && locality_ != NEXT_AND_REVERSE)
      throw Exception("Invalid locality '" + locality.ToString(locality_) + "' within '" + fname + "'");
    break;

  // the horizontal overhang is | ... | >= c^* and because of the >= we cannot bound as with slopes but need to smooth abs
  case OVERHANG_HOR:
    if(locality_ != MULT_DESIGNS_NEXT && locality_ != DEFAULT)
      throw Exception("Invalid locality '" + locality.ToString(locality_) + "' within '" + fname + "'");
    locality_ = MULT_DESIGNS_NEXT;
    break;

  // opposite to OVERHANG_VERT there is always reverse and no bounded stuff as in OVERHANG_HOR as in SLOPE
  case OVERHANG_VERT:
    if(locality_ != MULT_DESIGNS_NEXT_AND_REVERSE && locality_ != DEFAULT)
      throw Exception("Invalid locality '" + locality.ToString(locality_) + "' within '" + fname + "'");
    locality_ = MULT_DESIGNS_NEXT_AND_REVERSE;
    break;

  case DISTANCE:
    if(locality_ != FUNCTION_SPECIFIC && locality_ != DEFAULT)
      throw Exception("Invalid locality '" + locality.ToString(locality_) + "' within '" + fname + "'");
    locality_ = FUNCTION_SPECIFIC;
    break;

  case BENDING:
    if(user == DEFAULT)
      locality_ = db_opt ? FUNCTION_SPECIFIC : FUNCTION_SPECIFIC_TWO_SIGNS;
    if(locality_ != FUNCTION_SPECIFIC && locality_ != FUNCTION_SPECIFIC_TWO_SIGNS && locality_ != DEFAULT)
      throw Exception("Invalid locality '" + locality.ToString(locality_) + "' within '" + fname + "'");
    break;

  case CONES:
    if(locality_ != MULT_DESIGNS_PREV_NEXT && locality_ != DEFAULT)
      throw Exception("Invalid locality '" + locality.ToString(locality_) + "' within '" + fname + "'");
    locality_ = MULT_DESIGNS_PREV_NEXT;
    break;

  case GLOBAL_SLOPE:
    if (locality_ != NEXT && locality_ != DEFAULT)
      throw Exception("Invalid locality '" + locality.ToString(locality_) + "' within '" + fname + "'");
    locality_ = NEXT_AND_REVERSE;
    break;

  case CURVATURE:
    if(user == DEFAULT)
      locality_ = db_opt ? PREV_NEXT : PREV_NEXT_AND_REVERSE;
    if (locality_ != PREV_NEXT && locality_ != PREV_NEXT_AND_REVERSE && locality_ != DEFAULT)
      throw Exception("Invalid locality '" + locality.ToString(locality_) + "' within '" + fname + "'");
    break;

  case GLOBAL_CURVATURE:
    if (locality_ != PREV_NEXT && locality_ != DEFAULT)
      throw Exception("Invalid locality '" + locality.ToString(locality_) + "' within '" + fname + "'");
    locality_ = PREV_NEXT_AND_REVERSE;
    break;

  case PERIMETER:
    if (locality_ != NEXT && locality_ != DEFAULT)
      throw Exception("Invalid locality '" + locality.ToString(locality_) + "' within '" + fname + "'");
    locality_ = NEXT;
    break;

  case PERIODIC:
    if (locality_ != CYCLIC && locality_ != DEFAULT)
      throw Exception("Invalid locality '" + locality.ToString(locality_) + "' within '" + fname + "'");
    locality_ = CYCLIC;
    break;

  case OSCILLATION:
  case GLOBAL_OSCILLATION:
    if ((phase_ == BOTH && locality_ != DEG_45_STAR_AND_REVERSE
        && locality_ != PREV_NEXT_AND_REVERSE && locality_ != DEFAULT)
        || (phase_ != BOTH && locality_ != DEG_45_STAR && locality_ != PREV_NEXT  && locality_ != DEFAULT))
      throw Exception("Invalid locality '" + locality.ToString(locality_) + "' within '" + fname + "' and phase '" + phase.ToString(phase_) + "'");
    if (locality_ == DEFAULT)
      locality_ = phase_ == BOTH ? DEG_45_STAR_AND_REVERSE : DEG_45_STAR;
    break;

  case MOLE:
  case GLOBAL_MOLE:
    if (locality_ != DEG_45_STAR && locality_ != DEFAULT)
      throw Exception("Invalid locality '" + locality.ToString(locality_) + "' within '" + fname + "'");
    locality_ = DEG_45_STAR;
    break;

  case JUMP:
  case GLOBAL_JUMP:
    if (locality_ != BOUNDARY && locality_ != DEFAULT)
      throw Exception("Invalid locality '" + locality.ToString(locality_) + "' within '" + fname + "'");
    locality_ = BOUNDARY;
    break;

  case GLOBAL_STRESS:
  case LOCAL_STRESS:
  case LOCAL_BUCKLING_LOAD_FACTOR:
  case GLOBAL_BUCKLING_LOAD_FACTOR:
  case DESIGN:
  case GLOBAL_DESIGN:
    if (locality_ != ELEMENT && locality_ != DEFAULT)
      throw Exception("Invalid locality '" + locality.ToString(locality_) + "' within '" + fname + "'");
    locality_ = ELEMENT;
    break;

  case TENSOR_TRACE:
  case TENSOR_NORM:
  case GLOBAL_TENSOR_TRACE:
  case ORTHOTROPIC_TENSOR_TRACE:
  case GLOBAL_ORTHOTROPIC_TENSOR_TRACE:
  case SUM_MODULI:
  case GLOBAL_SUM_MODULI:
  case TWO_SCALE_VOL:
  case GLOBAL_TWO_SCALE_VOL:
  case PARAM_PS_POS_DEF:
  case POS_DEF_DET_MINOR_1:
  case POS_DEF_DET_MINOR_2:
  case POS_DEF_DET_MINOR_3:
  case BENSON_VANDERBEI_1:
  case BENSON_VANDERBEI_2:
  case BENSON_VANDERBEI_3:
  case MULTIMATERIAL_SUM:
    if (locality_ != MULT_DESIGNS_ELEMENT && locality_ != DEFAULT)
      throw Exception("Invalid locality '" + locality.ToString(locality_) + "' within '" + fname + "'");
    locality_ = MULT_DESIGNS_ELEMENT;
    break;

  case BUMP:
    if (locality_ != PREV_NEXT && locality_ != DEFAULT)
      throw Exception("Invalid locality '" + locality.ToString(locality_) + "' within '" + fname + "'");
    locality_ = PREV_NEXT;
    break;

  case LOCAL_PYTHON_FUNCTION:
    if (locality_ != EXTERNALLY_DEFINED && locality_ != DEFAULT)
      throw Exception("Invalid locality '" + locality.ToString(locality_) + "' within '" + fname + "'");
    locality_ = EXTERNALLY_DEFINED;
    break;

    
  case SHAPE_INF:
    locality_ = SHAPE;
    break;

  default: // no locality
    assert(false);
    break;
  }

  if(CouldDoubleBounded(ftype) && !db_opt && !Function::Local::IsReverse(locality_))
   throw Exception("The optimizer has no bounds for constraints: your choice for 'local' is invalid");

}

void Function::Local::PostInit()
{
  PtrParamNode pn = func_->pn->Get("local", ParamNode::PASS);
  Function::Type ftype = func_->GetType();
  string fname = Function::type.ToString(ftype);

  // this is actually pure constructor work, just extracted to handle function size
  FeaturedDesign* fd = dynamic_cast<FeaturedDesign*>(space); // only not null if we do not shape mapping or other feature mapping
  assert(!(BaseDesignElement::IsShapeMapType(func_->GetDesignType()) && space->GetNumberOfFeatureMappingVariables() == 0));

  switch (locality_)
  {
  case DEG_45_STAR:
  case DEG_45_STAR_AND_REVERSE:
  case BOUNDARY:
    if (!pn)
      throw Exception("sub element 'local' with neighborhood information mandatory for '" + fname + "'");
    structure_ = new NeighborhoodStructure(this, pn);
    if (locality_ == BOUNDARY)
      SetupBoundaryElementMap();
    else
      SetupStarLocalityElementMap(phase_);
    break;

  case ELEMENT:
    SetupSingularElementMap();
    break;

  case MULT_DESIGNS_ELEMENT:
    SetupMultDesignsElementMap(func_);
    break;
    
  case SHAPE:
    SetupShapeElementMap(func_, dynamic_cast<ShapeDesign*>(space));
    break;

  case MULT_DESIGNS_NEXT:
  case MULT_DESIGNS_PREV_NEXT:
  case MULT_DESIGNS_NEXT_AND_REVERSE:
  case MULT_DESIGNS_PREV_NEXT_AND_REVERSE:
    if(BaseDesignElement::IsFeatureMappingType(func_->GetDesignType()))
      fd->SetupVirtualMultiShapeElementMap(func_, virtual_elem_map, locality_);
    else
      SetupMultDesignsVirtualElementMap(func_);
    break;

  case NEXT_DIAG:
   // if(!pn)
   // throw Exception("sub element 'local' with neighborhood information mandatory for '" + fname + "'");
   //    structure_ = new NeighborhoodStructure(this, pn);
   SetupVirtualStarLocalElementMap(func_);
   break;

  case CYCLIC:
    if(BaseDesignElement::IsShapeMapType(func_->GetDesignType()))
      fd->SetupCyclicVirtualShapeElementMap(func_, virtual_elem_map, locality_);
    else
      throw Exception("the local function '" + func_->ToString() + "' is only mapping");
   break;

  case EXTERNALLY_DEFINED:
    assert(ftype == LOCAL_PYTHON_FUNCTION);
    func_->SetLocalPythonVirtualElementMap(virtual_elem_map, space);
    break;


  default:
    if(BaseDesignElement::IsFeatureMappingType(func_->GetDesignType()))
      fd->SetupVirtualShapeElementMap(func_, virtual_elem_map, locality_);
   else
      SetupVirtualElementMap(phase_);
    break;
  }

  if (virtual_elem_map.GetSize() == 0)
    throw Exception("mesh too small for locality of function '" + fname + "' or wrong design attribute");

  // needs to be set prior CalcSlopeConstraint() as the optimizers need the size
  local_values.Resize(virtual_elem_map.GetSize(), -1.0);
}


Function::Local::~Local() {
  if (structure_ != NULL) {
    delete structure_;
    structure_ = NULL;
  }
}

void Function::Local::SetupVirtualElementMap(Phase ph) {
  // we construct locality_ into reverse, prev and next
  // reverse means we have a REVERSE option which makes two constraints with different signs
  int dim = domain->GetDim();
  bool prev = locality_ == PREV_NEXT_AND_REVERSE || locality_ == PREV_NEXT;
  bool next = true; // always
  bool two_signs = locality_ == NEXT_AND_REVERSE || locality_ == PREV_NEXT_AND_REVERSE;

   //assert((ph == BOTH && two_signs) || (!two_signs && ph != BOTH));
  // assume ph is set correctly and Phase is in sync with the signs
  int sign_1 = ph != BOTH ? (int) ph : two_signs ? 1 : Identifier::NO_SIGN;
  int sign_2 = ph != BOTH ? (int) ph : -1;

  element_dimension_ = dim * (two_signs ? 2 : 1);

  virtual_elem_map.Reserve(element_dimension_ * space->GetNumberOfElements());

  // traverse all elements and check for full neighborhood
  for (int e = 0, ss = func_->elements.GetSize(); e < ss; ++e)
  {
    DesignElement* de = func_->elements[e];
    if (de->GetType()== ((func_->design_ == DesignElement::DEFAULT) ? DesignElement::DENSITY : func_->design_))
    {
      VicinityElement* ve = de->vicinity;

      for(int a = 0; a < dim; a++)
      {
        // do we have a full neighborhood? All or none as in the original slope paper
        bool full = true;
        if(prev)
          if(ve->design[2*a + 1] == NULL) full = false; // X_N=1, Y_N=3, Z_N=5
        if(next)
          if(ve->design[2*a] == NULL) full = false; // X_P=0, Y_P=2, Z_P=4
        LOG_DBG2(func) << "Local::Local e_num=" << de->elem->elemNum << " vicinity=" << ve->ToString() << " full=" << full;

        if(full)
        {
          DesignElement* prev_de = prev ? ve->GetNeighbour(VicinityElement::ToNeighbour(a, -1)) : NULL;
          DesignElement* next_de = ve->GetNeighbour(VicinityElement::ToNeighbour(a, 1));

          virtual_elem_map.Push_back(Identifier(de, prev_de, next_de, sign_1));
          if (two_signs)
            virtual_elem_map.Push_back(Identifier(de, prev_de, next_de, sign_2));
        }
      }
    }
  }
}

void Function::Local::SetupVirtualStarLocalElementMap(const Function* f)
{
  // we construct locality_ into reverse, prev and next
  // reverse means we have a REVERSE option which makes two constraints with different signs
  assert(locality_==NEXT_DIAG); // It is not clean, but I do it by hand to be sure of what I am doing

//  int  dim  = domain->GetDim();

  //So far works only for d=2

 // assert (dim == 2);

  // only this element!
  element_dimension_ = 1; // two boundary "stones" per dimension

  UInt elems = space->GetNumberOfElements();
  virtual_elem_map.Reserve(element_dimension_ * elems);

  assert(space->design.GetSize()== 2);
  if(f->GetDesignType() != DesignElement::ALL_DESIGNS)
    throw Exception("'tensor_norm' only defined 'ALL_DESIGNS' design");

  for(unsigned int e = 0; e < elems; e++)
   {
     DesignElement* de = func_->elements[e]; // -> GX_0
     assert((int) e == space->Find(de->elem, true)); // assert that we still are on the right finite element

       VicinityElement* ve = de->vicinity;

      // do we have a full neighborhood? All or none as in the original slope paper
      bool full = true;
        if(ve->design[VicinityElement::X_P] == NULL) full = false;
        if(ve->design[VicinityElement::Y_P] == NULL) full = false;
        if (full)
        {
        VicinityElement* ve_px = ve->design[VicinityElement::X_P]->vicinity;
          if(ve_px->design[VicinityElement::Y_P] == NULL) full = false;
        }


      LOG_DBG2(func) << "Local::Local e_num=" << de->elem->elemNum << " vicinity=" << ve->ToString() << " full=" << full;

      if(full)
      {

          DesignElement* de_px = ve->GetNeighbour(VicinityElement::X_P); //-GX_PX
          DesignElement* de_py = ve->GetNeighbour(VicinityElement::Y_P); //GX_PY
          VicinityElement* ve_px = ve->design[VicinityElement::X_P]->vicinity;
          DesignElement* de_pxy = ve_px->GetNeighbour(VicinityElement::Y_P); //GX_PXY

          DesignElement* other = &(space->data[elems * 1 + e]); //GY_0

          VicinityElement* vother = other->vicinity;

          DesignElement* other_px = vother->GetNeighbour(VicinityElement::X_P); //-GY_PX
          DesignElement* other_py = vother->GetNeighbour(VicinityElement::Y_P); //GY_PY
          VicinityElement* vother_px = vother->design[VicinityElement::X_P]->vicinity;
          DesignElement* other_pxy = vother_px->GetNeighbour(VicinityElement::Y_P); //GY_PXY
          //Add other

          StdVector<BaseDesignElement*> buddies;
          buddies.Push_back(de_px);
          buddies.Push_back(de_py);
          buddies.Push_back(de_pxy);

          buddies.Push_back(other);
          buddies.Push_back(other_px);
          buddies.Push_back(other_py);
          buddies.Push_back(other_pxy);
          // Identifier(BaseDesignElement* elem, StdVector<BaseDesignElement*> buddies, int si = NO_SIGN);
          virtual_elem_map.Push_back(Identifier(de, buddies, 1));

          // So in the order, we have:
          //GX_0(-1)  GX_PX(0)  GX_PY(1) GX_PXY(2) GY_0(3) GY_PX(4) GY_PY(5) GY_PXY(6)
      }
    }
}

void Function::Local::SetupStarLocalityElementMap(Phase ph) {
  unsigned int dim = domain->GetDim();
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
  for (int e = 0, ss = elems; e < ss; ++e)
  {
    DesignElement* de = &(space->data[e]);

    // do we have a full neighborhood? All or none
    bool full = true;

    // we only need to check orthogonal
    assert(VicinityElement::X_P == 0);
    for (int dir = VicinityElement::X_P; dir <= (dim == 2 ? VicinityElement::Y_N : VicinityElement::Z_N); dir++)
    {
      int a = VicinityElement::ToMainAxis((VicinityElement::Neighbour) dir);
      unsigned int n = struc->orthogonal[a];
      assert(n > 0);
      if (!VicinityElement::HasNeighbor(de, (VicinityElement::Neighbour) dir, n))
        full = false;
    }

    if (!full)
      continue;
    // orthogonal first
    StdVector<BaseDesignElement*> buddies;

    for (int dir = VicinityElement::X_P; dir <= (dim == 2 ? VicinityElement::Y_N : VicinityElement::Z_N); dir += 2)
    {
      VicinityElement::Neighbour pos = (VicinityElement::Neighbour) dir;
      VicinityElement::Neighbour neg = (VicinityElement::Neighbour) (dir + 1);
      unsigned int a = VicinityElement::ToMainAxis(pos);
      unsigned int n = struc->orthogonal[a];

      buddies.Resize(0);
      for (unsigned int i = n; i > 0; i--)
        buddies.Push_back(VicinityElement::GetNeighbour(de, neg, i));
      for (unsigned int i = 1; i <= n; i++)
        buddies.Push_back(VicinityElement::GetNeighbour(de, pos, i));

      LOG_DBG3(func) << "L:SSLEM: de=" << de->ToString() << " dir=" << dir << " pos=" << pos << " neg=" << neg << " a=" << a << " n=" << n << " buddies=" << BaseDesignElement::ToString(buddies);

      virtual_elem_map.Push_back(Identifier(de, buddies, sign_1));
      if (two_signs)
        virtual_elem_map.Push_back(Identifier(de, buddies, sign_2));
    }

    // In 2D we have 2 diagonals in the xy plane. In 3D also in the xz and the yz plane which makes
    // 6 diagonals in the planes plus 4 diagonals between the 8 corners. -> 10 diagonals in 3D

    // we make the assumption, that our elements are almost quadratic, such our diagonals are
    // simply alternating the directions. The total number of diagonals (we want to describe a
    // ball and not a cube) is determined by the ball radius and the element diagonals.
    // start with plane diagonals
    for (unsigned int d = 0; d < (dim == 2 ? 1 : 3); d++) {
      unsigned int n = struc->diagonal[d];
      // within a plane we go zig-zag. First the first axis then the second.
      // Plus first in negative direction then positive
      // the planes are xy, xy and yz
      int axis_first = d != 2 ? 0 : 1; // xy, xz, yz
      int axis_second = d == 0 ? 1 : 2; // xy, xz, yz

      // for the first  direction we have the X_N/Y_N elements and the X_P/Y_P elements.
      // for the second direction we have the X_N/Y_P elements and the X_P/Y_N elements.
      for (int dir = -1; dir <= 1; dir += 2) {
        buddies.Resize(0);
        for (int e = n; e > 0; e--) {
          DesignElement* tmp = VicinityElement::GetNeighbour(de, VicinityElement::ToNeighbour(axis_first, -1), e);
          buddies.Push_back(VicinityElement::GetNeighbour(tmp, VicinityElement::ToNeighbour(axis_second, dir), e));
        }
        for (unsigned int e = 1; e <= n; e++) {
          DesignElement* tmp = VicinityElement::GetNeighbour(de, VicinityElement::ToNeighbour(axis_first, 1), e);
          buddies.Push_back(VicinityElement::GetNeighbour(tmp, VicinityElement::ToNeighbour(axis_second, dir == 1 ? -1 : 1), e));
        }

        // LOG_DBG3(func) << "L:SSLEM: diag de=" << de->ToString() << " dir=" << dir << " buddies=" << BaseDesignElement::ToString(buddies);

        virtual_elem_map.Push_back(Identifier(de, buddies, sign_1));
        if (two_signs)
          virtual_elem_map.Push_back(Identifier(de, buddies, sign_2));
      }
    }
    if (dim == 3) {
      // now the four corners
      // -x -y -z -> +x +y +z
      // -x -y +z -> +x +y -z
      // -x +y -z -> +x -y +z
      // -x +y +z -> +x -y -z
      assert(struc->diagonal.GetSize() == 4);
      unsigned int n = struc->diagonal[3];
      for (int dir_y = -1; dir_y <= 1; dir_y += 2) {
        for (int dir_z = -1; dir_z <= 1; dir_z += 2) {
          buddies.Resize(0);
          for (int e = n; e > 0; e--) {
            DesignElement* tmp_x = VicinityElement::GetNeighbour(de, VicinityElement::X_N, e);
            DesignElement* tmp_y = VicinityElement::GetNeighbour(tmp_x, VicinityElement::ToNeighbour(1, dir_y), e);
            buddies.Push_back(VicinityElement::GetNeighbour(tmp_y, VicinityElement::ToNeighbour(2, dir_z), e));
          }
          for (unsigned int e = 1; e <= n; e++) {
            DesignElement* tmp_x = VicinityElement::GetNeighbour(de, VicinityElement::X_P, e);
            DesignElement* tmp_y = VicinityElement::GetNeighbour(tmp_x, VicinityElement::ToNeighbour(1, dir_y == 1 ? -1 : 1), e);
            buddies.Push_back(VicinityElement::GetNeighbour(tmp_y, VicinityElement::ToNeighbour(2, dir_z == 1 ? -1 : 1), e));
          }

          LOG_DBG3(func) << "L:SSLEM: corner de=" << de->ToString() << " dir_y=" << dir_y << " dir_z=" << dir_z << " buddies=" << BaseDesignElement::ToString(buddies);

          virtual_elem_map.Push_back(Identifier(de, buddies, sign_1));
          if (two_signs)
            virtual_elem_map.Push_back(Identifier(de, buddies, sign_2));
        }
      }
    }
  }
}

void Function::Local::SetupBoundaryElementMap() {
  unsigned int dim = domain->GetDim();
  // oscillation has with BOTH DEG_45_STAR_AND_REVERSE,
  // mole has always BOTH and DEG_45_STAR.
  // oscillation w/o BOTH needs to be DEG_45_STAR
  assert(structure_ != NULL);
  NeighborhoodStructure* struc = structure_;

  element_dimension_ = dim * 2; // two boundary "stones" per dimension
  virtual_elem_map.Reserve(element_dimension_ * space->GetNumberOfElements());

  space->AssertOneDesignOnly(); // can be extended we use the design from the condition
  int elems = space->GetNumberOfElements();
  for (int e = 0, ss = elems; e < ss; ++e) {
    DesignElement* de = &(space->data[e]);

    // do we have a full neighborhood? All or none
    bool full = true;

    for (int dir = VicinityElement::X_P; dir <= (dim == 2 ? VicinityElement::Y_N : VicinityElement::Z_N); dir++) {
      int a = VicinityElement::ToMainAxis((VicinityElement::Neighbour) dir);
      unsigned int n = struc->orthogonal[a];
      assert(n > 0);
      if (!VicinityElement::HasNeighbor(de, (VicinityElement::Neighbour) dir, n))
        full = false;
    }

    if (!full)
      continue;

    for (int dir = VicinityElement::X_P; dir <= (dim == 2 ? VicinityElement::Y_N : VicinityElement::Z_N); dir += 2) {
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

void Function::Local::SetupSingularElementMap() {
  // only this element!
  element_dimension_ = 1; // two boundary "stones" per dimension
  virtual_elem_map.Reserve(element_dimension_ * func_->elements.GetSize());

  StdVector<BaseDesignElement*> empty;

  for (int e = 0, en = func_->elements.GetSize(); e < en; e++) {
    DesignElement* de = func_->elements[e];
    virtual_elem_map.Push_back(Identifier(de, empty));
  }
}

void Function::Local::SetupMultDesignsElementMap(const Function* f) {
  // only this element!
  element_dimension_ = 1; // two boundary "stones" per dimension
  UInt elems = space->GetNumberOfElements();
  virtual_elem_map.Reserve(element_dimension_ * elems);

  // the neighbors are the design elements for the same FE-element but with other designs
  // one is not a neighbor of oneself
  StdVector<BaseDesignElement*> neighbours;

  StdVector<unsigned int> des_idx; // the design indices we consider here
  switch (f->GetType()) {
  case TENSOR_TRACE: {
    // we have the case of elast and piezo fmo with elec tensor
    // but also param mat where stiff1 and stiff2, ... form the tensor
    // The problem is, that the first *compatible* design shall not be within the neighborhood,
    // but for piezo-FMO the first element might not be a compatible one
    bool first_compatible_found = false;
    for (unsigned int i = 0; i < space->design.GetSize(); ++i) {
      if (DesignElement::IsCompatible(f->GetDesignType(), space->design[i].design)) {
        if (!first_compatible_found)
          first_compatible_found = true; // we do not add the first compatible designs to the neighbors
        else
          des_idx.Push_back(i);
      }
    }
    break;
  }
  case TENSOR_NORM:
    if (f->GetDesignType() != DesignElement::PIEZO_ALL)
      throw Exception("'tensor_norm' only defined for 'piezo_all' design");
    des_idx.Push_back(space->FindDesign(DesignElement::PIEZO_12));
    des_idx.Push_back(space->FindDesign(DesignElement::PIEZO_13));
    des_idx.Push_back(space->FindDesign(DesignElement::PIEZO_21));
    des_idx.Push_back(space->FindDesign(DesignElement::PIEZO_22));
    des_idx.Push_back(space->FindDesign(DesignElement::PIEZO_23));
    break;

  case POS_DEF_DET_MINOR_1:
  case BENSON_VANDERBEI_1:
    assert(space->design.GetSize() >= 6);
    if (space->design[0].design != DesignElement::MECH_11
        || space->design[0].design != DesignElement::DIELEC_11)
      throw Exception("'Expect first design to be 'tensor11' or 'dielec_11'");
    // only design TENSOR11 is not neighbor
    break;
  case POS_DEF_DET_MINOR_2:
  case BENSON_VANDERBEI_2:
    assert(space->design.GetSize() >= 6);
    if (f->GetDesignType() == DesignElement::DIELEC_ALL) {
      des_idx.Push_back(space->FindDesign(DesignElement::DIELEC_22));
      des_idx.Push_back(space->FindDesign(DesignElement::DIELEC_12));
    } else {
      des_idx.Push_back(space->FindDesign(DesignElement::MECH_22));
      des_idx.Push_back(space->FindDesign(DesignElement::MECH_12));
    }
    break;
  case POS_DEF_DET_MINOR_3:
  case BENSON_VANDERBEI_3:
    assert(space->design.GetSize() >= 6);
    // note, that the indices are sorted in sparse pattern
    des_idx.Push_back(space->FindDesign(DesignElement::MECH_22));
    des_idx.Push_back(space->FindDesign(DesignElement::MECH_33));
    des_idx.Push_back(space->FindDesign(DesignElement::MECH_23));
    des_idx.Push_back(space->FindDesign(DesignElement::MECH_13));
    des_idx.Push_back(space->FindDesign(DesignElement::MECH_12));
    break;
  default:
    // all designs but the first one
    des_idx.Reserve(space->design.GetSize() - 1);
    for (unsigned int i = 1; i < space->design.GetSize(); i++)
      des_idx.Push_back(i);
    break;
  }

  // LOG_DBG(func)<< "F:L:SMDEM des_idx=" << des_idx.ToString() << " total=" << space->design.ToString();


  if(elems > func_->elements.GetSize())
    EXCEPTION("for function " << Function::type.ToString(f->GetType()) << " with design " << DesignElement::type.ToString(f->GetDesignType())
              << " the design space is " << func_->elements.GetSize() << " but should be " << elems);

  for(unsigned int e = 0; e < elems; e++)
  {
    DesignElement* de = func_->elements[e];
    //assert((int ) e == space->Find(de->elem, true)); // assert that we still are on the right finite element

    neighbours.Resize(0);

    for (unsigned int d = 0; d < des_idx.GetSize(); d++) {
      // the first design = 0 is no neighbor!
      unsigned des = des_idx[d];
      DesignElement* other = &(space->data[elems * des + e]);
      neighbours.Push_back(other);
      LOG_DBG3(func) << "F:L:SMDEM e=" << e << " el=" << de->elem->elemNum << " d = " << d << " des=" << des << " design="
                     << DesignElement::type.ToString(space->design[des].design) << " idx=" << other->GetIndex() << " ed=" << de->GetType();
    }

    virtual_elem_map.Push_back(Identifier(de, neighbours));
  }
}

void Function::Local::SetupMultDesignsVirtualElementMap(const Function* f)//, const Phase ph)
{
  // only this element!
  //element_dimension_ = 1; // two boundary "stones" per dimension
  int  dim     = domain->GetDim();
  bool prev    = locality_ == MULT_DESIGNS_PREV_NEXT_AND_REVERSE || locality_ == MULT_DESIGNS_PREV_NEXT;
  bool next    = true; // always
  bool two_signs = locality_ == MULT_DESIGNS_NEXT_AND_REVERSE || locality_ == MULT_DESIGNS_PREV_NEXT_AND_REVERSE;

  int sign_1 =  two_signs ? 1 : Identifier::NO_SIGN;
  int sign_2 =  -1;

  element_dimension_ = 1.0* (two_signs ? 2 : 1);

  UInt elems = space->GetNumberOfElements();

  virtual_elem_map.Reserve(element_dimension_ * elems);

  // the neighbors are the design elements for the same FE-element but with other designs
  // one is not a neighbor of oneself
  StdVector<BaseDesignElement*> neighbours;

  StdVector<unsigned int> des_idx; // the design indices we consider here

  // all designs but the first one
  des_idx.Reserve(space->design.GetSize()-1);
  for(unsigned int i = 1; i < space->design.GetSize(); i++)
    des_idx.Push_back(i);

  // LOG_DBG(func) << "F:L:SMDEM des_idx=" << des_idx.ToString() << " total=" << space->design.ToString();

  for(unsigned int e = 0; e < elems; e++)
  {
    DesignElement* de = func_->elements[e];
    assert((int) e == space->Find(de->elem, true)); // assert that we still are on the right finite element

    //Check that the elements has full neighbours
    //if(de->GetType() == ( (func_->design_ == DesignElement::DEFAULT) ? DesignElement::DENSITY : func_->design_)){
    VicinityElement* ve = de->vicinity;

    // do we have a full neighborhood? All or none as in the original slope paper
    bool full = true;
    if(prev)
    {
      if(ve->design[VicinityElement::X_N] == NULL)
        full = false;
      if(ve->design[VicinityElement::Y_N] == NULL)
        full = false;
      if(dim == 3 && ve->design[VicinityElement::Z_N] == NULL)
        full = false;
    }
    if(next)
    {
      if(ve->design[VicinityElement::X_P] == NULL) full = false;
      if(ve->design[VicinityElement::Y_P] == NULL) full = false;
      if(dim == 3 && ve->design[VicinityElement::Z_P] == NULL) full = false;
    }

    //if this is the case, we can add the identifiers
    if(full)
    {
      neighbours.Resize(0);

      //We begin by adding the derivatives of the first design: same order as vicinity X_N =0, X_P =1; Y_N = 2; etc...
      //for the first design
      //assert((space->FindDesign(DesignElement::G11) == 0));

      for(int a = 0; a < dim; a++)
      {
        DesignElement* prev_de = prev ? ve->GetNeighbour(VicinityElement::ToNeighbour(a, -1)) : NULL;
        DesignElement* next_de = ve->GetNeighbour(VicinityElement::ToNeighbour(a, 1));

        if (prev)
          neighbours.Push_back(prev_de);

        neighbours.Push_back(next_de);
      }

      //We then add the values of the other designs
      for(unsigned int d = 0; d < des_idx.GetSize(); d++)
      {
        // the first design = 0 is no neighbor!
        unsigned des = des_idx[d];
        DesignElement* other = &(space->data[elems * des + e]);
        neighbours.Push_back(other);

        VicinityElement* veother = other->vicinity;

        for(int a = 0; a < dim; a++)
        {
          DesignElement* prev_other = prev ? veother->GetNeighbour(VicinityElement::ToNeighbour(a, -1)) : NULL;
          DesignElement* next_other = veother->GetNeighbour(VicinityElement::ToNeighbour(a, 1));

          if (prev)
            neighbours.Push_back(prev_other);

          neighbours.Push_back(next_other);
        }

        LOG_DBG3(func) << "F:L:SMDEM e=" << e << " el=" << de->elem->elemNum << " d = " << d << " des=" << des << " design="
                       << DesignElement::type.ToString(space->design[des].design) << " idx=" << other->GetIndex() << " ed=" << de->GetType();
      }

      virtual_elem_map.Push_back(Identifier(de, neighbours, sign_1));
      if(two_signs)
        virtual_elem_map.Push_back(Identifier(de, neighbours, sign_2));
    }
  }
}

void Function::Local::SetupShapeElementMap(const Function* func, ShapeDesign* design) {
  element_dimension_ = 0;
  StdVector<ShapeDesign::ShapeConstraint>& shapeconstraints = design->GetShapeConstraints();
  int n = shapeconstraints.GetSize();
  virtual_elem_map.Reserve(n);
  StdVector<BaseDesignElement*> neighbours;
  DesignElement* element(NULL);
  neighbours.Reserve(5);

  for(int e = 0; e < n; e++) {
    ShapeDesign::ShapeConstraint& c = shapeconstraints[e];
    neighbours.Resize(0);
    element = func_->elements[c.param[0]];      
    int t = c.param[1]; // TODO: pattern size may not vary in one constraint!!!, so constraints on a single parameters have added a second parameter with factor 0.0
    if(t >= 0){
      neighbours.Push_back(func_->elements[t]);
    }else{
      if(c.param[0] == 0){ // we must not have element = neighbours[0], this causes errors in sparsity
        neighbours.Push_back(func_->elements[func_->elements.GetSize()-1]);
      }else{
        neighbours.Push_back(func_->elements[0]);
      }
    }
    virtual_elem_map.Push_back(Identifier(element, neighbours));
  }
}

void Function::Local::ToInfo(PtrParamNode in) {
  Function::Type ft = func_->type_;
  in->Get("locality")->SetValue(locality.ToString(locality_));
  in->Get("local_size")->SetValue(virtual_elem_map.GetSize());

  if(IsGlobalized()) {
    in->Get("normalize")->SetValue(normalize_);
    in->Get("power")->SetValue(power_);
  }
  if(RequiresBeta(ft)) {
    in->Get("beta")->SetValue(beta_);
    in->Get("phase")->SetValue(phase.ToString(phase_));
  }

  if(RequiresEps(ft))
    in->Get("eps")->SetValue(eps_);

  // we simply handle periodic only in ShapeMapDesign, extend if you generalize!
  // the perdiodic attribute is for slope and curvature, for the function perdiodic it makes no sense
  if((func_->GetDesignType() == BaseDesignElement::NODE || func_->GetDesignType() == BaseDesignElement::PROFILE) && ft != Function::PERIODIC)
    in->Get("periodic")->SetValue(periodic);

  if(structure_ != NULL)
    structure_->ToInfo(in->Get("neighborhood"));
}

bool Function::Local::IsReverse(Locality loc)
{
  switch(loc)
  {
  case NEXT_AND_REVERSE:
  case PREV_NEXT_AND_REVERSE:
  case DEG_45_STAR_AND_REVERSE:
  case MULT_DESIGNS_NEXT_AND_REVERSE:
  case MULT_DESIGNS_PREV_NEXT_AND_REVERSE:
  case FUNCTION_SPECIFIC_TWO_SIGNS:
    return true;
  default:
    return false;
  }
  return false; // for the stupid compilers
}

bool Function::Local::RequiresEps(Type ft)
{
  switch(ft)
  {
  case MOLE:
  case GLOBAL_MOLE:
  case OVERHANG_HOR:
    return true;
  default:
    return false;
  }
}


bool Function::Local::RequiresBeta(Type ft)
{
  switch(ft)
  {
  case OSCILLATION:
  case GLOBAL_OSCILLATION:
    return true;
  default:
    return false;
  }
}


Function::Local::NeighborhoodStructure::NeighborhoodStructure(Local* local,
    PtrParamNode pn) {
  unsigned int dim = domain->GetDim();
  // sample design element -> assume regular grid
  DesignElement& de = local->space->data[0];

  value = pn->Get("neighbor_value")->As<double>();
  fs = Filter::filterSpace.Parse(
      pn->Get("neighbor_type")->As<string>());
  radius = DesignStructure::FindFilterRadius(fs, &de, value);

  // find the orthogonal dimensions based on radius
  StdVector<double> edges;
  domain->GetGrid()->GetElemShapeMap(de.elem, false)->GetEdgeLength(edges);

  assert(edges.GetSize() == dim);
  orthogonal.Resize(dim);
  for (unsigned int i = 0; i < edges.GetSize(); i++) {
    orthogonal[i] = (int) ((radius / edges[0]) + 0.5); // proper rounding
    LOG_DBG(func)<< "L:NS:NS orthogonal[" << i << "]: radius=" << radius << " edge=" << edges[i] << " -> " << orthogonal[i];
  }

  // validate
  for (unsigned int i = 0; i < orthogonal.GetSize(); i++)
    if (orthogonal[i] == 0)
      throw Exception("your local neighbor radius for '" + Function::type.ToString(local->func_->GetType()) + "' is too small");

  // diagonals are xy, xz and yz plane, only xy in 2D
  for (int i = 0; i < (dim == 2 ? 1 : 3); i++) {
    int x = i != 2 ? 0 : 1; // xy, xz, yz
    int y = i == 0 ? 1 : 2; // xy, xz, yz
    double diag = std::sqrt(edges[x] * edges[x] + edges[y] * edges[y]);
    diagonal.Push_back((int) ((radius / diag) + 0.5));
    LOG_DBG(func) << "L:NS:NS diagonal[0]: x=" << x << " y=" << y << " radius="
                  << radius << " diag=" << diag << " -> " << diagonal[0];
  }
  // all 4 "total" diagonals have the same size
  if (dim == 3) {
    double diag = std::sqrt(edges[0] * edges[0] + edges[1] * edges[1] + edges[2] * edges[2]);
    diagonal.Push_back((int) ((radius / diag) + 0.5));
  }

  assert((dim == 2 && diagonal.GetSize() == 1) || (dim == 3 && diagonal.GetSize() == 4));
}

void Function::Local::NeighborhoodStructure::ToInfo(PtrParamNode in) {
  in->Get("type")->SetValue(Filter::filterSpace.ToString(fs));
  in->Get("value")->SetValue(value);
  in->Get("radius")->SetValue(radius);

  StdVector<double> tmp;
  tmp.Push_back(orthogonal[0] * 2 + 1);
  tmp.Push_back(orthogonal[1] * 2 + 1);
  if (orthogonal.GetSize() == 3)
    tmp.Push_back(orthogonal[2] * 2 + 1);

  in->Get("orthogonal")->SetValue(tmp.ToString());

  tmp.Resize(0);
  for (unsigned int i = 0; i < diagonal.GetSize(); i++)
    tmp.Push_back(diagonal[i] * 2 + 1);

  in->Get("diagonal")->SetValue(tmp.ToString());
}

Function::Local::Identifier::Identifier(BaseDesignElement* elem, BaseDesignElement* prev, BaseDesignElement* next, int si)
{
  this->element = elem;

  assert(next != NULL);
  bool has_prev = prev != NULL;

  this->neighbor.Resize(has_prev ? 2 : 1);

  if (has_prev)
    this->neighbor[0] = prev;

  this->neighbor[has_prev ? 1 : 0] = next;

  this->sign = si;
}

Function::Local::Identifier::Identifier(BaseDesignElement* elem, StdVector<BaseDesignElement*> buddies, int si)
{
  this->element = elem;
  this->neighbor = buddies; // copy constructor
  assert(si == NO_SIGN || si == -1 || si == 1);
  this->sign = si;
}

Function::Local::Identifier::Identifier(BaseDesignElement* elem, StdVector<BaseDesignElement*> buddies, StdVector<BaseDesignElement*> sb_buddies, int si)
{
  this->element = elem;
  this->neighbor = buddies;
  this->sb_neighbor = sb_buddies;
  assert(sb_buddies.GetSize() == 2 * domain->GetDim());
  assert(elem->GetType() == DesignElement::Type::CP && si >= -12 && si <= 12);
  this->sign = si;
}

const BaseDesignElement* Function::Local::Identifier::GetElementByType(DesignElement::Type type) const
{
  for(int i = -1 ; i < (int) neighbor.GetSize(); i++)
  {
    const BaseDesignElement* de = GetElement(i);
    if (de->GetType() == type)
      return de; // do not check for non-uniqueness
  }
//  assert(false); (design of Type type may be param in ParamMat -> see GetDesign)
  return NULL;
}

double Function::Local::Identifier::GetDesign(BaseDesignElement::Type type, const Local* local, const DesignElement::Access access, bool get_parameter) const
{
  const BaseDesignElement* de = GetElementByType(type);
  if(de != NULL)
    return de->GetDesign(access);
  if(get_parameter)
    return Optimization::context->dm->GetParameter(type);
  throw Exception("Design type not found! If it is a ParamMat parameter make sure to query for parameters.");
}


double Function::Local::Identifier::EvalFunction(const Local* local,  bool grad_glob)
{
  // function value
  double fv = 0.0;
  Function* f = local->func_;
  DesignElement::Access access = f->isFiltered() ? DesignElement::SMART : DesignElement::PLAIN;

  // short cut for the gradient in the 1-norm
  if (grad_glob && local->power_ == 1.0) {
    double factor = 1.0;
    if (local->DoNormalizeGlobal())
      factor = (f->type_ == GLOBAL_TWO_SCALE_VOL && !local->space->IsRegular()) ? (1./local->total_vol_) : (1.0 / local->virtual_elem_map.GetSize());

    LOG_DBG(func)<< "L:I:EF: global! p=" << local->power_ << " gg=" << grad_glob << " normalize=" << local->DoNormalizeGlobal() << " -> " << factor;

    if (f->type_ == GLOBAL_BUCKLING_LOAD_FACTOR)
        factor *= (f->GetValue() < f->GetParameter()) ? -1.0 : 0.0;

    return factor;
  }

  switch (f->type_) {

  case GLOBAL_STRESS:
  case LOCAL_STRESS:
  case LOCAL_BUCKLING_LOAD_FACTOR:
  case GLOBAL_BUCKLING_LOAD_FACTOR:
    fv = f->GetValue(); // from local_values, set in ErsatzMaterial::Calc[Global/Local]VonMisesStress[OrLoadFactor]() or ErsatzMaterial::CalcGlobalFunction() for the global case, we might be objective function.
    assert(fv != -1);
    break;

  case SLOPE:
  case GLOBAL_SLOPE:
    fv = CalcSlope();
    break;

  case PERIMETER:
    fv = CalcPerimeter(local->func_->parameter_, 1. / local->virtual_elem_map.GetSize());
    break;

  case PERIODIC:
    fv = CalcPeriodic();
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

  case CURVATURE:
    // for symmetric shape mapping the curvature at the symmetry point reduces to slope.
    // Typically the curvature bound is much tighter than the slope bound,
    // therefore ShapeMapDesign::SetupVirtualShapeElementMap() cheats a single slope constraints
    if(this->neighbor.GetSize() == 1)
      fv = CalcSlope();
    else
      fv = CalcCurvature();
    break;

  case OVERHANG_VERT:
  case OVERHANG_HOR:
    fv = CalcOverhang(f->type_, local->eps_); // not GetEps() as we don't need it for VERT
    break;

  case DISTANCE:
    fv = CalcDistance(-1, false);
    break;

  case BENDING:
    fv = CalcBending(-1, false);
    break;

  case LOCAL_PYTHON_FUNCTION:
    fv = CalcLocalPythonFunc(local);
    break;

  case CONES:
    fv = CalcCones(local);
    break;

  case SUM_MODULI:
  case GLOBAL_SUM_MODULI:
    fv = CalcSumModuli(local, access);
    break;

  case TWO_SCALE_VOL:
  case GLOBAL_TWO_SCALE_VOL:
    fv = CalcTwoScaleVolume(local, access);
    break;

  case PARAM_PS_POS_DEF:
    fv = CalcParamPSPosDef(local, access, -1, false);
    break;

  case POS_DEF_DET_MINOR_1:
  case POS_DEF_DET_MINOR_2:
  case POS_DEF_DET_MINOR_3:
    fv = CalcPosDefDeterminant(-1, local, false, f->type_);
    break;

  case BENSON_VANDERBEI_1:
  case BENSON_VANDERBEI_2:
  case BENSON_VANDERBEI_3:
    fv = CalcBensonVanderbei(-1, local, false, f->type_);
    break;

  case ORTHOTROPIC_TENSOR_TRACE:
    fv = CalcOrthotropicTensorTrace(local, access);
    break;
  case GLOBAL_ORTHOTROPIC_TENSOR_TRACE:
    fv = CalcOrthotropicTensorTrace(local, access);
    break;

  case TENSOR_TRACE:
  case GLOBAL_TENSOR_TRACE:
    fv = CalcTensorTrace(-1, local, false);
    break;

  case TENSOR_NORM:
    fv = CalcTensorNorm(-1, local, false);
    break;

  case DESIGN:
  case GLOBAL_DESIGN:
    fv = CalcDesignBound(f, local, false);
    break;

  case MULTIMATERIAL_SUM:
    fv = CalcMultiMaterialSum(-1, local, false);
    break;
    
  case SHAPE_INF:
    fv = CalcShape(f, local);
    break;

  default:
    assert(false);
    break;
  }

  LOG_DBG2(func) << "L:I:EF: f=" << f->type.ToString(f->type_) << " loc=" << f->IsLocal()
                 << " de=" << element->GetIndex() << " sign=" << sign << " fv=" << fv;

  // handle globalization
  switch (f->type_) {
  case GLOBAL_SLOPE:
  case GLOBAL_OSCILLATION:
  case GLOBAL_MOLE:
  case GLOBAL_JUMP:
  case GLOBAL_DESIGN:
  case GLOBAL_STRESS:
  case GLOBAL_SUM_MODULI:
  case GLOBAL_TWO_SCALE_VOL:
  case GLOBAL_ORTHOTROPIC_TENSOR_TRACE:
  case GLOBAL_TENSOR_TRACE:
  case GLOBAL_BUCKLING_LOAD_FACTOR:
  {
    // we normalize all values by the number of "constraints". Note that it is
    // sufficient for the function value, the gradient is then also right
    double factor;
    if (local->DoNormalizeGlobal()) {
      factor = (f->type_ == GLOBAL_TWO_SCALE_VOL && !local->space->IsRegular()) ? (1./local->total_vol_) : (1.0 / local->virtual_elem_map.GetSize());
    } else {
      factor = 1.;
    }

    double fac = (f->type_ == GLOBAL_BUCKLING_LOAD_FACTOR) ? -1. : 1.;
    double v = std::max(0.0, fac * (fv - f->GetParameter()));

    double p = local->GetPower();

    double res = grad_glob ? fac * p * std::pow(v, p - 1.0) : std::pow(v, p);

    res *= factor;

    LOG_DBG(func)<< "L:I:EF: global! bound=" << f->GetParameter() << " fv=" << std::setprecision(12) << fv << " v=" << v << " p=" << p
       << " factor=" << factor << " gg=" << grad_glob << " power=" << local->GetPower() << " v^p="<< std::pow(v, local->GetPower()) << " -> " << res;

    return res;
  }
  default:
    return fv; // check is done before
  }
}

void Function::Local::Identifier::EvalGradient(const Local* local) {
  // TODO the dynamic_cast might be too slow, check! and do faster by IsObjective()
  // we need this pointers, note that C++ makes NULL for an invalid dynamic cast
  Function* funct = local->func_;
  Function::Type ft = funct->type_;
  Condition* g = dynamic_cast<Condition*>(funct);
  Objective* f = dynamic_cast<Objective*>(funct);
  DesignElement::Access access = funct->isFiltered() ? DesignElement::SMART : DesignElement::PLAIN;
  assert((f == NULL && g != NULL) || (f != NULL && g == NULL));

  LOG_DBG2(func) << "L:I:EvalGrad: f=" << funct->type.ToString(funct->type_) << " de="
                 << ( typeid(element) == typeid(DesignElement*) ? (int)dynamic_cast<DesignElement*>(element)->elem->elemNum : -1 ) << " sign=" << sign;

  // are we global? then we don't do anything if the globalization function gives zero
  // this applies the gradient of the globalization function (max(0, fv)^2)
  // EvalFunction() is very fast for power=1!
  double grad_glob_fv = local->IsGlobalized() ? EvalFunction(local, true) : -1.0; // if not global we don't need grad_glob_fv

  if(local->IsGlobalized() && grad_glob_fv == 0.0) {
    LOG_DBG2(func) << "L:I:EvalGrad: f=" << funct->type.ToString(funct->type_) << " de=" << element->GetIndex() << " sign=" << sign << " fv=0.0 -> return immediately";
    return;
  }
  assert(local->IsGlobalized() || ft == PERIMETER || g != NULL); // only constraints are local

  /** in the local python grad case we call python for a "full" local gradient */
  if(ft == LOCAL_PYTHON_FUNCTION) {
    CalcLocalPythonGrad(local->func_->py_local_grad_.Mine(), local);
  }

  for (int n = -1, nn = neighbor.GetSize(); n < nn; n++)
  {
    double gv = -5.0;

    switch (ft)
    {
    case SLOPE:
    case GLOBAL_SLOPE:
      gv = CalcSlopeGradient(n);
      break;

    case PERIMETER:
      gv = CalcPerimeterGradient(n, local->func_->parameter_, 1. / local->virtual_elem_map.GetSize());
      break;

    case PERIODIC:
      gv = CalcPeriodicGradient(n);
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

    case CURVATURE:
      // see CURVATUE in EvalFunction()
      if(this->neighbor.GetSize() == 1)
        gv = CalcSlopeGradient(n);
      else
        gv = CalcCurvatureGradient(n);
      break;

    case OVERHANG_VERT:
    case OVERHANG_HOR:
      gv = CalcOverhangGradient(n, g->type_, local->eps_); // no GetEps()!
      break;

    case DISTANCE:
      gv = CalcDistance(n, true);
      break;

    case BENDING:
      gv = CalcBending(n, true);
      break;

    case LOCAL_PYTHON_FUNCTION:
    {
      const Vector<double>& vec = local->func_->py_local_grad_.Mine();
      gv = vec[n + 1]; // n = -1 is this element, 0 is first neighbor
      break;
    }

    case CONES:
      gv = CalcConesGradient(n, local);
      break;

    case GLOBAL_STRESS:
    case LOCAL_STRESS:
    case LOCAL_BUCKLING_LOAD_FACTOR:
    case GLOBAL_BUCKLING_LOAD_FACTOR:
      assert(false); // in ErsatzMaterial::Calc[Local/Global]VonMisesStress[OrLoadFactor]() only!
      break;

    case SUM_MODULI:
    case GLOBAL_SUM_MODULI:
      gv = CalcSumModuli(local, access, n, true);
      break;

    case TWO_SCALE_VOL:
    case GLOBAL_TWO_SCALE_VOL:
      gv = CalcTwoScaleVolume(local, access, n, true);
      break;

    case PARAM_PS_POS_DEF:
      gv = CalcParamPSPosDef(local, access, n, true);
      break;

    case POS_DEF_DET_MINOR_1:
    case POS_DEF_DET_MINOR_2:
    case POS_DEF_DET_MINOR_3:
      gv = CalcPosDefDeterminant(n, local, true, ft);
      break;

    case BENSON_VANDERBEI_1:
    case BENSON_VANDERBEI_2:
    case BENSON_VANDERBEI_3:
      gv = CalcBensonVanderbei(n, local, true, ft);
      break;

    case ORTHOTROPIC_TENSOR_TRACE:
      gv = CalcOrthotropicTensorTrace(local, access, n, true);
      break;
    case GLOBAL_ORTHOTROPIC_TENSOR_TRACE:
      gv = CalcOrthotropicTensorTrace(local, access, n, true);
      break;

    case TENSOR_TRACE:
    case GLOBAL_TENSOR_TRACE:
      gv = CalcTensorTrace(n, local, true);
      break;

    case TENSOR_NORM:
      gv = CalcTensorNorm(n, local, true);
      break;

    case DESIGN:
    case GLOBAL_DESIGN:
      gv = CalcDesignBound(funct, local, true);
      break;

    case MULTIMATERIAL_SUM:
      gv = CalcMultiMaterialSum(n, local, true);
      break;
      
    case SHAPE_INF:
      gv = CalcShapeGradient(funct, local, n);
      break;

    default:
      assert(false);
      break;
    }

    LOG_DBG2(func) << "L:I:EvalGrad: f=" << funct->type.ToString(funct->type_) << " de="
                   << element->GetIndex() << " sign=" << sign << " n=" << n << " des=" << DesignElement::type.ToString(GetElement(n)->GetType())
                   << " curr=" << GetElement(n)->GetIndex() << " gv=" << gv;

    // post process the globalized functions. The perimeter is not globalized in that sense
    if (local->IsGlobalized())
    {
      gv *= grad_glob_fv;
      LOG_DBG2(func) << "L:I:EvalGrad: f=" << funct->type.ToString(funct->type_) << " de="
                     << element->GetIndex() << " sign=" << sign << " n=" << n
                     << " curr=" << GetElement(n)->GetIndex()
                     << " bound! grad_glob_gv=" << grad_glob_fv << " new gv=" << gv;
    }

    BaseDesignElement* bde = GetElement(n);
    assert(bde != NULL);


    // the perimeter is not globalized by sum max(g-g*, 0)^p but it is not local!
    if(!local->IsGlobalized() && ft != PERIMETER)
    {
      // reset the constraint data. Note, as we are local, there are no side effects by elements
      bde->Reset(DesignElement::CONSTRAINT_GRADIENT, g);
      if(g->isFiltered())
      {
        DesignElement* de = dynamic_cast<DesignElement*>(bde);
        if(de != NULL && !de->simp->filter.IsEmpty())
        {
          unsigned int fix = de->simp->DetermineFilterIndexNonInlined();
          const StdVector<Filter::NeighbourElement> neighborhood = de->simp->filter[fix].neighborhood;
          // for constraints using filtered design variables also reset the constraint data in the filter neighborhood
          for(unsigned int j = 0, nj = neighborhood.GetSize(); j < nj; j++)
          {
            DesignElement* de2 =  neighborhood[j].neighbour;
            de2->Reset(DesignElement::CONSTRAINT_GRADIENT, g);
            for(unsigned int k = 0, nk = de2->simp->filter[fix].neighborhood.GetSize(); k < nk; k++)
              de2->simp->filter[fix].neighborhood[k].neighbour->constraintGradient[g->GetIndex()] = 0.0;  // This is much faster than calling Reset()
          }
        }
      }
    }


    bde->AddGradient(f, g, gv);
    LOG_DBG2(func) << "L:I:EvalGrad: f=" << funct->type.ToString(funct->type_) << " de="
                   << element->GetIndex() << " sign=" << sign << " n=" << n
                   << " curr=" << GetElement(n)->GetIndex() << " gv=" << gv
                   << " stored_gv=" << bde->GetPlainGradient(funct)
                   << " current_position: " << (g != NULL ? ((LocalCondition*) g)->GetCurrentPosition()+1 : -1); //somehow only seems to work for constraints
  } // end loop over n
}

double Function::Local::Identifier::CalcSlope() const {
  double mine = element->GetDesign(DesignElement::SMART);
  assert(this->neighbor.GetSize() == 1);
  double other = neighbor[0]->GetDesign(DesignElement::SMART);

  assert(this->sign == 1 || this->sign == -1 || this->sign == BOTH); // 1/-1 for reverse (not snopt) and BOTH (-1000) for snopt
  double s = this->sign == -1 ? -1.0 : 1.0;

  LOG_DBG3(func) << "L:I:CS de=" << element->GetIndex() << " other=" << (typeid(neighbor[0]) == typeid(DesignElement*) ? (int)dynamic_cast<DesignElement*>(neighbor[0])->elem->elemNum : -1 )
                 << " sign=" << sign << " slope -> " << (s * (mine - other));
  return s * (mine - other);
}


double Function::Local::Identifier::CalcSlopeGradient(int neigh_idx) const
{
  assert(neigh_idx == -1 || neigh_idx == 0);
  // we have the cases sign=1, sign=-1, NO_SIGN. NO_SIGN is handled as sign=-1
  LOG_DBG3(func) << "L:I:CSG de=" << element->GetIndex() << " ni=" << neigh_idx << " sign=";

  if (neigh_idx == -1)
    return sign == -1 ? -1.0 : 1.0;
  else
    return sign == -1 ? 1.0 : -1.0;
}

double Function::Local::Identifier::CalcOverhang(Function::Type ft, double eps) const
{
  // this(node)->elem is implicit, then this(profile), then prev(node) and prev(profile) if exist, then next(node) and next(profile)
  //assert(this->sign == 1); // the other stuff is not considered yet
  assert(neighbor.GetSize() == 3);
  assert(element->GetType() == DesignElement::NODE);
  assert(neighbor[0]->GetType() == DesignElement::PROFILE);
  assert(neighbor[1]->GetType() == DesignElement::NODE);
  assert(neighbor[2]->GetType() == DesignElement::PROFILE);

  double a = element->GetPlainDesignValue(); // this
  double w = neighbor[0]->GetPlainDesignValue();
  double an = neighbor[1]->GetPlainDesignValue(); // next
  double wn = neighbor[2]->GetPlainDesignValue();

  assert(ft == OVERHANG_VERT || ft == OVERHANG_HOR);
  assert(!(ft == OVERHANG_HOR && eps <= 0.0)); // hor needs eps
  assert(!(ft == OVERHANG_HOR && sign != NO_SIGN)); // hor knows only one sign as it is real smooth abs() as >= cannot be double bounded

  // 1/-1 for reverse (not snopt) and BOTH (-1000) for snopt. The BOTH case is the positive case. See BaseOptimizer::GetBounds() and Function::IsDoubleBounded()
  assert(this->sign == 1 || this->sign == -1 || this->sign == BOTH);

  double s = sign == -1 ? -1.0 : 1.0; // 1 for 1 and BOTH

  double res = -1.0;

  if(ft == OVERHANG_HOR)
  {
    double tmp = ((an-wn) - (a-w)); // abs for the lower points only. For the upper part of a horizontal structure we can build everything. Note, we need >= c^*
    res = SmoothAbs(tmp, eps);
    LOG_DBG3(func) << "L:I:CO ft=" << Function::type.ToString(ft) << " a=" << element->GetIndex() << "(" << a << ") an=" << neighbor[1]->GetIndex() << "(" << an << ") w=" << w << " wn=" << wn
                   << " tmp=" << tmp << " eps=" << eps << " -> " << res;
  }
  else
  {
    res = s == 1 ? ((an+wn) - (a+w)) : ((a-w) - (an-wn)); // the right part is checked for right overhang only, the left part for a left overhang only
    LOG_DBG3(func) << "L:I:CO ft=" << Function::type.ToString(ft) + " sign=" << s
                   << " a=" << element->GetIndex() << "(" << a << ") an=" << neighbor[1]->GetIndex() << "(" << an << ") w=" << w << " wn=" << wn << " -> " << res;
  }
  return res;
}

double Function::Local::Identifier::CalcOverhangGradient(int neigh_idx, Function::Type ft, double eps) const
{
  assert(neighbor.GetSize() == 3);
  assert(element->GetType() == DesignElement::NODE);        // a
  assert(neighbor[0]->GetType() == DesignElement::PROFILE); // w
  assert(neighbor[1]->GetType() == DesignElement::NODE);    // an
  assert(neighbor[2]->GetType() == DesignElement::PROFILE); // wn

  assert(!(ft == OVERHANG_HOR && eps <= 0.0)); // hor needs eps
  assert(!(ft == OVERHANG_HOR && sign != NO_SIGN)); // hor knows only one sign as it is real smooth abs() as >= cannot be double bounded

  //assert(this->sign == 1); // the other stuff is not considered yet. n = i+1
  // OVERHANG_HOR:  |(an-wn) - (a-w)| >= c^* // real smoothed abs required!! no double signs
  // OVERHANG_VERT: ((an+wn) - (a+w)) <= c^* and ((a-w) - (an-wn)) <= c^*
  assert(sign == 1 || sign == -1 || sign == BOTH || sign == NO_SIGN);
  double s = sign == -1 ? -1.0 : 1.0; // 1 for 1 and BOTH
  // var        a    w   an  wn
  // hor_s=1   -1    1   +1   -1
  // vert_s=1  -1   -1   +1   +1
  // vert_s=-1 +1   -1   -1   +1

  // for the horizontal case, df is the factor for DerivSmoothAbs()
  double df = 0.0;
  switch(neigh_idx)
  {
  case -1: // a
    // same for hor and vert
    df = -1 * s;
    break;
  case 0:  // w
    df = (ft == OVERHANG_HOR) ? 1.0: -1.0;
    break;
  case 1: // an
     df = s;
     break;
  case 2: // wn
    df = (ft == OVERHANG_HOR) ? -1.0 : 1.0;
    break;
  default:
    break;
  }
  assert(df != 0.0);

  if(ft == OVERHANG_HOR)
    // we have to flip sign because we are lower bound
    return -1 * df * DerivSmoothAbs(CalcOverhang(ft, eps), eps);
  else
    return df;
}

double Function::Local::Identifier::CalcCones(const Local* local) const {
  assert(local->GetLocality() == MULT_DESIGNS_PREV_NEXT);
  unsigned int dim = domain->GetDim();

  assert(sb_neighbor.GetSize() == 2 * dim);
  assert(element->GetType() == DesignElement::CP);
  for(unsigned int i = 0; i < neighbor.GetSize(); ++i)
    assert(neighbor[i]->GetType() == DesignElement::CP);
  for(unsigned int i = 0; i < sb_neighbor.GetSize(); ++i)
    assert(sb_neighbor[i]->GetType() == DesignElement::CP);

  // this->x is implicit, then this->y, then this->z, then either prev->x ... or next->x ...
  double x, y, z = 0, xo, yo, zo = 0;
  x = sb_neighbor[0]->GetPlainDesignValue(); // this
  y = sb_neighbor[1]->GetPlainDesignValue();
  if(dim == 2) {
    xo = sb_neighbor[2]->GetPlainDesignValue(); // other
    yo = sb_neighbor[3]->GetPlainDesignValue();
  } else {
    z = sb_neighbor[2]->GetPlainDesignValue();
    xo = sb_neighbor[3]->GetPlainDesignValue(); // other
    yo = sb_neighbor[4]->GetPlainDesignValue();
    zo = sb_neighbor[5]->GetPlainDesignValue();
  }

  // add initial control point position to design elements
  SplineBoxDesign* sbd = dynamic_cast<SplineBoxDesign*>(local->space);
  // intentional integer division
  unsigned int cp_idx_e = sb_neighbor[0]->GetIndex() / dim;
  unsigned int cp_idx_n = sb_neighbor[3]->GetIndex() / dim;
  Point cp_e = sbd->GetInitialControlPoint(cp_idx_e);
  Point cp_n = sbd->GetInitialControlPoint(cp_idx_n);
  x += cp_e[0];
  y += cp_e[1];
  z += cp_e[2];
  xo += cp_n[0];
  yo += cp_n[1];
  zo += cp_n[2];


  assert(sign >= -12 && sign <= 12);
  assert((dim == 2 && sign >= -4 && sign <= 4) || dim == 3);

  double angle = local->func_->parameter_;
  assert(angle >= 0);
  double tana = std::tan(angle/180. * M_PI);

  double res = -1.0;

  // in the following comments we have comment(x) = code(xo) - code(x),
  // comment(x_0) = code(x), comment(x_1) = code(xo) and comment(a|b_*) = code(tana)
  // cases 1 and 2 describe cone in x direction,
  // cases 3 and 4 describe cone in y direction
  if(dim == 2) {
    switch(std::abs(sign)) {
    case 1:
      // a_1 * x >= y -> - a_1 * x + y <= 0 -> - a_1 * x_1 + a_1 * x_0 + y_1 - y_0 <= 0
      res =   tana * x - y - tana * xo + yo;
      break;
    case 2:
      // a_2 * x >= -y -> - a_2 * x - y <= 0 -> - a_2 * x_1 + a_2 * x_0 - y_1 + y_0 <= 0
      res =   tana * x + y - tana * xo - yo;
      break;
    case 3:
      // b_1 * y >= x -> x - b_1 * y <= 0 -> x_1 - x_0 - b_1 * y_1 + b_1 * y_0 <= 0
      res = - x + tana * y + xo - tana * yo;
      break;
    case 4:
      // b_2 * y >= - x -> - x - b_2 * y <= 0 -> - x_1 + x_0 - b_2 * y_1 + b_2 * y_0 <= 0
      res =   x + tana * y - xo - tana * yo;
      break;
    }
    LOG_DBG3(func) << "F:L:I:CC e=" << element->GetIndex() << " (" << x << "," << y << ") en="
        << sb_neighbor[3]->GetIndex() << " (" << x << "," << y << ") res=" << res;
  } else {
    // Todo implement angle
    // cases 1,2,3 and 4 describe cone in x direction,
    // cases 5,6,7 and 8 describe cone in y direction,
    // cases 9,10,11 and 12 describe cone in y direction
    switch(std::abs(sign)) {
    case 1:
      res = + x - y - z - xo + yo + zo;
      break;
    case 2:
    case 10:
      res = + x - y + z - xo + yo - zo;
      break;
    case 3:
    case 7:
      res = + x + y - z - xo - yo + zo;
      break;
    case 4:
    case 8:
    case 12:
      res = + x + y + z - xo - yo - zo;
      break;
    case 5:
    case 9:
      res = - x + y + z + xo - yo - zo;
      break;
    case 6:
      res = - x + y - z + xo - yo + zo;
      break;
    case 11:
      res = - x - y + z + xo + yo - zo;
      break;
    }
    LOG_DBG3(func) << "F:L:I:CC e=" << element->GetIndex() << " (" << x << "," << y << "," << z << ") en="
        << neighbor[3]->GetIndex() << " (" << x << "," << y << "," << z << ") res=" << res;
  }

  res *= sign < 0 ? -1.0 : 1.0;

  // Multiply with -1. Then the regularization parameter is a positive lower bound.
  res *= -1;

  return res;
}

double Function::Local::Identifier::CalcConesGradient(int neigh_idx, const Local* local) const {
  assert(local->GetLocality() == MULT_DESIGNS_PREV_NEXT);
  unsigned int dim = domain->GetDim();

  assert(sb_neighbor.GetSize() == 2 * dim);
  assert(element->GetType() == DesignElement::CP);
  for(unsigned int i = 0; i < neighbor.GetSize(); ++i)
    assert(neighbor[i]->GetType() == DesignElement::CP);
  for(unsigned int i = 0; i < sb_neighbor.GetSize(); ++i)
    assert(sb_neighbor[i]->GetType() == DesignElement::CP);

  assert(sign >= -12 && sign <= 12);
  assert((dim == 2 && sign >= -4 && sign <= 4) || dim == 3);

  double angle = local->func_->parameter_;
  assert(angle >= 0);
  double tana = std::tan(angle/180. * M_PI);

  const BaseDesignElement* bde = GetElement(neigh_idx);
  int sb_neigh_idx = sb_neighbor.Find(const_cast<BaseDesignElement*>(bde));
  LOG_DBG(func) << "F:L:I:CCG: e=" << bde->GetIndex() << " sb_neigh_idx=" << sb_neigh_idx;

  double df = 0.0;

  // this->x is implicit, then this->y, then this->z, then either prev->x ... or next->x ...
  if(dim == 2) {
    switch(std::abs(sign)) {
    case 1:
      // res =   tana * x - y - tana * xn + yn;
      df  = (sb_neigh_idx == 0 || sb_neigh_idx == 3) ? 1 : -1;
      df *= (sb_neigh_idx == 0 || sb_neigh_idx == 2) ? tana : 1;
      break;
    case 2:
      // res =   tana * x + y - tana * xn - yn;
      df  = (sb_neigh_idx == 0 || sb_neigh_idx == 1) ? 1 : -1;
      df *= (sb_neigh_idx == 0 || sb_neigh_idx == 2) ? tana : 1;
      break;
    case 3:
      // res = - x + tana * y + xn - tana * yn;
      df  = (sb_neigh_idx == 1 || sb_neigh_idx == 2) ? 1 : -1;
      df *= (sb_neigh_idx == 1 || sb_neigh_idx == 3) ? tana : 1;
      break;
    case 4:
      // res =   x + tana * y - xn - tana * yn;
      df  = (sb_neigh_idx == 0 || sb_neigh_idx == 1) ? 1 : -1;
      df *= (sb_neigh_idx == 1 || sb_neigh_idx == 3) ? tana : 1;
      break;
    }
  } else {
    switch(std::abs(sign)) {
    // Todo @see CalcCones
    case 1:
      // res = + x - y - z - xn + yn + zn;
      df = (sb_neigh_idx == 0 || sb_neigh_idx == 4 || sb_neigh_idx == 5) ? 1 : -1;
      break;
    case 2:
    case 10:
      // res = + x - y + z - xn + yn - zn;
      df = (sb_neigh_idx == 0 || sb_neigh_idx == 2 || sb_neigh_idx == 4) ? 1 : -1;
      break;
    case 3:
    case 7:
      // res = + x + y - z - xn - yn + zn;
      df = (sb_neigh_idx == 0 || sb_neigh_idx == 1 || sb_neigh_idx == 5) ? 1 : -1;
      break;
    case 4:
    case 8:
    case 12:
      // res = + x + y + z - xn - yn - zn;
      df = (sb_neigh_idx == 0 || sb_neigh_idx == 1 || sb_neigh_idx == 2) ? 1 : -1;
      break;
    case 5:
    case 9:
      // res = - x + y + z + xn - yn - zn;
      df = (sb_neigh_idx == 1 || sb_neigh_idx == 2 || sb_neigh_idx == 3) ? 1 : -1;
      break;
    case 6:
      // res = - x + y - z + xn - yn + zn;
      df = (sb_neigh_idx == 1 || sb_neigh_idx == 3 || sb_neigh_idx == 5) ? 1 : -1;
      break;
    case 11:
      // res = - x - y + z + xn + yn - zn;
      df = (sb_neigh_idx == 2 || sb_neigh_idx == 3 || sb_neigh_idx == 4) ? 1 : -1;
      break;
    }
  }

  df *= sign < 0 ? -1.0 : 1.0;

  // @see CalcCones
  df *= -1;

  return df;
}

double Function::Local::Identifier::CalcPerimeter(double eps, double l_k) const
{
  // P = sum_k^K l_k ( sqrt( (<p>_k)**2 + eps**2 ) - eps )
  // where K is the number of interfaces and l_k is the length of the interface
  // <p>_k is the jump at the interface rho_i - rho_i+1
  // We use the NEXT locality with right and upper neighbor. We don't want to count interfaces double
  // eps is param as this constraint does not look like a local one in the xml file
  // We ignore the interfaces to the boundaries
  assert(this->neighbor.GetSize() == 1);

  double mine  = element->GetDesign(DesignElement::SMART);
  double other = neighbor[0]->GetDesign(DesignElement::SMART);
  double res   = l_k * (sqrt( (mine-other) * (mine-other) + eps * eps ) - eps);

  LOG_DBG3(func) << "L:I:CP de=" << element->GetIndex() << " other=" << neighbor[0]->GetIndex()
                 << " mine=" << mine << " other=" << other << " eps=" << eps << " l_k=" << l_k << " -> " << res;
  return res;
}

double Function::Local::Identifier::CalcPerimeterGradient(int neigh_idx, double eps, double l_k) const
{
  // P = sum_k p_k
  // p_k' = l_k * 0.5 * ((rho_i - rho_i+1)**2 + eps**2)**-0.5 * 2 * (rho_i - rho_i+1) * s
  // s = 1 for d p_k / d rho_i and -1 for d p_k / d rho_i+1

  double mine  = element->GetDesign(DesignElement::SMART);
  double other = neighbor[0]->GetDesign(DesignElement::SMART);
  double s = neigh_idx == -1 ? 1.0 : -1.0;

  // using not the std::pow() gives wrong results!
  double res = l_k * std::pow((mine-other)*(mine-other) + eps*eps, -0.5) * (mine-other) * s;

  LOG_DBG3(func) << "L:I:CPG de=" << element->GetIndex() << " other=" << neighbor[0]->GetIndex()
                 << " mine=" << mine << " other=" << other << " eps=" << eps << " l_k=" << l_k << " neigh_idx=" << neigh_idx << " -> " << res;
  return res;
}

double Function::Local::Identifier::CalcPeriodic() const
{
  // this - neighbor
  assert(this->neighbor.GetSize() == 1);

  double mine  = element->GetDesign(DesignElement::SMART);
  double other = neighbor[0]->GetDesign(DesignElement::SMART);
  double res   = mine - other;

  LOG_DBG3(func) << "L:I:CP de=" << element->GetIndex() << " other=" << neighbor[0]->GetIndex()
                 << " mine=" << mine << " other=" << other << " -> " << res;
  return res;
}

double Function::Local::Identifier::CalcPeriodicGradient(int neigh_idx) const
{
  assert(neigh_idx == -1 || neigh_idx == 0);
  double s = neigh_idx == -1 ? 1.0 : -1.0;
  return s;
}

double Function::Local::Identifier::CalcOscillation(double beta) const {
  assert(sign == 1 || sign == -1);

  double own = element->GetDesign(DesignElement::SMART);
  // we divide the neighbors in lower and upper
  unsigned int half = neighbor.GetSize() / 2;

  tmp1.Resize(half);
  for (unsigned int i = 0; i < half; i++)
    tmp1[i] = neighbor[i]->GetDesign(DesignElement::SMART);
  double prev = sign == -1 ? SmoothMax(tmp1, beta) : SmoothMin(tmp1, beta);

  tmp2.Resize(half);
  for (unsigned int i = 0; i < half; i++)
    tmp2[i] = neighbor[i + half]->GetDesign(DesignElement::SMART);
  double next = sign == -1 ? SmoothMax(tmp2, beta) : SmoothMin(tmp2, beta);

  double min_max = 0.0; // min or max
  double res = 0.0;

  if (sign == 1) {
    // "Heaviside"(rho_i - max( rho_i-1, rho_i+1)
    // double smaller = std::max(0.0, own - std::max(left, right));
    min_max = beta < 0 ? std::max(prev, next) : SmoothMax(prev, next, beta);
    res = own - min_max;
  } else {
    // "Heaviside"(min( rho_i-1, rho_i+1) - rho_i)
    // double larger = std::max(0.0, std::min(left, right) - own);
    min_max = beta < 0 ? std::min(prev, next) : SmoothMin(prev, next, beta);
    res = min_max - own;
  }

  // LOG_DBG3(func) << "L:I:CO de=" << element->ToString() << " neigh=" << BaseDesignElement::ToString(neighbor)
  //               << " vals=" << tmp1.ToString() << "; " << tmp2.ToString() << " sign=" << sign << " own=" << own
  //               << " prev=" << prev << " next=" << next << " smooth=" << min_max << " hard="
  //               << (sign == 1 ? (own - std::max(prev, next)) : (std::min(prev, next) - own)) << " -> " << res;
  return res;
}

double Function::Local::Identifier::CalcOscillationGradient(int neigh_idx, double beta)
{
  assert(beta >= 0);

  // the own value is not within the min/max stuff
  if (neigh_idx == -1)
    return sign == 1 ? 1 : -1;

  // we divide the neighbors in lower and upper
  unsigned int half = neighbor.GetSize() / 2;

  tmp1.Resize(half);
  for (unsigned int i = 0; i < half; i++)
    tmp1[i] = neighbor[i]->GetDesign(DesignElement::SMART);
  double prev = sign == -1 ? SmoothMax(tmp1, beta) : SmoothMin(tmp1, beta);

  tmp2.Resize(half);
  for (unsigned int i = 0; i < half; i++)
    tmp2[i] = neighbor[i + half]->GetDesign(DesignElement::SMART);
  double next = sign == -1 ? SmoothMax(tmp2, beta) : SmoothMin(tmp2, beta);

  // sign = -1: min( max(x_0 .. x_half-1), max(x_half .. x_max) ) - own
  // sign =  1: own - max( min(x_0 .. x_half-1), min(x_half .. x_max) )

  // example derivative for sign = 1: - max'(min(x_0 .. x_half-1), min(x_half .. x_max) * min'(max(x_0 .. x_half-1) * 1

  // are we within prev (-1) or next (1)
  int side = neigh_idx < (int) half ? -1 : 1;

  // outer is simple
  double outer =  sign == 1 ? DerivSmoothMax(prev, next, beta, side) :  DerivSmoothMin(prev, next, beta, side);
  // we do not handle neighbor.GetSize() == 2 special as the Smooth tools are fast for this special case
  // inner depends on neigh_idx within the first or the next
  double inner = 0.0;
  if (side == -1)
    inner = sign == -1 ? DerivSmoothMax(tmp1, beta, neigh_idx) : DerivSmoothMin(tmp1, beta, neigh_idx);
  else
    inner = sign == -1 ? DerivSmoothMax(tmp2, beta, neigh_idx - half) : DerivSmoothMin(tmp2, beta, neigh_idx - half);
  assert(neighbor.GetSize() > 2 || close(inner, 1.0));

  return (sign == 1 ? -1.0 : 1.0) * outer * inner;
}

double Function::Local::Identifier::CalcMole(double eps) const {
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
  for (unsigned int i = 0; i < neighbor.GetSize() / 2; i++) {
    tmp1.Push_back(GetElement(i)->GetDesign(DesignElement::SMART));
    LOG_DBG3(func)<< "L:I:CalcMole de=" << element->ToString() << " other=" << GetElement(i)->ToString() << " v=" << tmp1.Last();
  }

  tmp1.Push_back(element->GetDesign(DesignElement::SMART));
  LOG_DBG3(func)<< "L:I:CalcMole de=" << element->ToString() << " other=" << element->ToString() << " v=" << tmp1.Last();

  for (unsigned int i = neighbor.GetSize() / 2; i < neighbor.GetSize(); i++) {
    tmp1.Push_back(GetElement(i)->GetDesign(DesignElement::SMART));
    LOG_DBG3(func)<< "L:I:CalcMole de=" << element->ToString() << " other=" << GetElement(i)->ToString() << " v=" << tmp1.Last();
  }

  assert(tmp1.GetSize() == neighbor.GetSize() + 1);

  double sum = 0.0;
  for (unsigned int i = 0; i < tmp1.GetSize() - 1; i++) {
    sum += SmoothAbs(tmp1[i + 1] - tmp1[i], eps);
    LOG_DBG3(func)<< "L:I:CalcMole de=" << element->ToString() << " i=" << i << "|" << tmp1[i+1] << " - " << tmp1[i] << "|="
    << (std::abs(tmp1[i+1] - tmp1[i])) << " smoothed=" << SmoothAbs(tmp1[i+1] - tmp1[i], eps) << " sum=" << sum;
  }
  double result = sum - SmoothAbs(tmp1.Last() - tmp1[0], eps);
  LOG_DBG3(func)<< "L:I:CalcMole de=" << element->ToString() << " bound=|" << tmp1.Last() << " - " << tmp1[0] << "| smoothed="
  << SmoothAbs(tmp1.Last() - tmp1[0], eps) << " -> " << result;
  return result;
}

double Function::Local::Identifier::CalcMoleGradient(int neigh_idx, double eps) {
  // see comments in the forward function implementation
  // three cases for the sensitivity analysis: first, intermediate, last

  // sort all values in a sequence -> TODO! optimize and do not copy and paste
  tmp1.Resize(0); // is static and keeps capacity
  for (unsigned int i = 0; i < neighbor.GetSize() / 2; i++)
    tmp1.Push_back(GetElement(i)->GetDesign(DesignElement::SMART));

  tmp1.Push_back(element->GetDesign(DesignElement::SMART));

  for (unsigned int i = neighbor.GetSize() / 2; i < neighbor.GetSize(); i++)
    tmp1.Push_back(GetElement(i)->GetDesign(DesignElement::SMART));

  assert(tmp1.GetSize() == neighbor.GetSize() + 1);

  // correct the index
  unsigned int half = neighbor.GetSize() / 2;
  int idx = neigh_idx == -1 ? half :
            neigh_idx < (int) half ? neigh_idx : neigh_idx + 1;
  assert(idx >= 0 && idx < (int ) tmp1.GetSize());

  double res = 0.0;

  if (idx == 0) {
    res = DerivSmoothAbs(tmp1.Last() - tmp1[0], eps)
        - DerivSmoothAbs(tmp1[1] - tmp1[0], eps);
    LOG_DBG3(func)<< "L:I:CalcMoleGrad de=" << element->ToString() << " neigh_idx=" << neigh_idx << " idx=" << idx << " tmp=" << tmp1.ToString()
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

double Function::Local::Identifier::CalcJump() const {
  assert(sign == NO_SIGN);
  assert(neighbor.GetSize() == 2);

  // sin(pi*(x_i-1 - x_+1))^2
  // no own value!
  double prev = neighbor[0]->GetDesign(DesignElement::SMART);
  double next = neighbor[1]->GetDesign(DesignElement::SMART);

  double sin = std::sin(M_PI * (prev - next));

  // LOG_DBG3(func) << "L:I:CJ de=" << element->ToString() << " prev=" << neighbor[0]->ToString() << "/" << prev
  //               << " next=" << neighbor[1]->ToString() << "/" << next << " slope=" << (prev-next) << " -> sin*sin";

  return sin * sin;
}

double Function::Local::Identifier::CalcJumpGradient(int neigh_idx) const {
  // g(x)=sin(pi*(x_i-1 - x_+1))^2
  // d g(x)/d x_i-1 = 2 * sin(pi*(x_i-1 - x_+1)) * cos (pi*(x_i-1 - x_+1)) * PI

  // no own value!
  if (neigh_idx == -1)
    return 0.0;

  double prev = neighbor[0]->GetDesign(DesignElement::SMART);
  double next = neighbor[1]->GetDesign(DesignElement::SMART);

  double slope = prev - next;

  assert(neigh_idx == 0 || neigh_idx == 1);
  double factor = neigh_idx == 0 ? 1.0 : -1.0;

  return 2.0 * std::sin(M_PI * slope) * std::cos(M_PI * slope) * M_PI * factor;
}

double Function::Local::Identifier::CalcBump() const {
  assert(sign == NO_SIGN);
  assert(neighbor.GetSize() == 2);

  // (x_i-1 - x_i)(x_i - x_i+1)
  // no own value!
  double prev = neighbor[0]->GetDesign(DesignElement::SMART);
  double mine = element->GetDesign(DesignElement::SMART);
  double next = neighbor[1]->GetDesign(DesignElement::SMART);

  double val = (prev - mine) * (mine - next);

  LOG_DBG3(func)<< "L:I:CB de=" << element->ToString()
  << " prev=" << neighbor[0]->ToString() << "/" << prev
  << " next=" << neighbor[1]->ToString() << "/" << next << " mine=" << mine
  << " -> " << (prev - mine) << " * " << (mine - next) << " = " << val;

  return val;
}

double Function::Local::Identifier::CalcBumpGradient(int neigh_idx) const {
  double prev = neighbor[0]->GetDesign(DesignElement::SMART);
  double mine = element->GetDesign(DesignElement::SMART);
  double next = neighbor[1]->GetDesign(DesignElement::SMART);

  double res = 0.0;

  switch (neigh_idx) {
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

double Function::Local::Identifier::CalcCurvature() const
{
  assert(neighbor.GetSize() == 2);
  double s = this->sign == -1 ? -1.0 : 1.0;

  // (x_i-1 -2*x_i + x_i+1) or the reverse
  double prev = neighbor[0]->GetDesign(DesignElement::SMART);
  double mine = element->GetDesign(DesignElement::SMART);
  double next = neighbor[1]->GetDesign(DesignElement::SMART);

  double val = s * (prev - 2 * mine + next);

  LOG_DBG3(func) << "L:I:CC de=" << element->ToString()
                 << " prev=" << neighbor[0]->ToString() << "/" << prev << " mine=" << mine
                 << " next=" << neighbor[1]->ToString() << "/" << next << " s=" << s
                 << " -> " << val;
  return val;
}

double Function::Local::Identifier::CalcCurvatureGradient(int neigh_idx) const
{
  assert(neigh_idx <= 1);
  double s = this->sign == -1 ? -1.0 : 1.0;

  return s * (neigh_idx == -1 ? -2 : 1);
}


double Function::Local::Identifier::CalcDistance(int neigh_idx, bool grad) const
{
  assert(neighbor.GetSize() >= 3); // to be also used by CalcBending()
  double px = GetElement(-1)->GetPlainDesignValue();
  double py = GetElement(0)->GetPlainDesignValue();
  double qx = GetElement(1)->GetPlainDesignValue();
  double qy = GetElement(2)->GetPlainDesignValue();

  double dist = std::sqrt((qx-px)*(qx-px) + (qy-py)*(qy-py));

  LOG_DBG(func) << "F:L:I:CD grad=" << grad << " ni=" << neigh_idx << " px=" << px << " py=" << py << " qx=" << qx << " qy=" << qy;

  if(grad)
    switch(neigh_idx)
    {
    case -1: // px
      return -(qx-px)/dist;
    case 0: // qy
      return -(qy-py)/dist;
    case 1: // qx
      return (qx-px)/dist;
    case 2: // qy
      return (qy-py)/dist;
    default:
      return 0; // when used by bending, this is the right answer
    }
  else
    return dist;
}

double Function::Local::Identifier::CalcBending(int neigh_idx, bool grad) const
{
  // see SpaghettiDesign::SetupVirtualShapeElementMap()
  double prev = 0;
  double mine = 0;
  double next = 0;

  int base = 3; // the first index for normals, distance is before. More complex for 3D and fixed for nodes
  switch(bending)
  {
  case ZNZ:
    mine = GetElement(base+0)->GetPlainDesignValue();
    break;
  case ZNN:
    mine = GetElement(base+0)->GetPlainDesignValue();
    next = GetElement(base+1)->GetPlainDesignValue();
    break;
  case NNN:
    prev = GetElement(base+0)->GetPlainDesignValue();
    mine = GetElement(base+1)->GetPlainDesignValue();
    next = GetElement(base+2)->GetPlainDesignValue();
    break;
  case NNZ:
    prev = GetElement(base+0)->GetPlainDesignValue();
    mine = GetElement(base+1)->GetPlainDesignValue();
    break;
  case NO_BENDING:
    assert(false);
    break;
  }

  // see CalcCurvature(Gradient)
  double s = this->sign == -1 ? -1.0 : 1.0;

  double curv = s * (prev - 2 * mine + next);
  double dist = CalcDistance(-1, false);
  assert(dist > 0);

  double res = -42;
  if(grad) {
    double d_curv = 0; // true for neigh_idx within distance nodes
    if(neigh_idx >= base)
    {
      switch(bending)
      {
      case ZNZ:
        assert(neigh_idx == base);
        d_curv = s * -2;
        break;
      case ZNN:
        assert(neigh_idx <= base+1);
        d_curv = s * (neigh_idx == base ? -2 : 1);
        break;
      case NNZ:
        assert(neigh_idx <= base+1);
        d_curv = s * (neigh_idx == base+1 ? -2 : 1);
        break;
      case NNN:
        d_curv = s * (neigh_idx == base+1 ? -2 : 1);
        break;
      case NO_BENDING:
        assert(false);
        break;
      }
    }

    double d_dist = CalcDistance(neigh_idx, true);

    res = (d_curv * dist - curv * d_dist) / (dist * dist);
    LOG_DBG2(func) << "F:L:CD(" << neigh_idx << ") gr=" << grad << " curv=" << curv << " d_curv=" << d_curv
                   << " dist=" << dist << " d_dist=" << d_dist << " -> " << res;
  }
  else
    res = curv / dist;

  return res;
}


double Function::Local::Identifier::CalcSumModuli(const Local* local, DesignElement::Access access, int neigh_idx, bool derivative) const
{
  double E1 = GetDesign(DesignElement::EMODULISO, local, access, true);
  double E3 = GetDesign(DesignElement::EMODUL, local, access, true);
  double G = GetDesign(DesignElement::GMODUL, local, access, true);
  double theta = GetDesign(DesignElement::POISSON, local, access, true);

  int dim = domain->GetDim();
  if(dim ==2){ //case PLANE_STRESS, reformulated theta version
    if(derivative)
    {
      switch(GetElement(neigh_idx)->GetType())
      {
      case DesignElement::EMODULISO:
      case DesignElement::EMODUL:
        return 1.0/(1.0-theta);
      case DesignElement::GMODUL:
        return 2.0;
      case DesignElement::POISSON:
      {
        return (E1+E3)/((1.0-theta)*(1.0-theta));
      }
      default:
        return 0.0;
      }
    }
    return (E1+E3)/(1.0-theta)+2*G;
  }
  else { // 3D case original version without theta, theta = nu_{oi}
    double nuiso = GetDesign(DesignElement::POISSONISO, local, access, true);
    double nuoisqrd = theta*theta;
    if(derivative)
    {
      switch(GetElement(neigh_idx)->GetType())
      {
      case DesignElement::EMODULISO:
      {
        double n = (2*E1*nuoisqrd - E3 + E3*nuiso);
        return (8*nuoisqrd*nuoisqrd)/(4*nuoisqrd*nuoisqrd*(1+nuiso)) - (E3*E3*(2*nuoisqrd + 1)*(nuiso - 1))/(n*n);
      }
      case DesignElement::EMODUL:
      {
        double n = (2*E1*nuoisqrd - E3 + E3*nuiso);
        return 1 - (2*E1*E1*nuoisqrd*(2*nuoisqrd + 1))/(n*n);
      }
      case DesignElement::GMODUL:
        return 4.0;
      case DesignElement::POISSON:
      {
        double n = (2*E1*nuoisqrd - E3 + E3*nuiso);
        return (4*E3*E1*theta*(E3 + E1 - E3*nuiso))/(n*n);
      }
      case DesignElement::POISSONISO:
      {
        double n = (nuiso + (2*E1*nuoisqrd)/E3 - 1);
        return (E1*(2*nuoisqrd + 1))/(n*n) - (2*E1)/((nuiso + 1)*(nuiso + 1));
      }

      default:
        return 0.0;
      }
    }
    return E3 + 4*G - (2*E1*nuoisqrd + E1)/(nuiso + (2*E1*nuoisqrd)/E3 - 1) + (2*E1)/(nuiso + 1);
  }
}


double Function::Local::Identifier::CalcOrthotropicTensorTrace(const Local* local, DesignElement::Access access, int neigh_idx, bool derivative) const
{
  double e11 = GetDesign(DesignElement::MECH_11, local, access, true);
  double e22 = GetDesign(DesignElement::MECH_22, local, access, true);
  double e33 = GetDesign(DesignElement::MECH_33, local, access, true);
  double e12 = GetDesign(DesignElement::MECH_12, local, access, true);
  double lowerEigBound = GetDesign(DesignElement::LOWER_EIG_BOUND, local, access, true);

  if(derivative)
  {
    DesignElement::Type type = GetElement(neigh_idx)->GetType();
    if (type == DesignElement::MECH_11)
      return 2.0*e11;
    else if(type == DesignElement::MECH_22)
      return 2.0*e22;
    else if(type == DesignElement::MECH_33)
      return 1.0;
    else if(type == DesignElement::MECH_12)
      return 4.0*e12;
    else if(type == DesignElement::LOWER_EIG_BOUND)
    {
      std::cout << "Warning: lowerEigBound is supposed to be a parameter only, not a design!" << std::endl;
      return 3.0;
    }
    else return 0.0;
  }
  else
    return 3.0*lowerEigBound+e11*e11+2.0*e12*e12+e22*e22+e33;
}

double Function::Local::Identifier::CalcLatticeVolume3D(const Local* local, DesignElement::Access access, int neigh_idx, bool derivative) const {
  // temporary data structure
  Vector<double> p(3);
  if (Optimization::context->dm->GetType() == DesignMaterial::HOM_ISO_C1) {
    p[0] = GetDesign(DesignElement::STIFF1, local, access, true);
    p[1] = p[0];
    p[2] = p[0];
  } else {
    p[0] = GetDesign(DesignElement::STIFF1, local, access, true);
    p[1] = GetDesign(DesignElement::STIFF2, local, access, true);
    p[2] = GetDesign(DesignElement::STIFF3, local, access, true);
  }

  // this is the case, where the volume is our parameter
  if (p[1] == 0 && p[2] == 0) {
    if (!derivative)
      return p[0];
    else
      return 1.0;
  }

  // this is the case, where the volume is interpolated with the values in the catalogue
  assert(local->volumeInterpolator_ != NULL);
  if (!derivative) {
    return local->volumeInterpolator_->EvaluateFunc(p);
  } else {
    switch (GetElement(neigh_idx)->GetType()) {
    case DesignElement::STIFF1:
      return local->volumeInterpolator_->EvaluateDeriv(p, 0);
    case DesignElement::STIFF2:
      return local->volumeInterpolator_->EvaluateDeriv(p, 1);
    case DesignElement::STIFF3:
      return local->volumeInterpolator_->EvaluateDeriv(p, 2);
    default:
      return 0.0;
    }
  }
  //should never be reached
  return -1.0;
}

double Function::Local::Identifier::CalcTwoScaleVolume(const Local* local, DesignElement::Access access, int neigh_idx, bool derivative) const {
  DesignElement* de = dynamic_cast<DesignElement*>(element);
  int dim = domain->GetDim();

  assert(Optimization::context->dm);
  if (Optimization::context->dm->GetType() == DesignMaterial::HOM_ISO_C1 && dim == 2) {
    throw Exception("CalcTwoScaleVolume is not implemented for dim = 2 and HOM_ISO_C1.");
  }

  double vol;
  bool regular = local->space->IsRegular();
  /** if grid is nonregular, the volume has to be scaled by element size */
  if (!regular) {
    assert(local->total_vol_ != 0);
  }
  /** svol is a scaling factor for unstructured, nonregular grids. */
  double svol = regular ? 1.0 : de->CalcVolume();
  LOG_DBG2(func) << "Element volume =  " << de->CalcVolume();

  if (Optimization::context->dm->GetType() == DesignMaterial::HOM_ISO_C1 && dim == 3) {
    return svol * CalcLatticeVolume3D(local, access, neigh_idx, derivative);
  }

  double stiff1 = GetDesign(DesignElement::STIFF1, local, access, true);
  double stiff2 = GetDesign(DesignElement::STIFF2, local, access, true);

  if (Optimization::context->dm->GetInterpolationMethod() == DesignMaterial::SG) {
    Vector<double> p(3);
    p[0] = stiff1;
    p[1] = stiff2;
    p[2] = GetDesign(DesignElement::SHEAR1, local, access, true);
    return svol * Optimization::context->dm->CalcHomVolume(p, GetElement(neigh_idx)->GetType(), derivative);
  }

  if (!derivative) {
    if (dim == 2) {
      // this is the volume formula for a cross or frame structure
      // but it also works for DesignMaterial with only one parameter (stiff1), if stiff2 = 0
      return svol * (stiff1 + stiff2 - stiff1 * stiff2);
    } else {
      return svol * CalcLatticeVolume3D(local, access, neigh_idx, derivative);
    }
  } else {
    if (dim == 2) {
      switch (GetElement(neigh_idx)->GetType()) {
      case DesignElement::STIFF1:
        return svol * (1.0 - stiff2);
      case DesignElement::STIFF2:
        return svol * (1.0 - stiff1);
      default:
        return 0.0;
      }
    } else {
      vol = svol * CalcLatticeVolume3D(local, access, neigh_idx, derivative);
      assert(vol != -1);
      return vol;
    }
  }
  //should never be reached
  assert(false);
  return -1.0;
}

double Function::Local::Identifier::CalcParamPSPosDef(const Local* local, DesignElement::Access access, int neigh_idx,
    bool derivative) const {
  double E1 = GetDesign(DesignElement::EMODULISO, local, access, true);
  double E3 = GetDesign(DesignElement::EMODUL, local, access, true);
  double nu31 = GetDesign(DesignElement::POISSON, local, access, true);

  if (derivative)
    switch (GetElement(neigh_idx)->GetType()) {
    case DesignElement::EMODULISO:
      return -nu31 * nu31;
    case DesignElement::EMODUL:
      return 1.0;
    case DesignElement::POISSON:
      return -2.0 * nu31 * E1;
    default:
      assert(false);
      return 0.0;
    }
  else {
    LOG_DBG3(func) << "Local::Local e_num=" << element->GetIndex() << ", E3-E1*nu31^2=" << E3-E1*nu31*nu31;
    return E3-E1*nu31*nu31;
  }
}


/* Condition: det(G-vId) >= eps*/

double Function::Local::Identifier::CalcPosDefDeterminant(int neigh_idx, const Local* local, bool derivative, Type type) const {
  const Condition* g = dynamic_cast<const Condition*>(local->func_);

  double v = g->GetParameter();
  double eps = 1.0 * g->GetBoundValue();

  MaterialTensor<double> tens(HILL_MANDEL);
  bool ok = Optimization::context->dm->GetTensor(tens, g->GetDesignType(), PLANE_STRAIN, dynamic_cast<DesignElement*>(element)->elem, DesignElement::NO_DERIVATIVE, HILL_MANDEL);
  Matrix<double>& E = tens.GetMatrix(HILL_MANDEL);
  // the sub-tensor-type does'nt matter
  // we need the HILL_MANDEL representation which is the plain design while it is transformed to Voigt for simulation (elasticity only)
  assert(ok);

  // LOG_DBG3(func) << "L::I::CPDD e_num=" << element->elem->elemNum << " v=" << v << " obs=" << obs << " eps=" <<  eps << " E=" << E.ToString(0, false);

  // see e.g. the phd thesis of Sonja Lehmann chap. 6.4 for implementation
  // for (E - v*I) >= 0 we have the determinants
  // d_1 can also be implemented by box constraints!
  // d_1 := e_11 - v >= 0
  // d_2 := (e_11-v)*(e_22-v) - e_12^2 >= 0
  // d_3 is done by Sonja using the Laplace formula, we use Sarrus
  // d_3 := ...

  double ret = -12345678.0 * (ok ? 1.0 : 1.0); // stupid stuff to use ok in release mode

  double e11, e22, e33, e23, e13, e12;

  switch (type) {
  case POS_DEF_DET_MINOR_1:
    e11 = E[0][0];
    // Standard: e11-v;
    // Standard with eps: e11-v-eps;

    if (!derivative)
      // ret = e11-v-eps;
      ret = e11 - v;
    else
      ret = neigh_idx == -1 ? 1.0 : 0.0;
    if (g->GetDesignType() == DesignElement::DIELEC_ALL)
      ret *= -1.0;
    break;

  case POS_DEF_DET_MINOR_2:
    e11 = E[0][0];
    e22 = E[1][1];
    e12 = E[0][1];
    // Standard: (e11-v)*(e22-v) - (e12*e12);
    // Standard with eps: (e11-v-eps)*(e22-v) - (e12*e12) - eps;
    if (!derivative)
      ret = (e11 - v - eps) * (e22 - v) - (e12 * e12) - eps;
    else {
      switch (GetElement(neigh_idx)->GetType()) {
      case DesignElement::MECH_11:
        ret = e22 - v;
        break;
      case DesignElement::DIELEC_11:
        ret = -1.0 * (e22 - v);
        break;
        // case DesignElement::TENSOR22: ret = e11-v-eps; break;
      case DesignElement::MECH_12:
        ret = -2.0 * e12;
        break;
      case DesignElement::DIELEC_12:
        ret = 2.0 * e12;
        break;
      case DesignElement::MECH_22:
        ret = e11 - v - eps;
        break;
      case DesignElement::DIELEC_22:
        ret = -1.0 * (e11 - v - eps);
        break;
      default:
        assert(false);
      }
    }
    break;

  case POS_DEF_DET_MINOR_3:
    e11 = E[0][0];
    e22 = E[1][1];
    e33 = E[2][2];
    e23 = E[1][2];
    e13 = E[0][2];
    e12 = E[0][1];

    // Sarrus: e11*e22*e33+e12*e23*e13+e13*e12*e23-e13*e22*e13-e12*e21*e33-e11*e23*e23;
    // Sarrus symmetric: (e11-v)*(e22-v)*(e33-v) + 2.0*e12*e23*e13 - e13*(e22-v)*e13 - e12*e12*(e33-v) - (e11-v)*e23*e23;
    // Sonja: (e33-v) *((e11-v)*(e22-v) - (e12*e12)) + 2.0* e12*e13 *e23 - e13**2 * (e22-v) - e23**2 * (e11-v)
    // Sonja with eps: ((e33-v) *((e11-v-eps)*(e22-v) - (e12*e12) - eps) + 2.0* e12*e13 *e23 - e13**2 * (e22-v) - e23**2 * (e11-v-eps)) - eps;
    if (!derivative) {
      ret = (e33 - v) * ((e11 - v - eps) * (e22 - v) - (e12 * e12) - eps)
          + 2.0 * e12 * e13 * e23 - e13 * e13 * (e22 - v)
          - e23 * e23 * (e11 - v - eps) - eps;
    } else {
      switch (GetElement(neigh_idx)->GetType()) {
      /** this are the nice gradients w.r.t Sarrus
       case DesignElement::TENSOR11: ret = (e22-v)*(e33-v) - e23*e23; break;
       case DesignElement::TENSOR22: ret = (e11-v)*(e33-v) - e13*e13; break;
       case DesignElement::TENSOR33: ret = (e11-v)*(e22-v) - e12*e12; break;
       case DesignElement::TENSOR23: ret = 2.0*e12*e13     - 2.0*e23*(e11-v); break;
       case DesignElement::TENSOR13: ret = 2.0*e12*e23     - 2.0*e13*(e22-v); break;
       case DesignElement::TENSOR12: ret = 2.0*e23*e13     - 2.0*e12*(e33-v); break;
       */
      case DesignElement::MECH_11:
        ret = (e22 - v) * (e33 - v) - e23 * e23;
        break;
      case DesignElement::MECH_12:
        ret = 2.0 * e23 * e13 - 2.0 * e12 * (e33 - v);
        break;
      case DesignElement::MECH_22:
        // Sarrus: ret = (e11-v)*(e33-v) - e13*e13;
        ret = (e33 - v) * (e11 - v - eps) - e13 * e13;
        // ret = (e11-v)*(e33-v) - e13*e13;
        break;
      case DesignElement::MECH_23:
        // Sarrus: ret = 2.0*e12*e13     - 2.0*e23*(e11-v);
        ret = 2.0 * e12 * e13 - 2.0 * e23 * (-v - eps + e11);
        // ret = 2.0*e12*e13     - 2.0*e23*(e11-v);
        break;
      case DesignElement::MECH_13:
        ret = 2.0 * e12 * e23 - 2.0 * e13 * (e22 - v);
        break;
      case DesignElement::MECH_33:
        // Sarrus: ret = (e11-v)*(e22-v) - e12*e12;
        ret = (e22 - v) * (-v - eps + e11) - eps - e12 * e12;
        // ret = (e11-v)*(e22-v) - e12*e12;
        break;
      default:
        assert(false);
      }
    }
    break;

  default:
    assert(false);
    break;
  } // end switch f->GetType()
  assert(ret != 12345678.0);
  LOG_DBG3(func) << "L::I::CPDD e_num=" << dynamic_cast<DesignElement*>(element)->elem->elemNum << " g=" << Function::type.ToString(g->GetType()) << " for " << Function::type.ToString(type)
                 << " ni=" << neigh_idx << " v=" << v << "  des=" << DesignElement::type.ToString(GetElement(neigh_idx)->GetType()) << " d=" << derivative << " -> " << ret;
  return ret;
}

double Function::Local::Identifier::CalcBensonVanderbei(int neigh_idx,
    const Local* local, bool derivative, Type type) const {
  const Condition* g = dynamic_cast<const Condition*>(local->func_);
  Matrix<double> E;
  // (E - v*I) >= gamma
  double v = g->GetParameter();
  // in case, we are used as "approximation" of the benson vanderbei constraints:
  //   det1(x) = bv1(x) - eps  <-- this eps is also on the left side!!
  double eps = g->GetBoundValue();

  // the sub-tensor-type does'nt matter
  // we need the HILL_MANDEL representation which is the plain design while it is transformed to Voigt for simulation
  assert(false);
  // local->space->GetErsatzMaterialTensor(E, PLANE_STRAIN, dynamic_cast<DesignElement*>(element)->elem, DesignElement::NO_DERIVATIVE, HILL_MANDEL);

  LOG_DBG3(func) << "L::I::CBV e_num=" << element->GetIndex() << " v=" << v << " E=" << E.ToString();

  double ret = -12345678.0;

  // soja: e1   e2   e3   e4   e5   e6
  double e11, e12, e22, e13, e23, e33;

  switch (type) {
  case BENSON_VANDERBEI_1:
    e11 = E[0][0];
    if (!derivative)
      // ret = e11-v;
      ret = e11 - v - eps;
    else
      ret = neigh_idx == -1 ? 1.0 : 0.0;
    break;

  case BENSON_VANDERBEI_2:
    e11 = E[0][0];
    e22 = E[1][1];
    e12 = E[0][1];
    if (!derivative)
      // (e3-v) - e2*e2/(e1-v)
      //ret = (e22-v) - (e12*e12)/(e11-v);
      ret = (e22 - v) - (e12 * e12) / (e11 - v - eps);
    else {
      switch (GetElement(neigh_idx)->GetType()) {
      case DesignElement::MECH_11:
        // ret = (e12*e12)/((e11-v)*(e11-v));
        ret = (e12 * e12) / ((e11 - v - eps) * (e11 - v - eps));
        break;
      case DesignElement::MECH_12:
        ret = -2.0 * e12 / (e11 - v - eps);
        break;
      case DesignElement::MECH_22:
        ret = 1.0;
        break;
      default:
        assert(false);
      }
    }
    break;

  case BENSON_VANDERBEI_3:
    e11 = E[0][0];
    e22 = E[1][1];
    e33 = E[2][2];
    e23 = E[1][2];
    e13 = E[0][2];
    e12 = E[0][1];

    if (!derivative) {
      // (e6-v) + (2.0*e2*e4*e5 - e4*e4*(e3-v) - e5*e5*(e1-v))/((e3-v)*(e1-v)-e2*e2);
      // ret = (e33-v) + (2.0*e12*e13*e23 - e13*e13*(e22-v) - e23*e23*(e11-v))/((e22-v)*(e11-v)-e12*e12);
      // from maxima: (d3 -eps) / (d2 - eps)
      ret = (-e23 * e23 * (-v - eps + e11)
          + ((e22 - v) * (-v - eps + e11) - eps - e12 * e12) * (e33 - v)
          - e13 * e13 * (e22 - v) - eps + 2.0 * e12 * e13 * e23)
          / ((e22 - v) * (-v - eps + e11) - eps - e12 * e12);

    } else {
      switch (GetElement(neigh_idx)->GetType()) {
      // case DesignElement::TENSOR11: ret = -((-e13*e13*(e22-v)-e23*e23*(e11-v)+2.0*e12*e13*e23)*(e22-v))
      //                                    / pow((e11-v)*(e22-v)-e12*e12, 2)
      //                                  - (e23*e23)
      //                                    / ((e11-v)*(e22-v)-e12*e12);
      case DesignElement::MECH_11:
        ret = ((e22 - v) * (e33 - v) - e23 * e23)
            / ((e22 - v) * (-v - eps + e11) - eps - e12 * e12)
            - ((-e23 * e23 * (-v - eps + e11)
                + ((e22 - v) * (-v - eps + e11) - eps - e12 * e12) * (e33 - v)
                - e13 * e13 * (e22 - v) - eps + 2.0 * e12 * e13 * e23)
                * (e22 - v))
                / pow(((e22 - v) * (-v - eps + e11) - eps - e12 * e12), 2);

        break;
        // case DesignElement::TENSOR12: ret = (2.0*e13*e23)
        //                                    / ((e11-v)*(e22-v)-e12*e12)
        //                                  + (2.0*e12*(-e13*e13*(e22-v)-e23*e23*(e11-v)+2.0*e12*e13*e23))
        //                                     / pow((e11-v)*(e22-v)-e12*e12, 2);
      case DesignElement::MECH_12:
        ret = (2.0 * e13 * e23 - 2.0 * e12 * (e33 - v))
            / ((e22 - v) * (-v - eps + e11) - eps - e12 * e12)
            + (2.0 * e12
                * (-e23 * e23 * (-v - eps + e11)
                    + ((e22 - v) * (-v - eps + e11) - eps - e12 * e12)
                        * (e33 - v) - e13 * e13 * (e22 - v) - eps
                    + 2.0 * e12 * e13 * e23))
                / pow(((e22 - v) * (-v - eps + e11) - eps - e12 * e12), 2);

        break;
        // case DesignElement::TENSOR22: ret = -((-e13*e13*(e22-v)-e23*e23*(e11-v)+2.0*e12*e13*e23)*(e11-v))
        //                                    / pow((e11-v)*(e22-v)-e12*e12,2)
        //                                  - (e13*e13)
        //                                    / ((e11-v)*(e22-v)-e12*e12);
      case DesignElement::MECH_22:
        ret = ((e33 - v) * (-v - eps + e11) - e13 * e13)
            / ((e22 - v) * (-v - eps + e11) - eps - e12 * e12)
            - ((-e23 * e23 * (-v - eps + e11)
                + ((e22 - v) * (-v - eps + e11) - eps - e12 * e12) * (e33 - v)
                - e13 * e13 * (e22 - v) - eps + 2.0 * e12 * e13 * e23)
                * (-v - eps + e11))
                / pow((e22 - v) * (-v - eps + e11) - eps - e12 * e12, 2);

        break;
        // case DesignElement::TENSOR13: ret = (2.0*e12*e23-2.0*e13*(e22-v))
        //                                    / ((e11-v)*(e22-v)-e12*e12);
      case DesignElement::MECH_13:
        ret = (2.0 * e12 * e23 - 2.0 * e13 * (e22 - v))
            / ((e22 - v) * (-v - eps + e11) - eps - e12 * e12);
        break;
        // case DesignElement::TENSOR23: ret = (2.0*e12*e13-2.0*e23*(e11-v))
        //                                      / ((e11-v)*(e22-v)-e12*e12);
      case DesignElement::MECH_23:
        ret = (2.0 * e12 * e13 - 2.0 * e23 * (-v - eps + e11))
            / ((e22 - v) * (-v - eps + e11) - eps - e12 * e12);
        break;
      case DesignElement::MECH_33:
        ret = 1.0;
        break;
      default:
        assert(false);
      }
    }
    break;

  default:
    assert(false);
    break;
  } // end switch f->GetType()
  assert(ret != 12345678.0);
  LOG_DBG3(func)<< "L::I::CBV e_num=" << element->GetIndex() << " g=" << Function::type.ToString(g->GetType())
  << " ni=" << neigh_idx << "  des=" << DesignElement::type.ToString(GetElement(neigh_idx)->GetType()) << " d=" << derivative << " -> " << ret;
  return ret;
}

double Function::Local::Identifier::CalcMultiMaterialSum(int neigh_idx, const Local* local, bool derivative) const
  {
    Matrix<double> E;

    double ret = 0.0;

    if(!derivative)
    {
      for(int i=-1; i < (int) neighbor.GetSize(); ++i)
      {
        ret += GetElement(i)->GetDesign(DesignElement::PLAIN);
        LOG_DBG3(func) << "L::I::CMMS e_num=" << element->GetIndex() << " i=" << i << " e=" <<  dynamic_cast<const DesignElement*>(GetElement(i))->elem->elemNum << " -> ret";
        // does not work in the mag opt case with density + rhsDensity
        // << " mi=" << (dynamic_cast<const DesignElement*>(GetElement(i))->multimaterial->index
      }
    }
    else
    {
      ret = 1.0;
    }
    LOG_DBG3(func) << "L::I::CMMS e_num=" << element->GetIndex() << " ni=" << neigh_idx << " d=" << derivative << " -> " << ret;
    return ret;
  }

double Function::Local::Identifier::CalcTensorTrace(int neigh_idx, const Local* local, bool derivative) const
{

  Function* f= local->func_;
  MaterialTensorNotation notation = f->notation_;
  SubTensorType stt = f->ctxt->stt;
  Elem* elem = dynamic_cast<DesignElement*>(element)->elem;
  const DesignElement* de = dynamic_cast<const DesignElement*>(GetElement(neigh_idx));
  DesignElement::Type der = derivative ? de->GetType() : DesignElement::NO_DERIVATIVE;

  MaterialTensor<double> tens(notation);
  bool ok = Optimization::context->dm->GetTensor(tens, f->GetDesignType(), stt, elem, der, notation); // the sub-tensor-type DOES matter)
  assert(ok);

  Matrix<double>& E = tens.GetMatrix(notation);
  assert((local->func_->GetDesignType() == DesignElement::DIELEC_TRACE && E.GetNumRows() == 2) || (local->func_->GetDesignType() != DesignElement::DIELEC_TRACE && (E.GetNumRows() == 3 || E.GetNumRows() == 6)));
  LOG_DBG3(func) << "L::I::CTT e_num=" << element->GetIndex() << " dt=" << de->type.ToString(local->func_->GetDesignType()) << " E=" << E.ToString();

  double ret = E.Trace() * (ok ? 1.0 : 1.0); // to use ok in assert
  assert(!(derivative && local->func_->GetDesignType() == DesignElement::DIELEC_TRACE && ret != -1.0));
  assert(!(derivative && notation == HILL_MANDEL && de->GetType() == DesignElement::MECH_33 && ret != 1.0));
  assert(!(derivative && notation == VOIGT && de->GetType() == DesignElement::MECH_33 && ret != 0.5));

  LOG_DBG3(func)<< "L::I::CTT e_num=" << de->elem->elemNum << " ni=" << neigh_idx << " nt=" << de->type.ToString(de->GetType()) << " d=" << derivative << " -> " << ret;
  return ret;
}

double Function::Local::Identifier::CalcTensorNorm(int neigh_idx, const Local* local, bool derivative) const
{
  const BaseDesignElement* de = GetElement(neigh_idx);
  assert(local->func_->GetDesignType() == DesignElement::PIEZO_ALL);
  MaterialTensor<double> tens(NO_NOTATION);
  // as we square we do not need the linear derivative
  Optimization::context->dm->GetPiezoCouplingTensor(tens, dynamic_cast<DesignElement*>(element)->elem, DesignElement::NO_DERIVATIVE);

  Matrix<double>& E = tens.GetMatrix(NO_NOTATION);
  LOG_DBG3(func) << "L::I::CTN e_num=" << element->GetIndex() << " E=" << E.ToString();

  double ret = 0.0;

  if(!derivative)
    ret = pow(E.NormL2(), 2);
  else {
    switch(de->GetType()) {
    case DesignElement::PIEZO_11:
      ret = 2.0 * E[0][0];
      break;
    case DesignElement::PIEZO_12:
      ret = 2.0 * E[0][1];
      break;
    case DesignElement::PIEZO_13:
      ret = 2.0 * E[0][2];
      break;
    case DesignElement::PIEZO_21:
      ret = 2.0 * E[1][0];
      break;
    case DesignElement::PIEZO_22:
      ret = 2.0 * E[1][1];
      break;
    case DesignElement::PIEZO_23:
      ret = 2.0 * E[1][2];
      break;
    default:
      assert(false);
      break;
    }
  }

  LOG_DBG3(func)<< "L::I::CTN e_num=" << element->GetIndex() << " ni=" << neigh_idx << " d=" << derivative << " -> " << ret;
  return ret;
}

double Function::Local::Identifier::CalcDesignBound(Function* f, const Local* l, bool derivative) const {
  assert(this->neighbor.GetSize() == 0);

  if(derivative)
  {
    return 1.0;
  }
  else
  {
    switch(f->GetAccess())
    {
    case Function::DEFAULT:
    case Function::PLAIN:
      return element->GetPlainDesignValue(); // will be replaced below in case

    case Function::FILTERED:
      return element->GetDesign(DesignElement::SMART); // Function::Access != DesignElement::Access!!

    case Function::PHYSICAL:
      // the gradient does not penalize :(

      // const Context& ctxt = Optimization::manager.GetContext(f->GetExcitation());
      // const TransferFunction* tf = l->space->GetTransferFunction(f->GetDesignType(), ctxt.ToApp());
      // DesignElement* de = dynamic_cast<DesignElement*>(element);
      // assert(element != NULL);
      // val = tf->Transform(de, DesignElement::SMART);
      // val = de->GetDesign(DesignElement::SMART);
      // LOG_DBG3(func) << "L::I::CDB e=" << element->GetIndex() << " de=" << de->ToString() << " plain=" << element->GetPlainDesignValue() << " smart=" << de->GetDesign(DesignElement::SMART) << " -> " << val;
      assert(false);
      EXCEPTION("not implemented");

    case Function::NO_ACCESS:
      assert(false);
    }
  }
  assert(false);
  return -1.0; // please compiler
}
  
double Function::Local::Identifier::CalcShape(Function* f, const Local* l) const {
  assert(f->type_ == SHAPE_INF);
  int idx = dynamic_cast<LocalCondition*>(f)->GetCurrentRelativePosition();
  ShapeDesign::ShapeConstraint& c = dynamic_cast<ShapeDesign*>(l->space)->GetShapeConstraints()[idx];
  // note that if neighbor[0] should not be given, it points to the first design element and c.factor[1] is 0.0
  double ret = this->element->GetDesign() * c.factor[0] - this->neighbor[0]->GetDesign() * c.factor[1];
  return(ret);
}

double Function::Local::Identifier::CalcShapeGradient(Function* f, const Local* l, int neigh_idx) const {
  assert(f->type_ == SHAPE_INF);
  int idx = dynamic_cast<LocalCondition*>(f)->GetCurrentRelativePosition();
  ShapeDesign::ShapeConstraint& c = dynamic_cast<ShapeDesign*>(l->space)->GetShapeConstraints()[idx];
  return(neigh_idx == -1 ? c.factor[0] : -c.factor[1]); // if no neighbor given c.factor[0][1] is 0.0
}

string Function::Local::Identifier::ToString() const
{
  std::stringstream ss;
  ss << "s=" << sign << " sv=" << signs.ToString() << " b=" << bending;
  ss << "e=" << element->ToString() << " n=";
  for(const auto e : neighbor)
    ss << e->ToString();

  return ss.str();
}

