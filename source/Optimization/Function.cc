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

DECLARE_LOG(func)
DEFINE_LOG(func, "opt_func")

// instantiation of the static elements is in Optimization::SetEnums()
Enum<Function::Type> Function::type;
Enum<Function::Local::Locality> Function::Local::locality;

// speed up by sharing
StdVector<double> Function::Local::Identifier::tmp;

const int Function::Local::Identifier::NO_SIGN = -1000;

using boost::lexical_cast;

Function::Function(PtrParamNode pn)
{
  this->harmonic_    = BasePDE::IsComplex(domain->GetDriver()->GetAnalysisType());

  this->pn = pn;
  this->local = NULL;

  this->type_ = type.Parse(pn->Get("type")->As<std::string>());

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

  // check parameter
  switch(type_)
  {
  case PENALIZED_VOLUME:
  case GAP:
  case GLOBAL_SLOPE:
  case GLOBAL_CHECKERBOARD:
    if(!pn->Has("parameter"))
      throw Exception("function '" + type.ToString(type_) + "' requires the 'parameter' attribute");
    break;

  default:
    if(pn->Has("parameter"))
      throw Exception("function '" + type.ToString(type_) + "' does not accept 'parameter' attribute");
  }


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
    case GLOBAL_CHECKERBOARD:
    case MOLE:
    case GLOBAL_MOLE:
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
  // we check for valid ocurence of paramter in the constructor
  if(pn->Has("parameter"))
    info->Get("parameter")->SetValue(parameter_);
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
  case GLOBAL_CHECKERBOARD:
  case MOLE:
  case GLOBAL_MOLE:
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
  case GLOBAL_CHECKERBOARD:
  case CHECKERBOARD:
  case GLOBAL_SLOPE:
  case SLOPE:
    assert(space->IsRegular()); // VicinityElements work only on a regular grid
    // the design elements require the vicinity element to be set which holds the direct
    // neighbors. Is save to call several times
    VicinityElement::Init(space, structure);
    InitLocal(space);

    if(type_ == SLOPE || type_ == GLOBAL_SLOPE)
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
  this->structure_ = NULL;

  // shortcuts
  Function::Type ftype = func->GetType();
  std::string    fname = Function::type.ToString(ftype);

  // read xml parameters -> might be null valued!
  PtrParamNode pn = func->pn->Get("local", ParamNode::PASS);

  this->beta_ = pn != NULL && pn->Has("beta") ? pn->Get("beta")->As<double>() : -3.14;
  this->beta_ = pn != NULL && pn->Has("eps") ? pn->Get("eps")->As<double>() : -3.14;

  this->normalize_ = pn != NULL ? pn->Get("normalize")->As<bool>() : true;

  // check beta
  switch(ftype)
  {
  case CHECKERBOARD:
  case GLOBAL_CHECKERBOARD:
    if(pn == NULL || !pn->Has("beta"))
      throw Exception("function '" + fname + "' requires the 'beta' attribute in a 'local' element");
    break;

  default:
    if(pn != NULL && pn->Has("beta"))
      throw Exception("function '" + fname + "' does not accept 'beta' attribute");
  }

  // check eps
  switch(ftype)
  {
  case MOLE:
  case GLOBAL_MOLE:
    if(pn == NULL || !pn->Has("eps"))
      throw Exception("function '" + fname + "' requires the 'eps' attribute in a 'local' element");
    break;

  default:
    if(pn != NULL && pn->Has("eps"))
      throw Exception("function '" + fname + "' does not accept 'eps' attribute");
  }

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
      throw Exception("Invalid choice for 'local' in slope constraint");
    break;

  case GLOBAL_SLOPE:
    if(locality_ != NEXT && locality_ != DEFAULT)
      throw Exception("Invalid choice for 'local' with " + fname);
    locality_ = NEXT_AND_REVERSE;
    break;

  case CHECKERBOARD:
  case GLOBAL_CHECKERBOARD:
    if(locality_ != PREV_NEXT_AND_REVERSE && locality_ != DEFAULT)
      throw Exception("Invalid choice for 'local' with " + fname);
    locality_ = PREV_NEXT_AND_REVERSE;
    break;

  case MOLE:
    if(locality_ != DEG_45_STAR && locality_ != DEFAULT)
      throw Exception("Invalid choice for 'local' with " + fname);
    locality_ = DEG_45_STAR;
    break;

  default: // no locality
    assert(false);
    break;
  }

  // this is actually pure constructor work, just extracted to handle function size
  if(locality_ == DEG_45_STAR)
  {
    if(!pn->Get("local"))
      throw Exception("subelement 'local' with neighborhood information mandatory for '" + fname + "'");
    structure_ = new NeighborhoodStructure(this, pn->Get("local"));
    SetupStarLocalityElementMap();
  }
  else
  {
    SetupVirtualElementMap();
  }

  // needs to be set prior CalcSlopeConstraint() as the optimizers need the size
  values.Resize(virtual_elem_map.GetSize(), -1.0);

  // Function::ToInfo() was already called, hence we hace the node
  ToInfo(func->info_);
}

Function::Local::~Local()
{
  if(structure_ != NULL) { delete structure_; structure_ = NULL; }
}

void Function::Local::SetupVirtualElementMap()
{
  // we construct locality_ into reverse, prev and next
  // reverse means we have a REVERSE option which makes two constraints with different signs
  bool reverse = locality_ == NEXT_AND_REVERSE || locality_ == PREV_NEXT_AND_REVERSE;
  bool prev    = locality_ == PREV_NEXT_AND_REVERSE;
  bool next    = true; // always

  unsigned int  dim  = domain->GetGrid()->GetDim();

  element_dimension_ = dim * (reverse ? 2.0 : 1.0);

  virtual_elem_map.Reserve(element_dimension_ * space->GetNumberOfElements());


  int no_sign = Identifier::NO_SIGN;

  // traverse all elements and check for full neighborhood
  space->AssertOneDesignOnly(); // can be extended we use the design from the conditon
  int elems = space->GetNumberOfElements();
  for(int e = 0, ss = elems; e < ss; ++e)
  {
    DesignElement* de = &(space->data[e]);
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
      assert(next);
      // for the slope constraint bounding box no sign is necessary!
      virtual_elem_map.Push_back(Identifier(de, prev ? ve->GetNeighbour(VicinityElement::X_N) : NULL, ve->GetNeighbour(VicinityElement::X_P), reverse ? 1 : no_sign));
      if(reverse)
        virtual_elem_map.Push_back(Identifier(de, prev ? ve->GetNeighbour(VicinityElement::X_N) : NULL, ve->GetNeighbour(VicinityElement::X_P), -1));

      virtual_elem_map.Push_back(Identifier(de, prev ? ve->GetNeighbour(VicinityElement::Y_N) : NULL, ve->GetNeighbour(VicinityElement::Y_P), reverse ? 1 : no_sign));
      if(reverse)
        virtual_elem_map.Push_back(Identifier(de, prev ? ve->GetNeighbour(VicinityElement::Y_N) : NULL, ve->GetNeighbour(VicinityElement::Y_P), -1));

      if(dim == 3)
        virtual_elem_map.Push_back(Identifier(de, prev ? ve->GetNeighbour(VicinityElement::Z_N) : NULL, ve->GetNeighbour(VicinityElement::Z_P), reverse ? 1 : no_sign));
      if(dim == 3 && reverse)
        virtual_elem_map.Push_back(Identifier(de, prev ? ve->GetNeighbour(VicinityElement::Z_N) : NULL, ve->GetNeighbour(VicinityElement::Z_P), -1));
    }
  }
}

void Function::Local::SetupStarLocalityElementMap()
{
  unsigned int dim  = domain->GetGrid()->GetDim();
  assert(locality_ == DEG_45_STAR);
  assert(structure_ != NULL);
  NeighborhoodStructure* struc = structure_;

  element_dimension_ = dim == 2 ? 4 : 13; // see paper
  virtual_elem_map.Reserve(element_dimension_ * space->GetNumberOfElements());

  space->AssertOneDesignOnly(); // can be extended we use the design from the conditon
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
      unsigned int n = struc->orthogonal[VicinityElement::ToMainAxis( (VicinityElement::Neighbour) dir )];
      if(!VicinityElement::HasNeighbor(de, (VicinityElement::Neighbour) dir, n));
        full = false;
    }

    if(!full) continue;
    // orthogonal first
    StdVector<DesignElement*> buddies;

    for(int dir = VicinityElement::X_P; dir <= (dim == 2 ? VicinityElement::Y_N : VicinityElement::Z_N); dir += 2)
    {
      VicinityElement::Neighbour pos = (VicinityElement::Neighbour) (dir +1);
      VicinityElement::Neighbour neg = (VicinityElement::Neighbour) dir;
      unsigned int a = VicinityElement::ToMainAxis(pos);
      unsigned int n = struc->orthogonal[a];

      for(unsigned int i = n; i > 0; i--)  buddies.Push_back(VicinityElement::GetNeighbour(de, neg, i));
      for(unsigned int i = 1; i <= n; i++) buddies.Push_back(VicinityElement::GetNeighbour(de, pos, i));
      virtual_elem_map.Push_back(Identifier(de, buddies));
      LOG_DBG3(func) << "L:SSLEM: de=" << de->ToString() << " dir=" << dir << " pos=" << pos << " neg=" << neg << " a=" << a << " n=" << n << " buddies=" << buddies.ToString();
    }

    // diagonal
    // from diagonal we know the number of element.
    // from orthogonal we know the ratio of dx and dy
    // -> evalutate how many dx and dy steps we need to identify a diagonal element
    unsigned int n = struc->diagonal[0];
    // double dx = struc->orthogonal[0] / (double) struc->orthogonal[1];
    // double dy = struc->orthogonal[1] / (double) struc->orthogonal[0];
    // TODO!! we assume dx = dy for simplicity!!! :( Do it better and smarter!!
    for(int i = n; i > 0; i--)
    {
      DesignElement* tmp = VicinityElement::GetNeighbour(de, VicinityElement::X_N, i);
      buddies.Push_back(VicinityElement::GetNeighbour(tmp, VicinityElement::Y_N, i));
    }
    for(unsigned int i = 1; i <= n; i++)
    {
      DesignElement* tmp = VicinityElement::GetNeighbour(de, VicinityElement::X_P, i);
      buddies.Push_back(VicinityElement::GetNeighbour(tmp, VicinityElement::Y_P, i));
    }
    // no 3D implemented yet
    assert(dim != 3);
  }

}


void Function::Local::ToInfo(PtrParamNode in)
{
  in->Get("locality")->SetValue(locality.ToString(locality_));
  in->Get("local_size")->SetValue(virtual_elem_map.GetSize());
  if(func_->type_ == GLOBAL_SLOPE || func_->type_ == GLOBAL_CHECKERBOARD)
    in->Get("normalize")->SetValue(normalize_);
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
    orthogonal[i] = (int) ((radius / edges[0]) + 0.5); // proper rounding

  // validate
  for(unsigned int i = 0; i < orthogonal.GetSize(); i++)
    if(orthogonal[i] == 0)
      throw Exception("your local neighbor radius for '" + Function::type.ToString(local->func_->GetType()) + "' is too small");

  // now diagonal
  if(dim == 2)
  {
    double diag = std::sqrt(edges[0] * edges[0] + edges[1] * edges[1]);
    diagonal.Resize(1);
    diagonal[0] = (int) ((radius / diag) + 0.5);
  }
  else
  {
    assert(false);
  }
}

void Function::Local::NeighborhoodStructure::ToInfo(PtrParamNode in)
{
  in->Get("type")->SetValue(DesignStructure::filterSpace.ToString(fs));
  in->Get("value")->SetValue(value);
  in->Get("radius")->SetValue(radius);
  in->Get("x_elems")->SetValue(orthogonal[0] * 2 + 1);
  in->Get("y_elems")->SetValue(orthogonal[1] * 2 + 1);
  if(orthogonal.GetSize() == 3)
    in->Get("z_elems")->SetValue(orthogonal[2] * 2 + 1);
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
  this->sign = si;
}


double Function::Local::Identifier::EvalFunction(const Local* local) const
{
  // function value
  double fv = 0.0;
  Function* f = local->func_;

  switch(f->type_)
  {
  case SLOPE:
  case GLOBAL_SLOPE:
    fv = CalcSlope();
    break;

  case CHECKERBOARD:
  case GLOBAL_CHECKERBOARD:
    fv = CalcCheckerboard(local->GetBeta());
    break;

  case MOLE:
  case GLOBAL_MOLE:
    fv = CalcMole(local->GetEps());
    break;

  default:
    assert(false);
  }

  LOG_DBG2(func) << "L:I:EF: f=" << f->type.ToString(f->type_)
                 << " de=" << element->elem->elemNum << " sign=" << sign << " fv=" << fv;

  // handle globalization
  switch(f->type_)
  {
  case GLOBAL_SLOPE:
  case GLOBAL_CHECKERBOARD:
  case GLOBAL_MOLE:
  {
    // we normalize all values by the number of "constraints". Note that it is
    // sufficient for the function value, the gradient is then also right
    double factor = local->DoNormalizeGlobal() ? 1.0/local->virtual_elem_map.GetSize() : 1.0;

    double v = std::max(0.0, fv - f->GetParameter());

    double res = factor * v*v;

    LOG_DBG2(func) << "L:I:EF: global! bound=" << f->GetParameter() << " factor=" << factor
                   << " power=" << (v*v) << " -> " << res;

    return res;
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

  // are we global? If so, we need the function value and don't need to to anything if
  // this function value is zero
  double fv = 0.0;
  bool global = true;
  // we need the real function values, EvalFunction would give the one for global!!
  switch(ft)
  {
  case GLOBAL_SLOPE:
    fv = std::max(0.0, CalcSlope() - funct->GetParameter());
    break;
  case GLOBAL_CHECKERBOARD:
    fv = std::max(0.0, CalcCheckerboard(local->beta_) - funct->GetParameter());
    break;
  case GLOBAL_MOLE:
    fv = std::max(0.0, CalcMole(local->eps_) - funct->GetParameter());
    break;
  default:
    global = false;
  }
  if(global && fv == 0.0)
  {
    LOG_DBG2(func) << "L:I:EvalGrad: f=" << funct->type.ToString(funct->type_) << " de="
                   << element->elem->elemNum << " sign=" << sign << " fv=0.0 -> return immediately";
    return;
  }

  for(int n = -1, nn = neighbor.GetSize(); n < nn; n++)
  {
    double gv = -5.0;

    switch(ft)
    {
    case SLOPE:
    case GLOBAL_SLOPE:
      gv = CalcSlopeGradient(n);
      break;

    case CHECKERBOARD:
    case GLOBAL_CHECKERBOARD:
      gv = CalcCheckerboardGradient(n, local->beta_);
      break;

    case MOLE:
    case GLOBAL_MOLE:
      gv = CalcMoleGradient(n, local->eps_);
      break;

    default:
      assert(false);
    }

    LOG_DBG2(func) << "L:I:EvalGrad: f=" << funct->type.ToString(funct->type_) << " de="
                   << element->elem->elemNum << " sign=" << sign << " n=" << n
                   << " curr=" << GetElement(n)->elem->elemNum << " gv=" << gv;

    // post process the globalized functions
    if(global)
    {
      double factor = local->DoNormalizeGlobal() ? 1.0/local->virtual_elem_map.GetSize() : 1.0;

      gv = factor * 2.0 * fv * gv;
      LOG_DBG2(func) << "L:I:EvalGrad: f=" << funct->type.ToString(funct->type_) << " de="
                     << element->elem->elemNum << " sign=" << sign << " n=" << n
                     << " curr=" << GetElement(n)->elem->elemNum << " gv=" << gv
                     << " bound! factor=" << factor << " fv=" << fv << " new gv=" << gv;
    }

    DesignElement* de = GetElement(n);
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


double Function::Local::Identifier::CalcCheckerboard(double beta) const
{
  double own  = element->GetDesign(DesignElement::SMART);
  assert(neighbor.GetSize() == 2);
  double prev = neighbor[0]->GetDesign(DesignElement::SMART);
  double next = neighbor[1]->GetDesign(DesignElement::SMART);

  assert(sign == 1 || sign == -1);
  if(sign == 1)
  {
    // "Heaviside"(rho_i - max( rho_i-1, rho_i+1)
    // double smaller = std::max(0.0, own - std::max(left, right));
    double max_left_right = beta < 0 ? std::max(prev, next) : SmoothMax(prev, next, beta);
    double smaller = std::max(0.0, own - max_left_right);
    return smaller;
  }
  else
  {
    // "Heaviside"(min( rho_i-1, rho_i+1) - rho_i)
    // double larger = std::max(0.0, std::min(left, right) - own);
    double min_left_right = beta < 0 ? std::min(prev, next) : SmoothMin(prev, next, beta);
    double larger = std::max(0.0, min_left_right - own);
    return larger;
  }
}

double Function::Local::Identifier::CalcCheckerboardGradient(int neigh_idx, double beta)
{
  double prev_val = GetElement(0)->GetDesign(DesignElement::SMART);
  double next_val = GetElement(1)->GetDesign(DesignElement::SMART);

  double exp_prev = 0.0;
  double exp_next = 0.0;

  double result = -125.56;

  if(sign == 1)
  {
    if(neigh_idx == -1) result = 1.0;
    else
    {
      exp_prev = std::exp(prev_val * beta);
      exp_next = std::exp(next_val * beta);
    }
  }
  else // sign = -1
  {
    if(neigh_idx == -1) result = -1.0;
    else
    {
      exp_prev = std::exp((1.0 - prev_val) * beta);
      exp_next = std::exp((1.0 - next_val) * beta);
    }
  }

  if(neigh_idx == 0)
    result = -1.0 * exp_prev / (exp_prev + exp_next);
  if(neigh_idx == 1)
    result = -1.0 * exp_next / (exp_prev + exp_next);

  LOG_DBG3(func) << "CCG el=" << element->elem->elemNum << " neigh_idx=" << neigh_idx
                     << " prev_val=" << prev_val << " next_val=" << next_val << " exp_prev=" << exp_prev
                     << " exp_next=" << exp_next << " result=" << result;

  return result;
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
  tmp.Resize(0); // is static and keeps capacity
  for(unsigned int i = 0; i < neighbor.GetSize() / 2; i++)
    tmp.Push_back(GetElement(i)->GetDesign(DesignElement::SMART));

  tmp.Push_back(element->GetDesign(DesignElement::SMART));

  for(unsigned int i = neighbor.GetSize() / 2; i < neighbor.GetSize(); i++)
      tmp.Push_back(GetElement(i)->GetDesign(DesignElement::SMART));

  assert(tmp.GetSize() == neighbor.GetSize() + 1);

  double sum = 0.0;
  for(unsigned int i = 0; i < tmp.GetSize() - 1; i++)
    sum += SmoothAbs(tmp[i+1] - tmp[i], eps);
  double result = sum - SmoothAbs(tmp[0] - tmp.Last(), eps);
  return result;
}

double Function::Local::Identifier::CalcMoleGradient(int neigh_idx, double eps)
{
  // see comments in the forward function implementation
  // three cases for the sensitivity analysis: first, intermediate, last

  // sort all values in a sequence -> TODO! optimize and do not copy and pase
  tmp.Resize(0); // is static and keeps capacity
  for(unsigned int i = 0; i < neighbor.GetSize() / 2; i++)
    tmp.Push_back(GetElement(i)->GetDesign(DesignElement::SMART));

  tmp.Push_back(element->GetDesign(DesignElement::SMART));

  for(unsigned int i = neighbor.GetSize() / 2; i < neighbor.GetSize(); i++)
      tmp.Push_back(GetElement(i)->GetDesign(DesignElement::SMART));

  assert(tmp.GetSize() == neighbor.GetSize() + 1);

  // correct the index
  unsigned int half = neighbor.GetSize() / 2;
  int idx = neigh_idx == -1 ? half : neigh_idx < (int) half ?  neigh_idx : neigh_idx + 1;
  assert(idx >= 0 && idx < (int) tmp.GetSize());

  if(idx == 0)
    return DerivSmoothAbs(tmp.Last() - tmp[0], eps) - DerivSmoothAbs(tmp[1] - tmp[0], eps);
  if(idx == (int) tmp.GetSize() - 1)
    return DerivSmoothAbs(tmp.Last() - tmp[idx-1], eps) - DerivSmoothAbs(tmp.Last() - tmp[0], eps);
  else
    return DerivSmoothAbs(tmp[idx] - tmp[idx-1], eps) - DerivSmoothAbs(tmp[idx+1] + tmp[idx], eps);
}
