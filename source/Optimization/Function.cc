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
Enum<Function::Local::Locality> Function::Local::locality;

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

  // shortcuts
  Function::Type ftype = func->GetType();
  std::string    fname = Function::type.ToString(ftype);

  // read xml parameters -> might be null valued!
  PtrParamNode pn = func->pn->Get("local", ParamNode::PASS);

  this->beta_ = pn != NULL && pn->Has("beta") ? pn->Get("beta")->As<double>() : -3.14;

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


  VicinityElement::Neighbour none = VicinityElement::NONE;
  int no_sign = Identifier::NO_SIGN;

  // traverse all elements and check for full neighborhood
  space->AssertOneDesignOnly(); // can be extended we use the design from the conditon
  int elems = space->GetNumberOfElements();
  for(int e = 0, ss = elems; e < ss; ++e)
  {
    DesignElement& de = space->data[e];

    // do we have a full neighborhood? All or none as in the original slope paper
    bool full = true;
    if(prev)
    {
      if(de.vicinity->design[VicinityElement::X_N] == NULL) full = false;
      if(de.vicinity->design[VicinityElement::Y_N] == NULL) full = false;
      if(dim == 3 && de.vicinity->design[VicinityElement::Z_N] == NULL) full = false;
    }
    if(next)
    {
      if(de.vicinity->design[VicinityElement::X_P] == NULL) full = false;
      if(de.vicinity->design[VicinityElement::Y_P] == NULL) full = false;
      if(dim == 3 && de.vicinity->design[VicinityElement::Z_P] == NULL) full = false;
    }

    LOG_DBG2(ofunc) << "Local::Local e_num=" << de.elem->elemNum << " vicinity=" << de.vicinity->ToString() << " full=" << full;

    if(full)
    {
      assert(next);
      // for the slope constraint bounding box no sign is necessary!
      virtual_elem_map.Push_back(Identifier(e, prev ? VicinityElement::X_N : none, VicinityElement::X_P, reverse ? 1 : no_sign));
      if(reverse)
        virtual_elem_map.Push_back(Identifier(e, prev ? VicinityElement::X_N : none, VicinityElement::X_P, -1));

      virtual_elem_map.Push_back(Identifier(e, prev ? VicinityElement::Y_N : none, VicinityElement::Y_P, reverse ? 1 : no_sign));
      if(reverse)
        virtual_elem_map.Push_back(Identifier(e, prev ? VicinityElement::Y_N : none, VicinityElement::Y_P, -1));

      if(dim == 3)
        virtual_elem_map.Push_back(Identifier(e, prev ? VicinityElement::Z_N : none, VicinityElement::Z_P, reverse ? 1 : no_sign));
      if(dim == 3 && reverse)
        virtual_elem_map.Push_back(Identifier(e, prev ? VicinityElement::Z_N : none, VicinityElement::Z_P, -1));
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
    DesignElement* curr = de;

    // do we have a full neighborhood? All or none
    bool full = true;

    // we only need to check orthogonal
    assert(VicinityElement::X_P == 0);
    for(int dir = VicinityElement::X_P; dir <= (dim == 2 ? VicinityElement::Y_N : VicinityElement::Z_N); dir++)
    {
      curr = de;
      for(unsigned int n = 0, nn = struc->orthogonal[VicinityElement::ToMainAxis((VicinityElement::Neighbour) dir)]; n < nn && full; n++)
      {
        curr = curr->vicinity->design[dir];
        if(curr == NULL) full = false;
      }
    }

    if(!full) continue;
    // orthogonal first
    StdVector<VicinityElement::Neighbour> buddies;

    unsigned int n = struc->orthogonal[0];
    for(unsigned int i = 0; i < n; i++)  buddies.Push_back(VicinityElement::X_N);
    for(unsigned int i = 0; i < n; i++)  buddies.Push_back(VicinityElement::X_P);
    virtual_elem_map.Push_back(Identifier(e, buddies));

    buddies.Resize(0);
    n = struc->orthogonal[1];
    for(unsigned int i = 0; i < n; i++)  buddies.Push_back(VicinityElement::Y_N);
    for(unsigned int i = 0; i < n; i++)  buddies.Push_back(VicinityElement::Y_P);
    virtual_elem_map.Push_back(Identifier(e, buddies));

    if(dim == 3)
    {
      buddies.Resize(0);
      n = struc->orthogonal[2];
      for(unsigned int i = 0; i < n; i++)  buddies.Push_back(VicinityElement::Z_N);
      for(unsigned int i = 0; i < n; i++)  buddies.Push_back(VicinityElement::Z_P);
      virtual_elem_map.Push_back(Identifier(e, buddies));
    }

    // first diagonal
    n = struc->diagonal[0];
    for(unsigned int i = 0; i < n; i++)  buddies.Push_back(VicinityElement::Y_N);
    for(unsigned int i = 0; i < n; i++)  buddies.Push_back(VicinityElement::Y_P);

    if(dim == 3)
    {
      assert(false);
    }

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

Function::Local::Identifier::Identifier(unsigned int el_idx, VicinityElement::Neighbour prev, VicinityElement::Neighbour next, int si)
{
  this->element_idx = el_idx;

  assert(next != VicinityElement::NONE);
  bool has_prev = prev != VicinityElement::NONE;

  this->neighbor.Resize(has_prev ? 2 : 1);

  if(has_prev)
    this->neighbor[0] = prev;

  this->neighbor[has_prev ? 1 : 0] = next;

  this->sign = si;
}

Function::Local::Identifier::Identifier(unsigned int el_idx, StdVector<VicinityElement::Neighbour> buddies, int si)
{
  this->element_idx = el_idx;
  this->neighbor = buddies;
  this->sign = si;
}

DesignElement* Function::Local::Identifier::GetElement(DesignSpace* design, int neigh_idx)
{
  if(neigh_idx == -1)
    return &(design->data[element_idx]);
  else
    return design->data[element_idx].vicinity->GetNeighbour(neighbor[neigh_idx]);
}

double Function::Local::Identifier::EvalFunction(const DesignSpace* design, const Local* local) const
{
  switch(local->func_->type_)
  {
  case SLOPE:
  case GLOBAL_SLOPE:
    return CalcSlope(design);

  case CHECKERBOARD:
  case GLOBAL_CHECKERBOARD:
    return CalcCheckerboard(design, local->beta_);

  default:
    assert(false);
  }
  return -1.0; // must not happen
}

void Function::Local::Identifier::EvalGradient(DesignSpace* design, const Local* local)
{
  // TODO the dynamic_cast might be to slow, check!
  Condition* g = dynamic_cast<Condition*>(local->func_);
  assert(g != NULL);

  for(int n = -1, nn = neighbor.GetSize(); n < nn; n++)
  {
    // reset the constraint data. Note, as we are local, there are no side effects by elements
    GetElement(design, n)->Reset(DesignElement::CONSTRAINT_GRADIENT, g);

    double gv = -5.0;

    switch(local->func_->type_)
    {
    case SLOPE:
      gv = CalcSlopeGradient(n);
      break;

    case CHECKERBOARD:
      gv = CalcCheckerboardGradient(n, design, local->beta_);
      break;

    default:
      assert(false);
    }

    DesignElement* de = GetElement(design, n);
    de->AddGradient(NULL, g, gv);
    LOG_DBG2(ofunc) << "FLI:EG: c=" <<  Function::type.ToString(local->func_->type_) << " de="
                    << design->data[element_idx].elem->elemNum << " sign=" << sign
                    << " n=" << n << " gv=" << gv;
  }
}


double Function::Local::Identifier::CalcSlope(const DesignSpace* design) const
{
  const DesignElement& de = design->data[this->element_idx];
  double mine  = de.GetDesign(DesignElement::SMART);
  assert(this->neighbor.GetSize() == 1);
  double other = de.vicinity->GetNeighbour(this->neighbor[0])->GetDesign(DesignElement::SMART);

  double s = this->sign == -1 ? -1.0 : 1.0;

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


double Function::Local::Identifier::CalcCheckerboard(const DesignSpace* design, double beta) const
{

  const DesignElement& de = design->data[element_idx];
  double own  = de.GetDesign(DesignElement::SMART);
  assert(neighbor.GetSize() == 2);
  double prev = de.vicinity->GetNeighbour(neighbor[0])->GetDesign(DesignElement::SMART);
  double next = de.vicinity->GetNeighbour(neighbor[1])->GetDesign(DesignElement::SMART);

  assert(sign == 1 || sign == -1);
  if(sign == 1)
  {
    // "Heaviside"(rho_i - max( rho_i-1, rho_i+1)
    // double smaller = std::max(0.0, own - std::max(left, right));
    double max_left_right = beta < 0 ? std::max(prev, next) : CalcMaxApproximation(prev, next, beta);
    double smaller = std::max(0.0, own - max_left_right);
    return smaller;
  }
  else
  {
    // "Heaviside"(min( rho_i-1, rho_i+1) - rho_i)
    // double larger = std::max(0.0, std::min(left, right) - own);
    double min_left_right = beta < 0 ? std::min(prev, next) : CalcMinApproximation(prev, next, beta);
    double larger = std::max(0.0, min_left_right - own);
    return larger;
  }
}

double Function::Local::Identifier::CalcCheckerboardGradient(int neigh_idx, DesignSpace* design,  double beta)
{
  double prev_val = GetElement(design, 0)->GetDesign(DesignElement::SMART);
  double next_val = GetElement(design, 1)->GetDesign(DesignElement::SMART);

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

  LOG_DBG3(ofunc) << "CCG el=" << GetElement(design, -1)->elem->elemNum << " neigh_idx=" << neigh_idx
                     << " prev_val=" << prev_val << " next_val=" << next_val << " exp_prev=" << exp_prev
                     << " exp_next=" << exp_next << " result=" << result;

  return result;
}


