#include "General/exception.hh"
#include "Optimization/DesignElement.hh"
#include "Optimization/DesignSpace.hh"
#include "Optimization/Condition.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"
#include "Elements/basefe.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/InfoNode.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "boost/lexical_cast.hpp"

using namespace std;
using namespace CoupledField;

DECLARE_LOG(del)
DEFINE_LOG(del, "designElement")

// the static enum
Enum<DesignElement::Filter>         DesignElement::filter;
Enum<DesignElement::Type>           DesignElement::type;
Enum<DesignElement::ValueSpecifier> DesignElement::valueSpecifier;
Enum<DesignElement::Access>         DesignElement::access;
Enum<DesignElement::Detail>         DesignElement::detail;

/** The default constructor for StdVector and ghost elements*/
DesignElement::DesignElement()
{
  Init();
}

DesignElement::DesignElement(ParamNode* pn, Elem* elem)
{
  Init();
  this->elem      = elem;

  // it is a little slow to perform this code for every DesignElement but the
  // implementations are rater fast and it should be not measurable in the end
  type_ = type.Parse(pn->Get("name"));

  upper_ = 1.0;
  // eventually overwrite
  pn->Get("upper", upper_, false);

  lower_ = type_ == POLARIZATION ? -1.0 : 0.001;
  pn->Get("lower", lower_, false);
}


DesignElement::~DesignElement()
{
  if(simp != NULL)     { delete simp; simp = NULL; }
  if(lsn != NULL)      { delete lsn; lsn = NULL; }
  if(vicinity != NULL) { delete vicinity; vicinity = NULL; }
  if(location_ != NULL) { delete location_; location_ = NULL; }
}

void DesignElement::Init()
{
  design          = 0.0;
  cost_gradient   = -1.0;
  simp            = NULL;
  vicinity        = NULL;
  lsn             = NULL;
  location_       = NULL;
  elem            = NULL;
  type_           = NO_TYPE;
  upper_          = 0.0;
  lower_          = 0.0;
}

void DesignElement::SetLocation(DesignElement* ref, VicinityElement::Neighbour pos1, VicinityElement::Neighbour pos2)
{
  this->location_ = new Point(*(ref->GetLocation()));
  
  // find dimension of reference
  Matrix<double>  coords; 
  StdVector<double> edges;
  // assume uniform grid!!
  domain->GetGrid()->GetElemNodesCoord(coords, ref->elem->connect, false );
  ref->elem->ptElem->GetEdgeLength(coords, edges);
  
  // a small switch is easier to read than a formula
  switch(pos1)
  {
    case VicinityElement::X_P: 
      location_->data[0] += edges[0];
      break;
      
    case VicinityElement::X_N:
      location_->data[0] -= edges[0];
      break;

    case VicinityElement::Y_P: 
      location_->data[1] += edges[1];
      break;
      
    case VicinityElement::Y_N:
      location_->data[1] -= edges[1];
      break;
      
    case VicinityElement::Z_P: 
      location_->data[2] += edges[2];
      break;
      
    case VicinityElement::Z_N:
      location_->data[2] -= edges[2];
      break;
      
    default: assert(false);
  }

  // todo: clean up dirty copy & pase
  switch(pos2)
  {
    case VicinityElement::NONE:
      break; // do nothinh
  
    case VicinityElement::X_P: 
      location_->data[0] += edges[0];
      break;
      
    case VicinityElement::X_N:
      location_->data[0] -= edges[0];
      break;

    case VicinityElement::Y_P: 
      location_->data[1] += edges[1];
      break;
      
    case VicinityElement::Y_N:
      location_->data[1] -= edges[1];
      break;
      
    case VicinityElement::Z_P: 
      location_->data[2] += edges[2];
      break;
      
    case VicinityElement::Z_N:
      location_->data[2] -= edges[2];
      break;
      
    default: assert(false);
  }

  
  LOG_DBG3(del) << "DesignElement::SetLocation " << ref->GetLocation()->ToString() << " pos1=" << pos1 << " pos2=" << pos2 << " -> " << location_->ToString();
}  


Point* DesignElement::GetLocation()
{
  if(location_ != NULL) return location_;
  
  location_ = new Point();
  // calc location
  
  // check for ghost elements
  if(elem == NULL) EXCEPTION("location_ not set but elem is NULL");
    
  Matrix<double>  coords; // we ignore the n times constructs

  StdVector<unsigned int>& connect = elem->connect;
  // do not use updated coordinates up to now!!
  domain->GetGrid()->GetElemNodesCoord(coords, connect, false );
  // a barycenter is simply the average of all coordinates -> location_ gets the Point out
  BaseFE::CalcBarycenter(coords, *location_);

  LOG_DBG3(del) << "DesignElement::GetLocation() find " << location_->ToString() << " for " << ToString();
  return location_;
}

void DesignElement::SetConstraintGradient(const Condition* condition, double value)
{
  this->constraintGradient[condition->GetIndex()] = value;
}

double DesignElement::GetConstraintGradient(const Condition* condition) const
{
  return this->constraintGradient[condition->GetIndex()];
}


void DesignElement::GetValue(ResultDescription& rd, StdVector<double>& out, unsigned int dofs) const
{
  // check for special result
  if(rd.value == OBJECTIVE || (rd.value == COST_GRADIENT && rd.detail != NONE))
  {
    if(dofs != 1) throw Exception("special results is only defined for scalar values");
    switch(rd.solutionType)
    {
    case OPT_RESULT_1:
      out[0] = specialResult[0];
      break;
      
    case OPT_RESULT_2:
      out[0] = specialResult[1];
      break;
      
    case OPT_RESULT_3:
      out[0] = specialResult[2];
      break;

    case OPT_RESULT_4:
      out[0] = specialResult[3];
      break;

    case OPT_RESULT_5:
      out[0] = specialResult[4];
      break;

    case OPT_RESULT_6:
      out[0] = specialResult[5];
      break;
      
    default: throw Exception("solution type not handled");
    }
  }
  else
  {
    if(dofs == 1)
      out[0] = GetValue(rd.value, rd.access);
    else
    {
      if(rd.access != PLAIN) throw Exception("non-scalar result '" + valueSpecifier.ToString(rd.value) + "' requicess access='plain'");
      GetValue(rd.value, out);
    }
  }
  return;
}

double DesignElement::GetValue(ValueSpecifier vs, Access access) const
{
  double val = 0.0;
  // we are silent if one wants filtering and it is not possible // todo: info
  if(simp == NULL || access == PLAIN)
    val = GetValue(vs);
  else
    val = simp->GetFilteredValue(vs, false);

  LOG_DBG2(del) << elem->elemNum << " GetValue(" << valueSpecifier.ToString(vs) << ", " 
  << access << " -> " << val;
  return val;
}



void DesignElement::GetValue(ValueSpecifier sp, StdVector<double>& out) const
{
  if(sp == LEVEL_SET_NORMAL)
    if(lsn == NULL) throw Exception("'" + valueSpecifier.ToString(sp) + "' is specific to Level-Set");
  
  switch(sp)
  {
  case LEVEL_SET_NORMAL:
    out.Import(lsn->normal.GetPointer(), lsn->normal.GetSize());
    break;
    
  default: 
    assert(false);
  }
}


double DesignElement::GetValue(ValueSpecifier sp) const
{
  // validate first:
  if(sp == LEVEL_SET_VALUE)
    if(lsn == NULL) throw Exception("'" + valueSpecifier.ToString(sp) + "' is specific to Level-Set");
  switch(sp)
  {
  case DESIGN:               return design;
  case DESIGN_COST_GRADIENT: return design * cost_gradient;
  case COST_GRADIENT:        return cost_gradient;
  case WEIGHT:
    if(simp == NULL) throw Exception("'" + valueSpecifier.ToString(sp) + "' is specific to SIMP");
    return simp->weight;
  case NUM_NEIGHBOURS:       
    if(simp == NULL) throw Exception("'" + valueSpecifier.ToString(sp) + "' is specific to SIMP");
    return simp->neighbourhood.GetSize();
  case CONSTRAINT_GRADIENT:  
    throw Exception("for constraint gradient we need an index!");
  case LEVEL_SET_VALUE:
    return lsn->GetValue();
  default: throw Exception(valueSpecifier.ToString(sp) + " is no scalar value"); 
  }
}

double DesignElement::GetObjectiveGradient(Access access) const
{
  double val = 0;

  if(simp == NULL || access == PLAIN)
    val = GetValue(COST_GRADIENT);
  else
    // see the filter definition for gradients in the
    // 99-lines paper why we use here "design TIMES cost_gradient"
    val = simp->GetFilteredValue(COST_GRADIENT, true);

  LOG_DBG2(del) << elem->elemNum << " GetObjectiveGradient() -> " << val;
  return val;
}

double DesignElement::GetDesign(Access access) const
{
  if(simp == NULL || access == PLAIN)
    return design;
  else
    return simp->GetFilteredValue(DESIGN, false);
}


void DesignElement::ToInfo(InfoNode* in) const
{
  in->Get("type")->SetValue(type.ToString(type_));
  in->Get("upperBound")->SetValue(upper_);
  in->Get("lowerBound")->SetValue(lower_);
}

std::string DesignElement::ToString(DesignElement* de)
{
  if(de == NULL) return "null";

  std::stringstream ss;
  if(de->elem == NULL)
  {
    ss << "ghost";
    if(de->vicinity != NULL) ss << de->vicinity->ToString();
  }
  else ss << boost::lexical_cast<std::string>(de->elem->elemNum);
  
  return ss.str();
}

void DesignElement::SetEnums()
{
  filter.SetName("DesignElement::Filter");
  filter.Add(RADIUS, "radius");
  filter.Add(VOLUME_RADIUS, "volumeRadius");
  filter.Add(MAX_EDGE, "maxEdge");

   type.SetName("DesignElement::Type");
   type.Add(TENSOR_TRACE, "tensor_trace");
   type.Add(DEFAULT, "default");
   type.Add(DENSITY, "density");
   type.Add(POLARIZATION, "polarization");
   type.Add(EMODUL, "emodul");
   type.Add(POISSON, "poisson");
   type.Add(LAMELAMBDA, "lamelambda");
   type.Add(LAMEMU, "lamemu");

  access.SetName("DesignElement::Access");
  access.Add(PLAIN, "plain");
  access.Add(SMART, "smart");

  valueSpecifier.SetName("DesignElement::ValueSpecifier");
  valueSpecifier.Add(DESIGN, "design");
  valueSpecifier.Add(DESIGN_COST_GRADIENT, "designTimesCostGradient");
  valueSpecifier.Add(COST_GRADIENT, "costGradient");
  valueSpecifier.Add(CONSTRAINT_GRADIENT, "constraintGradient");
  valueSpecifier.Add(WEIGHT, "weight");
  valueSpecifier.Add(OBJECTIVE, "objective");
  valueSpecifier.Add(NUM_NEIGHBOURS, "neighbours");
  valueSpecifier.Add(LEVEL_SET_VALUE, "levelSetValue");
  valueSpecifier.Add(LEVEL_SET_NORMAL, "levelSetNormal");

  detail.SetName("DesignElement::Detail");
  detail.Add(NONE, "none");
  detail.Add(SYMMETRY, "symmetry");
  detail.Add(FINITE_DIFF_COST_GRADIENT, "finiteDiffCostGrad");
  detail.Add(ERROR_COST_GRADIENT, "finiteDiffCostGradRelError");
  detail.Add(MECH_MECH, "mech_mech");
  detail.Add(ELEC_ELEC, "elec_elec");
  detail.Add(ELEC_ELEC_QUAD, "elec_elec_quad");
  detail.Add(ELEC_MECH, "elec_mech");
  detail.Add(MECH_ELEC, "mech_elec");
}


SIMPElement::SIMPElement(DesignElement* base)
{
  this->de_ = base;
  this->filter_ = false;
  this->weight  = 1.0;
}


void SIMPElement::InitFilter(StdVector<DesignElement>& data, ParamNode* pn)
{
  DesignElement::Filter filter = DesignElement::filter.Parse(pn->Get("type"));
  double value  = pn->Get("value")->AsDouble();

  if(value == 0.0) return; // todo: infoParam

  double avg_radius = 0;
  double avg_neighbours = 0;

  Grid* grid = domain->GetGrid();
  double dim = (double) grid->GetDim();
  Matrix<double>  coords; 
  // assume uniform grid!!
  grid->GetElemNodesCoord(coords, data[0].elem->connect, false );
  
  // look for the neighbours. This is n^2 but Manfred said, that there
  // will be a smarter structure soon in CFS to do it faster.
  // for the same reason GridCFS<DIM>::GetElemsNextToNodes() is not used.
  for(UInt element = 0; element < data.GetSize(); element++)
  {
    DesignElement* de = &data[element];

    de->simp->filter_ = value > 0.0;
    Point& loc = *de->GetLocation();

    double radius, tmp;

    switch(filter)
    {
    case DesignElement::RADIUS:
      radius = value;
      break;

    case DesignElement::VOLUME_RADIUS:
      // TODO really check for axis symmetry off
      tmp = de->elem->ptElem->CalcVolume(coords, false);
      // The radius is <value> times square/cube edge length where the
      // square/cube has the volume of the element
      radius = value * std::pow(tmp, 1.0/dim);
      LOG_DBG3(del) << "from volume " << tmp << " to radius " << radius;
      break;

    case DesignElement::MAX_EDGE:
      de->elem->ptElem->GetMaxMinEdgeLength(coords, radius, tmp);
      radius = value * radius;
      LOG_DBG3(del) << "Element " << de->elem->elemNum << " edge max=" << radius << " min=" << tmp;
      break;

    default:
      assert(false);
    }

    for(UInt n = 0; n < data.GetSize(); n++)
    {
      // an element is not a neighbour of itself!
      if(de->elem->elemNum == data[n].elem->elemNum) continue;

      double distance = Point::dist(loc, *(data[n].GetLocation()));
      if(distance < radius)
      {
        // value is here a double radius
        // this is the implementation from Bendsoe/ Sigmund
        NeighbourElement ne;
        ne.neighbour = &data[n];
        ne.weight    = value - distance;
        ne.distance  = distance;
        de->simp->neighbourhood.Push_back(ne); // let the default copy constructor work!
        // std::cout << "element " << data[element].elem->elemNum << " "
        //          << data[element].location_.ToString() << " has neighbour "
        //          << data[n].elem->elemNum << " " <<  data[n].location_.ToString()
        //          << " distance= " << distance << " weight=" << ne.weight << std::endl;
      }
    } // here the element is calculated

    // now normalize the weights. The weight of this element is by definition 1.0
    double sum = 1.0;
    for(unsigned int e = 0; e < de->simp->neighbourhood.GetSize(); e++)
      sum += de->simp->neighbourhood[e].weight;

    // now normalize all
    de->simp->weight /= sum;
    for(unsigned int e = 0; e < de->simp->neighbourhood.GetSize(); e++)
      de->simp->neighbourhood[e].weight /= sum;

    avg_radius += radius;
    avg_neighbours += de->simp->neighbourhood.GetSize();
  }

  // for direct debug output: determin neighbourhood statistics
  std::cout << "Filter: avg radius: " << (avg_radius / data.GetSize())
  << " avg neighbourhood: " << (avg_neighbours / data.GetSize()) << std::endl;
}

double SIMPElement::GetFilteredValue(DesignElement::ValueSpecifier sp, bool design_weighted) const
{
  // We filter over this element and the neighbours.

  double factor_sum = 0.0;
  double weight_sum = 0.0;

  for(UInt i = 0; i < neighbourhood.GetSize(); i++)
  {
    const NeighbourElement& ne = neighbourhood[i];
    double des_weight = design_weighted ?  ne.neighbour->GetDesign(DesignElement::PLAIN) : 1.0;
    double factor = ne.weight * ne.neighbour->GetValue(sp) * des_weight;
    factor_sum += factor;
    weight_sum += ne.weight;
    LOG_DBG3(del) << de_->elem->elemNum << ": neighbour " << ne.neighbour->elem->elemNum 
    << " weight = " << ne.weight
    << " des_weight " << des_weight
    << " value = " << ne.neighbour->GetValue(sp);
  }
  // add our part :)
  double des_weight = design_weighted ? de_->GetDesign(DesignElement::PLAIN) : 1.0;

  LOG_DBG2(del) << de_->elem->elemNum << ": factor_sum = " << factor_sum << " + " << de_->GetValue(sp)
  << " weight_sum = " << weight_sum << " + " << weight << " -> "
  << ((factor_sum+(de_->GetValue(sp)*des_weight))/((weight_sum+weight)*des_weight)) << " instead of " << de_->GetValue(sp);

  factor_sum += de_->GetValue(sp) * des_weight;
  weight_sum += this->weight;

  // in the 99-lines paper this is inverse of this_value
  double result = factor_sum / (weight_sum * des_weight);

  return result;
}  



void SIMPElement::Dump()
{
  double weight_sum = weight;
  double distance_avg = 0.0;
  for(unsigned int i = 0; i < neighbourhood.GetSize(); i++)
  {
    weight_sum   += neighbourhood[i].weight;
    distance_avg += neighbourhood[i].distance;
  }
  distance_avg /= neighbourhood.GetSize();

  std::cout << "\nelement: " << de_->elem->elemNum << " location " << de_->GetLocation()->ToString() 
            << " weight sum " << weight_sum << " this weight " << weight <<" distance avg " << distance_avg << std::endl;
  for(unsigned int i = 0; i < neighbourhood.GetSize(); i++)
    std::cout << "  n[" << i << "]: elem " << neighbourhood[i].neighbour->elem->elemNum << " location " 
              << neighbourhood[i].neighbour->GetLocation()->ToString() << " dist=" << neighbourhood[i].distance << " w=" << neighbourhood[i].weight << std::endl; 
}


VicinityElement::VicinityElement()
{
  design.Resize(domain->GetGrid()->GetDim() == 2 ? 4 : 6);
  design.Init(NULL);
}


void VicinityElement::Init(DesignSpace* design)
{
  Grid* grid = domain->GetGrid();
  StdVector<DesignElement>& data = design->data;
  
  // Our design region(s)
  StdVector<RegionIdType> regions;
  regions.Push_back(design->GetRegionId());
  
  // result of GetElemsNextToNodes()
  StdVector<Elem*> elems;
  
  // coordinatates of "this" element
  Matrix<double> reference;
  // coordinates of compare element
  Matrix<double> other;
  
  for(unsigned int element = 0; element < data.GetSize(); element++)
  {
    DesignElement& de = data[element];
    de.vicinity = new VicinityElement();
    
    // all nodes of this element
    StdVector<unsigned int>& nodes = de.elem->connect;

    // coordinates of current element, not updated lagrangian
    grid->GetElemNodesCoord(reference, nodes, false );
    
    // all elements that share "our" nodes, includes diagonal and "this" elements
    elems.Resize(0);
    grid->GetElemsNextToNodes(elems, nodes, regions);  
    
    // check the neighbour elements
    for(unsigned int n = 0; n < elems.GetSize(); n++)
    {
      Elem* other_elem = elems[n];
      // ignore "this" element
      if(de.elem->elemNum == other_elem->elemNum) continue;
      
      grid->GetElemNodesCoord(other, other_elem->connect, false );
      
      int orientation; // 0..2 = x..z neighbour
      bool positive;    // +/- x..z neighbour
      bool valid = VicinityElement::IdentifyNeighbor(reference, other, orientation, positive);
      int  idx   = orientation * 2 + (positive ? 0 : 1);
      LOG_DBG3(del) << "VicinityElement::Init: reference=" << de.elem->elemNum << " other=" << other_elem->elemNum 
                    << " vicinity=" << valid << " orientation=" << orientation << " postivie=" << positive << " idx=" << idx;
      
      // set the vicinity element
      if(!valid) continue;
      DesignElement* other_de = design->Find(other_elem->elemNum, DesignElement::DENSITY); // quick and dirty! 
      
      de.vicinity->design[idx] = other_de;
      
    }
  }
}
  

/** Compares the node coordinates of two elements and decides which orientation we have. only face/line (3D/2D) neighbours.
 * "diagonal" elments and identical elements are "false".
 * @param dimension 0, 1, 2 for x, y, z
 * @param positive if other is reference + dimension, false if reference - other
 * @return true if other is a true closest neighbour */
bool VicinityElement::IdentifyNeighbor(Matrix<double>& reference, Matrix<double>& other, int& dimension, bool& positive)
{
  // -------------------------
  // |       |       |       |
  // |   2   |  1    |   8   |
  // |       |  +y   |       |
  // -------------------------
  // |       |       |       |
  // |   3   |  0    |   7   |
  // |  -x   |       |  +x   |
  // -------------------------
  // |       |       |       |
  // |   4   |   5   |  6    |
  // |       |  -y   |       |
  // -------------------------
  // 
  // From the example above, for reference element 0, we return 1,3,5,7
  
  
  // the structure of the matrices is 
  // [x][i]
  // [y][i]
  // [z][i]
  unsigned int dim = reference.GetNumRows();
  unsigned int size  = reference.GetNumCols();

  assert(other.GetNumRows() == dim);
  assert(other.GetNumCols() == size);
  assert(dim == domain->GetGrid()->GetDim());
  if(dim != 2) EXCEPTION("no implementation for 3D yet")
  
  // defaults to make loggin in the false case easier
  dimension = -1;
  positive = false;
  
  // we first test for x, then for y and then for z
  for(unsigned int d = 0; d < dim; d++)
  {
    // assume 2D and d = 0 (=x), then all y-coordinates of 3 and 7 match 0
    
    // probe for the other dimensions
    for(unsigned int p = 0; p < dim; p++)
    {
      if(d == p) continue;

      bool good = true; // for above, only true for 0,3,7
      
      // check all elements
      for(unsigned int elem = 0; elem < size && good; elem++)
      {
        if(reference[p][elem] != other[p][elem]) good = false;
      }
      
      if(!good) continue;
      
      // this element is in the other dimension than equal to the reference element
      // check for positive or negative
      if(reference[d][0] == other[d][0]) return false; // identical version
      dimension = d;
      positive = other[d][0] > reference[d][0];
      return true;
    }
  }
  // nothing found
  return false;
}
    
    
int VicinityElement::GetNumberOfEntries()
{
  int count = 0;
  for(unsigned int i = 0; i < design.GetSize(); i++)
    if(design[i] != NULL) count++;
  
  return count;
}

DesignElement* VicinityElement::GetNeighbour(Neighbour first, Neighbour second)
{
  DesignElement* de = GetNeighbour(first);
  if(de == NULL) return NULL;
          return de->vicinity->GetNeighbour(second); 
}

void VicinityElement::GetFullNeighbourhood(StdVector<DesignElement*>& out)
{
  assert(domain->GetGrid()->GetDim() == 2);
  
  out.Resize(8);
  out[0] = GetNeighbour(X_P);
  out[1] = GetNeighbour(X_P, Y_N);
  out[2] = GetNeighbour(Y_N);
  out[3] = GetNeighbour(X_N, Y_N);
  out[4] = GetNeighbour(X_N);
  out[5] = GetNeighbour(X_N, Y_P);
  out[6] = GetNeighbour(Y_P);
  out[7] = GetNeighbour(X_P, Y_P);
}

string VicinityElement::ToString()
{
  stringstream ss;
  ss << "(";
  unsigned int max = domain->GetGrid()->GetDim() == 2 ? 4 : 6;
  for(unsigned int i = 0; i < max; i++)
  {
    if(design[i] == NULL) ss << "null";
    else
    {
      if(design[i]->elem == NULL) ss << "ghost";
                             else ss << design[i]->elem->elemNum;
    }
    if(i < max-1) ss << ",";
  }
  ss << ")";
  return ss.str();
}


LevelSetNode::LevelSetNode(DesignElement* de, double value)
{
  assert(de != NULL);  
  int dim = domain->GetGrid()->GetDim();
  this->de_     = de;
  this->value_   = value;
  this->domain_boundary = de->vicinity->GetNumberOfEntries() != (dim == 2 ? 4 : 6);
  this->normal.Resize(domain->GetGrid()->GetDim());
  this->first_order_grad.Resize(dim);
  this->second_order_grad.Resize(dim);

  UpdateDesign();
  LOG_DBG3(del) << "LevelSetNode::LevelSetNode: de=" << (de->elem != NULL ? (int) de->elem->elemNum : -1)<< " vicinity=" << de->vicinity->GetNumberOfEntries() 
                << " value=" << value << " design=" << de->GetDesign(DesignElement::PLAIN);
}

LevelSetNode::~LevelSetNode()
{
  // it is our duty to delete ghost-elements, the vicinity has just references.
  for(unsigned int i = 0; i < ghosts_.GetSize(); i++)
  {
    delete ghosts_[i]; ghosts_[i] = NULL; 
  }
}


int LevelSetNode::CreateGhostElements(double value)
{
  int dim = domain->GetGrid()->GetDim();
  // check if, and how many ghost elements we need!
  assert(dim == 2); // 3D not implemented yet
  // only domain boundary elements "contain" ghost elements
  if(!domain_boundary) return ghosts_.GetSize(); 
  // do we have an "edge" (2D: vicinity has 2 instead of 4 elements), then we have to set the edge element
  // as vicinity has only side neighbours
  VicinityElement& vicinity = *(de_->vicinity);
  
  LOG_DBG3(del) << "LevelSetNode::CreateGhostElements: elem: " << de_->ToString() << " vicinity: " << vicinity.ToString();
  
  // check for an corner element
  bool corner = false;
  if(de_->vicinity->GetNumberOfEntries() == 2) // 2D!
  {
    ghosts_.Push_back(new DesignElement());
    ghosts_.Last()->vicinity = new VicinityElement();
    ghosts_.Last()->lsn = new LevelSetNode(ghosts_.Last(), value);
    corner = true;
  }
  
  // check for regualar ghost elements
  for(int d = 0; d < dim; d++)
  {
    for(int p = 0; p <= 1; p++) // positive and negative component
    {
      unsigned int idx = 2 * d + p;
      if(vicinity.design[idx] == NULL)
      {
        DesignElement* de =  new DesignElement(); // is by default a ghost element
        de->vicinity = new VicinityElement();
        vicinity.design[idx] = de; // link to the other/new node
        de->vicinity->design[2 * d + (p == 0 ? 1 : 0)] = de_; // link the other node with "this"
        de->lsn = new LevelSetNode(de, value);
        de->SetLocation(de_, (VicinityElement::Neighbour) idx);
        ghosts_.Push_back(de);
      }
    }
  }

  // do we have a corner, then do the linking and set the corner
  if(corner)
  {
    assert(ghosts_.GetSize() == 3);
    assert(ghosts_[0]->vicinity->GetNumberOfEntries() == 0); // the corner has no entries yet
    // the corner is the first element!
    // in 2D we have a second and third ghost element
    
    // check lower left corner, order is X_P, X_N, Y_P, Y_N
    assert(ghosts_[1]->vicinity->GetNumberOfEntries() == 1); // the one reference to the design element
    if(IsGhost(vicinity.GetNeighbour(VicinityElement::X_N)) && IsGhost(vicinity.GetNeighbour(VicinityElement::Y_N)))
    {
      ghosts_[0]->SetLocation(de_, VicinityElement::X_N, VicinityElement::Y_N);
      ghosts_[1]->vicinity->design[VicinityElement::Y_N] = ghosts_[0];
      ghosts_[0]->vicinity->design[VicinityElement::Y_P] = ghosts_[1]; 
      ghosts_[2]->vicinity->design[VicinityElement::X_N] = ghosts_[0];
      ghosts_[0]->vicinity->design[VicinityElement::X_P] = ghosts_[2]; 
    }
    // check lower right corner
    if(IsGhost(vicinity.GetNeighbour(VicinityElement::X_P)) && IsGhost(vicinity.GetNeighbour(VicinityElement::Y_N)))
    {
      ghosts_[0]->SetLocation(de_, VicinityElement::X_P, VicinityElement::Y_N);      
      ghosts_[1]->vicinity->design[VicinityElement::Y_N] = ghosts_[0];
      ghosts_[0]->vicinity->design[VicinityElement::Y_P] = ghosts_[1]; 
      ghosts_[2]->vicinity->design[VicinityElement::X_P] = ghosts_[0];
      ghosts_[0]->vicinity->design[VicinityElement::X_N] = ghosts_[2]; 
    }
    // check upper right corner
    if(IsGhost(vicinity.GetNeighbour(VicinityElement::X_P)) && IsGhost(vicinity.GetNeighbour(VicinityElement::Y_P)))
    {
      ghosts_[0]->SetLocation(de_, VicinityElement::X_P, VicinityElement::Y_P);      
      ghosts_[1]->vicinity->design[VicinityElement::Y_P] = ghosts_[0];
      ghosts_[0]->vicinity->design[VicinityElement::Y_N] = ghosts_[1]; 
      ghosts_[2]->vicinity->design[VicinityElement::X_P] = ghosts_[0];
      ghosts_[0]->vicinity->design[VicinityElement::X_N] = ghosts_[2]; 
    }
    // check upper left corner
    if(IsGhost(vicinity.GetNeighbour(VicinityElement::X_N)) && IsGhost(vicinity.GetNeighbour(VicinityElement::Y_P)))
    {
      ghosts_[0]->SetLocation(de_, VicinityElement::X_N, VicinityElement::Y_P);      
      ghosts_[1]->vicinity->design[VicinityElement::Y_P] = ghosts_[0];
      ghosts_[0]->vicinity->design[VicinityElement::Y_N] = ghosts_[1]; 
      ghosts_[2]->vicinity->design[VicinityElement::X_N] = ghosts_[0];
      ghosts_[0]->vicinity->design[VicinityElement::X_P] = ghosts_[2]; 
    }
    // the corner shall have now 2 entries
    assert(ghosts_[0]->vicinity->GetNumberOfEntries() == 2);
    LOG_DBG2(del) << "CreateGhostElements: corner: " << ghosts_[0]->ToString() << " 1:" << ghosts_[1]->ToString() << " 2:" << ghosts_[2]->ToString();  
  }
  
  assert(de_->vicinity->GetNumberOfEntries() == (dim == 2 ? 4 : 6));
  return ghosts_.GetSize();
}

DesignElement* LevelSetNode::GetNeighbour(Neighbour idx)
{
  return de_->vicinity->design[idx];
}

DesignElement* LevelSetNode::GetNeighbour(Neighbour first, Neighbour second)
{
  DesignElement* de = de_->vicinity->design[first];
  if(de == NULL) return NULL;
  de = de->vicinity->design[second];
  if(de == NULL) return NULL;
  // else: found
  return de;
}


void LevelSetNode::SetGhostVicinity()
{
  // this method is only to be called on real elements, not on ghost elements itself
  assert(!IsGhost());
  
  assert(domain->GetGrid()->GetDim() == 2); // 3D not implemented yet

  // from SetGhostElement() every ghost element already links to its "parent" element
  // and the corners are set and linked!
  if(IsGhost(GetNeighbour(X_P)))
  {
    StdVector<DesignElement*>& design = GetNeighbour(X_P)->vicinity->design;
    if(design[Y_P] == NULL) design[Y_P] = GetNeighbour(Y_P, X_P);
    if(design[Y_N] == NULL) design[Y_N] = GetNeighbour(Y_N, X_P);
  }
  if(IsGhost(GetNeighbour(X_N)))
  {
    StdVector<DesignElement*>& design = GetNeighbour(X_N)->vicinity->design;
    if(design[Y_P] == NULL) design[Y_P] = GetNeighbour(Y_P, X_N);
    if(design[Y_N] == NULL) design[Y_N] = GetNeighbour(Y_N, X_N);
  }
  if(IsGhost(GetNeighbour(Y_P)))
  {
    StdVector<DesignElement*>& design = GetNeighbour(Y_P)->vicinity->design;
    if(design[X_P] == NULL) design[X_P] = GetNeighbour(X_P, Y_P);
    if(design[X_N] == NULL) design[X_N] = GetNeighbour(X_N, Y_P);
  }
  if(IsGhost(GetNeighbour(Y_N)))
  {
    StdVector<DesignElement*>& design = GetNeighbour(Y_N)->vicinity->design;
    if(design[X_P] == NULL) design[X_P] = GetNeighbour(X_P, Y_N);
    if(design[X_N] == NULL) design[X_N] = GetNeighbour(X_N, Y_N);
  }
}

void LevelSetNode::SetValue(double value, bool mirror_to_ghosts)
{
  value_ = value;
  
  // mirror to the ghosts, is kind of homogenoeous neumann
  for(unsigned int i = 0; mirror_to_ghosts == true && i < ghosts_.GetSize(); i++)
    ghosts_[i]->lsn->value_ = value;
}

bool LevelSetNode::IsGhost(DesignElement* de) 
{ 
  if(de == NULL) return false;
  if(de->lsn == NULL) return false; //  not all elements might be initialized!
  return de->lsn->IsGhost(); 
}




void LevelSetNode::CalcGradients()
{
  VicinityElement* vicinity = de_->vicinity;
  
  // dont' call this on ghost nodes itself, otherwise we have no neighbourhood!
  for(unsigned int d = 0; d < normal.GetSize(); d++)
  {
    double left  = vicinity->GetNeighbour((VicinityElement::Neighbour) (d * 2 + 1))->lsn->GetValue(); 
    double right = vicinity->GetNeighbour((VicinityElement::Neighbour) (d * 2 + 0))->lsn->GetValue();
    double own   = value_;
    
    // we ignore h and assume a constant h = 1
    first_order_grad[d]  = 0.5 * (right - left);
    second_order_grad[d] = right - 2 * own + left; 
  }
  
  LOG_DBG3(del) << "LevelSetNode::CalcGradients: " << de_->ToString() << " first_order_grad=" 
                << first_order_grad.ToString() << " second_order_grad=" << second_order_grad.ToString();
}



void LevelSetNode::UpdateNormal()
{
  // CalcGradients() has to be executed!
  double length = first_order_grad.NormL2();

  for(unsigned int d = 0; d < normal.GetSize(); d++)
    normal[d] = first_order_grad[d] / length;
}

void LevelSetNode::UpdateCurvature()
{
   assert(domain->GetGrid()->GetDim() == 2);
   
   // source: http://en.wikipedia.org/wiki/Curvature
   // 
   // curvature = |(x' y'' - y' x'') / (x'^2 + y'^2)^(3/2)|
   double numerator = first_order_grad[0] * second_order_grad[1] - first_order_grad[1] * second_order_grad[0];
   double tmp       = square(first_order_grad[0]) + square(first_order_grad[1]);
   double denominator = powl(tmp, 3.0/2.0);
   curvature = abs(numerator / denominator);
}

void LevelSetNode::CalcNextFirstOrderHamiltonJacobiStep(double delta_t)
{
  assert(domain->GetGrid()->GetDim() == 2);
  // this is based on forward and backward differences
  VicinityElement* vicinity = de_->vicinity;
  double own = value_;
  // again we assume h=1 it is possible to do it exact with GetEdgeLength() as in DesignElement::SetLocation()
  double d_p_x = vicinity->GetNeighbour(VicinityElement::X_P)->lsn->GetValue() - own;
  double d_n_x = own - vicinity->GetNeighbour(VicinityElement::X_N)->lsn->GetValue();
  double d_p_y = vicinity->GetNeighbour(VicinityElement::Y_P)->lsn->GetValue() - own;
  double d_n_y = own - vicinity->GetNeighbour(VicinityElement::Y_N)->lsn->GetValue();
  
  // the forward differences
  double forward  = square(max(d_n_x, 0.0)) + square(min(d_p_x, 0.0)) 
                  + square(max(d_n_y, 0.0)) + square(min(d_p_y, 0.0));
  forward = sqrt(forward);
  
  double backward = square(max(d_p_x, 0.0)) + square(min(d_n_x, 0.0))
                  + square(max(d_p_y, 0.0)) + square(min(d_n_y, 0.0));
  backward = sqrt(backward);
                  
  double velocity = de_->GetObjectiveGradient(DesignElement::SMART);                
  // it would be faster to evaluate only forward or backward depeneding on the sign of the velocity
  double next = own - delta_t * (max(velocity, 0.0) * forward + min(velocity, 0.0) * backward);       
  // we propagate this value to the ghosts!
  SetValue(next, true);
}


void LevelSetNode::UpdateDesign()
{
  de_->SetDesign(IsSolid() ? de_->GetUpperBound() : de_->GetLowerBound());
}

ResultDescription::ResultDescription()
{
  // does not set all!!
  access = DesignElement::PLAIN;
  value  = DesignElement::DESIGN;
  design = DesignElement::DEFAULT;
}

ResultDescription::ResultDescription(ParamNode* pn)
{
  // killme: do this nice when there is no more environment.hh
  SolutionType st;
  String2Enum(pn->Get("id")->AsString(), st);
  solutionType = st;

  design = DesignElement::DEFAULT;
  if(pn->Has("design"))
    design = DesignElement::type.Parse(pn->Get("design"));

  access = DesignElement::access.Parse(pn->Get("access"));

  value = DesignElement::valueSpecifier.Parse(pn->Get("value"));

  detail = DesignElement::detail.Parse(pn->Get("detail"));
}

std::string ResultDescription::ToString()
{
  // killme: do this nice when there is no more environment.hh
  std::string string;
  Enum2String((SolutionType) solutionType, string);
  return string;
}
