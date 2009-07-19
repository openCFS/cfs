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
#include "boost/date_time/posix_time/posix_time.hpp"
#include "Optimization/LevelSet.hh"
#include "Optimization/ShapeOptimizer.hh"

using namespace std;
using namespace CoupledField;

using boost::posix_time::ptime;
using boost::posix_time::second_clock;
using boost::posix_time::microsec_clock;

DECLARE_LOG(del)
DEFINE_LOG(del, "designElement")

BaseDesignElement::BaseDesignElement() {
  design          = 0.0;
  cost_gradient   = -1.0;
  upper_          = 0.0;
  lower_          = 0.0;
}

void BaseDesignElement::SetConstraintGradient(const Condition* condition, double value)
{
  if(condition)
    this->constraintGradient[condition->GetIndex()] = value;
  else
    this->cost_gradient = value;
}

void BaseDesignElement::SetConstraintGradientOrAddObjectiveGradient(const Condition* condition, double value){
  if(condition){
    this->constraintGradient[condition->GetIndex()] = value;
  }else{
    this->cost_gradient += value;
  }
}

double BaseDesignElement::GetConstraintGradient(const Condition* condition) const
{
  if(condition)
    return this->constraintGradient[condition->GetIndex()];
  else
    return this->cost_gradient;
}

void BaseDesignElement::PostInit(int constraints){
  // resize and init with 0.0 so constraint, which only act on one design variable, do not have to set all others explicitly to zero 
  constraintGradient.Resize(constraints, 0.0);
}



// the static enum
Enum<DesignElement::Filter>         DesignElement::filter;
Enum<DesignElement::Type>           DesignElement::type;
Enum<DesignElement::ValueSpecifier> DesignElement::valueSpecifier;
Enum<DesignElement::Access>         DesignElement::access;
Enum<DesignElement::Detail>         DesignElement::detail;

/** The default constructor for StdVector and ghost elements*/
DesignElement::DesignElement() : BaseDesignElement()
{
  Init();
}

DesignElement::DesignElement(ParamNode* pn, Elem* elem) : BaseDesignElement()
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
  if(simp != NULL)      { delete simp; simp = NULL; }
  if(vicinity_ != NULL) { delete vicinity_; vicinity_ = NULL; }
  if(location_ != NULL) { delete location_; location_ = NULL; }
}

void DesignElement::Init()
{
  simp            = NULL;
  vicinity_       = NULL;
  lse_            = 0;
  location_       = NULL;
  elem            = NULL;
  type_           = NO_TYPE;
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
      throw Exception("we only have one dof, only scalar values");
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


double DesignElement::GetValue(ValueSpecifier sp) const
{
  // validate first:
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
  case TOPGRAD_VALUE:
    if(lse_ == 0) throw Exception("'" + valueSpecifier.ToString(sp) + "' is specific to Levelset");
    return lse_->topGradValue;
  case SHAPEGRAD_VALUE:
    if(lse_ == 0) throw Exception("'" + valueSpecifier.ToString(sp) + "' is specific to Levelset");
    return lse_->shapeGradValue;
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
    if(de->vicinity_ != NULL) ss << de->vicinity_->ToString();
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
  type.Add(EMODULISO, "emodul-iso");
  type.Add(POISSONISO, "poisson-iso");
  type.Add(GMODUL, "gmodul");
  type.Add(UNITY, "unity");

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
  valueSpecifier.Add(LEVEL_SET_STATE, "levelSetState");
  valueSpecifier.Add(TOPGRAD_VALUE, "topGradValue");
  valueSpecifier.Add(SHAPEGRAD_VALUE, "shapeGradValue");
  valueSpecifier.Add(SHAPEGRAD_NODE_VALUE, "shapeGradNodeValue");
  valueSpecifier.Add(LEVEL_SET_GRAD_XP, "levelSetGradXP");
  valueSpecifier.Add(LEVEL_SET_GRAD_XN, "levelSetGradXN");
  valueSpecifier.Add(LEVEL_SET_GRAD_YP, "levelSetGradYP");
  valueSpecifier.Add(LEVEL_SET_GRAD_YN, "levelSetGradYN");
  valueSpecifier.Add(LEVEL_SET_GRAD_ZP, "levelSetGradZP");
  valueSpecifier.Add(LEVEL_SET_GRAD_ZN, "levelSetGradZN");

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
              << neighbourhood[i].neighbour->GetLocation()->ToString() 
              << " dist=" << neighbourhood[i].distance << " w=" << neighbourhood[i].weight << std::endl;
}


VicinityElement::VicinityElement()
{
  design.Resize(domain->GetGrid()->GetDim() == 2 ? 4 : 6, 0);
}


void VicinityElement::Init(DesignSpace* design)
{
  ptime before_time = second_clock::local_time();
  
  cout << "Vicinities " << flush;
  Grid* grid = domain->GetGrid();
  StdVector<DesignElement> &data = design->data;

  // Our design region(s)
  StdVector<RegionIdType> regions;
  regions.Push_back(design->GetRegionId());

  // result of GetElemsNextToNodes()
  StdVector<Elem*> elems;

  // coordinatates of "this" element
  Matrix<double> reference;
  // coordinates of compare element
  Matrix<double> other;

  unsigned int factor(1);
  for(unsigned int element = 0, data_size = data.GetSize(); element < data_size; ++element)
  {
    if(10*element == factor*data_size) //draw a dot every 1/10 elements to cout
    {
      cout << "." << flush;
      ++factor;
    }
    DesignElement& de = data[element];
    de.vicinity_ = new VicinityElement();

    // all nodes of this element
    StdVector<unsigned int> &nodes = de.elem->connect;

    // coordinates of current element, not updated lagrangian
    grid->GetElemNodesCoord(reference, nodes, false );

    // all elements that share "our" nodes, includes diagonal and "this" elements
    elems.Resize(0);
    grid->GetElemsNextToNodes(elems, nodes, regions);

    // check the neighbour elements
    for(unsigned int n = 0, elems_size = elems.GetSize(); n < elems_size; ++n)
    {
      Elem* other_elem = elems[n];
      // ignore "this" element
      if(de.elem->elemNum == other_elem->elemNum) continue;

      grid->GetElemNodesCoord(other, other_elem->connect, false );

      int orientation; // 0..2 = x..z neighbour
      bool positive;    // +/- x..z neighbour
      bool valid = VicinityElement::IdentifyNeighbor(reference, other, orientation, positive);
      
      // set the vicinity element
      if(!valid) continue;
      
      int  idx   = orientation * 2 + (positive ? 0 : 1);
      LOG_DBG3(del) << "VicinityElement::Init: reference=" << de.elem->elemNum << " other=" << other_elem->elemNum 
                    << " vicinity=" << valid << " orientation=" << orientation << " positive=" << positive << " idx=" << idx;

      DesignElement* other_de = design->Find(other_elem->elemNum, DesignElement::DENSITY); // quick and dirty!

      de.vicinity_->design[idx] = other_de;
    }
  }
  cout << "done (" << ShapeOptimizer::GetTimeString(second_clock::local_time() - before_time) << ") " << flush;
}


/** Compares the node coordinates of two elements and decides which orientation we have. only face/line (3D/2D) neighbours.
 * "diagonal" elements and identical elements are "false".
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
  // and i \in {0, ..., number_of_elem_nodes} (4 in 2D, 8 in 3D)
  const unsigned int dim(reference.GetNumRows());
  const unsigned int size(reference.GetNumCols());

  assert(other.GetNumRows() == dim);
  assert(other.GetNumCols() == size);
  assert(dim == domain->GetGrid()->GetDim());
  
  // defaults to make logging in the false case easier
  dimension = -1;
  positive = false;
  
  // we first test for x, then for y and then for z
  for(unsigned int d = 0; d < dim; ++d)
  {
    // assume 2D and d = 0 (=x), then all y-coordinates of 3 and 7 match 0
    bool good(true); // for above, only true for 0,3,7
    // probe for the other dimensions
    for(unsigned int p = 0; p < dim; ++p)
    {
      if(d == p) continue;
      // check all element nodes (size == 4 in 2D, size == 8 in 3D)
      for(unsigned int node = 0; node < size && good; ++node)
      {
        // for the neighbour in direction d, all coordinates of all points
        // in the other directions have to match
        // if we find a difference, we have no candidate!
        if(std::abs(reference[p][node] - other[p][node]) > 1e-6) good = false;
      }
    }

    if(!good) continue;

    // this element and the reference element differ only in direction d
    // check if we are comparing to self
    if(std::abs(reference[d][0] - other[d][0]) < 1e-6) return false;
    dimension = d;
    // check for positive or negative
    positive = (other[d][0] > reference[d][0]);
    return true;
  }
  // nothing found
  return false;
}


int VicinityElement::GetNumberOfEntries() const
{
  int count(0);
  for(unsigned int i = 0; i < design.GetSize(); ++i)
    if(design[i] != NULL) ++count;

  return count;
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

ResultDescription::ResultDescription()
{
  // does not set all!!
  access = DesignElement::PLAIN;
  value  = DesignElement::DESIGN;
  design = DesignElement::DEFAULT;
}

ResultDescription::ResultDescription(ParamNode* pn)
{
  solutionType = SolutionTypeEnum.Parse(pn->Get("id"));

  design = DesignElement::DEFAULT;
  if(pn->Has("design"))
    design = DesignElement::type.Parse(pn->Get("design"));

  access = DesignElement::access.Parse(pn->Get("access"));

  value = DesignElement::valueSpecifier.Parse(pn->Get("value"));

  detail = DesignElement::detail.Parse(pn->Get("detail"));
}
