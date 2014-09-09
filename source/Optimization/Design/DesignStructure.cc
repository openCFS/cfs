#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <limits>
#include <map>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/Logging/log.hpp"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Domain/BCs.hh"
#include "Domain/Domain.hh"
#include "Domain/ElemMapping/Elem.hh"
#include "Domain/ElemMapping/EntityLists.hh"
#include "Domain/Mesh/Grid.hh"
#include "FeBasis/BaseFE.hh"
#include "General/defs.hh"
#include "General/Exception.hh"
#include "MatVec/Matrix.hh"
#include "Optimization/Condition.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Design/DesignStructure.hh"
#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/Function.hh"
#include "Optimization/Objective.hh"
#include "Optimization/Optimization.hh"
#include "PDE/SinglePDE.hh"
#include "Utils/Timer.hh"
#include "Utils/tools.hh"

using std::string;
using std::map;
using namespace CoupledField;

DECLARE_LOG(ds)
DEFINE_LOG(ds, "designStructure")

Enum<DesignStructure::FilterSpace> DesignStructure::filterSpace;

DesignStructure::DesignStructure(DesignSpace* space, StdVector<RegionIdType>& regions)
{
  this->space = space;
  this->regions = regions;
  this->em = NULL;
  this->grid = domain->GetGrid();
  Constructor();
}

DesignStructure::DesignStructure(ErsatzMaterial* em)
{
  this->space = em->GetDesign();
  this->regions = em->GetDesign()->GetRegionIds();
  this->em = em;
  this->grid = domain->GetGrid();
  Constructor();
}

void DesignStructure::Constructor()
{
  initialized_ = false;

  this->dim  = grid->GetDim();

  periodic = false;
  if(em != NULL)
    for(map<Optimization::Application, SinglePDE*>::iterator it = em->pdes.begin(); it != em->pdes.end(); ++it)
      if(it->second->HasPeriodicBC()) periodic = true;

  filter_space_ = NO_FILTER;
  filter_ = Filter(); // defaults

  value  = -1.0;

}


void DesignStructure::Initialize()
{
  // save to be called multiple times. Has all neighbors and the number of common nodes
  grid->FindElementNeighorhood();

  // we will need the barycenters in FindNeibhborhood()
  for(unsigned int i = 0; i < regions.GetSize(); i++)
    grid->SetElementBarycenters(regions[i], false); // no updated coordinates
  // Handle also the off-design barycenters
  if(space->DoNonDesignVicinity())
    for(unsigned int i = 0; i < space->GetPseudoDesignRegions().GetSize(); i++)
      grid->SetElementBarycenters(space->GetPseudoDesignRegions()[i][0].elem->regionId, false);

  // can we assume all elements within all regions to be similar?
  regular = grid->IsRegionRegular(regions);

  if(periodic)
  {
    SetPeriodicConstraintMapping();
    SetNodeElemMapping();
    grid->CalcVolumeSpannedByNamedNodes(&dimension);
  }

  initialized_ = true;
}

void DesignStructure::SetFilters(PtrParamNode pn, PtrParamNode info, StdVector<DesignElement>* data_ptr)
{
  if(!initialized_)
    Initialize();

  StdVector<DesignElement>& data = data_ptr != NULL ? *data_ptr : space->data;

  filter_space_ = filterSpace.Parse(pn->Get("neighborhood")->As<string>());
  contribution_ = pn->Get("contribution")->As<string>() == "linear" ? LINEAR : CONSTANT;
  value  = pn->Get("value")->As<double>();
  if(pn->Has("design"))
    design = DesignElement::type.Parse(pn->Get("design")->As<std::string>());
  else
    design = DesignElement::ALL_DESIGNS;

  filter_.type_ = Filter::type.Parse(pn->Get("type")->As<std::string>());

  if(value <= 0.0)
    filter_.type_ = Filter::NO_FILTERING;

  if(filter_.type_ == Filter::SENSITIVITY && pn->Has("sensitivity"))
    filter_.sensitivity_ = Filter::sensitivity.Parse(pn->Get("sensitivity/type")->As<std::string>());

  if(filter_.type_ == Filter::DENSITY && pn->Has("density"))
    filter_.density_ = Filter::density.Parse(pn->Get("density/type")->As<string>());

  if(filter_.density_ != Filter::STANDARD)
  {
    if(!pn->Has("density/beta"))
      throw Exception("Attribute 'beta' required for '" + Filter::density.ToString(filter_.density_) + "' density filtering");
    filter_.SetBeta(pn->Get("density/beta")->As<double>(), space); // all relevant parameters set!

    if(pn->Has("density/force_lower_bound"))
      filter_.SetLowerBound(pn->Get("density/force_lower_bound")->As<double>());
  }

  if(filter_.density_ == Filter::TANH)
  {
    if(!pn->Has("density/eta"))
      throw Exception("Attribute 'eta' required for 'tanh' density filtering");
    filter_.eta = pn->Get("density/eta")->As<double>();
  }

  PtrParamNode in = info->Get(ParamNode::HEADER)->Get("regularization/filter", ParamNode::APPEND);

  in->Get("target")->SetValue(Filter::type.ToString(filter_.type_));

  // do we have to do something?
  if(filter_.type_ == Filter::NO_FILTERING)
    return;

  in->Get("type")->SetValue(filterSpace.ToString(filter_space_));

  in->Get("value")->SetValue(value);
  in->Get("contribution")->SetValue(contribution_ == LINEAR ? "linear" : "constant");

  if(filter_.type_ == Filter::SENSITIVITY)
    in->Get("sensitivity")->SetValue(Filter::sensitivity.ToString(filter_.sensitivity_));

  if(filter_.type_ == Filter::DENSITY)
  {
    in->Get("density")->SetValue(Filter::density.ToString(filter_.density_));
    if(filter_.density_ != Filter::STANDARD)
    {
      in->Get("beta")->SetValue(filter_.GetBeta());
      if(em != NULL && em->constraints.Has(Function::VOLUME) && em->constraints.Get(Function::VOLUME)->IsLinear())
        in->Get(ParamNode::WARNING)->SetValue("'volume' constraint shall be non-linear due to non-linear filter");
      if(filter_.density_ == Filter::HEAVISIDE)
        in->Get("heaviside_correction")->SetValue(filter_.heaviside_corr);
      if(pn->Has("density/force_lower_bound"))
        in->Get("force_lower_bound")->SetValue(filter_.GetLowerBound(NULL));
    }
  }

  in->Get("periodicBCs")->SetValue(periodic);

  // print about the function filtering
  for(unsigned int i = 0; em != NULL && i < em->objectives.data.GetSize(); i++)
  {
    PtrParamNode in_ = in->Get("functions")->Get("objective", ParamNode::APPEND);
    Objective* f = em->objectives.data[i];
    in_->Get("name")->SetValue(f->GetName());
    in_->Get("filtered")->SetValue(filter_.type_ == Filter::DENSITY ? f->ForDensityFiltering() : f->ForSensitivityFiltering());
  }
  for(unsigned int i = 0; em != NULL && i < em->constraints.active.GetSize(); i++)
  {
    PtrParamNode in_ = in->Get("functions")->Get("constraint", ParamNode::APPEND);
    Condition* g = em->constraints.active[i];
    in_->Get("name")->SetValue(g->ToString());
    in_->Get("filtered")->SetValue(filter_.type_ == Filter::DENSITY ? g->ForDensityFiltering() : g->ForSensitivityFiltering());
  }

  // the initialization was separated!
  boost::shared_ptr<Timer> timer(new Timer()); 
  in->Get("timer")->SetValue(timer);
  timer->Start();

  double avg_radius = 0;
  double avg_neighbours = 0;

  // find simp neighbors for all our elements
  double radius = -1.0; // for each element, set only once for regular.
  StdVector<SIMPElement::NeighbourElement> neighbors; // will become element neighborhood

  // for unstructured neighborhood search
  StdVector<unsigned int> too_far;   // element numbers too far away

  // our reference element dimensions for FindRegularNeighborhood()
  StdVector<double> edges;
  if(regular)
  {
    Matrix<double> coords; // temporary
    Elem* elem = data[0].elem;
    domain->GetGrid()->GetElemNodesCoord(coords, elem->connect);
    elem->ptElem->GetEdgeLength(coords, edges);

    // also initialize the vicinity elements!
    VicinityElement::Init(space, this);
  }

  unsigned int start = design == DesignElement::ALL_DESIGNS ? 0 : space->FindDesign(design)*space->GetNumberOfElements();
  unsigned int end = design == DesignElement::ALL_DESIGNS ? data.GetSize() : start + space->GetNumberOfElements();
  for(unsigned int e = start; e < end; e++)
  {
    DesignElement* de = &data[e];

    de->simp->filter = filter_;

    // independent of the filter type, radius determines the neighborhood
    // via barycenter distance.
    if(!regular || e == start)  // save calling if possible
      radius = FindFilterRadius(filter_space_, de, value);

    // set the filter neighborhood which is determined by radius
    // recursively via element neighbors.
    neighbors.Resize(0);
    too_far.Resize(0);

    LOG_DBG2(ds) << "SF: call FN for " << de->elem->ToString();
    if(regular)
      FindRegularNeighborhood(de, radius, edges, neighbors);
    else
      FindUnstructuredNeighborhood(de, radius, *(de->elem->neighborhood), neighbors, too_far); // works recursive
    // save neighborhood by copy constructor
    de->simp->neighborhood = neighbors;
    // set own weight
    assert(contribution_ == LINEAR || contribution_ == CONSTANT);
    de->simp->weight = (contribution_ == CONSTANT ? 1.0 : radius);

    // this is actually the re-implementation of a bug as it appeared to be not bad :)
    if(de->simp->filter.sensitivity_ == Filter::SHARP_SIGMUND || de->simp->filter.sensitivity_ == Filter::SHARP_PLAIN)
    {
      // normalize with a 'bug'
      double weight_sum = de->simp->CalcWeightSum(false) + 1.0;
      // assume 1.0 for this weight -> in the end it might be smaller! but in DesignElement::GetFilteredValue() we cheat 1.0 again
      de->simp->weight = 1.0 / weight_sum;
      for(unsigned int j = 0, n = de->simp->neighborhood.GetSize(); j < n; j++)
        de->simp->neighborhood[j].weight /= weight_sum;
    }

    avg_radius += radius;
    avg_neighbours += de->simp->neighborhood.GetSize();
    LOG_DBG2(ds) << "SF: final " << de->simp->ToString(0);
  }
  double normalized_avg_radius = avg_radius / (data.GetSize()/space->design.GetSize());
  double normalized_avg_neighbours = avg_neighbours / (data.GetSize()/space->design.GetSize());
  in->Get("avg_radius")->SetValue(normalized_avg_radius);
  in->Get("avg_neighbors")->SetValue(normalized_avg_neighbours);
//  in->Get("avg_radius")->SetValue(avg_radius / (end-start+1));
//  in->Get("avg_neighbors")->SetValue(avg_neighbours / (end-start+1));

  timer->Stop();
  
//  std::cout << "Filter: " << "avg radius: " << (avg_radius / (end-start+1))
//            << " avg neighbourhood: " << (avg_neighbours / (end-start+1)) << std::endl;

  std::cout << "Filter: " << "avg radius: " << normalized_avg_radius
            << " avg neighbourhood: " << normalized_avg_neighbours << std::endl;
}

void DesignStructure::FindRegularNeighborhood(DesignElement* base, double radius, const StdVector<double>& edges, StdVector<SIMPElement::NeighbourElement>& neighbors)
{
  assert(regular);
  // from the radius define a square/cube and check for every element. The corners are sorted out by distance

  // the legacy SHARP_PLAIN and SHARP_SIGMUND had the bug, that the weight was not
  // radius - distance but value - distance. To keep the legacy results we reproduce
  // the bug. Note, that another bug in SHARP_PLAIN and SHARP_SIGMUND is in normalization
  // and for SHARP_SIGMUND also in the filtering itself! :(
  double val_rad = filter_.sensitivity_ == Filter::SHARP_PLAIN || filter_.sensitivity_ == Filter::SHARP_SIGMUND ? value : radius;

  int x = static_cast<int>(ceil(radius / edges[0]));
  int y = static_cast<int>(ceil(radius / edges[1]));
  int z = dim < 3 ? 0 : static_cast<int>(ceil(radius / edges[2]));

  for(int i = -x; i <= x; i++)
  {
    for(int j = -y; j <= y; j++)
    {
      for(int k = -z; k <= z; k++) // ensure to enter one time in 2D
      {
        if(i == 0 && j == 0 && k == 0) // don't search for ourself!
          continue;

        DesignElement* other = GetNeighborElement(base, i, j, k);
        if(other != NULL)
        {
          // check the element
          double distance = RelaxedDistance(base->elem, other->elem);
          if(distance <= radius)
          {
            // value is here a double radius
            // this is the implementation from Bendsoe/ Sigmund
            SIMPElement::NeighbourElement ne;

            // map from element number to design
            ne.neighbour = other;

            // linear or constant weighting. will be normalized in the calling method!
            assert(contribution_ == LINEAR || contribution_ == CONSTANT);
            ne.weight = contribution_ == LINEAR ? val_rad  - distance : 1;
            ne.distance  = distance;
            neighbors.Push_back(ne); // cheap
          }
        }
      }
    }
  }
}

DesignElement* DesignStructure::GetNeighborElement(DesignElement* base, int i_steps, int j_steps, int k_steps)
{
  DesignElement* other = NULL;
  other = GetNeighborElement(base, abs(i_steps), i_steps < 0 ? VicinityElement::X_N : VicinityElement::X_P);
  if(other == NULL) return NULL;
  other = GetNeighborElement(other, abs(j_steps), j_steps < 0 ? VicinityElement::Y_N : VicinityElement::Y_P);
  if(other == NULL) return NULL;
  other = GetNeighborElement(other, abs(k_steps), k_steps < 0 ? VicinityElement::Z_N : VicinityElement::Z_P);
  return other;
}

DesignElement* DesignStructure::GetNeighborElement(DesignElement* base, unsigned int steps, VicinityElement::Neighbour dir)
{
  DesignElement* other = base;

  for(unsigned int i = 0; i < steps; i++)
  {
    // the periodic bcs are already done by vicinity element!!
    if(other->vicinity->HasNeighbor(dir)) {
      other = other->vicinity->GetNeighbour(dir);
    }
    else {
      assert(!periodic);
      return NULL;
    }
  }
  return other;
}


void DesignStructure::FindUnstructuredNeighborhood(DesignElement* base, double radius,
                                      StdVector<std::pair<Elem*, int> >& initial,
                                      StdVector<SIMPElement::NeighbourElement>& neighbors,
                                      StdVector<unsigned int>& too_far)
{
  LOG_DBG2(ds) << "FN: base= " << base->elem->elemNum << " initial=" << initial.ToString() << " n=" << ToString(neighbors) << " tf=" << too_far.ToString() << " ext=" << space->DoNonDesignVicinity();

  // the legacy SHARP_PLAIN and SHARP_SIGMUND had the bug, that the weight was not
  // radius - distance but value - distance. To keep the legacy results we reproduce
  // the bug. Note, that another bug in SHARP_PLAIN and SHARP_SIGMUND is in normalization
  // and for SHARP_SIGMUND also in the filtering itself! :(
  double val_rad = filter_.sensitivity_ == Filter::SHARP_PLAIN || filter_.sensitivity_ == Filter::SHARP_SIGMUND ? value : radius;

  assert(!periodic); // only regular may be periodic!!

    // the idea is as follows:
  // * We assume non regular grid.
  // * For an element t we check for all neighbors the distance to center
  // * If a neighbor is close enough we check also the neighbors recursively
  // * check means only, that the neighbors of check are checked!
  // * Hence buddies might grow (appending only) while traversing
  for(unsigned int e = 0, en = initial.GetSize(); e < en; e++)
  {
    // we ignore the grade of neighborhood (the int in the pair)
    const Elem* test_elem = initial[e].first;
    unsigned int test = test_elem->elemNum;

    if(test == base->elem->elemNum) continue; // we're not a neighbor of ourself

    // are we already a neighbor
    bool already = false;
    for(unsigned int n = 0; !already && n < neighbors.GetSize(); n++)
      if(neighbors[n].neighbour->elem->elemNum == test) already = true; // continue e loop!
    if(already) continue;

    // has it already been found that we are too far?
    if(too_far.Contains(test)) continue;

    // check the element if it is in the (possibly virtual) design space. If so we handle it as too far. May be NULL!
    DesignElement* test_de = space->Find(test, base->GetType(), false, space->DoNonDesignVicinity()); // silent

    // no need (and not possible!) to evaluate the distance for non-design elements
    double distance = test_de != NULL ? RelaxedDistance(base->elem, test_elem) : std::numeric_limits<double>::max();

    if(distance > radius || test_de == NULL)
    {
      too_far.Push_back(test);
    }
    else
    {
      // value is here a double radius
      // this is the implementation from Bendsoe/ Sigmund
      SIMPElement::NeighbourElement ne;

      // map from element number to design
      ne.neighbour = test_de;
      assert(ne.neighbour->elem->elemNum == test);

      // linear or constant weighting. will be normalized in the calling method!
      assert(contribution_ == LINEAR || contribution_ == CONSTANT);
      ne.weight = contribution_ == LINEAR ? val_rad  - distance : 1;
      ne.distance  = distance;
      neighbors.Push_back(ne); // cheap

      // now do the recursive call!!
      // test is in neighbors or too_far, hence the recursive call does't bounce back
      FindUnstructuredNeighborhood(base, radius, *test_elem->neighborhood, neighbors, too_far);
    }
  }
}

double DesignStructure::RelaxedDistance(const Elem* base, const Elem* test) const
{
  // default case
  const Point& bb = base->barycenter;
  const Point& tb = test->barycenter;

  assert(!(tb[0] == 0.0 && tb[1] == 0.0 && tb[2] == 0.0) && (test->ExpensiveCalcBarycenter()[0] != 0.0 ||  test->ExpensiveCalcBarycenter()[1] != 0.0));

  double dist = bb.Dist(tb);

  if(!periodic)
  {
    LOG_DBG3(ds) << "RD: dist " << base->elemNum << " <-> " << test->elemNum << " = " << dist << " : " << bb.ToString() << " <-> " << tb.ToString();
    return dist;
  }

  double preSqrt = dist * dist;

  double dx, dy, dz;

  // the idea is, that we map other via the periodic boundary condition.
  // we make assumptions for all possible cases and choose the shortest
  for(int z = -1; z <= 1; z += 1) // loop z only in 3D
  {
    for(int y = -1; y <= 1; y += 1)
    {
      for(int x = -1; x <= 1; x += 1)
      {
        dx = bb[0] - (tb[0] +  (double) x * dimension[0]);

        dy = bb[1] - (tb[1] +  (double) y * dimension[1]);

        dz = dim == 3 ? bb[2] - (tb[2] +  (double) z * dimension[2]) : 0.0;

        preSqrt = std::min(preSqrt, dx * dx + dy * dy + dz * dz);
      }
    }
  }

  LOG_DBG3(ds) << "RD: base=" << base->elemNum << " " << base->barycenter.ToString()
               << " test=" << test->elemNum << " " << test->barycenter.ToString()
               << " direct=" << dist << " relaxed=" << std::sqrt(preSqrt);

  return std::sqrt(preSqrt);
}

/** The is not performance tuned as for almost cases we have regular grids and then this method is
 * only called once. In the other cases - life with it */
double DesignStructure::FindFilterRadius(FilterSpace space, DesignElement* de, double value)
{
  Matrix<double>  coords;
  domain->GetGrid()->GetElemNodesCoord(coords, de->elem->connect, false );

  double radius, tmp;

  switch(space)
  {
    case RADIUS:
      radius = value;
      break;

    case VOLUME_RADIUS:
      // TODO really check for axis symmetry off
      tmp = de->elem->ptElem->CalcVolume(coords, false);
      // The radius is <value> times square/cube edge length where the
      // square/cube has the volume of the element
      radius = value * std::pow(tmp, 1.0/ (double) domain->GetGrid()->GetDim());
      LOG_DBG3(ds) << "FFR: de=" << de->ToString() << " from volume " << tmp << " to radius " << radius;
      break;

    case MAX_EDGE:
      de->elem->ptElem->GetMaxMinEdgeLength(coords, radius, tmp);
      radius = value * radius;
      LOG_DBG3(ds) << "FFR: de=" << de->ToString() << " edge max=" << radius << " min=" << tmp << " to radius " << radius;
      break;

    default:
      assert(false); // NO_FILTER
  }

  return radius;
}

void DesignStructure::SetPeriodicConstraintMapping()
{
  assert(em != NULL); // is not called otherwise!
   constraintMapping.Resize(grid->GetNumNodes() + 1); // 1-based

   ConstraintList glist = em->pde->GetConstraints();

   assert(em->pdes.size() == 1);
   assert(glist.GetSize() > 0);

   LOG_DBG(ds) << "SPCM: constraint list = " << glist.GetSize();

   StdVector<unsigned int> nlist;

   for(unsigned int i = 0, n = glist.GetSize(); i < n; i++)
   {
     shared_ptr<Constraint> g = glist[i];
     shared_ptr<EntityList> mlist = g->masterEntities;
     //EntityList* slist = g->slaveEntities.get();

     assert(mlist->GetType() == EntityList::NODE_LIST);
     assert(g->slaveEntities->GetType() == EntityList::NODE_LIST);

     // we assume two entries for the master and the same entries for the slave.
     assert(mlist->GetSize() == 2);
     assert(mlist->GetSize() == g->slaveEntities->GetSize());

     EntityIterator mit = mlist->GetIterator();
     nlist.Resize(0);

     for(mit.Begin(); !mit.IsEnd(); mit++)
       nlist.Push_back(mit.GetNode());

     assert(nlist.GetSize() == 2);

     // map against each other
     if(!constraintMapping[nlist[0]].Contains(nlist[1]))
       constraintMapping[nlist[0]].Push_back(nlist[1]);

     if(!constraintMapping[nlist[1]].Contains(nlist[0]))
       constraintMapping[nlist[1]].Push_back(nlist[0]);

     // 0 identifies not set
     assert(nlist[0] != 0 && nlist[1] != 0);

     LOG_DBG3(ds) << "stored pairs cn[" << nlist[0] << "]=" << constraintMapping[nlist[0]].ToString()
                  << " cn[" << nlist[1] << "]=" << constraintMapping[nlist[1]].ToString()
                  << " current pair: master=" << mlist->GetName() << " slave=" << g->slaveEntities->GetName();
   }

   // now we need a post processing as all connectesd 3D coner nodes can only be found by following
   // the chain. This is because the current constraint implemenation does not resolve multiple constraints
   // in the exact sense (but the computions are correct!)
   for(unsigned int n = 0, nn = constraintMapping.GetSize(); n < nn; n++)
   {
     StdVector<unsigned int>& cm = constraintMapping[n];
     if(cm.IsEmpty()) continue;
     // cm is extended by the recursive call, therefore it is important to start recursion only
     // for the pre existing entries.
     for(unsigned int o = 0, on = cm.GetSize(); o < on; o++)
       RecursiveCompletePeriodicity(n, cm, cm[o]);
     LOG_DBG3(ds) << "final cm[" << n << "]=" << cm.ToString();
   }

}

void DesignStructure::RecursiveCompletePeriodicity(unsigned int master, StdVector<unsigned int>& list, unsigned int test)
{
  StdVector<unsigned int>& other = constraintMapping[test];

  for(unsigned int o = 0; o < other.GetSize(); o++)
  {
    // anything new?
    if(other[o] == master || list.Contains(other[o]))
      continue;

    // yes, something is found!
    list.Push_back(other[o]);
    // now check if there are more references to come
    RecursiveCompletePeriodicity(master, list, other[o]);
  }
}

void DesignStructure::SetNodeElemMapping()
{
  assert(periodic);

  nodeToElem.Resize(grid->GetNumNodes() + 1,0); // 1-based

  StdVector<DesignElement>& data = space->data;
  // traverse all elements
  for(unsigned int e = 0, en = data.GetSize(); e < en; e++)
  {
    // we process all elements, as I don't see a cheap way to check
    // if the node has a periodic boundary node
    Elem* elem = data[e].elem;
    for(unsigned int n = 0, nn = elem->connect.GetSize(); n < nn; n++)
    {
      nodeToElem[elem->connect[n]] = elem; // we are only interested in one of the elements.
      LOG_DBG3(ds) << "SNEM: add elem " << elem->elemNum << " for node " << elem->connect[n];
    }
  }
}

bool DesignStructure::ExtendPeriodicNeighborhood(Elem* elem, int common, StdVector<std::pair<Elem*, int> >& neighbors)
{
  if(!initialized_)
    Initialize();

  assert(periodic);

  neighbors.Resize(0);

  // collect the constraint nodes of the main element
  constraintNodes_.Resize(0);

  for(unsigned int n = 0; n < elem->connect.GetSize(); n++)
  {
    unsigned int test = elem->connect[n];
    StdVector<unsigned int>& others = constraintMapping[test];

    for(unsigned int o = 0; o < others.GetSize(); o++)
      if(!constraintNodes_.Contains(others[o]))
        constraintNodes_.Push_back(others[o]);
  }

  // traverse over the elements of the pairing constraint nodes
  // Add only those elements with a minimum of 'common' nodes
  for(unsigned int i = 0; i < constraintNodes_.GetSize(); i++)
  {
    Elem* other_elem = nodeToElem[constraintNodes_[i]];
    AppendNeighbors(other_elem, constraintNodes_, common, neighbors);
  }

  // add the original neighborhood if there is a periodic case
  if(neighbors.GetSize() > 0)
    AppendNeighbors(*elem->neighborhood, neighbors);

  LOG_DBG3(ds) << "EPN_C: elem=" << elem->elemNum << " en=" << neighbors.ToString();

  return neighbors.GetSize() > 0;
}


void DesignStructure::AppendNeighbors(const StdVector<std::pair<Elem*, int> >& source, StdVector<std::pair<Elem*, int> >& out)
{
  for(unsigned int s = 0, sn = source.GetSize(); s < sn; s++)
  {
    bool found = false;
    for(unsigned int o = 0, on = out.GetSize(); o < on && !found; o++)
    {
      if(source[s].first == out[o].first)
        found = true;
    }
    if(!found)
      out.Push_back(source[s]);
  }
}

void DesignStructure::AppendNeighbors(Elem* check,
                                      const StdVector<unsigned int>& constraints, int min_common,
                                      StdVector<std::pair<Elem*, int> >& out)
{
  if(check == NULL) return;
  
  StdVector<std::pair<Elem*, int> >& source = *(check->neighborhood);

  // check all elements
  for(int s = -1, sn = (int) source.GetSize(); s < sn; s++)
  {
    Elem* test = s == -1 ? check : source[s].first;

    int common = 0;
    // check first if we have enough common
    for(unsigned int n = 0; n < test->connect.GetSize(); n++)
      if(constraints.Contains(test->connect[n])) common++;

    // is the element of interest?
    if(common < min_common)
      continue;

    bool found = false;
    // check if the element is known
    for(unsigned int o = 0, on = out.GetSize(); o < on && !found; o++)
    {
      if(test == out[o].first)
        found = true;
    }
    if(!found)
    {
      // fake the common information
      out.Push_back(std::make_pair(test, common));
    }
  }
}

string DesignStructure::ToString(StdVector<SIMPElement::NeighbourElement>& data)
{
  std::stringstream out;
  for(unsigned int i = 0, ni = data.GetSize(); i < ni; i++)
    out << data[i].neighbour->elem->elemNum
        << (i < ni-1 ? ", " : "");
  return out.str();
}

