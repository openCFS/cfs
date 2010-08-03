#include "Optimization/Design/DesignStructure.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/SIMP.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"
#include "PDE/SinglePDE.hh"
#include "Elements/basefe.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Utils/Timer.hh"
#include "DataInOut/Logging/cfslog.hh"

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
  Constructor();
}

DesignStructure::DesignStructure(ErsatzMaterial* em)
{
  this->space = em->GetDesign();
  this->regions = em->regionIds;
  this->em = em;
  Constructor();
}

void DesignStructure::Constructor()
{
  initialized_ = false;

  this->dim  = domain->GetGrid()->GetDim();

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
  Grid* grid = domain->GetGrid();

  // save to be called multiple times. Has all neighbors and the number of common nodes
  grid->FindElementNeighorhood();

  // we will need the barycenters in FindNeibhborhood()
  for(unsigned int i = 0; i < regions.GetSize(); i++)
    grid->SetElementBarycenters(regions[i], false); // no updated coordinates

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

void DesignStructure::SetFilters(PtrParamNode pn, PtrParamNode info)
{
  if(!initialized_)
    Initialize();

  filter_space_ = filterSpace.Parse(pn->Get("neighborhood")->As<string>());
  contribution_ = pn->Get("contribution")->As<string>() == "linear" ? LINEAR : CONSTANT;
  value  = pn->Get("value")->As<double>();

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
    filter_.beta_ = pn->Get("density/beta")->As<double>();
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
      in->Get("beta")->SetValue(filter_.beta_);
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
  Timer* timer = new Timer(); 
  in->Get("timer")->SetValue(timer);
  timer->Start();

  double avg_radius = 0;
  double avg_neighbours = 0;
  StdVector<DesignElement>& data = space->data;

  // find simp neighbors for all our elements
  double radius = -1.0; // for each element, set only once for regular.
  StdVector<SIMPElement::NeighbourElement> neighbors; // will become element neighborhood
  StdVector<unsigned int> too_far;   // element numbers too far away
  StdVector<std::pair<Elem*, int> > base_buddies; // starting neighbors extended by periodic stuff

  for(unsigned int e = 0; e < data.GetSize(); e++)
  {
    DesignElement* de = &data[e];

    de->simp->filter = filter_;

    // independent of the filter type, radius determines the neighborhood
    // via barycenter distance.
    if(!regular || e == 0)  // save calling if possible
      radius = FindFilterRadius(filter_space_, de, value);

    // set the filter neighborhood which is determined by radius
    // recursively via element neighbors.
    neighbors.Resize(0);
    too_far.Resize(0);
    base_buddies.Resize(0);

    // do we use the extended neighborhood?
    bool extend = periodic && ExtendPeriodicNeighborhood(de->elem, base_buddies);

    StdVector<std::pair<Elem*, int> >* start = extend ? &base_buddies : de->elem->neighborhood;

    LOG_DBG2(ds) << "SF: call FN for " << de->elem->ToString();
    FindNeighborhood(de, radius, *start, neighbors, too_far); // works recursive
    // save neighborhood
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
  }

  in->Get("avg_radius")->SetValue(avg_radius / data.GetSize());
  in->Get("avg_neighbors")->SetValue(avg_neighbours / data.GetSize());

  timer->Stop();

  std::cout << "Filter: avg radius: " << (avg_radius / data.GetSize())
            << " avg neighbourhood: " << (avg_neighbours / data.GetSize()) << std::endl;
}

void DesignStructure::FindNeighborhood(DesignElement* base, double radius,
                                      StdVector<std::pair<Elem*, int> >& buddies,
                                      StdVector<SIMPElement::NeighbourElement>& neighbors,
                                      StdVector<unsigned int>& too_far)
{
  LOG_DBG2(ds) << "FN: base= " << base->elem->elemNum << " buddies=" << buddies.ToString() << " n=" << ToString(neighbors) << " tf=" << too_far.ToString();

  // the legacy SHARP_PLAIN and SHARP_SIGMUND had the bug, that the weight was not
  // radius - distance but value - distance. To keep the legacy results we reproduce
  // the bug. Note, that anoter bug in SHARP_PLAIN and SHARP_SIGMUND is in normalization
  // and for SHARP_SIGMUND also in the filtering itself! :(
  double val_rad = filter_.sensitivity_ == Filter::SHARP_PLAIN || filter_.sensitivity_ == Filter::SHARP_SIGMUND ? value : radius;

    // the idea is as follows:
  // * We assume non regular grid.
  // * For an element t we check for all neighbors the distance to center
  // * If a neighbor is close enough we check also the neighbors recursively
  // * check means only, that the neighbors of check are checked!
  for(unsigned int e = 0; e < buddies.GetSize(); e++)
  {
    // we ignore the grade of neighborhood (the int in the pair)
    Elem* test_elem = buddies[e].first;
    unsigned int test = test_elem->elemNum;

    if(test == base->elem->elemNum) continue; // we're not a neighbor of ourself

    // are we already a neighbor
    bool already = false;
    for(unsigned int n = 0; !already && n < neighbors.GetSize(); n++)
      if(neighbors[n].neighbour->elem->elemNum == test) already = true; // continue e loop!
    if(already) continue;

    // has it already been found that we are too far?
    if(too_far.Contains(test)) continue;

    // check the element
    double distance = RelaxedDistance(base->elem, test_elem);

    if(distance > radius)
    {
      too_far.Push_back(test);
    }
    else
    {
      // value is here a double radius
      // this is the implementation from Bendsoe/ Sigmund
      SIMPElement::NeighbourElement ne;

      // map from element number to design
      ne.neighbour = space->Find(test, base->GetType());
      assert(ne.neighbour->elem->elemNum == test);

      // linear or constant weighting. will be normalized in the calling method!
      assert(contribution_ == LINEAR || contribution_ == CONSTANT);
      ne.weight = contribution_ == LINEAR ? val_rad  - distance : 1;
      ne.distance  = distance;
      neighbors.Push_back(ne); // cheap

      // now do the recursive call!!
      // test is in neighbors or too_far, hence the recursive call does't bounce back
      FindNeighborhood(base, radius, *test_elem->neighborhood, neighbors, too_far);
    }
  }
}

double DesignStructure::RelaxedDistance(const Elem* base, const Elem* test) const
{
  // default case
  const Point& bb = base->barycenter;
  const Point& tb = test->barycenter;

  double dist = bb.Dist(tb);

  if(!periodic)
  {
    LOG_DBG3(ds) << "RD: dist " << base->elemNum << " <-> " << test->elemNum << " = " << dist;
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
   constraintMapping.Resize(domain->GetGrid()->GetNumNodes() + 1); // 1-based

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

  nodeToElem.Resize(domain->GetGrid()->GetNumNodes() + 1,0); // 1-based

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


bool DesignStructure::ExtendPeriodicNeighborhood(Elem* elem, StdVector<std::pair<Elem*, int> >& neighbors)
{
  assert(periodic);

  neighbors.Resize(0);

  // check if any of the nodes is a periodic boundary node at all
  for(unsigned int n = 0; n < elem->connect.GetSize(); n++)
  {
    unsigned int test = elem->connect[n];
    StdVector<unsigned int>& others = constraintMapping[test];

    if(others.IsEmpty()) continue;

    // take one of the elements that share the node.
    for(unsigned int o = 0; o < others.GetSize(); o++)
    {
      Elem* other_elem = nodeToElem[others[o]];
      AppendNeighbors(*other_elem->neighborhood, neighbors);
    }
  }

  // add the original neighborhood if there is a periodic case
  if(neighbors.GetSize() > 0)
    AppendNeighbors(*elem->neighborhood, neighbors);

  LOG_DBG3(ds) << "EPN: elem=" << elem->elemNum << " en=" << neighbors.ToString();

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

