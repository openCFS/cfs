#include "Optimization/DesignFilter.hh"
#include "Optimization/DesignSpace.hh"
#include "Optimization/SIMP.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"
#include "PDE/SinglePDE.hh"
#include "Elements/basefe.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Utils/Timer.hh"
#include "DataInOut/Logging/cfslog.hh"

using std::string;
using namespace CoupledField;

DECLARE_LOG(df)
DEFINE_LOG(df, "designFilter")

DesignFilter::DesignFilter(SIMP* simp, PtrParamNode pn)
{
  info_ = info->Get("optimization")->Get(ParamNode::HEADER)->Get("regularization/filter", ParamNode::APPEND);
  timer_ = new Timer();
  info_->SetValue(timer_);
  timer_->Start();

  Grid* grid = domain->GetGrid();

  this->simp = simp;
  this->dim  = grid->GetDim();

  StdVector<RegionIdType>&  regions = simp->regionIds;
  StdVector<DesignElement>& data    = simp->GetDesign()->data;

  periodic = false;
  for(unsigned int i = 0; i < simp->pdes.GetSize(); i++)
    if(simp->pdes[i]->HasPeriodicBC()) periodic = true;
  info_->Get("periodicBCs")->SetValue(periodic);

  // save to be called multiple times. Has all neighbors and the number of common nodes
  grid->FindElementNeighorhood();

  // the element point to Elem* neighbors, but we need DesignElement* neighbors
  // fill partially for direct access by element number
  elemToDesign.Resize(grid->GetNumElems() + 1); // 1 based.
  for(unsigned int i = 0, in = data.GetSize(); i < in; i++)
    elemToDesign[data[i].elem->elemNum] = &data[i];

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

  filter = DesignElement::filter.Parse(pn->Get("type")->As<std::string>());
  value  = pn->Get("value")->As<Double>();

  info_->Get("type")->SetValue(DesignElement::filter.ToString(filter));
  info_->Get("value")->SetValue(value);

  info_->Get("enabled")->SetValue(value > 0.0);

  timer_->Stop(); // to be restarted by SetFilter
}

void DesignFilter::SetFilters()
{
  if(value <= 0.0) return; // most times, there will be no filter when no filter is wanted and not a filter with value 0

  timer_->Start();

  double avg_radius = 0;
  double avg_neighbours = 0;
  StdVector<DesignElement>& data = simp->GetDesign()->data;

  // find simp neighbors for all our elements
  double radius = -1.0; // for each element, set only once for regular.
  StdVector<SIMPElement::NeighbourElement> neighbors; // will become element neighborhood
  StdVector<unsigned int> too_far;   // element numbers too far away
  StdVector<std::pair<Elem*, int> > base_buddies; // starting neighbors extended by peridic stuff

  for(unsigned int e = 0; e < data.GetSize(); e++)
  {
    DesignElement* de = &data[e];

    de->simp->filter = true;

    // independent of the filter type, radius determines the neighborhood
    // via barycenter distance.
    if(!regular || e == 0)  // save calling if possible
      radius = FindFilterRadius(filter, de);

    // set the filter neighborhood which is determined by radius
    // recursively via element neighbors.
    neighbors.Resize(0);
    too_far.Resize(0);
    base_buddies.Resize(0);

    // do we use the extended neighborhood?
    bool extend = periodic && ExtendPeriodicNeighborhood(de->elem, base_buddies);

    StdVector<std::pair<Elem*, int> >* start = extend ? &base_buddies : de->elem->neighborhood;

    LOG_DBG2(df) << "SF: call FN for " << de->elem->ToString();
    FindNeighborhood(de, radius, *start, neighbors, too_far); // works recursive
    // save neighborhood
    de->simp->neighborhood = neighbors;

    // now normalize the weights. The weight of this element is by definition 1.0
    double sum = 1.0;
    for(unsigned int e = 0; e < de->simp->neighborhood.GetSize(); e++)
      sum += de->simp->neighborhood[e].weight;

    // now normalize all
    de->simp->weight /= sum;
    for(unsigned int e = 0; e < de->simp->neighborhood.GetSize(); e++)
      de->simp->neighborhood[e].weight /= sum;

    avg_radius += radius;
    avg_neighbours += de->simp->neighborhood.GetSize();
  }

  info_->Get("avg_radius")->SetValue(avg_radius / data.GetSize());
  info_->Get("avg_neighbors")->SetValue(avg_neighbours / data.GetSize());

  timer_->Stop();

  std::cout << "Filter: avg radius: " << (avg_radius / data.GetSize())
            << " avg neighbourhood: " << (avg_neighbours / data.GetSize()) << std::endl;
}

void DesignFilter::FindNeighborhood(DesignElement* base, double radius,
                                      StdVector<std::pair<Elem*, int> >& buddies,
                                      StdVector<SIMPElement::NeighbourElement>& neighbors,
                                      StdVector<unsigned int>& too_far)
{
  LOG_DBG2(df) << "FN: base= " << base->elem->elemNum << " buddies=" << buddies.ToString() << " n=" << ToString(neighbors) << " tf=" << too_far.ToString();

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

      ne.neighbour = elemToDesign[test]; // map from element number to design
      assert(ne.neighbour->elem->elemNum == test);

      ne.weight    = value - distance;
      ne.distance  = distance;
      neighbors.Push_back(ne); // cheap

      // now do the recursive call!!
      // test is in neighbors or too_far, hence the recursive call does't bounce back
      FindNeighborhood(base, radius, *test_elem->neighborhood, neighbors, too_far);
    }
  }
}

double DesignFilter::RelaxedDistance(const Elem* base, const Elem* test) const
{
  // default case
  const Point& bb = base->barycenter;
  const Point& tb = test->barycenter;

  double dist = bb.Dist(tb);

  if(!periodic) return dist;

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

  LOG_DBG3(df) << "RD: base=" << base->elemNum << " " << base->barycenter.ToString()
               << " test=" << test->elemNum << " " << test->barycenter.ToString()
               << " direct=" << dist << " relaxed=" << std::sqrt(preSqrt);

  return std::sqrt(preSqrt);
}

/** The is not performance tuned as for almost cases we have regular grids and then this method is
 * only called once. In the other cases - life with it */
double DesignFilter::FindFilterRadius(DesignElement::Filter filter, DesignElement* de) const
{
  Matrix<double>  coords;
  domain->GetGrid()->GetElemNodesCoord(coords, de->elem->connect, false );

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
      radius = value * std::pow(tmp, 1.0/ (double) domain->GetGrid()->GetDim());
      LOG_DBG3(df) << "FFR: de=" << de->ToString() << " from volume " << tmp << " to radius " << radius;
      break;

    case DesignElement::MAX_EDGE:
      de->elem->ptElem->GetMaxMinEdgeLength(coords, radius, tmp);
      radius = value * radius;
      LOG_DBG3(df) << "FFR: de=" << de->ToString() << " edge max=" << radius << " min=" << tmp << " to radius " << radius;
      break;
  }

  return radius;
}

void DesignFilter::SetPeriodicConstraintMapping()
{
   constraintMapping.Resize(domain->GetGrid()->GetNumNodes() + 1,0); // 1-based

   ConstraintList glist = simp->pde->GetConstraints();

   assert(simp->pdes.GetSize() == 1);
   assert(glist.GetSize() > 0);

   LOG_DBG(df) << "SPCM: constraint list = " << glist.GetSize();

   StdVector<unsigned int> nlist;

   for(unsigned int i = 0, n = glist.GetSize(); i < n; i++)
   {
     shared_ptr<Constraint> g = glist[i];
     EntityList* mlist = g->masterEntities.get();
     //EntityList* slist = g->slaveEntities.get();

     assert(mlist->GetType() == EntityList::NODE_LIST);
     assert(g->slaveEntities.get()->GetType() == EntityList::NODE_LIST);

     LOG_DBG3(df) << "master=" << mlist->GetName() << " slave=" << g->slaveEntities.get()->GetName();

     // we assume two entries for the master and the same entries for the slave.
     assert(mlist->GetSize() == 2);
     assert(mlist->GetSize() == g->slaveEntities.get()->GetSize());

     EntityIterator mit = mlist->GetIterator();

     nlist.Resize(0);

     for(mit.Begin(); !mit.IsEnd(); mit++)
       nlist.Push_back(mit.GetNode());

     assert(nlist.GetSize() == 2);

     // map against each other
     constraintMapping[nlist[0]] = nlist[1];
     constraintMapping[nlist[1]] = nlist[0];

     // 0 identifies not set
     assert(nlist[0] != 0 && nlist[1] != 0);

     LOG_DBG3(df) << " add pair " << nlist[0] << ", " << nlist[1];
   }
}

void DesignFilter::SetNodeElemMapping()
{
  assert(periodic);

  nodeToElem.Resize(domain->GetGrid()->GetNumNodes() + 1,0); // 1-based

  StdVector<DesignElement>& data = simp->GetDesign()->data;
  // traverse all elements
  for(unsigned int e = 0, en = data.GetSize(); e < en; e++)
  {
    // we process all elements, as I don't see a cheap way to check
    // if the node has a periodic boundary node
    Elem* elem = data[e].elem;
    for(unsigned int n = 0, nn = elem->connect.GetSize(); n < nn; n++)
    {
      nodeToElem[elem->connect[n]] = elem; // we are only interested in one of the elements.
      LOG_DBG3(df) << "SNEM: add elem " << elem->elemNum << " for node " << elem->connect[n];
    }
  }
}

bool DesignFilter::ExtendPeriodicNeighborhood(Elem* elem, StdVector<std::pair<Elem*, int> >& neighbors)
{
  assert(periodic);

  neighbors.Resize(0);

  // check if any of the nodes is a periodic boundary node at all
  for(unsigned int n = 0; n < elem->connect.GetSize(); n++)
  {
    unsigned int test = elem->connect[n];
    unsigned int other = constraintMapping[test];

    // 0 means not set
    if(other == 0) continue;

    // take one of the elemnts that share the node. Due to the
    // recursive FindNeighborhood() the other elements will be covered
    Elem* other_elem = nodeToElem[other];
    AppendNeighbors(*other_elem->neighborhood, neighbors);
  }

  // add the original neighborhood if there is a periodic case
  if(neighbors.GetSize() > 0)
    AppendNeighbors(*elem->neighborhood, neighbors);

  LOG_DBG3(df) << "EPN: elem=" << elem->elemNum << " en=" << neighbors.ToString();

  return neighbors.GetSize() > 0;
}


void DesignFilter::AppendNeighbors(const StdVector<std::pair<Elem*, int> >& source, StdVector<std::pair<Elem*, int> >& out)
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

string DesignFilter::ToString(StdVector<SIMPElement::NeighbourElement>& data)
{
  std::stringstream out;
  for(unsigned int i = 0, ni = data.GetSize(); i < ni; i++)
    out << data[i].neighbour->elem->elemNum
        << (i < ni-1 ? ", " : "");
  return out.str();
}

