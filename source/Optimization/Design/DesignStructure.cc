#include <assert.h>
#include <cmath>
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
#include "Domain/Mesh/GridCFS/GridCFS.hh"
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
#include <def_use_openmp.hh>

#ifdef USE_OPENMP
  #include <omp.h>
#endif

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
  this->regions = em->GetDesign()->GetRegionIds();
  this->em = em;

  Constructor();
}

void DesignStructure::Constructor()
{
  initialized_ = false;
  num_robust_  = 0;

  this->grid = domain->GetGrid();
  this->gridcfs = dynamic_cast<GridCFS*>(domain->GetGrid());
  assert(this->gridcfs != NULL);

  this->dim  = grid->GetDim();

  bool enable = domain->GetParamRoot()->Has("optimization/ersatzMaterial/filters/periodic") ?
      domain->GetParamRoot()->Get("optimization/ersatzMaterial/filters/periodic")->As<bool>() : true;
  periodic_ = enable & domain->HasPerdiodicBC();

  filter_space_ = NO_FILTER;

  value  = -1.0;
}


void DesignStructure::Initialize()
{
  // save to be called multiple times. Has all neighbors and the number of common nodes
  grid->FindElementNeighorhood();

  // we will need the barycenters in FindNeighborhood()
  for(unsigned int i = 0; i < regions.GetSize(); i++)
    grid->SetElementBarycenters(regions[i], false); // no updated coordinates
  // Handle also the off-design barycenters
  if(space->DoNonDesignVicinity())
    for(unsigned int i = 0; i < space->GetPseudoDesignRegions().GetSize(); i++)
      grid->SetElementBarycenters(space->GetPseudoDesignRegions()[i][0].elem->regionId, false);

  // can we assume all elements within all regions to be similar?
  regular = grid->IsRegionRegular(regions);

  if(periodic_)
  {
    SetPeriodicConstraintMapping();
    SetNodeElemMapping();
    Matrix<double>& bb = grid->CalcGridBoundingBox();
    for(unsigned int r = 0; r < bb.GetNumRows(); r++)
      dimension[r] = bb[r][1] - bb[r][0];
  }

  LOG_DBG2(ds) << "I: regular=" << regular << " periodic=" << periodic_ << " dimension=" << dimension.ToString();

  initialized_ = true;
}

void DesignStructure::SetFilter(PtrParamNode pn, PtrParamNode info)
{
  if(!initialized_)
    Initialize();

  StdVector<DesignElement>& data = space->data;

  if(pn->Has("design"))
    design = DesignElement::type.Parse(pn->Get("design")->As<std::string>());
  else
    design = DesignElement::ALL_DESIGNS;

  // This are the designs we deal with
  unsigned int start = space->FindDesign(design) * space->GetNumberOfElements(); // handles ALL_DESIGNS
  unsigned int end = design == DesignElement::ALL_DESIGNS ? data.GetSize() : start + space->GetNumberOfElements();

  LOG_DBG(ds) << "SF: design=" << DesignElement::type.ToString(design) << " start=" << start << " end=" << end << " total=" << space->data.GetSize() << " nel=" << space->GetNumberOfElements();

  // Here we store our reference filter or more, if we do robust.
  // If SetFilters() was already called for this design (element) we must be robust!
  DesignElement& ref_de = data[start];

  StdVector<Filter>& filter = ref_de.simp->filter;
  unsigned int rex = 0;
  if(pn->Has("robust_excitation"))
  {
    rex = pn->Get("robust_excitation")->As<unsigned int>();
    if(rex != filter.GetSize())
      EXCEPTION("'robust_excitation' in 'filter' is " << rex << " but " << filter.GetSize() << " is expected");
  }
  else
  {
    if(!filter.IsEmpty())
      throw Exception("expect 'robust_excitation' for 'filter' as more than one filter is defined for the design");
  }

  Filter ref;

  ref.SetType(Filter::type.Parse(pn->Get("type")->As<std::string>()));
  ref.region = space->GetRegionIds()[0];

  filter_space_ = filterSpace.Parse(pn->Get("neighborhood")->As<string>());
  contribution_ = pn->Get("contribution")->As<string>() == "linear" ? LINEAR : CONSTANT;
  value  = pn->Get("value")->As<double>();

  if(value <= 0.0)
    ref.SetType(Filter::NO_FILTERING);

  if(ref.GetType() == Filter::SENSITIVITY && pn->Has("sensitivity"))
    ref.sensitivity_ = Filter::sensitivity.Parse(pn->Get("sensitivity/type")->As<std::string>());

  if(ref.GetType() == Filter::DENSITY && pn->Has("density"))
    ref.density_ = Filter::density.Parse(pn->Get("density/type")->As<string>());

  if(ref.density_ != Filter::STANDARD)
  {
    if(!pn->Has("density/beta"))
      throw Exception("Attribute 'beta' required for '" + Filter::density.ToString(ref.density_) + "' density filtering");
    ref.beta = pn->Get("density/beta")->As<double>();

    if(pn->Has("density/force_lower_bound"))
      ref.SetLowerBound(pn->Get("density/force_lower_bound")->As<double>());

    assert(space->design.GetSize() == 1); // extend for multiple regions with lower bounds which are now stored in non_lin_*

    if(ref.density_ == Filter::TANH)
    {
      if(!pn->Has("density/eta"))
        throw Exception("Attribute 'eta' required for 'tanh' density filtering");
      ref.eta = pn->Get("density/eta")->As<double>();
    }

    // for any projection filter
    ref.SetNonLinCorrection(&data[space->FindDesign(design) * space->GetNumberOfElements()], 0); // further robust filter might already be set
  }


  info->Get(ParamNode::HEADER)->Get("filters/periodic")->SetValue(periodic_);
  PtrParamNode in = info->Get(ParamNode::HEADER)->Get("filters")->Get("filter", ParamNode::APPEND);

  // do we have to do something?
  if(ref.GetType() == Filter::NO_FILTERING)
    return;

  // the initialization was separated!
  shared_ptr<Timer> timer = in->Get("timer")->AsTimer();
  timer->Start();

  double avg_radius = 0;
  double avg_neighbours = 0;

  // avoids multiple filter radius warning
  bool done = true;

  // for unstructured neighborhood search
  StdVector<unsigned int> too_far;   // element numbers too far away

  // our reference element dimensions for FindRegularNeighborhood()
  StdVector<double> edges;
  if(regular)
  {
    domain->GetGrid()->GetElemShapeMap(data[0].elem, false)->GetEdgeLength(edges);
    // also initialize the vicinity elements!
    VicinityElement::Init(space, this);
  }
  assert(!(design == DesignElement::DEFAULT && space->design.GetSize() > 1));

  LOG_DBG(ds) << "SF: design=" << DesignElement::type.ToString(design) << " start=" << start << " end=" << end << " total=" << space->data.GetSize() << " nel=" << space->GetNumberOfElements() << " rex=" << rex;

  DesignElement::Type ref_design = data[start].GetType();

  // calculate radius for for first element
  // in case grid is regular, set only once and not in loop
  double radius = FindFilterRadius(filter_space_, &data[start], value);

  #pragma omp parallel shared(ref)
  {
    // don't do it in for-loop, thread local vector
    StdVector<Filter::NeighbourElement> neighbors;

    #pragma omp for schedule(dynamic) reduction(+:avg_radius,avg_neighbours) firstprivate(radius)
    for(unsigned int e = start; e < end; e++)
    {
      DesignElement* de = &data[e];

      // did we came across a new design or a new region? Then update ref
      if(de->elem->regionId != ref.region || de->GetType() != ref_design)
      {
        ref.region = de->elem->regionId;
        ref.SetNonLinCorrection(de,rex);
        ref_design = de->GetType();
      }
      de->simp->filter.Push_back(ref); // copy the reference data

      assert(de->simp->filter.GetSize() == rex + 1); // we always work on the last filter in the filter vector

      // independent of the filter type, radius determines the neighborhood
      // via barycenter distance.
      if(!regular)   // save calling if possible
        radius = FindFilterRadius(filter_space_, de, value);

      // set the filter neighborhood which is determined by radius
      // recursively via element neighbors.
      neighbors.Resize(0); // keeps capacity

      LOG_DBG2(ds) << "SF: call FN for " << de->elem->ToString();
      if(regular)
        FindRegularNeighborhood(de, radius, edges, neighbors);
      else
        FindUnstructuredNeighborhood(de, radius, neighbors);

      // set own weight
      assert(contribution_ == LINEAR || contribution_ == CONSTANT);
      de->simp->filter.Last().weight = (contribution_ == CONSTANT ? 1.0 : radius);

      // this is actually the re-implementation of a bug as it appeared to be not bad :)
      if(de->simp->filter.Last().sensitivity_ == Filter::SHARP_SIGMUND || de->simp->filter.Last().sensitivity_ == Filter::SHARP_PLAIN)
      {
        // normalize with a 'bug'
        double weight_sum = de->simp->filter.Last().CalcWeightSum(false) + 1.0;
        // assume 1.0 for this weight -> in the end it might be smaller! but in DesignElement::GetFilteredValue() we cheat 1.0 again
        de->simp->filter.Last().weight = 1.0 / weight_sum;
        for(unsigned int j = 0, n = neighbors.GetSize(); j < n; j++)
          neighbors[j].weight /= weight_sum;
      }

      // save neighborhood by copy constructor
      de->simp->filter.Last().neighborhood = neighbors;

      avg_radius += radius;
      avg_neighbours += neighbors.GetSize();
      #pragma single nowait
      if(done && neighbors.GetSize() > 1000) {
        in->SetWarning("Filter radius too large. Neighborhood is bigger than 1000!");
        done = false;
      }
      LOG_DBG2(ds) << "SF: final " << de->simp->ToString(0);
    } // end for loop
  } // end omp region

  WriteFilterInfo(pn, in, ref, avg_radius, avg_neighbours, rex == 0); // goes into the appended filters/filter

  timer->Stop();
}


void DesignStructure::WriteFilterInfo(PtrParamNode pn, PtrParamNode in, const Filter& ref, double avg_radius, double avg_neighbours, bool first)
{
  in->Get("type")->SetValue(filterSpace.ToString(filter_space_));

  in->Get("value")->SetValue(value);
  in->Get("contribution")->SetValue(contribution_ == LINEAR ? "linear" : "constant");

  if(ref.GetType() == Filter::SENSITIVITY)
    in->Get("sensitivity")->SetValue(Filter::sensitivity.ToString(ref.sensitivity_));

  if(first && ref.GetType() == Filter::DENSITY)  {
    // in->Get("density")->SetValue(Filter::density.ToString(ref.density_));
    if(ref.density_ != Filter::STANDARD && em != NULL && em->constraints.Has(Function::VOLUME) && em->constraints.Get(Function::VOLUME)->IsLinear())
      in->SetWarning("'volume' constraint shall be non-linear due to non-linear filter");
  }

  if(periodic_)
   in->Get("periodic")->SetValue(periodic_);

  in->Get("design")->SetValue(DesignElement::type.ToString(design));

  if(pn->Has("robust_excitation")) // here we trust that we did good checks in SetFilter()
    in->Get("robust_excitation")->SetValue(pn->Get("robust_excitation")->As<int>());

  if(ref.GetType() == Filter::DENSITY && ref.density_ != Filter::STANDARD)
  {
    PtrParamNode in_ = in->Get(ref.density.ToString(ref.density_));

    in_->Get("beta")->SetValue(ref.beta);

    if(ref.density_ == Filter::TANH)
      in_->Get("eta")->SetValue(ref.eta);

    in_->Get("scaling")->SetValue(ref.non_lin_scale);
    in_->Get("offset")->SetValue(ref.non_lin_offset);

    if(pn->Has("density/force_lower_bound"))
      in_->Get("force_lower_bound")->SetValue(ref.GetLowerBound(NULL));
  }

  double normalized_avg_radius = avg_radius / (space->data.GetSize()/space->design.GetSize());
  double normalized_avg_neighbours = avg_neighbours / (space->data.GetSize()/space->design.GetSize());

  in->Get("avg_radius")->SetValue(normalized_avg_radius);
  in->Get("avg_neighbors")->SetValue(normalized_avg_neighbours);

  std::cout << "Filter " << DesignElement::type.ToString(design) << ": avg radius=" << normalized_avg_radius
            << " avg neighbourhood=" << normalized_avg_neighbours << std::endl;

}

void DesignStructure::FindRegularNeighborhood(DesignElement* base, double radius, const StdVector<double>& edges, StdVector<Filter::NeighbourElement>& neighbors)
{
  assert(regular);
  // from the radius define a square/cube and check for every element. The corners are sorted out by distance

  // the legacy SHARP_PLAIN and SHARP_SIGMUND had the bug, that the weight was not
  // radius - distance but value - distance. To keep the legacy results we reproduce
  // the bug. Note, that another bug in SHARP_PLAIN and SHARP_SIGMUND is in normalization
  // and for SHARP_SIGMUND also in the filtering itself! :(

  Filter& filter = base->simp->filter[0];
  double val_rad = filter.sensitivity_ == Filter::SHARP_PLAIN || filter.sensitivity_ == Filter::SHARP_SIGMUND ? value : radius;

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
            Filter::NeighbourElement ne;

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
  if(other == NULL)
    return NULL;
  other = GetNeighborElement(other, abs(j_steps), j_steps < 0 ? VicinityElement::Y_N : VicinityElement::Y_P);
  if(other == NULL)
    return NULL;
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
      // with periodic b.c. we shall always have a neighbor! - but only if the full domain is design domain
      assert(!periodic_ || domain->GetGrid()->GetNumRegions() > space->design.GetSize());
      return NULL;
    }
  }
  return other;
}


void DesignStructure::FindUnstructuredNeighborhood(DesignElement* base, double radius,
                                                   StdVector<Filter::NeighbourElement>& neighbors)
{
  // for 3D ist shall be significantly faster (O(N)) to make a discrete grid (e.g. of size radius) and sort
  // the elements to this grid. Then we can linearly test all relevant grid cells for an element. It will help for
  // large meshes where the filter setup is otherwise hours and magnitudes slower than solving the system!

  // LOG_DBG2(ds) << "FN: base= " << base->elem->elemNum << " initial=" << ToString(initial) << " n=" << ToString(neighbors) << " tf=" << too_far.ToString() << " ext=" << space->DoNonDesignVicinity();

  // the legacy SHARP_PLAIN and SHARP_SIGMUND had the bug, that the weight was not
  // radius - distance but value - distance. To keep the legacy results we reproduce
  // the bug. Note, that another bug in SHARP_PLAIN and SHARP_SIGMUND is in normalization
  // and for SHARP_SIGMUND also in the filtering itself! :(
  Filter& filter = base->simp->filter[0];
  double val_rad = filter.sensitivity_ == Filter::SHARP_PLAIN || filter.sensitivity_ == Filter::SHARP_SIGMUND ? value : radius;

  assert(!periodic_); // only regular may be periodic!!

  // the idea is as follows:
  // * We assume non regular grid.
  // * For an element t we check for all neighbors the distance to center
  // * If a neighbor is close enough we check also the neighbors recursively
  // * check means only, that the neighbors of check are checked!
  // * Hence buddies might grow (appending only) while traversing
  std::set<unsigned int> checked;

  int dim = (int) grid->GetDim();

  // base is not a neighbor of itself and does not need to be checked
  checked.insert(base->elem->elemNum);

  // we start with the base neighborhood to be checked
  assert(base->elem->extended->neighborhood != NULL);
  std::set<unsigned int> to_check;

  const StdVector<std::pair<Elem*, int> >& initial = *(base->elem->extended->neighborhood);
  for(unsigned int j = 0; j < initial.GetSize(); j++) {
    // we reduce to surface neighbors. second is the number of sharing nodes
    // for 2d the surface is = 2 (edge) for 3D the surface is >= 3 (triangle)
    if(initial[j].second >= dim)
      to_check.insert(initial[j].first->elemNum);
    assert(initial[j].first->elemNum != base->elem->elemNum);
  }

  // The idea is as follows. We start with the Elem* neighbors of the base element. All these elements are candidate for the
  // filter (will be stored in neighbors) if they are not too far away. Then in principle we recursively check for all candidates
  // all Elem* neighbors up to all is either a neighbor or too_far. A direct recursive implementation of this is for large meshed
  // prohibitive expensive (> 33h for 3e6 element!). Revision 15237 of shared_opt still has this implementation
  // Therefore we have the life loop below which is more efficient than the recursive implementation
  while(to_check.size() > 0)
  {
    // note that his loop run my add several new items to to_check
    unsigned int test_num = *(to_check.begin());
    const Elem*  test     = gridcfs->GetElem(test_num); // shall be sped up!!
    // remove test, it will be added to checked in 1e-10 seconds
    to_check.erase(to_check.begin());

    // we can assume, that we were not checked before we were put in to_check
    assert(test_num != base->elem->elemNum);

    // test can already be considered as checked
    checked.insert(test_num);

    // check the element if it is in the (possibly virtual) design space. If so we handle it as too far. May be NULL!
    DesignElement* test_de = space->Find(test_num, base->GetType(), false, space->DoNonDesignVicinity()); // silent

    // no need (and not possible!) to evaluate the distance for non-design elements
    double distance = test_de != NULL ? RelaxedDistance(base->elem, test) : std::numeric_limits<double>::max();

    if(test_de != NULL && distance <= radius)
    {
      // value is here a double radius
      // this is the implementation from Bendsoe/ Sigmund
      Filter::NeighbourElement ne;

      // map from element number to design
      ne.neighbour = test_de;
      assert(ne.neighbour->elem->elemNum == test->elemNum);

      // linear or constant weighting. will be normalized in the calling method!
      assert(contribution_ == LINEAR || contribution_ == CONSTANT);
      ne.weight = contribution_ == LINEAR ? val_rad  - distance : 1;
      ne.distance  = distance;

      #ifndef NDEBUG
        for(unsigned int k=0; k < neighbors.GetSize(); k++)
          assert(neighbors[k].neighbour->elem->elemNum != test_num);
      #endif

      neighbors.Push_back(ne); // cheap

      // as test was a successful element we have to process the Elem* neighbors
      const StdVector<std::pair<Elem*, int> >& test_ne = *(test->extended->neighborhood);
      for(unsigned e = 0; e < test_ne.GetSize(); e++)
      {
        // see above!
        if(test_ne[e].second >= dim)
        {
          Elem* cand = test_ne[e].first;

          // only if the candidate is not already checked and also not already queued it will be processed
          // the sets are sorted and unique. insert() returns if insertion was necessary due to uniqueness or not
          assert(checked.find(base->elem->elemNum) != checked.end());
          if(checked.find(cand->elemNum) == checked.end())
          {
            assert(cand->elemNum != base->elem->elemNum);
            to_check.insert(cand->elemNum); // it it was already there it's no problem
          }
        }
      }
    }
  }
}

inline double DesignStructure::RelaxedDistance(const Elem* base, const Elem* test) const
{
  // default case
  const Point& bb = base->extended->barycenter;
  const Point& tb = test->extended->barycenter;

  assert(!(tb[0] == 0.0 && tb[1] == 0.0 && tb[2] == 0.0)/* && (test->ExpensiveCalcBarycenter()[0] != 0.0 ||  test->ExpensiveCalcBarycenter()[1] != 0.0)*/);

  double dist = bb.Dist(tb);

  if(!periodic_)
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

  LOG_DBG3(ds) << "RD: base=" << base->elemNum << " " << base->extended->barycenter.ToString()
               << " test=" << test->elemNum << " " << test->extended->barycenter.ToString()
               << " direct=" << dist << " relaxed=" << std::sqrt(preSqrt);

  return std::sqrt(preSqrt);
}

/** The is not performance tuned as for almost cases we have regular grids and then this method is
 * only called once. In the other cases - life with it */
double DesignStructure::FindFilterRadius(FilterSpace space, DesignElement* de, double value)
{
  Matrix<double>  coords;
  domain->GetGrid()->GetElemNodesCoord(coords, de->elem->connect, false );


  switch(space)
  {
    case RADIUS:
      return value;

    case VOLUME_RADIUS:
    {
      // TODO really check for axis symmetry off
      double tmp = domain->GetGrid()->GetElemShapeMap(de->elem, false)->CalcVolume();
      // The radius is <value> times square/cube edge length where the
      // square/cube has the volume of the element
      double radius = value * std::pow(tmp, 1.0/ (double) domain->GetGrid()->GetDim());
      LOG_DBG3(ds) << "FFR: de=" << de->ToString() << " from volume " << tmp << " to radius " << radius;
      return radius;
    }

    case MAX_EDGE:
    {
      double max, tmp;
      LagrangeElemShapeMap sm(domain->GetGrid());
      sm.SetElem(de->elem,false);
      sm.GetMaxMinEdgeLength(max, tmp);
      double radius = value * max;
      LOG_DBG3(ds) << "FFR: de=" << de->ToString() << " edge max=" << max << " min=" << tmp << " to radius " << radius;
      return radius;
    }

    default:
      assert(false); // NO_FILTER
      return -1.0;
  }
}

void DesignStructure::SetPeriodicConstraintMapping()
{
   assert(em != NULL); // is not called otherwise!
   constraintMapping.Resize(grid->GetNumNodes() + 1); // 1-based

   assert(Optimization::context->pde != NULL); // if(em->pde == NULL) em->SetPDEs(em->ParseSystem());

   ConstraintList glist = Optimization::context->pde->GetFeFunctions().begin()->second->GetConstraints();

   assert(Optimization::context->pdes.size() == 1);
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
  assert(periodic_);

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

  assert(periodic_);

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
    AppendNeighbors(*elem->extended->neighborhood, neighbors);

  LOG_DBG3(ds) << "EPN_C: elem=" << elem->elemNum << " en=" << ToString(neighbors);

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
  
  StdVector<std::pair<Elem*, int> >& source = *(check->extended->neighborhood);

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

std::string DesignStructure::ToString(const StdVector<Filter::NeighbourElement>& data)
{
  std::stringstream out;
  for(unsigned int i = 0, ni = data.GetSize(); i < ni; i++)
    out << data[i].neighbour->elem->elemNum << (i < ni-1 ? ", " : "");
  return out.str();
}

std::string DesignStructure::ToString(const StdVector<std::pair<Elem*, int> >& data)
{
  std::stringstream ss;
  ss << "(";
  for(unsigned int i = 0; i < data.GetSize(); i++)
    ss << i << ": " << data[i].first->elemNum << "," << data[i].second << (i < data.GetSize()-1 ? "; " :  ")");
  return ss.str();
}
