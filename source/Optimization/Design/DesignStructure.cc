#include <assert.h>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <limits>
#include <map>

#include "DataInOut/Logging/LogConfigurator.hh"
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
#include "Optimization/Tune.hh"
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

DEFINE_LOG(ds, "designStructure")

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
  timer = shared_ptr<Timer>(new Timer());

  this->grid = domain->GetGrid();
  this->gridcfs = dynamic_cast<GridCFS*>(domain->GetGrid());
  assert(this->gridcfs != NULL);

  this->dim  = grid->GetDim();

  bool enable = domain->GetParamRoot()->Has("optimization/ersatzMaterial/filters/periodic") ?
      domain->GetParamRoot()->Get("optimization/ersatzMaterial/filters/periodic")->As<bool>() : true;
  periodic_ = enable & domain->HasPerdiodicBC();

  assert(space);
  assert(space->filter.IsEmpty());
  // we only store pointers, so we can serve a lot!
  // space->region fails for constant_region test case
  space->filter.Reserve(domain->GetGrid()->regionData.GetSize() * space->design.GetSize() * 4); // hope 4 robust excitations are sufficient
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


void DesignStructure::SetFilter(PtrParamNode pn, bool skip_cfs_filtering)
{
  if(!initialized_)
    Initialize();

  StdVector<DesignElement>& data   = space->data;
  StdVector<GlobalFilter>&  filter = space->filter;

  // might be ALL_DESIGNS, in GlobalFilter, we set the real design
  DesignElement::Type design = pn->Has("design") ? DesignElement::type.Parse(pn->Get("design")->As<std::string>()) : DesignElement::ALL_DESIGNS;

  // This are the designs we deal with
  unsigned int start = space->FindDesign(design) * space->GetNumberOfElements(); // handles ALL_DESIGNS
  unsigned int end = design == DesignElement::ALL_DESIGNS ? data.GetSize() : start + space->GetNumberOfElements();

  LOG_DBG(ds) << "SF: design=" << DesignElement::type.ToString(design) << " start=" << start << " end=" << end << " total=" << space->data.GetSize() << " nel=" << space->GetNumberOfElements();

  assert(filter.GetCapacity() > filter.GetSize()); // e.g. in the robust case we might already have a filter

  filter.Push_back(Parse(pn, &data[start], skip_cfs_filtering));
  // while meeting other regions or designs we create new global filters based on this one.
  GlobalFilter* global = &filter.Last();
  LOG_DBG(ds) << "SF: initial global=" << global->ToString();
  global->PostCopy(true); // register tune

  // do we have to do something?
  // in case skip_cfs_filtering is true, we don't want cfs to apply filter but still assemble it for external (SGP) lib
  if(global->type == Filter::NO_FILTERING && !skip_cfs_filtering)
    return;

  // the initialization was separated!
  timer->Start();

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

  LOG_DBG(ds) << "SF: design=" << DesignElement::type.ToString(design) << " start=" << start << " end=" << end << " total=" << space->data.GetSize() << " nel=" << space->GetNumberOfElements() << " rex=" << global->robust;

  const DesignElement* initial_ref_de = &data[start]; // we later use a parallelized de
  global->design = initial_ref_de->GetType(); // no ALL_DESIGNS but be specific
  global->region = initial_ref_de->elem->regionId;
  global->SetNonLinCorrection(initial_ref_de);

  // calculate radius for for first element
  // in case grid is regular, set only once and not in loop
  double radius = FindFilterRadius(global->filterspace, initial_ref_de, global->value);

  int total_sum_neigbors = 0; // for the filter matrix

  // we call this function by filter element in xml but we might cross several regions or possibly designs
  // each region might have own bound, therefore each region needs a own GlobalFilter with own nonlinear filter scaling
  const StdVector<DesignSpace::DesignRegion*> units = space->GetRegions(design);
  for(unsigned int u = 0; u < units.GetSize(); u++)
  {
    const DesignSpace::DesignRegion* unit = units[u];
    assert(unit->base >= start && unit->base < end);

    double sum_radius = 0;
    double sum_neighbours = 0;

    if(u > 0) // we come along another GlobalFilter
    {
      assert(global->region != unit->regionId || global->design != unit->design); // there shall be a reason for change

      // append a new GlobalFilter which is a copy of the current one.
      if(filter.GetCapacity() <= filter.GetSize()) // resizing would destroy all pointer references
        throw Exception("DesignSpace::filter capacity too small: " + std::to_string(filter.GetCapacity()));
      filter.Push_back(*global); // copy constructor magic - but handle math parser!
      global = &(filter.Last());
      global->PostCopy(true); // register tune
      global->region = unit->regionId; // we need the region as one way to validate the change
      global->design = unit->design;
      global->SetNonLinCorrection(&(space->data[unit->base]));
    }
    // MSVC want's ints for parallel loops ?!
    int unit_start = (int) unit->base;
    int unit_end   = u < units.GetSize()-1 ? units[u+1]->base : end;

    LOG_DBG(ds) << "SF: unit=" << u << "/" << units.GetSize() << " us=" << unit_start << " ue=" << unit_end << " g=" << global->ToString();

    #pragma omp parallel num_threads(CFS_NUM_THREADS)
    {
      // don't do it in for-loop, thread local vector
      StdVector<Filter::NeighbourElement> neighbors;

      #pragma omp for reduction(+:sum_radius,sum_neighbours)
      for(int e = unit_start; e < unit_end; e++)
      {
        DesignElement* de = &data[e];

        assert(de->GetType() == global->design && de->elem->regionId == global->region);

        de->simp->filter.Push_back(Filter());
        Filter& elem_filter = de->simp->filter.Last();
        elem_filter.global = global;

        /* what does this assert do? deactivating for now */
        //assert(de->simp->filter.GetSize() == rex + 1); // we always work on the last filter in the filter vector
        // for unstructured grids only "radius" filter makes sense
        assert(regular || global->filterspace == Filter::RADIUS);

        // set the filter neighborhood which is determined by radius
        // recursively via element neighbors.
        neighbors.Resize(0); // keeps capacity

        // LOG_DBG2(ds) << "SF: call FN for " << de->elem->ToString();
        if(regular)
          FindRegularNeighborhood(global, de, radius, edges, neighbors);
        else
          FindUnstructuredNeighborhood(global, de, radius, neighbors);

        // set own weight
        assert(global->contribution == Filter::LINEAR || global->contribution == Filter::CONSTANT);
        elem_filter.weight = global->contribution == Filter::CONSTANT ? 1.0 : radius;

        // this is actually the re-implementation of a bug as it appeared to be not bad :)
        if(global->sensitivity == Filter::SHARP_SIGMUND || global->sensitivity == Filter::SHARP_PLAIN)
        {
          // normalize with a 'bug'
          double weight_sum = elem_filter.CalcWeightSum(false) + 1.0;
          // assume 1.0 for this weight -> in the end it might be smaller! but in DesignElement::GetFilteredValue() we cheat 1.0 again
          elem_filter.weight = 1.0 / weight_sum;
          for(unsigned int j = 0, n = neighbors.GetSize(); j < n; j++)
            neighbors[j].weight /= weight_sum;
        }

        // save neighborhood by copy constructor
        elem_filter.neighborhood = neighbors;

        sum_radius += radius;
        sum_neighbours += neighbors.GetSize();
        if(neighbors.GetSize() > 1000) {
          too_large_filter_++;
        }

        // have weight_sum ready for density filters
        elem_filter.weight_sum = elem_filter.CalcWeightSum(true);

        assert(elem_filter.global != NULL);

        LOG_DBG2(ds) << "SF: final " << de->simp->ToString(0);
      } // end for parallel element loop within unit
    } // end parallel bloc
    global->elements = unit_end - unit_start;
    assert(global->elements > 0);
    global->avg_radius = sum_radius / global->elements;
    global->avg_neigbor = sum_neighbours / global->elements;
    total_sum_neigbors += sum_neighbours;

  } // end loop units

  timer->Stop();

  if (space->is_matrix_filt)
  {
    if(filter.Last().density == Filter::MATERIAL)
      throw Exception("Disable use_mat_filt with material filter - it would only work in special cases");
    DensityFilterMat filter_mat;
    filter_mat.designType = design;
    space->density_filter.Push_back(filter_mat);
    // avoid robust case where we have density and multiple filters per element
    bool do_robust = space->density_filter.GetSize() > 1  && space->design.GetSize() == 1;
    int filter_index = do_robust ? space->density_filter.GetSize() - 1 : 0;
    space->density_filter.Last().AssembleFilterMatrix(data,total_sum_neigbors,filter_index,start,end);
  }

  if (space->write_matrix_filt && total_sum_neigbors > 0)
  {
    // writes filter matrix to .mtx file for first filter radius > 1
//    DensityFilterMat filter_mat;
//    space->density_filter.Push_back(filter_mat);
//    int filter_index =space->density_filter.GetSize() - 1 ;
//    // This are the designs we deal with
//    unsigned int start = space->FindDesign(design) * space->GetNumberOfElements(); // handles ALL_DESIGNS
//    unsigned int end = start + space->GetNumberOfElements();
//    //std::vector<DesignElement>   sub(&data[start],&data[end]);
//    //StdVector<DesignElement> data_red(sub);
//    //std::vector<DesignElement>::const_iterator first = data.Begin() + start;
//    //std::vector<DesignElement>::const_iterator last = data.Begin() + end;
//    //std::vector<DesignElement> sub(data.Begin() + start, data.Begin() + end);
//    //DesignElement *arrayOfT = &data[0] + start;
//    //size_t arrayOfTLength = (end - start);
//    //StdVector<DesignElement> data_red(arrayOfT);
//    space->density_filter.Last().AssembleFilterMatrix(data,sum_neighbours,filter_index,start,end);
    space->density_filter.Last().ExportDensityFilterMatrix("filterMat" + std::to_string(space->density_filter.GetSize()-1) + ".mtx");
    // makes sure that matrix is only written once for first filter radius > 1
    if (space->density_filter.GetSize() == space->density_filter.GetCapacity())
      space->write_matrix_filt = false;
  }
}


GlobalFilter DesignStructure::Parse(PtrParamNode pn, const DesignElement* ref_de, bool skip_cfs_filtering)
{
  GlobalFilter global;

  // valid values in xml are 0,1,2

  assert(ref_de != NULL && ref_de->simp != NULL);
  assert(global.robust == -1);

  if(pn->Has("robust_excitation"))
  {
    global.robust = pn->Get("robust_excitation")->As<unsigned int>();
    if(global.robust != (int) ref_de->simp->filter.GetSize())
      EXCEPTION("'robust_excitation' in 'filter' is " << global.robust << " but " << ref_de->simp->filter.GetSize() << " is expected");
  }
  else
  {
    if(!ref_de->simp->filter.IsEmpty())
      throw Exception("expect 'robust_excitation' for 'filter' as more than one filter is defined for the design");
  }
  // we don't set design as it might be allDesigns but set it by filling it in SetFilter()
  global.type = Filter::type.Parse(pn->Get("type")->As<std::string>());
  global.region = ref_de->elem->regionId; // was space->GetRegionIds()[0];

  // for unstructured grids only "radius" filter makes sense
  global.filterspace = Filter::filterSpace.Parse(pn->Get("neighborhood")->As<string>());
  if (!regular && global.filterspace != Filter::RADIUS)
    throw Exception("For non-regular grids filter has to be set with neighborhood='radius'.");

  global.contribution = pn->Get("contribution")->As<string>() == "linear" ? Filter::LINEAR : Filter::CONSTANT;

  global.value  = pn->Get("value")->As<double>();
  if(global.value <= 0.0 || skip_cfs_filtering)
    global.type = Filter::NO_FILTERING;

  if(global.type == Filter::SENSITIVITY && pn->Has("sensitivity"))
    global.sensitivity = Filter::sensitivity.Parse(pn->Get("sensitivity/type")->As<std::string>());

  if(global.type == Filter::DENSITY && pn->Has("density"))
    global.density = Filter::density.Parse(pn->Get("density/type")->As<string>());

  if(global.density == Filter::MATERIAL)
  {
    if(!pn->Has("density/material"))
      throw Exception("density filter 'material' needs 'material' child element");
    global.mat_filter = TransferFunction(pn->Get("density/material/filter"));
    global.mat_scale  = TransferFunction(pn->Get("density/material/scale"));
    global.mat_phase  = TransferFunction(pn->Get("density/material/phase"));
  }

  if(global.density != Filter::STANDARD && global.density != Filter::MATERIAL)
  {
    if(!pn->Has("density/beta"))
      throw Exception("Attribute 'beta' required for '" + Filter::density.ToString(global.density) + "' density filtering");
    string bv = pn->Get("density/beta")->As<string>();
    if(pn->Get("density/beta")->As<string>() == "tune")
    {
      global.tune.Init(pn->Get("density/tune", ParamNode::PASS), Tune::BETA);
      global.beta = global.tune.start;
    }
    else if(pn->Get("density/beta")->As<string>() == "tuned")
    {
      /** we can set global.ext_tune only to a Tune object when we do the final Registration in PostCopy() */
      global.ext_tune = (Tune*) GlobalFilter::UNSET_EXT_BETA_TUNE;
      global.beta = -1;
    }
    else
      global.beta = pn->Get("density/beta")->As<double>();

    assert(space->design.GetSize() == 1); // extend for multiple regions with lower bounds which are now stored in non_lin_*

    if(global.density == Filter::TANH || global.density == Filter::EXPRESSION)
    {
      if(!pn->Has("density/eta"))
        throw Exception("Attribute 'eta' required for '" + Filter::density.ToString(global.density) + "' density filtering");
      global.eta = pn->Get("density/eta")->As<double>();
    }
    if(global.density == Filter::EXPRESSION)
    {
      if(!pn->Has("density/expression"))
        throw Exception("density filter 'expression' requires element 'expression'");
      global.InitParser(pn->Get("density/expression/function")->As<string>(), pn->Get("density/expression/derivative")->As<string>());
    }
  }

  return global;
}


void DesignStructure::WriteFilterInfo(PtrParamNode in)
{
  PtrParamNode root = in->Get("filters");
  if(periodic_)
    root->Get("periodic")->SetValue(periodic_);

  int nlf = 0; // nonlinar filters?
  int mf = 0; // material filter?
  double ar = 0; // average radius
  double an = 0; // average neighbors
  int te = 0; // total elements

  for(GlobalFilter& gf : space->filter)
  {
    gf.ToInfo(root->Get("filter",ParamNode::APPEND));
    switch(gf.density)
    {
    case Filter::MATERIAL:
      mf++;
    case Filter::VOID_HEAVISIDE:
    case Filter::SOLID_HEAVISIDE:
    case Filter::TANH:
      nlf++;
      break;
    default:
      break;
    }
    ar += gf.avg_radius * gf.elements;
    an += gf.avg_neigbor * gf.elements;
    te += gf.elements;
  }

  if(em != NULL && em->constraints.Has(Function::VOLUME))
  {
    if(nlf && em->constraints.Get(Function::VOLUME)->IsLinear())
      root->SetWarning("'volume' constraint shall be non-linear due to non-linear filter");
    if(mf && em->constraints.Get(Function::VOLUME)->GetAccess() != Function::PLAIN)
      root->SetWarning("'volume' constraint for material filter shall have plain access");
  }
  if(too_large_filter_ > 0)
    root->SetWarning("filter radius probably too large");

  root->Get("timer")->SetValue(timer);
  std::cout << "Filter: avg radius=" << (ar/te) << " avg neighbourhood=" << (an/te) << std::endl;
}

void DesignStructure::FindRegularNeighborhood(const GlobalFilter* global, DesignElement* base, double radius, const StdVector<double>& edges, StdVector<Filter::NeighbourElement>& neighbors)
{
  assert(regular);
  // from the radius define a square/cube and check for every element. The corners are sorted out by distance

  // the legacy SHARP_PLAIN and SHARP_SIGMUND had the bug, that the weight was not
  // radius - distance but value - distance. To keep the legacy results we reproduce
  // the bug. Note, that another bug in SHARP_PLAIN and SHARP_SIGMUND is in normalization
  // and for SHARP_SIGMUND also in the filtering itself! :(

  double val_rad = global->sensitivity == Filter::SHARP_PLAIN || global->sensitivity == Filter::SHARP_SIGMUND ? global->value : radius;

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
            assert(global->contribution == Filter::LINEAR || global->contribution == Filter::CONSTANT);
            ne.weight = global->contribution == Filter::LINEAR ? val_rad  - distance : 1;
            ne.distance  = distance;
            neighbors.Push_back(ne); // cheap
          }
        }
      }
    }
  }
}


Filter::Type DesignStructure::GetCommonFilterType() const
{
  for(unsigned int i = 1; i < space->filter.GetSize(); i++)
    if(space->filter[i].type != space->filter[i-1].type)
      throw Exception("inconsistent filter types");

  return space->filter.IsEmpty() ? Filter::NO_FILTERING : space->filter.First().type;
}

Filter::Type DesignStructure::GuessFilterType()
{
  PtrParamNode pn = domain->GetParamRoot();
  if(!pn->Has("optimization/ersatzMaterial/filters"))
    return Filter::NO_FILTERING;

  ParamNodeList list = pn->Get("optimization/ersatzMaterial/filters")->GetList("filter");
  if(list.IsEmpty())
    return Filter::NO_FILTERING;

  assert(Filter::type.map.size() > 0);
  // later GetCommonFilterType() will be called and compared with this guess
  // it would be cooler to move some filter attributes to the common filters element
  return Filter::type.Parse(list.First()->Get("type")->As<string>());
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


void DesignStructure::FindUnstructuredNeighborhood(const GlobalFilter* global, DesignElement* base, double radius, StdVector<Filter::NeighbourElement>& neighbors)
{
  // for 3D ist shall be significantly faster (O(N)) to make a discrete grid (e.g. of size radius) and sort
  // the elements to this grid. Then we can linearly test all relevant grid cells for an element. It will help for
  // large meshes where the filter setup is otherwise hours and magnitudes slower than solving the system!

  // LOG_DBG2(ds) << "FN: base= " << base->elem->elemNum << " initial=" << ToString(initial) << " n=" << ToString(neighbors) << " tf=" << too_far.ToString() << " ext=" << space->DoNonDesignVicinity();

  // the legacy SHARP_PLAIN and SHARP_SIGMUND had the bug, that the weight was not
  // radius - distance but value - distance. To keep the legacy results we reproduce
  // the bug. Note, that another bug in SHARP_PLAIN and SHARP_SIGMUND is in normalization
  // and for SHARP_SIGMUND also in the filtering itself! :(
  double val_rad = global->sensitivity == Filter::SHARP_PLAIN || global->sensitivity == Filter::SHARP_SIGMUND ? global->value : radius;

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
      assert(global->contribution == Filter::LINEAR || global->contribution == Filter::CONSTANT);
      ne.weight = global->contribution == Filter::LINEAR ? val_rad  - distance : 1;
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

/** This is not performance tuned as for almost cases we have regular grids and then this method is
 * only called once. In the other cases - live with it */
double DesignStructure::FindFilterRadius(Filter::FilterSpace space, const DesignElement* de, double value)
{
  Matrix<double>  coords;
  domain->GetGrid()->GetElemNodesCoord(coords, de->elem->connect, false );

  switch(space)
  {
    case Filter::RADIUS:
      return value;

    case Filter::VOLUME_RADIUS:
    {
      // TODO really check for axis symmetry off
      double tmp = domain->GetGrid()->GetElemShapeMap(de->elem, false)->CalcVolume();
      // The radius is <value> times square/cube edge length where the
      // square/cube has the volume of the element
      double radius = value * std::pow(tmp, 1.0/ (double) domain->GetGrid()->GetDim());
      LOG_DBG3(ds) << "FFR: de=" << de->ToString() << " from volume " << tmp << " to radius " << radius;
      return radius;
    }

    case Filter::MAX_EDGE:
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
