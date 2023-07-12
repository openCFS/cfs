#ifndef DESIGNSTRUCTURE_HH_
#define DESIGNSTRUCTURE_HH_

#include <stddef.h>
#include <string>
#include <utility>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/Enum.hh"
#include "General/Environment.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/Filter.hh"
#include "Utils/Point.hh"
#include "Utils/StdVector.hh"

namespace CoupledField {
class DesignSpace;
class Grid;
class GridCFS;
struct Elem;
}  // namespace CoupledField

namespace CoupledField
{
class ErsatzMaterial;

/** This class has two purposes: it sets the filter for DesignElement (density and sensitivity filter)
 * and it provides the periodic b.c. for VinicityElement (slopes).
 * It has expensive arrays. You might delete them !*/
 class DesignStructure
{
public:

  /** Very cheap, the work is done by Initialize() only on SetFilter() or ExtendPeriodicNeighborhood() */
  DesignStructure(ErsatzMaterial* em);

  /** Alternative constructor when no ErsatzMaterial* is available (load ersatz material),
   * then periodic b.c.s are not set */
  DesignStructure(DesignSpace* space, StdVector<RegionIdType>& regions);

  /** This is the actual, expensive action. It sets the filters for all

   * relevant design elements of simp->design.data.
   * Might be called multiple times for multiple filter definitions.
   * @param pn our parameter section of the filer, includes design name (might be allDesigns) 
   * @param skip_cfs_filtering is a shortcut to prevent cfs from applying (only assembling) the declared filters in xml-file (e.g. for SGP)*/
  void SetFilter(PtrParamNode pn, bool skip_cfs_filtering = false);

  /** Service function to be called after the SetFilter() calls.
   * Searches common type in DesignSpace::filter */
  Filter::Type GetCommonFilterType() const;

  /** for ugly post init ordering chaos we need to know the filter type very early in the
   * Function constructors. Read this manually from param. To be confirmed later in Function::PostProc() */
  static Filter::Type GuessFilterType();

  /** Is a ToInfo() which creates a filters element but also prints avg data on command line */
  void WriteFilterInfo(PtrParamNode in);

  /** Do we have periodic BC. */
  bool IsPeriodic() const { return periodic_; }

  /** This is the relevant service function for VicinityElement::Init() -> check IsPeriodic().
   *  if elem has nodes in the constraintMapping, neighbors is extended by all neighbors of
   * the corresponding periodic elements.
   * Also the original elem neighbors are added, ONLY if periodic within elem!
   * @param min_common
   * @return false means not periodic and and neighbors is 0 size */
  bool ExtendPeriodicNeighborhood(Elem* elem, int min_common, StdVector<std::pair<Elem*, int> >& neighbors);


  /** Helper for DesignStructure::InitFilter() and Function::Local::SetupStarLocalityElementMap.
   * Call only once if regular as it is not performance tunes! */
  static double FindFilterRadius(Filter::FilterSpace filter, const DesignElement* de, double value);

  static std::string ToString(const StdVector<std::pair<Elem*, int> >& data);

  static std::string ToString(const StdVector<Filter::NeighbourElement>& data);



private:

  /** The common Constructor, does much less than Initialize() */
  void Constructor();

  /** The actual constructor, initializes for SetFilter() and ExtendPeriodicNeighborhood() on
   * the fly! The time is not recorded!
   * @see initialized_ */
  void Initialize();

  /** parses s single filter element. Service for SetFilter()
   * The return value is created and there shall be no copy operation necessary
   * @param skip_cfs_filtering: assemble filter but don't apply it to design, e.g. pass to external SGP solver*/
  GlobalFilter Parse(PtrParamNode pn, const DesignElement* ref_de, bool skip_cfs_filtering);

  /** finds quite efficiently the neighborhood with an regular grid.
   * The idea is that by radius and edge size we construct a 2D/3D cube and check every element for distance.
   * @param neighbors to be reused */
  void FindRegularNeighborhood(const GlobalFilter* global, DesignElement* base, double radius, const StdVector<double>& edges, StdVector<Filter::NeighbourElement>& neighbors);

  /** A helper for InitFilter() with a recursive like implementation. See implementation. */
  void FindUnstructuredNeighborhood(const GlobalFilter* global, DesignElement* base, double radius, StdVector<Filter::NeighbourElement>& neighbors);

  /** Helper for FindRegularNeighborhood().
   * Defines an element by the number of (+/-) steps in the main axes
   * @return NULL if nothing found. Must not be in the periodic case */
  DesignElement* GetNeighborElement(DesignElement* base, int i_steps, int j_steps, int k_steps);

  /** Helper for the other GetNeighborElement().
   * Is able to cross periodic boundaries */
  DesignElement* GetNeighborElement(DesignElement* base, unsigned int steps, VicinityElement::Neighbour dir);

  /** calc the distance between two points for the periodic case,
   * where periodic boundaries are considered.
   * @param base can on, close or far away from a periodic boundary node element
   * @param test element
   * @return a minimal distance with considered periodic boundaries */
  double RelaxedDistance(const Elem* base, const Elem* test) const;

  /** set the vector contraintMapping */
  void SetPeriodicConstraintMapping();

  /** Recursive helper for SetPeriodicConstraintMapping() */
  void RecursiveCompletePeriodicity(unsigned int master, StdVector<unsigned int>& list, unsigned int test);

  /** set vector nodeToElem for periodic only */
  void SetNodeElemMapping();

  /** append all new stuff */
  void AppendNeighbors(const StdVector<std::pair<Elem*, int> >& source,
                       StdVector<std::pair<Elem*, int> >& out);

  /** Potentially appends the element to the neighbor
   * Almost copy and paste to the other AppendNeighbors() */
  void AppendNeighbors(Elem* check,
                       const StdVector<unsigned int>& constraints,
                       int min_common,
                       StdVector<std::pair<Elem*, int> >& out);

  /** do we have periodic boundary conditions? */
  bool periodic_ = false;

  /** is the grid regular? */
  bool regular = false;

  /** This maps the periodic boundaries.
   * Nodes which are not part of a periodic b.c. have an empty vector.
   * Nodes with exactly one periodic partner are like: n_1 <-> n_2: constraintMapping[n_1][0] = n_2
   * In 2D the corners and in 3D also all edges have more than one corresponding node.
   * In 3D do the corners have 7 mapping nodes. */
  StdVector<StdVector<unsigned int> > constraintMapping;

  /** maps from the node number to one possible element.
   * Only Required for periodic and set by SetNodeElemMapping() */
  StdVector<Elem*> nodeToElem;

  /** Little helper to be reused in ExtendPeriodicNeighborhood() */
  StdVector<unsigned int> constraintNodes_;

  /** diameter of convex assumed domain */
  Point dimension;

  /** shortcut to the grid dimension */
  unsigned int dim = 0;

  Grid* grid = NULL;
  /** shortcut to grid cfs to save virtual table lookup for GetElem() */
  GridCFS* gridcfs = NULL;

  ErsatzMaterial* em = NULL;

  /** is Initialize() called yet? */
  bool initialized_ = false;

  unsigned int num_robust_ = 0;

  DesignSpace* space = NULL;
  StdVector<RegionIdType> regions;

  /** for warnings */
  int too_large_filter_ = 0;

  shared_ptr<Timer> timer;
};


} // end of namespace


#endif /* DESIGNSTRUCTURE_HH_ */
