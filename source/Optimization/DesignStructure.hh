#ifndef DESIGNSTRUCTURE_HH_
#define DESIGNSTRUCTURE_HH_

#include "Optimization/DesignElement.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField
{
class ErsatzMaterial;
class Timer;

/** This class has two purposes: it sets the the filter for DesignElement (sensitivity filter)
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
   * relevant design elements of simp->design.data
   * @param pn our parameter section of the filer, includes design name
   * @param info where to write (in HEADER) the information */
  void SetFilters(PtrParamNode pn, PtrParamNode info);

  /** Do we have periodic BC. */
  bool IsPeriodic() const { return periodic; }

  /** This is the relevant service function for VicinityElement::Init() -> check IsPeriodic().
   *  if elem has nodes in the constraintMapping, neighbors is extended by all neighbors of
   * the corresponding periodic elements.
   * Also the original elem neighbors are added, ONLY if periodic within elem!
   * @param min_common
   * @return false means not periodic and and neighbors is 0 size */
  bool ExtendPeriodicNeighborhood(Elem* elem, int min_common, StdVector<std::pair<Elem*, int> >& neighbors);

private:

  /** The common Constructor, does much less than Initialize() */
  void Constructor();

  /** Filter types we have
   * <ul>
   *   <li>RADIUS: this is the implementation following Sigmund in the 99lines paper.
   *               The drawback is the discretization dependency.</li>
   *   <li>VOLUME_RADIUS: The radius is *value* times square/cube edge length where the
   *               square/cube has the volume of the element</li>
   *   <li>MAX_EDGE: The largest edge size, discretization independent and preferable</li>
   * </ul>
   * This does not tell if we have design or sensitivity filtering! */
  typedef enum { NO_FILTER = -1, RADIUS, VOLUME_RADIUS, MAX_EDGE } FilterSpace;


  /** The actual constructor, initializes for SetFilter() and ExtendPeriodicNeighborhood() on
   * the fly! The time is not recorded!
   * @see initialized_ */
  void Initialize();

  /** Helper for InitFilter(). Call only once if regular! */
  double FindFilterRadius(FilterSpace filter, DesignElement* de) const;

  /** Almost copy and paste to the version with common */
  bool ExtendPeriodicNeighborhood(Elem* elem, StdVector<std::pair<Elem*, int> >& neighbors);

  /** This is a helper for InitFilter(). It is recursive!.
   * See implementation for docu. */
  void FindNeighborhood(DesignElement* base, double radius,
                          StdVector<std::pair<Elem*, int> >& buddies,
                          StdVector<SIMPElement::NeighbourElement>& neighbors,
                          StdVector<unsigned int>& too_far);

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


  /** adds source to out such that no double entries exist.
   * Helper for ExtendPeriodicNeighborhood() */
  void AppendNeighbors(const StdVector<std::pair<Elem*, int> >& source, StdVector<std::pair<Elem*, int> >& out);

  /** helper for the common stuff.
   * Almost copy and paste to the other AppendNeighbors() */
  void AppendNeighbors(Elem* check,
                         const StdVector<unsigned int>& constraints, int min_common,
                         StdVector<std::pair<Elem*, int> >& out);

  /** Helper for LOG_DBG() */
  std::string ToString(StdVector<SIMPElement::NeighbourElement>& data);

  /** the way of the weighting in the filter. CONSTANT e.g. for MAX filter */
  enum Contribution { LINEAR, CONSTANT };

  /** Parameter for filter */
  Contribution contribution_;

  Enum<FilterSpace> filterSpace;

  /** actual filter space setting */
  FilterSpace filter_space_;

  /** the filter value */
  double value;

  /** do we have periodic boundary conditions? */
  bool periodic;

  /** is the grid regular? */
  bool regular;

  /** The reference filter setting */
  Filter filter_;

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
  unsigned int dim;

  ErsatzMaterial* em;

  /** is Initialize() called yet? */
  bool initialized_;

  DesignSpace* space;
  StdVector<RegionIdType> regions;
};


} // end of namespace


#endif /* DESIGNSTRUCTURE_HH_ */
