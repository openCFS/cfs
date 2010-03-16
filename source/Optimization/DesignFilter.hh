#ifndef DESIGNFILTER_HH_
#define DESIGNFILTER_HH_

#include "Optimization/DesignElement.hh"

namespace CoupledField
{
class ParamNode;
class ParamNode;
class SIMP;
class Timer;

/** The only purpose of this tool class is, to set up the filter for DesignElement.
 * This is extracted here as it becomes more complex with periodic boundary conditions.
 * Give this element only a short live time as there are expensive arrays */
 class DesignFilter
{
public:

  /** Is already a little expensive. Give this object a short live time!
   * @param simp refers to the data and has pde reference for bcs.
   * @param pn our parameter section, includes design name */
  DesignFilter(SIMP* simp, PtrParamNode pn);

  /** This is the actual, expensive action. It sets the filters for all
   * relevant design elements of simp->design.data */
  void SetFilters();

private:

  /** Helper for InitFilter(). Call only once if regular! */
  double FindFilterRadius(DesignElement::Filter filter, DesignElement* de) const;

  /** This is a helper for InitFilter(). It is recursive!.
   * See implementation for docu. */
  void FindNeighborhood(DesignElement* base, double radius,
                          StdVector<std::pair<Elem*, int> >& buddies,
                          StdVector<SIMPElement::NeighbourElement>& neighbors,
                          StdVector<unsigned int>& too_far);

  /** if elem has nodes in the constraintMapping, neighbors is extended by all neighbors of
   * the correspinding periodic elements.
   * Also the originial elem neighbors are added, ONLY if periodic
   * @return false means not periodic and thend neighbors is 0 size */
  bool ExtendPeriodicNeighborhood(Elem* elem, StdVector<std::pair<Elem*, int> >& neighbors);

  /** calc the distance between two points for the periodic case,
   * where periodic boundaries are considered.
   * @param base can on, close or far away from a periodic boundary node element
   * @param test element
   * @return a minimal distance with considered periodic boundaries */
  double RelaxedDistance(const Elem* base, const Elem* test) const;

  /** set the vector contraintMapping */
  void SetPeriodicConstraintMapping();

  /** set vector nodeToElem for periodic only */
  void SetNodeElemMapping();


  /** adds source to out such that no double entries exist.
   * Helper for ExtendPeriodicNeighborhood() */
  void AppendNeighbors(const StdVector<std::pair<Elem*, int> >& source, StdVector<std::pair<Elem*, int> >& out);

  /** Helper for LOG_DBG() */
  std::string ToString(StdVector<SIMPElement::NeighbourElement>& data);

  /** the filter value */
  double value;

  /** do we have periodic boundary conditions? */
  bool periodic;

  /** is the grid regular? */
  bool regular;

  /** our filter mode */
  DesignElement::Filter filter;

  /** maps elem_num to design element pointer. Large! */
  StdVector<DesignElement*> elemToDesign;

  /** This maps the periodic boundies. For every mapped node position the date is the other node.
   * master -> slave is also slave -> master. Set by SetPeriodicConstraintMapping() if periodic.
   * Entry 0 means no content! */
  StdVector<unsigned int> constraintMapping;

  /** maps from the node number to one possible element.
   * Only Required for periodic and set by SetNodeElemMapping() */
  StdVector<Elem*> nodeToElem;

  /** diameter of convex assumed domain */
  Point dimension;

  /** shortcut to the grid dimension */
  unsigned int dim;

  SIMP* simp;

  PtrParamNode info_;

  /** this timer is stored within the info_ element */
  Timer* timer_;
};


} // end of namespace


#endif /* DESIGNFILTER_HH_ */
