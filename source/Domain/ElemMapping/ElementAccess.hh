#ifndef OPTIMIZATION_ELEMENTACCESS_HH_
#define OPTIMIZATION_ELEMENTACCESS_HH_

#include "Utils/StdVector.hh"
#include "Domain/ElemMapping/ElemShapeMap.hh"
#include "Forms/IntScheme.hh"

namespace CoupledField
{
struct Elem;
class BaseBDBInt;
class SinglePDE;
class BiLinFormContext;

/** This is a helper for dealing with FE-Elements, e.g. because the differential operator (BOp) which is usually to be evaluated at each integration point of an element is required.
 * Within the FE procedure itself this helper is probably not necessary but for other parts (e.g. optimization) it might be helpful to avoid code duplication, e.g. holding
 * the objects for integration points, weights, nodal equations, ... */
struct ElementAccess
{
public:
  /** @param comment is just an optional additional information for logging */
  ElementAccess(SinglePDE* pde, const std::string& bilinform, RegionIdType reg);

  /** extracts the first pde */
  ElementAccess(BiLinFormContext* context);

  /** Set all information from SetElem but also esm, stores the info in curr_elem. resets curr_ip.
   * TODO intPoints and weights could be cached!
   * @param elem if NULL it is marked as unset, resets also ip
   * @return @see CurrElemShapeMap() */
  ElemShapeMap* SetElem(const Elem* elem);  // SetElem()

  /** optionally set the integration point. Sets also esm content! */
  LocPointMapped& SetIP(unsigned int ip);

  /** last SetElem() or NULL */
  const Elem* CurrElem() const { return curr_elem; }

  /** Set via SetElem(). */
  BaseFE* CurrBaseFE() { return curr_base_fe; }

  /** @return note that is the content of the shared_ptr esm. Don't store the pointer!!! */
  ElemShapeMap* CurrElemShapeMap() { return esm.get(); }

  /** @return NULL if there was no SetIP after the last SetElem() */
  LocPoint* CurrIP() { return curr_ip; }

  BaseBDBInt* GetBDBInt() { return bdb; }

  /** helper for dumping information.
   * 1: about general state
   * 2: about elem and ip
   * 3: very verbose */
  std::string ToString(int level) const;

  /** TODO is not but could be specific to current element type */
  StdVector<LocPoint> intPoints; // Get integration Points
  StdVector<double> weights;
  IntegOrder order;

  /** specific to current element. the elem shape map contains the Jacobi determinant. */
  shared_ptr<ElemShapeMap> esm;

  /** 1-based element equation numbers, set by SetElem(). 0 = HDBC, < 0 is constrained. Not with empty constructor! */
  StdVector<int> elem_eqn_nr;

  /** 0-based element equation indices. < 0 is HDBC or constrained */
  StdVector<int> elem_eqn_idx;

  /* specific to current integration point */
  LocPointMapped lpm;

private:
  IntScheme::IntegMethod method = IntScheme::UNDEFINED;

  /** element specific */
  const Elem* curr_elem = NULL; // current element set by SetElem()

  /** set by SetElem() */
  BaseFE* curr_base_fe = NULL;

  /** ip specific. Reset to NULL for new elem */
  LocPoint* curr_ip = NULL;

  /** general stuff */
  BaseBDBInt* bdb = NULL;
  SinglePDE* pde = NULL;

};



} // end of namespace

#endif /* OPTIMIZATION_ELEMENTACCESS_HH_ */
