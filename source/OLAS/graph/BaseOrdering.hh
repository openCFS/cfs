#ifndef OLAS_BASEORDERING_HH
#define OLAS_BASEORDERING_HH

#include <iostream>
#include <vector>

#include "General/defs.hh"
#include "Utils/StdVector.hh"

namespace CoupledField {


  //! Base Class for ordering of a matrix graph 

  //! (notation: ordering of nodes although in combination with
  //! openCFS we order equation numbers).
  class BaseOrdering {
  public:
  //! Type of re-ordering for the graph

  //! This enumeration data type describes the type of re-ordering that is to
  //! be applied to the graph. The enumeration contains the following values:
  //! - NOREORDERING
  //! - METIS
  //! - SLOAN
  //! - MINIMUM_DEGREE
  //! - NESTED_DISSECTION
  typedef enum { NOREORDERING, METIS, SLOAN, NESTED_DISSECTION,
                 MINIMUM_DEGREE } ReorderingType;
  static Enum<ReorderingType> reorderingType;

  public:
    
    //! \name Constructors and Destructors
    //@{
  
    //! Constructor
    //! \param graph the graph to be re-ordered stored as array of linked
    //!              lists
    //! \param order contains for global node i the mapped index
    //! \param asize number of elements (nodes) in the graph
    BaseOrdering(StdVector<std::vector<UInt>>& graph, StdVector<UInt>& order );

    //! Destructor
    virtual ~BaseOrdering();

    //@}

    //! perform the reordering and compute the new node numbers
    virtual void LabelGraph(){;};

    //! computes the total profile (=sum of bandwidth of each row)

    //! This method can be used to compute both profiles, that of the old
    //! graph and that of the re-ordered graph. We return the profile sizes
    //! as floating point values to avoid problems with limited size of
    //! integral types.
    virtual void GetProfile( Double &oldprof, Double &newprof ) = 0;

  protected:

    StdVector<std::vector<UInt>>& graph_;  //!< Uncompressed matrix graph (STL list)
    StdVector<UInt>& order_;   //!< stores for global node i the mapped index

  };

} // namespace

#endif // OLAS_BASEORDERING_HH
