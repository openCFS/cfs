// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef OLAS_BASEORDERING_HH
#define OLAS_BASEORDERING_HH

#include <iostream>
#include <vector>
#include "utils/utils.hh"
#include "utils/plainAlloc.hh"


namespace OLAS {


  //! Base Class for ordering of a matrix graph 

  //! (notation: ordering of nodes although in combination with
  //! CFS++ we order equation numbers).
  class BaseOrdering {

  public:
    
    //! Shortcut for an STL vector of unsigned integers (UInt)
    typedef std::vector<UInt> NodeList;

    //! \name Constructors and Destructors
    //@{
  
    //! Constructor
    //! \param graph the graph to be re-ordered stored as array of linked
    //!              lists
    //! \param order contains for global node i the mapped index
    //! \param asize number of elements (nodes) in the graph
    BaseOrdering( NodeList *graph, Integer *order, Integer asize );

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

    NodeList* graph_;  //!< Uncompressed matrix graph (STL list)
    Integer *order_;   //!< stores for global node i the mapped index
    UInt size_;        //!< Total number of nodes in matrix graph

  };

} // namespace

#endif // OLAS_BASEORDERING_HH
