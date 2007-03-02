#include <iostream>
#include <cmath>

#include "graph/baseordering.hh"

namespace OLAS {


  // ***************
  //   Constructor
  // ***************
  BaseOrdering::BaseOrdering( NodeList *graph, Integer *order,
			      Integer size ) {

    ENTER_FCN( "BaseOrdering::BaseOrdering" );
    graph_ = graph;
    order_ = order;
    size_  = (UInt)size;

    // memory has been allocated by the calling routine set to zero
    for ( UInt i = 1; i <= size_; i++ ) {
      order[i] = 0;
    }
  }


  // **************
  //   Destructor
  // **************
  BaseOrdering::~BaseOrdering() {
    ENTER_FCN( "BaseOrdering::~BaseOrdering" );
  }


} // namespace
