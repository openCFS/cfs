// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <cmath>

#include "General/Enum.hh"
#include "baseordering.hh"

namespace CoupledField {

  static EnumTuple reorderingTypeTuples[] = 
  {
    EnumTuple( BaseOrdering::NOREORDERING, "noReordering" ),
    EnumTuple( BaseOrdering::SLOAN, "Sloan" ),
    EnumTuple( BaseOrdering::METIS, "Metis"),
    EnumTuple( BaseOrdering::MINIMUM_DEGREE, "minimumDegree"),
    EnumTuple( BaseOrdering::NESTED_DISSECTION, "nestedDissection" ),
  };

  Enum<BaseOrdering::ReorderingType> BaseOrdering::reorderingType = \
  Enum<BaseOrdering::ReorderingType>("Reordering Types",
      sizeof(reorderingTypeTuples) / sizeof(EnumTuple),
      reorderingTypeTuples); 
      
  // ***************
  //   Constructor
  // ***************
  BaseOrdering::BaseOrdering( NodeList *graph, StdVector<UInt>& order ) :
    graph_(graph),
    order_(order)
  {
    // take care: graph and order are zero based!

    // memory has been allocated by the calling routine set to zero
    for ( UInt i = 0, n=order.GetSize(); i < n; i++ ) {
      order_[i] = 0;
    }
  }


  // **************
  //   Destructor
  // **************
  BaseOrdering::~BaseOrdering() {
  }


} // namespace
