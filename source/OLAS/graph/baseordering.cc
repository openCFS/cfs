// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <cmath>

#include "graph/baseordering.hh"

namespace OLAS {


  // ***************
  //   Constructor
  // ***************
  BaseOrdering::BaseOrdering( NodeList *graph, Integer *order,
			      Integer size ) {

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
  }


} // namespace
