#include <iostream>
#include <cmath>

#include "General/Enum.hh"
#include "BaseOrdering.hh"

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
  BaseOrdering::BaseOrdering(StdVector<std::vector<UInt>>& graph, StdVector<UInt>& order ) :
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
