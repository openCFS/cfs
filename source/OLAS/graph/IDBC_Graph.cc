#include "OLAS/graph/IDBC_Graph.hh"

namespace CoupledField {

  void IDBC_Graph::FinaliseAssembly( StdVector<UInt>* newEqn ) {

    MapSetToVector();

    // Re-ordering strategy should be not to re-order
    if ( newOrder_ != BaseOrdering::NOREORDERING ) {
      std::string tmp;
      tmp = BaseOrdering::reorderingType.ToString( newOrder_ );

      WARN("IDBC_Graph::FinaliseAssembly: Re-ordering strategy = '"
           << tmp
           << "' makes no sense in my case! Ignoring it!");
    }

    // Determine number of non-zero entries
    CountNNE();

    // Check whether we got a NULL pointer. In this case there is
    // no re-ordering necessary
    bool doReorder = newEqn->GetSize() == 0 ? false : true;

    // Only do re-order, if graph is not empty
    doReorder = numNodes_ > 0 ? doReorder : false;

#ifdef DEBUG_IDBCGRAPH
    (*debug) << " IDBC_Graph:\n"
             << " graph contains " << numNodes_ << " vertices\n"
             << " with " << nne_ << " neighbours\n";
    if ( doReorder == true ) {
      (*debug) << " Going to sort vertices now" << std::endl;
    }
    else if ( newEqn == NULL ) {
      (*debug) << " newEqn = NULL, so no new vertex numbers!" << std::endl;
    }
    else {
      (*debug) << " No vertices -> no re-ordering!" << std::endl;
    }
#endif

    // Re-arrange neighbour list with respect to the new equation
    // numbers of the free dofs
    if ( doReorder == true ) {

      StdVector<std::vector<unsigned int>> newElement(numNodes_);

      for ( UInt i=0; i< numNodes_; i++ ) {
        UInt n = (*newEqn)[i];
        newElement[n-1] = element_[i];
      }

      for ( UInt i=0; i< numNodes_; i++ ) {
        element_[i] = newElement[i];
      }

#ifdef DEBUG_IDBCGRAPH
      (*debug) << " IDBC_Graph: New numbers of free dofs:\n";
      for ( UInt i = 0; i < numNodes_; i++ ) {
        (*debug) << i << " -> " << (*newEqn)[i] << std::endl;
      }
#endif

    }

    // Sort lists and remove duplicate entries
    SortLists();

    // Convert Storage format to CRS style
    ConvertToCRS();

    // The element vector is no longer required
    element_.Clear(false); // don't keep capcity

    // Now the graph object is fully assembled
    amAssembled_ = true;

  }

}
