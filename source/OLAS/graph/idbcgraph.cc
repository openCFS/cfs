#include "graph/idbcgraph.hh"

namespace OLAS {

  void IDBC_Graph::FinaliseAssembly( Integer *newEqn ) {

    ENTER_FCN( "IDBC_Graph::FinaliseAssembly" );

    // Re-ordering strategy should be not to re-order
    if ( newOrder_ != NOREORDERING ) {
      (*warning) << "IDBC_Graph::FinaliseAssembly: Re-ordering strategy = '"
                 << Enum2String( newOrder_ )
                 << "' makes no sense in my case! Ignoring it!";
      Warning( __FILE__, __LINE__ );
    }

    // Determine number of non-zero entries
    CountNNE();

    // Check whether we got a NULL pointer. In this case there is
    // no re-ordering necessary
    bool doReorder = newEqn == NULL ? false : true;

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
      NodeList *newElement;
      NewArray( newElement, NodeList, numNodes_ );
      for ( UInt i = 1; i <= numNodes_; i++ ) {
        newElement[ newEqn[i] ] = element_[i];
      }
      NodeList *tmpPtr = element_;
      element_ = newElement;
      newElement = NULL;
      DeleteArray( tmpPtr );

#ifdef DEBUG_IDBCGRAPH
      (*debug) << " IDBC_Graph: New numbers of free dofs:\n";
      for ( UInt i = 1; i <= numNodes_; i++ ) {
        (*debug) << i << " -> " << newEqn[i] << std::endl;
      }
#endif

    }

    // Sort lists and remove duplicate entries
    SortLists();

    // Convert Storage format to CRS style
    ConvertToCRS();

    // The element vector is no longer required
    DeleteArray( element_ );
    element_ = NULL;

    // Now the graph object is fully assembled
    amAssembled_ = true;

  }

}
