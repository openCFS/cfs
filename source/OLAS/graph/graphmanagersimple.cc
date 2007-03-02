// ===========================================================================
//
// For debugging:
//
// Definition of DEBUG_GRAPHMANAGERSIMPLE1 turns on some pieces of code for
// debugging which should be uncritical with respect to performance and do
// not log stuff to (*debug)
//
// Definition of DEBUG_GRAPHMANAGERSIMPLE2 turns on some pieces of code for
// debugging which are critical with respect to performance and should only
// be turned on for complete debugging
//
#ifdef DEBUG_GRAPHMANAGERSIMPLE
#define DEBUG_GRAPHMANAGERSIMPLE1
#define DEBUG_GRAPHMANAGERSIMPLE2
#else
#define DEBUG_GRAPHMANAGERSIMPLE1
#endif
//
// ===========================================================================


#include <iomanip>
#include "graph/graphmanagersimple.hh"


namespace OLAS {


  // ===============
  //   Constructor
  // ===============
  GraphManagerSimple::GraphManagerSimple() {

    ENTER_FCN( "GraphManagerSimple::GraphManager" );

    myPDEIdent_          = NO_PDE_ID;
    graph_               = NULL;
    graphIDBC_           = NULL;
    newOrdering_         = NULL;
    newOrderingPassedOn_ = false;
    reorderingDone_      = false;
    numLastFreeDof_      = 0;
    numEqns_             = 0;
  }


  // ==============
  //   Destructor
  // ==============
  GraphManagerSimple::~GraphManagerSimple() {
    ENTER_FCN( "GraphManagerSimple::~GraphManager" );

    // If no re-ordering was performed the pointer to the re-ordering vector
    // is still NULL. If it was claimed by the PDE it was re-set to NULL and
    // if it was not claimed, it's our responsibility to delete it now.
    if ( newOrdering_ != NULL ) {
      if ( newOrderingPassedOn_ == false ) {
	(*warning) << "GraphManagerSimple: I was told to perform a '"
                   << Enum2String(reorderType_) << "' re-ordering of the "
                   << "graph, but nobody ever claimed the permutation "
                   << "vector! Assuming it's my task to de-allocate the "
                   << "memory!";
	Warning( __FILE__, __LINE__ );
      }
      DeleteArray( newOrdering_ );
    }

    // Delete the graph objects
    delete graph_;
    delete graphIDBC_;

  }


  // =============
  //   SetupInit
  // =============
  void GraphManagerSimple::SetupInit( UInt numPDEs ) {

    ENTER_FCN( "GraphManagerSimple::SetupInit" );

    // Make sure that there is only one PDE, since we are the graph manager
    // for the simple case
    if ( numPDEs != 1 ) {
      (*error) << "GraphManagerSimple::SetupInit: This is the graph manager "
               << "for the simple case of one single PDE! However SetupInit "
               << "was informed that there are " << numPDEs << " PDEs!";
      Error( __FILE__, __LINE__ );
    }
  }


  // ===============
  //   RegisterPDE
  // ===============
  void GraphManagerSimple::RegisterPDE( const PdeIdType identifierPDE,  
                                        const UInt numEqns,
                                        const UInt numLastFreeDof,
					const ReorderingType reorder ) {

    ENTER_FCN( "GraphManagerSimple::RegisterPDE" );

#ifdef DEBUG_GRAPHMANAGERSIMPLE1

    // Avoid mis-use
    if ( myPDEIdent_ != NO_PDE_ID ) {
      (*error) << "GraphManagerSimple: A PDE with identifier '" << myPDEIdent_
	       << "' has already been registered with this instance! Refusing "
	       << "to register PDE " << identifierPDE;
      Error( __FILE__, __LINE__ );
    }

#endif

    // Make a sensibility check to avoid problems in the long-run
    if ( numLastFreeDof == 0 ) {
      (*error) << "GraphManagerSimple::RegisterPDE: You tried to register "
               << "a PDE with identifier '" << myPDEIdent_ << "' that has "
               << "numEqns = " << numEqns << ", but numLastFreeDof = "
               << numLastFreeDof << ". If this is really intended, please "
               << "set <setup idbcHandling=\"penalty\"/> in your xml-file.";
      Error( __FILE__, __LINE__ );
    }

    // Store identifier
    myPDEIdent_ = identifierPDE;

    // Store sizes
    numEqns_        = numEqns;
    numLastFreeDof_ = numLastFreeDof;

    // Generate graph object
    graph_ = New BaseGraph( numLastFreeDof_, numLastFreeDof_, reorder );
    AssertMem( graph_, sizeof(BaseGraph) );

    // Generate graph object for IDBC
    UInt graphSize = 0;
    if ( numLastFreeDof_ < numEqns_ ) {
      graphSize = numLastFreeDof_;
      graphIDBC_ = New IDBC_Graph( graphSize, numEqns_ - numLastFreeDof );
      AssertMem( graphIDBC_, sizeof(IDBC_Graph) );

#ifdef DEBUG_GRAPHMANAGERSIMPLE2
      (*debug) << " Generated IDBC_Graph with " << graphSize << " vertices"
               << std::endl;
#endif
    }

    // Allocate memory for the re-ordering and store re-ordering type
    reorderType_ = reorder;
    if ( reorderType_ != NOREORDERING ) {
      NewArray( newOrdering_, Integer, numEqns_ );
    }

  }


  // ================
  //   AssembleInit
  // ================
  void GraphManagerSimple::AssembleInit( const PdeIdType identifierPDE1,
					 const PdeIdType identifierPDE2,
                                         bool assemblingTranspose ) {

    ENTER_FCN( "GraphManagerSimple::AssembleInit" );

#ifdef DEBUG_GRAPHMANAGERSIMPLE1
    CheckConsistency( identifierPDE1, identifierPDE2, "AssembleInit" );
#endif

    // Nothing else to be done, since graph object itself does not
    // require a special call before the start of the assembly

  }


  // ================
  //   AssembleDone
  // ================
  void GraphManagerSimple::AssembleDone( const PdeIdType identifierPDE1,
					 const PdeIdType identifierPDE2,
                                         bool assemblingTranspose ) {

    ENTER_FCN( "GraphManagerSimple::AssembleDone" );

#ifdef DEBUG_GRAPHMANAGERSIMPLE1
    CheckConsistency( identifierPDE1, identifierPDE2, "AssembleDone" );
#endif

  }


  // ============
  //   GetGraph
  // ============
  BaseGraph* GraphManagerSimple::GetGraph( const PdeIdType identifierPDE1,
					   const PdeIdType identifierPDE2 ){

    ENTER_FCN( "GraphManagerSimple::GetGraph" );

    // Avoid trouble (keep defaults in mind)
    if ( identifierPDE1 == NO_PDE_ID ) {
      CheckConsistency( myPDEIdent_, identifierPDE2, "GetGraph" );
    }
    else {
      CheckConsistency( identifierPDE1, identifierPDE2, "GetGraph" );
    }

    // Return the graph object for real dofs
    return graph_;
  }


  // ================
  //   GetIDBCGraph
  // ================
  BaseGraph* GraphManagerSimple::GetIDBCGraph( const PdeIdType pdeID1,
                                               const PdeIdType pdeID2 ) const{

    ENTER_FCN( "GraphManagerSimple::GetIDBCGraph" );

#ifdef DEBUG_GRAPHMANAGERSIMPLE1
    if ( pdeID1 != NO_PDE_ID && pdeID2 != NO_PDE_ID ) {
      CheckConsistency( pdeID1, pdeID2, "GetIDBCGraph" );
    }
#endif

    // Return the graph object for dofs fixed by inhomogeneous Dirichlet
    // boundary conditions
    return graphIDBC_;
  }


  // =================
  //   SetElementPos
  // =================
  void GraphManagerSimple::SetElementPos( const PdeIdType identifierPDE1,
					  Integer *connect1,
					  Integer length1,
					  const PdeIdType identifierPDE2,
					  Integer *connect2,
					  Integer length2,
                                          bool setCounterPart ) {

    ENTER_FCN( "GraphManagerSimple::SetElementPos" );

#ifdef DEBUG_GRAPHMANAGERSIMPLE2
    CheckConsistency( identifierPDE1, identifierPDE2, "SetElementPos" );
#endif

    // -------------------------------------------------------
    // Transform and split the connect arrays to format usable
    // with the graph objects
    // -------------------------------------------------------

    // Clear the arrays
    vertexList1_.clear();
    vertexList2_.clear();
    edgeList1_.clear();
    edgeList2_.clear();

    // STEP 1: Generate vertex list from first connect array, dropping
    //         equation numbers for dofs fixed by inhomogeneous Dirichlet
    //         boundary conditions and changing the sign of those fixed by
    //         constraints
    UInt aux;
    for ( UInt i = 1; i <= length1; i++ ) {
      aux = std::abs( connect1[i] );
      if ( aux <= numLastFreeDof_ && aux > 0 ) {
        vertexList1_.push_back( aux );
      } else {
        vertexList2_.push_back( aux - numLastFreeDof_ );
      }
    }

    // STEP 2: Split the second connect array into two edge lists, one for
    //         the graph and one for the IDBCgraph (which handles the dofs
    //         fixed by inhomogeneous Dirichlet boundary conditions)
    for ( UInt i = 1; i <= length2; i++ ) {
      aux = std::abs( connect2[i] );
      if ( aux > 0 ) {
        if ( aux > numLastFreeDof_ ) {
          edgeList2_.push_back( aux - numLastFreeDof_ );
        }
        else {
          edgeList1_.push_back( aux );
        }
      }
    }

#ifdef DEBUG_GRAPHMANAGERSIMPLE2
    // output original connectivity
    (*debug) << "\n GraphManagerSimple::AdaptConnects\n"
             << " numLastFreeDof = " << numLastFreeDof_ << '\n';
    (*debug) << " connect1 ";
    for ( UInt i = 1; i <= length1; i++ ) {
      (*debug) << connect1[i] << " ";
    }
    (*debug) << std::endl;
    (*debug) << " connect2 ";
    for ( UInt i = 1; i <= length2; i++ ) {
      (*debug) << connect2[i] << " ";
    }
    (*debug) << std::endl;

    // output new connectivity
    (*debug) << " vertexList1: ";
    for ( UInt i = 0; i < vertexList1_.size(); i++ ) {
      (*debug) << vertexList1_[i] << " ";
    }
    (*debug) << std::endl;
    (*debug) << " vertexList2: ";
    for ( UInt i = 0; i < vertexList2_.size(); i++ ) {
      (*debug) << vertexList2_[i] << " ";
    }
    (*debug) << std::endl;
    (*debug) << " edgeList1: ";
    for ( UInt i = 0; i < edgeList1_.size(); i++ ) {
      (*debug) << edgeList1_[i] << " ";
    }
    (*debug) << std::endl;
    (*debug) << " edgeList2: ";
    for ( UInt i = 0; i < edgeList2_.size(); i++ ) {
      (*debug) << edgeList2_[i] << " ";
    }
    (*debug) << std::endl;
#endif

    // Check if we have any real dofs
    if ( vertexList1_.empty() == false ) {

      // Insert information into graph for real dofs
      graph_->AddVertexNeighbours( vertexList1_, edgeList1_ );

      // Check for assembly of counter part
      if( setCounterPart == true ) {

        // Insert reverse information into graph for real dofs
        graph_->AddVertexNeighbours( edgeList1_, vertexList1_ );
      }
      
      // Insert information into graph for fixed dofs
      if ( edgeList2_.empty() == false ) {
        graphIDBC_->AddVertexNeighbours( vertexList1_, edgeList2_ );

        if( setCounterPart == true ) {
          graphIDBC_->AddVertexNeighbours( edgeList1_, vertexList2_ );
        }
      }
    }
  }


  // =================
  //   GetReordering
  // =================
  Integer* GraphManagerSimple::GetReordering( const PdeIdType identifier ) {

    ENTER_FCN( "GraphManagerSimple::GetReordering" );

    Integer *retVal = NULL;

    // Test, whether we can return a re-ordering vector
    if ( reorderingDone_ == false ) {
      (*error) << "GraphManagerSimple::GetReordering: "
               << "No reordering vector available since the graph has not "
               << "been reordered, yet!";
      Error( __FILE__, __LINE__ );
    }
    else if ( newOrderingPassedOn_ == true ) {
      (*error) << "GraphManagerSimple::GetReordering: "
               << "No reordering vector available! The vector was already "
               << "queried by somebody else!";
      Error( __FILE__, __LINE__ );
    }

    // By passing the pointer to the array containing the re-ordering
    // information to the caller, this class forgets about the re-ordering
    retVal = newOrdering_;
    newOrdering_ = NULL;
    newOrderingPassedOn_ = true;

    return retVal;
  }


  // ====================
  //   CheckConsistency
  // ====================
  void GraphManagerSimple::CheckConsistency( const PdeIdType idPDE1,
                                             const PdeIdType idPDE2,
                                             std::string caller ) const {

    ENTER_FCN( "GraphManagerSimple::CheckConsistency" );

    // Check that first identifier is not nil
    if ( idPDE1 == NO_PDE_ID ) {
      (*error) << "GraphManagerSimple::"
               << caller << ": No PDE identifier was specified";
      Error( __FILE__, __LINE__ );
    }

    // Check that first identifier matches the one for this object
    if ( idPDE1 != myPDEIdent_ ) {
      (*error) << "GraphManagerSimple::"
               << caller << ": PDE identifier '" << idPDE1
               << "' does not match identifier '" << myPDEIdent_
               << "' which is managed by this object!";
      Error( __FILE__, __LINE__ );
    }

    // If second identifier is not nil, it must match first one
    if ( idPDE2 != NO_PDE_ID && idPDE2 != idPDE1 ) {
      (*error) << "GraphManagerSimple::"
               << caller << ": Received the two identifiers '"
	       << idPDE1 << "' and '" << idPDE2
               << "', but they do not match!";
      Error( __FILE__, __LINE__ );
    }

    // Check, if a graph is available
    if ( graph_ == NULL ) {
      (*error) << "GraphManagerSimple::"
               << caller << ": Pointer to graph object = NULL! "
	       << "Did you call RegisterPDE() once?";
      Error( __FILE__, __LINE__ );
    }
  }


  // ==========================
  //   BuildPermutationVector
  // ==========================
  void GraphManagerSimple::BuildPermutationVector() {

    ENTER_FCN( "GraphManagerSimple::BuildPermutationVector" );

    // CASE 1:
    //
    // This is the simple case of the penalty approach. In this case
    // we simply can pass newOrdering_ to CFS++ in the form graph_
    // returned it via FinaliseAssembly()
    if ( numLastFreeDof_ == numEqns_ ) {
      reorderingDone_ = true;
      return;
    }

    // CASE 2:
    //
    // If no re-ordering was performed by the graph_ object, then it is
    // also an easy case, since we simply can return the NULL pointer
    // to CFS++
    else if ( reorderType_ == NOREORDERING ) {
      if ( newOrdering_ != NULL ) {
        (*error) << "Internal error: Pointer should be NULL, but is not!";
        Error( __FILE__, __LINE__ );
      }
      reorderingDone_ = true;
      return;
    }

    // CASE 3:
    //
    // When the graph_ object has performed a re-ordering of the graph,
    // then we simply have to add an identity permutation for the fixed
    // dofs
    else {
      for ( UInt i = numLastFreeDof_ + 1; i <= numEqns_; i++ ) {
        newOrdering_[i] = i;
      }
      reorderingDone_ = true;
    }

    // Log new mapping to debug file
#ifdef DEBUG_GRAPHMANAGERSIMPLE2
    (*debug) << "\n GraphManagerSimple - new equation numbers:\n";
    if ( reorderType_ != NOREORDERING ) {
      for ( UInt i = 1; i <= numEqns_; i++ ) {
        (*debug) << i << " -> " << newOrdering_[i] << std::endl;
      }
    }
#endif

  }


  // =============
  //   SetupDone
  // =============
  void GraphManagerSimple::SetupDone() {

    ENTER_FCN( "GraphManagerSimple::SetupDone" );

    // Finish generation of primary graph
    graph_->FinaliseAssembly( newOrdering_ );

    // Finish generation of auxilliary graph
    if ( graphIDBC_ != NULL ) {
      graphIDBC_->FinaliseAssembly( newOrdering_ );
    }

    // Build full permutation vector including fixed dofs
    BuildPermutationVector();

    // Be verbose
    PrintStats( cla );
  }


  // ==============
  //   PrintStats
  // ==============
  void GraphManagerSimple::PrintStats( std::ostream *log ) {

    ENTER_FCN( "GraphManagerSimple::PrintStats" );

    // Determine width for pretty-printing
    UInt tw = numLastFreeDof_;
    tw = tw > 0 ? 47 + (UInt)std::log10( (float)tw ) + 1 : 48;

    // Begin report block
    (*log) << "\n " << std::setw(tw) << std::setfill( '=' ) << '=' << "\n"
           << std::setfill( ' ' )
           << " GRAPHMANAGER:\n\n";

    // Print standard graphmanager info
    (*log) << " Type of graph manager: GraphManagerSimple\n"
           << " Number of sub-graphs: 1\n"
           << " Number of attached PDEs: 1\n"
           << " PDE identifier: '" << myPDEIdent_ << "'"
           << "\n Total number of degrees of freedom: " << numEqns_
           << "\n Number of last free/unfixed degree of freedom: "
           << numLastFreeDof_ << '\n';

    // Close report block
    (*log) << ' ' << std::setw(tw) << std::setfill( '=' ) << '=' << "\n"
           << std::setfill( ' ' )  << std::endl;
  }

}
