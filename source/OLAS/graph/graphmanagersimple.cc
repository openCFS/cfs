// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

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

#include "Utils/tools.hh"
#include "OLAS/graph/graphmanagersimple.hh"

namespace CoupledField {

  //#define DEBUG_GRAPHMANAGERSIMPLE2

  // ===============
  //   Constructor
  // ===============
  GraphManagerSimple::GraphManagerSimple() :
    newOrderingPassedOn_(false),
    myPDEIdent_(NO_PDE_ID),
    graph_(NULL),
    graphIDBC_(NULL),
    reorderingDone_(false),
    numEqns_(0),
    numLastFreeDof_(0)
  {
  }


  // ==============
  //   Destructor
  // ==============
  GraphManagerSimple::~GraphManagerSimple() {

    // If no re-ordering was performed the pointer to the re-ordering vector
    // is still NULL. If it was claimed by the PDE it was re-set to NULL and
    // if it was not claimed, it's our responsibility to delete it now.
    if ( !newOrdering_.GetSize() ) {
      if ( newOrderingPassedOn_ == false ) {
        std::string tmp;
        Enum2String( reorderType_, tmp );
        
	(*warning) << "GraphManagerSimple: I was told to perform a '"
                   << tmp << "' re-ordering of the "
                   << "graph, but nobody ever claimed the permutation "
                   << "vector! Assuming it's my task to de-allocate the "
                   << "memory!";
	Warning( __FILE__, __LINE__ );
      }
      newOrdering_.Resize(0);
    }

    // Delete the graph objects
    delete graph_;
    delete graphIDBC_;

  }


  // =============
  //   SetupInit
  // =============
  void GraphManagerSimple::SetupInit( UInt numPDEs ) {


    // Make sure that there is only one PDE, since we are the graph manager
    // for the simple case
    if ( numPDEs != 1 ) {
      EXCEPTION("GraphManagerSimple::SetupInit: This is the graph manager "
               << "for the simple case of one single PDE! However SetupInit "
               << "was informed that there are " << numPDEs << " PDEs!");
    }
  }


  // ===============
  //   RegisterPDE
  // ===============
  void GraphManagerSimple::RegisterPDE( const FeFctIdType identifierPDE,  
                                        const UInt numEqns,
                                        const UInt numLastFreeDof,
					const ReorderingType reorder ) {


#ifdef DEBUG_GRAPHMANAGERSIMPLE1

    // Avoid mis-use
    if ( myPDEIdent_ != NO_PDE_ID ) {
      EXCEPTION("GraphManagerSimple: A PDE with identifier '" << myPDEIdent_
	       << "' has already been registered with this instance! Refusing "
	       << "to register PDE " << identifierPDE);
    }

#endif

    // Make a sensibility check to avoid problems in the long-run
    if ( numLastFreeDof == 0 ) {
      EXCEPTION("GraphManagerSimple::RegisterPDE: You tried to register "
               << "a PDE with identifier '" << myPDEIdent_ << "' that has "
               << "numEqns = " << numEqns << ", but numLastFreeDof = "
               << numLastFreeDof << ". If this is really intended, please "
               << "set <setup idbcHandling=\"penalty\"/> in your xml-file.");
    }

    // Store identifier
    myPDEIdent_ = identifierPDE;

    // Store sizes
    numEqns_        = numEqns;
    numLastFreeDof_ = numLastFreeDof;

    // Generate graph object
    graph_ = new BaseGraph( numLastFreeDof_, numLastFreeDof_, reorder );
    ASSERTMEM( graph_, sizeof(BaseGraph) );

    // Generate graph object for IDBC
    UInt graphSize = 0;
    if ( numLastFreeDof_ < numEqns_ ) {
      graphSize = numLastFreeDof_;
      graphIDBC_ = new IDBC_Graph( graphSize, numEqns_ - numLastFreeDof );
      ASSERTMEM( graphIDBC_, sizeof(IDBC_Graph) );

#ifdef DEBUG_GRAPHMANAGERSIMPLE2
      (*debug) << " Generated IDBC_Graph with " << graphSize << " vertices"
               << std::endl;
#endif
    }

    // Allocate memory for the re-ordering and store re-ordering type
    reorderType_ = reorder;
    if ( reorderType_ != NOREORDERING ) {
      newOrdering_.Resize( numEqns_ );
    }

  }


  // ================
  //   AssembleInit
  // ================
  void GraphManagerSimple::AssembleInit( const FeFctIdType identifierPDE1,
					 const FeFctIdType identifierPDE2,
                                         bool assemblingTranspose ) {


#ifdef DEBUG_GRAPHMANAGERSIMPLE1
    CheckConsistency( identifierPDE1, identifierPDE2, "AssembleInit" );
#endif

    // Nothing else to be done, since graph object itself does not
    // require a special call before the start of the assembly

  }


  // ================
  //   AssembleDone
  // ================
  void GraphManagerSimple::AssembleDone( const FeFctIdType identifierPDE1,
					 const FeFctIdType identifierPDE2,
                                         bool assemblingTranspose ) {


#ifdef DEBUG_GRAPHMANAGERSIMPLE1
    CheckConsistency( identifierPDE1, identifierPDE2, "AssembleDone" );
#endif

  }


  // ============
  //   GetGraph
  // ============
  BaseGraph* GraphManagerSimple::GetGraph( const FeFctIdType identifierPDE1,
					   const FeFctIdType identifierPDE2 ){


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
  BaseGraph* GraphManagerSimple::GetIDBCGraph( const FeFctIdType pdeID1,
                                               const FeFctIdType pdeID2 ) const{


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
  void GraphManagerSimple::SetElementPos( const FeFctIdType identifierPDE1,
                                          const StdVector<Integer>& eqnNrs1,
                                          const FeFctIdType identifierPDE2,
                                          const StdVector<Integer>& eqnNrs2,
                                          bool setCounterPart ) {

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
    for ( UInt i = 0, length1 = eqnNrs1.GetSize(); i < length1; i++ ) {
      aux = abs(eqnNrs1[i]);
      if ( aux > 0 ) {
        //since the graph will be 0-based, we substract from the 
        //equation number "aux" minus one
        if ( aux <= numLastFreeDof_  ) {
          vertexList1_.push_back( aux - 1 );
        } else {
          vertexList2_.push_back( aux - numLastFreeDof_ -1 );
        }
      }
    }

    // STEP 2: Split the second connect array into two edge lists, one for
    //         the graph and one for the IDBCgraph (which handles the dofs
    //         fixed by inhomogeneous Dirichlet boundary conditions)
    for ( UInt i = 0, length2 = eqnNrs2.GetSize(); i < length2; i++ ) {
      aux = abs(eqnNrs2[i]);
      if ( aux > 0 ) {
        //since the graph will be 0-based, we substract from the 
        //equation number "aux" minus one
        if ( aux > numLastFreeDof_ ) {
          edgeList2_.push_back( aux - numLastFreeDof_ - 1);
        }
        else {
          edgeList1_.push_back( aux - 1 );
        }
      }
    }

#ifdef DEBUG_GRAPHMANAGERSIMPLE2
    // output original connectivity
    (*debug) << "\n GraphManagerSimple::AdaptConnects\n"
             << " numLastFreeDof = " << numLastFreeDof_ << '\n';
    (*debug) << " connect1 ";
    for ( UInt i = 0; i < length1; i++ ) {
      (*debug) << connect1[i] << " ";
    }
    (*debug) << std::endl;
    (*debug) << " connect2 ";
    for ( UInt i = 0; i < length2; i++ ) {
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
  void GraphManagerSimple::GetReordering( const FeFctIdType identifier,
                                          StdVector<UInt>& order ) {

    // Test, whether we can return a re-ordering vector
    if ( reorderingDone_ == false ) {
      EXCEPTION("GraphManagerSimple::GetReordering: "
               << "No reordering vector available since the graph has not "
               << "been reordered, yet!");
    }
    else if ( newOrderingPassedOn_ == true ) {
      EXCEPTION("GraphManagerSimple::GetReordering: "
               << "No reordering vector available! The vector was already "
               << "queried by somebody else!");
    }

    // By passing the pointer to the array containing the re-ordering
    // information to the caller, this class forgets about the re-ordering
    order = newOrdering_;
    newOrdering_.Resize(0);
    newOrderingPassedOn_ = true;
  }


  // ====================
  //   CheckConsistency
  // ====================
  void GraphManagerSimple::CheckConsistency( const FeFctIdType idPDE1,
                                             const FeFctIdType idPDE2,
                                             std::string caller ) const {


    // Check that first identifier is not nil
    if ( idPDE1 == NO_PDE_ID ) {
      EXCEPTION("GraphManagerSimple::"
               << caller << ": No PDE identifier was specified");
    }

    // Check that first identifier matches the one for this object
    if ( idPDE1 != myPDEIdent_ ) {
      EXCEPTION("GraphManagerSimple::"
               << caller << ": PDE identifier '" << idPDE1
               << "' does not match identifier '" << myPDEIdent_
               << "' which is managed by this object!");
    }

    // If second identifier is not nil, it must match first one
    if ( idPDE2 != NO_PDE_ID && idPDE2 != idPDE1 ) {
      EXCEPTION("GraphManagerSimple::"
               << caller << ": Received the two identifiers '"
	       << idPDE1 << "' and '" << idPDE2
               << "', but they do not match!");
    }

    // Check, if a graph is available
    if ( graph_ == NULL ) {
      EXCEPTION("GraphManagerSimple::"
               << caller << ": Pointer to graph object = NULL! "
	       << "Did you call RegisterPDE() once?");
    }
  }


  // ==========================
  //   BuildPermutationVector
  // ==========================
  void GraphManagerSimple::BuildPermutationVector() {


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
      if ( newOrdering_.GetSize() ) {
        EXCEPTION("Internal error: Pointer should be NULL, but is not!");
      }
      reorderingDone_ = true;
      return;
    }

    // CASE 3:
    //
    // When the graph_ object has performed a re-ordering of the graph,
    // then we simply have to add an identity permutation for the fixed
    // dofs
    // note: the equation numbers are 1-based!!
    else {
      for ( UInt i = numLastFreeDof_; i < numEqns_; i++ ) {
        newOrdering_[i] = i+1;
      }
      reorderingDone_ = true;
    }

    // Log new mapping to debug file
#ifdef DEBUG_GRAPHMANAGERSIMPLE2
    (*debug) << "\n GraphManagerSimple - new equation numbers:\n";
    if ( reorderType_ != NOREORDERING ) {
      for ( UInt i = 0; i < numEqns_; i++ ) {
        (*debug) << i+1 << " -> " << newOrdering_[i] << std::endl;
      }
    }
#endif

  }


  // =============
  //   SetupDone
  // =============
  void GraphManagerSimple::SetupDone() {


    // Finish generation of primary graph
    graph_->FinaliseAssembly( &newOrdering_ );

    // Finish generation of auxilliary graph
    if ( graphIDBC_ != NULL ) {
      graphIDBC_->FinaliseAssembly( &newOrdering_ );
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
