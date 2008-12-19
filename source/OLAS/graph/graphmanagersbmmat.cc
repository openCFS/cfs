// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "graph/graphmanagersbmmat.hh"

#include <iterator>
#include <iomanip>

namespace OLAS {


  // ===============
  //   Constructor
  // ===============
  GraphManagerSBMMat::GraphManagerSBMMat() {


    graph_               = NULL;
    graphIDBC_           = NULL;
    isCoupled_           = NULL;
    numLastFreeDof_      = NULL;
    numEqn_              = NULL;
    newOrdering_         = NULL;

    reorderingDone_      = false;
    registrationDone_    = false;
    numPDEs_             = 0;
    numRegisteredPDEs_   = 0;
    numAssembledPDEs_    = 0;
  }


  // ==============
  //   Destructor
  // ==============
  GraphManagerSBMMat::~GraphManagerSBMMat() {


    // If no re-ordering was performed the pointers to the re-ordering
    // vectors are still NULL. If they were claimed by the PDEs they
    // were re-set to NULL and if they were not claimed, it is our
    // responsibility to delete then now.
    for ( UInt i = 1; i <= numPDEs_; i++ ) {
      if ( newOrdering_[i] != NULL ) {
        (*error) << "GraphManagerSBMMat::~GraphManagerSBMMat: "
                 << "Nobody ever claimed the permutation vector for the "
                 << "PDE with identifier '" << i
                 << "'! Assuming it's my task to de-allocate the memory!";
        Error( __FILE__, __LINE__ );
      }
      else {
        DeleteArray( newOrdering_[i] );
      }
    }
    DeleteArray( newOrdering_ );

    // Delete the graph objects
    for ( UInt i = 1; i <= numPDEs_ * numPDEs_; i++ ) {
      delete graph_[i];
    }
    DeleteArray( graph_ );

    // Delete the IDBC graph objects
    for ( UInt i = 1; i <= numPDEs_ * numPDEs_; i++ ) {
      delete graphIDBC_[i];
    }
    DeleteArray( graphIDBC_ );

    // Delete remaining arrays
    DeleteArray( numLastFreeDof_ );
    DeleteArray( numEqn_ );
    DeleteArray( isCoupled_ );
  }


  // =============
  //   SetupInit
  // =============
  void GraphManagerSBMMat::SetupInit( UInt numPDEs ) {


    // Now we now for how many PDEs we are responsible ...
    numPDEs_ = numPDEs;

    // ... and can build the empty graph pointer matrix
    NewArray( graph_, BaseGraph*, numPDEs_ * numPDEs_ );
    for ( UInt i = 1; i <= numPDEs_ * numPDEs_; i++ ) {
      graph_[i] = NULL;
    }

    // ... and the empty IDBC graph array
    NewArray( graphIDBC_, IDBC_Graph*, numPDEs_ * numPDEs_ );
    for ( UInt i = 1; i <= numPDEs_ * numPDEs_; i++ ) {
      graphIDBC_[i] = NULL;
    }

    // Allocate memory to store PDE specific information
    NewArray( numLastFreeDof_, UInt, numPDEs_ );
    NewArray( numEqn_        , UInt, numPDEs_ );
    for ( UInt i = 1; i <= numPDEs_; i++ ) {
      numLastFreeDof_[i] = 0;
      numEqn_[i] = 0;
    }

    // Allocate memory to store the permutation vectors for the reordering
    // of the unknowns of each PDE
    NewArray( newOrdering_, Integer*, numPDEs_ );
    for ( UInt i = 1; i <= numPDEs_; i++ ) {
      newOrdering_[i] = NULL;
    }

    // Setup empty array for coupling flags
    NewArray( isCoupled_, bool, numPDEs_ * numPDEs_ );
    for ( UInt i = 1; i <= numPDEs_ * numPDEs_; i++ ) {
      isCoupled_[i] = false;
    }

  }


  // =============
  //   SetupDone
  // =============
  void GraphManagerSBMMat::SetupDone() {


    // Print statistics to standard log stream
    PrintStats( cla );

  }


  // ===============
  //   RegisterPDE
  // ===============
  void GraphManagerSBMMat::RegisterPDE( const PdeIdType identifierPDE,
                                        const UInt numEqns,
                                        const UInt numLastFreeDof,
                                        const ReorderingType reorder ) {


    // Be cautious
    if ( registrationDone_ == true ) {
      (*error) << "Attempt to use RegisterPDE() after end of "
               << "registration phase, i.e. after the first call to "
               << "AssembleInit()!";
      Error( __FILE__, __LINE__ );
    }

    // Step counter for the number of registered PDEs and check number
    numRegisteredPDEs_++;
    if ( numRegisteredPDEs_ > numPDEs_ ) {
      (*error) << "GraphManagerSBMMat::RegisterPDE: You tried to "
               << "register a " << numRegisteredPDEs_ << "-th PDE "
               << "with id = '" << identifierPDE << "', but SetupInit "
               << "specified only " << numPDEs_ << " to be expected!";
      Error( __FILE__, __LINE__ );
    }
    if ( numRegisteredPDEs_ == NO_PDE_ID ) {
      (*error) << "GraphManagerSBMMat::RegisterPDE: You tried to "
               << "register a PDE, but identifier is NO_PDE_ID!";
      Error( __FILE__, __LINE__ );
    }

    // Store info on number of unknowns and euqation numbers of this PDE
    numLastFreeDof_[numRegisteredPDEs_] = numLastFreeDof;
    numEqn_[numRegisteredPDEs_] = numEqns;

    // Generate graph object for this PDE
    UInt idx = ComputeIndex( identifierPDE, NO_PDE_ID );
    graph_[idx] = New BaseGraph( numLastFreeDof, numLastFreeDof, reorder );
    if ( graph_[idx] == NULL ) {
      (*error) << "Generation of graph object for PDE with identifier "
               << identifierPDE << " failed!";
      Error( __FILE__, __LINE__ );
    }

#ifdef DEBUG_GRAPHMANAGERSBMMAT
    (*debug) << " GraphManagerSBMMat: Generated sub-graph ("
             << identifierPDE << ", " << identifierPDE << ")"
             << " for a " << numLastFreeDof << " x " << numLastFreeDof
             << " matrix" << std::endl;
#endif

    // If reordering is going to be performed for the current PDE then
    // we need to allocate memory to store the resulting permutation
    // vector
    if ( reorder != NOREORDERING ) {
      NewArray( newOrdering_[identifierPDE], Integer, numLastFreeDof );
    }
  }


  // ================
  //   AssembleInit
  // ================
  void GraphManagerSBMMat::AssembleInit( const PdeIdType idPDE1,
                                         const PdeIdType idPDE2,
                                         bool assemblingTranspose ) {


    // Perform a consistency check
    if ( idPDE1 == NO_PDE_ID || idPDE2 == NO_PDE_ID ) {
      (*error) << "GraphManagerSBMMat: Please specify two valid PDE "
               << "identifiers!";
      Error( __FILE__, __LINE__ );
    }

    // Generate coupling graph and also transpose if necessary
    if ( idPDE1 != idPDE2 ) {
      GenerateCouplingGraph( idPDE1, idPDE2 );
      if ( assemblingTranspose ) {
        GenerateCouplingGraph( idPDE2, idPDE1 );
      }
    }

    // Generate IDBC graph and its transpose if necessary
    GenerateIDBCGraph( idPDE1, idPDE2 );
    if ( assemblingTranspose ) {
      GenerateIDBCGraph( idPDE2, idPDE1 );
    }

    // Check that all PDEs have finished their assembly. This is important
    // since we do not want to do a reordering for the coupling objects!
    if ( idPDE2 != idPDE1 && numAssembledPDEs_ != numPDEs_ ) {
      (*error) << "GraphManagerSBMMat::AssembleInit: You are trying to "
               << "assemble a coupling object before assembling all PDEs! "
               << "Naughty you!";
      Error( __FILE__, __LINE__ );
    }
  }


  // ================
  //   AssembleDone
  // ================
  void GraphManagerSBMMat::AssembleDone( const PdeIdType idPDE1,
                                         const PdeIdType idPDE2,
                                         bool assemblingTranspose ) {


    // Compute index into graph pointer matrix
    UInt idx = ComputeIndex( idPDE1, idPDE2 );

    // Perform a consistency check
    if ( idPDE1 == NO_PDE_ID ) {
      (*error) << "GraphManagerSBMMat: First PDE identifier passed to "
               << "AssembleInit is empty!";
      Error( __FILE__, __LINE__ );
    }
    if ( graph_[idx] == NULL ) {
      (*error) << "GraphManagerSBMMat::AssembleInit: "
               << "Pointer to graph object = NULL! "
               << "Did you call RegisterPDE() for all " << numPDEs_
               << " PDEs?";
      Error( __FILE__, __LINE__ );
    }

#ifdef DEBUG_GRAPHMANAGERSBMMAT
      (*debug) << " GraphManagerSBMMat:\n"
               << " Finalising assembly of sub-graph ("
               << idPDE1 << ", " << idPDE2 << ")"
               << " with index " << idx << std::endl;
      if ( assemblingTranspose == true ) {
        (*debug) << " and of sub-graph ("
                 << idPDE2 << ", " << idPDE1 << ")"
                 << " with index " << ComputeIndex( idPDE2, idPDE1 )
                 << std::endl;
      }
      (*debug) << '\n';
#endif

    // If we assembled a coupling object, then we should not pass memory
    // for storing a permutation vector
    if ( idPDE2 != idPDE1 && idPDE2 != NO_PDE_ID ) {
      graph_[idx]->FinaliseAssembly( NULL );
      if ( assemblingTranspose == true ) {
        graph_[ComputeIndex(idPDE2, idPDE1)]->FinaliseAssembly( NULL );
      }
    }
    else {

      // Finalise assembly of graph
      graph_[idx]->FinaliseAssembly( newOrdering_[idPDE1] );

      // Increase counter of PDEs that have finalised their assembly
      numAssembledPDEs_++;
    }

    // Now finalise the assembly of the associated IDBC graph
    if ( graphIDBC_[idx] != NULL ) {
      graphIDBC_[idx]->FinaliseAssembly( newOrdering_[idPDE1] );

#ifdef DEBUG_GRAPHMANAGERSBMMAT
      (*debug) << " GraphManagerSBMMat:\n"
               << " Finalising assembly of IDBC graph ("
               << idPDE1 << ", " << idPDE2 << ")"
               << " with index " << idx << std::endl;
#endif
    }

    // Now finalise assembly of transpose IDBC graph
    idx = ComputeIndex( idPDE2, idPDE1 );
    if ( assemblingTranspose == true && graphIDBC_[idx] != NULL ) {
      graphIDBC_[idx]->FinaliseAssembly( newOrdering_[idPDE2] );

#ifdef DEBUG_GRAPHMANAGERSBMMAT
      (*debug) << " GraphManagerSBMMat:\n"
               << " Finalising assembly of IDBC graph ("
               << idPDE2 << ", " << idPDE1 << ")"
               << " with index " << idx << std::endl;
#endif
    }

    // Check, if all PDEs were assembled. If this is the case, then
    // all re-orderings are now available
    if ( numAssembledPDEs_ == numPDEs_ ) {
      reorderingDone_ = true;
    }
  }


  // =================
  //   SetElementPos
  // =================
  void GraphManagerSBMMat::SetElementPos( const PdeIdType pdeID1,
                                          Integer *connect1,
                                          Integer length1,
                                          const PdeIdType pdeID2,
                                          Integer *connect2,
                                          Integer length2,
                                          bool setCounterPart ) {


    // Compute index of graph in graph pointer matrix
    UInt idx = ComputeIndex( pdeID1, pdeID2 );

#ifdef DEBUG_GRAPHMANAGERSTDMAT

    // Perform some consistency checks in debug mode only
    if ( pdeID1 == NO_PDE_ID ) {
      (*error) << "GraphManagerSBMMat: First PDE identifier passed to "
               << "SetElementPos is empty!";
      Error( __FILE__, __LINE__ );
    }
    else if ( pdeID1 > numPDEs_ ) {
      (*error) << "GraphManagerSBMMat: First PDE identifier passed to "
               << "SetElementPos is " << pdeID1
               << " which exceeds the number of " << numPDEs_
               << " registered PDEs!";
      Error( __FILE__, __LINE__ );
    }
    if ( graph_[idx] == NULL ) {
      (*error) << "GraphManagerSBMMat::SetElementPos: "
               << "Pointer to graph object = NULL! "
               << "Did you call RegisterPDE() for all " << numPDEs_
               << " PDEs?";
      Error( __FILE__, __LINE__ );
    }

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
    for ( int i = 1; i <= length1; i++ ) {
      aux = std::abs( (double) connect1[i] );
      if ( aux > 0 ) {
        if ( aux <= numLastFreeDof_[pdeID1] ) {
          vertexList1_.push_back( aux );
        }
        else {
          vertexList2_.push_back( aux - numLastFreeDof_[pdeID1] );
        }
      }
    }

    // STEP 2: Split the second connect array into two edge lists, one for
    //         the graph and one for the IDBCgraph (which handles the dofs
    //         fixed by inhomogeneous Dirichlet boundary conditions)
    for ( int i = 1; i <= length2; i++ ) {
      aux = std::abs( (double) connect2[i] );
      if ( aux > 0 ) {
        if ( aux > numLastFreeDof_[pdeID2] ) {
          edgeList2_.push_back( aux - numLastFreeDof_[pdeID2] );
        }
        else {
          edgeList1_.push_back( aux );
        }
      }
    }

#ifdef DEBUG_GRAPHMANAGERSBMMAT
    // output original connectivity
    (*debug) << "\n GraphManagerSBMMat::AdaptConnects\n"
             << " idPDE1 = " << pdeID1 << '\n'
             << " idPDE2 = " << pdeID2 << '\n'
             << " numLastFreeDof1 = " << numLastFreeDof_[pdeID1] << '\n'
             << " numLastFreeDof2 = " << numLastFreeDof_[pdeID2] << '\n';
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

    idx = ComputeIndex( pdeID1, pdeID2 );

    // Insert information into graph for real dofs
    graph_[idx]->AddVertexNeighbours( vertexList1_, edgeList1_ );

    // Insert information into graph for fixed dofs
    graphIDBC_[idx]->AddVertexNeighbours( vertexList1_, edgeList2_ );

    // Check for assembly of counter part
    if ( (pdeID1 != pdeID2) && setCounterPart == true ) {

      idx = ComputeIndex( pdeID2, pdeID1 );

      // Insert information into (transpose) graph for real dofs
      graph_[idx]->AddVertexNeighbours( edgeList1_, vertexList1_);

      // Insert information into graph for fixed dofs
      graphIDBC_[idx]->AddVertexNeighbours( edgeList1_, vertexList2_ );
    }
  }


  // =================
  //   GetReordering
  // =================
  Integer* GraphManagerSBMMat::GetReordering( const PdeIdType identifier ) {


    Integer *retVal = NULL;

    // Small consistency check
    if ( identifier > numPDEs_ ) {
      (*error) << "GraphManagerSBMMat::GetReordering: "
               << "PDE with identifier '" << identifier << "' was not "
               << "registered using RegisterPDE()!";
      Error( __FILE__, __LINE__ );
    }

    // Test, whether we can return a re-ordering vector
    if ( reorderingDone_ == false ) {
      (*error) << "GraphManagerSBMMat::GetReordering: "
               << "No reordering vector available since the graphs have not "
               << "been reordered, yet!";
      Error( __FILE__, __LINE__ );
    }

    // By passing the pointer to the array containing the re-ordering
    // information to the caller, this class forgets about the re-ordering
    retVal = newOrdering_[identifier];
    newOrdering_[identifier] = NULL;
    return retVal;
  }


  // ============
  //   GetGraph
  // ============
  BaseGraph* GraphManagerSBMMat::GetGraph( const PdeIdType identifierPDE1,
                                           const PdeIdType identifierPDE2 ) {


    // Check consisteny
    if ( identifierPDE1 == NO_PDE_ID ) {
      (*error) << "GraphManagerSBMMat::GetGraph: Please specify at least one "
               << "PDE identifier!";
      Error( __FILE__, __LINE__ );
    }

    // Determine which graph we are looking for
    UInt idx = ComputeIndex( identifierPDE1, identifierPDE2 );

    // Check if a graph object was already created
    if ( graph_[idx] == NULL) {
      (*error) << "GraphManagerSBMMat: There exists no graph "
               << "for the PDE identifier pair ( " << identifierPDE1
               << " , " << identifierPDE2
               << ") which could be returned by the GetGraph() method.";
      Error( __FILE__, __LINE__);
    }

    // Return pointer to the graph object
    return graph_[idx];

  }


  // ================
  //   GetIDBCGraph
  // ================
  BaseGraph* GraphManagerSBMMat::GetIDBCGraph( const PdeIdType pdeID1,
                                               const PdeIdType pdeID2 ) const{


    // Check consisteny
    if ( pdeID1 == NO_PDE_ID || pdeID2 == NO_PDE_ID ) {
      (*error) << "GraphManagerSBMMat::GetIDBCGraph: "
               << "Please specify two valid PDE identifiers!";
      Error( __FILE__, __LINE__ );
    }

    UInt idx = ComputeIndex( pdeID1, pdeID2 );

    if ( graphIDBC_[idx] == NULL ) {
      (*error) << "GraphManagerSBMMat::GetIDBCGraph: "
               << "An IDBC graph object for PDE pair (" << pdeID1
               << " , " << pdeID2 << ") does not or not yet exist!";
      Error( __FILE__, __LINE__ );
    }

    // Return pointer to the IDBC graph object
    return graphIDBC_[idx];
  }


  // ==================
  //   SubGraphExists
  // ==================
  bool GraphManagerSBMMat::SubGraphExists( const PdeIdType idPDE1,
                                           const PdeIdType idPDE2 ) const {
    return graph_[ ComputeIndex( idPDE1, idPDE2 ) ] != NULL;
  }


  // ===================
  //   IDBCGraphExists
  // ===================
  bool GraphManagerSBMMat::IDBCGraphExists( const PdeIdType idPDE1,
                                            const PdeIdType idPDE2 ) const {
    return graphIDBC_[ ComputeIndex( idPDE1, idPDE2 ) ] != NULL;
  }


  // ==============
  //   PrintStats
  // ==============
  void GraphManagerSBMMat::PrintStats( std::ostream *log ) {



    // ***********************************
    //  Assemble info for pretty-printing
    // ***********************************

    // Compute maximal column widths (PDE table)
    UInt cw1 = 0, cw2 = 0, cw3 = 0, cw4 = 0, tw = 0, aux;
    for ( UInt i = 1; i <= numPDEs_; i++) {
      aux = numEqn_[i] > 0 ? (UInt)std::log10( (float)numEqn_[i] ) + 1 : 1;
      cw1 = cw1 < aux ? aux : cw1;

      aux = numLastFreeDof_[i] > 0 ?
        (UInt)std::log10( (float)numLastFreeDof_[i] ) + 1 : 1;
      cw2 = cw2 < aux ? aux : cw2;
    }

    // Compute maximal column widths (sub-graph table)
    UInt idx = 0;
    for ( UInt i = 1; i <= numPDEs_; i++ ) {
      for ( UInt j = 1; j <= numPDEs_; j++ ) {
        idx = ComputeIndex(i,j);
        if ( graph_[idx] != NULL ) {
          aux = graph_[idx]->GetSize() > 0 ?
            (UInt)std::log10( (float)graph_[idx]->GetSize() ) + 1 : 1;
          cw3 = cw3 < aux ? aux : cw3;
          aux = graph_[idx]->GetNNE() > 0 ?
            (UInt)std::log10( (float)graph_[idx]->GetNNE() ) + 1 : 1;
          cw4 = cw4 < aux ? aux : cw4;
        }
      }
    }

    // Correct field widths for column headers
    cw1 = cw1 >  6 ? cw1 :  6;
    cw2 = cw2 > 11 ? cw2 : 11;
    cw3 = cw3 >  9 ? cw3 :  9;
    cw4 = cw4 >  6 ? cw4 :  6;

    // Set total width
    tw = 7 + 12 + cw1 + cw2 + 5;


    // *******************
    //  Report statistics
    // *******************

    // Begin report block
    (*log) << "\n " << std::setw(tw) << std::setfill( '=' ) << '=' << "\n"
           << std::setfill( ' ' )
           << " GRAPHMANAGER:\n\n";

    // Determine number of sub-graphs we are holding
    UInt numGraphs = 0;
    for ( UInt i = 1; i <= numPDEs_ * numPDEs_; i++ ) {
      if ( graph_[i] != NULL ) {
        numGraphs++;
      }
    }

    // Print standard graphmanager info
    (*log) << " Type of graph manager: GraphManagerSBMMat\n"
           << " Number of attached PDEs: " << numPDEs_
           << std::endl;

    // Output table showing identifier and number of unknowns
    (*log) << ' ' << std::setw(tw) << std::setfill( '-' ) << '-' << "\n"
           << "  PDE ID "
           << "| numEqn | lastFreeDof\n"
           << ' ' << std::setw(tw) << std::setfill( '-' ) << '-' << "\n"
           << std::setfill( ' ' );

    // now the table rows
    log->setf( std::ios::right, std::ios::adjustfield );
    for ( UInt i = 1; i <= numPDEs_; i++) {
      (*log) << std::setw(8) << i
             << " | " << std::setw(cw1) << numEqn_[i]
             << " | " << std::setw(cw2) << numLastFreeDof_[i]
             << std::endl;
    }
    (*log) << ' ' << std::setw(tw) << std::setfill( '-' ) << '-' << "\n"
           << std::setfill( ' ' );

    // Output table showing sub-graph information
    (*log) << " Number of sub-graphs: " << numGraphs << std::endl;
    (*log) << ' ' << std::setw(tw) << std::setfill( '-' ) << '-' << "\n"
           << "  row | col | #vertices | #edges\n"
           << ' ' << std::setw(tw) << std::setfill( '-' ) << '-' << "\n"
           << std::setfill( ' ' );

    idx = 0;
    for ( UInt i = 1; i <= numPDEs_; i++ ) {
      for ( UInt j = 1; j <= numPDEs_; j++ ) {
        idx = ComputeIndex(i,j);
        if ( graph_[idx] != NULL ) {
          (*log) << std::setw(5) << i << " | "
                 << std::setw(3) << j << " | "
                 << std::setw(cw3) << graph_[idx]->GetSize() << " | "
                 << std::setw(cw4) << graph_[idx]->GetNNE()
                 << std::endl;
        }
      }
    }


    // Close report block
    (*log) << ' ' << std::setw(tw) << std::setfill( '=' ) << '=' << "\n"
           << std::setfill( ' ' )  << std::endl;
  }


  // =========================
  //   GenerateCouplingGraph
  // =========================
  void GraphManagerSBMMat::GenerateCouplingGraph( PdeIdType pdeID1,
                                                  PdeIdType pdeID2 ) {


    UInt idx = ComputeIndex( pdeID1, pdeID2 );

    // Safety check
    if ( graph_[idx] != NULL ) {
      (*error) << "GraphManagerSBMMat::GenerateCouplingGraph:\n "
               << "A graph object for PDE pair (" << pdeID1
               << " , " << pdeID2 << ") already exists!";
      Error( __FILE__, __LINE__ );
    }

    // Generate graph object
    graph_[idx] = New BaseGraph( numLastFreeDof_[pdeID1],
                                 numLastFreeDof_[pdeID2], NOREORDERING );

    if ( graph_[idx] == NULL ) {
      (*error) << " GraphManagerSBMMat: Generation of sub-graph "
               << "for PDE pair (" << pdeID1 << " , " << pdeID2
               << ") and a " << numLastFreeDof_[pdeID1]
               << " x " << numLastFreeDof_[pdeID2] << " matrix failed ";
      Error( __FILE__, __LINE__ );
    }

#ifdef DEBUG_GRAPHMANAGERSBMMAT
    (*debug) << " GraphManagerSBMMat: Generated sub-graph for PDE pair ("
             << pdeID1 << " , " << pdeID2 << ") and a "
             << numLastFreeDof_[pdeID1] << " x " << numLastFreeDof_[pdeID2]
             << " matrix " << std::endl;
#endif
  }


  // =====================
  //   GenerateIDBCGraph
  // =====================
  void GraphManagerSBMMat::GenerateIDBCGraph( PdeIdType pdeID1,
                                              PdeIdType pdeID2 ) {


    UInt idx = ComputeIndex( pdeID1, pdeID2 );

    // Compute number of fixed dofs in second PDE of pair
    UInt fixedDofs = numEqn_[pdeID2] - numLastFreeDof_[pdeID2];

    // If there are fixed dofs in one of the PDEs we need an IDBC
    // graph object for this pair
    if ( fixedDofs > 0 ) {

      // Safety check
      if ( graphIDBC_[idx] != NULL ) {
        (*error) << "GraphManagerSBMMat::GenerateIDBCGraph:\n "
                 << "An IDBC graph object for PDE pair (" << pdeID1
                 << " , " << pdeID2 << ") already exists!";
        Error( __FILE__, __LINE__ );
      }

      // Generate IDBC graph object
      graphIDBC_[idx] = New IDBC_Graph( numLastFreeDof_[pdeID1], fixedDofs );
      if ( graphIDBC_[idx] == NULL ) {
        (*error) << " GraphManagerSBMMat: Generation of IDBC sub-graph "
                 << "for PDE pair (" << pdeID1 << " , " << pdeID2
                 << ") and a " << numLastFreeDof_[pdeID1]
                 << " x " << fixedDofs << " matrix failed ";
        Error( __FILE__, __LINE__ );
      }
    }

#ifdef DEBUG_GRAPHMANAGERSBMMAT
    if ( fixedDofs > 0 ) {
      (*debug) << " GraphManagerSBMMat: Generated IDBC sub-graph for"
               << " PDE pair (" << pdeID1 << " , " << pdeID2
               << ") and a " << numLastFreeDof_[pdeID1]
               << " x " << fixedDofs << " matrix " << std::endl;
    }
    else {
      (*debug) << " GraphManagerSBMMat: PDE pair (" << pdeID1 << " , "
               << pdeID2 << ") needs no IDBC graph" << std::endl;
    }
#endif

  }

}
