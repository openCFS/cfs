// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "graph/graphmanagerstdmat.hh"

#include <iterator>
#include <iomanip>

namespace OLAS {


  // ===============
  //   Constructor
  // ===============
  GraphManagerStdMat::GraphManagerStdMat() {


    graph_               = NULL;
    graphIDBC_           = NULL;
    numEqn_              = NULL;
    numLastFreeDof_      = NULL;
    newOrdering_         = NULL;
    offsetIDBC_          = NULL;
    offset_              = NULL;

    reorderType_         = NOREORDERING;

    reorderingDone_      = false;

    numPDEs_             = 0;
    numRegisteredPDEs_   = 0;
    numFreeDofs_         = 0;
    numFixedDofs_        = 0;
    offset1_             = 0;
    offset2_             = 0;

  }


  // ==============
  //   Destructor
  // ==============
  GraphManagerStdMat::~GraphManagerStdMat() {


    // If no re-ordering was performed the pointers to the re-ordering
    // vectors are still NULL. If they were claimed by the PDEs they
    // were re-set to NULL and if they were not claimed, it is our
    // responsibility to delete then now.
    for ( UInt i = 1; i <= numPDEs_; i++ ) {
      if ( newOrdering_[i] != NULL ) {
	(*error) << "GraphManagerStdMat::~GraphManagerStdMat: "
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
    delete graph_;
    delete graphIDBC_;

    // Delete auxilliary information arrays
    DeleteArray( numLastFreeDof_ );
    DeleteArray( offsetIDBC_ );
    DeleteArray( offset_ );
    DeleteArray( numEqn_ );

  }


  // =============
  //   SetupInit
  // =============
  void GraphManagerStdMat::SetupInit( UInt numPDEs ) {


    // Now we now for how many PDEs we are repsonsible
    numPDEs_ = numPDEs;

    // Allocate memory to store the offsets and further PDE size information
    NewArray( numEqn_        , UInt, numPDEs_     );
    NewArray( numLastFreeDof_, UInt, numPDEs_     );
    NewArray( offset_        , UInt, numPDEs_ + 1 );
    NewArray( offsetIDBC_    , Integer, numPDEs_ + 1 );

    // Allocate memory for the re-ordering vectors and set all
    // pointers to NULL for the case that no re-ordering is performed
    NewArray( newOrdering_, Integer*, numPDEs_ );
    for ( UInt i = 1; i <= numPDEs; i++ ) {
      newOrdering_[i] = NULL;
    }

  }


  // =============
  //   SetupDone
  // =============
  void GraphManagerStdMat::SetupDone() {



    // *******************
    //  Memory allocation
    // *******************

    // If re-ordering is to be performed, allocate auxilliary
    // vector to store large permutation vector
    Integer *permutation = NULL;
    if ( reorderType_ != NOREORDERING ) {
      NewArray( permutation, Integer, numFreeDofs_ );
    }

    // Only need permutation vector for first PDE, if we are
    // using the elimination approach and there are fixed dofs
    // for the first PDE
    UInt firstPermVec = 1;
    if ( numLastFreeDof_[1] == numEqn_[1] && reorderType_ == NOREORDERING ) {
      firstPermVec = 2;
    }

    // Each PDE, besides maybe the first, needs an individual reordering
    // array, so allocate memory for them taking into accout fixed dofs, too.
    for ( UInt i = firstPermVec; i <= numPDEs_; i++ ) {
      NewArray( newOrdering_[i], Integer, numEqn_[i] );
    }


    // *****************
    //  Re-order graphs
    // *****************

    // Trigger final steps of assembly
    graph_->FinaliseAssembly( permutation );

    // Finish generation of auxilliary graph
    if ( graphIDBC_ != NULL ) {
      graphIDBC_->FinaliseAssembly( permutation );
    }


    // ******************************
    //  Generate permutation vectors
    // ******************************

    // If re-ordering was performed, split large permutation vector
    // up into small vectors for each individual PDE
    UInt i, k, index;
    if ( reorderType_ != NOREORDERING ) {

      // Loop over all registered PDEs
      for ( i = 1; i <= numPDEs_; i++ ) {

        // Copy new dof number from large vector into PDE individual vector
        index = 1;
        for ( k = offset_[i] + 1; k <= offset_[i+1]; k++ ) {
          newOrdering_[i][index] = permutation[k];
          index++;
        }
      }

      // Large vector is no longer needed
      DeleteArray( permutation );
    }

    // No "real" re-ordering was performed. However, CFS++ needs the new
    // equation numbers, that include the offsets for all but the first PDE.
    // Therefore we generate a simple reordering for the other PDEs
    else {

      // Loop over all registered PDEs
      for ( i = firstPermVec; i <= numPDEs_; i++ ) {
        index = 1;
        for ( k = offset_[i] + 1; k <= offset_[i+1]; k++ ) {
          newOrdering_[i][index] = k;
          index++;
        }
      }
    }

    // In the elimination case, we must add the new equation numbers for
    // the fixed dofs which include an offsetto the end of the permutation
    // vector
    if ( numFixedDofs_ > 0 ) {

      // We use a running index for all fixed dofs per PDE
      UInt fixedBlockInit = numFreeDofs_;
      if ( firstPermVec == 2 ) {
        fixedBlockInit += numEqn_[1] - numLastFreeDof_[1];
      }

      for ( i = firstPermVec; i <= numPDEs_; i++ ) {
        for ( UInt j = 1; j <= numEqn_[i] - numLastFreeDof_[i]; j++ ) {
          newOrdering_[i][ numLastFreeDof_[i] + j ] = fixedBlockInit + j;
        }
        fixedBlockInit += numEqn_[i] - numLastFreeDof_[i];
      }
    }

    // If reordering was demanded, newOrdering_ now contains
    // the respective permutation vector
    reorderingDone_ = true;

    // Print statistics to standard log stream
    PrintStats( cla );

  }


  // ===============
  //   RegisterPDE
  // ===============
  void GraphManagerStdMat::RegisterPDE( const PdeIdType identifierPDE,
					const UInt n,
                                        const UInt numLastFreeDof,
					const ReorderingType reorder ) {


    // Step counter for the number of registered PDEs and check number
    numRegisteredPDEs_++;
    if ( numRegisteredPDEs_ > numPDEs_ ) {
      (*error) << "GraphManagerStdMat::RegisterPDE: You tried to "
               << "register a " << numRegisteredPDEs_ << "-th PDE "
               << "with id = '" << identifierPDE << "', but SetupInit "
               << "specified only " << numPDEs_ << " to be expected!";
      Error( __FILE__, __LINE__ );
    }


    // **************
    //  Registration
    // **************

    // Store total size and number of last real dof
    numEqn_[identifierPDE] = n;
    numLastFreeDof_[identifierPDE] = numLastFreeDof;

    // Increase counters for dof types
    numFreeDofs_  += numLastFreeDof_[identifierPDE];
    numFixedDofs_ += numEqn_[identifierPDE] - numLastFreeDof_[identifierPDE];


    // ************
    //  Reordering
    // ************

    // We can only do one type of re-ordering for one graph, so we
    // adopt the policy to chose the re-ordering specified by the
    // first PDE that registers itself.
    if ( numRegisteredPDEs_ == 1 ) {
      reorderType_ = reorder;
    }


    // ********************************************
    //  Stuff to be done after end of registration
    // ********************************************

    if ( numRegisteredPDEs_ == numPDEs_ ) {

      // Compute the offsets for the free equation numbers to be inserted
      // into the main graph and the offsets for the fixed equation numbers
      // that will be inserted into the auxilliary IDBC graph
      offset_[1] = 0;
      for ( UInt i = 1; i <= numPDEs_; i++ ) {
        offset_[i+1] = offset_[i] + numLastFreeDof_[i];
      }

      UInt aux = 0;
      for ( UInt i = 1; i <= numPDEs_; i++ ) {
        offsetIDBC_[i] = aux - numLastFreeDof_[i];
        aux += numEqn_[i] - numLastFreeDof_[i];
      }

      // Once the last PDE has registered itself, we know the number of
      // vertices in the graph and can thus construct the graph object.
      graph_ = New BaseGraph( numFreeDofs_, numFreeDofs_, reorderType_ );
      AssertMem( graph_, sizeof(BaseGraph) );

      // Generate auxilliary graph object
      if ( numFixedDofs_ > 0 ) {
        graphIDBC_ = New IDBC_Graph( numFreeDofs_, numFixedDofs_ );
        AssertMem( graphIDBC_, sizeof(IDBC_Graph) );
      }
    }
  }


  // ================
  //   AssembleInit
  // ================
  void GraphManagerStdMat::AssembleInit( const PdeIdType identifierPDE1,
					 const PdeIdType identifierPDE2,
                                         bool assemblingTranspose ) {


    // Perform a consistency check
    if ( identifierPDE1 == NO_PDE_ID ) {
      (*error) << "GraphManagerStdMat: First PDE identifier passed to "
               << "AssembleInit is empty!";
      Error( __FILE__, __LINE__ );
    }
    if ( graph_ == NULL ) {
      (*error) << "GraphManagerStdMat::AssembleInit: "
               << "Pointer to graph object = NULL! "
	       << "Did you call RegisterPDE() for all " << numPDEs_
               << " PDEs?";
      Error( __FILE__, __LINE__ );
    }

    // Assign identifiers for checking in SetElementPos
    idPDE1_ = identifierPDE1;
    idPDE2_ = identifierPDE2;

    // Determine offsets
    if ( idPDE1_ > numPDEs_ ) {
      (*error) << "GraphManagerStdMat::AssembleInit: "
               << "PDE with identifier '" << idPDE1_ << "' was not "
               << "registered using RegisterPDE()!";
      Error( __FILE__, __LINE__ );
    }
    offset1_ = offset_[idPDE1_];
    offsetIDBC1_ = offsetIDBC_[idPDE1_];
    
    if ( idPDE2_ != NO_PDE_ID ) {
      if ( idPDE2_ > numPDEs_ ) {
        (*error) << "GraphManagerStdMat::AssembleInit: "
                 << "PDE with identifier '" << idPDE2_ << "' was not "
                 << "registered using RegisterPDE()!";
        Error( __FILE__, __LINE__ );
      }
      offset2_ = offset_[idPDE2_];
      offsetIDBC2_ = offsetIDBC_[idPDE2_];
    }

  }


  // =================
  //   SetElementPos
  // =================
  void GraphManagerStdMat::SetElementPos( const PdeIdType idPDE1,
					  Integer *connect1,
					  Integer length1,
					  const PdeIdType idPDE2,
					  Integer *connect2,
					  Integer length2,
                                          bool setCounterPart ) {


#ifdef DEBUG_GRAPHMANAGERSTDMAT

    if ( idPDE1 == NO_PDE_ID ) {
      (*error) << "GraphManagerStdMat: First PDE identifier passed to "
               << "SetElementPos is empty!";
      Error( __FILE__, __LINE__ );
    }
    if ( graph_ == NULL ) {
      (*error) << "GraphManagerStdMat::SetElementPos: "
               << "Pointer to graph object = NULL! "
	       << "Did you call RegisterPDE() for all " << numPDEs_
               << " PDEs?";
      Error( __FILE__, __LINE__ );
    }
    if ( idPDE1_ != idPDE1 || idPDE2_ != idPDE2 ) {
      (*error) << "Consistency error detected! You specified to AssembleInit "
               << " idPDE1 = '" << idPDE1_ << "',  idPDE2 = '" <<  idPDE2_
               << "', but SetElementPos received idPDE1 = '"
               << idPDE1 << "',  idPDE2 = '" <<  idPDE2
               << "'. Maybe you forgot to call AssembleDone & AssembleInit?";
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
    for ( UInt i = 1; i <= length1; i++ ) {
      aux = std::abs( (double) connect1[i] );
      if ( aux > 0 ) {
        if ( aux <= numLastFreeDof_[idPDE1] ) {
          vertexList1_.push_back( aux + offset1_ );
        }
        else {
          vertexList2_.push_back( aux + offsetIDBC1_ );
        }
      }
    }

    // STEP 2: Split the second connect array into two edge lists, one for
    //         the graph and one for the IDBCgraph (which handles the dofs
    //         fixed by inhomogeneous Dirichlet boundary conditions)
    for ( UInt i = 1; i <= length2; i++ ) {
      aux = std::abs( (double) connect2[i] );
      if ( aux > 0 ) {
        if ( aux > numLastFreeDof_[idPDE2] ) {
          edgeList2_.push_back( aux + offsetIDBC2_ );
        }
        else {
          edgeList1_.push_back( aux + offset2_ );
        }
      }
    }

#ifdef DEBUG_GRAPHMANAGERSTDMAT
    // output original connectivity
    (*debug) << "\n GraphManagerStdMat::AdaptConnects\n"
             << " idPDE1 = " << idPDE1 << '\n'
             << " idPDE2 = " << idPDE2 << '\n'
             << " offset1 = " << offset1_ << '\n'
             << " offset2 = " << offset2_ << '\n'
             << " offsetIDBC1 = " << offsetIDBC1_ << '\n'
             << " offsetIDBC2 = " << offsetIDBC2_ << '\n'
             << " numLastFreeDof1 = " << numLastFreeDof_[idPDE1] << '\n'
             << " numLastFreeDof2 = " << numLastFreeDof_[idPDE2] << '\n';
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

    // Insert information into graph for real dofs
    graph_->AddVertexNeighbours( vertexList1_, edgeList1_ );

    // Insert information into graph for fixed dofs
    if ( graphIDBC_ != NULL ) {
      graphIDBC_->AddVertexNeighbours( vertexList1_, edgeList2_ );
    }

    // Check for assembly of counter part
    if ( setCounterPart == true ) {

      // Insert information into graph for real dofs
      graph_->AddVertexNeighbours( edgeList1_, vertexList1_ );

      // Insert information into graph for fixed dofs
      if ( graphIDBC_ != NULL ) {
        graphIDBC_->AddVertexNeighbours( edgeList1_, vertexList2_ );
      }
    }
  }


  // =================
  //   GetReordering
  // =================
  Integer* GraphManagerStdMat::GetReordering( const PdeIdType identifier ) {


    Integer *retVal = NULL;

    // Get registration number of PDE
    if ( identifier > numPDEs_ ) {
      (*error) << "GraphManagerStdMat::GetReordering: "
               << "PDE with identifier '" << identifier << "' was not "
               << "registered using RegisterPDE()!";
      Error( __FILE__, __LINE__ );
    }

    // Test, whether we can return a re-ordering vector
    if ( reorderingDone_ == false ) {
      (*error) << "GraphManagerStdMat::GetReordering: "
               << "No reordering vector available since the graph has not "
               << "been reordered, yet!";
      Error( __FILE__, __LINE__ );
    }
    else if ( newOrdering_[identifier] == NULL &&
              reorderType_ != NOREORDERING ) {
      (*error) << "GraphManagerStdMat::GetReordering: "
               << "No reordering vector available for PDE with identifier '"
               << identifier << "' Maybe it was already claimed?";
      Error( __FILE__, __LINE__ );
    }

    // By passing the pointer to the array containing the re-ordering
    // information to the caller, this class forgets about the re-ordering
    retVal = newOrdering_[identifier];
    newOrdering_[identifier] = NULL;

#ifdef DEBUG_GRAPHMANAGERSTDMAT
    (*debug) << "\nGraphManagerStdMat::GetReordering:\n"
             << "Re-ordering for PDE with ID = '" << identifier
             << std::endl;
    if ( retVal == NULL ) {
      (*debug) << " NULL\n" << std::endl;
    }
    else {
      for ( UInt i = 1; i <= numEqn_[identifier]; i++ ) {
        (*debug) << " " << i << " -> " << retVal[i] << std::endl;
      }
    }
#endif

    return retVal;
  }


  // ============
  //   GetGraph
  // ============
  BaseGraph* GraphManagerStdMat::GetGraph( const PdeIdType identifierPDE1,
					   const PdeIdType identifierPDE2){


    // Check if a graph object was already created
    if ( graph_ == NULL ) {
      (*error) << "GraphManagerStdMat: There exists no graph yet, "
               << "which could be returned by the GetGraph() method.";
      Error( __FILE__, __LINE__);
    }

    // Return pointer to the one and only graph object
    return graph_;

  }


  // ================
  //   GetIDBCGraph
  // ================
  BaseGraph* GraphManagerStdMat::GetIDBCGraph( const PdeIdType pdeID1,
                                               const PdeIdType pdeID2 ) const{
    return graphIDBC_;
  }


  // ==============
  //   PrintStats
  // ==============
  void GraphManagerStdMat::PrintStats( std::ostream *log ) {


    std::map<std::string,UInt>::iterator it;


    // ***********************************
    //  Assemble info for pretty-printing
    // ***********************************

    // Compute maximal column widths
    UInt cw1 = 0, cw2 = 0, cw3 = 0, cw4 = 0, tw = 0, aux;
    for ( UInt i = 1; i <= numPDEs_; i++) {
      aux = numEqn_[i] > 0 ? (UInt)std::log10( (float)numEqn_[i] ) + 1 : 1;
      cw1 = cw1 < aux ? aux : cw1;

      aux = numLastFreeDof_[i] > 0 ?
        (UInt)std::log10( (float)numLastFreeDof_[i] ) + 1 : 1;
      cw2 = cw2 < aux ? aux : cw2;

      aux = offset_[i] > 0 ? (UInt)std::log10( (float)offset_[i] ) + 1 : 1;
      cw3 = cw3 < aux ? aux : cw3;

      aux = offsetIDBC_[i] > 0 ?
        (UInt)std::log10( (float)offsetIDBC_[i] ) + 1 : 1;
      cw4 = cw4 < aux ? aux : cw4;
    }

    // Correct field widths for column headers
    cw1 = cw1 >  6 ? cw1 :  6;
    cw2 = cw2 > 11 ? cw2 : 11;
    cw3 = cw3 >  6 ? cw3 :  6;
    cw4 = cw4 > 10 ? cw4 : 10;

    // Set total width
    tw = 7 + 12 + cw1 + cw2 + cw3 + cw4 + 1;


    // *******************
    //  Report statistics
    // *******************

    // Begin report block
    (*log) << "\n " << std::setw(tw) << std::setfill( '=' ) << '=' << "\n"
           << std::setfill( ' ' )
           << " GRAPHMANAGER:\n\n";

    // Print standard graphmanager info
    (*log) << " Type of graph manager: GraphManagerStdMat\n"
           << " Number of sub-graphs: 1\n"
           << " Number of attached PDEs: " << numPDEs_
           << std::endl;

    // Output table showing identifier and offset
    // First the header
    (*log) << ' ' << std::setw(tw) << std::setfill( '-' ) << '-' << "\n"
           << "  PDE ID "
           << "| numEqn | lastFreeDof | offset | offsetIDBC\n"
           << ' ' << std::setw(tw) << std::setfill( '-' ) << '-' << "\n"
           << std::setfill( ' ' );

    // now the table rows
    log->setf( std::ios::right, std::ios::adjustfield );
    for ( UInt i = 1; i <= numPDEs_; i++) {
      (*log) << std::setw(8) << i
             << " | " << std::setw(cw1) << numEqn_[i]
             << " | " << std::setw(cw2) << numLastFreeDof_[i]
             << " | " << std::setw(cw3) << offset_[i]
             << " | " << std::setw(cw4) << offsetIDBC_[i]
             << std::endl;
    }
    (*log) << ' ' << std::setw(tw) << std::setfill( '-' ) << '-' << "\n"
           << std::setfill( ' ' );

    // number of degrees of freedom
    (*log) << " number of free dofs : " << numFreeDofs_
           << "\n number of fixed dofs: " << numFixedDofs_
           << "\n total number of dofs: "
           << numFreeDofs_ + numFixedDofs_ << std::endl;

    // graph entries
    (*log) << " number of non-zero matrix entries: " << graph_->GetNNE()
           << "\n number of non-zeros in auxilliary matrix: ";
    if ( graphIDBC_ == NULL ) {
      (*log) << 0 << '\n';
    }
    else {
      (*log) << graphIDBC_->GetNNE() << '\n';
    }

    // Close report block
    (*log) << ' ' << std::setw(tw) << std::setfill( '=' ) << '=' << "\n"
           << std::setfill( ' ' )  << std::endl;

  }
}
