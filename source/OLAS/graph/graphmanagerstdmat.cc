// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iterator>
#include <iomanip>

#include "OLAS/graph/graphmanagerstdmat.hh"

namespace CoupledField {


  // ===============
  //   Constructor
  // ===============
  GraphManagerStdMat::GraphManagerStdMat()
    : 
      reorderType_(BaseOrdering::NOREORDERING),
      offset_(NULL),
      offsetIDBC_(NULL),
      offset1_(0),
      offset2_(0),
      offsetIDBC1_(0),
      offsetIDBC2_(0),
      idPDE1_(NO_FCT_ID),
      idPDE2_(NO_FCT_ID),
      graph_(NULL),
      graphIDBC_(NULL),
      numPDEs_(0),
      numRegisteredPDEs_(0),
      numFixedDofs_(0),
      numFreeDofs_(0),
      numLastFreeDof_(NULL),
      numEqn_(NULL),
      reorderingDone_(false)
  {
  }


  // ==============
  //   Destructor
  // ==============
  GraphManagerStdMat::~GraphManagerStdMat() {


    // If no re-ordering was performed the pointers to the re-ordering
    // vectors are still NULL. If they were claimed by the PDEs they
    // were re-set to NULL and if they were not claimed, it is our
    // responsibility to delete then now.
    for ( UInt i = 0; i < numPDEs_; i++ ) {
      if ( newOrdering_[i].GetSize() ) {
        newOrdering_[i].Clear();
      }
    }
    newOrdering_.Clear();

    // Delete the graph objects
    delete graph_;
    delete graphIDBC_;

    // Delete auxilliary information arrays
    delete [] ( numLastFreeDof_ );
    delete [] ( offsetIDBC_ );
    delete [] ( offset_ );
    delete [] ( numEqn_ );

  }


  // =============
  //   SetupInit
  // =============
  void GraphManagerStdMat::SetupInit( UInt numPDEs ) {


    // Now we now for how many PDEs we are repsonsible
    numPDEs_ = numPDEs;

    // Allocate memory to store the offsets and further PDE size information
    NEWARRAY( numEqn_        , UInt, numPDEs_     );
    NEWARRAY( numLastFreeDof_, UInt, numPDEs_     );
    NEWARRAY( offset_        , UInt, numPDEs_ + 1 );
    NEWARRAY( offsetIDBC_    , Integer, numPDEs_ + 1 );

    // Allocate memory for the re-ordering vectors and set all
    // pointers to NULL for the case that no re-ordering is performed
    newOrdering_.Resize( numPDEs_ );
    for ( UInt i = 0; i < numPDEs; i++ ) {
      newOrdering_[i].Clear();
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
    StdVector<UInt> permutation;
    if ( reorderType_ != BaseOrdering::NOREORDERING ) {
      permutation.Resize( numFreeDofs_ );
    }

    // Only need permutation vector for first PDE, if we are
    // using the elimination approach and there are fixed dofs
    // for the first PDE
    UInt firstPermVec = 0;
    if ( numLastFreeDof_[0] == numEqn_[0] && reorderType_ == BaseOrdering::NOREORDERING ) {
      firstPermVec = 1;
    }

    // Each PDE, besides maybe the first, needs an individual reordering
    // array, so allocate memory for them taking into accout fixed dofs, too.
    for ( UInt i = firstPermVec; i < numPDEs_; i++ ) {
      newOrdering_[i].Resize( numEqn_[i] );
    }


    // *****************
    //  Re-order graphs
    // *****************

    // Trigger final steps of assembly
    graph_->FinaliseAssembly( &permutation );

    // Finish generation of auxilliary graph
    if ( graphIDBC_ != NULL ) {
      graphIDBC_->FinaliseAssembly( &permutation );
    }


    // ******************************
    //  Generate permutation vectors
    // ******************************

    // If re-ordering was performed, split large permutation vector
    // up into small vectors for each individual PDE
    UInt i, k, index;
    if ( reorderType_ != BaseOrdering::NOREORDERING ) {

      // Loop over all registered PDEs
      for ( i = 0; i < numPDEs_; i++ ) {

        // Copy new dof number from large vector into PDE individual vector
        index = 0;
        for ( k = offset_[i]; k < offset_[i+1]; k++ ) {
          newOrdering_[i][index] = permutation[k];
          index++;
        }
      }

      // Large vector is no longer needed
      permutation.Clear();
    }

    // No "real" re-ordering was performed. However, CFS++ needs the new
    // equation numbers, that include the offsets for all but the first PDE.
    // Therefore we generate a simple reordering for the other PDEs
    else {

      // Loop over all registered PDEs
      for ( i = firstPermVec; i < numPDEs_; i++ ) {
        index = 0;
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
      UInt fixedBlockInit = numFreeDofs_+1;
      if ( firstPermVec == 1 ) {
        fixedBlockInit += numEqn_[0] - numLastFreeDof_[0];
      }

      for ( i = firstPermVec; i < numPDEs_; i++ ) {
        for ( UInt j = 0; j < numEqn_[i] - numLastFreeDof_[i]; j++ ) {
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
  void GraphManagerStdMat::RegisterPDE( const FeFctIdType identifierPDE,
					const UInt n,
                                        const UInt numLastFreeDof,
					const BaseOrdering::ReorderingType reorder ) {


    // Step counter for the number of registered PDEs and check number
    numRegisteredPDEs_++;
    if ( numRegisteredPDEs_ > numPDEs_ ) {
      EXCEPTION("GraphManagerStdMat::RegisterPDE: You tried to "
               << "register a " << numRegisteredPDEs_ << "-th PDE "
               << "with id = '" << identifierPDE << "', but SetupInit "
               << "specified only " << numPDEs_ << " to be expected!");
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
      offset_[0] = 0;
      for ( UInt i = 0; i < numPDEs_; i++ ) {
        offset_[i+1] = offset_[i] + numLastFreeDof_[i];
      }

      UInt aux = 0;
      for ( UInt i = 0; i < numPDEs_; i++ ) {
        offsetIDBC_[i] = aux - numLastFreeDof_[i];
        aux += numEqn_[i] - numLastFreeDof_[i];
      }

      // Once the last PDE has registered itself, we know the number of
      // vertices in the graph and can thus construct the graph object.
      graph_ = new BaseGraph( numFreeDofs_, numFreeDofs_, reorderType_ );
      ASSERTMEM( graph_, sizeof(BaseGraph) );

      // Generate auxilliary graph object
      if ( numFixedDofs_ > 0 ) {
        graphIDBC_ = new IDBC_Graph( numFreeDofs_, numFixedDofs_ );
        ASSERTMEM( graphIDBC_, sizeof(IDBC_Graph) );
      }
    }
  }


  // ================
  //   AssembleInit
  // ================
  void GraphManagerStdMat::AssembleInit( const FeFctIdType identifierPDE1,
					 const FeFctIdType identifierPDE2,
                                         bool assemblingTranspose ) {


    // Perform a consistency check
    if ( identifierPDE1 == NO_FCT_ID ) {
      EXCEPTION("GraphManagerStdMat: First PDE identifier passed to "
               << "AssembleInit is empty!");
    }
    if ( graph_ == NULL ) {
      EXCEPTION("GraphManagerStdMat::AssembleInit: "
          << "Pointer to graph object = NULL! "
          << "Did you call RegisterPDE() for all " << numPDEs_
          << " PDEs?");
    }

    // Assign identifiers for checking in SetElementPos
    idPDE1_ = identifierPDE1;
    idPDE2_ = identifierPDE2;

    // Determine offsets
    if ( idPDE1_ > static_cast<Integer>(numPDEs_) ) {
      EXCEPTION("GraphManagerStdMat::AssembleInit: "
               << "PDE with identifier '" << idPDE1_ << "' was not "
               << "registered using RegisterPDE()!");
    }
    offset1_ = offset_[idPDE1_];
    offsetIDBC1_ = offsetIDBC_[idPDE1_];
    
    if ( idPDE2_ != NO_FCT_ID ) {
      if ( idPDE2_ > static_cast<Integer>(numPDEs_) ) {
        EXCEPTION("GraphManagerStdMat::AssembleInit: "
                 << "PDE with identifier '" << idPDE2_ << "' was not "
                 << "registered using RegisterPDE()!");
      }
      offset2_ = offset_[idPDE2_];
      offsetIDBC2_ = offsetIDBC_[idPDE2_];
    }

  }


  // =================
  //   SetElementPos
  // =================
  void GraphManagerStdMat::SetElementPos( const FeFctIdType identifierPDE1,
                                          const StdVector<Integer>& eqnNrs1,
                                          const FeFctIdType identifierPDE2,
                                          const StdVector<Integer>& eqnNrs2,
                                          bool setCounterPart ) {  


#ifdef DEBUG_GRAPHMANAGERSTDMAT

    if ( identifierPDE1 == NO_FCT_ID ) {
      EXCEPTION("GraphManagerStdMat: First PDE identifier passed to "
                << "SetElementPos is empty!");
    }
    if ( graph_ == NULL ) {
      EXCEPTION("GraphManagerStdMat::SetElementPos: "
               << "Pointer to graph object = NULL! "
	       << "Did you call RegisterPDE() for all " << numPDEs_
                << " PDEs?");
    }
    if ( idPDE1_ != identifierPDE1 || idPDE2_ != identifierPDE2 ) {
      EXCEPTION("Consistency error detected! You specified to AssembleInit "
               << " idPDE1 = '" << idPDE1_ << "',  idPDE2 = '" <<  idPDE2_
               << "', but SetElementPos received idPDE1 = '"
               << identifierPDE1 << "',  idPDE2 = '" <<  identifierPDE2
                << "'. Maybe you forgot to call AssembleDone & AssembleInit?");
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
    UInt length1 = eqnNrs1.GetSize();
    UInt length2 = eqnNrs2.GetSize();

    // STEP 1: Generate vertex list from first connect array, dropping
    //         equation numbers for dofs fixed by inhomogeneous Dirichlet
    //         boundary conditions and changing the sign of those fixed by
    //         constraints
    UInt aux;
    for ( UInt i = 0; i < length1; i++ ) {
      aux = abs(eqnNrs1[i]);
      if ( aux > 0 ) {
        if ( aux <= numLastFreeDof_[identifierPDE1] ) {
          vertexList1_.push_back( aux + offset1_ - 1);
        }
        else {
          vertexList2_.push_back( aux + offsetIDBC1_ - 1);
        }
      }
    }

    // STEP 2: Split the second connect array into two edge lists, one for
    //         the graph and one for the IDBCgraph (which handles the dofs
    //         fixed by inhomogeneous Dirichlet boundary conditions)
    for ( UInt i = 0; i < length2; i++ ) {
      aux = abs(eqnNrs2[i]);
      if ( aux > 0 ) {
        if ( aux > numLastFreeDof_[identifierPDE2] ) {
          edgeList2_.push_back( aux + offsetIDBC2_ - 1);
        }
        else {
          edgeList1_.push_back( aux + offset2_ - 1);
        }
      }
    }
 
#ifdef DEBUG_GRAPHMANAGERSTDMAT
    
    // output original connectivity
    (std::cerr) << "\n GraphManagerStdMat::AdaptConnects\n"
             << " idPDE1 = " << identifierPDE1 << '\n'
             << " idPDE2 = " << identifierPDE2 << '\n'
             << " offset1 = " << offset1_ << '\n'
             << " offset2 = " << offset2_ << '\n'
             << " offsetIDBC1 = " << offsetIDBC1_ << '\n'
             << " offsetIDBC2 = " << offsetIDBC2_ << '\n'
             << " numLastFreeDof1 = " << numLastFreeDof_[identifierPDE1] << '\n'
             << " numLastFreeDof2 = " << numLastFreeDof_[identifierPDE2] << '\n';
    (std::cerr) << " connect1 ";
    for ( UInt i = 0; i < length1; i++ ) {
      (std::cerr) << eqnNrs1[i] << " ";
    }
    (std::cerr) << std::endl;
    (std::cerr) << " connect2 ";
    for ( UInt i = 0; i < length2; i++ ) {
      (std::cerr) << eqnNrs2[i] << " ";
    }
    (std::cerr) << std::endl;

    // output new connectivity
    (std::cerr) << " vertexList1: ";
    for ( UInt i = 0; i < vertexList1_.size(); i++ ) {
      (std::cerr) << vertexList1_[i] << " ";
    }
    (std::cerr) << std::endl;
    (std::cerr) << " vertexList2: ";
    for ( UInt i = 0; i < vertexList2_.size(); i++ ) {
      (std::cerr) << vertexList2_[i] << " ";
    }
    (std::cerr) << std::endl;
    (std::cerr) << " edgeList1: ";
    for ( UInt i = 0; i < edgeList1_.size(); i++ ) {
      (std::cerr) << edgeList1_[i] << " ";
    }
    (std::cerr) << std::endl;
    (std::cerr) << " edgeList2: ";
    for ( UInt i = 0; i < edgeList2_.size(); i++ ) {
      (std::cerr) << edgeList2_[i] << " ";
    }
    (std::cerr) << std::endl;
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
  void GraphManagerStdMat::GetReordering( const FeFctIdType identifier, 
                                          StdVector<UInt>& order ) {

    // Get registration number of PDE
    if ( identifier > static_cast<Integer>(numPDEs_) ) {
      EXCEPTION("GraphManagerStdMat::GetReordering: "
               << "PDE with identifier '" << identifier << "' was not "
               << "registered using RegisterPDE()!");
    }

    // Test, whether we can return a re-ordering vector
    if ( reorderingDone_ == false ) {
      EXCEPTION("GraphManagerStdMat::GetReordering: "
               << "No reordering vector available since the graph has not "
               << "been reordered, yet!");
    }
    else if ( !newOrdering_[identifier].GetSize() &&
              reorderType_ != BaseOrdering::NOREORDERING ) {
      EXCEPTION("GraphManagerStdMat::GetReordering: "
               << "No reordering vector available for PDE with identifier '"
               << identifier << "' Maybe it was already claimed?");
    }

    // By passing the pointer to the array containing the re-ordering
    // information to the caller, this class forgets about the re-ordering
    order = newOrdering_[identifier];
    newOrdering_[identifier].Clear();

#ifdef DEBUG_GRAPHMANAGERSTDMAT
    (*debug) << "\nGraphManagerStdMat::GetReordering:\n"
             << "Re-ordering for PDE with ID = '" << identifier
             << std::endl;
    if ( !order.GetSize() ) {
      (*debug) << " NULL\n" << std::endl;
    }
    else {
      for ( UInt i = 0; i < numEqn_[identifier]; i++ ) {
        (*debug) << " " << i << " -> " << order[i] << std::endl;
      }
    }
#endif
  }


  // ============
  //   GetGraph
  // ============
  BaseGraph* GraphManagerStdMat::GetGraph( const FeFctIdType identifierPDE1,
					   const FeFctIdType identifierPDE2){


    // Check if a graph object was already created
    if ( graph_ == NULL ) {
      EXCEPTION("GraphManagerStdMat: There exists no graph yet, "
               << "which could be returned by the GetGraph() method.");
    }

    // Return pointer to the one and only graph object
    return graph_;

  }


  // ================
  //   GetIDBCGraph
  // ================
  BaseGraph* GraphManagerStdMat::GetIDBCGraph( const FeFctIdType pdeID1,
                                               const FeFctIdType pdeID2 ) const{
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
    for ( UInt i = 0; i < numPDEs_; i++) {
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
    for ( UInt i = 0; i < numPDEs_; i++) {
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
