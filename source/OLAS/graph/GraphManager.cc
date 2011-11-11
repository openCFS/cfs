#include <iterator>
#include <iomanip>

#include "OLAS/graph/GraphManager.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "OLAS/graph/BaseGraph.hh"
#include "OLAS/graph/IDBC_Graph.hh"


DECLARE_LOG(graphMan)
DEFINE_LOG(graphMan, "graphManager")

namespace CoupledField {


  // ===============
  //   Constructor
  // ===============
  GraphManager::GraphManager() {

    numBlocks_           = 0;
    numRegisteredBlocks_ = 0;

    reorderingDone_      = false;
    registrationDone_    = false;
    WARN("GraphManager: Adapt documentation!")
  }


  // ==============
  //   Destructor
  // ==============
  GraphManager::~GraphManager() {

    // Delete the graph objects
    for ( UInt i = 0; i < numBlocks_ * numBlocks_; i++ ) {
      delete graph_[i];
    }
    graph_.Clear();

    // Delete the IDBC graph objects
    for ( UInt i = 0; i < numBlocks_ * numBlocks_; i++ ) {
      delete graphIDBC_[i];
    }
    graphIDBC_.Clear();

  }


  // =============
  //   SetupInit
  // =============
  void GraphManager::SetupInit( UInt numBlocks,
                                bool useDistinctGraphs ) {
    
    LOG_TRACE(graphMan) << "Initializing GraphManager with " 
                         << numBlocks << " blocks";
    LOG_DBG(graphMan) << "\tuseDistinctGraphs:" 
        << (useDistinctGraphs ? "yes" : "no" );
    
    // Note: Currently we do not support distinct graphs
    // for the system. This is a feature yet to come
    if( useDistinctGraphs ) {
      EXCEPTION("Definition of distinct graphs not supported yet");
    }

    // Now we now for how many blocks we are responsible
    numBlocks_ = numBlocks;

    // ... and can build the empty graph pointer matrix
    graph_.Resize( numBlocks_ * numBlocks_ );
    graph_.Init( NULL );


    // ... and the empty IDBC graph array
    graphIDBC_.Resize( numBlocks_ * numBlocks_ );
    graphIDBC_.Init( NULL );

    // Resize also the vertex and edge lists
    vertexList1_.resize( numBlocks_ );
    vertexList2_.resize( numBlocks_ );
    edgeList1_.resize( numBlocks_ );
    edgeList2_.resize( numBlocks_ );
    
    // Allocate memory to store the permutation vectors for the reordering
    // of the unknowns of each block
    newOrdering_.Resize( numBlocks_ );
    for ( UInt i = 0; i < numBlocks_; i++ ) {
      newOrdering_[i].Resize(0);
    }

    // Setup empty array for coupling flags
    isCoupled_.Resize( numBlocks_ * numBlocks_ );
    isCoupled_.Init( false );
    
    // Resize array with block information
    blockInfo_.Resize(numBlocks);
  }


  // =============
  //   SetupDone
  // =============
  void GraphManager::SetupDone() {
    
    LOG_TRACE(graphMan) << "Finalize Graphmanager (SetupDone)";
    
    // ----------------------------------
    //   D I A G O N A L    B L O C K S 
    // ----------------------------------
    for( UInt iBlock = 0; iBlock < numBlocks_; ++iBlock ) {
      UInt idx = ComputeIndex( iBlock, iBlock );

      // Before we finalize the graph (i.e. make column entries unique,
      // sort entries, perform reordering, convert to CRS-structure),
      // we have to provide the sub-matrix block information, so the
      // graph can reorder internally the entries in a way that entries
      // in one submatrix block are numbered consecutively.
      if( blockInfo_[iBlock]->hasSubBlocks) {
        graph_[idx]->SetBlockInfo(&(blockInfo_[iBlock]->indexBlocks));
      }
      
      // Finalize assembly of graph
      LOG_DBG(graphMan) << "Finalize diagonal graph (" << iBlock
                        << ", " << iBlock << ")";
      graph_[idx]->FinaliseAssembly( false, &newOrdering_[iBlock] );
      LOG_DBG3(graphMan) << "Reordering array is "  
                         << newOrdering_[iBlock].ToString(); 
      
      // Now finalize the assembly of the associated IDBC graph
      if ( graphIDBC_[idx] != NULL ) {
        LOG_DBG(graphMan) << "Finalize IDBC graph (" << iBlock
                          << ", " << iBlock << ")";
        graphIDBC_[idx]->FinaliseAssembly( &newOrdering_[iBlock] );
      }
    }

    // -------------------------------------------
    //   O F F  - D I A G O N A L    B L O C K S 
    // -------------------------------------------
    for( UInt iRow = 0; iRow < numBlocks_; ++iRow ) {
      for( UInt iCol = 0; iCol < numBlocks_; ++iCol ) {
        
        // leave at diagonal blocks
        if( iRow == iCol)
          continue;

        UInt idx = ComputeIndex( iRow, iCol );

        //  Finalize assembly of graph (sorting, re-ordering, conversion to 
        // CRS-structure)
        // Note: For the off-diagonal entries we use the re-ordering arrays
        //       of the related row/col diagonal blocks
        LOG_DBG(graphMan) << "Finalize off-diagonal graph (" << iRow
                          << ", " << iCol << ")";
        graph_[idx]->FinaliseAssembly( true, &newOrdering_[iRow],
                                       &newOrdering_[iCol] );

        // Set also the corresponding diagonal graph objects of the row/col.
        // This is needed, as the matrix sub block definition of the coupling
        // graph has to be taken from the corresponding diagonal blocks.
        BaseGraph * rowDiagGraph = graph_[ComputeIndex(iRow,iRow)];
        BaseGraph * colDiagGraph = graph_[ComputeIndex(iCol,iCol)];
        graph_[idx]->SetRowColDiagGraphs( rowDiagGraph, colDiagGraph );


        // Now finalize the assembly of the associated IDBC graph
        // Note: Also here we use the external re-ordering provided by
        //       the diagonal block of the same column
        if ( graphIDBC_[idx] != NULL ) {
          LOG_DBG(graphMan) << "Finalize off-diagonal IDBC-graph (" 
                            << iRow << ", " << iCol << ")";
          graphIDBC_[idx]->FinaliseAssembly( &newOrdering_[iRow] );
        }
      }
    }
    
    // set flag for successful reordering
    reorderingDone_ = true;

    // Print statistics to standard log stream
    PrintStats();
  }


  // ===============
  //   RegisterBlock
  // ===============
  void GraphManager::RegisterBlock( const UInt blockNum,
                                    SBMBlockInfo* blockInfo,
                                    const BaseOrdering::ReorderingType reorder ) {
    
    LOG_TRACE(graphMan) << "Registering block " << blockNum;
    
    if (IS_LOG_ENABLED(graphMan, dbg2) ){
    LOG_DBG(graphMan) << "Detailed block information:";
    LOG_DBG(graphMan) << "\ttotal size: " << blockInfo->size;
    LOG_DBG(graphMan) << "\tlastFreeIndex: " << blockInfo->numLastFreeIndex;
    LOG_DBG(graphMan) << "\thasSubBlocks: " << blockInfo->hasSubBlocks;
    }

    // Be cautious
    if ( registrationDone_ == true ) {
      EXCEPTION("Attempt to use RegisterBlock() after end of "
               << "registration phase, i.e. after the first call to "
               << "AssembleInit()!");
    }
    
    
    // Store SBM block information
    blockInfo_[blockNum] = blockInfo; 

    // Step counter for the number of registered blocks and check number
    numRegisteredBlocks_++;
    if ( numRegisteredBlocks_ > numBlocks_ ) {
      EXCEPTION("GraphManager::RegisterBlock: You tried to "
               << "register a " << numRegisteredBlocks_ << "-th block "
               << "but SetupInit specified only " << numBlocks_ 
               << " to be expected!");
    }
    
    if( numRegisteredBlocks_ == numBlocks_ )
      registrationDone_ = true;

    // Generate graph object for this block
    UInt idx = ComputeIndex( blockNum, blockNum );
    graph_[idx] = new BaseGraph( blockInfo->numLastFreeIndex, blockInfo->numLastFreeIndex,
                                 reorder );
    if ( graph_[idx] == NULL ) {
      EXCEPTION("Generation of graph object for block #" << blockNum 
                << " failed!");
    }

    LOG_DBG(graphMan) << " GraphManager: Generated sub-graph ("
             << blockNum << ", " << blockNum << ")"
             << " for a " << blockInfo->numLastFreeIndex << " x " 
             << blockInfo->numLastFreeIndex
             << " matrix" << std::endl;
    
    // Generate IDBC graph object for this block
    GenerateIDBCGraph( blockNum, blockNum );
    

    // If reordering is going to be performed for the current block then
    // we need to allocate memory to store the resulting permutation
    // vector
    if ( reorder != BaseOrdering::NOREORDERING ) {
      newOrdering_[blockNum].Resize( blockInfo->numLastFreeIndex );
    }
  }

  // =================
  //   SetElementPos
  // =================
  void GraphManager::SetElementPos( const StdVector<UInt>& rowBlocks,
                                    const StdVector<UInt>& rowNums,
                                    const StdVector<UInt>& colBlocks,
                                    const StdVector<UInt>& colNums,
                                    FEMatrixType matrixType,
                                    bool setCounterPart ) {

    // Just some logging for debugging
    LOG_TRACE(graphMan) << "Setting element connectivity";
    if( IS_LOG_ENABLED(graphMan, dbg3) ) {
      LOG_DBG3(graphMan) << "setCounterPart: " << (setCounterPart ? "yes" : "no");
      LOG_DBG3(graphMan) << "\t(rowBlock, rowNum)";
      for( UInt i = 0; i < rowBlocks.GetSize(); ++i ) {
        LOG_DBG3(graphMan) << "\t(" << rowBlocks[i] << ", " << rowNums[i] << ")"; 
      }
      LOG_DBG3(graphMan) << "\t(colBlock, colNum)";
      for( UInt i = 0; i < colBlocks.GetSize(); ++i ) {
        LOG_DBG3(graphMan) << "\t(" << colBlocks[i] << ", " << colNums[i] << ")"; 
      }
    }
    
    // Clear the arrays
    for( UInt i = 0; i < numBlocks_; ++i ) {
      vertexList1_[i].clear();
      vertexList2_[i].clear();
      edgeList1_[i].clear();
      edgeList2_[i].clear();
    }

    UInt numRows = rowBlocks.GetSize();
    UInt numCols = colBlocks.GetSize();

    // Loop over all rows
    for( UInt iRow = 0; iRow < numRows; ++iRow ) {
      // get hold of block numbers and indices
      const UInt & rowBlock = rowBlocks[iRow];
      const UInt & rowNum = rowNums[iRow];

      // Compute index of graph in graph pointer matrix
      // get hold of vertex and edgelists
      std::vector<UInt> & vList1 = vertexList1_[rowBlock];
      std::vector<UInt> & vList2 = vertexList2_[rowBlock];
      // get limits of free indices
      const UInt & lastFreeRowIndex = blockInfo_[rowBlock]->numLastFreeIndex;

      // STEP 1: Generate vertex list from first connect array, dropping
      //         equation numbers for dofs fixed by inhomogeneous Dirichlet
      //         boundary conditions
      if ( rowNum > 0 ) {
        if ( rowNum <= lastFreeRowIndex ) {
          vList1.push_back( rowNum - 1);
        }
        else {
          vList2.push_back( rowNum - lastFreeRowIndex -1 );
        }
      }
    }

    // Loop over all columns
    for( UInt iCol = 0; iCol < numCols; ++iCol ) {

      // get hold of block numbers and indices
      const UInt & colBlock = colBlocks[iCol];
      const UInt & colNum = colNums[iCol];

      // get hold of vertex and edgelists
      std::vector<UInt> & eList1 = edgeList1_[colBlock];
      std::vector<UInt> & eList2 = edgeList2_[colBlock];

      // get limits of free indices
      const UInt & lastFreeColIndex = blockInfo_[colBlock]->numLastFreeIndex;

      // STEP 2: Split the second connect array into two edge lists, one for
      //         the graph and one for the IDBCgraph (which handles the indices
      //         fixed by inhomogeneous Dirichlet boundary conditions)
      if( colNum > 0 ) {
        if ( colNum > lastFreeColIndex ) {
          eList2.push_back( colNum - lastFreeColIndex - 1);
        }
        else {
          eList1.push_back( colNum - 1);
        }
      }
    } // loop over cols

    
    // Now we have the edge/vertex list for every block available, so
    // we can loop explictly over vertices (=rows) and edges (=cols)
    // of each block.
    
    // loop over all blocks and pass for every block the information to
    // the corresponding graph / IDBC graph
    for( UInt row = 0; row < numBlocks_; ++row ) {
      for( UInt col = 0; col < numBlocks_; ++col ) {

        // Compute index of graph in graph pointer matrix
        UInt idx = ComputeIndex( row, col );
        
        // Generate coupling graph and also transpose if necessary
        if ( row != col ) {
          GenerateCouplingGraph( row, col );
          if ( setCounterPart ) {
            GenerateCouplingGraph( col, row );
          }

          // Generate IDBC graph and its transpose if necessary
          GenerateIDBCGraph( row, col );
          if ( setCounterPart ) {
            GenerateIDBCGraph( col, row );
          }
        }

        // --- logging output ---
        if( IS_LOG_ENABLED(graphMan, dbg3) ) {
          LOG_DBG3(graphMan) << "IDBC/Graph insertion for block (" 
              << row << ", " << col << ")";
          
          LOG_DBG3(graphMan) << "vertexList1: ";
          for(UInt i=0; i <vertexList1_[row].size(); ++i ) 
            LOG_DBG3(graphMan) << "\t" << vertexList1_[row][i];
          
          LOG_DBG3(graphMan) << "vertexList2: ";
          for(UInt i=0; i <vertexList2_[row].size(); ++i ) 
            LOG_DBG3(graphMan) << "\t" << vertexList2_[row][i];
          
          LOG_DBG3(graphMan) << "edgeList1: ";
          for(UInt i=0; i <edgeList1_[col].size(); ++i ) 
            LOG_DBG3(graphMan) << "\t" << edgeList1_[col][i];
          
          LOG_DBG3(graphMan) << "edgeList2: ";
          for(UInt i=0; i <edgeList2_[col].size(); ++i ) 
            LOG_DBG3(graphMan) << "\t" << edgeList2_[col][i];
        }
        
        // --- logging output ---
        // Insert information into graph for real dofs
        graph_[idx]->AddVertexNeighbours( vertexList1_[row], edgeList1_[col] );

        // Insert information into graph for fixed dofs
        graphIDBC_[idx]->AddVertexNeighbours( vertexList1_[row], edgeList2_[col] );
//        if( row!= col ) {
//          std::cerr << "IDBC (" << row << ", " << col << ": Inserting ";
//          for( UInt i = 0; i < vertexList1_[row].size(); ++i ) 
//            std::cerr << vertexList1_[row][i] << ", ";
//          std::cerr << " AND ";
//          for( UInt i = 0; i < edgeList2_[col].size(); ++i ) 
//            std::cerr << edgeList2_[col][i] << ", ";
//          std::cerr << std::endl;
//        }

        // Check for assembly of counter part
        if ( (row != col) && setCounterPart == true ) {

          idx = ComputeIndex( col, row );
          LOG_DBG3(graphMan) << "IDBC: Inserting into (" << col << ", " << row << ")" << std::endl;

          // Insert information into (transpose) graph for real dofs
          graph_[idx]->AddVertexNeighbours( edgeList1_[col], vertexList1_[row]);
          //graph_[idx]->AddVertexNeighbours( vertexList1_[col], edgeList1_[row]);

          // Insert information into graph for fixed dofs

          graphIDBC_[idx]->AddVertexNeighbours( vertexList1_[col], edgeList2_[row] );
//          std::cerr << "IDBC (" << col << ", " << row << "): Inserting ";
//          for( UInt i = 0; i < vertexList1_[col].size(); ++i ) 
//            std::cerr << vertexList1_[col][i] << ", ";
//          std::cerr << " AND ";
//          for( UInt i = 0; i < edgeList2_[row].size(); ++i ) 
//            std::cerr << edgeList2_[row][i] << ", ";
//          std::cerr << std::endl;
        }
      } // row loop
    } // col loop
  }


  // =================
  //   GetReordering
  // =================
  void GraphManager::GetReordering( UInt blockNum,
                                    StdVector<UInt>& order) {
    
    LOG_TRACE(graphMan) << "Returning reordering for block #" << blockNum;
    
    // Small consistency check
    if ( blockNum > numBlocks_ ) {
      EXCEPTION("GraphManager::GetReordering: "
               << "block with number '" << blockNum << "' was not "
               << "registered using RegisterBlock()!");
    }

    // Test, whether we can return a re-ordering vector
    if ( reorderingDone_ == false ) {
      EXCEPTION("GraphManager::GetReordering: "
               << "No reordering vector available since the graphs have not "
               << "been reordered, yet!");
    }

    // By passing the pointer to the array containing the re-ordering
    // information to the caller, this class forgets about the re-ordering
    order = newOrdering_[blockNum];
    newOrdering_[blockNum].Clear();
  }


  // ============
  //   GetGraph
  // ============
  BaseGraph* GraphManager::GetGraph( UInt rowNum,
                                     UInt colNum ) {

    LOG_TRACE(graphMan) << "Returning graph (" << rowNum << ", " << colNum << ")";

    // Determine which graph we are looking for
    UInt idx = ComputeIndex( rowNum, colNum );

    // Check if a graph object was already created
    if ( graph_[idx] == NULL) {
      EXCEPTION("GraphManager: There exists no graph "
               << "for the block index  pair ( " << rowNum
               << " , " << colNum
               << ") which could be returned by the GetGraph() method.");
    }

    // Return pointer to the graph object
    return graph_[idx];
  }


  // ================
  //   GetIDBCGraph
  // ================
  BaseGraph* GraphManager::GetIDBCGraph( UInt rowNum,
                                         UInt colNum ) const{
    LOG_TRACE(graphMan) << "Returning IDBC-graph (" 
                        << rowNum << ", " << colNum << ")";

    UInt idx = ComputeIndex( rowNum, colNum );
    if ( graphIDBC_[idx] == NULL ) {
      EXCEPTION("GraphManager::GetIDBCGraph: "
               << "An IDBC graph object for block index pair (" << rowNum
               << " , " << colNum << ") does not or not yet exist!");
    }
    // Return pointer to the IDBC graph object
    return graphIDBC_[idx];
  }

  // ==============
  //   PrintStats
  // ==============
  void GraphManager::PrintStats() {

//    std::ostream * out = &std::cout;

//    // ***********************************
//    //  Assemble info for pretty-printing
//    // ***********************************
//
//    // Compute maximal column widths (block table)
//    UInt cw1 = 0, cw2 = 0, cw3 = 0, cw4 = 0, tw = 0, aux;
//    for ( UInt i = 1; i <= numBlocks_; i++) {
//      aux = numEqn_[i] > 0 ? (UInt)std::log10( (float)numEqn_[i] ) + 1 : 1;
//      cw1 = cw1 < aux ? aux : cw1;
//
//      aux = numLastFreeDof_[i] > 0 ?
//        (UInt)std::log10( (float)numLastFreeDof_[i] ) + 1 : 1;
//      cw2 = cw2 < aux ? aux : cw2;
//    }
//
//    // Compute maximal column widths (sub-graph table)
//    UInt idx = 0;
//    for ( UInt i = 1; i <= numBlocks_; i++ ) {
//      for ( UInt j = 1; j <= numBlocks_; j++ ) {
//        idx = ComputeIndex(i,j);
//        if ( graph_[idx] != NULL ) {
//          aux = graph_[idx]->GetSize() > 0 ?
//            (UInt)std::log10( (float)graph_[idx]->GetSize() ) + 1 : 1;
//          cw3 = cw3 < aux ? aux : cw3;
//          aux = graph_[idx]->GetNNE() > 0 ?
//            (UInt)std::log10( (float)graph_[idx]->GetNNE() ) + 1 : 1;
//          cw4 = cw4 < aux ? aux : cw4;
//        }
//      }
//    }
//
//    // Correct field widths for column headers
//    cw1 = cw1 >  6 ? cw1 :  6;
//    cw2 = cw2 > 11 ? cw2 : 11;
//    cw3 = cw3 >  9 ? cw3 :  9;
//    cw4 = cw4 >  6 ? cw4 :  6;
//
//    // Set total width
//    tw = 7 + 12 + cw1 + cw2 + 5;
//
//
//    // *******************
//    //  Report statistics
//    // *******************
//
//    // Begin report block
//    (*log) << "\n " << std::setw(tw) << std::setfill( '=' ) << '=' << "\n"
//           << std::setfill( ' ' )
//           << " GRAPHMANAGER:\n\n";
//
//    // Determine number of sub-graphs we are holding
//    UInt numGraphs = 0;
//    for ( UInt i = 1; i <= numBlocks_ * numBlocks_; i++ ) {
//      if ( graph_[i] != NULL ) {
//        numGraphs++;
//      }
//    }
//
//    // Print standard graphmanager info
//    (*log) << " Type of graph manager: GraphManager\n"
//           << " Number of attached blcocks: " << numBlocks_
//           << std::endl;
//
//    // Output table showing identifier and number of unknowns
//    (*log) << ' ' << std::setw(tw) << std::setfill( '-' ) << '-' << "\n"
//           << "  Block  "
//           << "| numEqn | lastFreeDof\n"
//           << ' ' << std::setw(tw) << std::setfill( '-' ) << '-' << "\n"
//           << std::setfill( ' ' );
//
//    // now the table rows
//    log->setf( std::ios::right, std::ios::adjustfield );
//    for ( UInt i = 1; i <= numBlocks_; i++) {
//      (*log) << std::setw(8) << i
//             << " | " << std::setw(cw1) << numEqn_[i]
//             << " | " << std::setw(cw2) << numLastFreeDof_[i]
//             << std::endl;
//    }
//    (*log) << ' ' << std::setw(tw) << std::setfill( '-' ) << '-' << "\n"
//           << std::setfill( ' ' );
//
//    // Output table showing sub-graph information
//    (*log) << " Number of sub-graphs: " << numGraphs << std::endl;
//    (*log) << ' ' << std::setw(tw) << std::setfill( '-' ) << '-' << "\n"
//           << "  row | col | #vertices | #edges\n"
//           << ' ' << std::setw(tw) << std::setfill( '-' ) << '-' << "\n"
//           << std::setfill( ' ' );
//
//    idx = 0;
//    for ( UInt i = 1; i <= numBlocks_; i++ ) {
//      for ( UInt j = 1; j <= numBlocks_; j++ ) {
//        idx = ComputeIndex(i,j);
//        if ( graph_[idx] != NULL ) {
//          (*log) << std::setw(5) << i << " | "
//                 << std::setw(3) << j << " | "
//                 << std::setw(cw3) << graph_[idx]->GetSize() << " | "
//                 << std::setw(cw4) << graph_[idx]->GetNNE()
//                 << std::endl;
//        }
//      }
//    }
//
//
//    // Close report block
//    (*log) << ' ' << std::setw(tw) << std::setfill( '=' ) << '=' << "\n"
//           << std::setfill( ' ' )  << std::endl;
  }


  // =========================
  //   GenerateCouplingGraph
  // =========================
  void GraphManager::GenerateCouplingGraph( UInt rowNum,
                                            UInt colNum ) {

    UInt idx = ComputeIndex( rowNum, colNum );
    // Safety check: As this method can get called several times,
    // we silently leave, if the coupling graph already exists
    if ( graph_[idx] != NULL ) {
      return;
    }
    

    // Generate graph object
    graph_[idx] = new BaseGraph( blockInfo_[rowNum]->numLastFreeIndex,
                                 blockInfo_[colNum]->numLastFreeIndex, 
                                 BaseOrdering::NOREORDERING );

    if ( graph_[idx] == NULL ) {
      EXCEPTION("GraphManager: Generation of sub-graph "
               << "for pair (" << rowNum << " , " << colNum
               << ") and a " << blockInfo_[rowNum]->numLastFreeIndex
               << " x " << blockInfo_[colNum]->numLastFreeIndex << " matrix failed ");
    }
    
    // log message
    LOG_DBG(graphMan) << " GraphManager: Generated sub-graph for  pair ("
        << rowNum << " , " << colNum << ") and a "
        << blockInfo_[rowNum]->numLastFreeIndex << " x " 
        << blockInfo_[colNum]->numLastFreeIndex << " matrix " << std::endl;
  }


  // =====================
  //   GenerateIDBCGraph
  // =====================
  void GraphManager::GenerateIDBCGraph( UInt rowNum,
                                        UInt colNum ) {

    UInt idx = ComputeIndex( rowNum, colNum );

    // Safety check: As this method can get called several times,
    // we silently leave, if the IDBC graph already exists
    if( graphIDBC_[idx] != NULL ) 
      return;

    // Compute number of fixed column indices
    UInt fixedDofs = blockInfo_[colNum]->size - 
                     blockInfo_[colNum]->numLastFreeIndex;

    // If there are fixed indices in the column block, we have to
    // generate an IDBC graph object for this pair
    if ( fixedDofs > 0 ) {

      // Generate IDBC graph object
      graphIDBC_[idx] = new IDBC_Graph( blockInfo_[rowNum]->numLastFreeIndex, 
                                        fixedDofs );
      if ( graphIDBC_[idx] == NULL ) {
        EXCEPTION(" GraphManager: Generation of IDBC sub-graph " 
                 << "for index pair (" << rowNum << " , " << colNum
                 << ") and a " << blockInfo_[rowNum]->numLastFreeIndex
                 << " x " << fixedDofs << " matrix failed ");
      }
    }

    if ( fixedDofs > 0 ) {
      // log message
      LOG_DBG(graphMan) << " GraphManager: Generated IDBC sub-graph "
          << "for index pair (" << rowNum << " , " << colNum
          << ") and a " << blockInfo_[rowNum]->numLastFreeIndex
          << " x " << fixedDofs << " matrix";
    }
  }

} // namespace
