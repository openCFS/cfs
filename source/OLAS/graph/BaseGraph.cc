#include "General/defs.hh"
#include <iostream>
#include <cmath>
#include <vector>
#include <set>
#include <algorithm>

#include "Utils/Timer.hh"

// use boost accumulators to gather usage statistics of the graph
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/sum.hpp>
#include <boost/accumulators/statistics/moment.hpp>

#include "DataInOut/Logging/LogConfigurator.hh"

using namespace boost::accumulators;

#include <def_use_metis.hh>

#ifdef USE_METIS

// this namespace serves to hide metis names
extern "C"{
#include "metis.h"
}
#endif


#include "Domain/Domain.hh"
#include "OLAS/graph/BaseGraph.hh"
#include "OLAS/graph/Sloan.hh"


namespace CoupledField {
  DEFINE_LOG(graphLogger, "graph")

  // ***************
  //   Constructor
  // ***************
  BaseGraph::BaseGraph( UInt nRows, UInt nCols ) {

    // Avoid problems with partially empty graphs
    if ( nRows == 0 || nCols ==0 ) {
      EXCEPTION("Refusing to generate a graph for a "
               << nRows << " x " << nCols << " matrix!");
    }

    // Initialize attributes
    csEdges_        = NULL;
    csNodes_        = NULL;
    amAssembled_    = false;
    amReordered_    = false;
    newOrder_       = BaseOrdering::NOREORDERING;
    numNodes_       = nRows;
    numNonDiagEntries_=0;
    numRowsMat_     = nRows;
    numColsMat_     = nCols;
    bwlower_        = 0;
    bwupper_        = 0;
    unsortedBlocks_ = NULL;
    rowDiagGraph_   = NULL;
    colDiagGraph_   = NULL;
    setToElemDone_  = false;

    // Allocate memory for linked lists
    if ( numNodes_ > 0 ) {
      element_ = new std::vector<UInt>[numNodes_];
      setElements_ = new std::unordered_set<UInt>[numNodes_];
    }
    else {
      element_ = NULL;
    }
  }


  // *************************
  //   Alternate Constructor
  // *************************
  BaseGraph::BaseGraph( UInt nRows, UInt nCols, UInt nne, UInt *cs_nodes,
                        UInt *cs_edges ) :
    amReordered_(false),
    amAssembled_(false) {

    numNodes_    = nRows;
    numRowsMat_  = nRows;
    numColsMat_  = nCols;
    nne_         = nne;
    numNonDiagEntries_ = 0;
    newOrder_    = BaseOrdering::NOREORDERING;
    csNodes_     = cs_nodes;
    csEdges_     = cs_edges;
    element_     = NULL;
    bwlower_     = 0;
    bwupper_     = 0;
    rowDiagGraph_   = NULL;
    colDiagGraph_   = NULL;

    amAssembled_ = true;
    setToElemDone_  = false;
  }


  // *******************
  //   Deep Destructor
  // *******************
  BaseGraph::~BaseGraph() {


    delete [] ( csNodes_ );  csNodes_  = NULL;
    delete [] ( csEdges_ );  csEdges_  = NULL;

    // Note: element_ should already have been deleted
    // in FinaliseAssembly()
    delete [] ( element_ );  element_  = NULL;

  }


  // ***********************
  //   AddVertexNeighbours
  // ***********************
  void BaseGraph::AddVertexNeighbours( std::vector<UInt>& vertexList,
                                       std::vector<UInt>& neighbourList ) {

    // For each vertex in vertexList add all vertices in the neighbourList
    // as neighbours / edges
    UInt i, j;
    for ( i = 0; i < vertexList.size(); i++ ) {
      for ( j = 0; j < neighbourList.size(); j++ ) {
        setElements_[ vertexList[i] ].insert(neighbourList[j]);
//        // Insert edge only if not inserted yet
//        if( std::find(element_[ vertexList[i] ].begin(),
//                      element_[ vertexList[i] ].end(), neighbourList[j] )
//        == element_[ vertexList[i] ].end() ) {
//        element_[ vertexList[i] ].push_back( neighbourList[j] );
//      }
    }
  }
  }

  // ****************
  //   SetBlockInfo
  // ****************
  void BaseGraph::
  SetBlockInfo( const StdVector<StdVector<UInt> >* indexBlocks ) {
    
    unsortedBlocks_ = indexBlocks;
  }
  
  // **********************
  //   SetRowColDigaGraph
  // **********************
  void BaseGraph::SetRowColDiagGraphs( BaseGraph* rowDiagGraph,
                                       BaseGraph* colDiagGraph ) {
    
    // Check, if graph are defined at all
    if( rowDiagGraph == 0 || colDiagGraph == NULL ) {
      EXCEPTION("row/colDiagGraph object is not set.")
    }
    // check coupling graphs for size
    if( rowDiagGraph->numNodes_ != numNodes_ )  {
      EXCEPTION( "Node mismatch: Row diagonal graph has " << 
                 rowDiagGraph->numNodes_ << " nodes, own graph has " <<
                 numNodes_ );
    }
    
    rowDiagGraph_ = rowDiagGraph;
    colDiagGraph_ = colDiagGraph;
  }
  
  // ********************
  //   ReorderForBlocks  
  // ********************
  void BaseGraph::ReorderForBlocks( StdVector<UInt>& order ) {
    
    LOG_TRACE(graphLogger) << "Calling ReorderForBlocks";
    
    // Check, if we are not yet finalized
    if( amAssembled_ ) {
      EXCEPTION("Can not reorder an already finalised graph.");
    }
    
    // Check, if this graph has set diagonal block graphs, i.e.
    // if this graph is on a non-diagonal position
    if( rowDiagGraph_ != NULL || colDiagGraph_ != NULL ) {
      EXCEPTION("This graph object seems to be a coupling graph. "
                "Thus, no internal reordering according to blocks "
                "can be performed ");
    }
    
    // Resize order
    order.Resize(numNodes_);
    order.Init();
    
    // temporary set for already used indices
    std::set<UInt> usedIndices;
    
    // Loop over all unsorted blocks
    UInt newIndex = 0;
    UInt numBlocks = unsortedBlocks_->GetSize();
    UInt blockStart = 0, blockEnd = 0;
    LOG_DBG(graphLogger) << "graph has " << numBlocks << " blocks";
    sortedBlocks_.Reserve(numBlocks);
    
    for( UInt i = 0; i < numBlocks; ++i ) {
      LOG_DBG3(graphLogger) << "treating block " << i;
      const StdVector<UInt> & actBlock = (*unsortedBlocks_)[i]; 
      UInt blockSize = actBlock.GetSize();
      
      // Note: we can have empty blocks
      if (blockSize == 0 ) 
        continue;
      
      blockStart =  newIndex;
      for( UInt j = 0; j < blockSize; ++j ) {

        LOG_DBG3(graphLogger) << "\tblock-size is " << blockSize;
        // check, if index was already used in different block
        const UInt & index = actBlock[j]; 
        if( usedIndices.find(index) != usedIndices.end() ) {
          EXCEPTION("Index " << index << " was already used in a block.");
        }
        usedIndices.insert(index);
        order[index] = ++newIndex;
        LOG_DBG3(graphLogger) << "\torder[" << index << "] = " << newIndex;
      } // entries in one block
      
      blockEnd = newIndex-1;
      sortedBlocks_.Push_back(std::pair<UInt,UInt>(blockStart,blockEnd));
    } // loop over blocks
    
    // Consistency check: Make sure, that all indices are reordered
    if( numNodes_ != newIndex ) {
      EXCEPTION("BaseGraph::ReorderForBlocks: Graph has " << numNodes_
                << " nodes but there was only a block definition for " 
                << newIndex << " nodes! In case blocks are defined,"
                    "all nodes have to be referenced." );
    }
  }
  

  // *******************
  //   Count non-zeros
  // *******************
  void BaseGraph::CountNNE() {
    UInt i;
    nne_ = 0;
    numNonDiagEntries_ = 0;
    std::vector<unsigned int>::iterator iter;
    for ( i = 0; i < numNodes_; i++ ) {
      nne_ += element_[i].size();
      bool diagFound = false;
      for ( iter = element_[i].begin(); iter != element_[i].end(); iter++ ) {
        if ( *iter == i ) {
          diagFound = true;
          break;
        }
      }
      
      if( !diagFound) {
        numNonDiagEntries_++;
      }

    }
  }

  std::string BaseGraph::ToString() const
  {
    std::stringstream os;

    if (amAssembled_ == false) // uncompressed graph
    {
      os << "Graph is in uncompressed state" << std::endl;
      for(unsigned int i=0; i<numNodes_; i++) {
        os << "line " << i << ", " << element_[i].size() << " entries:\t";
        for (std::vector<unsigned int>::iterator j=element_[i].begin();j!=element_[i].end();j++)
          os << *j << "\t";
        os << std::endl;
      }
    }
    else if (csNodes_ != nullptr) // compressed graph
    {
      os << "Graph is in compressed state" << std::endl;
      for(unsigned int i=0;i<numNodes_;i++) {
        os << "line " << i << ", " << csNodes_[i+1]-csNodes_[i] << " entries:\t";
        for (UInt k=csNodes_[i];k<csNodes_[i+1];k++)
          os << csEdges_[k] << "\t";
        os << std::endl;
      }
    }
    return os.str();
  }

  std::string BaseGraph::PrintCRS(unsigned int numNodes, unsigned int* rows, unsigned int* cols, int level) const
  {
    std::stringstream os;

    switch(level)
    {
    case 0:
      for(unsigned int i = 0; i < numNodes; i++ )
      {
        os << " row " << i << ": ";
        for(unsigned int j = rows[i]; j < rows[i+1]; j++)
          os << cols[j] << " ";
        os << std::endl;
      }
      break;
    case 1:
      os << "rows: " << ::ToString<unsigned int>(rows, numNodes + 1) << " (" << numNodes << "+1)" << std::endl;
      os << "cols: " << ::ToString<unsigned int>(cols, rows[numNodes_]) << " (" << rows[numNodes] << ")";
      break;
    case 2: // metis graph file format (1-based)
      os << numNodes << std::endl;
      for(unsigned int i = 0; i < numNodes; i++ )
      {
        for(unsigned int j = rows[i]; j < rows[i+1]; j++)
          os << cols[j]+1 << " "; // 1-based
        os << std::endl;
      }
      break;
    default:
      EXCEPTION("invalid output level");
    }
    return os.str();
  }


  // ********************
  //   FinaliseAssembly
  // ********************
  void BaseGraph::FinaliseAssembly( BaseOrdering::ReorderingType reorder,  
                                    bool useExternalOrdering,
                                    StdVector<UInt>* vertexOrder,
                                    StdVector<UInt>* edgeOrder  ) {
    this->MapSetToVector();

    newOrder_ = reorder;
    assert(vertexOrder != NULL);
    std::vector<unsigned int>::iterator iter;
    
    // Store requested ordering scheme
    newOrder_ = reorder;
    
    
    // Sort lists and remove duplicate entries
    SortLists();

    // Determine number of non-zero entries
    CountNNE();

    if ( vertexOrder == NULL ) {
      std::string tmp = BaseOrdering::reorderingType.ToString( reorder );

      EXCEPTION("BaseGraph::FinaliseAssembly: I was told to do a '"
          << tmp
          << "' reordering, but you provided a NULL pointer "
          << "for storing the permutation array!"
          << "Who you think you are? I'm mortally offended! ;-)");
    }

    // ==============================
    //  Distinguish reordering cases
    // ==============================
    
    // 1) no external ordering provided : We use the internal reordering strategy
    //    and store the reordering mapping back in the vector vertexOrder
    if( !useExternalOrdering ) {
      
      // check, a block-reordering was set
      if( unsortedBlocks_ != NULL &&
          reorder != BaseOrdering::NOREORDERING ) {
        EXCEPTION( "Can not reorder the graph, as there are sub-blocks defined!");
      }
      
      // resize ordering array
      vertexOrder->Resize(numNodes_);
      edgeOrder = vertexOrder;

      // If user wants, try to reorder the graph, otherwise return identity array
      if ( reorder != BaseOrdering::NOREORDERING ) {

        // compute new ordering
        // in the array "order" we still have vertices starting with one
        // since we pass it back to eqnMap in CFS
        Reorder( reorder, *vertexOrder );
      } else {
        
        // check, if sub-blocks are defined. If yes, we have to reorder
        // the graph to make the blocks continuously numbered. Otherwise we
        // just return the identity mapping.
        if( unsortedBlocks_ != NULL) {
          ReorderForBlocks(*vertexOrder);
        } else {

          // set vertexOrder and edgeOrder accordingly
          vertexOrder->Resize(numNodes_);
          for( UInt i = 0; i < numNodes_; ++i ) {
            (*vertexOrder)[i] = i+1;
          }
        }
      }
      
    }  else {
      // 2) An external ordering (maybe different for vertices and edges)
      //    was provided, so we just take these.
      
      // check if  block-reordering was set
      if( unsortedBlocks_ != NULL ) {
        EXCEPTION( "Can not reorder the graph according to externally provided "
                   "vector, as there are sub-blocks defined!" );
      }
      
      if( edgeOrder == NULL )
        edgeOrder = vertexOrder;
    }
    
    // sort graph according to new ordering
    // step 1: re-number every list
    // step 2: re-arrange lists in vector
    for ( UInt i = 0; i < numNodes_; i++ ) {
      for ( iter = element_[i].begin(); iter != element_[i].end(); iter++ ){
        // minus one, since in OLAS the vertices starts with zero
        *iter = (*edgeOrder)[ *iter ] - 1;
      }
    }

    StdVector<std::vector<unsigned int>> newElement(numNodes_);
    for ( UInt i=0; i< numNodes_; i++ ) {
      UInt n = (*vertexOrder)[i];
      newElement[n-1] = element_[i];
    }

    // Sort the list again (but no need to make unique)
    for ( UInt i=0; i< numNodes_; i++ ) {
      element_[i] = newElement[i];
      std::sort(element_[i].begin(), element_[i].end());
    }


    // Convert Storage format to CRS style
    ConvertToCRS();
    
    // The element vector is no longer required
    delete [] ( element_ );  element_  = NULL;
    element_ = NULL;

    // Now the graph object is fully assembled
    amAssembled_ = true;
 
  }


  // **************
  //   Sort Lists
  // **************
  void BaseGraph::SortLists() {

#ifndef NDEBUG
    UInt totalSize = 0;     // number of entries, including doubles
    UInt totalCapacity = 0; // total capacity 
    boost::accumulators::accumulator_set<Double, stats<boost::accumulators::tag::mean, boost::accumulators::tag::moment<2>, boost::accumulators::tag::sum > > size;
#endif
    // Sort lists and remove duplicate entries 
    for ( UInt i = 0; i < numNodes_; i++ ) {

      // collect information
#ifndef NDEBUG
      totalSize += element_[i].size();
      totalCapacity += element_[i].capacity();
      size(element_[i].size());
#endif
      std::sort( element_[i].begin(), element_[i].end() );
    }

#ifndef NDEBUG
    LOG_DBG(graphLogger) << "SL: statistics";
    LOG_DBG(graphLogger) << "\ttotal number of entries (unsorted): " << totalSize;
    LOG_DBG(graphLogger) << "\ttotal capacity: " << totalCapacity;
    LOG_DBG(graphLogger) << "\tmemory over-estimation: "<< double(totalCapacity/totalSize*100) << " %\n";
    LOG_DBG(graphLogger) << "\taverage row size: " << mean(size) << "\n";
#endif
    
  }


  // ****************
  //   ConvertToCRS
  // ****************
  void BaseGraph::ConvertToCRS()
  {
    // Allocate memory for CRS structure
    if(csNodes_ == nullptr) { // for Metis already called and overwritten reordered
      csNodes_ = new unsigned int[numNodes_+1];
      csEdges_ = new unsigned int[nne_];
    }

    csNodes_[0] = 0;

    std::vector<unsigned int>::iterator iter;
    for(unsigned int i = 0; i < numNodes_; i++ )
    {
      // set number of edges for current node
      csNodes_[i+1] = csNodes_[i] + element_[i].size();

      // set adjacent nodes
      unsigned int j = csNodes_[i];
      for ( iter = element_[i].begin(); iter != element_[i].end(); iter++ ) {
        csEdges_[j] = (*iter);
        j++;
      }
    }

    // LOG_DBG3(graphLogger) << "C2CRS\n" << PrintCRS(numNodes_, csNodes_, csEdges_);
    // LOG_DBG3(graphLogger) << "C2CRS\n" << PrintCRS(numNodes_, csNodes_, csEdges_, 1);
  }



  // *************************
  //   Compute a re-ordering
  // *************************
  void BaseGraph::Reorder( BaseOrdering::ReorderingType newOrder, StdVector<UInt>& order )
  {
    std::string name = BaseOrdering::reorderingType.ToString(newOrder);

    LOG_DBG(graphLogger) << "R: type=" << name << " order capacity=" << order.GetCapacity();

    PtrParamNode pn = domain->GetInfoRoot()->GetList("sequenceStep").Last();
    shared_ptr<Timer> timer = pn->Get("OLAS")->Get("reorder_"  + name + "/timer")->AsTimer();
    timer->SetSub();
    timer->Start();

    UInt i;
    
    switch( newOrder ) {

    case BaseOrdering::METIS:
      {
#ifdef USE_METIS
        // metis can be compiled with configurable layout
        assert(IDXTYPEWIDTH == 32);
        assert(sizeof(int) == sizeof(idx_t)); // idx_t is signed int
        assert(sizeof(int*) == sizeof(idx_t*));

        assert(csNodes_ == nullptr && csEdges_ == nullptr);
        ConvertToCRS();

        // convert CSR mesh to graph
        idx_t* xadj = nullptr; // graph rows
        idx_t* adjncy = nullptr; // graph data
        idx_t numflag = 0; // c-style
        int ret = METIS_MeshToNodal((idx_t*) &numNodes_, (idx_t*) &nne_, (idx_t*) csNodes_, (idx_t*) csEdges_, &numflag, &xadj, &adjncy);
        if(ret != METIS_OK)
          EXCEPTION("cannot convert mesh to metis graph: " << ret);

        // LOG_DBG3(graphLogger) << "R: metis graph\n" << PrintCRS(numNodes_, (unsigned int*) xadj, (unsigned int*) adjncy, 1);

        // Set parameters for Metis
        idx_t options[METIS_NOPTIONS];
        METIS_SetDefaultOptions(options);
        options[METIS_OPTION_NUMBERING] = 0; // C-style
        options[METIS_OPTION_DBGLVL] = 0; // 1 for standard status output

        // Allocate memory for the ordering
        idx_t* perm = new idx_t[numNodes_];
        // for iperm we use the argument order.
        // Note that what for metis is the inverse permutation vector is for us the re-ordering vector.
        assert(order.GetCapacity() == numNodes_);

        ret = METIS_NodeND((idx_t*) &numNodes_, xadj, adjncy, nullptr, options, perm, (idx_t*) order.GetPointer());
        //ret = METIS_NodeND((idx_t*) &numNodes_, (idx_t*) csNodes_, (idx_t*) csEdges_, nullptr, options, perm, (idx_t*) order.GetPointer());
        if(ret != METIS_OK)
          EXCEPTION("cannot reorder metis graph: " << ret);

        // TODO cfs seems to assume 1-based reodering and corrects this. So this shall be set to 0-based also for Sloan!
        for(unsigned int i = 0; i < order.GetSize(); i++)
          order[i]++; // would be nicer if order would be a Vector

        // We don't need the inverse mapping and the graph structure any more
        delete[] perm;
        METIS_Free(xadj); // just a C free
        METIS_Free(adjncy);
        // Generate log message
        // LOG_DBG3(graphLogger) << "R: metis reordering: " << order.ToString();
#else
        EXCEPTION("Re-compile with USE_METIS = ON to get Metis support");
#endif
      }
      break;

    case BaseOrdering::SLOAN:
      {
        MapSetToVector();

        // Generate log message
        LOG_DBG(graphLogger) << " Sloan Reordering";

        // Remove form the neighboring nodes the node itself 
        // ararys are zero based, the equation number however starts by one
        std::vector<unsigned int>::iterator it;
        for ( i = 0; i < numNodes_; i++ ) {

          for ( it = element_[i].begin(); it != element_[i].end(); it++ ) {
            if ( *it == i ) {
              element_[i].erase(it);
              break;
            }
          }
        }

        // put it to one based
        for ( i = 0; i < numNodes_; i++ ) {
          for ( it = element_[i].begin(); it != element_[i].end(); it++ ) {
            *it += 1;
          }
        }

        // Create re-ordering object
        //        Sloan reorder( element_, order, numNodes_ );
        Sloan reorder( element_, order );
        reorder.LabelGraph();

        // Use old graph, if bandwidth was better as after reordering
        Double profOld, profNew;
        reorder.GetProfile( profOld, profNew );

        // Start report
        LOG_DBG(graphLogger) << "\n --> Old Profile: " << profOld  << "\n --> New Profile: " << profNew;

        // If new ordering does not improve the profile, discard it
        if ( profOld <= profNew ) {
          LOG_DBG(graphLogger) << " Discarded re-ordering, since new profile >=" << " old one\n";
          for ( i = 0; i < numNodes_; i++ )
            order[i] = i+1;
        }
        else {
          LOG_DBG(graphLogger) << " --> Profile reduced to " << profNew/profOld * 100.0 << "%\n\n" << " Accepted re-ordering\n";
        }

        // put it to zero based
        for ( i = 0; i < numNodes_; i++ ) {
          for ( it = element_[i].begin(); it != element_[i].end(); it++ ) {
            *it -= 1;
          }
        }
        // re-add the node itself (with the new number)
        for ( i = 0; i < numNodes_; i++ ) {
          element_[i].push_back( i );
          //element_[i].insert( i );
        }
      }
      break;

    default:
      EXCEPTION("Request for unknown re-ordering type");
    }

    // now we are re-ordered
    amReordered_ = true;

    timer->Stop();

  }


  // *********************************************************************
  //   Compute bandwidth (upper and lower) used e.g. for LAPACK Matrices
  // *********************************************************************
  void BaseGraph::GetBandwidth( UInt& bwlower, UInt& bwupper, UInt& bwAvg ) {

                
#ifdef DEBUG_BASEGRAPH
//     if ( amAssembled_ == false ) {
//       EXCEPTION("Attempt to obtain information from graph object, "
//                << "before assembly was completed by calling "
//                << "FinaliseAssembly()");
//     }
#endif
                
    MapSetToVector();

    // Note: If we have a rectangular matrix, this algorithm would fail,
    // as it relies on the presence of a diagonal element.
    if( numRowsMat_ != numColsMat_ ) {
      WARN("GetBandWith() produces wrong results for rectangular matrices.");
      return;
    }
    // If both values are still 0, they have not been computed yet,
    // so we do it now (or we have a diagonal matrix, which is unlikely)
    if ( (bwlower_ == 0) && (bwupper_ == 0) ) {

      UInt bwu, bwl; // distance of current edge
      Integer bwtmp;
      bwavg_ = 0;

      // Loop over all nodes, i.e. all matrix rows
      for ( UInt node = 0; node < numNodes_; node++ ) {

        // since the row is ordered
        bwl = bwu = 0;

        // The row entries are stored with the diagonal entry first,
        // followed by the remaining entries in ascending lexicographic
        // ordering, so the entry with smallest column index is the second
        // one in this row. If the second entry in the row is right of the
        // diagonal, bwl will get < 0, which is okay, but we risk segfault,
        // if we do not test for the case that the row contains only one
        // entry.
        if ( (Integer) csNodes_[node+1] - (Integer) csNodes_[node] > 1 ) {
          // In case of general rectangular matrices, we can not rely on the
          // strict ordering, so we have to check the sign.
          bwtmp = Integer(node - csEdges_[ csNodes_[node] ]);
          if ( bwtmp > 0 )
            bwl = (UInt) bwtmp;
        }

        // The row entries are stored with the diagonal entry first,
        // followed by the remaining entries in ascending lexicographic
        // ordering, so the entry with largest column index is the last
        // one in this row. If the last entry in the row is left of the
        // diagonal, or the diagonal itself, bwu will get <= 0, which is
        // okay.
//std::cout<<"csNodes_[node+1]"<<csNodes_[node+1]<<std::endl;
//std::cout<<"csEdges_[csNodes_[node+1] -1  ]"<<csEdges_[csNodes_[node+1] -1  ]<<std::endl;
//std::cout<<"node"<<node<<std::endl;
        bwtmp = Integer(csEdges_[csNodes_[node+1] -1  ] - node);
        if( bwtmp > 0) {
          bwu = (UInt) bwtmp;
        }
        
        // Sum up for average total bandwidth
        bwavg_ += bwl + bwu + 1;
        
        // Has the lower bandwidth grown?
        bwlower_ = bwl > (UInt)bwlower_ ? (UInt)bwl : bwlower_;

        // Ha the upper bandwidth grown?
        bwupper_ = bwu > (UInt)bwupper_ ? (UInt)bwu : bwupper_;
      }
    }
    bwlower = bwlower_;
    bwupper = bwupper_;
    bwavg_  = UInt( bwavg_ / Double(numNodes_) );
    bwAvg   = bwavg_;

#ifdef DEBUG_BASEGRAPH
    (*debug) << "Bandwidth of the graph: l = " << bwlower_
             << ", u = " << bwupper_ << std::endl;
#endif

  }
  
  // ***************************
  //   Return block definitions
  // ***************************
  StdVector<std::pair<UInt,UInt> >& BaseGraph::GetRowBlockDefinition() {
    
    // Check if graph is already assembled at all
    if( !amAssembled_ ) {
      EXCEPTION("")
    }
    // check, if we have coupling graphs set
    if( rowDiagGraph_ != NULL) {
      return rowDiagGraph_->GetRowBlockDefinition();
    } else {
      return sortedBlocks_;
    }
  }
  
  StdVector<std::pair<UInt,UInt> >& BaseGraph::GetColBlockDefinition() {
    
    // check, if we have coupling graphs set
    if( colDiagGraph_ != NULL) {
      return colDiagGraph_->GetColBlockDefinition();
    } else {
      return sortedBlocks_;
    }
  }

  void BaseGraph::MapSetToVector(){
    //convert set to vector
    if(!setToElemDone_)
    {
      // dynamic scheduling fails for gcc9.
      // it was static,10 before but a quit test showed slight advantages without the option
      #pragma omp parallel for schedule(static) num_threads(CFS_NUM_THREADS)
      for(int i=0;i< (int) numNodes_;i++){
        element_[i].resize(setElements_[i].size()); // TODO: This touches every element - combine with copy.
        std::copy(setElements_[i].begin(), setElements_[i].end(), element_[i].begin());
        setElements_[i].clear();
      }
      delete [] ( setElements_ );  setElements_  = NULL;
      setToElemDone_ = true;
    }
  }


}
