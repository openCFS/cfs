// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;


// ===========================================================================
//
// Variation of the GRAPH-Implementation
//
// 1) Each edge is inserted, regardless if it is already contained in the
//    matrix graph. Fast, but not memory efficient
// 2) For each edge it is checked, if it is already contained in the graph.
//    Slow, but memory efficient
//
// Second variant is used, if GRAPH_VECTOR_SORT macro is enabled. This is
// done by setting FAST_EDGE_INSERTION = no in Makefile.option
//
// ===========================================================================


#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

#include <def_use_metis.hh>

#ifdef USE_METIS

// this namespace serves to hide metis names
namespace olas_metis{
  extern "C"{
#include "metis.h"
  }
}
#endif

#include "graph/basegraph.hh"
#include "graph/sloan.hh"


namespace OLAS {


  // ***************
  //   Constructor
  // ***************
  BaseGraph::BaseGraph( UInt nRows, UInt nCols, ReorderingType reorder ) {

    ENTER_FCN("BaseGraph::BaseGraph");

    // Avoid problems with partially empty graphs
    if ( nRows == 0 || nCols ==0 ) {
      (*error) << "Refusing to generate a graph for a "
               << nRows << " x " << nCols << " matrix!";
      Error( __FILE__, __LINE__ );
    }

    // Initialise attributes
    csEdges_     = NULL;
    csNodes_     = NULL;
    amAssembled_ = false;
    amReordered_ = false;
    newOrder_    = reorder;
    numNodes_    = nRows;
    numRowsMat_  = nRows;
    numColsMat_  = nCols;
    bwlower_     = 0;
    bwupper_     = 0;

    // Allocate memory for linked lists
    if ( numNodes_ > 0 ) {
      NewArray( element_, NodeList, numNodes_ );
    }
    else {
      element_ = NULL;
    }
  }


  // *************************
  //   Alternate Constructor
  // *************************
  BaseGraph::BaseGraph( UInt nRows, UInt nCols, UInt nne, Integer *cs_nodes,
                        Integer *cs_edges, ReorderingType reorder ) :
    amReordered_(false),
    amAssembled_(false) {

    ENTER_FCN( "BaseGraph::BaseGraph" );

    newOrder_    = reorder;
    numNodes_    = nRows;
    numRowsMat_  = nRows;
    numColsMat_  = nCols;
    nne_         = nne;
    csNodes_     = cs_nodes;
    csEdges_     = cs_edges;
    element_     = NULL;
    bwlower_     = 0;
    bwupper_     = 0;
    amAssembled_ = true;
  }


  // *******************
  //   Deep Destructor
  // *******************
  BaseGraph::~BaseGraph() {

    ENTER_FCN( "BaseGraph::~BaseGraph" );

    DeleteArray( csNodes_ );
    DeleteArray( csEdges_ );

    // Note: element_ should already have been deleted
    // in FinaliseAssembly()
    DeleteArray( element_ );

  }


  // ***********************
  //   AddVertexNeighbours
  // ***********************
  void BaseGraph::AddVertexNeighbours( std::vector<UInt> vertexList,
                                       std::vector<UInt> neighbourList ) {

    ENTER_IFCN( "BaseGraph::AddVertexNeighbours" );

    Integer i, j, nodeIndex;

#ifdef DEBUG_BASEGRAPH
    (*debug) << "\n BaseGraph::AddVertexNeighbours:\n";
    (*debug) << " vertexList: ";
    for ( i = 0; i < vertexList.size(); i++ ) {
      (*debug) << vertexList[i] << " ";
    }
    (*debug) << std::endl;
    (*debug) << " neighbourList: ";
    for ( i = 0; i < neighbourList.size(); i++ ) {
      (*debug) << neighbourList[i] << " ";
    }
    (*debug) << '\n' << std::endl;
#endif

#ifdef GRAPH_VECTOR_SORT
    NodeList curList;
    bool edgeFound = false;
    UInt k = 0;
#endif

    // For each vertex in vertexList add all vertices in the neighbourList
    // as neighbours / edges
    for ( i = 0; i < vertexList.size(); i++ ) {
      for ( j = 0; j < neighbourList.size(); j++ ) {

#ifdef GRAPH_VECTOR_SORT

        // Check whether the edge is already contained in the graph
        edgeFound = false;
        curList = element_[ vertexList[i] ];
        for ( k = 0; k < curList.size(); k++ ) {
          if ( curList[k] == neighbourList[j] ) {
            edgeFound = true;
            break;
          }
        }
        if ( edgeFound == false ) {
          element_[ vertexList[i] ].push_back( neighbourList[j] );
        }
#else
        // Insert edge in any case
        element_[ vertexList[i] ].push_back( neighbourList[j] );
#endif
      }
    }
  }


  // *******************
  //   Count non-zeros
  // *******************
  void BaseGraph::CountNNE() {
    ENTER_FCN( "BaseGraph::CountNNE" );
    Integer i;
    nne_ = 0;
    for ( i = 1; i <= numNodes_; i++ ) {
      nne_ += element_[i].size();
    }
  }


  // ***********************************
  //   Print graph to an output stream
  // ***********************************
  void BaseGraph :: Print(std::ostream &os) const {
    ENTER_FCN("BaseGraph::Print");
    Integer i;
    NodeListIterator j;
        
    os << "printing auxiliary graph" << std::endl;

    // uncompressed graph
    if (amAssembled_ == false) {
      os << "Graph is in uncompressed state" << std::endl;
      for (i=1; i<=numNodes_; i++) {
        os << "line " << i << ", " << element_[i].size() << " entries:\t";
        for (j=element_[i].begin();j!=element_[i].end();j++)
          os << *j << "\t";
        os << std::endl;
      }
    }

    // compressed graph
    else if (csNodes_){
      os << "Graph is in compressed state" << std::endl;
      for (i=1;i<=numNodes_;i++) {
        os << "line " << i << ", " << csNodes_[i+1]-csNodes_[i]
           << " entries:\t";
        for (Integer k=csNodes_[i];k<csNodes_[i+1];k++)
          os << csEdges_[k] << "\t";
        os << std::endl;
      }
    }
  }


  // ********************
  //   FinaliseAssembly
  // ********************
  void BaseGraph::FinaliseAssembly( Integer *order ) {

    ENTER_FCN( "BaseGraph::FinaliseAssembly" );

    NodeListIterator iter;

    // Sort lists and remove duplicate entries
    SortLists();

    // Determine number of non-zero entries
    CountNNE();

    // If user wants, try to reorder the graph
    if ( newOrder_ != NOREORDERING ) {

      if ( order == NULL ) {
        (*error) << "BaseGraph::FinaliseAssembly: I was told to do a '"
                 << Enum2String( newOrder_ )
                 << "' reordering, but you provided a NULL pointer "
                 << "for storing the permutation array!"
                 << "Who you think you are? I'm mortally offended! ;-)";
        Error( __FILE__, __LINE__ );
      }

      // compute new ordering
      Reorder( newOrder_, order );

      // sort graph according to new ordering
      // step 1: re-number every list
      // step 2: re-arrange lists in vector
      for ( Integer i = 1; i <= numNodes_; i++ ) {
        for ( iter = element_[i].begin(); iter != element_[i].end(); iter++ ){
          *iter = order[ *iter ];
        }
      }

      NodeList *newElement;
      NewArray( newElement, NodeList, numNodes_ );

      for ( Integer i=1; i<= numNodes_; i++ ) {
        newElement[ order[i] ] = element_[i];
      }

      NodeList *tmpList = NULL;
      tmpList = element_;
      element_ = newElement;
      newElement = NULL;
      DeleteArray( tmpList );

      // Sort the list again (but no need to make unique)
      for ( Integer i = 1; i <= numNodes_; i++ ) {
        std::sort(element_[i].begin(), element_[i].end());
      }

#ifdef DEBUG_BASEGRAPH
      (*debug) << "Reordering: New mapping" << std::endl;
      for ( Integer i = 1; i <= numNodes_; i++ ) {
        (*debug) << i << " -> " << order[i] << std::endl;
      }
#endif

    }

    // Convert Storage format to CRS style
    ConvertToCRS();
    
    // The element vector is no longer required
    DeleteArray( element_ );
    element_ = NULL;

    // Now the graph object is fully assembled
    amAssembled_ = true;
 
  }


  // **************
  //   Sort Lists
  // **************
  void BaseGraph::SortLists() {

    ENTER_FCN( "BaseGraph::SortLists" );

    // Sort lists and remove duplicate entries 
    for ( Integer i = 1; i <= numNodes_; i++ ) {
      std::sort( element_[i].begin(), element_[i].end() );

      NodeList::iterator new_end = std::unique( element_[i].begin(), 
						element_[i].end() );
      NodeList temp( element_[i].begin(), new_end );
      element_[i] = temp;
    }

  }


  // ****************
  //   ConvertToCRS
  // ****************
  void BaseGraph::ConvertToCRS() {

    ENTER_FCN( "BaseGraph::ConvertToCRS" );

    UInt i, j;

    // Allocate memory for CRS structure
    if ( csNodes_ == NULL ) {
      NewArray( csNodes_, Integer, numNodes_ + 1 );
      NewArray( csEdges_, Integer, nne_ );
    }

    csNodes_[1] = 1;

    NodeListIterator iter;
    UInt pos;

    for ( i = 1; i <= numNodes_; i++ ) {

      // set number of edges for current node
      csNodes_[i+1] = csNodes_[i] + element_[i].size();

      // set adjacent nodes
      j = csNodes_[i];
      for ( iter = element_[i].begin(); iter != element_[i].end(); iter++ ) {
        csEdges_[j] = (*iter);
        j++;
      }
    }

#ifdef DEBUG_BASEGRAPH
    (*debug) << "Graph in CRS format: " << std::endl;
    for ( i = 1; i <= numNodes_; i++ ) {
      (*debug) << " Row " << i << ": ";
      for ( j = csNodes_[i]; j < csNodes_[i+1]; j++ ) {
        (*debug) << csEdges_[j] << " ";
      }
      (*debug) << std::endl;
    }
#endif

  }


  // *********************
  //   ConvertToMetisCRS
  // *********************
  void BaseGraph::ConvertToMetisCRS( Integer **rptr, Integer **cidx ) {

    ENTER_FCN( "BaseGraph::ConvertToMetisCRS" );

    UInt i, j;

    // Make sure nne_ is up-to-date
    CountNNE();

    // Allocate memory for CRS structure
    NewArray( (*rptr), Integer, numNodes_ + 1 );
    NewArray( (*cidx), Integer, nne_ - numNodes_ );

    (*rptr)[1] = 1;

    NodeListIterator iter;
    Integer pos;

    for ( i = 1; i <= numNodes_; i++ ) {

      // set number of edges for current node (minus self-reference of node)
        
      // insert edge information (minus self-reference of node)
      j = (*rptr)[i];
      for ( iter = element_[i].begin(); iter != element_[i].end(); iter++ ) {
        pos = (*iter);
        if ( pos != i ) {
          (*cidx)[j] = pos;
          j++;
        }
      }

      // Set start of next CRS row
      (*rptr)[i+1] = j;

    }
        
#ifdef DEBUG_BASEGRAPH
    (*debug) << "Graph in CRS format (minus self references): " << std::endl;
    for ( i = 1; i <= numNodes_; i++ ) {
      for ( j = (*rptr)[i]; j < (*rptr)[i+1]; j++ ) {
        (*debug) << (*cidx)[j] << " ";
      }
      (*debug) << std::endl;
    }
#endif

  }


  // *************************
  //   Compute a re-ordering
  // *************************
  void BaseGraph::Reorder( ReorderingType newOrder, Integer *order ) {

    ENTER_FCN( "BaseGraph::Reorder" );
    Integer i;
    
    // Test, if there is memory available for storing new ordering
    if ( order == NULL ) {
      (*error) << "No memory provided for storing re-ordering";
      Error( __FILE__, __LINE__ );
    }

    switch( newOrder ) {

    case METIS:
      {
#ifdef USE_METIS

        //if ( !amAssembled_ ) {
        //  ConvertToCRS();
        //}

        // Generate auxilliary CRS representation of the graph for
        // Metis. We cannot use ConvertToCRS() here, since this
        // contains self-references of the nodes
        Integer *rptr = NULL;
        Integer *cidx = NULL;
        ConvertToMetisCRS( &rptr, &cidx );

        // Allocate memory for the ordering
        Integer *iorder = NULL;
        NewArray( iorder, Integer, numNodes_ );

        // Set parameters for Metis
        Integer options[8];
        options[0] = 0; // (use defaults)

        // options[1] = 3;  // default
        // options[2] = 1;  // default
        // options[3] = 2;  // default
        // options[4] = 0;  // default
        // options[5] = 1;  // default
        // options[6] = 0;  // default
        // options[7] = 5;

        Integer numflag = 1;
        int nNodes = (int)numNodes_;

        // Note that what for metis is the inverse permutation vector,
        // for us is the re-ordering vector
        olas_metis::METIS_NodeND( &nNodes, rptr + 1, cidx + 1, &numflag,
                                  options, iorder + 1, order + 1 );

        // We don't need the inverse mapping and the auxilliary CRS structure
        // anymore, so free the memory
        DeleteArray( iorder );
        DeleteArray( rptr );
        DeleteArray( cidx );

        // Generate log message
        (*cla) << " Reordering: Re-ordered graph using Metis" << std::endl;

#else
        (*error) << "Re-compile with USE_METIS = yes to get Metis support";
        Error( __FILE__, __LINE__ );
#endif
      }
      break;

    case SLOAN:
      {

        // Generate log message
        (*cla) << " -----------------------------------------------\n"
               << " Sloan Reordering:"
               << std::endl;

        // Remove form the neighboring nodes the node itself 
        NodeList::iterator it;
        for ( i = 1; i <= numNodes_; i++ ) {

          for ( it = element_[i].begin(); it != element_[i].end(); it++ ) {
            
            if ( *it == i ) {
              element_[i].erase(it);
              break;
            }
          }

        }

        // Create re-ordering object
        Sloan reorder( element_, order, numNodes_ );
        reorder.LabelGraph();

        // Use old graph, if bandwidth was better as after reordering
        Double profOld, profNew;
        reorder.GetProfile( profOld, profNew );

        // Start report
        (*cla) << "\n --> Old Profile: " << profOld
               << "\n --> New Profile: " << profNew
               << std::endl;

        // If new ordering does not improve the profile, discard it
        if ( profOld <= profNew ) {

          (*cla) << " Discarded re-ordering, since new profile >="
                 << " old one\n";

          for ( i = 1; i <= numNodes_; i++ ) {
            order[i] = i;
          }
        }

        else {
          (*cla) << " --> Profile reduced to "
                 << profNew/profOld * 100.0 << "%\n\n"
                 << " Accepted re-ordering\n";
        }

        // Finish report
        (*cla) << " -----------------------------------------------"
               << std::endl;

        // re-add the node itself (with the new number)
        for ( i = 1; i <= numNodes_; i++ ) {
          element_[i].push_back( i );
          //element_[i].insert( i );
        }
      }
      break;

    case NOREORDERING:
      // nothing to do
      break;

    default:
      (*error) << "Request for unknown re-ordering type";
      Error( __FILE__, __LINE__ );
    }

    // now we are re-ordered
    amReordered_ = true;

  }


  // *********************************************************************
  //   Compute bandwidth (upper and lower) used e.g. for LAPACK Matrices
  // *********************************************************************
  void BaseGraph::GetBandwidth( UInt& bwlower, UInt& bwupper ) {

    ENTER_FCN( "BaseGraph::GetBandwidth" );
                
#ifdef DEBUG_BASEGRAPH
    if ( amAssembled_ == false ) {
      (*error) << "Attempt to obtain information from graph object, "
               << "before assembly was completed by calling "
               << "FinaliseAssembly()";
      Error( __FILE__, __LINE__ );
    }
#endif
                
    // If both values are still 0, they have not been computed yet,
    // so we do it now (or we have a diagonal matrix, which is unlikely)
    if ( (bwlower_ == 0) && (bwupper_ == 0) ) {

      Integer bwu, bwl; // distance of current edge

      // Loop over all nodes, i.e. all matrix rows
      for ( UInt node = 1; node <= numNodes_; node++ ) {

        // since the row is ordered
        bwl = bwu = 0;

        // The row entries are stored with the diagonal entry first,
        // followed by the remaining entries in ascending lexicographic
        // ordering, so the entry with smalles column index is the second
        // one in this row. If the second entry in the row is right of the
        // diagonal, bwl will get < 0, which is okay, but we risk segfault,
        // if we do not test for the case that the row contains only one
        // entry.
        if ( csNodes_[node+1] - csNodes_[node] > 1 ) {
          bwl = (Integer)node - csEdges_[ csNodes_[node] + 1 ];
        }

        // The row entries are stored with the diagonal entry first,
        // followed by the remaining entries in ascending lexicographic
        // ordering, so the entry with largest column index is the last
        // one in this row. If the last entry in the row is left of the
        // diagonal, or the diagonal itself, bwu will get <= 0, which is
        // okay.
        bwu = csEdges_[ csNodes_[node+1] - 1 ] - (Integer)node;
           
        // Has the lower bandwidth grown?
        bwlower_ = bwl > (Integer)bwlower_ ? (UInt)bwl : bwlower_;

        // Ha the upper bandwidth grown?
        bwupper_ = bwu > (Integer)bwupper_ ? (UInt)bwu : bwupper_;
      }
    }

    bwlower = bwlower_;
    bwupper = bwupper_;

#ifdef DEBUG_BASEGRAPH
    (*debug) << "Bandwidth of the graph: l = " << bwlower_
             << ", u = " << bwupper_ << std::endl;
#endif

  }

}
