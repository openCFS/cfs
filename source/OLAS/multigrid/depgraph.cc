// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

/* $Id$ */

#include <string.h>

#include "multigrid/depgraph.hh"
#include "DataInOut/coloredConsole.hh"

/**********************************************************/
#ifdef DEBUG_TO_CERR
#ifndef DEBUG_DEPENDENCYGRAPH
#define DEBUG_DEPENDENCYGRAPH
#endif // DEBUG_DEPENDENCYGRAPH
#define  debug  &std::cerr
#endif // DEBUG_TO_CERR
/**********************************************************/

namespace CoupledField {
/**********************************************************/

template <typename T>
DependencyGraph<T>::DependencyGraph()
                  : NumNodes_(0),
                    NumEdges_(0),
                    StartIndex_(NULL),
                    NodeSize_(NULL),
                    Edges_(NULL),
                    ownStartIndex_(false)
{
}

/**********************************************************/

template <typename T>
DependencyGraph<T>::~DependencyGraph()
{

    Reset();
}

/**********************************************************/

template <typename T>
bool DependencyGraph<T>::Create( const Integer num_nodes,
                                 const Integer num_edges_per_node )
{

    // create arrays, in particular also start indices
    if( !CreateArrays(num_nodes,
                      num_edges_per_node * num_nodes,
                      true) )  return false;

    // initialize the start indices and the node sizes
    for( Integer i = 1; i <= num_nodes; i++ ) {
        StartIndex_[i] = (i-1)*num_edges_per_node + 1;
        NodeSize_[i]   = 0;
    }
    // set the last start index for the virtual node [n+1]
    StartIndex_[num_nodes+1] = num_nodes*num_edges_per_node + 1;
    
    return true;
}

/**********************************************************/

template <typename T>
bool DependencyGraph<T>::
CreateWithDiagonals( const Integer num_nodes,
                     const Integer num_edges_per_node )
{

    // create arrays, in particular also start indices
    if( !CreateArrays(num_nodes,
                      num_edges_per_node * num_nodes,
                      true) )  return false;

    // initialize the start indices and the node sizes
    Integer startindex = 1;
    for( Integer i = 1; i <= num_nodes; i++ ) {
        StartIndex_[i]     = startindex;
        Edges_[startindex] = i;
        NodeSize_[i]       = 1;
        startindex        += num_edges_per_node;
    }
    // set the last start index for the virtual node [n+1]
    StartIndex_[num_nodes+1] = startindex;
    
    return true;
}

/**********************************************************/

template <typename T>
bool DependencyGraph<T>::Create( const Integer        num_nodes,
                                 const Integer *const node_size,
                                       Integer        num_edges )
{

    // get total number of edges, if not passed
    if( num_edges < 0 ) {
        num_edges = 0;
        for( Integer i = 1; i <= num_nodes; i++ ) {
            num_edges += node_size[i];
        }
    }
    // in debug mode check the correctness of the passed
    // maximal total number of edges
#ifdef DEBUG_DEPENDENCYGRAPH
    else {
        int calculated_num_edges = 0;
        for( Integer i = 1; i <= num_nodes; i++ ) {
            calculated_num_edges += node_size[i];
        }
        if( calculated_num_edges != num_edges ) {
            Warning( "passed wrong maximal number of edges, proceeding "
                     "with the correct value", __FILE__, __LINE__ );
            num_edges = calculated_num_edges;
        }
    }
#endif

    // create arrays, including start indices
    if( !CreateArrays(num_nodes, num_edges, true) )  return false;
    // fill array with start indices
    num_edges = StartIndex_[1] = 1;
    for( Integer i = 1; i <= num_nodes; i++ ) {
        StartIndex_[i+1] = (num_edges += node_size[i]);
    }    
    for( Integer i = 1; i <= num_nodes; i++ )  NodeSize_[i] = 0;

    return true;
}

/**********************************************************/

template <typename T>
bool DependencyGraph<T>::
Create( const Integer        num_nodes,
        const Integer *const startarray,
        const bool           accroach_startarray,
        const bool           copy_startarray )
{

    if( num_nodes <= 0 ) { Reset(); return true; }

#ifdef DEBUG_DEPENDENCYGRAPH
    if( accroach_startarray && copy_startarray ) {
        Warning( "The combination of accroach_startarray == true "
                 "and copy_startarray == true does not make sense",
                 __FILE__, __LINE__ );
    }
#endif

    // create arrays
    if( !CreateArrays(num_nodes,
                      startarray[num_nodes+1],
                      copy_startarray) ) {
        return false;
    }
    
    // insert startarray
    if( !InsertStartArray(startarray,          // array with start indices
                          num_nodes,           // number of nodes
                          accroach_startarray, // pipe parameter
                          copy_startarray,     // pipe parameter
                          true) ) {            // force insertion of start array
        return false;
    }
    
    // initialize node sizes (#edges / node[i])
    for( Integer i = 1; i <= num_nodes; i++ )  NodeSize_[i] = 0;

    return true;
}

/**********************************************************/

template <typename T>
bool DependencyGraph<T>::Create( const CRS_Matrix<T>& matrix,
                                 const bool buildemptygraph,
                                 const bool copystartarray  )
{

    if( matrix.GetNnz() == 0 ) { Reset(); return true; }

    // create the arrays
    if( !CreateArrays(matrix.GetNumRows(),
                      matrix.GetNnz(),
                      copystartarray) ) {
        return false;
    }    
    // insert the array of start indices
    if( !InsertStartArray(matrix.GetRowPointer(), // pointer to start indices
                          matrix.GetNumRows(), // #rows == #nodes
                          false, // never accroach an array from a matrix
                          copystartarray, // just pipe the parameter
                          true) ) { // important, if copystartarray == true
        return false;
    }

    // initialize number of edges
    if( buildemptygraph ) {
        for( int i = 1; i <= NumNodes_; i++ )  NodeSize_[i] = 0;
    // copy matrix graph
    } else {
        for( int i = 1; i <= NumNodes_; i++ ) {
            NodeSize_[i] = matrix.GetRowSize(i);
        }
        memcpy( Edges_+1, matrix.GetColPointer()+1, NumEdges_*sizeof(Integer) );
    }

    return true;
}

/**********************************************************/

template <typename T>
bool DependencyGraph<T>::
InsertStartArray( const Integer *const startarray,
                  const Integer        size,
                  const bool           accroach,
                  const bool           copy,
                  const bool           force )
{
    
    // if there is already an array with start indices
    if( StartIndex_ ) {
        // only overwrite, if explicitly demanded
        if( false == force )  return false;
        // if this object owns the start array
        if( ownStartIndex_ ) {
            // check, if the sizes match
            if( size != NumNodes_ ) {
                // cannot change number of nodes, if connections already exist
                if( NumEdges_ )  return false;
                delete [] ( StartIndex_ );  StartIndex_  = NULL;
                NEWARRAY( StartIndex_, Integer, size+1 );
            }
        }
    }
    // insert new array
    if( copy )  memcpy( StartIndex_+1, startarray+1, (size+1)*sizeof(Integer) );
    else        StartIndex_ = ((Integer*)startarray);

    NumNodes_      = size;
    ownStartIndex_ = accroach || copy;

    return true;    
}

/**********************************************************/

template <typename T>
inline void DependencyGraph<T>::AddEdge( const Integer i,
                                         const Integer j )
{
    
#ifdef DEBUG_DEPENDENCYGRAPH
    // check range of the parameters
    if( i < 1 || i > NumNodes_  ) {
        std::cerr
          << fg_red << "DependencyGraph<T>::AddEdge" << fg_reset << "\n"
           "  >> edge ("<<i<<","<<j<<") is out of range [1,"
        << NumNodes_ << "] for node index" << std::endl;
        return;
    }
    // check, if this edge is already present
    for( Integer ij = 0; ij < NodeSize_[i]; ij++ ) {
        if( Edges_[StartIndex_[i]+ij] == j ) {
            std::cerr
            << fg_red << "DependencyGraph<T>::AddEdge" << fg_reset
            << "  >> edge ("<<i<<","<<j<<") already exists.\n     "
               "To add edges savely in non-debug mode use "
               "AddEdgeSavely\n";
            return;                                                                  
        }
    }
    // check, if the node can take further edges
    if( NodeSize_[i] >= StartIndex_[i+1] - StartIndex_[i] ) {
        std::cerr
        << fg_red << "DependencyGraph<T>::AddEdge" << fg_reset << "\n  >> "
           "node " << i << " cannot take the additional edge ("
        << i<<","<<j<<") (max " << NodeSize_[i] << ")";
        if( GetNumEdges(i) > 0 ) {
            for( Integer k = 0; k < GetNumEdges(i); k++ ) {
                std::cerr << GetEdges(i)[k];
                if( k < GetNumEdges(i)-1 )  std::cerr << ", ";
            }
        }
        std::cerr << std::endl;
        return;
    }
#endif

    // add the new edge
    Edges_[StartIndex_[i] + NodeSize_[i]++] = j;
}

/**********************************************************/

template <typename T>
inline void DependencyGraph<T>::AddEdgeSavely( const Integer i,
                                               const Integer j )
{
    
#ifdef DEBUG_DEPENDENCYGRAPH
    // check range of the parameters
    if( i < 1 || i > NumNodes_ ) {
        std::cerr << fg_red << "DependencyGraph<T>::AddEdgeSavely" << fg_reset
                     "  >> edge ("<<i<<","<<j<<") is out of range "
                  << NumNodes_ << std::endl;
        return;
    }
#endif

    const Integer start = StartIndex_[i];
    const Integer end   = start + NodeSize_[i];
    // search for the edge in the graph
    for( Integer ij = start; ij < end; ij++ ) {
        // edge already present?, job done
        if( Edges_[ij] == j )  return;
    }

#ifdef DEBUG_DEPENDENCYGRAPH
    // check, if the node can take further edges
    if( NodeSize_[i] >= StartIndex_[i+1] - StartIndex_[i] ) {
        std::cerr << fg_red << "DependencyGraph<T>::AddEdgeSavely" << fg_reset << "\n  >> "
                     "node " << i << " cannot take further edges (maximal "
                  << NodeSize_ << ")\n";
        return;
    }
#endif

    // add the new edge
    Edges_[start + NodeSize_[i]++] = j;
}

/**********************************************************/

template <typename T>
inline void DependencyGraph<T>::
AddEdgeSorted( const Integer i,
               const Integer j,
               const Integer first_position )
{

    // simply append the edge
    if( first_position >= NodeSize_[i] ) {
#ifdef DEBUG_DEPENDENCYGRAPH
        // check, if the node can take further edges
        if( NodeSize_[i] >= StartIndex_[i+1] - StartIndex_[i] ) {
            std::cerr << fg_red << "DependencyGraph<T>::AddEdgeSorted" << fg_reset << "\n  >> "
                         "node " << i << " cannot take further edges (maximal "
                      << NodeSize_ << ")\n";
            return;
        }
#endif
        Edges_[StartIndex_[i]] = j;
        NodeSize_[i]++;
        return;
    }

    Integer        first = first_position,
                   last  = NodeSize_[i] - 1;
    Integer *const edges = Edges_ + StartIndex_[i];

    // search the position of the new edge via bisection
    while( true ) {
        // found edge?
        if( edges[first] == j || edges[last] == j )  return;
        // add edge (i,j) BEFORE position [first]
        if( j < edges[first] ) {
            for( Integer ij = NodeSize_[i]; ij > first; ij-- ) {
                edges[ij] = edges[ij-1];
            }
#ifdef DEBUG_DEPENDENCYGRAPH
            // check, if the node can take further edges
            if( NodeSize_[i] >= StartIndex_[i+1] - StartIndex_[i] ) {
                std::cerr << fg_red << "DependencyGraph<T>::AddEdgeSorted" << fg_reset << "\n  >> "
                             "node " << i << " cannot take further edges (maximal "
                          << NodeSize_ << ")\n";
                return;
            }
#endif
            edges[first] = j;
            NodeSize_[i]++;
            return;
        }
        // add edge (i,j) AFTER position [last]
        if( edges[last] < j ) {
            for( Integer ij = NodeSize_[i]-1; ij > last; ij-- ) {
                edges[ij+1] = edges[ij];
            }
#ifdef DEBUG_DEPENDENCYGRAPH
            // check, if the node can take further edges
            if( NodeSize_[i] >= StartIndex_[i+1] - StartIndex_[i] ) {
                std::cerr << fg_red << "DependencyGraph<T>::AddEdgeSorted" << fg_reset << "\n  >> "
                             "node " << i << " cannot take further edges (maximal "
                          << NodeSize_ << ")\n";
                return;
            }
#endif
            edges[last+1] = j;
            NodeSize_[i]++;
            return;
        }
        // NOTE: at this position we cannot get first == last, because
        // this would mean with j_edge := edge[last] = edge[first],
        //                  ==> j_edge <= j <= j_edge
        //                  ==> j_edge == j
        // but this is checked on the top of the loop body
        // ==> bisection
        if( first+1 == last ) {
            if( j < edges[last] )  first = last;
            else                   last  = first;
        } else {
            const Integer middle = (first + last) >> 1;
            if( j < edges[middle] )  last  = middle;
            else                     first = middle;
        }
    }
}
    
/**********************************************************/

template <typename T>
inline void DependencyGraph<T>::
AddEdgeSortedAfterDiag( const Integer i, const Integer j )
{
    if( i == j )  return;
    AddEdgeSorted( i, j, 1 );
}

/**********************************************************/

template <typename T>
inline void DependencyGraph<T>::SetEdgeAtPosition( const Integer i,
                                                   const Integer j,
                                                   const Integer position )
{
    
#ifdef DEBUG_DEPENDENCYGRAPH
    const char fn[] = "DependencyGraph::AddEdgeAtPosition\n";
    // check range of the node number
    if( i < 1 || i > NumNodes_ ) {
        std::cerr << fg_red << fn << fg_reset << "  >> edge ("<<i<<","<<j<<") is out of range "
                  << NumNodes_ << std::endl;
        return;
    }
    // check range of position
    if( StartIndex_[i] + position >= StartIndex_[i+1] ) {
        std::cerr << fg_red << fn << fg_reset << "  >> position " << position << " ist out of "
                     "range [0,"<<GetMaxNumEdges(i)<<"] for node "<<i
                  << std::endl;
        return;
    }
    // check if the edge is already present
    for( int ij = StartIndex_[i]; ij < StartIndex_[i+1]; i++ ) {
        if( Edges_[ij] == j && StartIndex_[i] + position != ij ) {
            std::cerr << fg_red << fn << fg_reset << "  >> edge ("<<i<<","<<j<<") is already"
                         " present at position "<<(ij - StartIndex_[i])
                      << " at node "<<i<<". Is this intended?";
        }
    }
#endif

    Edges_[StartIndex_[i] + position] = j;
    if( position >= NodeSize_[i] )  NodeSize_[i] = position + 1;
}

/**********************************************************/

template <typename T>
inline void DependencyGraph<T>::RemoveEdge( const Integer i,
                                            const Integer j )
{
    
#ifdef DEBUG_DEPENDENCYGRAPH
    // check range of the parameters
    if( i < 1 || i > NumNodes_ ) {
        std::cerr << fg_red << "DependencyGraph<T>::RemoveEdge" << fg_reset << "\n"
                     "  >> edge ("<<i<<","<<j<<") is out of range "
                  << NumNodes_ << std::endl;
        return;
    }
#endif

    const Integer start = StartIndex_[i];
    const Integer end   = start + NodeSize_[i];

    // seek the edge
    for( Integer ij = start; ij < end; ij++ ) {
        // if edge is present, ...
        if( Edges_[ij] == j ) {
            // ... remove it by overwriting it with the last edge
            Edges_[ij] = Edges_[end-1];
            NodeSize_[i]--;
            return;
        }
    }
}

/**********************************************************/

template <typename T>
inline void DependencyGraph<T>::RemoveAllEdges( const Integer i )
{
    
    NodeSize_[i] = 0;
}

/**********************************************************/

template <typename T>
inline Integer DependencyGraph<T>::
ShiftEdges( const Integer offset )
{

    Integer max_j = 0;
    
    for( Integer i = 1; i <= NumNodes_; i++ ) {
        Integer* const edges = Edges_ + StartIndex_[i];
        for( Integer ij = 0; ij < NodeSize_[i]; ij++ ) {
            if( max_j < (edges[ij] += offset) ) {
                max_j = edges[ij];
            }
        }
    }

    return max_j;
}

/**********************************************************/

template <typename T>
inline void DependencyGraph<T>::Sort()
{
    for( Integer i = 1; i <= NumNodes_; i++ ) {
        QuickSort( Edges_ + StartIndex_[i], NodeSize_[i] );
    }
}

/**********************************************************/

template <typename T>
inline void DependencyGraph<T>::SortDiagFirst()
{
    for( Integer i = 1; i <= NumNodes_; i++ ) {
        Integer* const edges = Edges_ + StartIndex_[i];
        // first place diagonal edge at first position
        for( Integer j = 0; j < NodeSize_[i]; j++ ) {
            if( edges[j] == i ) {
                edges[j] = edges[0];
                edges[j] = i;
                break;
            }
        }
        // sort the rest
        QuickSort( edges+1, NodeSize_[i]-1 );
    }
}

/**********************************************************/

template <typename T>
bool DependencyGraph<T>::
AssignTransposed( const DependencyGraph<T>& graph,
                  bool    use_start_array,
                  bool    treat_as_square,
                  Integer overlap_per_node )
{

    if( use_start_array == true ) {
#ifdef DEBUG_DEPENDENCYGRAPH
        if( treat_as_square == false ) {
            Warning( "DependencyGraph::AssignTransposed: parameter "
                     "\"treat_as_square\" == false does not make sense "
                     "in combination with \"use_start_array\" == true "
                     "-> \"treat_as_square\" is set to true." );
        }
        if( overlap_per_node > 1 ) {
            Warning( "DependencyGraph::AssignTransposed: parameter "
                     "\"overlap_per_node\" > 0 does not make sense in "
                     "combination with \"use_start_array\" == true -> "
                     "\"overlap_per_node\" will be ignored." );
        }
#endif // DEBUG_DEPENDENCYGRAPH
        treat_as_square = true;
    }    

    // NOTE: Although this assignment might cause a total rebuild
    //       of the graph, we need not call Reset() here, since
    //       the call of CreateArrays(..) will check if we can
    //       reuse the old memory.

    // get number of nodes and number of edges from the source graph
    // NOTE that we will create arrays, in particular the edge array,
    // in that the edges fit exactly, so that there cannot be added
    // further edges in a transposed graph.
    const Integer newNumEdges = graph.GetNumEdges();
          Integer newNumNodes = graph.GetNumNodes();

    // calculate the new number of nodes, if not explicitely
    // specified that the graph's edges can be represented
    // by a square matrix
    if( treat_as_square == false ) {
        for( Integer i = 1; i <= graph.GetNumNodes(); i++ ) {
            const Integer* const  edges = graph.GetEdges( i );
            const Integer        nedges = graph.GetNumEdges( i );
            for( Integer ij = 0; ij < nedges; ij++ ) {
                if( edges[ij] > newNumNodes ) {
                    newNumNodes = edges[ij];
                }
            }
        }
    }
#ifdef DEBUG_DEPENDENCYGRAPH
    // in debug mode always check the edge range
    else {
        Integer debug_highestEdgeIndex = graph.GetNumNodes(),
                debug_criticalNode     = 0;
        for( Integer i = 1; i <= graph.GetNumNodes(); i++ ) {
            const Integer* const edges = graph.GetEdges(i);
            for( Integer ij = 0; ij < graph.GetNumEdges(i); ij++ ) {
                if( edges[ij] > debug_highestEdgeIndex ) {
                    debug_highestEdgeIndex = edges[ij];
                    debug_criticalNode     = i;
                }
            }
        }
        if( debug_highestEdgeIndex > graph.GetNumNodes() ) {
            const char msgf[] = "DependencyGraph<T>::AssignTransposed: "
                "edge checking suppressed and predicted graph as square"
                ", but there are %d nodes in the graph, and the edge "
                "with highest j is (i,j) = (%d,%d) -> (this check is "
                "only performed in debug mode)\n";
            char *msg = NULL;
            NEWARRAY( msg, char, strlen(msgf)+50 )
            sprintf( msg+1, msgf, graph.GetNumNodes(),
                     debug_criticalNode, debug_highestEdgeIndex );
            Warning( msg+1, __FILE__, __LINE__ );
            return false;
        }
    }
#endif // DEBUG_DEPENDENCYGRAPH

    // create arrays
    // If we share the start array with the original graph, ...
    if( use_start_array ) {
        // ... we create the arrays ...
        if( CreateArrays(newNumNodes,            // for newNumNodes nodes
                         graph.GetMaxNumEdges(), // maximal newNumEdges edges,
                         false) // but without the start index array.
            &&
            // Insert the start index array of the original graph:
            InsertStartArray(graph.StartIndex_,
                             graph.GetNumNodes(),
                             false,    // do not accroach the array
                             false,    // do not copy the array
                             true) ) { // force the insertion
            // initialize the node size array
            for( Integer i = 1; i <= GetNumNodes(); i++ ) {
                NodeSize_[i] = 0;
            }
        } else {
            return false;
        }
    // If we should create our own start index array ...
    } else {
        // create the arrays including the start index array.
        if( CreateArrays(newNumNodes,
                         newNumEdges + newNumNodes*overlap_per_node,
                         true) ) {
            // initialize arrays
            for( Integer i = 1; i <= graph.GetNumNodes(); i++ ) {
                NodeSize_[i] = 0;
            }
            // get the node sizes of the transposed graph
            for( Integer i = 1; i <= graph.NumNodes_; i++ ) {
                const Integer* const  edges = graph.GetEdges( i );
                const Integer        nedges = graph.GetNumEdges( i );
                for( Integer ij = 0; ij < nedges; ij++ ) {
                    NodeSize_[edges[ij]]++;
                }
            }
            // adapt the start index array
            StartIndex_[1] = 1;
            for( Integer i = 1; i <= NumNodes_; i++ ) {
                StartIndex_[i+1] =   StartIndex_[i]
                                   + NodeSize_[i]
                                   + overlap_per_node;
                NodeSize_[i] = 0;
            }
        } else {
            return false;
        }
    }

    // insert edges as transposed edges
    for( Integer i = 1; i <= graph.GetNumNodes(); i++ ) {
        const Integer* const  edges = graph.GetEdges( i );
        const Integer        nedges = graph.GetNumEdges( i );
        for( Integer ij = 0; ij < nedges; ij++ ) {
            AddEdge( edges[ij], i );
        }
    }

    return true;
}

/**********************************************************/

template <typename T>
bool DependencyGraph<T>::
CheckTransposition( const DependencyGraph<T>& graph ) const
{

    bool result = true;

    // check number of nodes
    if( GetNumNodes() != graph.GetNumNodes() ) {
#ifdef DEBUG_DEPENDENCYGRAPH
        (*debug) << "DependencyGraph::CheckTransposition: "
                 << "number of nodes did not match " << GetNumNodes()
                 << " <-> " << graph.GetNumNodes() << std::endl;
#endif
        return false;
    }

    // check edges in passed graph to be represented in this one
    for( Integer i = 1; i <= graph.GetNumNodes(); i++ ) {
        const Integer        nedges = graph.GetNumEdges( i );
        const Integer *const  edges = graph.GetEdges( i );
        for( Integer ij = 0; ij < nedges; ij++ ) {
            if( false == IsElement(edges[ij], i) ) {
#ifdef DEBUG_DEPENDENCYGRAPH
                (*debug) << "DependencyGraph::CheckTransposition: " << std::endl
                         << "edge ("<<i<<","<<edges[ij]<<") in the passed "
                            "graph has no corresponding edge in this graph"
                         << std::endl;
#endif
                result = false;
            }
        }
    }
    // and vice versa
    for( Integer i = 1; i <= GetNumNodes(); i++ ) {
        const Integer        nedges = GetNumEdges( i );
        const Integer *const  edges = GetEdges( i );
        for( Integer ij = 0; ij < nedges; ij++ ) {
            if( false == graph.IsElement(edges[ij], i) ) {
#ifdef DEBUG_DEPENDENCYGRAPH
                (*debug) << "DependencyGraph::CheckTransposition: " << std::endl
                         << "edge ("<<i<<","<<edges[ij]<<") in this graph "
                            "has no corresponding edge in the passed graph"
                         << std::endl;
#endif
                result = false;
            }
        }
    }

    return result;
}

/**********************************************************/

template <typename T>
inline bool DependencyGraph<T>::IsElement( const Integer i,
                                           const Integer j ) const
{
    
    const Integer end = StartIndex_[i] + NodeSize_[i];
    
    for( Integer ij = StartIndex_[i]; ij < end; ij++ ) {
        if( Edges_[ij] == j )  return true;
    }
    
    return false;
}

/**********************************************************/

template <typename T>
inline Integer DependencyGraph<T>::GetNumNodes() const
{
    return NumNodes_;
}

/**********************************************************/

template <typename T>
inline Integer DependencyGraph<T>::GetNumEdges() const
{
    Integer ne = 0;
    
    for( Integer i = 1; i <= NumNodes_; i++ )  ne += NodeSize_[i];

    return ne;
}

/**********************************************************/

template <typename T>
inline Integer DependencyGraph<T>::
GetNumEdges( const Integer node ) const
{
#ifdef  DEBUG_DEPENDENCYGRAPH
    if( node < 1 || node > GetNumNodes() ) {
        Warning( __FILE__, __LINE__, "DependencyGraph::GetNumEdges(%d) : "
                 "node %d is out of range [1,%d] -> returning 0. In "
                 "non-debug mode this call might cause a segmentation "
                 "fault\n", node, node, GetNumNodes() );
        return 0;
    }
#endif

    return NodeSize_[node];
}

/**********************************************************/

template <typename T>
inline Integer DependencyGraph<T>::GetMaxNumEdges() const
{
    return NumEdges_;
}

/**********************************************************/

template <typename T>
inline Integer DependencyGraph<T>::
GetMaxNumEdges( const Integer node ) const
{
#ifdef  DEBUG_DEPENDENCYGRAPH
    if( node < 1 || node > GetNumNodes() ) {
        Warning( __FILE__, __LINE__, "DependencyGraph::GetMaxNumEdges(%d) : "
                 "node %d is out of range [1,%d] -> returning 0. In "
                 "non-debug mode this call might cause a segmentation "
                 "fault\n", node, node, GetNumNodes() );
        return 0;
    }
#endif

    return StartIndex_[node+1] - StartIndex_[node];
}

/**********************************************************/

template <typename T>
inline const Integer* DependencyGraph<T>::
GetEdges( const Integer node ) const
{
#ifdef  DEBUG_DEPENDENCYGRAPH
    if( node < 1 || node > GetNumNodes() ) {
        Warning( __FILE__, __LINE__, "DependencyGraph::GetEdges(%d) : "
                 "node %d is out of range [1,%d] -> expect a crash\n",
                 node, node, GetNumNodes() );
    }
#endif

    return (const Integer*)(Edges_ + StartIndex_[node]);
}

/**********************************************************/

template <typename T>
inline Integer DependencyGraph<T>::
GetEdges( const Integer        node,
                Integer *const edges ) const
{
    memcpy( edges, Edges_ + StartIndex_[node],
            NodeSize_[node] * sizeof(Integer) );
    return NodeSize_[node];
}

/**********************************************************/

template <typename T>
void DependencyGraph<T>::Reset()
{
    
    // destroy start array only if it is mine
    if( ownStartIndex_ ) {
			delete [] ( StartIndex_ );  StartIndex_  = NULL;
		}
    // destroy other arrays
    delete [] ( NodeSize_ );  NodeSize_  = NULL;
    delete [] ( Edges_ );  Edges_  = NULL;

    NumNodes_      = 0;
    NumEdges_      = 0;
    ownStartIndex_ = false;
}

/**********************************************************/

template <typename T>
std::ostream& DependencyGraph<T>::Print( std::ostream& out,
                                         const bool color ) const
{

    if( NumNodes_ == 0 ) {
        out << "graph is empty\n";
    } else {
        out
        << fg_green << "DependencyGraph"
        << fg_reset << std::endl
        << "  "<<NumNodes_ << "  nodes"<<std::endl
        << "  "<<GetNumEdges()<<':'<<NumEdges_<<" edges"
        << std::endl;
        for( Integer i = 1; i <= NumNodes_; i++ ) {
            out << fg_cyan <<'['<<i
                << ']'<< fg_reset <<": ";
            for( Integer j = 0; j < NodeSize_[i]; j++ ) {
                out << "[" << Edges_[StartIndex_[i]+j] << "] ";
            }
            out
            <<  fg_cyan <<'<'<< NodeSize_[i]
            << ':'<<(StartIndex_[i+1] - StartIndex_[i])
            << '>'<< fg_reset <<std::endl;
        }
    }

    return out;
}

/**********************************************************/

template <typename T>
bool DependencyGraph<T>::CreateArrays( const int  numnodes,
                                       const int  numedges,
                                       const bool createstart )
{

    // create start array, if demanded
    if( createstart ) {
        // if there is already a start array ...
        if( StartIndex_ ) {
            // if this class owns the start array ...
            if( ownStartIndex_ ) {
                // ... check if we can reuse it.
                if( numnodes != NumNodes_ ) {
                    delete [] ( StartIndex_ );  StartIndex_  = NULL;
                }
            // if the present start array is not owned by this
            // class remove by simply removing the access
            } else {
                StartIndex_ = NULL;
            }
        }
        // (StartIndex_ != NULL) means, that we have found an
        // array, we can reuse
        if( !StartIndex_ )  NEWARRAY( StartIndex_, Integer, numnodes+1 );
        ownStartIndex_ = true;
    }
    // delete old arrays, if their sizes do not match
    if( NodeSize_  && numnodes != NumNodes_ ) {
        delete [] ( NodeSize_ );  NodeSize_  = NULL;
		}
    if( Edges_ && numedges != NumEdges_ ) {
        delete [] ( Edges_ );  Edges_  = NULL;
    }
    // create new arrays
    if( NULL == NodeSize_ )  NEWARRAY( NodeSize_, Integer, numnodes );
    if( NULL == Edges_    )  NEWARRAY( Edges_,    Integer, numedges );
    // NOTE that there is NO INITIALISATION, for example of NodeSize_, here
//    for( Integer i = 1; i <= numnodes; i++ )  NodeSize_[i] = 0;
    
    NumNodes_ = numnodes;
    NumEdges_ = numedges;    

    return true;
}

/**********************************************************/

template <typename T>
void DependencyGraph<T>::QuickSort(       Integer *const array,
                                    const Integer        length )
{
    if( length <= 1 )  return;

    const Integer last     = length-1,
                  splitter = array[last];
    Integer temp, i = 0;

    // splitt array
    for( Integer j = 0; j < last; j++ ) {
        if( array[j] < splitter ) {
            temp       = array[j];
            array[j]   = array[i];
            array[i++] = temp;
        }
    }
    array[last] = array[i];
    array[i]    = splitter;

    // quicksort the two parts
    QuickSort( array,     i++        );
    QuickSort( array + i, length - i );
}

/**********************************************************/
} // namespace CoupledField

/**********************************************************
 * print-out operator
 **********************************************************/

namespace std {
template <typename T> std::ostream& operator <<
( std::ostream& out, const OLAS::DependencyGraph<T>& S ) {
    return S.Print( out ); }
}

/**********************************************************/
#ifdef DEBUG_TO_CERR
#undef debug
#endif // DEBUG_TO_CERR
/**********************************************************/
