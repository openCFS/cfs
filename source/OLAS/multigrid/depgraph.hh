/* $Id$ */

#ifndef OLAS_DEPGRAPH_HH
#define OLAS_DEPGRAPH_HH

#include "matvec/crs_matrix.hh"
#include "multigrid/ppflags.hh"

namespace OLAS {
/**********************************************************/

//! dependency graph for the AMG topology

/*! A special graph for the dependencies needed in the AMG setup phase.
 *  This graph does not inherit any of the other graph classes, since
 *  it needs a totally different flexibility and therefore storage
 *  structure, that could be added to the existing graph classes only
 *  with reduced performance. As this graph is heavily used and even
 *  might be modified (i.e. augmented and reduced) in the AMG setup
 *  phase, a reusage of one of the other graph classes is likely to
 *  spoil the setup's performance.
 *  <br>
 *  The graph's data structure is similar to the one of a CRS matrix,
 *  with the exception that the array containing the start indices does
 *  not tell us that there are StartIndex_[i+1] - StartIndex_[i] edges
 *  at node [i], but maximal StartIndex_[i+1] - StartIndex_[i] edges
 *  at node [i]. The real number of edges at node [i] is stored in a
 *  separat array NodeSize_. This allows the efficient implementation
 *  of a dynamic graph, of which the maximal numbers of edges at each
 *  node are known.
 *  <br>
 *  Some notes:
 *   - the graph is designed to keep edges (i,j) with i,j <= 0,
 *     otherwise for example some methods' return values would lose
 *     their meaning (e.g. GetMaxEdge()).
 */
template <typename T>
class DependencyGraph
{
    public:

        DependencyGraph();
        ~DependencyGraph();

        //! creates an empty graph

        /*! Creates an empty graph with a predicted maximal size, given
         *  as number of nodes and maximal number of edges per node.
         *  Here, of course, <b>all</b> arrays are created internally.
         *  \param num_nodes number of nodes
         *  \param num_edges_per_node expected maximal number of edges
         *         per node
         *  \return true, if successful, false, otherwise
         */
        bool Create( const Integer num_nodes,
                     const Integer num_edges_per_node );

        //! creates graph with already set "diagonal" edges at first position

        /*! Creates graph with a predicted maximal size, given as number
         *  of nodes and maximal number of edges per node. All edges (i,i).
         *  for i in [1,num_nodes] are inserted in this creation routine.
         *  Here, of course, <b>all</b> arrays are created internally.
         *  \param num_nodes number of nodes
         *  \param num_edges_per_node expected maximal number of edges
         *         per node
         *  \return true, if successful, false, otherwise
         */
        bool CreateWithDiagonals( const Integer num_nodes,
                                  const Integer num_edges_per_node );

        //! creates an empty graph

        /*! Creates an empty graph with a predicted maximal size, given
         *  as the maximal number of edges for each node. The total
         *  maximal number of edges can be calculated from these data.
         *  This is done, whenever the parameter <code>num_edges</code> is
         *  not passed (i.e. a value < 0). But often this total number
         *  is already known in the calling scope of this function. So
         *  it is only a tiny feature for slight acceleration, to pass
         *  the correct number of maximal total number of edges.
         *  \param num_nodes number of nodes in the graph
         *  \param node_size array containing the maximal number of edges
         *         for each node
         *  \param num_edges maximal total number of edges in the graph
         *  \return true, if successful, false, otherwise
         */
        bool Create( const Integer        num_nodes,
                     const Integer *const node_size,
                           Integer        num_edges = -1 );

        //! creates an empty graph

        /*! Creates an empty graph from a specified number of nodes and
         *  a passed array of start indices. This array of start indices
         *  is built like the one in a compressed row matrix (CRS_Matrix)
         *  and has length [num_nodes+1] (one-based indexed). The entry
         *  startarray[num_nodes+1] contains the maximal number of edges
         *  +1. Typically such a start array is received from a CRS_Matrix
         *  or another graph object. Already present buffers are reused,
         *  if possible, but the content (which means in particular the
         *  the array NodeSize_ is reset to all node sizes 0).
         *  \param num_nodes number of nodes
         *  \param startarray array of size [num_nodes+1], which contains
         *         start indices for all nodes. In particular position
         *         [num_nodes+1] contains the start index for a virtual
         *         node [num_nodes+1].
         *  \param accroach_startarray if true, the startarray will be
         *         administrated from now on from this class, in
         *         particular deleted, if necessary
         *  \param copy_startarray if true, the startarray will be copied
         *         in a buffer that is created by this class. Therefore
         *         accraoch_startarray == true and copy_startarray ==
         *         true does not make sense and will cause a warning in
         *         debug mode
         *  \return true, if the creation was successful
         */
        bool Create( const Integer        num_nodes,
                     const Integer *const startarray,
                     const bool           accroach_startarray,
                     const bool           copy_startarray );
        
        //! creates the whole dependency graph from a crs-matrix

        /*! Creates a graph from a passed matrix in compressed row
         *  storage format. The graph either can be built as an exact
         *  copy of the matrix graph or only take the matrix graph as
         *  the maximal graph, itself can grow up to. In the last case
         *  this means, that a graph will be prepared, that <b>could</b>
         *  take the whole matrix graph, but none of the edges is copied.
         *  This is implemented for example for the graph of strong
         *  dependencies in a matrix. We cannot predict in detail, how
         *  the graph will look like, but we know, that it will be a
         *  subset of the matrix graph. To create such an empty graph
         *  set buildemptygraph to true. If the graph should be independent
         *  of the further existence of the matrix, the array containing
         *  the start indices should be copied (=> copystartarray = true)
         *  \param matrix matrix the graph should be built form
         *  \param buildemptygraph if true, an empty graph is created,
         *         which is limited exactly by the matrix graph
         *  \param copystartarray if true, the array containing the start
         *         indices will be copied
         */
        bool Create( const CRS_Matrix<T>& matrix,
                     const bool buildemptygraph = false,
                     const bool copystartarray  = false );

        //! insert an array with start positions
        
        /*! Inserts an array containing start positions of the nodes.
         *  Per default DependencyGraph only <b>reads</b> and uses the
         *  passed array, but does <b>not delete</b> it, e.g. in Reset()
         *  or the destructor. If accroach is passed as true (default
         *  is false), the class will treat the passed array as its own and
         *  also delete it, if necessary (therefore the array must have
         *  been created with the OLAS makro CreateArray). The array can
         *  also be copied in one, that is owned by this object. Therefor
         *  copy must be true (default false). The new array is only
         *  accepted if there is none already present, but by setting
         *  force to true (default false) the insertion can be forced,
         *  and an eventually old one will be replaced. In principle this
         *  is a function, that should only be used as a public function,
         *  if you really know, what you are doing.
         *  \param startarray the array containing the start indices
         *         <b>NOTE</b>: startarray is expected to be one based
         *         indexed, e.g. as returned by CRS_Matrix::GetRowPointer(),
         *         and must also contain a start index for a virtual node
         *         with index [size+1]!
         *  \param size size-1 of startarray (== number of nodes)
         *  \param accroach if set to true (default false), the control
         *         of the passed startarray is passed to the graph object
         *  \param copy if true, a new class-owned array is created and the
         *         passed array is copied (accroach = false will not have
         *         any effect then)
         *  \param force  if set to true (default false), an eventually
         *         already present array will be replaced
         */
        bool InsertStartArray( const Integer *const startarray,
                               const Integer size,
                               const bool    accroach = false,
                               const bool    copy     = false,
                               const bool    force    = false  );

        //! resets the object to the state after creation with standard constructor
        void Reset();

        //! adds the edge (i,j) to the graph

        /*! Adds the edge (i,j) to the graph. This insertion is done without
         *  any savety checks. Use this function <b>only</b>, if you know for
         *  sure, that the edge (i,j) is <b>not already element</b> of the
         *  graph. If the makro DEBUG_DEPENDENCYGRAPH is defined, the following
         *  checks are applied, and an error message is printed, if
         *   - node number i is not in the range of valid nodes
         *   - the edge (i,j) already exists (see also AddEdgeSavely())
         *   - the node i cannot add further nodes, because it has reached
         *     its maximal number of edges
         */
        inline void AddEdge( const Integer i, const Integer j );
        
        //! adds the edge (i,j) to the graph, but ommitting double edges

        /*! Adds the edge (i,j) to the graph. This function is similar to
         *  AddEdge(), except that the edge (i,j) is only added, if it is
         *  not already present. Use this function instead of AddEdge(), if
         *  you cannot be sure, if (i,j) might already exist in this graph.
         *  If the makro DEBUG_DEPENDENCYGRAPH is defined, the following
         *  checks are applied, and an error message is printed, if
         *   - node number i is not in the range of valid nodes
         *   - the node i cannot add further nodes, because it has reached
         *     its maximal number of edge.
         */
        inline void AddEdgeSavely( const Integer i, const Integer j );

        //! adds edge (i,j) to the graph, that is assumed to be sorted

        /*! Adds the edge (i,j) to the graph. This function assumes the
         *  graph to be sorted, and will insert the edge sorted, too.
         *  Therefore do not mix this insertion routine with the other
         *  ones. You may optionally specifiy the lower most position in
         *  the edge array, the new edge should be inserted at. This
         *  might be usefull for example, if you want to exclude the
         *  first position from sorting, since it should hold the edge
         *  (i,i) for all nodes i (see also AddEdgeSortedAfterDiag).
         *  Also note, that this function avoids redundant insertion
         *  of an already present edge.
         */
        inline void AddEdgeSorted( const Integer i, const Integer j,
                                   const Integer first_position = 0 );

        //! adds the edge (i,j) to a sorted graph, but after the diagonal position

        /*! This function works similar to AddEdgeSorted, except that
         *  it will insert the edge after the leading position, that
         *  might be preserved for a special edge, e.g. the diagonal edge
         *  (i,i).
         */
        inline void AddEdgeSortedAfterDiag( const Integer i, const Integer j );
        
        //! adds the edge (i,j) to the graph as edge [position] at node i

        /*! This function sets the edge (i,j) to the graph at position
         *  [position] in the edge array (range [0,GetMaxNumEdges(i)-1]).
         *  This function modifies the graph data directly without any
         *  check in non-debug mode. The edge is not "added" to the other
         *  edges, it is "set". If there is another edge present at position
         *  [position], it is overwrittten. Use carefully. In debug mode
         *  (i.e DEBUG_DEPENDENCYGRAPH is defined) the following checks
         *  are applied, and an error message is printed, if
         *   - i is out of range (edge not set)
         *   - position is out of range (edge not set)
         *   - the edge (i,j) is already present <b>at a different</b>
         *     position at node i. (edge is set nevertheless)
         */
        inline void SetEdgeAtPosition( const Integer i, const Integer j,
                                       const Integer position );
        
        //! removes edge (i,j) from the graph

        /*! Removes edge (i,j) from the graph. The absence of (i,j) in
         *  the graph is no problem. <b>Note</b> that the ordering of
         *  the edges will be changed by this function (e.g. if the edge
         *  (i,4) is removed from node i, that has the edges (1,4,7,9),
         *  the ordering of the edges will remain as (1,9,7)). Only if
         *  the makro flag DEBUG_DEPENDENCYGRAPH is defined, there will
         *  be a boundary check of the node numbers i and j.
         */
        inline void RemoveEdge( const Integer i, const Integer j );

        //! removes all edges of node [i] in the graph
        inline void RemoveAllEdges( const Integer i );

        //! shift all edges: (i,j) \f$\longrightarrow\f$ (i,j + offset)

        /*! \param offset (i,j) \f$\longrightarrow\f$ (i,j + offset)
         *  \result \f$ \max_{i} \{ j | (i,j) \in G \} \f$ after the
         *          shift, 0 if the graph is empty
         */
        inline Integer ShiftEdges( const Integer offset );

        //! sort the edges at all nodes
        inline void Sort();
        //! like Sort(), but positions "diagonal edges" at first position
        inline void SortDiagFirst();

        //! assigns the transposed of the passed graph

        /*! This function assigns the transposed graph of the passed
         *  one to this graph. Transposed graph means: \f$ (i,j) \in
         *  G \Leftrightarrow (j,i) \in G' \f$. If G has got n nodes
         *  and contains the edge \f$ (i,J), J = \max_{1 \le i \le G}
         *  \{ j | (i,j) \in G \} \f$, then G' has got \f$J\f$ nodes.
         *  Therefore the function first must find this \f$J\f$. If
         *  the caller of this method already knows that \f$J \le n
         *  \f$, this search can be suppressed by specifying
         *  treat_as_square as true. The transposed graph is created
         *  in a compact format, which means that no further nodes can
         *  be added in this graph. But by specifying overlap_per_node
         *  space for exactly overlap_per_node edges is reservated.
         *  Also note, that the transposed graph is automatically
         *  sorted due to the loop ordering in the algorithm, but the
         *  diagonal element is not placed at leading position.
         *  \param graph the source graph, from that the transposed
         *         graph is taken and assigned
         *  \param use_start_array if true, the transposed graph will
         *         share the start index array with the original graph.
         *         For example in AMG we know, that the original graph
         *         S of strong dependencies and its transpose ST both
         *         are subgraphs of a matrix graph. Therefore we can
         *         create S from the matrix without copying the start
         *         index array and then share the start index array also
         *         in the transpose ST.
         *  \param treat_as_square if true, the maximal j for edges
         *         (i,j) is not searched and assumed to be equal to
         *         or less than the number of nodes. In debug mode
         *         the validity of a "true" is verified (if true
         *         although the maximal j is greater than n, in
         *         non-debug mode the program may crash)
         *  \param overlap_per_node number of empty space per node,
         *         that should be preserved.
         */
        bool AssignTransposed( const DependencyGraph<T>& graph,
                               bool use_start_array = false,
                               bool treat_as_square = false,
                               Integer overlap_per_node = 0 );

        //! checks, if this graph is the transposed of the passed one

        /*! Checks, if this graph is the transposed of the passed graph.
         *  This function is quite slow and should be used for debugging
         *  purposes only.
         *  \param graph the graph, that should be compared with this one
         *  \return true, if the passed graph and this one are transposed
         *          to each other, otherwise false.
         */
        bool CheckTransposition( const DependencyGraph<T>& graph ) const;

        //! returns true, if the edge (i,j) is element of the graph
        inline bool IsElement( const Integer i, const Integer j ) const;
        //! returns the number of nodes
        inline Integer GetNumNodes() const;
        //! returns the total number of edges in the graph (calculated at calling time)
        inline Integer GetNumEdges() const;
        //! returns the number of edges at node "node"
        inline Integer GetNumEdges( const Integer node ) const;
        //! returns the maximal number of edges in the whole graph
        inline Integer GetMaxNumEdges() const;
        //! returns the maximal number of edges at node "node"
        inline Integer GetMaxNumEdges( const Integer node ) const;

        /*! \brief returns a pointer to the edges of node "node" (valid
         *  range: [0, GetNumEdges(node)-1]).
         */
        inline const Integer* GetEdges( const Integer node ) const;
        //! writes the edges of node "node" into the buffer "edges" and
        //! returns the number of written edges
        inline Integer GetEdges( const Integer        node,
                                       Integer *const edges ) const;
        //! returns the array with start indices
        const Integer* GetStartArray() const { return StartIndex_; }

        //! returns \f$ \max_{i} \{ j | (i,j) \in G \} \f$ (0, if graph is empty)
        inline Integer GetMaxEdge() const {
            Integer maxEdge = 0;
            for( Integer i = 1; i <= NumNodes_; i++ ) {
                const Integer end = StartIndex_[i] + NodeSize_[i];
                for( Integer ij = StartIndex_[i]; ij < end; ij++ ) {
                    if( maxEdge < Edges_[ij] )  maxEdge = Edges_[ij];
                }
            }
            return maxEdge;
        }

        //! prints out the graph, identical to operator <<
        std::ostream& Print( std::ostream& out,
                             const bool color = false ) const;

    protected:

        //! creates the needed arrays

        /*! Creates the arrays, that are needed in this class, i.e.
         *  StartIndex_[numnodes+1], NumEdges_[numnodes], and Edges_[numedges].
         *  As the array containing the nodes' start indices may be set
         *  separately, its creation can be suppressed by setting the
         *  parameter createstart to false. This method checks, if
         *  already present arrays can be reused, so this needs not to
         *  be checked in the calling function. The arrays are <b>not</b>
         *  initialized.
         *  \param numnodes number of nodes in the graph
         *  \param numedges total number of edges in the graph
         *  \param createstart if true, the array with start indices is
         *         also created
         */
        bool CreateArrays( const int  numnodes,
                           const int  numedges,
                           const bool createstart = true );

        //! a quick sort implementations for Integer arrays

        /*! An implementation of Quick Sort for Integers, needed for
         *  sorting the edges at a node.
         *  \param array pointer to the Integer array, that should be
         *         sorted. array is assumed to be zero based indexed,
         *         i.e. array[0] is the first valid entry, array[length-1]
         *         its last.
         *  \param length number of intergers in the array.
         */
        static void QuickSort(       Integer *const array,
                               const Integer        length );
        

        Integer  NumNodes_, //!< number of nodes in the graph
                 NumEdges_; //!< maximal number of <b>possible</b> edges

        Integer *StartIndex_; //!< start indices of nodes in Edges_
        Integer *NodeSize_;   //!< number of connections for each node
        Integer *Edges_;      //!< array containing the edges

        bool ownStartIndex_;
};

/**********************************************************/
} // namespace OLAS

#endif // OLAS_DEPGRAPH_HH
