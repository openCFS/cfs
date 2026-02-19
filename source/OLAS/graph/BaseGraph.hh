#ifndef OLAS_BASEGRAPH_HH
#define OLAS_BASEGRAPH_HH

#include <iostream>
#include <vector>
#include <boost/unordered/unordered_flat_set.hpp>

#include "General/Environment.hh"
#include "BaseOrdering.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField 
{
  //! Base Class for handling the graph associated with a matrix
  
  //! This class represents the base class for all classes related to the
  //! handling of graphs associated with a matrix.
  //! It basically implements a node graph based on linked lists which can be
  //! transformed to compressed storage after Setup
  class BaseGraph {

  public:

    // =======================================================================
    // CONSTRUCTORS
    // =======================================================================
  
    //! \name Methods for construction, assembly and destruction
    //@{

    //! Constructor

    //! This constructor initialises some attributes to default values,
    //! stores the given re-ordering type and allocates memory to hold the
    //! neighbour lists for each node.
    //! \param nRows    number of rows of associated matrix
    //! \param nCols    number of columns of associated matrix
    //! \note In the current implementation we use nRows to set the value of
    //!       the numNodes_ attribute, i.e. we generate a vertex for each
    //!       row of the matrix.
    BaseGraph( UInt nRows, UInt nCols, unsigned int estimated_row_size);

    //! Insert data into a graph

    //! This constructor can be used to create a compressed graph immediately,
    //! if the number of edges and the connectivities are already known.
    //! \param nRows    number of rows of associated matrix
    //! \param nCols    number of columns of associated matrix
    //! \param numEdge  number of edges (including self-edges, if present)
    //! \param cs_node  pointers to the position of the connectivity for node
    //!                 i in cs_edge
    //! \param cs_edge  stores the connectivity information for all nodes
    BaseGraph( UInt nRows, UInt nCols, UInt numEdge, UInt* cs_node, UInt* cs_edge );

    //! Default destructor
    virtual ~BaseGraph();

    //! Method to be called after all vertices were inserted into the graph

    //! This method must be called after all vertices were inserted into the
    //! graph. This finalizes the assembly phase of the graph and triggers
    //! the re-ordering and the conversion into CRS storage format.
    //! This method two possible strategies:
    //! - 1) If useExternalOrdering = false, the internal reordering strategy
    //!      will be used (if set) and the new ordering array will be stored
    //!      in vertexOrder
    //! - 2) If useExternalOrdering = true, the two vector vertexOrder and 
    //!      edgeOrder will be used to re-order the graph. If only the vertex-
    //!      reordering array is provided, it will be also applied to the 
    //!      edges.
    //! \param reorder  Specifies the re-ordering strategy to be applied to
    //!                 to the graph once it was completely assembled.
    //! \param useExternalOrdering use the ordering provided in the array(s)
    //!                            vertexOrder (optionally: also edgeOrder)
    //! \param vertexOrder new ordering for vertices. If useExternalOrdering
    //!                    = false, the internal reordering array will
    //!                    be taken. If useExternalOrdering = true, this vector
    //!                    will be used as input.
    //! \param edgeOrder new ordering for edges as input. Can be omitted. 
    //! \param order One-based array for storing the re-ordering vector. If no
    //!              re-ordering is performed, this may be a NULL pointer
    void FinaliseAssembly( BaseOrdering::ReorderingType reorder, 
                           bool useExternalOrdering,
                           StdVector<UInt>* vertexOrder,
                           StdVector<UInt>* edgeOrder = nullptr );

    //! Add edges between vertices and their neighbours

    //! This method can be used to add neighbours to the neighbour list of
    //! vertices. For each vertex contained in the vertexList vector we add
    //! an edge connecting the vertex to each vertex contained in the
    //! neighbourList.
    //! \note Both input parameters are zero-based STL vectors!
    //! \param vertexList     a vector containing vertex numbers of vertices
    //!                       to which neighbours/edges are to be added
    //! \param neighbourList  a vector containing vertex numbers of the
    //!                       neighbours of the vertices in vertexList
    void AddVertexNeighbours( std::vector<UInt>& vertexList,
                              std::vector<UInt>& neighbourList );

    
    //! Set block definition
    
    //! This method can be used to define non-overlapping blocks within a
    //! graph / matrix. This can be later used to define block matrices.
    //! \param indexBlock contains the blocks. First index is the blockindex
    void SetBlockInfo( const StdVector<StdVector<UInt> >* indexBlocks );
    
    
    //! Set related diagonal graphs of SBM row / col
    
    //! If this graph object represents not a diagonal block in the 
    //! SBM-structure, we can set here the related diagonal graphs of
    //! the SBM structure. This is especially useful, if the graph has a 
    //! block definition. In this case, the non-diagonal graph is not 
    //! allowed to perform any reordering to match the block information,
    //! but rather uses the block information provided by the diagonal
    //! graphs.
    //! \param rowDiaGraph graph object of the corresponding row diagonal block
    //! \param colDiagGraph graph object of the corresponding col diagonal block
    void SetRowColDiagGraphs( BaseGraph* rowDiagGraph,
                              BaseGraph* colDiagGraph );
    
    //@}

    // =======================================================================
    // QUERY METHODS
    // =======================================================================

    //@{ \name Query methods

    //! get the node neighbors for node i

    //! This method returns a pointer to an integer array containing the
    //! indices of the nodes that are connected to node i including possibly
    //! the self-reference for node i, if it exists. The list of indices is
    //! lexicographically sorted with increasing index number.
    inline UInt* GetGraphRow( UInt i ) {
      return csEdges_ + csNodes_[i];
    }

    //! Get total number of vertices in the graph

    //! This method can be used to query the number of vertices in the graph.
    //! \note In the case that the graph corresponds to the sparsity pattern
    //! of a rectangular or unsymmetric matrix the graph is a directed one.
    //! In this case we understand by the number of vertices in the graph
    //! the number of vertices that are (potential) start points of an edge.
    UInt GetSize() const {
      return numNodes_;
    }

    //! Get number of edges

    //! This method can be used to query the number of edges in the graph,
    //! which corresponds to the number of non-zero entries in its CRS
    //! representation
    UInt GetNNE() const {
      return nne_;
    }

    //! Get number of columns in matrix associated with the graph

    //! When the graph corresponds to the sparsity pattern of a symmetric
    //! matrix we can use the number of vertices in order to determine the
    //! dimension of this matrix.
    //! In the case of an unsymmetric graph the number of vertices corresponds
    //! only to the number of rows of the associated matrix.
    //! This method allows to obtain the number of columns in this case.
    UInt GetNumColsMat() {
      return numColsMat_;
    }
    
    //! Query reordering type
    
    //! This method returns the requested reordering strategy.
    BaseOrdering::ReorderingType GetReorderType() const {
      return newOrder_;
    }

    //! Query the state of the graph

    //! This method can be used to query the status of the graph object. If
    //! the assembly of the graph was finalised, then this method will return
    //! true and it is safe to use other query methods on the graph.
    bool IsAssembled() {
      return amAssembled_;
    }

    //! Query wether the graph has been reordered
    bool IsReordered() {
      return isReordered_;
    }

    //! Return lower and upper bandwidth of the graph-connectivity matrix

    //! This method determines the lower and upper bandwidth of the matrix
    //! pattern corresponding to the graph.
    //! The lower bandwidth \f$w_l\f$ and the upper bandwidth \f$w_u\f$ of a
    //! matrix are defined as follows
    //! \f[
    //! \begin{array}{l@{\,=\,}l}
    //! w_l & \max\limits_{i=1,\ldots,n}\left( \max \left\{ |i-j|\,\Bigm|\,j<i
    //!       \mbox{ and } a_{ij}\neq 0 \right\}\right)\\[2ex]
    //! w_u & \max\limits_{i=1,\ldots,n}\left( \max \left\{ |i-j|\,\Bigm|\,j>i
    //!       \mbox{ and } a_{ij}\neq 0 \right\}\right)\enspace.
    //! \end{array}
    //! \f]
    //! In addition the average total bandwidth is returned
    void GetBandwidth( UInt &bwlower, UInt& bwupper, UInt& bwAvg );

    //! Query the number of neighbours of node i

    //! This method returns the number of entries in the i-th row of the
    //! companion matrix of the graph. This number corresponds to the number
    //! of edges leaving vertex i.
    //! \note The return value includes the self-reference/loop for vertex i,
    //!       if this is contained in the graph.
    inline UInt GetRowSize( UInt i ) const {
      return csNodes_[i+1] - csNodes_[i];
    }
    
    //! Return block definition of rows
    
    //! This method returns the row block definition. If this graph represents
    //! a diagonal one, this method just return the internal 
    //! #sortedBlocks-array, otherwise it returns the row-block-definition of 
    //! the related block row-diagonal.
    StdVector<std::pair<UInt,UInt> >& GetRowBlockDefinition();
    
    //! Return block definition of columns
    
    //! This method returns the col block definition. If this graph represents
    //! a diagonal one, this method just return the internal 
    //! #sortedBlocks-array, otherwise it returns the col-block-definition of 
    //! the related block col-diagonal.
    StdVector<std::pair<UInt,UInt> >& GetColBlockDefinition();
    //@}

    /** this is the full graph timer (a sub-timer) 
     * Started by BaseGraph constructor and ended by AlgebraicSystem::StopGraphTimers()  */
    boost::shared_ptr<Timer> timer;

  protected:

    //! Counts the number of edges

    //! While the graph is uncompressed this function can be used
    //! to count how many edges are in the node lists
    //! the result is stored in the nne_ member
    void CountNNE();

    //! Auxilliary methods used by FinaliseAssembly

    //! This is an auxilliary methods used by FinaliseAssembly(). It uses
    //! STL algorithms for sorting the edge lists of each vertex and making
    //! the indices in each list unique.
    void SortLists();

    /** returns a debug string of the system with the comple graph for debugging
     * @see PrintCRS() */
    std::string ToString() const;
  
    /** print a crs structure to a debug string
     * @param structure row by row, otherwise one line with rows and one with column
     * @param level 0 = each rows by line, 1 = one line for row and col, 2 = metis 5 graph file formant for ndmetis */
    std::string PrintCRS(unsigned int numNodes, unsigned int* rows, unsigned int* cols, int level = 0) const;

    /** Compute a fill-reducing ordering of the graph.

    * The result is stored in order_, the inverse mapping in iorder_.
    * If the graph is not compressed, it will be compressed now.
    * @newOrder METIS, SLOAN  */
    void Reorder(BaseOrdering::ReorderingType newOrder, StdVector<UInt>& order);


    //! Create reordering to get consecutive arranges blocks
    
    //! This method tries to generate a reordering such that all defined
    //! subBlocks contain indices in sequential order. 
    //! If the method can not find a ordering (e.g. due to overlapping subBlocks)
    //! it will throw an exception.
    //! \param order 1-based reordering array
    void ReorderForBlocks( StdVector<UInt>& order );
    
    //! function for converting set contents to vector
    void MapSetToVector();

    // =======================================================================
    // CRS FORMAT stuff
    // =======================================================================

    //! Convert the linked list data structure to CRS format

    //! This method is used to transform the storage format of the graph from
    //! the vector of linked lists to a sparse matrix in compressed row
    //! storage (CRS) format. Here a node corresponds to an unknown and each
    //! non-zero entry \f$a_{ij}\f$ in the matrix corresponds to an edge from
    //! node \f$n_i\f$ to node \f$n_j\f$. The entries in each row are sorted
    //! lexicographically with increasing column index.
    void ConvertToCRS();
 
    /** common part of constructors */
    void Init();

    /** estimated row size. If too small we need re-hashing of setElements_ */
    unsigned int estimated_row_size_ = 0;

    PtrParamNode info_;

    //! Uncompressed matrix graph, consisting of STL lists
    StdVector<std::vector<unsigned int>> element_;

    //! Number of nodes/vertices in the matrix graph
    UInt numNodes_ = 0;

    //! Number of nodes without self-reference, i.e. without diagonal element
    UInt numNonDiagEntries_ = 0;
    
    //! Number of edges in the graph
    UInt nne_ = 0;

    //! store the lower bandwidth of the graph
    UInt bwlower_ = 0;

    //! store the upper bandwidth of the graph
    UInt bwupper_ = 0;
    
    //! store the average total bandwidth of the graph
    UInt bwavg_ = 0;


    // =======================================================================
    // Neighbouring graphs
    // =======================================================================
    
    //! Row diagonal graph
    BaseGraph* rowDiagGraph_ = nullptr;
    
    //! Column diagonal graph
    BaseGraph* colDiagGraph_ = nullptr;

    // =======================================================================
    // REORDERING stuff
    // =======================================================================

    //! Strategy for re-ordering the graph
    BaseOrdering::ReorderingType newOrder_ = BaseOrdering::NOREORDERING;

    //! Keep track of wether the graph has been reordered
    bool isReordered_ = false;


    //! This array contains the indices of connected nodes (column array)

    //! For every node we have a list of node indices of the nodes to which
    //! it is connected (including itself, if a 'diagonal' node). These lists
    //! are stored one after the other in this data array. This corresponds
    //! to the column index vector in compressed row storage format.
    UInt* csEdges_ = nullptr;

    //! This array points to the start of the rows in cs_edges

    //! This array contains for every node the starting index of its neighbor
    //! list in the cs_edgdes_ array. This corresponds to the row pointer
    //! vector in compressed row storage format.
    UInt* csNodes_ = nullptr;

    //! Keep track of whether assembly of the graph was finalized

    //! This attribute keeps track of the status of the graph. It is initially
    //! set to false and switched to true, once the assembly of the graph was
    //! finalized and it was converted to the CRS structure.
    bool amAssembled_ = false;

    //! Default Constructor

    //! The default constructor is dis-allowed since the graph needs
    //! a parameter object, which is passed to the constructor
    BaseGraph() {};

    //! Number of columns in the matrix associated with the graph

    //! This graph class is used to handle a matrix graph, i.e. the graph
    //! represents the sparsity pattern of a matrix. This attribute stores
    //! the number of columns in this associated matrix.
    UInt numColsMat_ = 0;

    //! Number of rows in the matrix associated with the graph

    //! This graph class is used to handle a matrix graph, i.e. the graph
    //! represents the sparsity pattern of a matrix. This attribute stores
    //! the number of rows in this associated matrix.
    UInt numRowsMat_ = 0;
    
    //! Original definition of blocks (unsorted)
    //! \note The graph object will not take ownership of this pointer, i.e.
    //! the user has to delete this pointer from outside after the graph
    //! object is deleted.
    const StdVector<StdVector<UInt> >* unsortedBlocks_ = nullptr;
    
    //! Final representation of blocks
    StdVector<std::pair<UInt,UInt> > sortedBlocks_;

    /** set for faster add of element neighbors - will be finalized to elements_
     * the rows are sorted in SortLists(), so unordered makes no numerical issues */
    StdVector<boost::unordered_flat_set<UInt>> setElements_;

    //! flag to check if element_ pointer is ready
    bool setToElemDone_ = false;

  };

} // namespace

#endif // OLAS_BASEGRAPH_HH
