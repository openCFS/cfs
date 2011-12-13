// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef OLAS_BASEGRAPH_HH
#define OLAS_BASEGRAPH_HH

#include <iostream>
#include <string>
#include <vector>

#include "General/defs.hh"
#include "General/environment.hh"
#include "baseordering.hh"


// Variation of the GRAPH-Implementation
//
// 1) Each edge is inserted, regardless if it is already contained in the
//    matrix graph. Fast, but not memory efficient
// 2) For each edge it is checked, if it is already contained in the graph.
//    Slow, but memory efficient
//
// Second variant is used, if GRAPH_VECTOR_SORT macro is enabled. This is
// done by setting FAST_EDGE_INSERTION = no in Makefile.option


namespace CoupledField {

  //! Base Class for handling the graph associated with a matrix
  //! This class represents the base class for all classes related to the
  //! handling of graphs associated with a matrix.

  // It basically implements a node graph based on linked lists which can be
  // transformed to compressed storage after Setup
template <class TYPE> class StdVector;

  class BaseGraph {

  public:

    //! Shortcut for an STL vector of unsigned integers (UInt)
    typedef std::vector<UInt> NodeList;

    //! Shortcut for an iterator over an STL vector of unsigned integers (UInt)
    typedef NodeList::iterator NodeListIterator;


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
    //! \param reorder  Specifies the re-ordering strategy to be applied to
    //!                 to the graph once it was completely assembled.
    //! \note In the current implementation we use nRows to set the value of
    //!       the numNodes_ attribute, i.e. we generate a vertex for each
    //!       row of the matrix.
    BaseGraph( UInt nRows, UInt nCols, BaseOrdering::ReorderingType reorder );

    //! Insert data into a graph

    //! This constructor can be used to create a compressed graph immediately,
    //! if the number of edges and the connectivities are already known.
    //! \param nRows    number of rows of associated matrix
    //! \param nCols    number of columns of associated matrix
    //! \param numEdge  number of edges (including self-edges, if present)
    //! \param cs_node  pointers to the position of the connectivity for node
    //!                 i in cs_edge
    //! \param cs_edge  stores the connectivity information for all nodes
    //! \param reorder  Specifies the re-ordering strategy to be applied to
    //!                 to the graph once it was completely assembled.
    BaseGraph( UInt nRows, UInt nCols, UInt numEdge, UInt *cs_node,
               UInt *cs_edge, BaseOrdering::ReorderingType reorder );

    //! Default destructor
    virtual ~BaseGraph();

    //! Method to be called after all vertices were inserted into the graph

    //! This method must be called after all vertices were inserted into the
    //! graph. This finalises the assembly phase of the graph and triggers
    //! the re-ordering and the conversion into CRS storage format.
    //! \param order One-based array for storing the re-ordering vector. If no
    //!              re-ordering is performed, this may be a NULL pointer
    void FinaliseAssembly( StdVector<UInt>* order );

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
    void AddVertexNeighbours( std::vector<UInt> vertexList,
                              std::vector<UInt> neighbourList );

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


#ifdef DEBUG_BASEGRAPH
      if ( amAssembled_ == false ) {
        EXCEPTION("Attempt to obtain information from graph object, "
                 << "before assembly was completed by calling "
                 << "FinaliseAssembly()");
      }
#endif

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

    //! Query the state of the graph

    //! This method can be used to query the status of the graph object. If
    //! the assembly of the graph was finalised, then this method will return
    //! true and it is safe to use other query methods on the graph.
    bool IsAssembled() {
      return amAssembled_;
    }

    //! Query wether the graph has been reordered
    bool IsReordered() {
      return amReordered_;
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
    void GetBandwidth( UInt &bwlower, UInt& bwupper );

    //! Query the number of neighbours of node i

    //! This method returns the number of entries in the i-th row of the
    //! companion matrix of the graph. This number corresponds to the number
    //! of edges leaving vertex i.
    //! \note The return value includes the self-reference/loop for vertex i,
    //!       if this is contained in the graph.
    inline UInt GetRowSize( UInt i ) const {


#ifdef DEBUG_BASEGRAPH
      if ( amAssembled_ == false ) {
        EXCEPTION("Attempt to obtain information from graph object, "
                 << "before assembly was completed by calling "
                 << "FinaliseAssembly()");
      }
#endif

      return csNodes_[i+1] - csNodes_[i];
    }

    //@}


    // =======================================================================
    // DEBUGGING methods
    // =======================================================================

    //@{
    //! \name Methods for debugging

    std::string ToString() const
    {
      std::ostringstream os;
      Print(os);
      return os.str();
    }
    
    //! prints the complete Graph for Debugging
    void Print(std::ostream &os) const;
  
    //! print the Graph to the .las file
    void Print() const {
      Print(*cla);
    }

    //@}

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

    //! Uncompressed matrix graph, consisting of STL lists
    NodeList *element_;

    //! Number of nodes/vertices in the matrix graph
    UInt numNodes_;

    //! Number of edges in the graph
    UInt nne_;

    //! store the lower bandwidth of the graph
    UInt bwlower_;

    //! store the upper bandwidth of the graph
    UInt bwupper_;


    // =======================================================================
    // REORDERING stuff
    // =======================================================================

    //! Strategy for re-ordering the graph
    BaseOrdering::ReorderingType newOrder_;

    //! Keep track of wether the graph has been reordered
    bool amReordered_;

    //! Compute a fill-reducing ordering of the graph. 

    //! The result is stored in order_, the inverse mapping 
    //! in iorder_. Currently, the ordering is obtained using Metis.
    //! If METIS is not defined, the ordering will be identity. 
    //! If the graph is not compressed, it will be compressed now.
    void Reorder( BaseOrdering::ReorderingType newOrder, StdVector<UInt>& order );


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

    //! Convert the linked list data structure to CRS format

    //! This method can be used to generate a representation of the graph as a
    //! CRS data structure. It is different from ConvertToCRS in the respect
    //! that the self-reference of a node (i.e. the edge connecting the node
    //! to itself / the diagonal matrix entry) is not stored. This is the
    //! format expected e.g. by Metis. The method will handle dynamic
    //! allocation of memory for the CRS structure.
    //! \param rptr on return stores adress of row pointer array of CRS
    //!             structure
    //! \param cidx on return stores adress of column index array of CRS
    //!             structure
    void ConvertToMetisCRS( UInt **rptr, UInt **cidx );

    //! This array contains the indices of connected nodes (column array)

    //! For every node we have a list of node indices of the nodes to which
    //! it is connected (including itself, if a 'diagonal' node). These lists
    //! are stored one after the other in this data array. This corresponds
    //! to the column index vector in compressed row storage format.
    UInt *csEdges_;

    //! This array points to the start of the rows in cs_edges

    //! This array contains for every node the starting index of its neighbour
    //! list in the cs_edgdes_ array. This corresponds to the row pointer
    //! vector in compressed row storage format.
    UInt *csNodes_;

    //! Keep track of whether assembly of the graph was finalised

    //! This attribute keeps track of the status of the graph. It is initialy
    //! set to false and switched to true, once the assembly of the graph was
    //! finalised and it was converted to the CRS structure.
    bool amAssembled_;

    //! Default Constructor

    //! The default constructor is dis-allowed since the graph needs
    //! a parameter object, which is passed to the constructor
    BaseGraph() {};

    //! Number of columns in the matrix associated with the graph

    //! This graph class is used to handle a matrix graph, i.e. the graph
    //! represents the sparsity pattern of a matrix. This attribute stores
    //! the number of columns in this associated matrix.
    UInt numColsMat_;

    //! Number of rows in the matrix associated with the graph

    //! This graph class is used to handle a matrix graph, i.e. the graph
    //! represents the sparsity pattern of a matrix. This attribute stores
    //! the number of rows in this associated matrix.
    UInt numRowsMat_;

  };

} // namespace

#endif // OLAS_BASEGRAPH_HH
