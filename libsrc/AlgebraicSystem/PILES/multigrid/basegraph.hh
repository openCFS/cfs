#ifndef FILE_BASEGRAPH_CLA
#define FILE_BASEGRAPH_CLA

namespace CoupledField
{

//! Base Class for matrix graph
class BaseGraph
{
public:
  //! Constructor
  /* \param asize number of nodes (edges)
  */
  BaseGraph(Integer asize);

  // Destructor
  virtual ~BaseGraph();

  //! Set the eolemnt position in the matzrix graph
  /*! \param connect contains the global node (edge) numbers
      \param elemsize number of nodes (edges) in array connect
  */
  void SetElementPos(Integer * connect, Integer elemsize);

  //! perform sorting of graph array
  void Create();

  //! get the node (edge) neighbors for node(edge) i
  Integer * GetGraphRow(Integer i) const {return graph[i-1].pos;};

  //! get numbers of neighbors for node (edge) i
  Integer GetRowSize(Integer i) const {return graph[i-1].actsize;};

  //! get total numbers of node (edge) in graph
  Integer GetSize() const {return size;};

  //! get number of nonzero entries 
  Integer GetNumNNE() const {return nne;};

  void Print() const;

private:
  //! Structure for one element in matrix graph
  /*!
    \param pos stores the neigboring elements
    \param actsize actual size for neighbors of one element after ordering
    \param size predefined value for neighbors of one element (=50)
  */
  struct Graph
  {
    Integer * pos;
    Integer actsize;
    Integer size;
  };

  //! matrix graph
  Graph * graph;
  //! total number of elements (nodes, edges) in matrix graph
  Integer size;
  //! number of nonzero entries
  Integer nne;
  //! array for storing constraints-data?
  Integer * cm;
};

}

#endif // FILE_BASEGRAPH_CLA
