#ifndef FILE_BASEGRAPH
#define FILE_BASEGRAPH

namespace CoupledField
{

class BaseGraph
{
public:
  ///
  BaseGraph(Integer asize);

  ///
  virtual ~BaseGraph();

  ///
  void SetElementPos(Integer * connect, Integer elemsize);

  ///
  void Create();

  ///
  Integer * GetGraphRow(Integer i) const {return graph[i-1].pos;};

  ///
  Integer GetRowSize(Integer i) const {return graph[i-1].actsize;};

  ///
  Integer GetSize() const {return size;};

  ///
  Integer GetNumNNE() const {return nne;};

  void Print() const;

private:
  ///
  struct Graph
  {
    Integer * pos;
    Integer actsize;
    Integer size;
  };

  Graph * graph;
  Integer size;
  Integer nne;
  Integer * cm;
};

}

#endif // FILE_BASEGRAPH
