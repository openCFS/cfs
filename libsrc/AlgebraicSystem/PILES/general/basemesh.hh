#ifndef FILE_BASEMESH_PILES
#define FILE_BASEMESH_PILES

namespace CoupledField
{
class BaseMesh
{
public:
  ///
  BaseMesh();

  ///
  virtual ~BaseMesh();

  ///
  virtual void Read() = 0;

  ///
  Integer GetNumNode() const {return numnode;};

  ///
  Integer GetNumElem() const {return numelem;};

  ///
  Integer GetElemSize() const {return elemsize;};

  ///
  Integer GetNumDir() const {return numdir;};

  ///
  Integer GetNumMat() const {return nummat;};

  ///
  Integer GetDim() const {return dim;};

  ///
  Integer * GetConnect(Integer i) const {return &connect[(i-1)*elemsize];};

  ///
  Double * GetCoord(Integer i) const {return &coord[(i-1)*dim];};

  ///
  Integer GetDirNode(Integer i) const {return bnode[i-1];};

  ///
  Double GetDirVal(Integer i) const {return bval[i-1];};

  ///
  Integer GetMaterial(Integer i) const {return mat[i-1];};

  ///
  Integer GetNNE() const {return nne;};
  
  ///
  Integer * GetGraphRow(Integer i) const {return &graph[(i-1)*offset];};

  ///
  Integer GetRowSize(Integer i) const {return rowsize[i-1];};

  ///
  void Print() const;

  ///
  void PrintMETIS() const;

protected:
  ///
  Char * filename, * pathname;

  ///
  Double * coord, * bval;

  ///
  Integer * connect, * mat, * bnode, * graph, * rowsize;

  ///
  Integer numnode, numelem, numdir, numneu, nummat, dim, elemsize, offset, nne;

  void ConstructGraph();

};

}

#endif // FILE_BASEMESH_PILES
