#ifndef FILE_BASEMATRIX_PILES
#define FILE_BASEMATRIX_PILES

namespace CoupledField
{

class BaseMatrix
{
public:
  ///
  BaseMatrix();

  ///
  virtual ~BaseMatrix();

  ///
  virtual void Mult(BaseVector & vec1, BaseVector & vec2, Double factor = 1) const  = 0;

  ///
  virtual void Mult(Double * vec1, Double * vec2, Double factor = 1) const = 0;
  
  ///
  virtual void Assemble(Double * v, Integer * p, Integer elemsize) = 0;

  ///
  virtual void Galerkin(BaseTransfer * ares, BaseTransfer * apro, BaseMatrix & amat) = 0;

  ///
  virtual void SetAuxMatrix(BaseMatrix & amat) = 0;

  ///
  Integer GetNumDoF() const {return dof;};

  ///
  Integer GetSize() const {return size;};

  ///
  Integer GetLength() const {return length;};

  ///
  Integer GetNNE() const {return nne;};

  ///
  Integer GetRowSize(Integer row) const {return (start[row]-start[row-1]);};

  ///
  Integer GetGraphPos(Integer row, Integer col) const;

  ///
  Integer GetMatrixPos(Integer row, Integer k) const {return pos[start[row-1]+k-1];};

  ///
  void SetStart(Integer row, Integer k) {start[row-1] = k;};

  ///
  void SetGraphRow(Integer end, Integer * apos, Integer row);

  ///
  virtual void Print() const = 0;

  ///
  virtual void BuildInDirichlet() = 0;
  
  ///
  virtual void BuildInDirichlet(BaseVector & rhs) = 0;

  ///
  virtual void Factor() = 0;

  ///
  virtual void Solve(BaseVector & rhs, BaseVector & sol) = 0;

  ///
  void SetDirichlet(Integer i, Integer num, Double v, Integer comp) 
    {numdir[i-1]=num; valdir[i-1]=v; compdir[i-1]=comp;};

  ///
  void UpdateDirichlet(Integer num, Double v)
    {valdir[num-1] = v;};

  ///
  void SetConstraint(Integer * cm, Double * cf){;};

  ///
  Integer GetNumDirichlet() const {return dir;};

  ///
  Integer GetIndDirichlet(Integer p) const {return numdir[p-1];};

  ///
  Integer GetValDirichlet(Integer p) const {return valdir[p-1];};

  ///
  Integer GetCmpDirichlet(Integer p) const {return compdir[p-1];};

  ///
  void SetMatrixFactor(Double amatfac){matfac = amatfac;};

  ///
  Integer * GetPosPointer() const {return pos;};

  ///
  Integer * GetStartPointer() const {return start;};

  ///
  void Calculated() {calculated = TRUE;};

  ///
  Boolean IsCalculated() const {return calculated;};

  ///
  void SetGraph() {setgraph = TRUE;};

  ///
  Boolean IsSetGraph() const {return setgraph;};

  ///
  virtual void SetCoarseGraph(BaseTopology * topology);

  ///
  Integer MaxRowSize() const;

  ///
  void InitMatrix();
  
  ///
  virtual void Copy(BaseMatrix * mat) = 0;

  ///
  void CopyGraph(Integer * astart, Integer * apos);

  ///
  void SetPointer(Integer * astart, Integer * apos, Double * aval)
    {start = astart; pos = apos; val = aval;};

protected:

  ///
  Double * val, * valdir, matfac;

  ///
  Integer * pos, * start, * numdir, * compdir;

  ///
  Integer * cm;

  ///
  Integer size, nne, dof, length, dir;

  ///
  Boolean calculated, setgraph, buildindir,outback;

  /// SuperLU
  SuperMatrix A, B, L, U;

  int * perm_c, * perm_r;
  double * f;
};

}

#endif // FILE_BASEMATRIX_PILES
