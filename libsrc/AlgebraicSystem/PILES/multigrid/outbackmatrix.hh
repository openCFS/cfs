#ifndef FILE_OUTBACKMATRIX_PILES
#define FILE_OUTBACKMATRIX_PILES

namespace CoupledField
{

class OutbackMatrix : public BaseMatrix
{
public:
  ///
  OutbackMatrix(Integer asize, Integer afsize, Integer adir, Integer alevel);

  ///
  virtual ~OutbackMatrix();

  ///
  virtual void Mult(BaseVector & vec1, BaseVector & vec2, Double factor) const;

  ///
  virtual void Mult(Double * vec1, Double * vec2, Double factor) const {;};

  ///
  virtual void Assemble(Double * v,Integer * p,Integer elemsize){;};

  ///
  virtual void Print() const;

  ///
  virtual void BuildInDirichlet();

  ///
  virtual void BuildInDirichlet(BaseVector & rhs);

  ///
  virtual void Factor();

  ///
  virtual void Solve(BaseVector & rhs, BaseVector & sol);
  
  ///
  virtual void Galerkin(BaseTransfer * ares, BaseTransfer * apro, BaseMatrix & amat);

  ///
  virtual void SetAuxMatrix(BaseMatrix & amat);

  ///
  virtual void Copy(BaseMatrix * mat);

  ///
  virtual void SetCoarseGraph(BaseTopology * topology){;};

  ///
  Double & Diag(Integer p) {return diag[p-1];};

  ///
  Double GetDiag(Integer p) {return diag[p-1];};

  ///
  void CreateTransfer(Integer asize, Integer * arsw)
    {accu = new RScalarTransfer(asize, arsw);};

  ///
  void SetTransfer(BaseTransfer * trans)
  {accu = trans;};

  ///
  BaseTransfer * GetTransfer() {return accu;};

  ///
  virtual void SetMatrixRow(Integer p, Integer q, Integer * apos, Double * aval){;};

private:
  ///
  Integer level, fsize, * rsw;

  long hsize;
  
  ///
  Double * diag;

  ///
  BaseTransfer * accu;

  ///
  BaseVector * r, * s;
};
}

#endif // FILE_OUTBACKMATRIX_PILES
