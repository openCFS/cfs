#ifndef FILE_BLOCKMATRIX_PILES
#define FILE_BLOCKMATRIX_PILES

namespace CoupledField
{

class RBlockMatrix : public BaseMatrix
{
public:
  ///
  RBlockMatrix(Integer asize, Integer anne, Integer adof, Integer adir);

  ///
  virtual ~RBlockMatrix();

  ///
  virtual void Mult(BaseVector & vec1, BaseVector & vec2, Double factor) const;

  ///
  virtual void Mult(Double * vec1, Double * vec2, Double factor) const {;};

  ///
  virtual void Assemble(Double * v,Integer * p, Integer elemsize);

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
  void SetDiagInv();

  ///
  Double & Elem(Integer p, Integer q, Integer r, Integer s) 
    {return val[(start[p-1]+q-1)*dof*dof+(r-1)*dof+s-1];};

  ///
  DenseMatrix * Get(Integer i, Integer j) const 
    {*dm = &val[(start[i-1]+j-1)*dof*dof]; return dm;};

  ///
  DenseMatrix * GetDiag(Integer i) const
    {*dm = &val[start[i-1]*dof*dof]; return dm;};

  ///
  DenseMatrix * GetDiagInv(Integer i) const
    {*dm = &diaginv[(i-1)*dof*dof]; return dm;};

private:
  ///
  DenseMatrix * dm;

  ///
  Double * diaginv, * sum;
};


  /////////////////////////////////////////////////////////////////////////////////////


class CBlockMatrix : public BaseMatrix
{
public:
  ///
  CBlockMatrix();

  ///
  virtual ~CBlockMatrix();

  ///
  virtual void Mult(BaseVector & vec1, BaseVector & vec2, Double factor) const;

  ///
  virtual void Mult(Double * vec1, Double * vec2, Double factor) const {;};

  ///
  virtual void Assemble(Double * v,Integer * p, Integer elemsize);

  ///
  virtual void Print() const;

  ///
  virtual void BuildInDirichlet(){;};

  ///
  virtual void BuildInDirichlet(BaseVector & rhs){;};

  ///
  virtual void Factor(){;};

  ///
  virtual void Solve(BaseVector & rhs, BaseVector & sol){;};
  
  ///
  virtual void Galerkin(BaseTransfer * ares, BaseTransfer * apro, BaseMatrix & amat) {;};
  
  ///
  virtual void SetAuxMatrix(BaseMatrix & amat){;};

  ///
  virtual void Copy(BaseMatrix * mat) {;};

};
}

#endif // FILE_BLOCKMATRIX_PILES
