#ifndef FILE_SCALARMATRIX_PILES
#define FILE_SCALARMATRIX_PILES

namespace CoupledField
{

class RScalarMatrix : public BaseMatrix
{
public:
  ///
  RScalarMatrix(Integer asize, Integer anne, Integer adir);

  ///
  virtual ~RScalarMatrix();

  ///
  virtual void Mult(BaseVector & vec1, BaseVector & vec2, Double factor) const;

  ///
  virtual void Mult(Double * vec1, Double * vec2, Double factor) const {;};

  ///
  virtual void MultAdd(Double * vec1, BaseVector &vec2) const;

  ///
  virtual void Assemble(Double * v,Integer * p,Integer elemsize);

  ///
  virtual void Print() const;
  
  ///
  Double GetDiag(Integer i) const {return val[start[i-1]];};

  ///
  Double Get(Integer i, Integer j) const {return val[start[i-1]+j-1];};

  ///
  Double MinRowVal(Integer i) const;

  ///
  Double & Elem(Integer p, Integer q) {return val[start[p-1]+q-1];};

  ///
  virtual void BuildInDirichlet();

  ///
  virtual void BuildInDirichlet(BaseVector & rhs);

  ///
  virtual void UpdateDirichletRHS(BaseVector & rhs);

  ///
  virtual void Factor();

  ///
  virtual void Solve(BaseVector & rhs, BaseVector & sol);

  ///
  virtual void Galerkin(BaseTransfer * ares, BaseTransfer * apro, BaseMatrix & amat);

  ///
  virtual void SetAuxMatrix(BaseMatrix & amat){;};

  ///
  virtual void Copy(BaseMatrix * mat);

  ///
  virtual void ConstructEffectiveMatrix(BaseMatrix ** amat, Double * matrix_fac);

};


  /////////////////////////////////////////////////////////////////////////////////////


class CScalarMatrix : public BaseMatrix
{
public:
  ///
  CScalarMatrix(Integer asize, Integer anne, Integer adir);

  ///
  virtual ~CScalarMatrix();

  ///
  virtual void Mult(BaseVector & vec1, BaseVector & vec2, Double factor) const;

  ///
  virtual void Mult(Double * vec1, Double * vec2, Double factor) const {;};

  ///
  virtual void Assemble(Double * v,Integer * p,Integer elemsize);

  ///
  virtual void Print() const;

  ///
  virtual void BuildInDirichlet();

  ///
  virtual void BuildInDirichlet(BaseVector & rhs);

  ///
  virtual void UpdateDirichletRHS(BaseVector & rhs){;};

  ///
  Double * GetDiag(Integer i) const {return &val[2*start[i-1]];};  

  ///
  virtual void Factor();

  ///
  virtual void Solve(BaseVector & rhs, BaseVector & sol);

  ///
  Double * Get(Integer i, Integer j) const {return &val[2*start[i-1]+2*(j-1)];};

  ///
  virtual void Galerkin(BaseTransfer * ares, BaseTransfer * apro, BaseMatrix & amat);

  ///
  virtual void SetAuxMatrix(BaseMatrix & amat);

  ///
  virtual void Copy(BaseMatrix * mat);

  ///
  virtual void ConstructEffectiveMatrix(BaseMatrix ** amat, Double * matrix_fac) {;};
};

}

#endif // FILE_SCALARMATRIX_PILES
