#ifndef FILE_MIXEDMATRIX_CLA
#define FILE_MIXEDMATRIX_CLA

namespace CoupledField
{

class MixedMatrix : public BaseMatrix
{
public:
  ///
  MixedMatrix();

  ///
  virtual ~MixedMatrix();

  ///
  virtual void Mult(BaseVector & vec1, BaseVector & vec2, Double factor) const;

  ///
  virtual void Mult(Double * vec1, Double * vec2, Double factor) const {;};

  ///
  virtual void Assemble(Double * v,Integer * p, Integer elemsize);

  ///
  virtual void Print() const;

  ///
  virtual void BuildInDirichlet() {;};

  ///
  virtual void BuildInDirichlet(BaseVector & rhs){;};

  ///
  virtual void Factor(){;};

  ///
  virtual void Solve(BaseVector & rhs, BaseVector & sol){;};

  ///
  virtual void Galerkin(BaseTransfer * ares, BaseTransfer * apro, BaseMatrix & amat){;};
  
  ///
  virtual void SetAuxMatrix(BaseMatrix & amat){;};

  ///
  virtual void Copy(BaseMatrix * mat) {;};

  ///
  virtual void ConstructEffectiveMatrix(BaseMatrix ** amat, Double * matrix_fac) {;};

private:
  ///
  BaseMatrix * sysmat[4];

};

}

#endif // FILE_MIXEDMATRIX_CLA
