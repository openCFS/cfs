#ifndef FILE_FULLMATRIX_PILES
#define FILE_FULLMATRIX_PILES

namespace CoupledField
{

class RFullMatrix : public BaseMatrix
{
public:
  ///
  RFullMatrix(Integer asize, Integer adir, Boolean mem=TRUE);

  ///
  virtual ~RFullMatrix();

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
  Double Get(Integer i, Integer j) const {return val[(i-1)*size+j-1];};

  ///
  Double GetDiag(Integer i) const {return val[(i-1)*size+i-1];};

  ///
  Double & Elem(Integer i, Integer j) {return val[(i-1)*size+j-1];};

  ///
  virtual void SetCoarseGraph(BaseTopology * topology){;};

private:
  ///
  Double * fac, * y;
};



class CFullMatrix : public BaseMatrix
{
public:
  ///
  CFullMatrix();

  ///
  virtual ~CFullMatrix();

  ///
  virtual void Mult(BaseVector & vec1, BaseVector & vec2, Double factor) const;

  ///
  virtual void Mult(Double * vec1, Double * vec2, Double factor) const {;};

  ///
  virtual void Assemble(Double * v,Integer * p,Integer elemsize);

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

#endif // FILE_FULLMATRIX_PILES
