#ifndef FILE_DENSEMATRIX_PILES
#define FILE_DENSEMATRIX_PILES

namespace CoupledField
{

class DenseMatrix
{
public:
  ///
  DenseMatrix(Integer asize, Integer awidth, Boolean mem = FALSE);

  ///
  virtual ~DenseMatrix();

  ///
  void Mult(DenseMatrix * a, DenseMatrix * b);

  ///
  void MultT(DenseMatrix * a, DenseMatrix * b);

  ///
  Double MaxNorm() const;

  ///
  Double Get(Integer i, Integer j) const {return val[(i-1)*size+j-1];};

  ///
  Double & Elem(Integer i, Integer j) {return val[(i-1)*size+j-1];};

  ///
  void Inverse(DenseMatrix * a);

  ///
  void Print() const;

  ///
  void InitMatrix();

  ///
  DenseMatrix & operator=(Double * x);

private:
  ///
  Integer size, width;

  ///
  Double * val, * fac, * u, * v, * y;
};

}

#endif // FILE_DENSEMATRIX_PILES

