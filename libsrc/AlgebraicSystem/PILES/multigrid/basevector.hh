#ifndef FILE_BASEVECTOR_PILES
#define FILE_BASEVECTOR_PILES

namespace CoupledField
{

class BaseVector
{
public:
  ///
  BaseVector();

  ///
  virtual ~BaseVector();

  ///
  Integer GetNumRHS() const {return numrhs;};

  ///
  Integer GetSize() const {return size;};

  ///
  Integer GetNumDoF() const {return dof;};

  ///
  virtual void Assemble(Double * v, Integer * p, Integer elemsize) = 0;

  ///
  virtual void Add(DenseVector & a, BaseVector &vec1, DenseVector & b, BaseVector &vec2) = 0;

  ///
  virtual void Add(DenseVector & a, BaseVector &vec) = 0;
  
  ///
  virtual void Add(BaseVector &vec) = 0;

  ///
  virtual void Scal(DenseVector &sum) = 0;

  ///
  virtual void Inner(BaseVector &vec, DenseVector &sum) const = 0;

  ///
  virtual DenseVector & L2Norm() const = 0;

  ///
  virtual DenseVector & Max() const = 0;

  ///
  virtual void Rand() = 0;

  ///
  virtual void Print() const = 0;

  ///
  virtual void Init() = 0;

  ///
  virtual void InitOne() = 0;

  ///
  Double * GetPointer() {return val;};

  ///
  void SetPointer(Double * x) {val = x;};

protected:
  ///
  Double * val;

  ///
  Integer size, numrhs, dof, length;

  ///
  Boolean outback;
};

}

#endif // FILE_BASEVECTOR_PILES
