#ifndef FILE_FULLVECTOR
#define FILE_FULLVECTOR

namespace CoupledField
{

class RealVector : public BaseVector
{
public:
  ///
  RealVector(Integer asize, Integer anumrhs, Integer adof, Boolean mem=TRUE);

  ///
  virtual ~RealVector();

  ///
  virtual void Add(DenseVector & a, BaseVector &vec1, DenseVector & b, BaseVector &vec2);

  ///
  virtual void Add(DenseVector & a, BaseVector &vec);
  
  ///
  virtual void Add(BaseVector &vec);

  ///
  virtual void Inner(BaseVector &vec, DenseVector &sum) const;

  ///
  virtual void Scal(DenseVector &sum);

  ///
  virtual DenseVector & L2Norm() const;

  ///
  virtual DenseVector & Max() const;

  ///
  virtual void Rand();

  ///
  virtual void Print() const;

  ///
  virtual void Assemble(Double * v, Integer * p, Integer elemsize);

  ///
  virtual void Init();

  ///
  virtual void InitOne();

  ///
  Double Get(Integer i, Integer j) {return val[(i-1)*numrhs+j-1];};

  ///
  Double Get(Integer i) {return val[i-1];};

  ///
  Double & Elem(Integer i) {return val[i-1];};
};



class ComplexVector : public BaseVector
{
public:
  ///
  ComplexVector(Integer asize, Integer anumrhs, Integer adof, Boolean mem=TRUE);

  ///
  virtual ~ComplexVector();

  ///
  virtual void Add(DenseVector & a, BaseVector &vec1, DenseVector & b, BaseVector &vec2);

  ///
  virtual void Add(DenseVector & a, BaseVector &vec);
  
  ///
  virtual void Add(BaseVector &vec);

  ///
  virtual void Scal(DenseVector &sum);

  ///
  virtual void Inner(BaseVector &vec, DenseVector &sum) const;

  ///
  virtual DenseVector & L2Norm() const;

  ///
  virtual DenseVector & Max() const;

  ///
  virtual void Rand();

  ///
  virtual void Print() const;

  ///
  virtual void Assemble(Double * v, Integer * p, Integer elemsize);

  ///
  virtual void Init();

  ///
  virtual void InitOne() {;};

  ///
  Double Get(Integer i, Integer j) {return val[(i-1)*numrhs+j-1];};

  ///
  Double Get(Integer i) {return val[i-1];};

  ///
  Double & Elem(Integer i) {return val[i-1];};

};

}

#endif // FILE_FULLVECTOR
