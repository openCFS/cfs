#ifndef FILE_DENSEVECTOR_CLA
#define FILE_DENSEVECTOR_CLA

namespace CoupledField
{

// class DenseVector
// {
// public:
//   ///
//   DenseVector();

//   ///
//   virtual ~DenseVector();

//   ///
//   virtual void Print() const = 0;

// protected:
//   ///
//   Integer size;

//   ///
//   Double * val;
// };

class DenseVector
{
public:
  ///
  DenseVector(Integer asize);

  ///
  virtual ~DenseVector();

  ///
  Double Get(Integer i) const {return val[i-1];};

  ///
  Double & Elem(Integer i) {return val[i-1];};

  ///
  void Sqrt();

  ///
  void Abs();

  ///
  DenseVector & operator=(const DenseVector &src);

  ///
  DenseVector & operator=(const Double x);

  ///
  DenseVector & operator*=(const Double x);

  ///
  DenseVector operator-() const;

  ///
  DenseVector operator/(const DenseVector &v) const;

  ///
  void Div(DenseVector &vec1, DenseVector &vec2);

  ///
  Boolean operator<=(const DenseVector &x) const;

  ///
  Boolean operator<=(const Double x) const;

  ///
  virtual void Print() const;

private:
  ///
  Integer size;

  ///
  Double * val;
};


// class CDenseVector : public DenseVector
// {
// public:
//   ///
//   CDenseVector(Integer asize);

//   ///
//   virtual ~CDenseVector();

//   ///
//   Double * Get(Integer i) const {return val[i-1];};

//   ///
//   Double & Elem(Integer i, Integer j) {return val[(i-1)*2+j-1];};

//   ///
//   void Sqrt();

//   ///
//   void Abs();

//   ///
//   CDenseVector & operator=(const CDenseVector &src);

//   ///
//   CDenseVector & operator=(const Double * x);

//   ///
//   CDenseVector & operator*=(const Double * x);

//   ///
//   CDenseVector operator-() const;

//   ///
//   CDenseVector operator/(const CDenseVector &v) const;

//   ///
//   Boolean operator<=(const CDenseVector &x) const;

//   ///
//   Boolean operator<=(const Double x) const;

//   ///
//   virtual void Print() const;
// };

}

#endif // FILE_DENSEVECTOR_CLA

