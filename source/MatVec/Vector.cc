#include <def_use_embedded_python.hh>

#include "Vector.hh"
#include "opdefs.hh"
#include "Matrix.hh"

#include <boost/type_traits/is_complex.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

#include <cstdio>
#include <cmath>
#include <cfloat>
#include <fstream>

namespace CoupledField {

  // **********************************************
  //   Constructor for vector of specified length
  // **********************************************
  template<typename T>
  Vector<T>::Vector(const unsigned int size, const T entry) :
    SingleVector(), data_(new T[size]), capacity_(size), memBelongsToMe_(true)
  {
    size_ = size;
    for(unsigned int i = 0; i < size; ++i)
      data_ [i] = entry;
  }

  // *******************
  //   Deep destructor
  // *******************
  template<typename T>
  Vector<T>::~Vector()
  {
    if(memBelongsToMe_)
      delete[] data_;
  }


  // **********************
  //   Replace data array
  // **********************
  template<typename T>
  void Vector<T>::Replace( unsigned int length, T* entries, bool transferMem )
  {
    // De-allocate old array, if required
    if(memBelongsToMe_)
      delete[] data_;

    // Re-set internal attributes
    data_           = entries;
    size_           = length;
    capacity_       = length;
    memBelongsToMe_ = transferMem;
  }


  // *********************************
  //   Change the size of the vector
  // *********************************
  template<typename T>
  void Vector<T>::Resize( unsigned int newSize ) {
    if(!memBelongsToMe_ && newSize != size_ )
      EXCEPTION("Refusing to re-size the data_ array, since the memory does not belong to me!");
    
    // If new vector is not longer
    // than the currently allocated maximum
    // only adapt size_
    if(newSize <= capacity_)
    {
      size_ = newSize;
      return;
    }

    // New vector is longer than the currently
    // allocated data array, so we re-allocate
    delete[] data_;
    data_ = new T[newSize];
    size_ = newSize;
    capacity_ = newSize;
  }
  
  template<typename T>
  void Vector<T>::Resize(const unsigned int newSize, const T val)
  {
    Vector<T>::Resize(newSize);
    Vector<T>::Init(val);
  }

  template<typename T>
  void Vector<T>::Fill(const T* source, unsigned int new_size)
  {
    Resize(new_size);
    std::copy_n(source, new_size, data_);
  }

  template<typename T>
  void Vector<T>::Fill(const T& start, const T& increment)
  {
    if(size_ > 0) data_[0] = start;
    for(unsigned int i = 1; i < size_; i++)
      data_[i] = data_[i-1] + increment;
  }

  // *********************************
  //   Add vec to this vector object
  // *********************************
  template <typename T>
  void Vector<T>::Add(const SingleVector &vec)
  {
    const Vector<T>& idvec = dynamic_cast<const Vector<T>&>(vec);

    for(Integer i = 0; i < (Integer) size_; ++i)
      data_[i] += idvec[i];
  }


  // ***************************************************************
  //   Replace this vector object by the sum of two scaled vectors
  // ***************************************************************
  template <typename T>
  void Vector<T>::Add(T a, const SingleVector &vec1, 
		      T b, const SingleVector &vec2)
  {
    const Vector<T>& idvec1 = dynamic_cast<const Vector<T>&>(vec1);
    const Vector<T>& idvec2 = dynamic_cast<const Vector<T>&>(vec2);

    for(Integer i = 0; i < (Integer) size_; ++i)
      data_[i] = a * idvec1[i] + b * idvec2[i];	
  }
  
  
  // **********************************************************
  //   Add a scaled version of a vector to this vector object
  // **********************************************************
  template <typename T>
  void Vector<T>::Add(T a, const SingleVector &vec)
  {
    assert(vec.GetSize() == this->GetSize());
    const Vector<T>& idvec = dynamic_cast<const Vector<T>&>(vec);

    for(Integer i = 0; i < (Integer) size_; ++i)
      data_[i] += a * idvec[i];
  }
  
  template <typename T>
  void Vector<T>::Set(T a, const SingleVector &vec)
  {
    assert(vec.GetSize() == this->GetSize());
    const Vector<T>* ptr = dynamic_cast<const Vector<T>*>(&vec);
    assert(ptr != NULL);
    const Vector<T>& idvec = *ptr;

    for(unsigned int i = 0; i < size_; ++i)
      data_[i] = a * idvec[i];
  }

  template <typename T>
  void Vector<T>::Hadamard(const Vector<T>& v1, const Vector<T>& v2)
  {
    assert((v1.GetSize() == v2.GetSize()) && (v1.GetSize() == GetSize()));

    for(unsigned int i = 0; i < v1.GetSize(); i++)
      data_[i] = v1[i] * v2[i];
  }


  // *********
  //   Inner
  // *********

  template <typename T>
  void Vector<T>::Inner(const SingleVector &vec, T &sum) const
  {
    const Vector<T>& secVec = dynamic_cast<const Vector<T>&>(vec);
    sum = 0;

    for(unsigned int i = 0; i < size_; ++i)
      sum += OpType<T>::dotProduct(data_[i], secVec[i]);
  }

  template <typename T>
  T Vector<T>::Inner(const SingleVector &vec) const 
  {
    T sum(0);

    const Vector<T>& secVec = dynamic_cast<const Vector<T>&>(vec);

    for(unsigned int i = 0; i < size_; i++)  
      sum += OpType<T>::dotProduct(data_[i], secVec[i]);

    return sum;
  }

  template <typename T>
  T Vector<T>::Inner(const SingleVector& vec, unsigned int start, unsigned int end) const
  {
    T sum(0);

    const Vector<T>& secVec = dynamic_cast<const Vector<T>&>(vec);

    for(unsigned int i = start; i < end; i++)
      sum += OpType<T>::dotProduct(data_[i], secVec[i]);

    return sum;
  }

  template <typename T>
  T Vector<T>::Inner() const
  {
    T sum(0);

    for(unsigned int i = 0; i < size_; i++)
      sum += OpType<T>::dotProduct(data_[i], data_[i]);

    return sum;
  }
  
  template<typename T>
  Vector<T> Vector<T>::Conj() const {
    return *this;
  }

  template<>
  Vector<Complex> Vector<Complex>::Conj() const {
    Vector<Complex> ret(size_);
    for( UInt i = 0; i < size_; ++i )
      ret[i] = std::conj(data_[i]);

    return ret;
  }

  /*
  template <typename T>
  void Vector<T>::Inner(const SingleVector& vec,Double& s) const {
    if(typeid(Double) != typeid(T))
      EXCEPTION("Wrong typeids in Vector<T>::Inner ("
                << typeid(Double).name() << " != "
                << typeid(T).name() << ")");
  }

  template <typename T>
  void Vector<T>::Inner(const SingleVector& vec,Complex& s) const {
    if(typeid(Complex) != typeid(T))
      EXCEPTION("Wrong typeids in Vector<T>::Inner ("
                << typeid(Complex).name() << " != "
                << typeid(T).name() << ")");
  }
  */

  // ***********************************
  //   Inner angle between two vectors
  // ***********************************
  template<typename T>
  Double Vector<T>::InnerAngle(const Vector<T>& other, bool preferPositiveZ) const
  {
    EXCEPTION("Inner angle only for Vector<Double>.")
  }
  
  template<>
  Double Vector<Double>::InnerAngle(const Vector<Double>& other, bool preferPositiveZ) const
  {
    /*
     * Simple helper function computing the angle between two vectors
     * 
     * > no optimization, 2d and 3d only, basically following:
     * https://stackoverflow.com/questions/14066933/direct-way-of-computing-clockwise-angle-between-2-vectors
     *
     * > in 2d (x,y) the returned angle is measured in a counterclockwise fashion from "this" towards "other" rotating
     *    around the positve z-axis
     *   i.e., the returned angle is negative if "this" has a larger counterclockwise angle than "other"
     * > in 3d the returned angle is always positive as the rotation axis is not fix; per definition the one
     *    that leads to a positive angle is utilized (see source above)
     * > by setting the newly introduced flag "preferPositiveZ" (default = false) to true,
     *    the actual rotation axis is evaluated from the cross product of "this" and "other" and 
     *    its z-component will be checked; if it is negative, the rotation axis will get reverted
     *    and consequntly the returned angle will be multiplied by -1;
     *    Obviously, this does not work if the normal has a zero value z-component; however it ensures
     *    that the angle between two vectors inside the x-y-plane has the same sign whether it is
     *    evaluated in 2d or 3d
     * 
     * Returned angle in radians
     */
    assert(size_ == other.GetSize());
    assert(size_ >= 2);
    assert(size_ <= 3);
    
    Double len1Sq = this->NormL2_squared();
    Double len2Sq = other.NormL2_squared();
    
    Double angleRad = 0.0;
    if((len1Sq == 0) || (len2Sq == 0)){
      return angleRad;
    }
    
    Double innerProduct = 0.0;
    this->Inner(other,innerProduct);
    
    if(size_ == 2){
      Double crossProduct2d = data_[0]*other[1] - data_[1]*other[0];   
      angleRad = std::atan2(crossProduct2d,innerProduct);
    } else {
      // floating point precision might cause argument to exceed bounds of allowed input for acos
      // furthermore, formula was wrong; had + instead of *
      Double argument = innerProduct/std::sqrt(len1Sq * len2Sq);
      if(argument>1.0){ argument = 1.0;}
      if(argument< -1.0){ argument = -1.0;}
      angleRad = std::acos(argument);
      
      if(preferPositiveZ){
        Double crossProduct3d_z = data_[0]*other[1] - data_[1]*other[0];
        if(crossProduct3d_z<0){
          angleRad = -angleRad;
        }
      }
    }

    return angleRad;
  }

  
  // **********************************
  //   Set all vector entries to zero
  // **********************************
  template<typename T>
  void Vector<T>::Init( T entry)
  {
    std::fill(data_, data_+size_, entry);
  }


  // **************************
  //   Compute Euclidean Norm
  // **************************
  template <typename T>
  inline Double Vector<T>::NormL2() const
  {				
    double sum = 0;

    //#pragma omp parallel for reduction(+:sum)
    for(Integer i = 0; i < (Integer) size_; ++i)
      sum += OpType<T>::zConjz(data_[i]);

    return sqrt(sum);
  }


  template <typename T>
  inline double Vector<T>::NormL2(const Vector<T>& other) const
  {
    if(size_ != other.GetSize()) EXCEPTION("incompatible sizes");

    double sum = 0;

    //#pragma omp parallel for reduction(+:sum)
    for(Integer i = 0; i < (Integer) size_; ++i)
      sum += OpType<T>::zConjz(data_[i] - other[i]);

    return sqrt(sum);
  }

  /*
   * sometimes NormL2^2 is needed (e.g. for trust region checking for
   * hysteresis; in these cases it makes no sense to take the sqrt and
   * then square it; just skip the sqrt
   */
  template <typename T>
  inline Double Vector<T>::NormL2_squared() const
  {
    double sum = 0;

    //#pragma omp parallel for reduction(+:sum)
    for(Integer i = 0; i < (Integer) size_; ++i)
      sum += OpType<T>::zConjz(data_[i]);

    return sum;
  }

  template <typename T>
  inline double Vector<T>::MAC(const Vector<T>& other) const
    {
      if(size_ != other.GetSize()) EXCEPTION("incompatible sizes");
      Vector<T> v1c =  this->Conj();
      Double up = std::abs(v1c.Inner(other));
      T down1 = (this->Conj()).Inner(other);
      T down2 = this->Inner(other.Conj());
      Complex down = Complex(down1*down2);
      return up*up/down.real();
    }

  template <>
  inline double Vector<UInt>::MAC(const Vector<UInt>& other) const {
    EXCEPTION("Not useful for unsigned integer");
  }

  template <typename T>
  inline T Vector<T>::Sum() const
  {
    T s(0);
    for(unsigned int i = 0; i < size_; ++i)
      s+=data_[i];

    return s;
  }

  template <typename T>
  inline T Vector<T>::Product() const
  {
    T s(1);
    for(unsigned int i = 0; i < size_; ++i)
      s*=data_[i];

    return s;
  }


  template <typename T>
  inline T Vector<T>::Avg() const
  {
    assert(size_ > 0);

    return Sum() * (1.0/size_);
  }


  template <typename T>
  inline T Vector<T>::Min() const
  {
    assert(size_ > 0);

    T m = data_[0];
    //unsigned int max_size = *(std::max_element(&n_.First(), &n_.Last()));
    for(unsigned int i = 1; i < size_; ++i)
      m = std::min(m, data_[i]);

    return m;
  }

  template <>
  Complex Vector<Complex>::Min() const
  {
    assert(size_ > 0);
    Complex m = data_[0];

    for(unsigned int i = 1; i < size_; ++i) {
      if(data_[i].real() < m.real())
        m.real(data_[i].real());
      if(data_[i].imag() < m.imag())
        m.imag(data_[i].imag());
    }

    return m;
  }


  template <typename T>
  inline T Vector<T>::Max() const
  {
    assert(size_ > 0);

    T m = data_[0];

    for(unsigned int i = 1; i < size_; ++i)
      m = std::max(m, data_[i]);

    return m;
  }

  template <>
  Complex Vector<Complex>::Max() const
  {
    assert(size_ > 0);
    Complex m = data_[0];

    for(unsigned int i = 1; i < size_; ++i) {
      if(data_[i].real() > m.real())
        m.real(data_[i].real());
      if(data_[i].imag() > m.imag())
        m.imag(data_[i].imag());
    }

    return m;
  }

  template <typename T>
  inline T Vector<T>::MaxAbs(int & loc) const {
    T m = data_[0];
    double abs, m_abs = 0;
    for(int i = 0; i < (int)size_; ++i) {
      abs = (double) std::abs(data_[i]);
      if ( abs > m_abs ) {
        m_abs = abs;
        m = data_[i];
        loc = i;
      }
    }
    return m;
  }

  template <>
  inline UInt Vector<UInt>::MaxAbs(int & loc) const {
    EXCEPTION("Not useful for unsigned integer");
  }

  template <typename T>
  inline void Vector<T>::MinMax(T& min, T& max) const
  {
    assert(size_ > 0);

    min = data_[0];
    max = data_[0];

    for(unsigned int i = 1; i < size_; ++i) {
      min = std::min(min, data_[i]);
      max = std::max(max, data_[i]);
    }
  }

  template <>
  void Vector<Complex>::MinMax(Complex& min, Complex& max) const
  {
    min = Min();
    max = Max();
  }


  template<class TYPE> 
  Double Vector<TYPE>::NormMax() const 
  { 
    EXCEPTION( "Vector<TPYE>::NormMax only defined for TYPE=Complex/Double" ); 
    return TYPE(); 
  } 

  // this functions localized the maximal component (absolute value) and returns it with its original sign
  // example: SignedMax([1,0,0]) = 1; SignedMax([-1,0,0]) = -1
  template<class TYPE> 
  Double Vector<TYPE>::SignedMax() const 
  { 
    EXCEPTION( "Vector<TPYE>::SignedMax only defined for TYPE=/Double" ); 
    return TYPE(); 
  } 


  template<> 
  double Vector<Double>::SignedMax() const 
  { 
    double ret(0.0),ret_new(0.0); 
    int idx = 0;
    if(size_ == 0) EXCEPTION("empty vector"); 

    for(unsigned int i = 0; i < size_; ++i){
      ret_new = std::max(ret, std::abs(data_[i]));
      if(ret != ret_new){
        // std::abs(data_[i]) > ret -> save index
        idx = i;
      }
      ret = ret_new;
    }

    return data_[idx]; 
  }

  template<> 
  double Vector<Complex>::SignedMax() const 
  { 
    double ret(0.0),ret_new(0.0); 
    int idx = 0;
    if(size_ == 0) EXCEPTION("empty vector"); 

    for(unsigned int i = 0; i < size_; ++i) {
      ret_new = std::max(std::abs(ret), std::abs(data_[i].real())); 
      if(std::abs(ret) != std::abs(ret_new)){
        // std::abs(data_[i]) > ret -> save index
        idx = i;
        ret = ret_new;
      }
    }
    return data_[idx].real(); 
  } 

  template<> 
  double Vector<Double>::NormMax() const 
  { 
    double ret(0.0); 
    if(size_ == 0) EXCEPTION("empty vector"); 

    for(unsigned int i = 0; i < size_; ++i) 
      ret = std::max(ret, std::abs(data_[i]));

    return ret; 
  }
  
  template<> 
  double Vector<Complex>::NormMax() const 
  { 
    double ret(0.0); 
    if(size_ == 0) EXCEPTION("empty vector"); 

    for(unsigned int i = 0; i < size_; ++i) 
      ret = std::max(std::abs(ret), std::abs(data_[i].real())); 

    return ret; 
  }  

  template<class TYPE>
  Double Vector<TYPE>::NormMax(const SingleVector& other) const
  {
    EXCEPTION( "Vector<TPYE>::NormMax only defined for TYPE=Complex/Double" );
    return TYPE();
  }

  template<>
  double Vector<Double>::NormMax(const SingleVector& other) const
  {
    double ret(0.0);
    if(size_ == 0) EXCEPTION("empty vector");
    if(size_ != other.GetSize()) EXCEPTION("incompatible sizes");
    const Vector<Double>& o = dynamic_cast<const Vector<Double>&>(other);
    for(unsigned int i = 0; i < size_; ++i)
      ret = std::max(ret, std::abs(data_[i] - o.data_[i]));

    return ret;
  }

  template<>
  double Vector<Complex>::NormMax(const SingleVector& other) const
  {
    double ret(0.0);
    if(size_ == 0) EXCEPTION("empty vector");
    if(size_ != other.GetSize()) EXCEPTION("incompatible sizes");
    const Vector<Complex>& o = dynamic_cast<const Vector<Complex>&>(other);
    for(unsigned int i = 0; i < size_; ++i)
      ret = std::max(std::abs(ret), std::abs(data_[i].real() - o.data_[i].real()));

    return ret;
  }


  template<typename T>
  Double Vector<T>::Normalize() {
    EXCEPTION("Vector<TYPE>::Normalize() is only defined for TYPE=Float,Double,Complex");
  }

  template<>
  Double Vector<Float>::Normalize() {
    Double norm = NormL2();
    ScalarDiv(norm);
    return norm;
  }

  template<>
  Double Vector<Double>::Normalize() {
    Double norm = NormL2();
    ScalarDiv(norm);
    return norm;
  }

  template<>
  Double Vector<Complex>::Normalize() {
    Double norm = NormL2();
    ScalarDiv(norm);
    return norm;
  }

  template<typename T>
  Integer Vector<T>::CountNonZero() const
  {
    Integer count = 0;

    for(unsigned int i = 0; i < size_; i++)
      if(data_[i] != T())
        count++;

    return count;
  }


  // *******************************************
  //   Same as the BLAS functions of that name
  // *******************************************
  template <typename T>
  void Vector<T>::Axpy(const T alpha, const SingleVector &y)
  {
    const Vector<T>& vec = dynamic_cast<const Vector<T>&>(y);

#pragma omp parallel for 
    for ( Integer i = 0; i < (Integer) size_; ++i )
      data_[i] = alpha * data_[i] + vec[i];
  }


  // ******************
  //   SetVectorEntry
  // ******************
  template<typename T>
  void Vector<T>::SetEntry( const unsigned int i, const T &val )
  {
#ifdef DEBUG_VECTOR
    if ( i <= 0 || i > size_ ) {
      EXCEPTION( "Vector<>::SetVectorEntry: "
               << "Detected index error:"
               << "\n index = " << i
               << "\n size  = " << size_ );
    }
#endif
    data_[i] = val;
  }


  // ********************
  //   AddToVectorEntry
  // ********************
  template<typename T>
  void Vector<T>::AddToEntry( const unsigned int i, const T &val )
  {
    data_[i] += val;
  }


  // ********************
  //   Get vector entry
  // ********************
  template<typename T>
  void Vector<T>::GetEntry( const unsigned int i, T &val ) const
  {
    //assert(i<size_);
    val = data_[i];
  }

  template<typename T>
  T Vector<T>::Last() const
  {
    assert(size_ > 0);
    return data_[size_-1];
  }


  // ********************
  //   Get several entries
  // ********************
  template<typename T>
  const Vector<T> Vector<T>::GetEntries( const StdVector<UInt>& in) const
  {
    Vector<T> vals;
    vals.Resize(in.GetSize());
    for(UInt i = 0; i < in.GetSize(); ++i){
      vals[i] = data_[in[i]];
    }
    return vals;
  }

  template<typename T>
  void  Vector<T>::Push_back(const T & y)
  {
    // Check whether capacity is sufficiently large to perform
    // a push-back operation. If not allocate memory according
    // to the following simply scheme: Each time capacity is
    // exceeded allocate twice as much memory as there was before.
    if(size_ >= capacity_)
    {
      T *help;

      // perform memory allocation
      capacity_ = size_ == 0 ? 1 : 2 * size_;
      help = new T[capacity_];

      // copy old entries into new buffer
      for(unsigned int i = 0; i < size_; ++i)
        help[i] = data_[i];

      // delete old buffer and re-set pointer
      delete[] data_;
      data_ = help;
    }

    // Perform push-back and increase size
    data_[size_++] = y; 
  }


  template<typename T>
  Vector<Double> Vector<T>::GetPart(Global::ComplexPart part) const
  {
    EXCEPTION( "Matrix::GetPart: Only Implemented for Real and "
               << " Complex matrices!" );
    Vector<Double> temp;
    return temp;
  }


  template<>
  Vector<Double> Vector<Double>::GetPart(Global::ComplexPart part) const
  {
    if(part != Global::REAL)
      EXCEPTION("Vector<Double>::GetPart: Only possible for REAL part." );
    
    return *this;
  }

  template<>
  Vector<Double> Vector<Complex>::GetPart(Global::ComplexPart part) const
  {
    Vector<Double> ret(size_);
    
    switch(part)
    {
    case Global::REAL:
      for ( unsigned int i = 0; i < size_; ++i )
        ret[i]  = data_[i].real();
      break;
    case Global::IMAG:
      for(unsigned int i = 0; i < size_; ++i)
        ret[i]  = data_[i].imag();
      break;
    default:
      EXCEPTION("Vector<Complex>::GetPart: Only possible for REAL or IMAG part!" );
      break;
    }
    
    return ret;
  }
  
  template<typename T>
  void Vector<T>::SetPart( Global::ComplexPart part, const Vector<Double> & partVector,
                           bool zeroOtherPart )
  {
    EXCEPTION( "Vector::SetPart: Only Implemented for Real and Complex vectors!" );
  }

  template<>
  void Vector<Double>::SetPart( Global::ComplexPart part, const Vector<Double> & partVector,
                                bool zeroOtherPart )
  {
    if(size_ != partVector.GetSize())
      EXCEPTION( "Vector<Double>::SetPart: Dimension of vectors do not match!" );

    if(part != Global::REAL)
      EXCEPTION( "Vector<Double>::SetPart: Only possible for REAL part." );
    
    *this = partVector;
  }

  template<>
  void Vector<Complex>::SetPart( Global::ComplexPart part, const Vector<Double> & partVector, bool zeroOtherPart ) {
    if(size_ != partVector.GetSize())
      EXCEPTION( "Vector<Complex>::SetPart: Dimension of vectors do not match!" );

    if( zeroOtherPart) {
      // ------------------
      //  Zero other part
      // ------------------
      switch(part)  {
        case Global::REAL:
          for(unsigned int i = 0; i < size_; ++i)
            data_[i] = Complex(partVector[i], 0.0);
          break;
        case Global::IMAG:
          for(unsigned int i = 0; i < size_; ++i)
            data_[i] = Complex(0.0, partVector[i]);
          break;
        default:
          EXCEPTION( "Vector<Complex>::SetPart: Only possible for REAL or IMAG part!" );
          break;
      }
    } else {
      // ------------------
      //  Keep other part
      // ------------------
      switch(part)  {
        case Global::REAL:
          for(unsigned int i = 0; i < size_; ++i)
            data_[i] = Complex(partVector[i], data_[i].imag());
          break;
        case Global::IMAG:
          for(unsigned int i = 0; i < size_; ++i)
            data_[i] = Complex(data_[i].real(), partVector[i]);
          break;
        default:
          EXCEPTION( "Vector<Complex>::SetPart: Only possible for REAL or IMAG part!" );
          break;
      }
    }
  }

  template<typename T>

  std::string Vector<T>::ToString(ToStringFormat format, const std::string& sep_in, int digits) const
  {
    std::ostringstream os;

    std::string sep = sep_in != "" ? sep_in : (format == TS_PLAIN ? " " : ", ");

    // ignored for int case, we cannot fill ints easily as it makes no sense for neg and >= not possible for complex
    if(digits > 0)
      os << std::scientific << std::setprecision(digits);

    char cplx = format == TS_PYTHON ? 'j' : 'i';

    switch(format)
    {
    case TS_INFO: {
      int nnz=0;
      for(unsigned int i = 0; i < size_; ++i)
        if(Abs(data_[i]) != 0)
          nnz++;
      os << "size=" << size_ << " nnz=" << nnz;
      }
      break;

    case TS_MATLAB:
    case TS_PYTHON:
      os << "[";
      // intentionally no break;
    case TS_PLAIN: {
      for(unsigned int i = 0; i < size_; i++)
      {
        if(boost::is_complex<T>::value)
        {
          Complex cval = (Complex) data_[i];
          os << cval.real() << "+" << cval.imag() << cplx;
        }
        else
          os << data_[i];
        if(i < size_-1)
          os << sep;
      }
      if(format != TS_PLAIN)
        os << "]";
    }
    break;

    case TS_NONZEROS:
      for(unsigned int i = 0; i < size_; ++i)
        if(Abs(data_[i]) != 0)
          os << i << ":" << data_[i] << (i < size_-1 ? sep : "");
      break;
    } // end switch
    return os.str();
  }
  

  // *****************
  //   Export vector
  // *****************
  template<typename T>
  void Vector<T>::Export(const std::string& fname, BaseMatrix::OutputFormat format) const
  {
    std::stringstream sstr;
    BaseMatrix::EntryType eType = GetEntryType();
    UInt nnz;

    // Assemble file name depending on output format and entry type.
    sstr << fname;
    switch(format)
    {
    case BaseMatrix::MATRIX_MARKET: // MatrixMarket
      sstr << ".mtx";
      break;
    case BaseMatrix::HARWELL_BOEING: // Harwell-Boeing
      switch(eType) 
      {
      case BaseMatrix::DOUBLE:
        sstr << ".rra";
        break;
      case BaseMatrix::COMPLEX:
        sstr << ".cra";
        break;
      default:
        break;
      }
      break;
    case BaseMatrix::PLAIN:
      sstr << ".vec";
      break;
    }
    
    // open output file and check for errors
    FILE *fp = fopen( sstr.str().c_str(), "wb" );
    if(fp == NULL)
      EXCEPTION( "Cannot open file '" << fname << "' for writing!" );

    // Determine number of non-zero entries and non-zero rows.
    StdVector<UInt> nzRows;
    for(unsigned int k = 0; k < size_; ++k)
    {
      if( ! IsZero<T>( data_[k] ) ) {
        nzRows.Push_back(k);
      }
    }
    nnz = nzRows.GetSize();

    // Set number of rows and columns
    UInt nrows = size_;
    UInt ncols = 1;
    
    switch(format)
    {
    case BaseMatrix::MATRIX_MARKET: // MatrixMarket
    {
      switch(eType) 
      {
      case BaseMatrix::DOUBLE:
        fprintf( fp, "%%%%MatrixMarket matrix coordinate real general\n" );
        break;
      case BaseMatrix::COMPLEX:
        fprintf( fp, "%%%%MatrixMarket matrix coordinate complex general\n" );
        break;
      default:
        break;
      }
      
      // Print comment
      fprintf( fp, "%%\n%% Vector exported by openCFS\n%%\n" );
      
      // Information on number of rows, columns and entries
      fprintf( fp, "%d\t%d\t%d\n", nrows, ncols, nnz );
      for(unsigned int k = 0; k < nnz; ++k)
      {
        // store row and column index
        fprintf( fp, "%6d\t%6d\t", nzRows[k] + 1, 1);
        
        // store non-zero entry
        OpType<T>::ExportEntry( data_[nzRows[k]], 0, fp );
        fprintf( fp, "\n" );
      }
    } 
    break;
    case BaseMatrix::HARWELL_BOEING: // Harwell-Boeing
    {
      std::string code;
      std::string fmt;
      
      switch(eType)
      {
      case BaseMatrix::DOUBLE:
        code = "RRA           ";
        fmt  = "(1E24.16)           ";
        break;
      case BaseMatrix::COMPLEX:
        code = "CRA           ";
        fmt  = "(2E24.16)           ";
        break;
      default:
        break;
      }

      UInt ptrcrd = 1;
      UInt indcrd = nnz/8 + UInt(nnz % 8 != 0);
      UInt valcrd = nnz;
      UInt totcrd = 4 + ptrcrd + indcrd + valcrd;
      //fprintf( fp, "No label                                                                No key  \n");
      fprintf( fp, "Harwell-Boeing vector exported from openCFS                              No key  \n");
      fprintf( fp, "% 14d% 14d% 14d% 14d% 14d\n", totcrd, ptrcrd, indcrd, valcrd, 0 );
      fprintf( fp, "%s% 14d% 14d% 14d\n", code.c_str(), size_, 1, nnz );
      fprintf( fp, "(8I10)          (8I10)          %s\n", fmt.c_str());
      fprintf( fp, "% 10d% 10d\n", 1, (nnz+1) );
      
      UInt k;
      for(k = 0; k < nnz; )
      {
        // store row index
        fprintf( fp, "% 10d", (nzRows[k] + 1));
        
        k++;

        if( (k % 8) == 0) {
          fprintf( fp, "\n" );
        }
      }
      
      if( (k % 8) != 0) {
        fprintf( fp, "\n" );
      }

      for(k = 0; k < nnz; k++)
      {
        OpType<T>::ExportEntry( data_[nzRows[k]], 0, fp );
        fprintf( fp, "\n" );
      }
    }  
    break;
    case BaseMatrix::PLAIN:
    for(unsigned int k = 0; k < size_; ++k)
    {
      OpType<T>::ExportEntry(data_[k], 0, fp);
      fprintf(fp, "\n");
    }
    break;
    }
    
    // close output file
    if(fclose( fp ) == EOF)
    {
      WARN("Could not close file '" << fname << "' after writing!");
    }
  }

  // ********************
  //   Import vector
  // ********************
  template<typename T>
  void Vector<T>::Import(const std::string& fname, bool checkSize)
  {
    // Determine expected entry Type of the imported file
    std::string matrixType;
    if(GetEntryType() == BaseMatrix::DOUBLE)
      matrixType = "real";
    else if(GetEntryType() == BaseMatrix::COMPLEX)
      matrixType = "complex";
    else
      EXCEPTION( "Entry Type of the Matrix is not supported by the Import Method" );

    //adding the ending to the fileName
    std::stringstream sstr;
    sstr << fname << ".mtx";

    // open input file
    std::ifstream file(sstr.str());
    if(!file)
    {
      EXCEPTION( "Cannot open file '" << fname << "' for reading!" );
    }

    // Determine the matrix format
    int lineLength = 1024;
    std::string line;
    std::string fileType;
    std::string objectType;
    std::string formatType;
    std::string firstQualifier;

    //reading the header
    if(!(file >> fileType >> objectType >> formatType >> firstQualifier))
    {
      EXCEPTION("The header of the input file '" << fname << "'could not be read properly");
    }

    //jump to the next line
    file.ignore(lineLength, '\n');

    //check if the file is of type MatrixMarket format
    if(std::string(fileType) != "%%MatrixMarket")
    {
      EXCEPTION("The Input File '" << fname << "'has not the MatrixMarket format, first line must start with %%MatrixMarket");
    }

    if(objectType != "matrix")
    {
      EXCEPTION("Only object type 'matrix' of the MatrixMarket format is supported");
    }

    //Check if entry Type of the Matrix match with the entry type of the imported matrix
    if(firstQualifier != matrixType)
    {
      EXCEPTION("Entry Type does not match the type of the imported file");
    }

    // Ignore comments
    while (file.peek() == '%')
    {
      file.ignore(lineLength, '\n');
    }

    // Read in the line containing size info
    UInt rows;
    UInt columns;
    std::string lineInfo;
    if (file.peek() == EOF)
    {
      EXCEPTION("Imported File " << fname << " contains not enough entries");
    }
    std::getline(file, lineInfo);
    std::stringstream ssInfo(lineInfo);

    // Import the matrix entries depending on format type
    if(formatType == "array")
    {
      // Determine number of rows and columns
      if (!(ssInfo >> rows >> columns))
      {
        EXCEPTION("Unable to read the amount of rows and columns in the input file '" << fname << "'");
      }

      // Handle the case when the size of the vector is initially unknown
      if (checkSize)
      {
        if ((columns != 1) || (rows != size_))
        {
          EXCEPTION("Input Matrix should be a [" << size_ << " x 1] - matrix");
        }
      }
      else
      {
        if (columns != 1)
        {
          EXCEPTION("Input Matrix should have only 1 column");
        }
        Resize(rows);
      }

      T entry;
      for ( UInt i = 0; i < rows; i++)
      {
        if (file.peek() == EOF)
        {
          EXCEPTION("Imported File " << fname << " contains not enough entries");
        }
        std::getline(file, line);
        std::stringstream ssline(line);
        ReadSingleEntry(&entry,ssline);
        SetEntry(i, entry);
      }

    }
    else if(formatType == "coordinate")
    {
      //Determine number of rows, columns and non-zero entries
      UInt nonZero;
      if (!(ssInfo >> rows >> columns >> nonZero))
      {
        EXCEPTION("Unable to read the amount of rows, columns and non-zero entries in the input file '" << fname << "'");
      }

      // Handle the case when the size of the vector is initially unknown
      if (checkSize)
      {
        if ((columns != 1) || (rows != size_))
        {
          EXCEPTION("Input Matrix should be a [" << size_ << " x 1] - matrix");
        }
      }
      else
      {
        if (columns != 1)
        {
          EXCEPTION("Input Matrix should have only 1 column");
        }
        Resize(rows);
      }

      Init();
      UInt row;
      UInt column;
      T entry;
      for ( UInt i = 0; i < nonZero; i++)
      {
        if (file.peek() == EOF)
        {
          EXCEPTION("Imported File " << fname << " contains not enough entries");
        }
        std::getline(file, line);
        std::stringstream ssline(line);
        ssline >> row >> column;
        ReadSingleEntry(&entry,ssline);
        SetEntry(row-1, entry);
      }
    }
    else
    {
      EXCEPTION("Matrix format is not specified in the Import file");
    }

    // close input file
    file.close();
  }

  // ********************
  //   ScalarDiv (real)
  // ********************
  template<typename T>
  void Vector<T>::ScalarDiv(const Double factor)
  {
#ifdef DEBUG_VECTOR
    if(factor == 0)
    {
      EXCEPTION( "Vector::ScalarDiv: Division by Zero!" );
    }
#endif
    T fac(static_cast<T>(factor)); // for compilers
    for(unsigned int i = 0; i < size_; ++i)
      data_[i] /= fac;
  }


  // *********************
  //   ScalarMult (real)
  // *********************
  // NOTE: this function is used only once in whole cfs!!
  template<typename T>
  void Vector<T>::ScalarMult(const Double factor)
  {
    T fac(static_cast<T>(factor)); // for compilers
    for(unsigned int i = 0; i < size_; ++i)
      data_[i] *= fac;
  }


  // ***********************
  //   ScalarDiv (complex)
  // ***********************
  template<typename T>
  void Vector<T>::ScalarDiv(const Complex factor)
  {
    for(unsigned int i = 0; i < size_; ++i)
      OpType<T>::DivByComplex(data_[i], factor);
  }


  // ************************
  //   ScalarMult (complex)
  // ************************
  template<typename T>
  void Vector<T>::ScalarMult(const Complex factor)
  {
    for(unsigned int i = 0; i < size_; ++i)
      OpType<T>::MultWithComplex( data_[i], factor ); 
  }

  template<class TYPE>
  bool Vector<TYPE>::Contains(const TYPE val) const
  {
    for(UInt k = 0, s = size_; k < s; ++k)
      if(data_[k] == val) return true;

    return false;
  }


  template<class TYPE>
  bool Vector<TYPE>::ContainsNaN() const
  {
    for(UInt k = 0, s = size_; k < s; ++k)
      if(std::isnan(data_[k]))
        return true;

    return false;
  }

  template<>
  bool Vector<int>::ContainsNaN() const
  {
    return false;
  }

  template<>
  bool Vector<unsigned int>::ContainsNaN() const
  {
    return false;
  }

  template<>
  bool Vector<Complex>::ContainsNaN() const
  {
    for(UInt k = 0, s = size_; k < s; ++k)
    {
      if(std::isnan(data_[k].real()))
        return true;
      if(std::isnan(data_[k].imag()))
        return true;
    }
    return false;
  }


  template<class TYPE>
  bool Vector<TYPE>::ContainsInf() const
  {
    for(UInt k = 0, s = size_; k < s; ++k)
      if(std::isinf(data_[k]))
        return true;

    return false;
  }

  template<>
  bool Vector<int>::ContainsInf() const
  {
    return false;
  }

  template<>
  bool Vector<unsigned int>::ContainsInf() const
  {
    return false;
  }

  template<>
  bool Vector<Complex>::ContainsInf() const
  {
    for(UInt k = 0, s = size_; k < s; ++k)
    {
      if(std::isinf(data_[k].real()))
        return true;
      if(std::isinf(data_[k].imag()))
        return true;
    }
    return false;
  }


  template<typename T>
  inline bool Vector<T>::operator==(const Vector<T> &x) const {
    if ( this == &x ) return true; // we compare pointers, not references, therefore not recusively
    if ( size_ != x.size_ ) return false;

    // memcmp is significantly faster than looping manually:
    
    // vector compare: size=3 n=333333 opt=loop dt=00:00:00.008847
    // vector compare: size=3 n=333333 opt=memcmp dt=00:00:00.005975
    // vector compare: size=300 n=3333 opt=loop dt=00:00:00.002672
    // vector compare: size=300 n=3333 opt=memcmp dt=00:00:00.000158
    // vector compare: size=30000 n=33 opt=loop dt=00:00:00.002576
    // vector compare: size=30000 n=33 opt=memcmp dt=00:00:00.000252

    return memcmp(data_, x.data_, size_ * sizeof(T)) == 0 ? true : false;
  }

  template<typename T>
  inline bool Vector<T>::operator!=( const Vector<T>& x) const
  {
    if(this == &x) // we compare pointers, not references, therefore not recusively
      return false;
    if(size_ != x.size_)
      return true;

    // we assume memcmp is compiler optimized and not based on a function call
    return memcmp(data_, x.data_, size_ * sizeof(T)) == 0 ? false : true;
  }



  // ********************************
  //   Overload Assignment Operator
  // ********************************
  template<typename T>
  SingleVector& Vector<T>::operator=(const SingleVector &stdvec)
  {
    // Down-cast base vector
    const Vector<T>& vec = dynamic_cast<const Vector<T>&>(stdvec);

    // Assign vector entries avoiding self-assignments
    if(this != &vec)
    {
      if(size_ != vec.size_)
        Resize(vec.size_);

      std::copy(vec.data_, vec.data_+size_, data_);
    }

    return *this;
  }
  
  // **********************************************
  //   Cross Product Calculation for vector in 3D
  // **********************************************  
  template <typename T> 
  void Vector<T>::CrossProduct(const Vector<T>& b, Vector<T>& v) const
  {
    if( size_ != 3 && b.size_ != 3 )
      EXCEPTION("CrossProduct can only be calculated for vector of size 3!");

    v.Resize(3);
    v[0] = data_[1] * b.data_[2] - data_[2] * b.data_[1];
    v[1] = data_[2] * b.data_[0] - data_[0] * b.data_[2];
    v[2] = data_[0] * b.data_[1] - data_[1] * b.data_[0];
  }

  template<typename T>
  bool Vector<T>::Collinear( const Vector<T>& vec) {
    Vector<T> a;
    CrossProduct(vec, a);
    return ( fabs(a.NormL2()) < 1e-12 );
  }

// the USE_EMBEDDED_PYTHON implementation is VectorPython.cc
#ifndef USE_EMBEDDED_PYTHON

  template<typename T>
  Vector<T>::Vector(PyObject* obj, bool decref)
  {
    this->capacity_ = 0;
    this->size_ = 0;
    this->memBelongsToMe_= true;
    this->data_= NULL;
    EXCEPTION("Compile with USE_EMBEDDED_PYTHON");
  }

  template<typename T>
  void Vector<T>::Fill(PyObject* obj, bool decref)
  {
    EXCEPTION("Compile with USE_EMBEDDED_PYTHON");
  }

  template<typename T>
  void Vector<T>::Export(PyObject* obj)
  {
    EXCEPTION("Compile with USE_EMBEDDED_PYTHON");
  }
#endif


  
  // ***********************************************************************
  //   Operator implementation for debug case without expression templates
  // ***********************************************************************  
   
#ifndef USE_EXPRESSION_TEMPLATES

  template<typename T>
  Vector<T> &Vector<T>::operator=(const Vector<T> &x)
  {
    // Assign vector entries avoiding self-assignments
    if(this != &x)
    {
      if(size_ != x.size_)
      {
        assert(memBelongsToMe_);
        Resize(x.size_);
      }
      
      std::copy(x.data_, x.data_+size_, data_);
    }
    
    return *this;
  }


  template<typename T>
  Vector<T> &Vector<T>::operator+=(const Vector<T> &x)
  {
#ifdef CHECK_INITIALIZED
    if ((size_ == 0) || (x.size_ == 0))
      EXCEPTION("Vector: undefined Vector in operator +=(vector)" );
#endif  

#ifdef CHECK_INDEX
    if (size_ != x.size_)
      EXCEPTION("Vector: incompatible dimension for operator +=(vector)" );
#endif

    for(unsigned int i = 0; i < size_; ++i)
      data_[i] += x.data_[i];

    return *this;
  }

  template<typename T>
  Vector<T> Vector<T>::operator-() const
  {
#ifdef CHECK_INITIALIZED
    if (size_ == 0)
      EXCEPTION("Vector: undefined Vector in operator-()" );
#endif

    Vector ret(size_);

    for(unsigned int i = 0; i < size_; ++i)
      ret.data_ [i] = -data_ [i];

    return ret;
  }

  template<typename T>
  Vector<T> Vector<T>::operator+() const
  {
#ifdef CHECK_INITIALIZED
    if (size_ == 0)
      EXCEPTION("Vector: undefined Vector in operator +()" );
#endif
    

    // this is not the addition operator but the positive sign, therefore nothing changes!
    return *this;
  }

 

  template<typename T>
  Vector<T> &Vector<T>::operator-=(const Vector<T> &x)
  {
#ifdef CHECK_INITIALIZED
    if ((size_ == 0) || (x.size_ == 0))
      EXCEPTION("Vector: undefined Vector in operator -=(vector)" );
#endif  
  
#ifdef CHECK_INDEX
    if (size_ != x.size_)
      EXCEPTION("Vector: incompatible dimension for operator -=(vector)" );
#endif

    for(unsigned int i = 0; i < size_; ++i)
      data_[i] -= x.data_[i];
  
    return *this;
  }


  template<typename T>
  Vector<T> &Vector<T>::operator/= (T x)
  {
#ifdef CHECK_INITIALIZED
    if (size_ == 0)
      EXCEPTION("Vector: undefined Vector in operator /=(number)" );
#endif
  
    for(unsigned int i = 0; i < size_; ++i)
      data_[i] /= x;
  
    return *this;
  }


  template<typename T>
  Vector<T> &Vector<T>::operator*= (T x)
  {
#ifdef CHECK_INITIALIZED
    if (size_ == 0)
      EXCEPTION("Vector: undefined Vector in operator*=" );
#endif
  
    for(unsigned int i = 0; i < size_; ++i)
      data_[i] *= x;
  
    return *this;
  }
#endif

  // Overloading << for Vector
   template<typename T>
   std::ostream & operator<<(std::ostream &out, const Vector<T> &vc)
   {
     out << vc.ToString();
     return out;
   }

// Explicit template instantiation
  template class Vector<Double>;
  template class Vector<Complex>;
  template class Vector<Integer>;
  template class Vector<UInt>;
  template std::ostream & operator<<<Double> (std::ostream & , const Vector<Double> &);
  template std::ostream & operator<<<Complex> (std::ostream & , const Vector<Complex> & );
  template std::ostream & operator<<<UInt> (std::ostream & , const Vector<UInt> &);
  template std::ostream & operator<<<Integer> (std::ostream & , const Vector<Integer> &);


}
