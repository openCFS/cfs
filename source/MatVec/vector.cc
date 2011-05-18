// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

// The following headers are required for Export()
#include <cstdio>

// OLAS headers
#include "vector.hh"
#include "opdefs.hh"

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


  // *******************************
  //   Constructor from Point class
  // *******************************
  template<typename T>
  Vector<T>::Vector(const Point & p) : 
    SingleVector(), data_(new T[3]), capacity_(3), memBelongsToMe_(true)
  {
    size_ = 3;
    for (unsigned int i = 0; i < 3; ++i)
      data_ [i] = static_cast<T>(p[i]);
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

    for(unsigned int i = 0; i < new_size; i++)
      data_[i] = source[i];
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

#pragma omp parallel for 
    for(unsigned int i = 0; i < size_; ++i)
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

#pragma omp parallel for 
    for(unsigned int i = 0; i < size_; ++i)
      data_[i] = a * idvec1[i] + b * idvec2[i];	
  }
  
  /*
  template <typename T>
  void Vector<T>::Add(Double a, const SingleVector& vec1,
                      Double b, const SingleVector& vec2) {
    if(typeid(Double) != typeid(T))
      EXCEPTION("Wrong typeids in Vector<T>::Add ("
                << typeid(Double).name() << " != "
                << typeid(T).name() << ")");

    //    unsigned int result = Add(a, vec1, b, vec2);
  }

  template <typename T>
  void Vector<T>::Add(Complex a, const SingleVector& vec1,
                      Complex b, const SingleVector& vec2) {
    if(typeid(Double) != typeid(T))
      EXCEPTION("Wrong typeids in Vector<T>::Add ("
                << typeid(Double).name() << " != "
                << typeid(T).name() << ")");
  }
*/
  
  // **********************************************************
  //   Add a scaled version of a vector to this vector object
  // **********************************************************
  template <typename T>
  void Vector<T>::Add(T a, const SingleVector &vec)
  {
    const Vector<T>& idvec = dynamic_cast<const Vector<T>&>(vec);

#pragma omp parallel for 
    for(unsigned int i = 0; i < size_; ++i)
      data_[i] += a * idvec[i];
  }
  
  /*
  template <typename T>
  void Vector<T>::Add(Double a,const SingleVector &vec) {
    if(typeid(Double) != typeid(T))
      EXCEPTION("Wrong typeids in Vector<T>::Add ("
                << typeid(Double).name() << " != "
                << typeid(T).name() << ")");
  }
  
  template <typename T>
  void Vector<T>::Add(Complex a,const SingleVector &vec) {
    if(typeid(Complex) != typeid(T))
      EXCEPTION("Wrong typeids in Vector<T>::Add ("
                << typeid(Complex).name() << " != "
                << typeid(T).name() << ")");
  }
*/

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

  // **********************************
  //   Set all vector entries to zero
  // **********************************
  template<typename T>
  void Vector<T>::Init( T entry)
  {
#pragma omp parallel for 
    for(unsigned int i = 0; i < size_; ++i)
      data_[i] = entry;
  }


  // **************************
  //   Compute Euclidean Norm
  // **************************
  template <typename T>
  double Vector<T>::NormL2() const
  {				
    double sum(0.0);

#pragma omp parallel for reduction(+:sum)
    for(unsigned int i = 0; i < size_; ++i)
      sum += OpType<T>::zConjz(data_[i]);
    
    return sqrt(sum);
  }

  template<class TYPE> 
  Double Vector<TYPE>::NormMax() const 
  { 
    EXCEPTION( "Vector<TPYE>::NormMax only defined for TYPE=Complex/Double" ); 
    return TYPE(); 
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
    for ( unsigned int i = 0; i < size_; ++i )
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
      EXCEPTION( "Vector<" << AssocType<T>::tagV << ">::SetVectorEntry: "
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
    val = data_[i];
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
    }
    
    return ret;
  }
  
  template<typename T>
  void Vector<T>::SetPart( Global::ComplexPart part, const Vector<Double> & partVector )
  {
    EXCEPTION( "Vector::SetPart: Only Implemented for Real and Complex vectors!" );
  }

  template<>
  void Vector<Double>::SetPart( Global::ComplexPart part, const Vector<Double> & partVector )
  {
    if(size_ != partVector.GetSize())
      EXCEPTION( "Vector<Double>::SetPart: Dimension of vectors do not match!" );

    if(part != Global::REAL)
      EXCEPTION( "Vector<Double>::SetPart: Only possible for REAL part." );
    
    *this = partVector;
  }

  template<>
  void Vector<Complex>::SetPart( Global::ComplexPart part, const Vector<Double> & partVector )
  {
    if(size_ != partVector.GetSize())
      EXCEPTION( "Vector<Complex>::SetPart: Dimension of vectors do not match!" );

    switch(part)
    {
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
    }
  }

  template<typename T>
  std::string Vector<T>::ToString(const int level, const char separator) const
  {
    assert(level == 0 || level == 1);
    std::ostringstream os;
    int nnz = 0;
    for(unsigned int i = 0; i < size_; ++i)
    {
      if(level == 1 && Abs(data_[i]) == 0) continue;
      if(level == 1) os << " " << i << ":";
      os << data_[i];
      if(i < size_-1) os << separator;
      nnz++;
    }

    if(level > 0)
      os << " size=" << size_ << " nnz=" << nnz;

    return os.str();
  }
  


  // *****************
  //   Export vector
  // *****************
  template<typename T>
  void Vector<T>::Export(const char *fname) const
  {
    // open output file and check for errors
    FILE *fp = fopen( fname, "w" );
    if(fp == NULL)
      EXCEPTION( "Cannot open file '" << fname << "' for writing!" );

    // Print dimension and entries
    for(unsigned int k = 0; k < size_; ++k)
    {
      OpType<T>::ExportEntry(data_[k], 0, fp);
      fprintf(fp, "\n");
    }

    // close output file
    if(fclose( fp ) == EOF)
    {
      WARN("Could not close file '" << fname << "' after writing!");
    }
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
  bool Vector<TYPE>::ContainsNaN() const
  {
    for(UInt k = 0, s = size_; k < s; ++k)
      if(std::isnan(data_[k])) return true;

    return false;
  }

  template<>
  bool Vector<Complex>::ContainsNaN() const
  {
    for(UInt k = 0, s = size_; k < s; ++k)
    {
      if(std::isnan(data_[k].real())) return true;
      if(std::isnan(data_[k].imag())) return true;
    }
    return false;
  }


  template<class TYPE>
  bool Vector<TYPE>::ContainsInf() const
  {
    for(UInt k = 0, s = size_; k < s; ++k)
      if(std::isinf(data_[k])) return true;

    return false;
  }

  template<>
  bool Vector<Complex>::ContainsInf() const
  {
    for(UInt k = 0, s = size_; k < s; ++k)
    {
      if(std::isinf(data_[k].real())) return true;
      if(std::isinf(data_[k].imag())) return true;
    }
    return false;
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

      for(unsigned int i = 0; i < size_; ++i)
        data_[i] = vec.data_[i];
    }

    return *this;
  }
  
  // **********************************************
  //   Cross Product Calculation for vector in 3D
  // **********************************************  
  template <typename T> 
  void Vector<T>::CrossProduct(const Vector<T>& b, Vector<T>& v)
  {
    if( size_ != 3 || b.size_ != 3 )
      EXCEPTION("CrossProduct can only be calculated for vector of size 3!");

    v.Resize(3);
    v[0] = data_[1] * b.data_[2] - data_[2] * b.data_[1];
    v[1] = data_[2] * b.data_[0] - data_[0] * b.data_[2];
    v[2] = data_[0] * b.data_[1] - data_[1] * b.data_[0];
  }
   
  // ***********************************************************************
  //   Operator implementation for debug case without expression templates
  // ***********************************************************************  
   
#ifndef EXPR_TEMPLATES

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
      
      for(unsigned int i = 0; i < size_; ++i)
        data_[i] = x.data_[i];
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
     out << vc.ToString(0, '\n' );
     return out;
   }


// Explicit template instantiation
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template class Vector<Double>;
  template class Vector<Complex>;
  template class Vector<Integer>;
  template class Vector<unsigned int>;
  template std::ostream & operator<<<Double> (std::ostream & , const Vector<Double> &);
  template std::ostream & operator<<<Complex> (std::ostream & , const Vector<Complex> & );
  template std::ostream & operator<<<unsigned int> (std::ostream & , const Vector<unsigned int> &);
  template std::ostream & operator<<<Integer> (std::ostream & , const Vector<Integer> &);
#endif



}
