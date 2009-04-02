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
  Vector<T>::Vector( const UInt size, const T entry ) {


    size_ = size;
    capacity_ = 0;
    data_ = new T[size];
    memBelongsToMe_ = true;
    
    for (UInt i = 0; i < size; i++)
      data_ [i] = entry;
     
  }


  // *******************************
  //   Constructor from Point class
  // *******************************
  template<typename T>
  Vector<T>::Vector(const Point & p)
  {
    size_ = 3;
    capacity_ = 3;
    data_ = new T[3];
    memBelongsToMe_ = true;
    data_[0] = (T) p[0];
    data_[1] = (T) p[1];
    data_[2] = (T) p[2];
  }

  // *******************
  //   Deep destructor
  // *******************
  template<typename T>
  Vector<T>::~Vector() {
    if ( memBelongsToMe_ == true ) {
      delete[] data_;
    }
  }


  // **********************
  //   Replace data array
  // **********************
  template<typename T>
  void Vector<T>::Replace( UInt length, T* entries, bool transferMem ) {


    // De-allocate old array, if required
    if ( memBelongsToMe_ == true ) {
      delete[] data_;
    }

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
  void Vector<T>::Resize( UInt newSize, bool init ) {
    if ( memBelongsToMe_ == false &&
         newSize != size_ ) {
         EXCEPTION( "Refusing to re-size the data_ array, since "
             << "the memory does not belong to me!" );
       }
    // If new vector is shorter, only adapt size_
    // but do not re-allocate memory
    if ( newSize <= size_ ) {
      size_ = newSize;
    }

    // If new vector is longer, but not longer
    // than the currently allocated maximum
    // only adapt size_
    else if ( newSize <= capacity_ ) {
      size_ = newSize;
    }

    // New vector is longer than the currently
    // allocated data array, so we re-allocate
    else {
      delete[] data_;
      data_ = new T[newSize];
      size_ = newSize;
      capacity_ = newSize;
    }

    // Zero vector entries if desired
    if ( init == true ) {
      Init();
    }
  }


  // *********************************
  //   Add vec to this vector object
  // *********************************
  template <typename T>
  void Vector<T>::Add(const SingleVector &vec) {

    PROFILE("Vector::Add (1)",size_*BlockSize<T>::size);
    TRY_CAST
    CONSTREFCAST(vec,Vector<T>,idvec);

#pragma omp parallel for 

    for ( unsigned int i = 0; i < size_; i++ ) {
      data_[i] += idvec[i];
    }
    CATCH_CAST
  }


  // ***************************************************************
  //   Replace this vector object by the sum of two scaled vectors
  // ***************************************************************
  template <typename T>
  void Vector<T>::Add(T a, const SingleVector &vec1, 
		      T b, const SingleVector &vec2) {

    PROFILE("Vector::Add (2)",3*size_*BlockSize<T>::size);
    TRY_CAST
    CONSTREFCAST(vec1,Vector<T>,idvec1);
    CONSTREFCAST(vec2,Vector<T>,idvec2);

#pragma omp parallel for 

    for ( unsigned int i = 0; i < size_; i++ ) {
      data_[i] = a * idvec1[i] + b * idvec2[i];	
    }
    CATCH_CAST;
    
    //    return 0;
  }
  /*
  template <typename T>
  void Vector<T>::Add(Double a, const SingleVector& vec1,
                      Double b, const SingleVector& vec2) {
    if(typeid(Double) != typeid(T))
      EXCEPTION("Wrong typeids in Vector<T>::Add ("
                << typeid(Double).name() << " != "
                << typeid(T).name() << ")");

    //    UInt result = Add(a, vec1, b, vec2);
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
  void Vector<T>::Add(T a, const SingleVector &vec) {
    PROFILE("Vector::Add (3)",2*size_*BlockSize<T>::size);

    TRY_CAST
      CONSTREFCAST(vec,Vector<T>,idvec);

#pragma omp parallel for 

    for ( unsigned int i = 0; i < size_; i++ ) {
      data_[i] += a * idvec[i];
    }
    CATCH_CAST;

    //    return 0;
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
  void Vector<T>::Inner( const SingleVector &vec, T &sum ) const {

    PROFILE( "Vector::Inner", size_ * 2 * BlockSize<T>::size );

    TRY_CAST {
      CONSTREFCAST( vec, Vector<T>, secVec );
      sum = 0;

      for ( UInt i = 0; i < size_; i++ ) {	
        sum += OpType<T>::dotProduct( data_[i], secVec[i] );
      }
    } CATCH_CAST;

    //    return 0;
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
  void Vector<T>::Init( T entry) {

#pragma omp parallel for 

    for ( unsigned int i = 0; i < size_; i++ ) {
      data_[i] = entry;
    }
  }


  // **************************
  //   Compute Euclidean Norm
  // **************************
  template <typename T>
  Double Vector<T>::NormL2() const {				
    PROFILE("Vector::NormEuclid",size_*2*BlockSize<T>::size);
    Double sum = 0.0;

#pragma omp parallel for reduction(+:sum)

    for ( unsigned int i = 0; i < size_; i++ ){
      sum += OpType<T>::zConjz(data_[i]);
    }
    return sqrt(sum);
  }

  template<class TYPE> 
  Double Vector<TYPE>::NormMax() const 
  { 
    EXCEPTION( "Vector<TPYE>::NormMax only defined for TYPE=Complex/Double" ); 
    return TYPE(); 
  } 


  template<> 
  Double Vector<Double>::NormMax() const 
  { 
    Double ret = 0.0; 
    if(size_ == 0) EXCEPTION("empty vector"); 

    for(UInt i = 0; i < size_; i++) 
      ret = std::max(std::abs(ret), std::abs(data_[i])); 

    return ret; 
  } 
  template<> 
  Double Vector<Complex>::NormMax() const 
  { 
    Double ret = 0.0; 
    if(size_ == 0) EXCEPTION("empty vector"); 

    for(UInt i = 0; i < size_; i++) 
      ret = std::max(std::abs(ret), std::abs(data_[i].real())); 

    return ret; 
  }  


  // *******************************************
  //   Same as the BLAS functions of that name
  // *******************************************
  template <typename T>
  void Vector<T>::Axpy( const T alpha, const SingleVector &y ) {
    PROFILE("Vector::Axpy",size_*2*BlockSize<T>::size);
    TRY_CAST
    CONSTREFCAST( y, Vector<T>, vec );

#pragma omp parallel for 

    for ( unsigned int i = 0; i < size_; i++ ) {
      data_[i] = alpha * data_[i] + vec[i];
    }
    CATCH_CAST
  }


  // ******************
  //   SetVectorEntry
  // ******************
  template<typename T>
  void Vector<T>::SetEntry( const UInt i, const T &val ) {
#ifdef DEBUG_VECTOR
    if ( i <= 0 || i > size_ ) {
      (*error) << "Vector<" << AssocType<T>::tagV << ">::SetVectorEntry: "
               << "Detected index error:"
               << "\n index = " << i
               << "\n size  = " << size_;
      Error( __FILE__, __LINE__ );
    }
#endif
    data_[i] = val;
  }


  // ********************
  //   AddToVectorEntry
  // ********************
  template<typename T>
  void Vector<T>::AddToEntry( const UInt i, const T &val ) {
    data_[i] += val;
  }


  // ********************
  //   Get vector entry
  // ********************
  template<typename T>
  void Vector<T>::GetEntry( const UInt i, T &val ) const {
    val = data_[i];
  }
  
  template<typename T>
   void  Vector<T>::Push_back (const T & y ) {
    // Check whether capacity is sufficiently large to perform
      // a push-back operation. If not allocate memory according
      // to the following simply scheme: Each time capacity is
      // exceeded allocate twice as much memory as there was before.
      if ( size_ >= capacity_ ) {

        T *help;

        // perform memory allocation
        capacity_ = size_ == 0 ? 1 : 2 * size_;
        help = new T[ capacity_ ];

        // copy old entries into new buffer
        for ( unsigned int i = 0; i < size_; i++ ) {
          help[i] = data_[i];
        }

        // Perform push-back and increase size
        help[size_] = y;
        size_++;

        // delete old buffer and re-set pointer
        delete[] data_;
        data_ = help;

      } else {
        
        // Perform push-back and increase size
        data_[size_] = y;
        size_++;
      }
  }


  template<typename T>
  Vector<Double> Vector<T>::GetPart( Global::ComplexPart part ) const {
    EXCEPTION( "Matrix::GetPart: Only Implemented for Real and "
               << " Complex matrices!" );
    Vector<Double> temp;
    return temp;
  }


  template<>
  Vector<Double> Vector<Double>::GetPart(  Global::ComplexPart part ) const {

    if ( part != Global::REAL ) {
      EXCEPTION("Vector<Double>::GetPart: Only possible for REAL part." );
    }
    return *this;
  }

  template<>
  Vector<Double> Vector<Complex>::GetPart(  Global::ComplexPart part ) const {

    Vector<Double> ret;
    if ( part == Global::REAL ) {
      ret.Resize( size_ );
      for ( UInt i = 0; i < size_; i++ ) {
        ret[i]  = data_[i].real();
      }
    }
    else if ( part == Global::IMAG ) {
      ret.Resize( size_ );
      for ( UInt i = 0; i < size_; i++ ) {
        ret[i]  = data_[i].imag();
      }
    }
    else {
      EXCEPTION("Vector<Complex>::GetPart: Only possible for REAL or IMAG part!" );
    }

    return ret;
  }
  
  template<typename T>
  void Vector<T>::SetPart( Global::ComplexPart part, const Vector<Double> & partVector ) {
    EXCEPTION( "Vector::SetPart:" <<
               "Only Implemented for Real and Complex vectors!" );
  }

  template<>
  void Vector<Double>::SetPart( Global::ComplexPart part, const Vector<Double> & partVector ) {

    if ( size_ != partVector.GetSize() ) {
      EXCEPTION( "Vector<Double>::SetPart: Dimension of vectors do not match!" );
    }

    if ( part != Global::REAL ) {
      EXCEPTION( "Vector<Double>::SetPart: Only possible for REAL part." );
    }
    *this = partVector;
  }

  template<>
  void Vector<Complex>::SetPart( Global::ComplexPart part, const Vector<Double> & partVector ) {

    if ( size_ != partVector.GetSize() ) {
      EXCEPTION( "Vector<Complex>::SetPart: Dimension of vectors do not match!" );
    }

    if ( part == Global::REAL ) {
      for ( UInt i = 0; i < size_; i++ ) {
        data_[i]  = Complex( partVector[i], data_[i].imag() );
      }
    }
    else if ( part == Global::IMAG ) {
      for ( UInt i = 0; i < size_; i++ ) {
        data_[i]  = Complex( data_[i].real(), partVector[i] );
      }
    }
    else {
      EXCEPTION( "Vector<Complex>::SetPart: Only possible for REAL or IMAG part!" );
    }
  }

  template<typename T>
  std::string Vector<T>::ToString( char separator ) const
  {
    std::ostringstream os;
    
    for(unsigned int i = 0; i < size_; i++)
    {
      os << data_[i];
      if(i < size_) os << separator;
    }

    return os.str();
  }
  


  // *****************
  //   Export vector
  // *****************
  template<typename T>
  void Vector<T>::Export( const char *fname ) const {

    // open output file and check for errors
    FILE *fp = fopen( fname, "w" );
    if ( fp == NULL ) {
      EXCEPTION( "Cannot open file '" << fname << "' for writing!" );
    }

    // Print dimension and entries
    for ( UInt k = 0; k < size_; k++ ) {
      OpType<T>::ExportEntry( data_[k], 0, fp );
      fprintf( fp, "\n" );
    }

    // close output file
    if ( fclose( fp ) == EOF ) {
      (*warning) << "Could not close file '" << fname << "' after writing!";
      Warning( __FILE__, __LINE__ );
    }
  }


  // ********************
  //   ScalarDiv (real)
  // ********************
  template<typename T>
  void Vector<T>::ScalarDiv( const Double factor ) {


#ifdef DEBUG_VECTOR
    if ( factor == 0 ) {
      (*error) << "Vector::ScalarDiv: Division by Zero!";
      Error( __FILE__, __LINE__ );
    }
#endif

    for ( UInt i = 0; i < size_; i++ ) {
      data_[i] = static_cast<T>( data_[i] / factor );
    }
  }


  // *********************
  //   ScalarMult (real)
  // *********************
  template<typename T>
  void Vector<T>::ScalarMult( const Double factor ) {


    for ( UInt i = 0; i < size_; i++ ) {
      data_[i] = static_cast<T>( data_[i] * factor );
    }
  }


  // ***********************
  //   ScalarDiv (complex)
  // ***********************
  template<typename T>
  void Vector<T>::ScalarDiv( const Complex factor ) {


    for ( UInt i = 0; i < size_; i++ ) {
      OpType<T>::DivByComplex( data_[i], factor );
    }
  }


  // ************************
  //   ScalarMult (complex)
  // ************************
  template<typename T>
  void Vector<T>::ScalarMult( const Complex factor ) {


    for ( UInt i = 0; i < size_; i++ ) {
      OpType<T>::MultWithComplex( data_[i], factor );
    }
  }


  // ********************************
  //   Overload Assignment Operator
  // ********************************
  template<typename T>
  SingleVector& Vector<T>::operator=( const SingleVector &stdvec ) {

    // Down-cast base vector
#ifdef DEBUG_VECTOR
    try {
#endif

      CONSTREFCAST( stdvec, Vector<T>, vec );

      if ( size_ != vec.size_ ) {
        Resize(vec.size_);
      }

      // Assign vector entries avoiding self-assignments
      if ( this != &vec ) {
        for ( UInt i = 0; i < size_; i++ ) {
          data_[i] = vec.data_[i];
        }
      }

#ifdef DEBUG_VECTOR
    } catch( std::bad_cast e ) {
      (*error) << "Vector::operator=\n"
      << " could not down-cast right hand side vector to Vector<T> "
      << "with\n"
      << " tagM = " << this->GetTagM()
      << ", tagV = " << this->GetTagV()
      << ", tagS = " << this->GetTagS();
      Error( __FILE__, __LINE__ );
    }
#endif

    return *this;
  }
  
  // **********************************************
  //   Cross Product Calculation for vector in 3D
  // **********************************************  
  template <typename T> 
   void Vector<T>::CrossProduct( const Vector<T>& b, 
                                 Vector<T>& v ) {
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
  
    if (this == &x)
      return *this;
  
    if (size_ != x.size_)
      {   
        if (memBelongsToMe_ == false ) {
          EXCEPTION( "Refusing to resize vector, since memory does not " 
                     << "belong to me!" );
        }
        if (data_)
          delete [] data_;
      
        size_ = x.size_;
        data_ = new T [size_];
      }
  
    for (UInt i = 0; i < size_; i++)
      data_ [i] = x.data_[i];
  
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

    for (UInt i = 0; i < size_; i++)
      data_[i] += x.data_[i];
  
    return *this; 
  }

  template<typename T>
  Vector<T> Vector<T>::operator- () const
  {
#ifdef CHECK_INITIALIZED
    if (size_ == 0)
      EXCEPTION("Vector: undefined Vector in oprator -()" );
#endif

    Vector ret(size_);
  
    for (UInt i = 0; i < size_; i++)
      ret.data_ [i] = -data_ [i];
  
    return ret;
  }

  template<typename T>
  Vector<T> Vector<T>::operator+ () const
  {
#ifdef CHECK_INITIALIZED
    if (size_ == 0)
      EXCEPTION("Vector: undefined Vector in oprator +()" );
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

    for (UInt i = 0; i < size_; i++)
      data_[i] -= x.data_[i];
  
    return *this;
  }


 

  template<typename T>
  Vector<T> &Vector<T>::operator/= (const T &x)
  {
#ifdef CHECK_INITIALIZED
    if (size_ == 0)
      EXCEPTION("Vector: undefined Vector in operator /=(number)" );
#endif

    T y = x;
  
    for (UInt i = 0; i < size_; i++)
      data_[i] /= y;
  
    return *this;
  }


  template<typename T>
  Vector<T> &Vector<T>::operator*= (const T &x)
  {
#ifdef CHECK_INITIALIZED
    if (size_ == 0)
      EXCEPTION("Vector: undefined Vector in operator*=" );
#endif
  
    T y = x;
  
    for (UInt i = 0; i < size_; i++)
      data_[i] *= y;
  
    return *this;
  }
#endif

  // Overloading << for Vector
   template<typename T>
   std::ostream & operator << ( std::ostream & out, const Vector<T> & vc)
   {
     out << vc.ToString( '\n' );
     return out;
   }


// Explicit template instantiation
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template class Vector<Double>;
  template class Vector<Complex>;
  template class Vector<Integer>;
  template class Vector<UInt>;
  template std::ostream & operator<<<Double> (std::ostream & , const Vector<Double> &);
  template std::ostream & operator<<<Complex> (std::ostream & , const Vector<Complex> & );
  template std::ostream & operator<<<UInt> (std::ostream & , const Vector<UInt> &);
  template std::ostream & operator<<<Integer> (std::ostream & , const Vector<Integer> &);
#endif



}
