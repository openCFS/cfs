// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

// The following headers are required for Export()
#include <cstdio>

// OLAS headers
#include "matvec/vector.hh"

namespace OLAS {

  // **********************************************
  //   Constructor for vector of specified length
  // **********************************************
  template<typename T>
  Vector<T>::Vector( Integer size ) {

    size_ = size;

    // Allocate memory and initialise entries to zero
    dataSize_ = size;
    NewArray( data_, T_Vtype, size_ );
    for ( UInt i = 1; i <= dataSize_; i++ ) {
      data_[i] = 0.0;
    }
    memBelongsToMe_ = true;

    // Set number of unknowns per entry
    dof_ = BlockSize<T>::size;
  }


  // *******************
  //   Deep destructor
  // *******************
  template<typename T>
  Vector<T>::~Vector() {
    if ( memBelongsToMe_ == true ) {
      DeleteArray( data_ );
    }
  }


  // **********************
  //   Replace data array
  // **********************
  template<typename T>
  void Vector<T>::Replace( UInt length, T_Vtype* entries, bool transferMem ) {


    // De-allocate old array, if required
    if ( memBelongsToMe_ == true ) {
      DeleteArray( data_ );
    }

    // Re-set internal attributes
    data_           = entries;
    size_           = length;
    dataSize_       = length;
    memBelongsToMe_ = transferMem;

  }


  // *********************************
  //   Change the size of the vector
  // *********************************
  template<typename T>
  void Vector<T>::Resize( UInt newSize, bool init ) {


    if ( memBelongsToMe_ == false ) {
      (*error) << "I'm cowardly refusing to re-size the data_ array, since "
	       << "the memory does not belong to me!";
      Error( __FILE__, __LINE__ );
    }

    // If new vector is shorter, only adapt size_
    // but do not re-allocate memory
    if ( newSize <= size_ ) {
      size_ = newSize;
    }

    // If new vector is longer, but not longer
    // than the currently allocated maximum
    // only adapt size_
    else if ( newSize <= dataSize_ ) {
      size_ = newSize;
    }

    // New vector is longer than the currently
    // allocated data array, so we re-allocate
    else {
      DeleteArray( data_ );
      NewArray( data_, T_Vtype, newSize );
      size_ = newSize;
      dataSize_ = newSize;
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
  void Vector<T>::Add(const SparseVector &vec) {

    PROFILE("Vector::Add (1)",size_*BlockSize<T>::size);
    TRY_CAST
    ConstRefCast(vec,Vector<T>,idvec);

#pragma omp parallel for 

    for ( unsigned int i = 1; i <= size_; i++ ) {
      data_[i] += idvec[i];
    }
    CATCH_CAST
  }


  // ***************************************************************
  //   Replace this vector object by the sum of two scaled vectors
  // ***************************************************************
  template <typename T>
  void Vector<T>::Add(T_Stype a, const SparseVector &vec1, 
		      T_Stype b, const SparseVector &vec2) {

    PROFILE("Vector::Add (2)",3*size_*BlockSize<T>::size);
    TRY_CAST
    ConstRefCast(vec1,Vector<T>,idvec1);
    ConstRefCast(vec2,Vector<T>,idvec2);

#pragma omp parallel for 

    for ( unsigned int i = 1; i <= size_; i++ ) {
      data_[i] = a * idvec1[i] + b * idvec2[i];	
    }
    CATCH_CAST	
  }


  // **********************************************************
  //   Add a scaled version of a vector to this vector object
  // **********************************************************
  template <typename T>
  void Vector<T>::Add(T_Stype a, const SparseVector &vec) {
    PROFILE("Vector::Add (3)",2*size_*BlockSize<T>::size);

    TRY_CAST
      ConstRefCast(vec,Vector<T>,idvec);

#pragma omp parallel for 

    for ( unsigned int i = 1; i <= size_; i++ ) {
      data_[i] += a * idvec[i];
    }
    CATCH_CAST
  }


  // *********
  //   Inner
  // *********
  template <typename T>
  void Vector<T>::Inner( const SparseVector &vec, T_Stype &sum ) const {

    PROFILE( "Vector::Inner", size_ * 2 * BlockSize<T>::size );

    TRY_CAST {
      ConstRefCast( vec, Vector<T>, secVec );
      sum = 0;

      for ( UInt i = 1; i <= size_; i++ ) {	
        sum += opType<T_Vtype>::dotProduct( data_[i], secVec[i] );
      }
    } CATCH_CAST;
  }


  // **********************************
  //   Set all vector entries to zero
  // **********************************
  template<typename T>
  void Vector<T>::Init() {

#pragma omp parallel for 

    for ( unsigned int i = 1; i <= size_; i++ ) {
      // data_[i] = T_Vtype();
      data_[i] = 0.0;
    }
  }


  // **************************
  //   Compute Euclidean Norm
  // **************************
  template <typename T>
  Double Vector<T>::NormEuclid() const {				
    PROFILE("Vector::NormEuclid",size_*2*BlockSize<T>::size);
    Double sum = 0.0;

#pragma omp parallel for reduction(+:sum)

    for ( unsigned int i = 1; i <= size_; i++ ){
      sum += opType<T_Vtype>::zConjz(data_[i]);
    }
    return sqrt(sum);
  }


  // *******************************************
  //   Same as the BLAS functions of that name
  // *******************************************
  template <typename T>
  void Vector<T>::Axpy( const T_Stype alpha, const SparseVector &y ) {
    PROFILE("Vector::Axpy",size_*2*BlockSize<T>::size);
    TRY_CAST
    ConstRefCast( y, Vector<T>, vec );

#pragma omp parallel for 

    for ( unsigned int i = 1; i <= size_; i++ ) {
      data_[i] = alpha * data_[i] + vec[i];
    }
    CATCH_CAST
  }


  // ******************
  //   SetVectorEntry
  // ******************
  template<typename T>
  void Vector<T>::SetVectorEntry( const Integer i, const T_Vtype &val ) {
#ifdef DEBUG_VECTOR
    if ( i <= 0 || i > size_ ) {
      (*error) << "Vector<" << assocType<T>::tagV << ">::SetVectorEntry: "
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
  void Vector<T>::AddToVectorEntry( const Integer i, const T_Vtype &val ) {
    data_[i] += val;
  }


  // ********************
  //   Get vector entry
  // ********************
  template<typename T>
  void Vector<T>::GetVectorEntry( const Integer i, T_Vtype &val ) const {
    val = data_[i];
  }


  template<typename T>
  std::string Vector<T>::ToString() const
  {
    std::ostringstream os;
    for(unsigned int i = 1; i <= size_; i++)
    {
      os << data_[i];
      if(i < size_) os << ", ";
    }
    
    return os.str();
  }
  
  // ****************
  //   Print vector
  // ****************
  template<typename T>
  void Vector<T>::Print(std::ostream& os) const {
    for ( UInt i = 1; i <= size_; i++ ) {
      (os) << data_[i] << " ";
    }
    (os) << std::endl;
  }


  // *****************
  //   Export vector
  // *****************
  template<typename T>
  void Vector<T>::Export( const Char *fname ) const {

    // open output file and check for errors
    FILE *fp = fopen( fname, "w" );
    if ( fp == NULL ) {
      (*error) << "Cannot open file '" << fname << "' for writing!";
      Error( __FILE__, __LINE__ );
    }

    // Print dimension and entries
    fprintf( fp, "%s%d\n", "%", size_ * dof_ );
    for ( UInt k = 1; k <= size_; k++ ) {
      for ( UInt j = 0; j < dof_; j++ ) {
	opType<T_Vtype>::ExportEntry( data_[k], j, fp );
	fprintf( fp, "\n" );
      }
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

    for ( UInt i = 1; i <= size_; i++ ) {
      data_[i] /= factor;
    }
  }


  // *********************
  //   ScalarMult (real)
  // *********************
  template<typename T>
  void Vector<T>::ScalarMult( const Double factor ) {


    for ( UInt i = 1; i <= size_; i++ ) {
      data_[i] *= factor;
    }
  }


  // ***********************
  //   ScalarDiv (complex)
  // ***********************
  template<typename T>
  void Vector<T>::ScalarDiv( const Complex factor ) {


    for ( UInt i = 1; i <= size_; i++ ) {
      opType<T_Vtype>::DivByComplex( data_[i], factor );
    }
  }


  // ************************
  //   ScalarMult (complex)
  // ************************
  template<typename T>
  void Vector<T>::ScalarMult( const Complex factor ) {


    for ( UInt i = 1; i <= size_; i++ ) {
      opType<T_Vtype>::MultWithComplex( data_[i], factor );
    }
  }


  // ********************************
  //   Overload Assignment Operator
  // ********************************
  template<typename T>
  SparseVector& Vector<T>::operator=( const SparseVector &stdvec ) {

    // Down-cast base vector
#ifdef DEBUG_VECTOR
    try {
#endif

      ConstRefCast( stdvec, Vector<T>, vec );

      if ( size_ != vec.size_ ) {
	Resize(vec.size_);
      }

      // Assign vector entries avoiding self-assignments
      if ( this != &vec ) {
	for ( UInt i = 1; i <= size_; i++ ) {
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


  // ***********************
  //  Forced Instantiation
  // ***********************
  template <typename T>
  void Vector<T>::InstantiatePublicMethods() {


    Error( "This function should never be called", __FILE__, __LINE__ );

    T_Vtype aux;

    Resize( 0, true );
    SetVectorEntry( 0, aux );
  }

}
