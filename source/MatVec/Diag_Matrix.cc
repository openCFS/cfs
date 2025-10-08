// the following headers are required for Export()
#include <cstdio>


#include "Diag_Matrix.hh"

#include "opdefs.hh"


// Implementation of methods for a diagonal matrix

namespace CoupledField {

  // *******************************
  //   Add value to a matrix entry
  // *******************************
  template<typename T>
  inline void Diag_Matrix<T>::AddToMatrixEntry( UInt i, UInt j, const T& v ) {


    //we just add the diagonal entries
    // Otherwise, we add the value to the existing entry

    if ( i == j ) {
      data_[i] += v;
    }
  }


  // ************************
  //   Add (another matrix)
  // ************************
  template<typename T>
  inline void Diag_Matrix<T>::Add( const Double alpha, const StdMatrix& mat ) {
    // Obtain pointer to data array of other matrix
    const T *data = dynamic_cast<const Diag_Matrix<T>&>(mat).GetDataPointer();

    // We now assume that the matrix have matching
    // dimensions and sparsity patterns
    for ( UInt i = 0; i < this->nrows_; i++ ) {
      data_[i] += alpha * data[i];
    }
  }

  template<typename T>
  inline void Diag_Matrix<T>::Add( const Complex alpha, const StdMatrix& mat ) {
    EXCEPTION("Complex Add is not implemented for Diag_Matrix.");
  }
  
  // ******************************************
  //   Add (another matrix, only index subset)
  // ******************************************
  template<typename T>
  void Diag_Matrix<T>::Add( const Double alpha, const StdMatrix& mat,
                            const std::set<UInt>& rowIndices,
                            const std::set<UInt>& colIndices ) {
    EXCEPTION("Not implemented");
  }

  template<typename T>
  void Diag_Matrix<T>::Add( const Complex alpha, const StdMatrix& mat, const std::set<UInt>& rowIndices, const std::set<UInt>& colIndices ) {
    EXCEPTION("Complex Add is not implemented for Diag_Matrix.");
  }


  // ********************************
  //   Matrix-Vector multiplication
  // ********************************
  template<typename T>
  inline void Diag_Matrix<T>::Mult( const Vector<T> &mvec,
                                   Vector<T> &rvec ) const {

    Integer i = 0;

#pragma omp parallel for 
    for ( i = 0; i < (Integer) this->nrows_; i++ ) {
      rvec[i] =  data_[i] * mvec[i];
    }
  }


  // ***********
  //   MultAdd
  // ***********
  template<typename T>
  inline void Diag_Matrix<T>::MultAdd( const Vector<T> & mvec,
                                      Vector<T> & rvec ) const {


    Integer i = 0;

#pragma omp parallel for
    for ( i = 0; i < (Integer) this->nrows_; i++ ) {
      rvec[i] += data_[i] * mvec[i];
    }
  }


  // ***********
  //   MultSub
  // ***********
  template<typename T>
  inline void Diag_Matrix<T>::MultSub( const Vector<T> &mvec,
                                      Vector<T> &rvec ) const {


    Integer i = 0;

#pragma omp parallel for
    for ( i = 0; i < (Integer) this->nrows_; i++ ) {
      rvec[i] -= data_[i] * mvec[i];
    }
  }


  // ***********
  //   CompRes
  // ***********
  template<typename T>
  inline void Diag_Matrix<T>::CompRes( Vector<T> &r, const Vector<T> &x,
                                      const Vector<T> &b ) const {

    Integer i = 0;

#pragma omp parallel for
    for ( i = 0; i < (Integer) this->nrows_; i++ ) {
      r[i] = b[i] - data_[i] * x[i];
    }
  }


  // **********
  //   MultT
  // **********
  template<typename T>
  inline void Diag_Matrix<T>::MultT( const Vector<T> & mvec,
                                    Vector<T> & rvec ) const {


    UInt i;

    // loop over matrix rows
    for ( i = 0; i < this->nrows_; i++ ) {
      rvec[i] = data_[i] * mvec[i];
    }
  }


  // ************
  //   MultTAdd
  // ************
  template<typename T>
  inline void Diag_Matrix<T>::MultTAdd( const Vector<T> &mvec,
                                       Vector<T> &rvec ) const {

    UInt i;

    // loop over matrix rows
    for( i = 0; i < this->nrows_; i++ ) {
      rvec[i] += data_[i] * mvec[i];
    }
  }
  
  // ************
   //   MultTSub
   // ************
   template<typename T>
   inline void Diag_Matrix<T>::MultTSub( const Vector<T> &mvec,
                                        Vector<T> &rvec ) const {

     UInt i;

     // loop over matrix rows
     for( i = 0; i < this->nrows_; i++ ) {
       rvec[i] -= data_[i] * mvec[i];
     }
   }


  // ****************************
  //   Print matrix to a stream
  // ****************************
  template<typename T>
  std::string Diag_Matrix<T>::ToString( char colSep,
                                        char rowSep )  const {


    UInt i;
    std::stringstream os;

    for ( i = 0; i < this->nrows_; i++ ) {
      os << std::scientific << i << colSep
          << data_[i] << rowSep;
    }
    
    os << std::endl;
    return os.str();
  }


  // *********
  //   Scale
  // *********
  template<typename T>
  void Diag_Matrix<T>::Scale( Double factor ) {
    for ( UInt i = 0; i < this->nnz_; i++ ) {
      data_[i] *= factor;
    }
  }
  
  // ************************
  //   Scale on index subset
  // ************************
  template<typename T>
  void Diag_Matrix<T>::Scale( Double factor,
                              const std::set<UInt>& rowIndices,
                              const std::set<UInt>& colIndices ) {
    EXCEPTION("Implement me");
  }

  // ******************
  //   Re-size matrix 
  // ******************
  template<typename T>
  void Diag_Matrix<T>::SetSize( UInt nrows, UInt ncols, UInt nnz ) {

    
#ifdef DEBUG_Diag_MATRIX
    if ( nrows < 0 || nnz < 0 ) {
      EXCEPTION( "invalid dimensions" );
    }
#endif

    if ( nrows != ncols  || nrows != nnz ) {
      EXCEPTION("In a pure diagonal matrix nnz = nrows = ncols!");
    }

    this->ncols_ = ncols;
    if ( this->nrows_ != nrows ) {
      this->nrows_ = nrows;
    }

    if ( this->nnz_ != nnz ) {
      this->nnz_ = nnz;
      delete [] ( data_ );  data_  = NULL;
      NEWARRAY( data_, T, this->nnz_ );
    }
  }
  
  // ******************
  //   GetMemoryUsage
  // ******************
  template <typename T>
  Double Diag_Matrix<T>::GetMemoryUsage() const {
    return this->nnz_ * sizeof(T);
  }

// Explicit template instantiation
  template class Diag_Matrix<Double>;
  template class Diag_Matrix<Complex>;
  
}
