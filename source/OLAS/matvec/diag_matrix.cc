// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

// the following headers are required for Export()
#include <cstdio>


#include "matvec/diag_matrix.hh"
#include "multigrid/prematrix.hh"


// Implementation of methods for a diagonal matrix

namespace OLAS {

  // *******************************
  //   Add value to a matrix entry
  // *******************************
  template<typename T>
  inline void Diag_Matrix<T>::AddToMatrixEntry( Integer i, Integer j, T& v ) {


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
     TRY_CAST {

      // Down-cast input matrix
      ConstRefCast( mat, Diag_Matrix<T>, diagMat );

      // Obtain pointer to data array of other matrix
      const T *data = diagMat.GetDataPointer();
  
      // We now assume that the matrix have matching
      // dimensions and sparsity patterns
      for ( Integer i = 1; i <= this->nrows_; i++ ) {
	data_[i] += alpha * data[i];
      }

    } CATCH_CAST;
  }


  // ********************************
  //   Matrix-Vector multiplication
  // ********************************
  template<typename T>
  inline void Diag_Matrix<T>::Mult( const Vector<T> &mvec,
                                   Vector<T> &rvec ) const {


    PROFILE( (char*)"Diag_Matrix::Mult", 2*this->nnz_*BlockSize<T>::size *
             BlockSize<T>::size );

    Integer i;

#pragma omp parallel for 
    for ( i = 1; i <= this->nrows_; i++ ) {
      rvec[i] =  data_[i] * mvec[i];
    }
  }


  // ***********
  //   MultAdd
  // ***********
  template<typename T>
  inline void Diag_Matrix<T>::MultAdd( const Vector<T> & mvec,
                                      Vector<T> & rvec ) const {


    Integer i;

#pragma omp parallel for
    for ( i = 1; i <= this->nrows_; i++ ) {
      rvec[i] += data_[i] * mvec[i];
    }
  }


  // ***********
  //   MultSub
  // ***********
  template<typename T>
  inline void Diag_Matrix<T>::MultSub( const Vector<T> &mvec,
                                      Vector<T> &rvec ) const {


    Integer i;

#pragma omp parallel for
    for ( i = 1; i <= this->nrows_; i++ ) {
      rvec[i] -= data_[i] * mvec[i];
    }
  }


  // ***********
  //   CompRes
  // ***********
  template<typename T>
  inline void Diag_Matrix<T>::CompRes( Vector<T> &r, const Vector<T> &x,
                                      const Vector<T> &b ) const {


    PROFILE((char*)"Diag_Matrix::CompRes",
            (2*this->nnz_*BlockSize<T>::size + this->nrows_)
            *BlockSize<T>::size);

    Integer i;

#pragma omp parallel for
    for ( i = 1; i <= this->nrows_; i++ ) {
      r[i] = b[i] - data_[i] * x[i];
    }
  }


  // **********
  //   MultT
  // **********
  template<typename T>
  inline void Diag_Matrix<T>::MultT( const Vector<T> & mvec,
                                    Vector<T> & rvec ) const {


    Integer i;

    // loop over matrix rows
    for ( i = 1; i <= this->nrows_; i++ ) {
      rvec[i] = data_[i] * mvec[i];
    }
  }


  // ************
  //   MultTAdd
  // ************
  template<typename T>
  inline void Diag_Matrix<T>::MultTAdd( const Vector<T> &mvec,
                                       Vector<T> &rvec ) const {

    Integer i;

    // loop over matrix rows
    for( i = 1; i <= this->nrows_; i++ ) {
      rvec[i] += data_[i] * mvec[i];
    }
  }


  // ****************************
  //   Print matrix to a stream
  // ****************************
  template<typename T>
  void Diag_Matrix<T>::Print( std::ostream& os ) const {


    Integer i; // COMPWARNING: unused variables , j, k;

    for ( i = 1; i <= this->nrows_; i++ ) {
      os << std::scientific << i << "\t" << data_[i] << '\n';
    }
    
    os << std::endl;
  }


  // ***********************
  //  Forced Instantiation
  // ***********************
  template <typename T>
  void Diag_Matrix<T>::InstantiatePublicMethods() {


    Error( "This function should never be called", __FILE__, __LINE__ );
  }


  // *********
  //   Scale
  // *********
  template<typename T>
  void Diag_Matrix<T>::Scale( Double factor ) {
    for ( Integer i = 1; i <= this->nnz_; i++ ) {
      data_[i] *= factor;
    }
  }


  // **************
  //   GetMaxDiag
  // **************
  template<typename T>
  Double Diag_Matrix<T>::GetMaxDiag() const {


    double maxDiag = 0;
    double current = 0;
    Integer i;

    for ( i = 1; i <= this->nrows_; i++ ) {

      // use an opType to ensure that tiny matrices
      // are treated correctly
      current = opType<T>::MaxDiag( data_[i] );
      maxDiag = maxDiag > current ? maxDiag : current;
    }

    return maxDiag;
  }


  // ******************
  //   Re-size matrix 
  // ******************
  template<typename T>
  void Diag_Matrix<T>::SetSize( Integer nrows, Integer ncols, Integer nnz ) {

    
#ifdef DEBUG_Diag_MATRIX
    if ( nrows < 0 || nnz < 0 ) {
      Error( "invalid dimensions", __FILE__, __LINE__ );
    }
#endif

    if ( nrows != ncols  || nrows != nnz ) {
      Error("In a pure diagonal matrix nnz = nrows = ncols!",__FILE__,__LINE__);
    }

    this->ncols_ = ncols;
    if ( this->nrows_ != nrows ) {
      this->nrows_ = nrows;
    }

    if ( this->nnz_ != nnz ) {
      this->nnz_ = nnz;
      DeleteArray( data_ );
      NewArray( data_, T, this->nnz_ );
    }
  }

}
