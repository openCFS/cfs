// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "matvec/sbmmatrix.hh"


namespace OLAS {


  // ***************
  //   Constructor
  // ***************
  SBM_Matrix::SBM_Matrix( UInt nrows, UInt ncols, bool amSymm ) :
    nrows_(nrows), ncols_(ncols) {


    // Set symmetry flag and check for matrix being square
    amSymm_ = amSymm;
    if ( amSymm_ == true && (ncols_ != nrows_) ) {
      (*error) << "SBM_Matrix::SBM_Matrix: Request for a symmetric "
               << "SBM_Matrix of dimension " << nrows << " x " << ncols
               << "! However, we only support square symmetric matrices!";
      Error( __FILE__, __LINE__ );
    }

    // No entry type or block size known so far
    myEntryType_ = NOENTRYTYPE;
    myBlockSize_ = 0;

    // Generate array for sub-matrices and set pointers to NULL
    NewArray( subMat_, StdMatrix*, nrows_ * ncols_ );
    for( Integer i = 1; i <= nrows_ * ncols_; i++ ) {
      subMat_[i] = NULL;
    }

  }


  // *****************************
  //   Initialise matrix to zero
  // *****************************
  void SBM_Matrix::Init() {
    for( Integer i = 1; i <= nrows_ * ncols_; i++ ) {
      if ( subMat_[i] != NULL ) {
        subMat_[i]->Init();
      }
    }
  }


  // ***********************
  //   Insert a sub-matrix
  // ***********************
  void SBM_Matrix::SetSubMatrix( Integer i, Integer j,
                                 MatrixEntryType eType,
                                 MatrixStorageType sType,
                                 Integer blockSize,
                                 Integer nrows, Integer ncols,
                                 Integer nnz ) {


    // If matrix is symmetric check that the submatrix lies
    // in the upper triangular part
    if ( amSymm_ == true ) {
      if ( j < i ) {
        (*warning) << "SBM_Matrix::SetSubMatrix: Not setting sub-matrix at "
                   << "position ( " << i << " , " << j << " )"
                   << ", since SBM_Matrix is symmetric we only store its "
                   << "triangle!";
        Warning( __FILE__, __LINE__ );
        return;
      }
    }

    // Check that there is no sub-matrix at this position yet
    if ( subMat_[ ComputeIndex(i,j) ] != NULL ) {
      (*error) << "SBM_Matrix::SetSubMatrix: Cowardly refusing the over-write"
               << " sub-matrix at position ( " << i << " , " << j << " )";
      Error( __FILE__, __LINE__ );
    }

    // Check that entry type of sub-matrix is valid, i.e. matches
    // that of the existing sub-matrices
    if ( myEntryType_ == NOENTRYTYPE ) {
      myEntryType_ = eType;
    }
    else if ( eType != myEntryType_ ) {
      (*error) << "SBM_Matrix::SetSubMatrix: Entry type '"
               << Enum2String( eType )
               << "' does not match entry type '"
               << Enum2String( myEntryType_ )
               << "' of other sub-matrices";
      Error( __FILE__, __LINE__ );
    }

    // Check that block size of sub-matrix is valid, i.e. matches
    // that of the existing sub-matrices
    if ( myBlockSize_ == 0 ) {
      myBlockSize_ = blockSize;
    }
    else if ( blockSize != myBlockSize_ ) {
      (*error) << "SBM_Matrix::SetSubMatrix: Block size " << blockSize
               << " does not match block size " << myBlockSize_
               << " of other sub-matrices";
      Error( __FILE__, __LINE__ );
    }

    // Check that the dimensions of the new matrix fit to its neighbours
    // in the super-block matrix
    UInt auxDim;

    // left neighbour
    if ( j > 1 && subMat_[ComputeIndex(i,j-1)] != NULL ) {
      auxDim = subMat_[ComputeIndex(i,j-1)]->GetNrows();
      if ( auxDim != nrows ) {
      (*error) << "SBM_Matrix::SetSubMatrix: "
               << "Cannot insert matrix at (" << i << "," << j << ") with "
               << nrows << " rows, since matrix at ("
               << i << "," << j-1 << ") "
               << "has " << auxDim << " rows!";
      Error( __FILE__, __LINE__ );
      }
    }

    // right neighbour
    if ( j < ncols_ && subMat_[ComputeIndex(i,j+1)] != NULL ) {
      auxDim = subMat_[ComputeIndex(i,j+1)]->GetNrows();
      if ( auxDim != nrows ) {
      (*error) << "SBM_Matrix::SetSubMatrix: "
               << "Cannot insert matrix at (" << i << "," << j << ") with "
               << nrows << "rows, since matrix at ("
               << i << "," << j+1 << ") "
               << "has " << auxDim << " rows!";
      Error( __FILE__, __LINE__ );
      }
    }

    // top neighbour
    if ( i > 1 && subMat_[ComputeIndex(i-1,j)] != NULL ) {
      auxDim = subMat_[ComputeIndex(i-1,j)]->GetNcols();
      if ( auxDim != ncols ) {
      (*error) << "SBM_Matrix::SetSubMatrix: "
               << "Cannot insert matrix at (" << i << " , " << j << ") with "
               << ncols << " columns, since matrix at ("
               << i-1 << "," << j << ") "
               << "has " << auxDim << " column!";
      Error( __FILE__, __LINE__ );
      }
    }

    // bottom neighbour
    if ( i < nrows_ && subMat_[ComputeIndex(i+1,j)] != NULL ) {
      auxDim = subMat_[ComputeIndex(i+1,j)]->GetNcols();
      if ( auxDim != ncols ) {
      (*error) << "SBM_Matrix::SetSubMatrix: "
               << "Cannot insert matrix at (" << i << "," << j << ") with "
               << ncols << "columns, since matrix at ("
               << i+1 << "," << j << ") "
               << "has " << auxDim << " column!";
      Error( __FILE__, __LINE__ );
      }
    }

    // No problems detected, so generate and insert sub-matrix
    StdMatrix *stdMat = GenerateStdMatrixObject( eType, sType, blockSize,
                                                 nrows, ncols, nnz );
    subMat_[ ComputeIndex(i,j) ] = stdMat;
  }


  // *******************************************
  //   Perform a scaled matrix-matrix addition
  // *******************************************
  void SBM_Matrix::Add( const Double fac, const BaseMatrix &mat ) {


    // Downcast BaseMatrix to SBM_Matrix
    try {

      const SBM_Matrix& sbmMat = dynamic_cast<const SBM_Matrix&>(mat);

      // Make sure that both matrices are either symmetric
      // or unsymmetric
      if ( amSymm_ != sbmMat.amSymm_ ) {
        (*error) << "SBM_Matrix::Add: Cannot mix symmetric and "
                 << "unsymmetric matrices in Add() so far!";
        Error( __FILE__, __LINE__ );
      }

      StdMatrix *auxMat = NULL;

      // Loop over all sub-blocks and trigger addition
      for ( UInt k = 1; k <= nrows_ * ncols_; k++ ) {

        // Get hold of submatrix
        auxMat = sbmMat.subMat_[k];

        // Only perform addition, if submatrix is not NULL
        if ( auxMat != NULL ) {

          // Check whether our submatrix is NULL, in this case
          // we copy the submatrix and scale it
          if ( subMat_[k] == NULL ) {
            subMat_[k] = CopyStdMatrixObject( *auxMat );
            subMat_[k]->Scale( fac );
          }

          // Both submatrices are not NULL so do a simple scaled addition
	  subMat_[k]->Add( fac, *auxMat );
	}
      }

    }

    // Test whether downcast was successful
    catch(std::bad_alloc &e) {
      Error( WRONG_CAST_MSG, __FILE__, __LINE__ );
    }

  }


  // **********************************
  //   Perform a residual computation
  // **********************************
  void SBM_Matrix::CompRes( BaseVector &r, const BaseVector &x,
			    const BaseVector& b ) const {


    // Downcast BaseVectors to SBM_Vectors
    try {
      SBM_Vector& sbm_r = dynamic_cast<SBM_Vector&>(r);
      const SBM_Vector& sbm_x = dynamic_cast<const SBM_Vector&>(x);
      const SBM_Vector& sbm_b = dynamic_cast<const SBM_Vector&>(b);

      // Set result vector to zero
      sbm_r.Init();

      // ---------------------
      //  Non-symmetric case:
      // ---------------------

      if ( amSymm_ == false ) {

        // Loop over all block-rows
        for ( Integer i = 1; i <= nrows_; i++ ) {

          // Loop over all blocks in this row
          for ( Integer j = 1; j <= ncols_; j++ ) {

            // Compute i-th component of A*x
            if ( subMat_[ComputeIndex(i,j)] != NULL ) {
              subMat_[ComputeIndex(i,j)]->MultAdd( sbm_x(j), sbm_r(i) );
            }
          }

          // Subtract A*x from b
          sbm_r.Axpy( -1.0, sbm_b );
        }
      }

      // -----------------
      //  Symmetric case:
      // -----------------

      else {

        // Loop over all block-rows
        for ( Integer i = 1; i <= nrows_; i++ ) {

          // Loop over all sub-diagonal blocks in this row
          // Add contribution of sub-diagonal blocks to the i-th component
          // of A * x using S(i,j) * x = S(j,i)^T * x
          for ( Integer j = 1; j < i; j++ ) {
            if ( subMat_[ComputeIndex(j,i)] != NULL ) {
              subMat_[ComputeIndex(j,i)]->MultTAdd( sbm_x(j), sbm_r(i) );
            }
          }

          // Add contribution of diagonal block
          if ( subMat_[ComputeIndex(i,i)] != NULL ) {
            subMat_[ComputeIndex(i,i)]->MultAdd( sbm_x(i), sbm_r(i) );
          }

          // Add contribution of super-diagonal blocks
          for ( Integer j = i+1; j <= ncols_; j++ ) {
            if ( subMat_[ComputeIndex(i,j)] != NULL ) {
              subMat_[ComputeIndex(i,j)]->MultAdd( sbm_x(j), sbm_r(i) );
            }
          }

          // Subtract A*x from b
          sbm_r.Axpy( -1.0, sbm_b );
        }
      }
    }

    // Test whether downcast was successful
    catch(std::bad_alloc &e) {
      Error( WRONG_CAST_MSG, __FILE__, __LINE__ );
    }

  }


  // ******************************************
  //   Perform a matrix-vector multiplication
  // ******************************************
  void SBM_Matrix::Mult( const BaseVector& mvec, BaseVector& rvec ) const {
    rvec.Init();
    MultAdd( mvec, rvec );
  }


  // *********************************************************
  //   Perform a matrix-vector multiplication with transpose
  // *********************************************************
  void SBM_Matrix::MultT( const BaseVector& mvec, BaseVector& rvec ) const {
    rvec.Init();
    this->MultTAdd( mvec, rvec );
  }


  // ********************************************************
  //   Perform a matrix-vector multiplication plus addition
  // ********************************************************
  void SBM_Matrix::MultAdd( const BaseVector& mvec, BaseVector& rvec ) const {


    // Downcast BaseVectors to SBM_Vectors
    try {
      const SBM_Vector& vec_m = dynamic_cast<const SBM_Vector&>(mvec);
      SBM_Vector& vec_r = dynamic_cast<SBM_Vector&>(rvec);

      // ---------------------
      //  Non-symmetric case:
      // ---------------------

      if ( amSymm_ == false ) {

        // perform the multiplication
        for ( Integer i = 1; i <= nrows_; i++ ) {
          for ( Integer j = 1; j <= ncols_; j++ ) {
            if ( subMat_[ComputeIndex(i,j)] != NULL ) {
              subMat_[ComputeIndex(i,j)]->MultAdd( vec_m(j), vec_r(i) );
            }
          }
        }
      }

      // -----------------
      //  Symmetric case:
      // -----------------

      else {

        // Loop over all block-rows
        for ( Integer i = 1; i <= nrows_; i++ ) {

          // Add contribution of sub-diagonal blocks to the i-th component
          // of A * x using S(i,j) * x = S(j,i)^T * x
          for ( Integer j = 1; j < i; j++ ) {
            if ( subMat_[ComputeIndex(j,i)] != NULL ) {
              subMat_[ComputeIndex(j,i)]->MultTAdd( vec_m(j), vec_r(i) );
            }
          }

          // Add contribution of diagonal block
          if ( subMat_[ComputeIndex(i,i)] != NULL ) {
            subMat_[ComputeIndex(i,i)]->MultAdd( vec_m(i), vec_r(i) );
          }

          // Add contribution of super-diagonal blocks
          for ( Integer j = i+1; j <= ncols_; j++ ) {
            if ( subMat_[ComputeIndex(i,j)] != NULL ) {
              subMat_[ComputeIndex(i,j)]->MultAdd( vec_m(j), vec_r(i) );
            }
          }
        }
      }
    }

    // Test whether downcast was successful
    catch(std::bad_alloc &e) {
      Error( WRONG_CAST_MSG, __FILE__, __LINE__ );
    }

  }


  // **********************************************************
  //   Perform a matrix-vector multiplication with transposed
  //   of matrix and add result to vector
  // **********************************************************
  void SBM_Matrix::MultTAdd( const BaseVector& mvec, BaseVector& rvec ) const{


    Warning( "SBM_Matrix::MultTAdd has not yet been tested, my contain bugs");

    // Downcast BaseVectors to SBM_Vectors
    try {
      const SBM_Vector& vec_m = dynamic_cast<const SBM_Vector&>(mvec);
      SBM_Vector& vec_r = dynamic_cast      <SBM_Vector&>(rvec);
	    
      // ---------------------
      //  Non-symmetric case:
      // ---------------------

      if ( amSymm_ == false ) {

        // Work on one column after the other
        for ( Integer j = 1; j <= ncols_; j++ ) {
          for ( Integer i = 1; i <= nrows_; i++ ) {
            if ( subMat_[ComputeIndex(j,i)] != NULL ) {
              subMat_[ComputeIndex(j,i)]->MultTAdd( vec_m(i), vec_r(j) );
            }
          }
        }
      }

      // -----------------
      //  Symmetric case:
      // -----------------

      else {
        this->MultAdd( mvec, rvec );
      }
    }

    // Test whether downcast was successful
    catch(std::bad_alloc &e) {
      Error( WRONG_CAST_MSG, __FILE__, __LINE__ );
    }
  }


  // ***********************************************************
  //   Perform a matrix-vector multiplication plus subtraction
  // ***********************************************************
  void SBM_Matrix::MultSub( const BaseVector& mvec, BaseVector& rvec ) const {


    // Downcast BaseVectors to SBM_Vectors
    try {
      const SBM_Vector& vec_m = dynamic_cast<const SBM_Vector&>(mvec);
      SBM_Vector& vec_r = dynamic_cast<SBM_Vector&>(rvec);

      // ---------------------
      //  Non-symmetric case:
      // ---------------------

      if ( amSymm_ == false ) {

        for ( Integer i = 1; i <= nrows_; i++ ) {
          for ( Integer j = 1; j <= ncols_; j++ ) {
            if ( subMat_[ComputeIndex(i,j)] != NULL ) {

              // Sub-vector of vec_m exists
              if ( vec_m.GetPointer( j ) != NULL ) {

                // Sub-vector of vec_r exists
                if ( vec_r.GetPointer( j ) != NULL ) {
                  subMat_[ComputeIndex(i,j)]->MultSub( vec_m(j), vec_r(i) );
                }

                // Generate error until we can dynamically generate
                // the sub-vector of vec_r here
                else {
                  (*error) << "Pattern mis-match in SBM_Vectors!";
                  Error( __FILE__, __LINE__ );
                }
              }
            }
          }
        }
      }

      // -----------------
      //  Symmetric case:
      // -----------------

      else {
        (*error) << "MultSub not implemented for symmetric SBM_Matrices, "
                 << "since MultTSub is not for normal matrices!";
        Error( __FILE__, __LINE__ );
      }
    }

    // Test whether downcast was successful
    catch(std::bad_alloc &e) {
      Error( WRONG_CAST_MSG, __FILE__, __LINE__ );
    }
  }


  // **********
  //   Export
  // **********
  void SBM_Matrix::Export( const Char *fname, const Char *comment ) const {


    std::stringstream fileName;
    std::string outFile;

    for ( Integer j = 1; j <= ncols_; j++ ) {
      for ( Integer i = 1; i <= nrows_; i++ ) {

        // construct file name
        fileName.str( "" );
        fileName << fname << '_' << i << '_' << j << ".mtx";
        outFile = fileName.str();

        // export sub-matrix
        if ( subMat_[ComputeIndex(i,j)] != NULL ) {
          subMat_[ComputeIndex(i,j)]->Export( outFile.c_str(), comment );
        }
      }
    }
  }


  // **************
  //   GetMaxDiag
  // **************
  Double SBM_Matrix::GetMaxDiag() const {


    Double retVal = 0.0;
    Double auxVal = 0.0;
    StdMatrix *stdMat = NULL;

    for ( Integer i = 1; i <= nrows_; i++ ) {
      stdMat = subMat_[ComputeIndex(i,i)];
      if ( stdMat != NULL ) {
        auxVal = stdMat->GetMaxDiag();
        retVal = auxVal > retVal ? auxVal : retVal;
      }
    }

    return retVal;
  }


  // ******************
  //   GetNrowsScalar
  // ******************
  Integer SBM_Matrix::GetNrowsScalar() const {


    Integer sum = 0;
    StdMatrix *stdMat = NULL;

    for ( UInt i = 1; i <= nrows_; i++ ) {
      stdMat =  subMat_[ ComputeIndex(i,i) ];
      if ( stdMat != NULL ) {
        sum += stdMat->GetNrowsScalar();
      }
      else {
        (*error) << "SBM_Matrix::GetNrowsScalar:\n "
                 << "Current implementation assumes that all diagonal "
                 << "matrices are set/non-zero!";
        Error( __FILE__, __LINE__ );
      }
    }

    return sum;
  }


  // ******************
  //   GetNcolsScalar
  // ******************
  Integer SBM_Matrix::GetNcolsScalar() const {


    Integer sum = 0;
    StdMatrix *stdMat = NULL;

    for ( UInt i = 1; i <= ncols_; i++ ) {
      stdMat =  subMat_[ ComputeIndex(i,i) ];
      if ( stdMat != NULL ) {
        sum += stdMat->GetNcolsScalar();
      }
      else {
        (*error) << "SBM_Matrix::GetNcolsScalar:\n "
                 << "Current implementation assumes that all diagonal "
                 << "matrices are set/non-zero!";
        Error( __FILE__, __LINE__ );
      }
    }

    return sum;
  }

}
