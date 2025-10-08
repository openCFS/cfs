#include "SBM_Matrix.hh"


namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
  SBM_Matrix::SBM_Matrix( UInt nrows, UInt ncols, bool amSymm ) :
    nrows_(nrows), ncols_(ncols) {
    ownSubMatrices_ = true;


    // Set symmetry flag and check for matrix being square
    amSymm_ = amSymm;
    if ( amSymm_ == true && (ncols_ != nrows_) ) {
      EXCEPTION( "SBM_Matrix::SBM_Matrix: Request for a symmetric "
               << "SBM_Matrix of dimension " << nrows << " x " << ncols
               << "! However, we only support square symmetric matrices!" );
    }

    // No entry type or block size known so far
    myEntryType_ = NOENTRYTYPE;

    // Generate array for sub-matrices and set pointers to NULL
    subMat_.Resize( nrows_ * ncols_ );
    subMat_.Init( NULL );

  }

  // ********************
  //   Copy Constructor
  // ********************
  SBM_Matrix::SBM_Matrix( const SBM_Matrix& origMat ) {

    //  initialize from origMat
    nrows_ = origMat.nrows_;
    ncols_ = origMat.ncols_;
    amSymm_ = origMat.amSymm_;
    ownSubMatrices_ = true;
    
    myEntryType_  = origMat.myEntryType_;
    subMat_.Resize( nrows_ * ncols_ );
    subMat_.Init( NULL );
    for ( UInt k = 0; k < subMat_.GetSize(); k++ ) {
      if(origMat.subMat_[k] != NULL)
        subMat_[k] = CopyStdMatrixObject( *(origMat.subMat_[k]));
    }
  }
  
  // *************************
  //   Weak Copy Constructor
  // *************************
  SBM_Matrix::SBM_Matrix( const SBM_Matrix& origMat,
                          UInt numRows, UInt numCols ) {

    //  initialize from origMat
    nrows_ = numRows;
    ncols_ = numCols;
    amSymm_ = origMat.amSymm_;
    // This is a weak copy, so no ownership of the pointers
    ownSubMatrices_ = false;
    
    myEntryType_  = origMat.myEntryType_;
    subMat_.Resize( nrows_ *  ncols_ );
    subMat_.Init( NULL );
    for( UInt iRow = 0; iRow < nrows_; ++iRow ) {
      for( UInt iCol = 0; iCol < nrows_; ++iCol ) {
      subMat_[ComputeIndex(iRow,iCol)] = 
          origMat.subMat_[origMat.ComputeIndex(iRow,iCol)];
      }
    }
  }
  

  // *****************************
  //   Initialise matrix to zero
  // *****************************
  void SBM_Matrix::Init() {
    for( UInt i = 0; i < subMat_.GetSize(); i++ ) {
      if ( subMat_[i] != NULL ) {
        subMat_[i]->Init();
      }
    }
  }


  // ***********************
  //   Insert a sub-matrix
  // ***********************
  void SBM_Matrix::SetSubMatrix( UInt i, UInt j,
                                 BaseMatrix::EntryType eType,
                                 BaseMatrix::StorageType sType,
                                 UInt nrows, UInt ncols,
                                 UInt nnz, bool forceOverwrite, 
                                 bool silentOverWrite ) {
    
    // Only set submatrix, if this is not a cheap copy
    if( !ownSubMatrices_ ) {
      EXCEPTION( "As this is a weak copy of a SBM-matrix, I refuse to "
                 "set a submatrix" );
    }


    // If matrix is symmetric check that the submatrix lies
    // in the upper triangular part
    if ( amSymm_ == true ) {
      if ( j < i ) {
        WARN("SBM_Matrix::SetSubMatrix: Not setting sub-matrix at "
             << "position ( " << i << " , " << j << " )"
             << ", since SBM_Matrix is symmetric we only store its "
             << "triangle!");
        return;
      }
    }

    // Check that there is no sub-matrix at this position yet
    if ( subMat_[ ComputeIndex(i,j) ] != NULL ) {
      if ( forceOverwrite ) {
        if( ownSubMatrices_ ) {
          delete subMat_[ ComputeIndex(i,j) ];
          if ( !silentOverWrite ) {
            WARN("SBM_Matrix got overwritten");
          }
        } else {
          EXCEPTION("I do not own this matrix and therefore can't overwrite it!");
        }
      } else {
        EXCEPTION( "SBM_Matrix::SetSubMatrix: Cowardly refusing the over-write"
                 << " sub-matrix at position ( " << i << " , " << j << " )");
      }
    }

    // Check that entry type of sub-matrix is valid, i.e. matches
    // that of the existing sub-matrices
    if ( myEntryType_ == NOENTRYTYPE ) {
      myEntryType_ = eType;
    }
    else if ( eType != myEntryType_ ) {
      EXCEPTION( "SBM_Matrix::SetSubMatrix: Entry type '"
               << BaseMatrix::entryType.ToString( eType )
               << "' does not match entry type '"
               << BaseMatrix::entryType.ToString( myEntryType_ )
               << "' of other sub-matrices" );
    }

    // Check that the dimensions of the new matrix fit to its neighbours
    // in the super-block matrix
    UInt auxDim;

    // left neighbour
    if ( j > 0 && subMat_[ComputeIndex(i,j-1)] != NULL ) {
      auxDim = subMat_[ComputeIndex(i,j-1)]->GetNumRows();
      if ( auxDim != nrows ) {
        EXCEPTION( "SBM_Matrix::SetSubMatrix: "
               << "Cannot insert matrix at (" << i << "," << j << ") with "
               << nrows << " rows, since matrix at ("
               << i << "," << j-1 << ") "
               << "has " << auxDim << " rows!" );
      }
    }

    // right neighbour
    if ( j < ncols_-1 && subMat_[ComputeIndex(i,j+1)] != NULL ) {
      auxDim = subMat_[ComputeIndex(i,j+1)]->GetNumRows();
      if ( auxDim != nrows ) {
        EXCEPTION( "SBM_Matrix::SetSubMatrix: "
               << "Cannot insert matrix at (" << i << "," << j << ") with "
               << nrows << "rows, since matrix at ("
               << i << "," << j+1 << ") "
               << "has " << auxDim << " rows!");
      }
    }

    // top neighbour
    if ( i > 1 && subMat_[ComputeIndex(i-1,j)] != NULL ) {
      auxDim = subMat_[ComputeIndex(i-1,j)]->GetNumCols();
      if ( auxDim != ncols ) {
        EXCEPTION( "SBM_Matrix::SetSubMatrix: "
               << "Cannot insert matrix at (" << i << " , " << j << ") with "
               << ncols << " columns, since matrix at ("
               << i-1 << "," << j << ") "
               << "has " << auxDim << " column!" );
      }
    }

    // bottom neighbour
    if ( i < nrows_-1 && subMat_[ComputeIndex(i+1,j)] != NULL ) {
      auxDim = subMat_[ComputeIndex(i+1,j)]->GetNumCols();
      if ( auxDim != ncols ) {
        EXCEPTION( "SBM_Matrix::SetSubMatrix: "
               << "Cannot insert matrix at (" << i << "," << j << ") with "
               << ncols << "columns, since matrix at ("
               << i+1 << "," << j << ") "
               << "has " << auxDim << " column!");
      }
    }

    // No problems detected, so generate and insert sub-matrix
    StdMatrix *stdMat = GenerateStdMatrixObject( eType, sType, 
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
        EXCEPTION( "SBM_Matrix::Add: Cannot mix symmetric and "
                 << "unsymmetric matrices in Add() so far!" );
      }

      StdMatrix *auxMat = NULL;

      // Loop over all sub-blocks and trigger addition
      for ( UInt k = 0; k < nrows_ * ncols_; k++ ) {

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
      EXCEPTION( WRONG_CAST_MSG );
    }

  }
  
  // ******************************************************
  //   Perform a scaled matrix-matrix addition on a subset
  // ******************************************************
  void SBM_Matrix::Add( const Double fac, const BaseMatrix& mat,
                        std::map<UInt, std::set<UInt> >& rowIndPerBlock,
                        std::map<UInt, std::set<UInt> >& colIndPerBlock ) {
    
    // Downcast BaseMatrix to SBM_Matrix
    try {

      const SBM_Matrix& sbmMat = dynamic_cast<const SBM_Matrix&>(mat);

      // Make sure that both matrices are either symmetric
      // or unsymmetric
      if ( amSymm_ != sbmMat.amSymm_ ) {
        EXCEPTION( "SBM_Matrix::Add: Cannot mix symmetric and "
            << "unsymmetric matrices in Add() so far!" );
      }

      StdMatrix *auxMat = NULL;
      
      // fill the row / col sets. If a set is empty,
      // all rows / cols are taken
      std::set<UInt> rowBlocks, colBlocks;
      if( rowIndPerBlock.size() == 0 )  {
        for( UInt i = 0; i < this->nrows_; ++i ) 
          rowBlocks.insert(i);
      } else {
        std::map<UInt, std::set<UInt> >::iterator rowIt;
        rowIt = rowIndPerBlock.begin();
        for( ; rowIt != rowIndPerBlock.end(); rowIt++ ) 
          rowBlocks.insert( rowIt->first);
      }
      if( colIndPerBlock.size() == 0 )  {
        for( UInt i = 0; i < this->ncols_; ++i ) 
          colBlocks.insert(i);
      } else {
        std::map<UInt, std::set<UInt> >::iterator colIt;
        colIt = colIndPerBlock.begin();
        for( ; colIt != colIndPerBlock.end(); colIt++ ) 
          colBlocks.insert( colIt->first);
      }
      
      // Iterate over all SBM rows / cols in given index sets
      std::set<UInt>::iterator rowIt, colIt;
      
      for( rowIt = rowBlocks.begin(); 
          rowIt != rowBlocks.end(); ++rowIt ) {
        
        for( colIt = colBlocks.begin(); 
            colIt != colBlocks.end(); ++colIt ) {
          
          const UInt rowInd = *rowIt;
          const UInt colInd = *colIt;
          UInt index = ComputeIndex( rowInd, colInd);
          
          // Get hold of submatrix
          auxMat = sbmMat.subMat_[index];
          
          // Only perform addition, if submatrix is not NULL
          if ( auxMat != NULL ) {

            
            // Check whether our submatrix is NULL, in this case
            // we copy the submatrix and scale it
            if ( subMat_[index] == NULL ) {
              subMat_[index] = CopyStdMatrixObject( *auxMat );
              subMat_[index]->Scale( fac,
                                     rowIndPerBlock[rowInd],
                                     colIndPerBlock[colInd] );
            }

            // Both submatrices are not NULL so do a simple scaled addition
            subMat_[index]->Add( fac, *auxMat,
                                 rowIndPerBlock[rowInd],
                                 colIndPerBlock[colInd] );
          }
        }
      }
    }
    // Test whether downcast was successful
    catch(std::bad_alloc &e) {
      EXCEPTION( WRONG_CAST_MSG );
    }
  }

  // Overloaded function for complex-valued factor
  void SBM_Matrix::Add( const Complex fac, const BaseMatrix& mat,
                        std::map<UInt, std::set<UInt> >& rowIndPerBlock,
                        std::map<UInt, std::set<UInt> >& colIndPerBlock ) {
    
    // Downcast BaseMatrix to SBM_Matrix
    try {

      const SBM_Matrix& sbmMat = dynamic_cast<const SBM_Matrix&>(mat);

      // Make sure that both matrices are either symmetric
      // or unsymmetric
      if ( amSymm_ != sbmMat.amSymm_ ) {
        EXCEPTION( "SBM_Matrix::Add: Cannot mix symmetric and "
            << "unsymmetric matrices in Add() so far!" );
      }

      StdMatrix *auxMat = NULL;
      
      // fill the row / col sets. If a set is empty,
      // all rows / cols are taken
      std::set<UInt> rowBlocks, colBlocks;
      if( rowIndPerBlock.size() == 0 )  {
        for( UInt i = 0; i < this->nrows_; ++i ) 
          rowBlocks.insert(i);
      } else {
        std::map<UInt, std::set<UInt> >::iterator rowIt;
        rowIt = rowIndPerBlock.begin();
        for( ; rowIt != rowIndPerBlock.end(); rowIt++ ) 
          rowBlocks.insert( rowIt->first);
      }
      if( colIndPerBlock.size() == 0 )  {
        for( UInt i = 0; i < this->ncols_; ++i ) 
          colBlocks.insert(i);
      } else {
        std::map<UInt, std::set<UInt> >::iterator colIt;
        colIt = colIndPerBlock.begin();
        for( ; colIt != colIndPerBlock.end(); colIt++ ) 
          colBlocks.insert( colIt->first);
      }
      
      // Iterate over all SBM rows / cols in given index sets
      std::set<UInt>::iterator rowIt, colIt;
      
      for( rowIt = rowBlocks.begin(); 
          rowIt != rowBlocks.end(); ++rowIt ) {
        
        for( colIt = colBlocks.begin(); 
            colIt != colBlocks.end(); ++colIt ) {
          
          const UInt rowInd = *rowIt;
          const UInt colInd = *colIt;
          UInt index = ComputeIndex( rowInd, colInd);
          
          // Get hold of submatrix
          auxMat = sbmMat.subMat_[index];
          
          // Only perform addition, if submatrix is not NULL
          if ( auxMat != NULL ) {

            
            // Check whether our submatrix is NULL, in this case
            // we copy the submatrix and scale it
            if ( subMat_[index] == NULL ) {
              subMat_[index] = CopyStdMatrixObject( *auxMat );
              // This now calls the Complex version of Scale
              subMat_[index]->Scale( fac,
                                     rowIndPerBlock[rowInd],
                                     colIndPerBlock[colInd] );
            }

            // Both submatrices are not NULL so do a simple scaled addition
            // This now calls the Complex version of Scale
            subMat_[index]->Add( fac, *auxMat,
                                 rowIndPerBlock[rowInd],
                                 colIndPerBlock[colInd] );
          }
        }
      }
    }
    // Test whether downcast was successful
    catch(std::bad_alloc &e) {
      EXCEPTION( WRONG_CAST_MSG );
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
        for ( UInt i = 0; i < nrows_; i++ ) {

          // Loop over all blocks in this row
          for ( UInt j = 0; j < ncols_; j++ ) {

            // Compute i-th component of A*x
            if ( subMat_[ComputeIndex(i,j)] != NULL ) {
              subMat_[ComputeIndex(i,j)]->MultAdd( sbm_x(j), sbm_r(i) );
            }
          }
        }

        // Subtract A*x from b
        if(sbm_r.GetEntryType() == BaseMatrix::EntryType::COMPLEX) sbm_r.Axpy( (Complex)-1.0, sbm_b );
        else sbm_r.Axpy( -1.0, sbm_b );
      }

      // -----------------
      //  Symmetric case:
      // -----------------

      else {
        // Loop over all block-rows
        for ( UInt i = 0; i < nrows_; i++ ) {

          // Loop over all sub-diagonal blocks in this row
          // Add contribution of sub-diagonal blocks to the i-th component
          // of A * x using S(i,j) * x = S(j,i)^T * x
          for ( UInt j = 0; j < i; j++ ) {
            if ( subMat_[ComputeIndex(j,i)] != NULL ) {
              subMat_[ComputeIndex(j,i)]->MultTAdd( sbm_x(j), sbm_r(i) );
            }
          }

          // Add contribution of diagonal block
          if ( subMat_[ComputeIndex(i,i)] != NULL ) {
            subMat_[ComputeIndex(i,i)]->MultAdd( sbm_x(i), sbm_r(i) );
          }

          // Add contribution of super-diagonal blocks
          for ( UInt j = i+1; j < ncols_; j++ ) {
            if ( subMat_[ComputeIndex(i,j)] != NULL ) {
              subMat_[ComputeIndex(i,j)]->MultAdd( sbm_x(j), sbm_r(i) );
            }
          }
        }

        // Subtract A*x from b
        sbm_r.Axpy( -1.0, sbm_b );
      }
    }

    // Test whether downcast was successful
    catch(std::bad_alloc &e) {
      EXCEPTION( WRONG_CAST_MSG );
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
        for ( UInt i = 0; i < nrows_; i++ ) {
          for ( UInt j = 0; j < ncols_; j++ ) {
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
        for ( UInt i = 0; i < nrows_; i++ ) {

          // Add contribution of sub-diagonal blocks to the i-th component
          // of A * x using S(i,j) * x = S(j,i)^T * x
          for ( UInt j = 0; j < i; j++ ) {
            if ( subMat_[ComputeIndex(j,i)] != NULL ) {
              subMat_[ComputeIndex(j,i)]->MultTAdd( vec_m(j), vec_r(i) );
            }
          }

          // Add contribution of diagonal block
          if ( subMat_[ComputeIndex(i,i)] != NULL ) {
            subMat_[ComputeIndex(i,i)]->MultAdd( vec_m(i), vec_r(i) );
          }

          // Add contribution of super-diagonal blocks
          for ( UInt j = i+1; j < ncols_; j++ ) {
            if ( subMat_[ComputeIndex(i,j)] != NULL ) {
              subMat_[ComputeIndex(i,j)]->MultAdd( vec_m(j), vec_r(i) );
            }
          }
        }
      }
    }

    // Test whether downcast was successful
    catch(std::bad_alloc &e) {
      EXCEPTION( WRONG_CAST_MSG );
    }

  }


  // **********************************************************
  //   Perform a matrix-vector multiplication with transposed
  //   of matrix and add result to vector
  // **********************************************************
  void SBM_Matrix::MultTAdd( const BaseVector& mvec, BaseVector& rvec ) const{


    WARN("SBM_Matrix::MultTAdd has not yet been tested, my contain bugs");

    // Downcast BaseVectors to SBM_Vectors
    try {
      const SBM_Vector& vec_m = dynamic_cast<const SBM_Vector&>(mvec);
      SBM_Vector& vec_r = dynamic_cast<SBM_Vector&>(rvec);

      // ---------------------
      //  Non-symmetric case:
      // ---------------------

      if ( amSymm_ == false ) {

        // Work on one column after the other
        for ( UInt j = 0; j < ncols_; j++ ) {
          for ( UInt i = 0; i < nrows_; i++ ) {
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
      EXCEPTION( WRONG_CAST_MSG );
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

        for ( UInt i = 0; i < nrows_; i++ ) {
          for ( UInt j = 0; j < ncols_; j++ ) {
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
                  EXCEPTION( "Pattern mis-match in SBM_Vectors!" );
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
        EXCEPTION( "MultSub not implemented for symmetric SBM_Matrices, "
                 << "since MultTSub is not for normal matrices!" );
      }
    }

    // Test whether downcast was successful
    catch(std::bad_alloc &e) {
      EXCEPTION( WRONG_CAST_MSG );
    }
  }


  // **********
  //   Export
  // **********
  void SBM_Matrix::Export( const std::string& fname,
                           BaseMatrix::OutputFormat format,
                           const std::string& comment) const {


    std::stringstream ss;
    std::string outFile;

    for ( UInt j = 0; j < ncols_; j++ ) {
      for ( UInt i = 0; i < nrows_; i++ ) {
        // construct file name
        //if(ncols_ > 1 || nrows_ > 1)
        ss << fname << '_' << i << '_' << j;
        outFile =  ss.str();
        ss.str(std::string());
        // export sub-matrix
        if ( subMat_[ComputeIndex(i,j)] != NULL ) {
          subMat_[ComputeIndex(i,j)]->Export( outFile.c_str(), format, comment );
        }
      }
    }
  }

  void SBM_Matrix::Export_MultHarm( const std::string& fname,
                                    BaseMatrix::OutputFormat format,
                                    const UInt& N,
                                    const StdVector<UInt>& sbmInd,
                                    const std::string& comment) const {


    std::stringstream ss;
    std::string outFile;
    UInt sbmRow, sbmCol;

    for(auto i : sbmInd){
      StdVector<UInt> c = DeflattenIndex(i);
      sbmRow = c[0];
      sbmCol = c[1];
      ss << fname << '_' << sbmRow << '_' << sbmCol;
      outFile =  ss.str();
      ss.str(std::string());
      // export sub-matrix
      if ( subMat_[ComputeIndex(sbmRow, sbmCol)] != NULL ) {
        subMat_[ComputeIndex(sbmRow, sbmCol)]->Export( outFile.c_str(), format, comment );
      }
    }


  }

  // ******************
  //   GetMemoryUsage
  // ******************
  Double SBM_Matrix::GetMemoryUsage() const {

    StdMatrix *stdMat = NULL;
    Double sum = 0.0;
    
    // distinguish symmetry case
    if( amSymm_ ) {
      for ( UInt i = 0; i < nrows_; i++ ) {
        for ( UInt j = 0; i < ncols_; j++ ) {
          stdMat = subMat_[ComputeIndex(i,j)];
          if ( stdMat != NULL ) {
            sum += stdMat->GetMemoryUsage();
          }
        }
      }
    } else {
      // only upper-diagonal and diagonal
      for ( UInt i = 0; i < nrows_; i++ ) {
        for ( UInt j = i; i < ncols_; j++ ) {
          stdMat = subMat_[ComputeIndex(i,j)];
          if ( stdMat != NULL ) {
            sum += stdMat->GetMemoryUsage();
          }
        }
      }
    }
    return sum;
  }

  std::string SBM_Matrix::Dump() const {
    std::stringstream ss;
    ss << "number of sub matrices: " << subMat_.GetSize();

    return ss.str();
  }


  // **************
  //   GetMaxDiag
  // **************
  Double SBM_Matrix::GetMaxDiag() const {


    Double retVal = 0.0;
    Double auxVal = 0.0;
    StdMatrix *stdMat = NULL;

    for ( UInt i = 0; i < nrows_; i++ ) {
      stdMat = subMat_[ComputeIndex(i,i)];
      if ( stdMat != NULL ) {
        auxVal = stdMat->GetMaxDiag();
        retVal = auxVal > retVal ? auxVal : retVal;
      }
    }

    return retVal;
  }

  // explicit template instantiation
/*
  template void SBM_Matrix::Add(const Double fac, const BaseMatrix& mat,
                                std::map<UInt, std::set<UInt> >& rowIndPerBlock,
                                std::map<UInt, std::set<UInt> >& colIndPerBlock);
  template void SBM_Matrix::Add(const Complex fac, const BaseMatrix& mat,
                                std::map<UInt, std::set<UInt> >& rowIndPerBlock,
                                std::map<UInt, std::set<UInt> >& colIndPerBlock);
  template void SBM_Matrix::Add( const Double fac, const BaseMatrix &mat );
  template void SBM_Matrix::Add( const Complex fac, const BaseMatrix &mat );
*/
}
