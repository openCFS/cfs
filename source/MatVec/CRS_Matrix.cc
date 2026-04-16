// the following headers are required for ExportMatrixMarket()
#include <cstdio>
#include <sstream>


#include "CRS_Matrix.hh"
#include "SCRS_Matrix.hh"
#include "opdefs.hh"
#include "Utils/SyncAccess.hh"
#include "Utils/ToolsFull.hh"

#ifdef USE_MKL
#include <mkl_spblas.h>
#endif

// Implementation of methods for the compressed row storage matrix class


namespace CoupledField {


  // *********************************
  //   Constructor using CoordFormat
  // *********************************
  template<typename T>
  CRS_Matrix<T>::CRS_Matrix( CoordFormat<T> &sparseMat, bool sort )
    : diagPtr_( NULL ) {

    // Test, if the matrix is stored in symmetric format.
    // If yes, issue a warning, since we due not expand
    // it to a full matrix currently.
    bool isSymm = sparseMat.GetSymmStorageFlag();
    if ( isSymm == true ) {
      WARN("Input matrix is stored in symmetric format!\n" <<
           "We cannot expand this! Creating a matrix containing " <<
           "only the upper triangle.\n Maybe you should use " << 
           "a SCRS_Matrix?");
    }

    // Set general matrix pattern and dimension info
    this->nrows_ = sparseMat.GetNumRows();
    this->ncols_ = sparseMat.GetNumCols();
    this->nnz_   = sparseMat.GetNumEntries();

    // Allocate memory
    NEWARRAY( data_, T, this->nnz_ );
    NEWARRAY( rowPtr_, UInt, this->nrows_ + 1 );
    NEWARRAY( colInd_, UInt, this->nnz_ );

    // Generate auxilliary vector
    UInt *auxVec = NULL;
    NEWARRAY( auxVec, UInt, this->nrows_ );

    // Determine the number of non-zero entries in the different rows
    UInt *entriesPerRow = auxVec;
    for ( UInt k = 0; k < this->nrows_; k++ ) {
      entriesPerRow[k] = 0;
    }
    for ( UInt i = 0; i < this->nnz_; i++ ) {
      //check if we have some special equation numbers...
      entriesPerRow[sparseMat.ridx(i)]++;
    }

    // Now prepare the row pointer array
    rowPtr_[0] = 0;
    for ( UInt k = 0; k < this->nrows_; k++ ) {
      rowPtr_[k+1] = rowPtr_[k] + entriesPerRow[k];
    }

    // Make a copy of values and column indices
    UInt *pos = auxVec;
    int curRow;
    for ( UInt i = 0; i < (UInt)this->nrows_; i++ ) {
      pos[i] = 0;
    }
    for ( UInt j = 0; j < this->nnz_; j++ ) {
      curRow = sparseMat.ridx(j);
      data_  [ rowPtr_[curRow] + pos[curRow] ] = sparseMat.val (j);
      colInd_[ rowPtr_[curRow] + pos[curRow] ] = sparseMat.cidx(j);
      pos[ curRow ]++;
    }

    // Destroy auxilliary vector
    delete [] ( auxVec );  auxVec  = NULL;

    // We do not expect the input matrix to be sorted in anyway
    currentLayout_ = UNSORTED;

    // By default we generate a matrix in LEX sub-format
    // Note: This will also setup the diagPtr_ array
    NEWARRAY( diagPtr_, UInt, this->nrows_ );
    if(sort){
      ChangeLayout( CRS_Matrix<T>::LEX );
    }
    
    // Set pattern pool pointer to NULL, since we allocated pattern
    // ourselves
    patternPool_ = NULL;
    patternID_ = NO_PATTERN_ID;
  }


  // ********************
  //   Copy Constructor
  // ********************
  template<typename T>
  CRS_Matrix<T>::CRS_Matrix( const CRS_Matrix<T> &origMat ) {

    colInd_           = NULL;
    rowPtr_           = NULL;
    diagPtr_          = NULL;

    currentLayout_ = UNSORTED;

    patternPool_ = NULL;
    patternID_ = NO_PATTERN_ID;


    // Set basic size informations
    this->nnz_        = origMat.nnz_;
    this->nrows_      = origMat.nrows_;
    this->ncols_      = origMat.ncols_;

    // ---------------------
    // Pattern is not shared
    // ---------------------
    if ( origMat.patternPool_ == NULL ) {

      // Allocate memory for internal arrays
      NEWARRAY( colInd_ , UInt, this->nnz_        );
      NEWARRAY( rowPtr_ , UInt, this->nrows_ + 1  );
      NEWARRAY( diagPtr_, UInt, this->nrows_      );
      NEWARRAY( data_   , T      , this->nnz_        );

      // Copy information
      for (UInt i = 0; i < this->nnz_; i++ ) {
        data_[i]   = origMat.data_[i];
        colInd_[i] = origMat.colInd_[i];
      }

      for (UInt i = 0; i < this->nrows_; i++ ) {
        rowPtr_[i]  = origMat.rowPtr_[i];
        diagPtr_[i] = origMat.diagPtr_[i];
      }
      rowPtr_[this->nrows_] = origMat.rowPtr_[this->nrows_];

      // Copy layout flag
      currentLayout_ = origMat.currentLayout_;
      
      // Set pattern pool pointer to NULL, since we allocated pattern
      // ourselves
      patternPool_ = NULL;
      patternID_ = NO_PATTERN_ID;
      
    } else {
      
      // -----------------
      // Pattern is shared
      // -----------------
      // Copy sparsity pattern from poool
      this->SetSparsityPattern( origMat.patternPool_, origMat.patternID_ );
      
      // Generate copy of data array
      NEWARRAY( data_, T, this->nnz_ );
      for ( UInt i = 0; i < this->nnz_; i++ ) {
        data_[i]   = origMat.data_[i];
      }
      // Copy layout flag
      currentLayout_ = origMat.currentLayout_;
    }
  }

  // Create CRS matrix from given input data
  template<typename T>
      CRS_Matrix<T>::CRS_Matrix(UInt nr, UInt nc, UInt nnz, UInt* row_ptr, UInt* col_ptr, T* data_ptr) {

        colInd_           = NULL;
        rowPtr_           = NULL;
        diagPtr_          = NULL;

        currentLayout_ = UNSORTED;

        patternPool_ = NULL;
        patternID_ = NO_PATTERN_ID;


        // Set basic size informations
        this->nnz_        = nnz;
        this->nrows_      = nr;
        this->ncols_      = nr;

        // Allocate memory for internal arrays
        NEWARRAY( colInd_ , UInt, this->nnz_        );
        NEWARRAY( rowPtr_ , UInt, this->nrows_ + 1  );
        NEWARRAY( diagPtr_, UInt, this->nrows_      );
        NEWARRAY( data_   , T      , this->nnz_        );

        // Copy information
        for (UInt i = 0; i < this->nnz_; i++ ) {
          data_[i]   = data_ptr[i];
          colInd_[i] = col_ptr[i];
        }
        for (UInt i = 0; i < this->nrows_; i++ ) {
          rowPtr_[i]  = row_ptr[i];
        }
        rowPtr_[this->nrows_] = row_ptr[this->nrows_];

        // Copy layout flag
        //currentLayout_ = origMat.currentLayout_;

        // Set pattern pool pointer to NULL, since we allocated pattern
        // ourselves
        patternPool_ = NULL;
        patternID_ = NO_PATTERN_ID;
      }

  // Create CRS matrix from given input data
  template<typename T>
      CRS_Matrix<T>::CRS_Matrix(UInt nr, UInt nc, UInt nnz, const UInt* row_ptr, const UInt* col_ptr, T* data_ptr) {

        colInd_           = NULL;
        rowPtr_           = NULL;
        diagPtr_          = NULL;

        currentLayout_ = UNSORTED;

        patternPool_ = NULL;
        patternID_ = NO_PATTERN_ID;


        // Set basic size informations
        this->nnz_        = nnz;
        this->nrows_      = nr;
        this->ncols_      = nr;

        // Allocate memory for internal arrays
        NEWARRAY( colInd_ , UInt, this->nnz_        );
        NEWARRAY( rowPtr_ , UInt, this->nrows_ + 1  );
        NEWARRAY( diagPtr_, UInt, this->nrows_      );
        NEWARRAY( data_   , T      , this->nnz_        );

        // Copy information
        for (UInt i = 0; i < this->nnz_; i++ ) {
          data_[i]   = data_ptr[i];
          colInd_[i] = col_ptr[i];
        }
        for (UInt i = 0; i < this->nrows_; i++ ) {
          rowPtr_[i]  = row_ptr[i];
        }
        rowPtr_[this->nrows_] = row_ptr[this->nrows_];

        // Copy layout flag
        //currentLayout_ = origMat.currentLayout_;

        // Set pattern pool pointer to NULL, since we allocated pattern
        // ourselves
        patternPool_ = NULL;
        patternID_ = NO_PATTERN_ID;
      }

  // *************************
    //   Copy SCRS Matrix to CRS Matrix
    // *************************
  template<typename T>
    CRS_Matrix<T>::CRS_Matrix( const SCRS_Matrix<T> &origMat ) {

      colInd_           = NULL;
      rowPtr_           = NULL;
      diagPtr_          = NULL;

      currentLayout_ = UNSORTED;

      patternPool_ = NULL;
      patternID_ = NO_PATTERN_ID;


      // Set basic size informations
      this->nrows_      = origMat.GetNumRows();
      this->ncols_      = origMat.GetNumCols();

      if (this->nrows_ != this->ncols_) {
        EXCEPTION("SCRS matrix is not a square matrix!");
      }
      const UInt * row_ptr = origMat.GetRowPointer();
      const UInt * col_ptr = origMat.GetColPointer();
      const T * data_ptr = origMat.GetDataPointer();

      //setup temporary vector
      std::vector<UInt> tmp_vec (this->nrows_ + 1,0);

      // Compute number of elements in strict right upper part and save it in tmp_vec
      UInt j;
      for (UInt i = 0; i < this->nrows_;i++) {
        for (UInt k = row_ptr[i]; k < row_ptr[i+1];k++) {
          j = col_ptr[k];
          if (i<j) {
            tmp_vec[j+1]++;
          }
        }
      }

      // save the last nrows entries of tmp_vec in ent_col
      std::vector<UInt> ent_col(this->nrows_,0);
      UInt sum_ent_col = 0;
      for (UInt i = 0; i < this->nrows_;i++) {
        ent_col[i] = tmp_vec[i+1];
        sum_ent_col += tmp_vec[i+1];
      }


      // Find addresses of first elements of output matrix
     for (UInt i = 0; i < this->nrows_;i++) {
       j = row_ptr[i+1] - row_ptr[i];
       tmp_vec[i+1] = tmp_vec[i] + tmp_vec[i+1] + j;
     }
     this->nnz_        = tmp_vec[this->nrows_];

     // Allocate memory for internal arrays
     NEWARRAY( colInd_ , UInt, this->nnz_        );
     NEWARRAY( rowPtr_ , UInt, this->nrows_ + 1  );
     NEWARRAY( diagPtr_, UInt, this->nrows_      );
     NEWARRAY( data_   , T      , this->nnz_        );

     std::vector<Integer> col_tmp (this->nnz_,-1);


     // Copy upper part of the matrix in reverse order tmp_ is the diagonal pointer after this loop
     j = tmp_vec[this->nrows_];
     UInt tmp2;
     UInt last, first;
     for (UInt i = this->nrows_; i-- > 0;) {
       last = row_ptr[i+1];
       first = row_ptr[i];
       rowPtr_[i+1] = j;
       tmp2 = tmp_vec[i+1];
       j =  tmp_vec[i];
       for (int k = last-1; k >= (int) first; k--) {
         tmp2--;
         data_[tmp2] = data_ptr[(UInt) k];
         col_tmp[tmp2] = (Integer) col_ptr[(UInt) k];
       }
       diagPtr_[i] = tmp2;
       tmp_vec[i+1] = tmp2;
     }
     rowPtr_[0] = 0;

     // Determine upper off diagonal element data and column index in columnwise order and
     // save them in dat_index and col_index
     std::vector<UInt> col_index(sum_ent_col,0), count_vec(this->nrows_,0);
     std::vector<T> dat_index(sum_ent_col,0.);
     UInt count = 0, row_count =0;
     // setup index counter vector
     for (UInt i = 0; i <this->nrows_;i++) {
       count_vec[i] = count;
       count += ent_col[i];
     }
     for (UInt i = 0;i < origMat.GetNumEntries();i++) {
         if (i >= row_ptr[row_count])
           row_count++;
         if (row_count - 1 == col_ptr[i])
           count ++;
         else {
           j = col_ptr[i];
           dat_index[count_vec[j]] = data_ptr[i];
           col_index[count_vec[j]] = row_count -1;
           count_vec[j]++;
         }
     }

     //Copy strict lower part to crs matrix
     count = 0;
     Integer jj;
     for (UInt i = 0; i < this->nrows_;i++) {
       for (UInt k = rowPtr_[i]; k < rowPtr_[i+1];k++) {
         jj = col_tmp[k];
         if (jj != -1) {
           break;
         }
         data_[k] = dat_index[count];
         count++;
       }
     }

     // set column indices for remaining elements in lower matrix
     count = 0;
     for (UInt i = 0; i < this->nnz_;i++) {
       if (col_tmp[i] == -1) {
         colInd_[i] = col_index[count];
         count++;
       } else {
         colInd_[i] = col_tmp[i];
       }
     }
}

#ifdef USE_MULTIGRID

  // ******************************
  //   Construct from a prematrix
  // ******************************
  template <typename T>
  CRS_Matrix<T>::CRS_Matrix( UInt nrows, UInt ncols,
			     const PreMatrix<T>& premat )
    : diagPtr_( NULL ) {


    colInd_           = NULL;
    rowPtr_           = NULL;
    diagPtr_          = NULL;
    data_             = NULL;
    this->nnz_              = 0;
    this->nrows_            = 0;
    this->ncols_            = 0;
    currentLayout_    = CRS_Matrix<T>::UNSORTED;

    // Resize matrix
    SetSize( nrows, ncols, premat.GetNumNonzeros() );

    // Make a consistency check
    // Why do we need to specify nrows and ncols anyway?
    if ( nrows != premat.GetNumRows() ) {
      WARN( "CRS_Matrix: PreMatrix has mismatched size " );
    }

    rowPtr_[0] = 0;

    UInt rowlength;
    UInt position = 0;

    for( Integer i = 0; i < this->nrows_; i++ ) {

      rowlength = premat.GetRowSize( i );

      // copy column indices
      memcpy( colInd_ + position,
	      premat.GetRowCols(i) + 1,
	      rowlength * sizeof(Integer) );

      // copy entries
      memcpy( this->data_ + position,
	      premat.GetRowData(i)+1, rowlength * sizeof(T) );

      // set row length
      rowPtr_[i+1] = rowPtr_[i] + rowlength;

      // update the position
      position += rowlength;
    }

    // Try to determine sub-format (will not detect LEX!)
    if ( premat.GetLayoutFlag() == true ) {
      currentLayout_ = CRS_Matrix<T>::LEX_DIAG_FIRST;
    }
    else {
      currentLayout_ = CRS_Matrix<T>::UNSORTED;
    }

    // Generate array with diagonal entries
    NEWARRAY( diagPtr_, UInt, this->nrows_ );
    FindDiagonalEntries();
  }

#endif

  // ************************
  //   Set Sparsity Pattern
  // ************************
  template<typename T>
  void CRS_Matrix<T>::SetSparsityPattern( BaseGraph &graph ) {

    // NOTE:
    //
    // The implementation of this routine expects the GetGraphMethod
    // of the BaseGraph to return the list of nodes (including the
    // node itself, in case of a self-reference) in lexicographic
    // ordering, i.e. in increasing number of node indices.

    // Check that no pattern was allocated
     if ( rowPtr_ != NULL || colInd_ != NULL || patternPool_ != NULL ) {
       EXCEPTION( "There seems to already be a sparsity pattern!" );
     }
     
    UInt rs;
    UInt *index;

#ifdef DEBUG_CRS_MATRIX
    (*debug) << "The Matrix Graph has size " << graph.GetSize() << std::endl;
    graph.Print( *debug );
#endif
    
    // Allocate memory for row and diagonal pointers
    NEWARRAY( rowPtr_ , UInt, this->nrows_ + 1 );
    NEWARRAY( diagPtr_, UInt, this->nrows_     );
    rowPtr_[0]  = 0;
    diagPtr_[0] = 0;

    // Allocate memory for column indices
    NEWARRAY( colInd_, UInt, this->nnz_ );


    // loop over all rows
    for ( UInt i = 0; i < this->nrows_; i++ ) {

      // get neighbours (indices and numbers) of i-th node
      index = graph.GetGraphRow(i);
      rs = graph.GetRowSize(i);

      // set next row pointer (beginning of next row)
      rowPtr_[i+1] = rowPtr_[i] + rs;

      // Assume that there is no diagonal entry (and correct
      // it below, if we find one)
      diagPtr_[i] = 0;

      // set column indices of non-zero entries in this row
      for ( UInt k = 0; k < rs; k++ ) {
        colInd_[rowPtr_[i]+k] = index[k];

        // look out for diagonal entry
        if ( index[k] == i ) {
          diagPtr_[i] = rowPtr_[i]+k;
        }
      }
    }

    // setting the sparsity pattern using the information provided by
    // the BaseGraph leads to a matrix in LEX layout
    currentLayout_ = CRS_Matrix<T>::LEX;

  }
  
  // **********************
  //   SetSparsityPattern
  // **********************
  template<typename T>
  void CRS_Matrix<T>::SetSparsityPattern( PatternPool *patternPool,
                                          PatternIdType patternID ) {

    // Safety checks
    if ( patternPool == NULL ) {
      EXCEPTION( "CRS_Matrix<T>::SetSparsityPattern: patternPool = NULL! "
               << "This cannot work ;-)" );
    }
    if ( rowPtr_ != NULL || colInd_ != NULL || patternPool_ != NULL ) {
      EXCEPTION( "There seems to already be a sparsity pattern!" );
    }

    // Set information
    patternPool_ = patternPool;
    patternID_   = patternID;

    // Get hold of pattern
    BaseSparsityPattern *basePattern = NULL;
    basePattern = patternPool_->RegisterUser( patternID );

    CRS_Pattern *myPattern = NULL;
    myPattern = reinterpret_cast<CRS_Pattern*>( basePattern );
    if ( myPattern == NULL ) {
      EXCEPTION( WRONG_CAST_MSG );
    }

    // Extract pointers
    colInd_  = myPattern->cidx_;
    rowPtr_  = myPattern->rptr_;
    diagPtr_ = myPattern->diagPtr_;
  }



  // **********************
  //   SetSparsityPatternData
  // **********************
  template<typename T>
  void CRS_Matrix<T>::SetSparsityPatternData( const StdVector<UInt>& rowP,
                                               const StdVector<UInt>& colI,
                                               const StdVector<T>& data){
    // Check that no pattern was allocated
    if ( rowPtr_ != NULL || colInd_ != NULL || patternPool_ != NULL ||
        this->nrows_ != (rowP.GetSize() - 1) ) {
      EXCEPTION( "There seems to already be a sparsity pattern!" );
    }
    //if(this->nrows_ != (rowP.GetSize() - 1) ) EXCEPTION("CRS_Matrix: rowPointer-1 has other size than number of rows!!")
    this->nrows_ = rowP.GetSize() - 1;
    this->nnz_ = colI.GetSize();

    // Allocate memory for row pointers and initialise first one
    NEWARRAY( rowPtr_, UInt, rowP.GetSize() );

    // Allocate memory for column indices
    NEWARRAY( colInd_, UInt, colI.GetSize() );

    // Allocate memory for the data vector
    delete[] data_;
    NEWARRAY( data_ , T, colI.GetSize() );

    UInt maxCol = 0;

    for (UInt i = 0; i < rowP.GetSize(); ++i ) {
      rowPtr_[i] = rowP[i];
    }
    for(UInt i = 0; i < colI.GetSize(); ++i ){
//      std::cout << colI.GetPointer();
//      std::cout << "\n";
//      std::cout << data_;
//      std::cout << "\n";
    	//std::cout << *(colI.GetPointer()+i);//TODO: for debugging, remove after
      //std::cout << "\n";
    	//std::cout << colI[i];
//      std::cout << "\n";
      colInd_[i] = colI[i];
      data_[i] = data[i];

      if(colI[i] > maxCol) maxCol = colI[i];
    }

    // Diagonal pointer only if the matrix has square shape, otherwise what's the diagonal?
    if( maxCol == rowP.GetSize() - 2){
      // Allocate memory for diagonal indices
      NEWARRAY( diagPtr_, UInt, this->nrows_ );

      // fill diagonal indices
      for(UInt i = 0; i < this->nrows_; ++i){
        for(UInt j = rowPtr_[i]; j < rowPtr_[i + 1]; ++j){
          if(colInd_[j] == i){
            diagPtr_[i] = j;
            break;
          }
        }
      }
    }
  }

  // *************************
   //   TransferPatternToPool
   // *************************
   template<typename T> PatternIdType
   CRS_Matrix<T>::TransferPatternToPool( PatternPool *patternPool ) {


     // Safety checks
     if ( patternPool == NULL ) {
       EXCEPTION( "CRS_Matrix<T>::TransferPatternToPool: "
                << "patternPool = NULL! This cannot work ;-)" );
     }
     if ( patternPool_ != NULL ) {
       EXCEPTION( "CRS_Matrix<T>::TransferPatternToPool: "
                << "Do not own pattern, so I refuse to transfer it!" );
     }

     // Generate pattern object for putting into the pool
     CRS_Pattern *myPattern = new CRS_Pattern;
     myPattern->cidx_    = colInd_;
     myPattern->rptr_    = rowPtr_;
     myPattern->diagPtr_ = diagPtr_;

     // Put pattern into the pool
     patternID_   = patternPool->InsertPattern( myPattern );
     patternPool_ = patternPool;
     myPattern = NULL;

     // We still use the pattern, so register ourselves with pool
     patternPool->RegisterUser( patternID_ );

     // Return identifier
     return patternID_;
   }
   
  // *************************
  //   CopySparstiyPattern
  // *************************
  template<typename T>
  void CRS_Matrix<T>::SetSparsityPattern( const StdMatrix &srcMat ) {
    if ( patternPool_ != NULL ) {
      EXCEPTION("Setting of sparsity pattern not allowed if a pattern pool is present.");
    }
    
    if( srcMat.GetEntryType() == BaseMatrix::DOUBLE ) 
    {
      const CRS_Matrix<Double>& mat = dynamic_cast< const CRS_Matrix<Double>& >(srcMat);
      const UInt* srcColInd = mat.GetColPointer();
      const UInt* srcRowPtr = mat.GetRowPointer();
      const UInt* srcDiagPtr = mat.GetDiagPointer();
      
      if ( colInd_ == NULL )
    	  NEWARRAY( colInd_ , UInt, this->nnz_        );
      if ( rowPtr_ == NULL )
            NEWARRAY( rowPtr_ , UInt, this->nrows_ + 1  );
      if ( diagPtr_ == NULL )
            NEWARRAY( diagPtr_, UInt, this->nrows_      );
      if ( data_ == NULL)
            NEWARRAY( data_   , T      , this->nnz_        );

      // Copy information
      for (UInt i = 0; i < this->nnz_; i++ ) {
        colInd_[i] = srcColInd[i];
      }
      
      for (UInt i = 0; i < this->nrows_; i++ ) {
        rowPtr_[i]  = srcRowPtr[i];
        diagPtr_[i] = srcDiagPtr[i];
      }
      rowPtr_[this->nrows_] = srcRowPtr[this->nrows_];
      
      // Copy layout flag
      //      currentLayout_ = mat.currentLayout_;
    } else 
    {
      EXCEPTION("CRS_Matrix<T>::SetSparsityPattern not yet implemented for complex matrices.");
    }
    
  }

  // ********************************
  //   Matrix-Vector multiplication
  // ********************************
  template<typename T>
  inline void CRS_Matrix<T>::Mult( const Vector<T> &mvec,
                                   Vector<T> &rvec ) const {
    UInt rs, j, k;
    Integer i = 0;
    T sum;

#pragma omp parallel for private(sum,k,rs,j)
    for ( i = 0; i < (Integer) this->nrows_; i++ ) {
      sum = 0;
      k  = rowPtr_[i];
      rs = rowPtr_[i+1] - k;
      for ( j = 0; j < rs; j++ ) {
        sum += data_[k] * mvec[colInd_[k]];
        k++;
      }
      rvec[i] = sum;
    }
  }

  template<>
  void CRS_Matrix<Double>::Mult_type( const Vector<Complex> &mvec,
                                    Vector<Complex> &rvec ) {

    UInt i, j, rs, k;
    Double sum_re, sum_im;

    for ( i = 0; i < this->nrows_; i++ ) {
      sum_re = 0;
      sum_im = 0;
      k  = rowPtr_[i];
      rs = rowPtr_[i+1] - k;
      for ( j = 0; j < rs; j++ ) {
        sum_re += data_[k] * mvec[colInd_[k]].real();
        sum_im += data_[k] * mvec[colInd_[k]].imag();
        k++;
      }
      rvec[i].real(sum_re);
      rvec[i].imag(sum_im);
    }

  }


  // ***********
  //   MultAdd
  // ***********
  template<typename T>
  inline void CRS_Matrix<T>::MultAdd( const Vector<T> & mvec,
                                      Vector<T> & rvec ) const {


    UInt rs, j, k;
    Integer i = 0;
    T sum;

#pragma omp parallel for private(sum,k,rs,j)
    for ( i = 0; i < (Integer) this->nrows_; i++ ) {
      sum = 0;
      k  = rowPtr_[i];
      rs = rowPtr_[i+1]-k;
      for ( j = 0; j < rs; j++ ) {
        sum += data_[k] * mvec[colInd_[k]];
        k++;
      }
      rvec[i] += sum;
    }
  }


  template<>
  void CRS_Matrix<Double>::MultAdd_type( const Vector<Complex> & mvec,
                                      Vector<Complex> & rvec ) {
    UInt i, j, rs, k;
    Double sum_re, sum_im;

    for ( i = 0; i < this->nrows_; i++ ) {
      sum_re = 0;
      sum_im = 0;
      k  = rowPtr_[i];
      rs = rowPtr_[i+1]-k;
      for ( j = 0; j < rs; j++ ) {
        sum_re += data_[k] * mvec[colInd_[k]].real();
        sum_im += data_[k] * mvec[colInd_[k]].imag();
        k++;
      }
      rvec[i].real(rvec[i].real() + sum_re);
      rvec[i].imag(rvec[i].imag() + sum_im);
    }
  }

  template<typename T>
  T CRS_Matrix<T>::MultColumnWithVec(const UInt & r, const Vector<T>& vec) const{
    T sum = 0.0;
    UInt i, b;
    for( i = rowPtr_[r]; i < rowPtr_[r+1]; ++i){
      b = colInd_[i];
      sum += data_[i] * vec[b];
    }
    return sum;
  }

  // ***********
  //   MultSub
  // ***********
  template<typename T>
  inline void CRS_Matrix<T>::MultSub( const Vector<T> &mvec,
                                      Vector<T> &rvec ) const {


    UInt rs, j, k;
    Integer i = 0;
    T sum;

    // since this function is mainly used for the off-diagonal entries
    // in the parallel case we disable profiling to avoid inconsistent data

#pragma omp parallel for private(sum,k,rs,j)
    for ( i = 0; i < (Integer) this->nrows_; i++ ) {
      sum = 0;
      k  = rowPtr_[i];
      rs = rowPtr_[i+1]-k;
      for ( j = 0; j < rs; j++ ) {
        sum += data_[k] * mvec[colInd_[k]];
        k++;
      }
      rvec[i] -= sum;
    }
  }


  // ***********
  //   CompRes
  // ***********
  template<typename T>
  inline void CRS_Matrix<T>::CompRes( Vector<T> &r, const Vector<T> &x,
                                      const Vector<T> &b ) const {

    UInt rs, j, k;
    Integer i = 0;
    T sum;

#pragma omp parallel for private(sum,k,rs,j)
    for ( i = 0; i < (Integer) this->nrows_; i++ ) {
      sum = 0;
      k  = rowPtr_[i];
      rs = rowPtr_[i+1]-k;
      for ( j = 0; j < rs; j++ ) {
        sum += data_[k] * x[colInd_[k]];
        k++;
      }
      r[i] = b[i] - sum;
    }
  }


  // **********
  //   MultT
  // **********
  template<typename T>
  inline void CRS_Matrix<T>::MultT( const Vector<T> & mvec,
                                    Vector<T> & rvec ) const {


    UInt i, j, rs, k;

    // set result vector to zero
    rvec.Init();

    // loop over matrix rows
    for ( i = 0; i < this->nrows_; i++ ) {

      // get row start and row size
      k  = rowPtr_[i];
      rs = rowPtr_[i+1] - k;

      // loop over this row
      for ( j = 0; j < rs; j++ ) {
        rvec[colInd_[k]] += data_[k] * mvec[i];
        k++;
      }
    }
  }


  template<>
  void CRS_Matrix<Double>::MultT_type( const Vector<Complex> & mvec,
                                    Vector<Complex> & rvec ) const{


    UInt i, j, rs, k;

    // set result vector to zero
    rvec.Init();

    // loop over matrix rows
    for ( i = 0; i < this->nrows_; i++ ) {

      // get row start and row size
      k  = rowPtr_[i];
      rs = rowPtr_[i+1] - k;

      // loop over this row
      for ( j = 0; j < rs; j++ ) {
        rvec[colInd_[k]].real(rvec[colInd_[k]].real() + data_[k] * mvec[i].real());
        rvec[colInd_[k]].imag(rvec[colInd_[k]].imag() + data_[k] * mvec[i].imag());
        k++;
      }
    }
  }


  // ************
  //   MultTAdd
  // ************
  template<typename T>
  inline void CRS_Matrix<T>::MultTAdd( const Vector<T> &mvec,
                                       Vector<T> &rvec ) const {

    UInt i, j, end;

    // loop over matrix rows
    for( i = 0; i < this->nrows_; i++ ) {

      // get end of row
      end = rowPtr_[i+1];

      // loop over this row
      for( j = rowPtr_[i]; j < end; j++ ) {
        rvec[colInd_[j]] += this->data_[j] * mvec[i];
      }
    }
  }
  

  template<>
  void CRS_Matrix<Double>::MultTAdd_type( const Vector<Complex> &mvec,
                                       Vector<Complex> &rvec ) const {

    UInt i, j, end;

    // loop over matrix rows
    for( i = 0; i < this->nrows_; i++ ) {

      // get end of row
      end = rowPtr_[i+1];

      // loop over this row
      for( j = rowPtr_[i]; j < end; j++ ) {
        rvec[colInd_[j]].real(rvec[colInd_[j]].real() + this->data_[j] * mvec[i].real());
        rvec[colInd_[j]].imag(rvec[colInd_[j]].imag() + this->data_[j] * mvec[i].imag());
      }
    }
  }





  // ************
  //   MultTSub
  // ************
  template<typename T>
  inline void CRS_Matrix<T>::MultTSub( const Vector<T> &mvec,
                                       Vector<T> &rvec ) const {

    UInt i, j, end;

    // loop over matrix rows
    for( i = 0; i < this->nrows_; i++ ) {

      // get end of row
      end = rowPtr_[i+1];

      // loop over this row
      for( j = rowPtr_[i]; j < end; j++ ) {
        rvec[colInd_[j]] -= this->data_[j] * mvec[i];
      }
    }
   }

  template<typename T>
  void CRS_Matrix<T>::Transpose(UInt* col_ptr, UInt* row_ptr, T* data_ptr) const
  {
     // first set the explicit row-index temporaly along the colInd_.
     UInt* tmpInd;
     NEWARRAY(tmpInd, UInt, this->nnz_);

     for(UInt i = 0; i < this->nrows_; i++)
     {
        for(UInt k = rowPtr_[i]; k < rowPtr_[i+1]; k++)
        {
          // int j = colInd_[k]; // TODO: Check if this is still needed
           tmpInd[k] = i;
        }
     }

     // now set up the column structure
     int p = 0;
     for(UInt j = 0; j < this->ncols_; j++)
     {
        // store the column index
        col_ptr[j] = p;
        // search in colInd_ for this column
        for(UInt k = 0; k < this->nnz_; k++)
        {
            // do we have our column index?
            if(colInd_[k] == j)
            {
               // copy the value and store the row-index
               data_ptr[p] = data_[k];
               row_ptr[p]  = tmpInd[k];

               p++;
            }
        }
     }

     // write the tail
     col_ptr[this->ncols_] = this->nnz_;

     delete [] (tmpInd); tmpInd = NULL;
  }


  // ****************************
  //   Print matrix to a stream
  // ****************************
  template<typename T>
  std::string CRS_Matrix<T>::ToString( char colSeparator,
                                       char rowSeparator) const {

    std::stringstream os;
    UInt i, j, k;

    for ( i = 0; i < this->nrows_; i++ ) {
      for ( k = rowPtr_[i]; k < rowPtr_[i+1]; k++ ) {
        j = colInd_[k];
        os << std::scientific << i << colSeparator 
        << j << colSeparator << data_[k] << rowSeparator;;
      }
    }

    os << std::endl;

#ifdef DEBUG_CRS_MATRIX
    os << " diagPtr_: ";
    for ( i = 0; i < this->nrows_; i++ ) {
      os << diagPtr_[i] << " ";
    }
    os << '\n' << std::endl;
#endif
    return os.str(); 
  }

  template<typename T>
  std::string CRS_Matrix<T>::Dump() const
  {
    std::stringstream ss;
    // don't use ToString() from this class but the glocal ToString() from ToolsFull.hh
    ss << " row=" << ::ToString<unsigned int>(rowPtr_, this->nrows_+1) << std::endl;
    ss << " col=" << ::ToString<unsigned int>(colInd_, this->nnz_)  << std::endl;
    ss << " val=" << ::ToString<T>(data_, this->nnz_);
    return ss.str();
  }


  // *************************
  //   Export matrix to file
  // *************************
  template<typename T>
  void CRS_Matrix<T>::ExportMatrixMarket( const std::string& fname,
                                          const std::string& comment ) const{


    // Open output file and check for errors
    FILE *fp = fopen( fname.c_str(), "w" );
    if ( fp == NULL ) {
      EXCEPTION( "Cannot open file " << fname << " for writing!" );
    }

    // ---------------------
    //   Write file header
    // ---------------------

    // Matrix Market Format Specification
    if ( this->GetEntryType() == BaseMatrix::DOUBLE ) {
      fprintf( fp, "%%%%MatrixMarket matrix coordinate real general\n" );
    }
    else if ( this->GetEntryType() == BaseMatrix::COMPLEX ) {
      fprintf( fp, "%%%%MatrixMarket matrix coordinate complex general\n" );
    }
    else {
      EXCEPTION( "Expected either DOUBLE or COMPLEX as matrix entry type" );
    }

    // User-supplied private comment
    if ( comment != "" ) {
      fprintf( fp, "%%\n%% %s\n%%\n", comment.c_str() );
    }
    else {
      fprintf( fp, "%%\n%% Matrix exported by openCFS\n%%\n" );
    }

    // Information on number of rows, columns and entries
    fprintf( fp, "%d\t%d\t%d\n", this->nrows_ , this->ncols_ ,
             this->nnz_ );


    // ------------------------
    //   Write matrix entries
    // ------------------------

    UInt i, j, nblocks;
    // UInt i, j, ib, jb, nblocks; // TODO: Check if this is still needed

    // loop over all block rows
    for ( i = 0; i < this->nrows_; i++ ) {

      // get number of blocks in i-th row
      nblocks = rowPtr_[i+1] - rowPtr_[i];

      // loop over all blocks in this row
      for ( j = 0; j < nblocks; j++ ) {

        // store row and column index
        fprintf( fp, "%6d\t%6d\t", i  + 1,
                 ( colInd_[rowPtr_[i]+j] )  + 1);
        
        // store non-zero entry
        OpType<T>::ExportEntry( this->data_[rowPtr_[i]+j],
                                0, 0, fp );
        fprintf( fp, "\n" );
      }
    }

    // close output file
    if ( fclose( fp ) == EOF ) {
      WARN("Could not close file " << fname << " after writing!");
    }
  }


  // *******************************
  //   Add value to a matrix entry
  // *******************************
  template<typename T>
  void CRS_Matrix<T>::AddToMatrixEntry( UInt i, UInt j, const T& v ){


    bool found = false;

    // Try to determine index for matrix entry at position (i,j)
    // in the data_ array and add value
    UInt l = rowPtr_[i];
    UInt u = rowPtr_[i+1];
    switch(currentLayout_){
    case UNSORTED: // we have to search linearly here
      for ( UInt k = l; k < u; k++ ) {
        if ( colInd_[k] == j ) {
          found = true;
          SyncAccess<SYNC_DATA>::AddTo(data_[k],v);
          break;
        }
      }
      break;
    case LEX_DIAG_FIRST:
      if(colInd_[l] == j){ // the diagonal exists
        if(i == j){ // we want the diagonal
          SyncAccess<SYNC_DATA>::AddTo(data_[j],v);
          found = true;
          break;
        }else{
          l++; // we do not want the diagonal element, search all others
        }
      }
      // Note: no break
    case LEX:
      // logarithmic search (this has complexity O(log(n)), n=u-l)
      u--; //instead of UInt u = rowPtr_[i+1]-1;
      while(l <= u){
        UInt k = (l+u) >> 1;
        if(colInd_[k] > j){
          u = k-1;
        }else if(colInd_[k] < j){
          l = k+1;
        }else{
	        SyncAccess<SYNC_DATA>::AddTo(data_[k],v);
          found = true;
          break;
        }
      }
    }

    // Its a serious error, if the index pair was not found
    if ( found == false ) {
      std::cerr << ToString();
      EXCEPTION( "AddToMatrixEntry: Index pair = (" << i << " , "
               << j << ") not found\n" );
    }
  }


  // *****************************
  //   Get specific matrix entry
  // *****************************
  template<typename T>
  void CRS_Matrix<T>::GetMatrixEntry( UInt i, UInt j,
                                      T &v ) const {

    bool found = false;

    // Try to determine index for matrix entry at position (i,j)
    // in the data_ array and get the entry
    UInt l = rowPtr_[i];
    UInt u = rowPtr_[i+1];
    switch(currentLayout_){
    case UNSORTED: // we have to search linearly here
      for ( UInt k = l; k < u; k++ ) {
        if ( colInd_[k] == j ) {
          found = true;
          v = data_[k];
          break;
        }
      }
      break;
    case LEX_DIAG_FIRST:
      if(colInd_[l] == j){ // the diagonal exists
        if(i == j){ // we want the diagonal
          v = data_[j];
          found = true;
          break;
        }else{
          l++; // we do not want the diagonal element, search all others
        }
      }
      // Note: no break
    case LEX:
      // logarithmic search (this has complexity O(log(n)), n=u-l)
      u--; //instead of UInt u = rowPtr_[i+1]-1;
      while(l <= u){
        UInt k = (l+u) >> 1;
        if(colInd_[k] > j){
          u = k-1;
        }else if(colInd_[k] < j){
          l = k+1;
        }else{
          v = data_[k];
          found = true;
          break;
        }
      }
    }

    // Its a serious error, if the index pair was not found
    if ( found == false ) {
      std::cerr << ToString();
      EXCEPTION( "GetMatrixEntry: Index pair = (" << i << " , "
               << j << ") not found\n" );
    }
  }

  // ***************************************************
  //   Get the row, column and data vectors of a matrix
  // ***************************************************
  template<typename T>
  void CRS_Matrix<T>::GetVectors(StdVector<UInt> *Vec_col, StdVector<UInt> *Vec_row, StdVector<T> *Vec_val ) const {
  	UInt i,j=0, k;
  	*(Vec_row->GetPointer()) =0;
  	    for ( i = 0; i < this->nrows_; i++ ) {
  	      for ( k = rowPtr_[i]; k < rowPtr_[i+1]; k++ ) {
  	        *(Vec_row->GetPointer()+i+1) =j+1;
  	        *(Vec_col->GetPointer()+j)=colInd_[k];
  	        *(Vec_val->GetPointer()+j)=data_[k];
  	        j++;
  	      }
  	    }
  }


  // *****************************
  //   Has specific matrix entry ?
  // *****************************
  template<typename T>
  bool CRS_Matrix<T>::HasMatrixEntry( UInt i, UInt j, T& v) const {
    bool found = false;
    // Try to determine index for matrix entry at position (i,j)
    UInt l = rowPtr_[i];
    UInt u = rowPtr_[i+1];
    switch(currentLayout_){
    case UNSORTED: // we have to search linearly here
      for ( UInt k = l; k < u; k++ ) {
        if ( colInd_[k] == j ) {
          found = true;
          v = data_[k];
          break;
        }
      }
      break;
    case LEX_DIAG_FIRST:
      if(colInd_[l] == j){ // the diagonal exists
        if(i == j){ // we want the diagonal
          v = data_[j];
          found = true;
          break;
        }else{
          l++; // we do not want the diagonal element, search all others
        }
      }
      // Note: no break
    case LEX:
      // logarithmic search (this has complexity O(log(n)), n=u-l)
      u--; //instead of UInt u = rowPtr_[i+1]-1;
      while(l <= u){
        UInt k = (l+u) >> 1;
        if(colInd_[k] > j){
          u = k-1;
        }else if(colInd_[k] < j){
          l = k+1;
        }else{
          v = data_[k];
          found = true;
          break;
        }
      }
    }
    return found;
  }

  // *****************************
  //   Set specific matrix entry
  // *****************************
  template<typename T>
  void CRS_Matrix<T>::SetMatrixEntry( UInt i, UInt j, const T &v,
                                      bool setCounterPart ) {

    bool found = false;

    // Try to determine index for matrix entry at position (i,j)
    // in the data_ array and set the value
    UInt l = rowPtr_[i];
    UInt u = rowPtr_[i+1];
    switch(currentLayout_){
    case UNSORTED: // we have to search linearly here
      for ( UInt k = l; k < u; k++ ) {
        if ( colInd_[k] == j ) {
          found = true;
          SyncAccess<SYNC_DATA>::Set(data_[k],v);
          break;
        }
      }
      break;
    case LEX_DIAG_FIRST:
      if(colInd_[l] == j){ // the diagonal exists
        if(i == j){ // we want the diagonal
          SyncAccess<SYNC_DATA>::Set(data_[j],v);
          found = true;
          break;
        }else{
          l++; // we do not want the diagonal element, search all others
        }
      }
      // Note: no break
    case LEX:
      // logarithmic search (this has complexity O(log(n)), n=u-l)
      u--; //instead of UInt u = rowPtr_[i+1]-1;
      while(l <= u){
        UInt k = (l+u) >> 1;
        if(colInd_[k] > j){
          u = k-1;
        }else if(colInd_[k] < j){
          l = k+1;
        }else{
          SyncAccess<SYNC_DATA>::Set(data_[k],v);
          found = true;
          break;
        }
      }
    }

    // Its a serious error, if the index pair was not found
    if ( found == false ) {
      std::cerr << ToString();
      EXCEPTION( "SetMatrixEntry: Index pair = (" << i << " , "
               << j << ") not found\n" );
    }

    // If setCounterPart is true and the index-pair denotes a non-diagonal
    // position, we do also insert the counterpart
    if ( i != j && setCounterPart == true ) {
      found = false;
      UInt l = rowPtr_[j];
      UInt u = rowPtr_[j+1];
      switch(currentLayout_){
      case UNSORTED: // we have to search linearly here
        for ( UInt k = l; k < u; k++ ) {
          if ( colInd_[k] == i ) {
            found = true;
            data_[k] = v;
            break;
          }
        }
        break;
      case LEX_DIAG_FIRST:
        if(colInd_[l] == i){ // the diagonal exists
          l++; // we do not want the diagonal element, search all others (i != j from above)
        }
        // Note: no break
      case LEX:
        // logarithmic search (this has complexity O(log(n)), n=u-l)
        u--; //instead of UInt u = rowPtr_[i+1]-1;
        while(l <= u){
          UInt k = (l+u) >> 1;
          if(colInd_[k] > i){
            u = k-1;
          }else if(colInd_[k] < i){
            l = k+1;
          }else{
            data_[k] = v;
            found = true;
            break;
          }
        }
      }

      // Its a serious error, if the index pair was not found
      if ( found == false ) {
        std::cerr << ToString();
        EXCEPTION( "SetMatrixEntry: Index pair = (" << j << " , "
                 << i << ") not found\n");
      }
    }
  }

  // **********
  //   Insert
  // **********
  template<typename T>
  void CRS_Matrix<T>::Insert(UInt row, UInt col,  UInt pos){
#ifdef DEBUG_CRS_MATRIX
    (*debug) << "CRS_Matrix: inserting ("<<row<<","<<col<<") as Entry ";
    (*debug) << pos << std::endl;
    if( rowPtr_[row+1] - rowPtr_[row] < pos ) {
      (*debug) << "ERROR: row " << row << " has got only "
               << (rowPtr_[row+1] - rowPtr_[row]) << " entries"
               << std::endl;
    }
#endif

    colInd_[rowPtr_[row]+pos-1] = col;
  }

  // *********
  //   Scale
  // *********
  template<>
  void CRS_Matrix<Double>::Scale( Double factor ) {
    for ( UInt i = 0; i < this->nnz_; i++ ) {
      data_[i] *= factor;
    }
  }
  
  template<>
  void CRS_Matrix<Complex>::Scale( Double factor ) {
    for ( UInt i = 0; i < this->nnz_; i++ ) {
      data_[i] *= factor;
    }
  }

  template<>
  void CRS_Matrix<Double>::Scale( Complex factor ) {
	  EXCEPTION("CRS: Matrix is Double; you can't multiply be a complex value");
  }


  // *********
  //   Scale
  // *********
  template<>
  void CRS_Matrix<Complex>::Scale( Complex factor ) {
    for ( UInt i = 0; i < this->nnz_; i++ ) {
      data_[i] *= factor;
    }
  }

  // ************************
  //   Scale on index subset
  // ************************
  template<typename T>
  void CRS_Matrix<T>::Scale( Double factor,
                             const std::set<UInt>& rowIndices,
                             const std::set<UInt>& colIndices ) {
    EXCEPTION("Implement me");
  }

  // ************************
  //   Add (another matrix)
  // ************************
  template<typename T>
  void CRS_Matrix<T>::Add( const Double alpha, const StdMatrix& mat ) {

      // Obtain pointer to data array of other matrix
      const T *data = dynamic_cast<const CRS_Matrix<T>&>(mat).GetDataPointer();

      // We now assume that the matrices have matching
      // dimensions and sparsity patterns
      for ( UInt i = 0; i < this->nnz_; i++ ) {
        data_[i] += alpha * data[i];
      }
  }
  template<typename T>
  void CRS_Matrix<T>::Add( const Complex alpha, const StdMatrix& mat ) {
    EXCEPTION("Complex Add is not implemented for CRS_Matrix.");
  }

  template<>
  void CRS_Matrix<Complex>::Add( const Complex alpha, const StdMatrix& mat ) {
    EXCEPTION("Complex Add is not implemented for CRS_Matrix.");
  }

  template<>
  void CRS_Matrix<Complex>::Add( const Double alpha, const StdMatrix& mat ) {

    // Check for entry type of mat
	BaseMatrix::EntryType eType = mat.GetEntryType();

	if( eType == BaseMatrix::DOUBLE ) {
	  // Obtain pointer to data array of other matrix
	  const Double *data = dynamic_cast<const CRS_Matrix<Double>&>(mat).GetDataPointer();

	  // We now assume that the matrices have matching
	  // dimensions and sparsity patterns
	  for ( UInt i = 0; i < this->nnz_; i++ ) {
	    data_[i] += alpha * Complex(data[i], 0.0 );
	  }
	} else {
	  // Obtain pointer to data array of other matrix
	  const Complex * data = dynamic_cast<const CRS_Matrix<Complex>&>(mat).GetDataPointer();

	  // We now assume that the matrices have matching
	  // dimensions and sparsity patterns
	  for ( UInt i = 0; i < this->nnz_; i++ ) {
	    data_[i] += alpha * data[i];
	  }
	}
  }
  
  // ******************************************
  //   Add (another matrix, only index subset)
  // ******************************************
  template<typename T>
  void CRS_Matrix<T>::Add( const Double alpha, const StdMatrix& mat,
                           const std::set<UInt>& rowIndices,
                           const std::set<UInt>& colIndices ) {
    // Down-cast input matrix
    const CRS_Matrix<T>& crsMat = dynamic_cast<const CRS_Matrix<T>&>(mat);

    // Obtain pointer to data array of other matrix
    const T *data = crsMat.GetDataPointer();

    // Distinguish 4 cases:
    // 1) Neither row- nor col-indices are set (i.e. take all incides)
    //    -> use standard Add methods
    // 2) Row and col-indices are set and contain all rows / cols
    //    -> use standard Add methods
    if( (rowIndices.size() == 0 && colIndices.size() == 0) ||
        (rowIndices.size() == this->nrows_ &&
         colIndices.size() == this->ncols_ ) ) {
      this->Add(alpha, mat);
      return;
    }

    // 3) Both sets are non-empty sub-sets of all indices
    //    --> loop over row / column indices to be set
    if( rowIndices.size() > 0 && colIndices.size() > 0 ) {
      std::set<UInt>::const_iterator rowIt, colIt;
      UInt k, rs;
      UInt j;
      rowIt = rowIndices.begin();

      // loop over all rows in the rowIndex set
      for( ; rowIt != rowIndices.end(); ++rowIt ) {
        k = rowPtr_[*rowIt];
        rs = rowPtr_[*rowIt+1] - k;

        // loop over all columns of the current row
        colIt = colIndices.begin();
        for ( j = 0; j < rs; j++ ) {

          // iterate over desired column indices until the
          // current column index is smaller
          while( *colIt < colInd_[k] && colIt != colIndices.end() ) 
            colIt++;
          // only perform operation, if column indices match
          if( *colIt == colInd_[k] )
            data_[k] += alpha * data[k];
          k++;
        }
      }
    } else {
      // 4) Only one set is empty: rowIndices or colIndices
      if ( rowIndices.size() > 0 ) {
        //    -> either loop over selected rows and take into account
        //       all columns
        std::set<UInt>::const_iterator rowIt;
        UInt k, rs;
        UInt j;
        rowIt = rowIndices.begin();

        // loop over all rows in the rowIndex set
        for( ; rowIt != rowIndices.end(); ++rowIt ) {
          k = rowPtr_[*rowIt];
          rs = rowPtr_[*rowIt+1] - k;

          // loop over all columns of the current row
          for ( j = 0; j < rs; j++ ) {
            data_[k] += alpha * data[k];
            k++;
          }
        }
      } else {
        //  -> or loop over all rows and take into account only
        //     selected columns
        //this is simpler right now we can assume that all matrices have the same graph
        // Obtain pointer to data array of other matrix
        const T *data = dynamic_cast<const CRS_Matrix<T>&>(mat).GetDataPointer();
        //another optimization check if data in set is continuous
        std::set<UInt>::iterator starting  = colIndices.begin();
        std::set<UInt>::iterator ending = colIndices.end();
        UInt min = *starting;
        UInt max = *colIndices.rbegin();
        UInt size = colIndices.size();
        if( (max - min) == (size-1)){
          //data is continuous
#pragma omp parallel for
          for ( Integer i = 0; i < (Integer) this->nnz_; i++ ) {
            if (colInd_[i] >= min && colInd_[i] <= max)
              data_[i] += alpha * data[i];
          }
        }else{
          // we have to searach the set for every entry (nnz_ * O(log sizeof(colIndices)))
#pragma omp parallel for
          for ( Integer i = 0; i < (Integer) this->nnz_; i++ ) {
            if (colIndices.find(colInd_[i]) != ending)
              data_[i] += alpha * data[i];
          }
        }
      } // if column set is nonzero
    } // if one set is nonzero
  }

  template<typename T>
  void CRS_Matrix<T>::Add( const Complex alpha, const StdMatrix& mat, const std::set<UInt>& rowIndices, const std::set<UInt>& colIndices ) {
    EXCEPTION("Complex Add is not implemented for CRS_Matrix.");
  }

  // *************
  //   QuickSort
  // *************
  template <typename T>
  void CRS_Matrix<T>::QuickSort( UInt *const cols, T *const data,
                                 const UInt length ) {

    if( length > 1 ) {

      const UInt last         = length-1;
      const UInt splitter     = cols[last];
      const T splitterData = data[last];

      UInt tcol, i = 0;
      T tdat;

      // split array
      for( UInt j = 0; j < last; j++ ) {

        if( cols[j] < splitter ) {

          // swap the column index
          tcol    = cols[j];
          cols[j] = cols[i];
          cols[i] = tcol;

          // swap the data
          tdat    = data[j];
          data[j] = data[i];
          data[i] = tdat;
          i++;
        }
      }

      // swap the splitter
      cols[last] = cols[i];
      cols[i]    = splitter;

      // and the corresponding data
      data[last] = data[i];
      data[i]    = splitterData;

      // quicksort the two parts
      QuickSort( cols,   data,   i++        );
      QuickSort( cols+i, data+i, length - i );
    }
  }


  // ****************
  //   ChangeLayout
  // ****************
  template <typename T>
  void CRS_Matrix<T>::ChangeLayout( const subFormat newLayout ) {


    // Check, if we must do anything at all
    if ( newLayout != currentLayout_ ) {

#ifdef DEBUG_CRS_MATRIX
      (*debug) << " Old format:\n";
      Print( *debug );
#endif

      T tdat;

      // ----------------
      //  Convert to LEX
      // ----------------
      if ( newLayout == CRS_Matrix<T>::LEX ) {

        // Use specialised routine, if we currently are in
        // LEX_DIAG_FIRST sub-format
        if ( currentLayout_ == CRS_Matrix<T>::LEX_DIAG_FIRST ) {
          SortLexDiagFirst2Lex();
          currentLayout_ = newLayout;
        }

        // Use quick-sort otherwise
        else {

          for( UInt i = 0; i < this->GetNumRows(); i++ ) {

            // sort the row entries lexicographically
            QuickSort( colInd_ +  rowPtr_[i], data_ + rowPtr_[i],
                       rowPtr_[i+1] - rowPtr_[i] );
          }

          // Now we conform to the new layout
          currentLayout_ = newLayout;

          // Re-built the diagonal entry array
          FindDiagonalEntries();
        }

#ifdef DEBUG_CRS_MATRIX
        (*debug) << " New format:\n";
        Print( *debug );
#endif
      }


      // ---------------------------
      //  Convert to LEX_DIAG_FIRST
      // ---------------------------
      else if ( newLayout == CRS_Matrix<T>::LEX_DIAG_FIRST ) {

        // Use specialised routine, if we currently are in LEX sub-format
        if ( currentLayout_ == CRS_Matrix<T>::LEX ) {
          SortLex2LexDiagFirst();
          currentLayout_ = newLayout;
        }

        // Use quick-sort otherwise
        else {
          for( UInt i = 0; i < this->GetNumRows(); i++ ) {

            // search the diagonal entry
            for( UInt ij = rowPtr_[i]; ij < rowPtr_[i+1]; ij++ ) {

              // place the diagonal entry at leading position
              if( colInd_[ij] == i ) {
                const UInt idiag = rowPtr_[i];

                // swap the column indices
                colInd_[ij]    = colInd_[idiag];
                colInd_[idiag] = i;

                // swap the data
                tdat               = this->data_[ij];
                this->data_[ij]    = this->data_[idiag];
                this->data_[idiag] = tdat;

                // sort the remaining row using quick sort
                QuickSort( colInd_ + idiag + 1,
                           this->data_   + idiag + 1,
                           rowPtr_[i+1] - idiag - 1 );

                // next row
                break;
              }
            }
          }

          // Now we conform to the new layout
          currentLayout_ = newLayout;

          // Re-built the diagonal entry array
          FindDiagonalEntries();
        }

#ifdef DEBUG_CRS_MATRIX
        (*debug) << " New format:\n";
        Print( *debug );
#endif
      }


      // ---------------------
      //  Convert to UNSORTED
      // ---------------------
      else if ( newLayout == CRS_Matrix<T>::UNSORTED ) {
        EXCEPTION( " Cannot convert CRS_Matrix to UNSORTED sub-format!" );
      }


      // ---------
      //  Default
      // ---------
      else {
        EXCEPTION( "Congratulations! You have found a missing case "
                 << "implementation!" );
      }
    }
  }


  // ******************
  //   GetMemoryUsage
  // ******************
  template <typename T>
  Double CRS_Matrix<T>::GetMemoryUsage() const {
    
    Double sum = 0.0;
    if( !patternPool_ ) {
      sum += ( (this->nrows_ + 1) // rowPtr
          + this->nrows_       // diagPtr
          + this->nnz_ )       // colInd
          * sizeof(UInt);
    }
    sum += this->nnz_ * sizeof(T); // data
    return sum;
  }
  
  // ******************
  //   Re-size matrix
  // ******************
  template<typename T>
  void CRS_Matrix<T>::SetSize( UInt nrows, UInt ncols, UInt nnz ) {


#ifdef DEBUG_CRS_MATRIX
    if ( nrows < 0 || nnz < 0 ) {
      EXCEPTION( "invalid dimensions" );
    }
#endif

    if(patternPool_ && ((this->nrows_ != nrows) || (this->nnz_ != nnz))){ // If we don't own the pattern, deregister and clear pattern
      patternPool_->DeRegisterUser(patternID_);
      rowPtr_ = nullptr;
      diagPtr_ = nullptr;
      colInd_ = nullptr;
      patternID_ = NO_PATTERN_ID;
      patternPool_ = nullptr;
    }

    this->ncols_ = ncols;
    if ( this->nrows_ != nrows ) {
      this->nrows_ = nrows;
      delete [] ( rowPtr_  );
      rowPtr_ = nullptr;
      delete [] ( diagPtr_ );
      diagPtr_ = nullptr;
      currentLayout_ = CRS_Matrix<T>::UNSORTED;
    }

    if ( this->nnz_ != nnz ) {
      this->nnz_ = nnz;
      delete [] ( colInd_ );
      colInd_ = nullptr;
      delete [] ( data_ );
      data_ = nullptr;
    }
    
    if(!rowPtr_){
      rowPtr_  = new UInt[this->nrows_ + 1];
      rowPtr_[0]  = 0;
    }
    if(!diagPtr_){
      diagPtr_ = new UInt[this->nrows_];
      diagPtr_[0] = 0;
    }
    if(!colInd_){
      colInd_ = new UInt[this->nnz_];
    }
    if(!data_){
      data_ = new T[this->nnz_];
    }
  }


  // ***********************
  //   FindDiagonalEntries
  // ***********************
  template <typename T>
  void CRS_Matrix<T>::FindDiagonalEntries() {

    // If the matrix is in LEX_DIAG_FIRST format, setting up
    // the array is pretty easy
    if ( currentLayout_ == CRS_Matrix<T>::LEX_DIAG_FIRST ) {
      for ( UInt i = 0; i < this->nrows_; i++ ) {
        diagPtr_[i] = 0;
        if ( colInd_[rowPtr_[i]] == i ) {
          diagPtr_[i] = rowPtr_[i];
        }
      }
    }

    // For other formats we must search the diagonal entries
    else {

      // loop over all rows
      for ( UInt i = 0; i < this->nrows_; i++ ) {

        // assume that there is no diagonal entry (and correct
        // it below, if we find one)
        diagPtr_[i] = 0;

        // loop over column
        for ( UInt k = rowPtr_[i]; k < rowPtr_[i+1]; k++ ) {

          // look out for diagonal entry
          if ( colInd_[k] == i ) {
            diagPtr_[i] = k;
          }
        }
      }
    }
  }


  // ************************
  //   SortLex2LexDiagFirst
  // ************************
  template<typename T>
  void CRS_Matrix<T>::SortLex2LexDiagFirst() {

    // Since the row is lexicographically ordered, we just have to move the
    // diagonal entry to the front and move the part of the row in front of
    // it to the right.

    T auxVal;

    // loop over all rows
    for ( UInt i = 0; i < this->nrows_; i++ ) {

      // Check if row has a diagonal entry (otherwise this might
      // not be sensible)
      //TODO had to disable this for AMG (zero-based), we have to check
      // if this has further consequences; maybe we can check
      // if data_[diagPtr_[i]] == 0.0 ?
      /*
      if ( diagPtr_[i] == 0 ) {
        EXCEPTION( "CRS_Matrix<T>::SortLex2LexDiagFirst: There is no "
                 << "diagonal entry in row " << i );
      }
    */

      // Store diagonal entry
      auxVal = data_[ diagPtr_[i] ];

      // shift left part of row to the right
      for ( UInt k = diagPtr_[i]; k > rowPtr_[i]; k-- ) {
        colInd_[k] = colInd_[k-1];
        data_[k] = data_[k-1];
      }

      // Put diagonal entry to the front
      data_[ rowPtr_[i] ] = auxVal;
      colInd_[ rowPtr_[i] ] = i;
      diagPtr_[i] = rowPtr_[i];
    }
  }


  // ************************
  //   SortLexDiagFirst2Lex
  // ************************
  template<typename T>
  void CRS_Matrix<T>::SortLexDiagFirst2Lex() {

    // Since the rest of the row is lexicographically ordered, we just have
    // to insert the diagonal entry at the correct position and move the
    // left part of the row to the left.

    T auxVal;
    UInt auxInd;

    // loop over all rows
    for ( UInt i = 0; i < this->nrows_; i++ ) {

      // Only, if there is a diagonal entry, do we have something to do
      if ( diagPtr_[i] == rowPtr_[i] ) {

        for ( UInt j = rowPtr_[i] + 1; j < rowPtr_[i+1]; j++ ) {

          // If next entry lies left of diagonal, then permute
          if ( colInd_[j] < i ) {
            auxVal = data_[j];
            auxInd = colInd_[j];

            colInd_[j] = colInd_[j-1];
            data_[j] = data_[j-1];

            colInd_[j-1] = auxInd;
            data_[j-1] = auxVal;

            diagPtr_[i] = j;
          }
          else {
            break;
          }
        }
      }
    }
  }


/*  template<>
  double* CRS_Matrix<Complex>::GetDataPointerReal(){
    double *r = new double[nnz_];
    for(UInt i = 0; i < nnz_; ++i) r[i] = data_[i].real();
    double *ret = r;
    delete [] r; r = NULL;
    return ret;
  }*/


  // ************************
  // Matrix-Matrix and Matrix-Matrix-Matrix Multiplication
  // ************************
#ifdef USE_MKL
  /* To avoid constantly repeating the part of code that checks inbound SparseBLAS functions' status,
       use macro CALL_AND_CHECK_STATUS */
  #define CALL_AND_CHECK_STATUS(function, error_message) do { \
      if(function != SPARSE_STATUS_SUCCESS)             \
      {                                                 \
        printf(error_message); fflush(0);                 \
      }                                                 \
  } while(0)
#endif

  template<typename T>
  void CRS_Matrix<T>::MultTriple_MKL(CRS_Matrix<T> B,
                                CRS_Matrix<Double>& A,
                                StdVector<UInt>& rPC,
                                StdVector<UInt>& cPC,
                                StdVector<T>& dPC,
                                UInt version,
                                bool setSparsity){
    EXCEPTION("CRS_Matrix<T>::MultTriple_MKL Not implemented for the general case");
  }


  template<>
  void CRS_Matrix<Double>::MultTriple_MKL(CRS_Matrix<Double> B,
                                    CRS_Matrix<Double>& A,
                                    StdVector<UInt>& rPC,
                                    StdVector<UInt>& cPC,
                                    StdVector<Double>& dPC,
                                    UInt version,
                                    bool setSparsity){

#ifdef USE_MKL
    /******************** Some dimension checks *********************/
    UInt numRowB = B.GetNumRows();
    UInt numColB = B.GetNumCols();
    UInt numRowA = A.GetNumRows();
    UInt numColA = A.GetNumCols();
    if(numColB != numRowB){
      EXCEPTION("MultTriple_MKL: Matrix B assumed to be square");
    }
    if(version==1 && numRowA!=numColB){ EXCEPTION("MultTriple_MKL: Matrices have wrong dimension");}
    if(version==2 && numColA!=numRowB){ EXCEPTION("MultTriple_MKL: Matrices have wrong dimension");}

    sparse_matrix_t A_MKL = NULL;
    sparse_index_base_t tB = SPARSE_INDEX_BASE_ZERO;
    int rowsA = A.GetNumRows();
    int colsA = A.GetNumCols();
    int *rPA = (Integer*) A.GetRowPointer();
    int *cPA = (Integer*) A.GetColPointer();
    double *dPA = A.GetDataPointer();

    sparse_matrix_t B_MKL = NULL;
    int rowsB = B.GetNumRows();
    int colsB = B.GetNumCols();
    int *rPB = (Integer*) B.GetRowPointer();
    int *cPB = (Integer*) B.GetColPointer();
    double *dPB = (double*) B.GetDataPointer();



    /********* Convert our CRS_Matrices into MKL-sparse csr-matrix **********/
    CALL_AND_CHECK_STATUS(mkl_sparse_d_create_csr( &A_MKL, tB, rowsA, colsA, rPA, rPA+1, cPA, dPA ),
        "Error after MKL_SPARSE_D_CREATE_CSR \n");

    CALL_AND_CHECK_STATUS(mkl_sparse_d_create_csr( &B_MKL, tB, rowsB, colsB, rPB, rPB+1, cPB, dPB),
        "Error after MKL_SPARSE_D_CREATE_CSR \n");

    /**************** Perform tempC1 = A^T * B ************************/
    sparse_matrix_t tmpC1 = NULL;
    sparse_operation_t t1 = (version==1)? SPARSE_OPERATION_TRANSPOSE:SPARSE_OPERATION_NON_TRANSPOSE;
    CALL_AND_CHECK_STATUS(mkl_sparse_spmm(t1, A_MKL, B_MKL, &tmpC1),
        "Error after MKL_SPARSE_SPMM \n");

    /**************** Perform C = (A^T * B) * A *********************/
    sparse_matrix_t C_MKL = NULL;
    sparse_matrix_t tmpC_MKL = NULL;
    if(version == 2){
      sparse_operation_t t = SPARSE_OPERATION_TRANSPOSE;
      CALL_AND_CHECK_STATUS(mkl_sparse_convert_csr(tmpC1, t, &tmpC_MKL),
            "Error after MKL_SPARSE_CONVERT_CSR, tempBHT = tempBH^T \n");
      sparse_operation_t t3 = SPARSE_OPERATION_NON_TRANSPOSE;
      CALL_AND_CHECK_STATUS(mkl_sparse_spmm(t3, A_MKL, tmpC_MKL, &C_MKL),"Error after MKL_SPARSE_SPMM \n");
    }else{
      sparse_operation_t t3 = SPARSE_OPERATION_NON_TRANSPOSE;
      CALL_AND_CHECK_STATUS(mkl_sparse_spmm(t3, tmpC1, A_MKL, &C_MKL),"Error after MKL_SPARSE_SPMM \n");
    }

    /**************** Export C in MKL-csr format ********************/
    sparse_index_base_t sbC;
    int rowsC, colsC;
    int *solColsC;
    double *valuesC = NULL;
    int *pCB, *pCE;
    CALL_AND_CHECK_STATUS(mkl_sparse_d_export_csr( C_MKL, &sbC, &rowsC, &colsC, &pCB, &pCE, &solColsC, &valuesC),
        "Error after MKL_SPARSE_D_EXPORT_CSR \n");

    /*********** Set row-, column- and data-pointers *****************/
    rPC.Resize( rowsC + 1, 0);
    //TODO I'm sure there's a method in MKL to get nnz, but I couldn't find it
    UInt nnz = 0;
    for(Integer i = 0; i < rowsC; i++ ){
      nnz += pCE[i] - pCB[i];
      rPC[i + 1] = nnz;
    }
    cPC.Resize(nnz, 0);
    dPC.Resize(nnz, 0.0);
    for(UInt i = 0; i < (UInt)rowsC; ++i){
      for(UInt j = rPC[i]; j < rPC[i+1]; ++j){
        cPC[j] = (UInt)solColsC[j];
        dPC[j] = valuesC[j];
      }
    }

    if(setSparsity) this->SetSparsityPatternData(rPC, cPC, dPC);

    /************ Release handles and deallocate memory *************/
    if( mkl_sparse_destroy( A_MKL ) != SPARSE_STATUS_SUCCESS)
    { printf(" Error after MKL_SPARSE_DESTROY, P \n");fflush(0); }

    if( mkl_sparse_destroy( B_MKL ) != SPARSE_STATUS_SUCCESS)
    { printf(" Error after MKL_SPARSE_DESTROY, Bh \n");fflush(0); }

    if( mkl_sparse_destroy( tmpC1 ) != SPARSE_STATUS_SUCCESS)
    { printf(" Error after MKL_SPARSE_DESTROY, tempBH \n");fflush(0); }

    if(version == 2){
      if( mkl_sparse_destroy( tmpC_MKL ) != SPARSE_STATUS_SUCCESS)
    { printf(" Error after MKL_SPARSE_DESTROY, tmpC_MKL \n");fflush(0); }
    }

    if( mkl_sparse_destroy( C_MKL ) != SPARSE_STATUS_SUCCESS)
    { EXCEPTION("Error in the MKL sparse matrix-matrix multiplication\n"
        "try to increase the size of the coarse system"); }

#else
    EXCEPTION("Compile with USE_MKL = ON in order to use the AMG-framework!")
#endif
  }




// Explicit template instantiation
  template class CRS_Matrix<Double>;
  template class CRS_Matrix<Complex>;
#ifdef MSVC
  template void CRS_Matrix<Double>::Mult(const Vector<Double>&, Vector<Double>&) const;
  template void CRS_Matrix<Complex>::Mult(const Vector<Complex>&, Vector<Complex>&) const;
  template void CRS_Matrix<Double>::MultT(const Vector<Double>&, Vector<Double>&) const;
  template void CRS_Matrix<Complex>::MultT(const Vector<Complex>&, Vector<Complex>&) const;
  template void CRS_Matrix<Double>::MultAdd(const Vector<Double>&, Vector<Double>&) const;
  template void CRS_Matrix<Complex>::MultAdd(const Vector<Complex>&, Vector<Complex>&) const;
  template void CRS_Matrix<Double>::MultTAdd(const Vector<Double>&, Vector<Double>&) const;
  template void CRS_Matrix<Complex>::MultTAdd(const Vector<Complex>&, Vector<Complex>&) const;
  template void CRS_Matrix<Double>::MultSub(const Vector<Double>&, Vector<Double>&) const;
  template void CRS_Matrix<Complex>::MultSub(const Vector<Complex>&, Vector<Complex>&) const;
  template void CRS_Matrix<Double>::MultTSub(const Vector<Double>&, Vector<Double>&) const;
  template void CRS_Matrix<Complex>::MultTSub(const Vector<Complex>&, Vector<Complex>&) const;
#endif

}
