// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

// the following headers are required for Export()
#include <cstdio>


#include "matvec/crs_matrix.hh"
#include "multigrid/prematrix.hh"


// Implementation of methods for the compressed row storage matrix class


namespace OLAS {


  // *********************************
  //   Constructor using CoordFormat
  // *********************************
  template<typename T>
  CRS_Matrix<T>::CRS_Matrix( CoordFormat<T> &sparseMat )
    : diagPtr_( NULL ) {


    // Test, if the matrix is stored in symmetric format.
    // If yes, issue a warning, since we due not expand
    // it to a full matrix currently.
    bool isSymm = sparseMat.GetSymmStorageFlag();
    if ( isSymm == true ) {
      (*warning) << "Input matrix is stored in symmetric format!\n"
		 << "We cannot expand this! Creating a matrix containing "
		 << "only the upper triangle.\n Maybe you should use "
		 << "an SCRS_Matrix?";
      Warning( __FILE__, __LINE__ );
    }

    // Set general matrix pattern and dimension info
    this->nrows_ = sparseMat.GetNrows();
    this->ncols_ = sparseMat.GetNcols();
    this->nnz_   = sparseMat.GetNumEntries();

    // Allocate memory
    NewArray( data_, T, this->nnz_ );
    NewArray( rowPtr_, Integer, this->nrows_ + 1 );
    NewArray( colInd_, Integer, this->nnz_ );

    // Generate auxilliary vector
    Integer *auxVec = NULL;
    NewArray( auxVec, Integer, this->nrows_ );

    // Determine the number of non-zero entries in the different rows
    Integer *entriesPerRow = auxVec;
    for ( Integer k = 1; k <= this->nrows_; k++ ) {
      entriesPerRow[k] = 0;
    }
    for ( Integer i = 1; i <= this->nnz_; i++ ) {
      entriesPerRow[sparseMat.ridx(i)]++;
    }

    // Now prepare the row pointer array
    rowPtr_[1] = 1;
    for ( Integer k = 1; k <= this->nrows_; k++ ) {
      rowPtr_[k+1] = rowPtr_[k] + entriesPerRow[k];
    }

    // Make a copy of values and column indices
    Integer *pos = auxVec;
    int curRow;
    for ( UInt i = 1; i <= (UInt)this->nrows_; i++ ) {
      pos[i] = 0;
    }
    for ( Integer j = 1; j <= this->nnz_; j++ ) {
      curRow = sparseMat.ridx(j);
      data_  [ rowPtr_[curRow] + pos[curRow] ] = sparseMat.val (j);
      colInd_[ rowPtr_[curRow] + pos[curRow] ] = sparseMat.cidx(j);
      pos[ curRow ]++;
    }

    // Destroy auxilliary vector
    DeleteArray( auxVec );

    // We do not expect the input matrix to be sorted in anyway
    currentLayout_ = UNSORTED;

    // By default we generate a matrix in LEX sub-format
    // Note: This will also setup the diagPtr_ array
    NewArray( diagPtr_, Integer, this->nrows_ );
    ChangeLayout( CRS_Matrix<T>::LEX );

  }


  // ********************
  //   Copy Constructor                    
  // ********************
  template<typename T>
  CRS_Matrix<T>::CRS_Matrix( const CRS_Matrix<T> &origMat )
    : diagPtr_( NULL ) {


    // Set basic size informations
    this->nnz_        = origMat.nnz_;
    this->nrows_      = origMat.nrows_;
    this->ncols_      = origMat.ncols_;

    // Allocate memory for internal arrays
    NewArray( colInd_ , Integer, this->nnz_        );
    NewArray( rowPtr_ , Integer, this->nrows_ + 1  );
    NewArray( diagPtr_, Integer, this->nrows_      );
    NewArray( data_   , T      , this->nnz_        );

    // Copy information
    for (int i = 1; i <= this->nnz_; i++ ) {
      data_[i]   = origMat.data_[i];
      colInd_[i] = origMat.colInd_[i];
    }

    for (int i = 1; i <= this->nrows_; i++ ) {
      rowPtr_[i]  = origMat.rowPtr_[i];
      diagPtr_[i] = origMat.diagPtr_[i];
    }
    rowPtr_[this->nrows_+1] = origMat.rowPtr_[this->nrows_+1];

    // Copy layout flag
    currentLayout_ = origMat.currentLayout_;
  }


  // ******************************
  //   Construct from a prematrix 
  // ******************************
  template <typename T>
  CRS_Matrix<T>::CRS_Matrix( Integer nrows, Integer ncols, 
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
    if ( nrows != premat.GetNrows() ) {
      Warning( "CRS_Matrix: PreMatrix has mismatched size ", __FILE__,
	       __LINE__ );
    }

    rowPtr_[1] = 1;

    Integer rowlength;
    Integer position = 1;

    for( Integer i = 1; i <= this->nrows_; i++ ) {

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
    NewArray( diagPtr_, Integer, this->nrows_ );
    FindDiagonalEntries();
  }


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


    Integer rs;
    Integer *index;

    rowPtr_[1] = 1;

#ifdef DEBUG_CRS_MATRIX
    (*debug) << "The Matrix Graph has size " << graph.GetSize() << std::endl;
    graph.Print( *debug );
#endif

    // loop over all rows
    for ( Integer i = 1; i <= this->nrows_; i++ ) {

      // get neighbours (indices and numbers) of i-th node
      index = graph.GetGraphRow(i);
      rs = graph.GetRowSize(i);

      // set next row pointer (beginning of next row)
      rowPtr_[i+1] = rowPtr_[i] + rs;

      // Assume that there is no diagonal entry (and correct
      // it below, if we find one)
      diagPtr_[i] = 0;

      // set column indices of non-zero entries in this row
      for ( Integer k = 1; k <= rs; k++ ) {
        colInd_[rowPtr_[i]+k-1] = index[k];

        // look out for diagonal entry
        if ( index[k] == i ) {
          diagPtr_[i] = rowPtr_[i]+k-1;
        }
      }
    }

    // setting the sparsity pattern using the information provided by
    // the BaseGraph leads to a matrix in LEX layout
    currentLayout_ = CRS_Matrix<T>::LEX;

  }

  // *************************
  //   CopySparstiyPattern
  // *************************
  template<typename T>
  void CRS_Matrix<T>::CopySparsityPattern( StdMatrix &mat ) const {

    (*error) << "CRS_Matrix::CopySparsityPattern: Not implemented";
    Error( __FILE__, __LINE__ );
   
  }

  // ********************************
  //   Matrix-Vector multiplication
  // ********************************
  template<typename T>
  inline void CRS_Matrix<T>::Mult( const Vector<T> &mvec,
                                   Vector<T> &rvec ) const {


    PROFILE( "CRS_Matrix::Mult", 2*this->nnz_*BlockSize<T>::size *
             BlockSize<T>::size );

    Integer i, j, rs, k;
    T sum;

#pragma omp parallel for private(sum,k,rs,j)
    for ( i = 1; i <= this->nrows_; i++ ) {
      sum = 0;
      k  = rowPtr_[i];
      rs = rowPtr_[i+1] - k;                                
      for ( j = 1; j <= rs; j++ ) {
        sum += data_[k] * mvec[colInd_[k]];
        k++;
      }
      rvec[i] = sum;
    }
  }


  // ***********
  //   MultAdd
  // ***********
  template<typename T>
  inline void CRS_Matrix<T>::MultAdd( const Vector<T> & mvec,
                                      Vector<T> & rvec ) const {


    Integer i, j, rs, k;
    T sum;

    // since this function is mainly used for the off-diagonal entries
    // in the parallel case we disable profiling to avoid inconsistent data
#ifndef PARALLEL
    PROFILE( "CRS_Matrix::MultAdd", (2*this->nnz_*BlockSize<T>::size +
                                       this->nrows_)*BlockSize<T>::size);
#endif

#pragma omp parallel for private(sum,k,rs,j)
    for ( i = 1; i <= this->nrows_; i++ ) {
      sum = 0;  
      k  = rowPtr_[i];
      rs = rowPtr_[i+1]-k;                                
      for ( j = 1; j <= rs; j++ ) {
        sum += data_[k] * mvec[colInd_[k]];
        k++;
      }
      rvec[i] += sum;
    }
  }


  // ***********
  //   MultSub
  // ***********
  template<typename T>
  inline void CRS_Matrix<T>::MultSub( const Vector<T> &mvec,
                                      Vector<T> &rvec ) const {


    Integer i, j, rs, k;
    T sum;

    // since this function is mainly used for the off-diagonal entries
    // in the parallel case we disable profiling to avoid inconsistent data
#ifndef PARALLEL
    PROFILE( "CRS_Matrix::MultAdd", (2*this->nnz_*BlockSize<T>::size +
                                       this->nrows_)*BlockSize<T>::size);
#endif

#pragma omp parallel for private(sum,k,rs,j)
    for ( i = 1; i <= this->nrows_; i++ ) {
      sum = 0;  
      k  = rowPtr_[i];
      rs = rowPtr_[i+1]-k;                                
      for ( j = 1; j <= rs; j++ ) {
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


    PROFILE("CRS_Matrix::CompRes",
            (2*this->nnz_*BlockSize<T>::size + this->nrows_)
            *BlockSize<T>::size);

    Integer i, j, rs, k;
    T sum;

#pragma omp parallel for private(sum,k,rs,j)
    for ( i = 1; i <= this->nrows_; i++ ) {
      sum = 0;
      k  = rowPtr_[i];
      rs = rowPtr_[i+1]-k;
      for ( j = 1; j <= rs; j++ ) {
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


    Integer i, j, rs, k;

    // set result vector to zero
    rvec.Init();

    // loop over matrix rows
    for ( i = 1; i <= this->nrows_; i++ ) {

      // get row start and row size
      k  = rowPtr_[i];
      rs = rowPtr_[i+1] - k;

      // loop over this row
      for ( j = 1; j <= rs; j++ ) {
        rvec[colInd_[k]] += data_[k] * mvec[i];
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

    Integer i, j, end;

    // loop over matrix rows
    for( i = 1; i <= this->nrows_; i++ ) {

      // get end of row
      end = rowPtr_[i+1];

      // loop over this row
      for( j = rowPtr_[i]; j < end; j++ ) {
        rvec[colInd_[j]] += this->data_[j] * mvec[i];
      }
    }
  }

  template<typename T>
  void CRS_Matrix<T>::Transpose(int* col_ptr, int* row_ptr, T* data_ptr) const
  {
     // first set the explicit row-index temporaly along the colInd_.
     Integer* tmpInd;
     NewArray(tmpInd, int, this->nnz_);

     for(int i = 1; i <= this->nrows_; i++) 
     {
        for(int k = rowPtr_[i]; k < rowPtr_[i+1]; k++) 
        {
           int j = colInd_[k];
           tmpInd[k] = i;
        }
     }

     // now set up the column structure
     int p = 1;
     for(int j = 1; j <= this->ncols_; j++)
     {
        // store the column indexx
        col_ptr[j] = p;
        // search in colInd_ for this column
        for(int k = 1; k <= this->nnz_; k++)
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
     col_ptr[this->ncols_+1] = this->nnz_+1;    

     DeleteArray(tmpInd);
  }


  // ****************************
  //   Print matrix to a stream
  // ****************************
  template<typename T>
  void CRS_Matrix<T>::Print( std::ostream& os ) const {


    Integer i, j, k;

    for ( i = 1; i <= this->nrows_; i++ ) {
      for ( k = rowPtr_[i]; k < rowPtr_[i+1]; k++ ) {
        j = colInd_[k];
        os << std::scientific << i << "\t" << j << "\t" << data_[k] << '\n';
      }
    }

    os << std::endl;

#ifdef DEBUG_CRS_MATRIX
    os << " diagPtr_: ";
    for ( i = 1; i <= this->nrows_; i++ ) {
      os << diagPtr_[i] << " ";
    }    
    os << '\n' << std::endl;
#endif

  }


  // *************************
  //   Export matrix to file
  // *************************
  template<typename T>
  void CRS_Matrix<T>::Export( const Char *fname,
                              const Char *comment ) const{


    // Obtain blocksize of matrix
    Integer dof = this->GetBlockSize();

    // Open output file and check for errors
    FILE *fp = fopen( fname, "w" );
    if ( fp == NULL ) {
      Char *errmsg = NULL;
      NewArray( errmsg, Char, strlen(fname)+40 );
      sprintf( errmsg, "Cannot open file %s for writing!", fname );
      Error( errmsg, __FILE__, __LINE__ );
    }

    // ---------------------
    //   Write file header
    // ---------------------

    // Matrix Market Format Specification
    if ( this->GetEntryType() == DOUBLE ) {
      fprintf( fp, "%%%%MatrixMarket matrix coordinate real general\n" );
    }
    else if ( this->GetEntryType() == COMPLEX ) {
      fprintf( fp, "%%%%MatrixMarket matrix coordinate complex general\n" );
    }
    else {
      Error( "Expected either DOUBLE or COMPLEX as matrix entry type",
             __FILE__, __LINE__ );
    }

    // User-supplied private comment
    if ( comment != NULL ) {
      fprintf( fp, "%%\n%% %s\n%%\n", comment );
    }
    else {
      fprintf( fp, "%%\n%% Matrix exported by OLAS\n%%\n" );
    }

    // Information on number of rows, columns and entries
    fprintf( fp, "%d\t%d\t%d\n", this->nrows_ * dof, this->ncols_ * dof,
             this->nnz_ * dof * dof );


    // ------------------------
    //   Write matrix entries
    // ------------------------

    Integer i, j, ib, jb, nblocks;

    // loop over all block rows
    for ( i = 1; i <= this->nrows_; i++ ) {

      // get number of blocks in i-th row
      nblocks = rowPtr_[i+1] - rowPtr_[i];

      // loop over all blocks in this row
      for ( j = 1; j <= nblocks; j++ ) {

        // loop over block entries
        for ( ib = 0; ib < dof; ib++ ) {
          for ( jb = 0; jb < dof; jb++ ) {

            // store row and column index
            fprintf( fp, "%6d\t%6d\t", (i-1) * dof + ib + 1,
                     ( colInd_[rowPtr_[i]+j-1] - 1 ) * dof + jb + 1);

            // store non-zero entry
            opType<T>::ExportEntry( this->data_[rowPtr_[i]+j-1],
                                          ib, jb, fp );
            fprintf( fp, "\n" );
          }
        }
      }
    }

    // close output file
    if ( fclose( fp ) == EOF ) {
      Char *errmsg = NULL;
      NewArray( errmsg, Char, strlen(fname)+40 );
      sprintf( errmsg, "Could not close file %s after writing!", fname );
      Warning( errmsg, __FILE__, __LINE__ );
    }  
  }


  // *******************************
  //   Add value to a matrix entry
  // *******************************
  template<typename T>
  void CRS_Matrix<T>::AddToMatrixEntry( Integer i, Integer j, T& v ){


    bool found = false;

    // Try to determine index for matrix entry at position (i,j)
    // in the data_ array and add value
    for ( Integer k = rowPtr_[i]; k < rowPtr_[i+1]; k++ ) {
      if ( colInd_[k] == j ) {
        found = true;
        data_[k] += v;
        break;
      }
    }

    // Its a serious error, if the index pair was not found
    if ( found == false ) {
      Print( std::cerr );
      (*error) << "AddToMatrixEntry: Index pair = (" << i << " , "
               << j << ") not found\n";
      Error( __FILE__, __LINE__ );
    }
  }


  // *****************************
  //   Get specific matrix entry
  // *****************************
  template<typename T>
  void CRS_Matrix<T>::GetMatrixEntry( Integer i, Integer j, 
                                      T &v ) const {

    bool found = false;

    // Try to determine index for matrix entry at position (i,j)
    // in the data_ array and get the entry
    for ( Integer k = rowPtr_[i]; k < rowPtr_[i+1]; k++ ) {
      if ( colInd_[k] == j ) {
        found = true;
        v = data_[k];
        break;
      }
    }

    // Its a serious error, if the index pair was not found
    if ( found == false ) {
      Print( std::cerr );
      (*error) << "GetMatrixEntry: Index pair = (" << i << " , "
               << j << ") not found\n";
      Error( __FILE__, __LINE__ );
    }
  }

  // *****************************
  //   Set specific matrix entry
  // *****************************
  template<typename T>
  void CRS_Matrix<T>::SetMatrixEntry( Integer i, Integer j, T &v,
                                      bool setCounterPart ) {

    bool found = false;
    
    // Try to determine index for matrix entry at position (i,j)
    // in the data_ array and set the value
    for ( Integer k = rowPtr_[i]; k < rowPtr_[i+1]; k++ ) {
      if ( colInd_[k] == j ) {
        found = true;
        data_[k] = v;
        break;
      }
    }
    
    // Its a serious error, if the index pair was not found
    if ( found == false ) {
      Print( std::cerr );
      (*error) << "SetMatrixEntry: Index pair = (" << i << " , "
               << j << ") not found\n";
      Error( __FILE__, __LINE__ );
    }

    // If setCounterPart is true and the index-pair denotes a non-diagonal
    // position, we do also insert the counterpart
    if ( i != j && setCounterPart == TRUE ) {
      found = false;
      for ( Integer k = rowPtr_[j]; k < rowPtr_[j+1]; k++ ) {
        if ( colInd_[k] == i ) {
          found = true;
          data_[k] = v;
          break;
        }
      }
      
      // Its a serious error, if the index pair was not found
      if ( found == false ) {
        Print( std::cerr );
        (*error) << "SetMatrixEntry: Index pair = (" << j << " , "
                 << i << ") not found\n";
        Error( __FILE__, __LINE__ );
      }
    }
  }

  // **********
  //   Insert
  // **********
  template<typename T>
  void CRS_Matrix<T>::Insert(Integer row, Integer col,  Integer pos){
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


  // ***********************
  //  Forced Instantiation
  // ***********************
  template <typename T>
  void CRS_Matrix<T>::InstantiatePublicMethods() {


    Error( "This function should never be called", __FILE__, __LINE__ );
    
    Integer dummyInt = 0;
    Integer *dummyPInt = NULL;
    BaseGraph dummyGraph( 0, 0, NOREORDERING );
    bool dummyBool = false;
    MatrixStorageType dummyType = NOSTORAGETYPE;
    Vector<T> dummyVec;
    T  dummyEntry;
    T* dummyPVal = NULL;
    double dummyVal;

    // Constructors
    PreMatrix<T> preMat;
    CoordFormat<T>  coordMat( 0, 0, 0, true );
    // CoordFormat<Double>  coordMatD( 0, 0, 0, true );
    // CoordFormat<Complex> coordMatC( 0, 0, 0, true );

    CRS_Matrix<T> dummyMat( 10, 10, 10 );
    CRS_Matrix<T> dummyPMat( 10, 10, preMat );
    CRS_Matrix<T>  dummyDMat( coordMat );
    // CRS_Matrix<Double>  dummyDMat( coordMatD );
    // CRS_Matrix<Complex> dummyCMat( coordMatC );

    dummyMat.Init();

    // Layout stuff
    subFormat myLayout = CRS_Matrix<T>::UNSORTED;
    ChangeLayout( myLayout );
    SetCurrentLayout( myLayout );
    myLayout = GetCurrentLayout();

    SetSparsityPattern( dummyGraph );
    SetSize( dummyInt, dummyInt, dummyInt );
    Print( std::cout );
    Export( "" );
    SetRowSize( dummyInt, dummyInt );
    Insert( dummyInt, dummyInt, dummyInt );
    AddToMatrixEntry( dummyInt, dummyInt, dummyEntry );

    dummyType = GetStorageType();

    Mult    ( dummyVec, dummyVec );
    MultAdd ( dummyVec, dummyVec );
    MultSub ( dummyVec, dummyVec );
    CompRes ( dummyVec, dummyVec, dummyVec );
    MultT   ( dummyVec, dummyVec );
    MultTAdd( dummyVec, dummyVec );
    Transpose(NULL, NULL, NULL);

    Scale( 1.0 );
    Add( 1.0, dummyMat );

    T v;
    GetDiagEntry( 14, v );
    SetDiagEntry( 27, v );

    // AccountForPenalty( dummyVec, dummyVec );

    dummyType = GetStorageType();
    dummyInt  = GetRowSize( dummyInt );
    dummyInt  = this->GetBlockSize();
    dummyVal  = GetMaxDiag();

    dummyPInt = GetRowPointer();
    dummyPInt = GetColPointer();
    dummyPInt = GetDiagPointer();
    dummyPVal = GetDataPointer();

    MatrixEntryType met = this->GetEntryType();
    met = NOENTRYTYPE;

    const Integer* cInt1 = GetRowPointer();
    const Integer* cInt2 = GetColPointer();

    // use the stuff to satisfy gcc
    std::cout << cInt1 << cInt2;
    dummyBool = dummyBool == true ? false : true;

    // For deprecated methods
    SortConformingToLayout();

  }


  // *********
  //   Scale
  // *********
  template<typename T>
  void CRS_Matrix<T>::Scale( Double factor ) {
    for ( Integer i = 1; i <= this->nnz_; i++ ) {
      data_[i] *= factor;
    }
  }


  // **************
  //   GetMaxDiag
  // **************
  template<typename T>
  Double CRS_Matrix<T>::GetMaxDiag() const {


    double maxDiag = 0;
    double current = 0;
    Integer i;

    for ( i = 1; i <= this->nrows_; i++ ) {

      // use an opType to ensure that tiny matrices
      // are treated correctly
      current = opType<T>::MaxDiag( data_[ diagPtr_[i] ] );
      maxDiag = maxDiag > current ? maxDiag : current;
    }

    return maxDiag;
  }


  // ************************
  //   Add (another matrix)
  // ************************
  template<typename T>
  void CRS_Matrix<T>::Add( const Double alpha, const StdMatrix& mat ) {

    TRY_CAST {

      // Down-cast input matrix
      ConstRefCast( mat, CRS_Matrix<T>, crsMat );

      // Obtain pointer to data array of other matrix
      const T *data = crsMat.GetDataPointer();

      // We now assume that the matrices have matching
      // dimensions and sparsity patterns
      for ( Integer i = 1; i <= this->nnz_; i++ ) {
        data_[i] += alpha * data[i];
      }

    } CATCH_CAST;

  }

  template<>
  void CRS_Matrix<Complex>::Add( const Double alpha, const StdMatrix& mat ) {
    
    // Check for entry type of mat
    MatrixEntryType eType = mat.GetEntryType();
    
    if( eType == DOUBLE ) { 
      TRY_CAST {
        
        // Down-cast input matrix
        ConstRefCast( mat, CRS_Matrix<Double>, crsMat );
        
        // Obtain pointer to data array of other matrix
        const Double *data = crsMat.GetDataPointer();
        
        // We now assume that the matrices have matching
        // dimensions and sparsity patterns
        for ( Integer i = 1; i <= this->nnz_; i++ ) {
          data_[i] += alpha * Complex(data[i], 0.0 );
        }
        
      } CATCH_CAST;

    } else {
      TRY_CAST {
        
        // Down-cast input matrix
        ConstRefCast( mat, CRS_Matrix<Complex>, crsMat );
        
        // Obtain pointer to data array of other matrix
        const Complex * data = crsMat.GetDataPointer();
        
        // We now assume that the matrices have matching
        // dimensions and sparsity patterns
        for ( Integer i = 1; i <= this->nnz_; i++ ) {
          data_[i] += alpha * data[i];
        }
        
      } CATCH_CAST;
    }

  }


  // *************
  //   QuickSort
  // *************
  template <typename T>
  void CRS_Matrix<T>::QuickSort( Integer *const cols, T *const data,
                                 const Integer length ) {

    if( length > 1 ) {

      const Integer last         = length-1;
      const Integer splitter     = cols[last];
      const T splitterData = data[last];

      Integer tcol, i = 0;
      T tdat;

      // split array
      for( Integer j = 0; j < last; j++ ) {

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

      // Store current layout for preparing final report
      subFormat oldLayout = currentLayout_;

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

          for( Integer i = 1; i <= this->GetNrows(); i++ ) {

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
          for( Integer i = 1; i <= this->GetNrows(); i++ ) {

            // search the diagonal entry
            for( Integer ij = rowPtr_[i]; ij < rowPtr_[i+1]; ij++ ) {

              // place the diagonal entry at leading position
              if( colInd_[ij] == i ) {
                const Integer idiag = rowPtr_[i];

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
        (*error) << " Cannot convert CRS_Matrix to UNSORTED sub-format!";
        Error( __FILE__, __LINE__ );
      }


      // ---------
      //  Default
      // ---------
      else {
        (*error) << "Congratulations! You have found a missing case "
                 << "implementation!";
        Error( __FILE__, __LINE__ );
      }

      // Report to standard report file
      (*cla) << " CRS_Matrix: Changed sub-format from '"
             << Enum2String( oldLayout )
             << "' to '"
             << Enum2String( currentLayout_ ) << "'" << std::endl;
    }
  }


  // ******************
  //   Re-size matrix 
  // ******************
  template<typename T>
  void CRS_Matrix<T>::SetSize( Integer nrows, Integer ncols, Integer nnz ) {

    
#ifdef DEBUG_CRS_MATRIX
    if ( nrows < 0 || nnz < 0 ) {
      Error( "invalid dimensions", __FILE__, __LINE__ );
    }
#endif

    this->ncols_ = ncols;
    if ( this->nrows_ != nrows ) {
      this->nrows_ = nrows;
      DeleteArray( rowPtr_  );
      DeleteArray( diagPtr_ );
      NewArray( rowPtr_ , Integer, this->nrows_ + 1 );
      NewArray( diagPtr_, Integer, this->nrows_     );
      rowPtr_[1]  = 1;
      diagPtr_[1] = 1;
      currentLayout_ = CRS_Matrix<T>::UNSORTED;
    }

    if ( this->nnz_ != nnz ) {
      this->nnz_ = nnz;
      DeleteArray( colInd_ );
      DeleteArray( data_ );
      NewArray( colInd_, Integer, this->nnz_ );
      NewArray( data_, T, this->nnz_ );
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
      for ( Integer i = 1; i <= this->nrows_; i++ ) {
        diagPtr_[i] = 0;
        if ( colInd_[rowPtr_[i]] == i ) {
          diagPtr_[i] = rowPtr_[i];
        }
      }
    }
    
    // For other formats we must search the diagonal entries
    else {

      // loop over all rows
      for ( Integer i = 1; i <= this->nrows_; i++ ) {

        // assume that there is no diagonal entry (and correct
        // it below, if we find one)
        diagPtr_[i] = 0;
      
        // loop over column
        for ( Integer k = rowPtr_[i]; k < rowPtr_[i+1]; k++ ) {

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
    for ( Integer i = 1; i <= this->nrows_; i++ ) {

      // Check if row has a diagonal entry (otherwise this might
      // not be sensible)
      if ( diagPtr_[i] == 0 ) {
        (*error) << "CRS_Matrix<T>::SortLex2LexDiagFirst: There is no "
                 << "diagonal entry in row " << i;
        Error( __FILE__, __LINE__ );
      }

      // Store diagonal entry
      auxVal = data_[ diagPtr_[i] ];

      // shift left part of row to the right
      for ( Integer k = diagPtr_[i]; k > rowPtr_[i]; k-- ) {
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
    Integer auxInd;

    // loop over all rows
    for ( Integer i = 1; i <= this->nrows_; i++ ) {

      // Only, if there is a diagonal entry, do we have something to do
      if ( diagPtr_[i] == rowPtr_[i] ) {

        for ( Integer j = rowPtr_[i] + 1; j < rowPtr_[i+1]; j++ ) {

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

}
