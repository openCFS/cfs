// the following headers are required for Export()
#include <cstdio>
#include <sstream>


#include "CRS_Matrix.hh"
#include "opdefs.hh"


// Implementation of methods for the compressed row storage matrix class


namespace CoupledField {


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
      WARN("Input matrix is stored in symmetric format!\n"
		       << "We cannot expand this! Creating a matrix containing "
		       << "only the upper triangle.\n Maybe you should use "
		       << "an SCRS_Matrix?");
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


    UInt rs;
    UInt *index;

    rowPtr_[0] = 0;

#ifdef DEBUG_CRS_MATRIX
    (*debug) << "The Matrix Graph has size " << graph.GetSize() << std::endl;
    graph.Print( *debug );
#endif

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

  // *************************
  //   CopySparstiyPattern
  // *************************
  template<typename T>
  void CRS_Matrix<T>::CopySparsityPattern( StdMatrix &mat ) const {

    EXCEPTION( "CRS_Matrix::CopySparsityPattern: Not implemented" );

  }

  // ********************************
  //   Matrix-Vector multiplication
  // ********************************
  template<typename T>
  inline void CRS_Matrix<T>::Mult( const Vector<T> &mvec,
                                   Vector<T> &rvec ) const {
    UInt i, j, rs, k;
    T sum;

#pragma omp parallel for private(sum,k,rs,j)
    for ( i = 0; i < this->nrows_; i++ ) {
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


  // ***********
  //   MultAdd
  // ***********
  template<typename T>
  inline void CRS_Matrix<T>::MultAdd( const Vector<T> & mvec,
                                      Vector<T> & rvec ) const {


    UInt i, j, rs, k;
    T sum;


#pragma omp parallel for private(sum,k,rs,j)
    for ( i = 0; i < this->nrows_; i++ ) {
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


  // ***********
  //   MultSub
  // ***********
  template<typename T>
  inline void CRS_Matrix<T>::MultSub( const Vector<T> &mvec,
                                      Vector<T> &rvec ) const {


    UInt i, j, rs, k;
    T sum;

    // since this function is mainly used for the off-diagonal entries
    // in the parallel case we disable profiling to avoid inconsistent data

#pragma omp parallel for private(sum,k,rs,j)
    for ( i = 0; i < this->nrows_; i++ ) {
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

    UInt i, j, rs, k;
    T sum;

#pragma omp parallel for private(sum,k,rs,j)
    for ( i = 0; i < this->nrows_; i++ ) {
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
        // store the column indexx
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


  // *************************
  //   Export matrix to file
  // *************************
  template<typename T>
  void CRS_Matrix<T>::Export( const char *fname,
                              const char *comment ) const{


    // Open output file and check for errors
    FILE *fp = fopen( fname, "w" );
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
    if ( comment != NULL ) {
      fprintf( fp, "%%\n%% %s\n%%\n", comment );
    }
    else {
      fprintf( fp, "%%\n%% Matrix exported by OLAS\n%%\n" );
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
          data_[k] += v;
          break;
        }
      }
      break;
    case LEX_DIAG_FIRST:
      if(colInd_[l] == j){ // the diagonal exists
        if(i == j){ // we want the diagonal
          data_[j] += v;
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
          data_[k] += v;
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
          data_[k] = v;
          break;
        }
      }
      break;
    case LEX_DIAG_FIRST:
      if(colInd_[l] == j){ // the diagonal exists
        if(i == j){ // we want the diagonal
          data_[j] = v;
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
          data_[k] = v;
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
  template<typename T>
  void CRS_Matrix<T>::Scale( Double factor ) {
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


  // **************
  //   GetMaxDiag
  // **************
  template<typename T>
  Double CRS_Matrix<T>::GetMaxDiag() const {


    double maxDiag = 0;
    double current = 0;
    UInt i;

    for ( i = 0; i < this->nrows_; i++ ) {

      // use an opType to ensure that tiny matrices
      // are treated correctly
      current = OpType<T>::MaxDiag( data_[ diagPtr_[i] ] );
      maxDiag = maxDiag > current ? maxDiag : current;
    }

    return maxDiag;
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
      register UInt k, rs;
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
        register UInt k, rs;
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
        std::set<UInt>::const_iterator colIt;
        register UInt k, rs;
        UInt i, j;

        // loop over all rows
        for ( i = 0; i < this->nrows_; i++ ) {
          k = rowPtr_[i];
          rs = rowPtr_[i+1] - k;

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
          } // cols
        } // rows
      } // if column set is nonzero
    } // if one set is nonzero
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
  void CRS_Matrix<T>::SetSize( UInt nrows, UInt ncols, UInt nnz ) {


#ifdef DEBUG_CRS_MATRIX
    if ( nrows < 0 || nnz < 0 ) {
      EXCEPTION( "invalid dimensions" );
    }
#endif

    this->ncols_ = ncols;
    if ( this->nrows_ != nrows ) {
      this->nrows_ = nrows;
      delete [] ( rowPtr_  );  rowPtr_   = NULL;
      delete [] ( diagPtr_ );  diagPtr_  = NULL;
      NEWARRAY( rowPtr_ , UInt, this->nrows_ + 1 );
      NEWARRAY( diagPtr_, UInt, this->nrows_     );
      rowPtr_[0]  = 0;
      diagPtr_[0] = 0;
      currentLayout_ = CRS_Matrix<T>::UNSORTED;
    }

    if ( this->nnz_ != nnz ) {
      this->nnz_ = nnz;
      delete [] ( colInd_ );  colInd_  = NULL;
      delete [] ( data_ );  data_  = NULL;
      NEWARRAY( colInd_, UInt, this->nnz_ );
      NEWARRAY( data_, T, this->nnz_ );
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
      if ( diagPtr_[i] == 0 ) {
        EXCEPTION( "CRS_Matrix<T>::SortLex2LexDiagFirst: There is no "
                 << "diagonal entry in row " << i );
      }

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

// Explicit template instantiation
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template class CRS_Matrix<Double>;
  template class CRS_Matrix<Complex>;
#endif

}
