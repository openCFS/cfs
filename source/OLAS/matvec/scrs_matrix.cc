// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <cassert>

#include "matvec/scrs_matrix.hh"
// Implementation of methods for the symmetric compressed row storage SCRS
// matrix class

namespace OLAS {


  // ********************
  //   Copy Constructor                    
  // ********************

  template<typename T>
  SCRS_Matrix<T>::SCRS_Matrix( const SCRS_Matrix<T> &origMat ) {


    // Set basic size informations
    numEntries_ = origMat.numEntries_;
    this->nnz_        = origMat.nnz_;
    this->nrows_      = origMat.nrows_;
    this->ncols_      = origMat.ncols_;

    // ---------------------
    // Pattern is not shared
    // ---------------------

    if ( origMat.patternPool_ == NULL ) {

      // Allocate memory for internal arrays
      NewArray( colInd_, Integer, numEntries_ );
      NewArray( rowPtr_, Integer, this->nrows_ + 1 );
      NewArray( data_, T, numEntries_ );

      // Copy information
      for ( UInt i = 1; i <= (UInt)numEntries_; i++ ) {
        data_[i]   = origMat.data_[i];
        colInd_[i] = origMat.colInd_[i];
      }

      for ( UInt i = 1; i <= (UInt)this->nrows_ + 1; i++ ) {
        rowPtr_[i] = origMat.rowPtr_[i];
      }
    }

    // -----------------
    // Pattern is shared
    // -----------------
    else {

      // Get hold of pool, pattern identifier and patttern
      patternPool_ = origMat.patternPool_;
      patternID_   = origMat.patternID_;
      colInd_      = origMat.colInd_;
      rowPtr_      = origMat.rowPtr_;

      // Generate copy of data array
      NewArray( data_, T, numEntries_ );
      for ( UInt i = 1; i <= (UInt)numEntries_; i++ ) {
        data_[i]   = origMat.data_[i];
      }
    }


  }


  // *************************
  //   Alternate Constructor                    
  // *************************

  // For all cases, but RMatrix1
  template<typename T>
  SCRS_Matrix<T>::SCRS_Matrix( CoordFormat<T> &sparseMat ) {
    Error( "Alternate Constructor for SCRS_Matrix only works for RMatrix1",
           __FILE__, __LINE__ );
  }


  // Specialisation for the RMatrix1 case
  template<>
  SCRS_Matrix<double>::SCRS_Matrix( CoordFormat<double> &sparseMat ) {


    // Set general matrix pattern and dimension info
    this->nrows_ = sparseMat.GetNrows();
    this->ncols_ = sparseMat.GetNcols();
    this->nnz_   = sparseMat.GetNNZ();
    numEntries_ = (this->nnz_ + this->nrows_) / 2;
    UInt storedEntries = sparseMat.GetNumEntries();

    // Test that matrix is square
    if ( this->nrows_ != this->ncols_ ) {
      (*error) << "The input matrix is of dimension " << this->nrows_ << " x "
	       << this->ncols_ << ", but should be square!";
      Error( __FILE__, __LINE__ );
    }

    // Make a simply consistency check
    if ( sparseMat.GetSymmStorageFlag() == true ) {

      // In this case the number of entries stored in
      // sparseMat should be identical to the number
      // of entries we are going to store
      if ( storedEntries != (UInt)numEntries_ ) {
	(*error) << "Inconsistency detected in input matrix in CoordFormat! "
		 << "Matrix claims to be in symmetric storage format, but "
		 << "fails test on number of entries! "
                 << "storedEntries = " << storedEntries << " , "
                 << "numEntries_ = " << numEntries_;
	Error( __FILE__, __LINE__ );
      }
    }

    else {

      // In this case the number of stored entries should be identical
      // to the number of non-zero matrix entries
      if ( storedEntries != (UInt)this->nnz_ ) {
	(*error) << "Inconsistency detected in input matrix in CoordFormat!\n"
		 << "Matrix claims to be in general storage format, but "
		 << "fails test on number of entries!";
	Error( __FILE__, __LINE__ );
      }
    }

    // Allocate memory
    NewArray( rowPtr_, Integer, this->nrows_ + 1  );
    NewArray( colInd_, Integer, numEntries_ );
    NewArray( data_  , Double , numEntries_ );

    // Sort input matrix, such that it conforms to our layout requirements
    sparseMat.Sort( CoordFormat<double>::PURELEX );

    // Generate auxilliary vector
    Integer *auxVec = NULL;
    NewArray( auxVec, Integer, this->nrows_ );

    // Determine the number of non-zero entries in the different rows.
    // In doing this we ignore entries in the lower triangle.
    Integer *entriesPerRow = auxVec;
    for ( Integer k = 1; k <= this->nrows_; k++ ) {
      entriesPerRow[k] = 0;
    }
    for ( UInt i = 1; i <= storedEntries; i++ ) {
      if ( sparseMat.cidx(i) >= sparseMat.ridx(i) ) {
	entriesPerRow[sparseMat.ridx(i)]++;
      }
    }

    // Now prepare the row pointer array
    rowPtr_[1] = 1;
    for ( Integer k = 1; k <= this->nrows_; k++ ) {
      rowPtr_[k+1] = rowPtr_[k] + entriesPerRow[k];
    }

    // Make a copy of values and column indices
    Integer *pos = auxVec;
    UInt curRow;
    for ( UInt i = 1; i <= (UInt)this->nrows_; i++ ) {
      pos[i] = 0;
    }
    for ( UInt j = 1; j <= storedEntries; j++ ) {
      curRow = sparseMat.ridx(j);
      if ( curRow <= sparseMat.cidx(j) ) {
	data_  [ rowPtr_[curRow] + pos[curRow] ] = sparseMat.val (j);
	colInd_[ rowPtr_[curRow] + pos[curRow] ] = sparseMat.cidx(j);
	pos[ curRow ]++;
      }
    }

    // Destroy auxilliary vector
    DeleteArray( auxVec );

    // Set pattern pool pointer to NULL, since we allocated pattern
    // ourselves
    patternPool_ = NULL;
    patternID_ = NO_PATTERN_ID;
  }


  // Specialisation for the CMatrix1 case
  template<>
  SCRS_Matrix<Complex>::SCRS_Matrix( CoordFormat<Complex> &sparseMat ) {


    // Set general matrix pattern and dimension info
    this->nrows_ = sparseMat.GetNrows();
    this->ncols_ = sparseMat.GetNcols();
    this->nnz_   = sparseMat.GetNNZ();
    numEntries_ = (this->nnz_ + this->nrows_) / 2;
    UInt storedEntries = sparseMat.GetNumEntries();

    // Test that matrix is square
    if ( this->nrows_ != this->ncols_ ) {
      (*error) << "The input matrix is of dimension " << this->nrows_ << " x "
	       << this->ncols_ << ", but should be square!";
      Error( __FILE__, __LINE__ );
    }

    // Make a simply consistency check
    if ( sparseMat.GetSymmStorageFlag() == true ) {

      // In this case the number of entries stored in
      // sparseMat should be identical to the number
      // of entries we are going to store
      if ( storedEntries != (UInt)numEntries_ ) {
	(*error) << "Inconsistency detected in input matrix in CoordFormat!\n"
		 << "Matrix claims to be in symmetric storage format, but "
		 << "fails test on number of entries!";
	Error( __FILE__, __LINE__ );
      }
    }

    else {

      // In this case the number of stored entries should be identical
      // to the number of non-zero matrix entries
      if ( storedEntries != (UInt)this->nnz_ ) {
	(*error) << "Inconsistency detected in input matrix in CoordFormat!\n"
		 << "Matrix claims to be in general storage format, but "
		 << "fails test on number of entries!";
	Error( __FILE__, __LINE__ );
      }
    }

    // Allocate memory
    NewArray( rowPtr_, Integer, this->nrows_ + 1  );
    NewArray( colInd_, Integer, numEntries_ );
    NewArray( data_  , Complex, numEntries_ );

    // Sort input matrix, such that it conforms to our layout requirements
    sparseMat.Sort( CoordFormat<Complex>::PURELEX );

    // Generate auxilliary vector
    Integer *auxVec = NULL;
    NewArray( auxVec, Integer, this->nrows_ );

    // Determine the number of non-zero entries in the different rows.
    // In doing this we ignore entries in the lower triangle.
    Integer *entriesPerRow = auxVec;
    for ( Integer k = 1; k <= this->nrows_; k++ ) {
      entriesPerRow[k] = 0;
    }
    for ( UInt i = 1; i <= storedEntries; i++ ) {
      if ( sparseMat.cidx(i) >= sparseMat.ridx(i) ) {
	entriesPerRow[sparseMat.ridx(i)]++;
      }
    }

    // Now prepare the row pointer array
    rowPtr_[1] = 1;
    for ( Integer k = 1; k <= this->nrows_; k++ ) {
      rowPtr_[k+1] = rowPtr_[k] + entriesPerRow[k];
    }

    // Make a copy of values and column indices
    Integer *pos = auxVec;
    UInt curRow;
    for ( UInt i = 1; i <= (UInt)this->nrows_; i++ ) {
      pos[i] = 0;
    }
    for ( UInt j = 1; j <= storedEntries; j++ ) {
      curRow = sparseMat.ridx(j);
      if ( curRow <= sparseMat.cidx(j) ) {
	data_  [ rowPtr_[curRow] + pos[curRow] ] = sparseMat.val (j);
	colInd_[ rowPtr_[curRow] + pos[curRow] ] = sparseMat.cidx(j);
	pos[ curRow ]++;
      }
    }

    // Destroy auxilliary vector
    DeleteArray( auxVec );

    // Set pattern pool pointer to NULL, since we allocated pattern
    // ourselves
    patternPool_ = NULL;
    patternID_ = NO_PATTERN_ID;
  }


  // **********************
  //   SetSparsityPattern
  // **********************
  template<typename T>
  void SCRS_Matrix<T>::SetSparsityPattern( BaseGraph &graph ) {

    // NOTE:
    //
    // The implementation of this routine expects the graph to be complete,
    // i.e. it stores all neighbourhood relations. It does not expect,
    // however, that the neighbours for a node are ordered in the array
    // obtained with GetGraphRow(). But, it is also expected that the array
    // contains an entry for the node itself.


    Integer rs;
    Integer *index;
    Integer i, k, actPos;


    // Check that no pattern was allocated
    if ( rowPtr_ != NULL || colInd_ != NULL || patternPool_ != NULL ) {
      (*error) << "There seems to already be a sparsity pattern!";
      Error( __FILE__, __LINE__ );
    }

    // Allocate memory for row pointers and initialise first one
    NewArray( rowPtr_, Integer, this->nrows_ + 1 );

    // Allocate memory for column indices
    NewArray( colInd_, Integer, numEntries_ );

    // loop over all rows
    rowPtr_[1] = 1;
    actPos = 1;
    for ( i = 1; i <= this->nrows_; i++ ) {

      // get neighbours (indices and numbers) of i-th node
      index = graph.GetGraphRow(i);
      rs = graph.GetRowSize(i);

      // set column indices of non-zero entries in the right
      // part of this row and count their number
      for ( k = 1; k <= rs; k++ ) {

	// test whether entry lies in the upper triangle
	if ( index[k] >= i ) {

	  // insert column index
	  colInd_[actPos] = index[k];
	  actPos++;

#ifdef DEBUG_SCRS_MATRIX
	  // Report new index
	  (*debug) << " colInd_[" << actPos-1 << "] = " << index[k]
		   << std::endl;
#endif
	}
      }

      // Now we know where the next row starts in data_ and colInd_
      rowPtr_[i+1] = actPos;

    }

    // Check that we are in agreement with the information
    // supplied previously
    if ( numEntries_ != rowPtr_[this->nrows_+1] - 1 ) {
      (*error) << "SCRS_Matrix: Mismatch in number of stored entries!\n"
	       << " numEntries_ = " << numEntries_
	       << "\n rowPtr_[this->nrows_+1] = " << rowPtr_[this->nrows_+1] << '\n';
      Error( __FILE__, __LINE__ );
    }

  }


  // **********************
  //   SetSparsityPattern
  // **********************
  template<typename T>
  void SCRS_Matrix<T>::SetSparsityPattern( PatternPool *patternPool,
                                           PatternIdType patternID ) {


    // Safety checks
    if ( patternPool == NULL ) {
      (*error) << "SCRS_Matrix<T>::SetSparsityPattern: patternPool = NULL! "
               << "This cannot work ;-)";
      Error( __FILE__, __LINE__ );
    }
    if ( rowPtr_ != NULL || colInd_ != NULL || patternPool_ != NULL ) {
      (*error) << "There seems to already be a sparsity pattern!";
      Error( __FILE__, __LINE__ );
    }

    // Set information
    patternPool_ = patternPool;
    patternID_   = patternID;

    // Get hold of pattern
    BaseSparsityPattern *basePattern = NULL;
    basePattern = patternPool_->RegisterUser( patternID );

    SCRS_Pattern *myPattern = NULL;
    myPattern = dynamic_cast<SCRS_Pattern*>( basePattern );
    if ( myPattern == NULL ) {
      Error( WRONG_CAST_MSG, __FILE__, __LINE__);
    }

    // Extract pointers
    colInd_ = myPattern->cidx_;
    rowPtr_ = myPattern->rptr_;
  }


  // *************************
  //   TransferPatternToPool
  // *************************
  template<typename T> PatternIdType
  SCRS_Matrix<T>::TransferPatternToPool( PatternPool *patternPool ) {


    // Safety checks
    if ( patternPool == NULL ) {
      (*error) << "SCRS_Matrix<T>::TransferPatternToPool: "
               << "patternPool = NULL! This cannot work ;-)";
      Error( __FILE__, __LINE__ );
    }
    if ( patternPool_ != NULL ) {
      (*error) << "SCRS_Matrix<T>::TransferPatternToPool: "
               << "Do not own pattern, so I refuse to transfer it!";
      Error( __FILE__, __LINE__ );
    }

    // Generate pattern object for putting into the pool
    SCRS_Pattern *myPattern = New SCRS_Pattern;
    myPattern->cidx_ = colInd_;
    myPattern->rptr_ = rowPtr_;

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
  void SCRS_Matrix<T>::CopySparsityPattern( StdMatrix &mat ) const {

    // Check, if pointer to pattern pool exists
    if( patternPool_ == NULL ) {
      (*error) << "SCRS_Matrix::CopySparsityPattern: only implemented "
               << "for matrices sharing a common pattern!\n";
      Error( __FILE__, __LINE__ );
    }
    
    mat.SetSparsityPattern( patternPool_, patternID_ );

  }
      

  // ***********************************
  //  Matrix-vector multiplication
  //  (specialised for symmetric case)
  // ***********************************
  template<typename T>
  inline void SCRS_Matrix<T>::Mult( const Vector<T> &mvec,
				    Vector<T> &rvec ) const {

    assert(mvec.GetSize() == rvec.GetSize());
    assert(this->ncols_ == (int) mvec.GetSize());
     

    register Integer k, rs;
    register Integer c;
    Integer i, j;

    rvec.Init();

    for ( i = 1; i <= this->nrows_; i++ ) {
      k = rowPtr_[i];
      rs = rowPtr_[i+1] - k;
      rvec[i] += this->data_[k] * mvec[i];
      k++;
      for ( j = 2; j <= rs; j++ ){
        c = colInd_[k];
        rvec[i] += this->data_[k] * mvec[c];
        rvec[c] += opType<T>::multT( this->data_[k], mvec[i] );
        k++;
      }
    }
  }


  // ***********
  //   MultAdd
  // ***********
  template<typename T>
  inline void SCRS_Matrix<T>::MultAdd( const Vector<T> & mvec,
				       Vector<T> & rvec ) const {


    register Integer k,rs;
    register Integer c;
    Integer i, j;

    for ( i = 1; i <= this->nrows_; i++ ) {
      k = rowPtr_[i];
      rs = rowPtr_[i+1] - k;
      rvec[i] += this->data_[k] * mvec[i];
      k++;
      for ( j = 2; j <= rs; j++ ) {
	c = colInd_[k];
	rvec[i] += this->data_[k] * mvec[c];
	rvec[c] += opType<T>::multT( this->data_[k], mvec[i] );
	k++;
      }
    }
  }


  // ***********
  //   MultSub
  // ***********
  template<typename T>
  inline void SCRS_Matrix<T>::MultSub( const Vector<T> &mvec,
				       Vector<T> &rvec ) const {


    register Integer k,rs;
    register Integer c;
    Integer i, j;

    for ( i = 1; i <= this->nrows_; i++ ) {
      k = rowPtr_[i];
      rs = rowPtr_[i+1] - k;
      rvec[i] -= this->data_[k] * mvec[i];
      k++;
      for ( j = 2; j <= rs; j++ ) {
	c = colInd_[k];
	rvec[i] -= this->data_[k] * mvec[c];
	rvec[c] -= opType<T>::multT( this->data_[k], mvec[i] );
	k++;
      }
    }
  }


  // ***********
  //   CompRes
  // ***********
  template<typename T>
  inline void SCRS_Matrix<T>::CompRes( Vector<T> &r, const Vector<T> &x,
				       const Vector<T> &b ) const {


    Integer i;

    // Copy right hand side b to r
    for ( i = 1; i <= this->nrows_; i++ ) {
      r[i] = b[i];
    }

    // Use MultSub (r = r - A*x) to compute the residual
    MultSub( x, r );

  }


  //  ========================================================================

  // In the symmetric case multiplication with the transpose of this matrix
  // is the same as multiplication with the matrix itself, so we can simply
  // fallback to Mult() resp. MultAdd()

  //  ========================================================================


  // **********
  //   MultT
  // **********
  template<typename T>
  inline void SCRS_Matrix<T>::MultT( const Vector<T> & mvec,
				     Vector<T> & rvec ) const {
    Mult( mvec, rvec );
  }


  // ************
  //   MultTAdd
  // ************
  template<typename T>
  inline void SCRS_Matrix<T>::MultTAdd( const Vector<T> & mvec,
					Vector<T> & rvec ) const {
    MultAdd( mvec, rvec );
  }


  // ****************************
  //   Print matrix to a stream
  // ****************************
  template<typename T>
  void SCRS_Matrix<T>::Print(std::ostream& os) const {


    register Integer k,rs;
    register Integer c;
    Integer i, j;

    for ( i = 1; i <= this->nrows_; i++ ) {
      k = rowPtr_[i];
      rs = rowPtr_[i+1] - k;
      os << std::scientific << i << "\t" << i << "\t" << data_[rowPtr_[i]]
	 << '\n';
      k++;
      for ( j = 2; j <= rs; j++ ) {
	c = colInd_[k];
	os << std::scientific << i << "\t" << c << "\t" << data_[k] << '\n';
	os << std::scientific << c << "\t" << i << "\t" << data_[k] << '\n';
	k++;
      }
    }

    os << std::endl;
  }


  // *************************
  //   Export matrix to file
  // *************************
  template<typename T>
  void SCRS_Matrix<T>::Export( const Char *fname, const Char *comment ) const{


    // Obtain blocksize of matrix
    Integer dof = this->GetBlockSize();

    // open output file and check for errors
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
      fprintf( fp, "%%%%MatrixMarket matrix coordinate real symmetric\n" );
    }
    else if ( this->GetEntryType() == COMPLEX ) {
      fprintf( fp, "%%%%MatrixMarket matrix coordinate complex symmetric\n" );
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
	     this->numEntries_ * dof );


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

	// loop over scalar entries in block entry
	for ( ib = 0; ib < dof; ib++ ) {
	  for ( jb = 0; jb < dof; jb++ ) {

	    // store row and column index (note that the MatrixMarket format
            // specifies to store the lower triangle, so we swap the two
            // indices)
	    fprintf( fp, "%6d\t%6d\t",
                     ( colInd_[rowPtr_[i]+j-1] - 1 ) * dof + jb + 1,
                     (i-1) * dof + ib + 1 );

	    // store non-zero entry
	    opType<T>::ExportEntry( this->data_[rowPtr_[i]+j-1], ib, jb,
					  fp );
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



  // *********************************************************************
  //   Copy system matrix to vector containing all upper triangle entries
  // *********************************************************************

  template<typename T>
  void SCRS_Matrix<T>::CopySCRSMatrix2Vec(Complex* &dataVec){
    
    (*error) << "CopySCRSMatrix2Vec is only implemented for complex valued "
             << "scalar matrices!";
    Error( __FILE__, __LINE__ );
  }

  // explicit specialization for complex valued matrices
  template<>
  void SCRS_Matrix<Complex>::CopySCRSMatrix2Vec(Complex* &dataVec){


    UInt k=0;

    dataVec = new Complex[nrows_*nrows_];

    // ------------------------
    //   Write matrix entries
    // ------------------------
    
    Integer i, j, nblocks ;
      
    // loop over all  rows
    for ( i = 1; i <= this->nrows_; i++ ) {
      
      // get number of blocks in i-th row
      nblocks = rowPtr_[i+1] - rowPtr_[i];
      
      //       // loop over all blocks in this row
      for ( j = 1; j <= nblocks; j++ )
        dataVec[(i-1)*nrows_+colInd_[rowPtr_[i]+j-1]-1] =
          Complex(this->data_[rowPtr_[i]+j-1]);
    }    
  }
  

  // *******************************
  //   Add value to a matrix entry
  // *******************************
  template<typename T>
  void SCRS_Matrix<T>::AddToMatrixEntry( Integer i, Integer j, T& v ) {


    // If the entry lies in the lower triangular part, we simply ignore it
    if ( j < i ) {
      return;
    }

    // Otherwise, we add the value to the existing entry
    else {

      bool found = false;
      for ( Integer k = rowPtr_[i]; k < rowPtr_[i+1]; k++ ) {
	if ( colInd_[k] == j ) {
	  found = true;
	  data_[k] += v;
	  break;
	}
      }

      if ( found == false ) {
        std::cerr << " Index pair = (" << i << " , " << j << ")\n";
	Error( "Index pair not found!", __FILE__, __LINE__ );
      }
    }
  }

  // *****************************
  //   Get specific matrix entry
  // *****************************
  template<typename T>
  void SCRS_Matrix<T>::GetMatrixEntry( Integer i, Integer j, 
                                       T &v ) const {
    
    // Determine, which of i and j is row and column index for upper
    // triangular part
    Integer row, col;
    if ( j > i ) {
      row = i;
      col = j;
    } else {
      row = j;
      col = i;    
    }

    bool found = false;
    for ( Integer k = rowPtr_[row]; k < rowPtr_[row+1]; k++ ) {
      if ( colInd_[k] == col ) {
        found = true;
        v = data_[k];
        break;
      }
    }
    if ( found == false ) {
      std::cerr << " Index pair = (" << i << " , " << j << ")\n";
      Error( "Index pair not found!", __FILE__, __LINE__ );
    }
    
  }
  
  // *****************************
  //   Set specific matrix entry
  // *****************************
  template<typename T>
  void SCRS_Matrix<T>::SetMatrixEntry( Integer i, Integer j, T &v, 
                                       bool setCounterPart ) {

    // Check, if a non-diagonal entry should be set without its counterpart
    if ( i!= j && setCounterPart == false ) {
      (*error) << "SCRS_Matrix: Can not set entry at position (" 
               << i <<", " << j << ") without its counterPart!";
      Error( __FILE__, __LINE__ );
    }

    // Determine, which of i and j is row and column index for upper
    // triangular part
    Integer row, col;
    if ( j > i ) {
      row = i;
      col = j;
    } else {
      row = j;
      col = i;    
    }

    bool found = false;
    for ( Integer k = rowPtr_[row]; k < rowPtr_[row+1]; k++ ) {
      if ( colInd_[k] == col ) {
        found = true;
        data_[k] = v;
        break;
      }
    }
    if ( found == false ) {
      std::cerr << " Index pair = (" << i << " , " << j << ")\n";
      Error( "Index pair not found!", __FILE__, __LINE__ );
    }
  }

  // ***********
  //   Insert
  // ***********
  template<typename T>
  void SCRS_Matrix<T>::Insert(Integer row, Integer col,  Integer pos){
    Error( "not implemented", __FILE__, __LINE__ );
  }


  // ***********************
  //  Forced Instantiation
  // ***********************
  template <typename T>
  void SCRS_Matrix<T>::InstantiatePublicMethods() {


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
    SCRS_Matrix<T> dummyMat( 10, 10, 10 );

    dummyMat.Init();

    SetSparsityPattern( dummyGraph );
    SetSize( dummyInt, dummyInt, dummyInt );
    Print( std::cout );
    Export( "" );
    SetRowSize( dummyInt, dummyInt );
    Insert( dummyInt, dummyInt, dummyInt );
    AddToMatrixEntry( dummyInt, dummyInt, dummyEntry );

    Mult    ( dummyVec, dummyVec );
    MultAdd ( dummyVec, dummyVec );
    MultSub ( dummyVec, dummyVec );
    CompRes ( dummyVec, dummyVec, dummyVec );
    MultT   ( dummyVec, dummyVec );
    MultTAdd( dummyVec, dummyVec );
    AccountForPenalty( dummyVec, dummyVec );

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
    dummyPVal = GetDataPointer();

    MatrixEntryType met = this->GetEntryType();
    met = NOENTRYTYPE;

    const Integer* cInt1 = GetRowPointer();
    const Integer* cInt2 = GetColPointer();
    const T* aux = GetDataPointer();

    // use the stuff to satisfy gcc
    std::cout << cInt1 << cInt2;
    dummyBool = dummyBool == true ? false : true;

  }


  // *********
  //   Scale
  // *********
  template<typename T>
  void SCRS_Matrix<T>::Scale( Double factor ) {
    for ( Integer i = 1; i <= numEntries_; i++ ) {
      data_[i] *= factor;
    }
  }


  // **************
  //   GetMaxDiag
  // **************
  template<typename T>
  Double SCRS_Matrix<T>::GetMaxDiag() const {


    double maxDiag = 0;
    double current = 0;
    Integer i;

    for ( i = 1; i <= this->nrows_; i++ ) {

      // use an opType to ensure that tiny matrices
      // are treated correctly
      current = opType<T>::MaxDiag( data_[ rowPtr_[i] ] );
      maxDiag = maxDiag > current ? maxDiag : current;
    }

    return maxDiag;
  }


  // ************************
  //   Add (another matrix)
  // ************************
  template<typename T>
  void SCRS_Matrix<T>::Add( const Double alpha, const StdMatrix& mat ) {
     TRY_CAST {

      // Down-cast input matrix
      ConstRefCast( mat, SCRS_Matrix<T>, scrsMat );

      // Obtain pointer to data array of other matrix
      const T *data = scrsMat.GetDataPointer();

      // We now assume that the matrix have matching
      // dimensions and sparsity patterns
      for ( Integer i = 1; i <= numEntries_; i++ ) {
	data_[i] += alpha * data[i];
      }

    } CATCH_CAST;
  }

  
  // Sepcialisation for complex matrices
  template<>
  void SCRS_Matrix<Complex>::Add( const Double alpha, const StdMatrix& mat ) {

    // Check for entry type of mat
    MatrixEntryType eType = mat.GetEntryType();

    if( eType == DOUBLE ) { 

      TRY_CAST {
        
        // Down-cast input matrix
        ConstRefCast( mat, SCRS_Matrix<Double>, scrsMat );
        
        // Obtain pointer to data array of other matrix
        const Double *data = scrsMat.GetDataPointer();
        
        // We now assume that the matrix have matching
        // dimensions and sparsity patterns
        for ( Integer i = 1; i <= numEntries_; i++ ) {
          data_[i] += alpha * Complex(data[i],0.0);
        }
        
      } CATCH_CAST;
      
    } else {
      
      TRY_CAST {
        
        // Down-cast input matrix
        ConstRefCast( mat, SCRS_Matrix<Complex>, scrsMat );
        
        // Obtain pointer to data array of other matrix
        const Complex *data = scrsMat.GetDataPointer();
        
        // We now assume that the matrix have matching
        // dimensions and sparsity patterns
        for ( Integer i = 1; i <= numEntries_; i++ ) {
          data_[i] += alpha * data[i];
        }
        
      } CATCH_CAST;
      
    }
  }
  


  // *********************
  //   AccountForPenalty
  // *********************
  template<typename T>
  void SCRS_Matrix<T>::AccountForPenalty( Vector<T> &vec,
                                          const Vector<T> &rhs ) const {

      Error( "Method not defined for this kind of template!",
             __FILE__,__LINE__ );
  }


  // *********************
  //   AccountForPenalty
  // *********************
  template<>
  void SCRS_Matrix<Double>::AccountForPenalty( Vector<Double> &vec,
                                               const Vector<Double> &rhs )
    const {


    // Determine value of penalty approach
    Double penalty = GetMaxDiag();

    // We need a safety factor, since at some place, we have some
    // rounding errors in determining GetMaxDiag()
    Double safety = 0.1;

    // Loop over diagonal and look for Dirichlet cases
    Double diag, entry;
    for ( Integer k = 1; k <= this->nrows_; k++ ) {
      GetDiagEntry( k, diag );
      if ( diag > penalty * safety ) {
	rhs.GetVectorEntry( k, entry );
	vec.SetVectorEntry( k, entry / penalty );
      }
    }
  }

}
