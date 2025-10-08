#include <cassert>
#include <sstream>

#include "CoordFormat.hh"
#include "SCRS_Matrix.hh"
#include "opdefs.hh"
#include "Matrix.hh"
#include "Utils/SyncAccess.hh"

// Implementation of methods for the symmetric compressed row storage SCRS
// matrix class

namespace CoupledField {

//#define DEBUG_SCRS_MATRIX

  // ********************
  //   Copy Constructor
  // ********************

  template<typename T>
  SCRS_Matrix<T>::SCRS_Matrix( const SCRS_Matrix<T> &origMat ) {
    
    colInd_      = NULL;
    rowPtr_      = NULL;
    diagPtr_     = NULL;
    data_        = NULL;
    patternPool_ = NULL;
    patternID_   = NO_PATTERN_ID;

    // Set basic size informations
    numEntries_ = origMat.numEntries_;
    this->nnz_        = origMat.nnz_;
    this->nrows_      = origMat.nrows_;
    this->ncols_      = origMat.ncols_;

    // ---------------------
    // Pattern is not shared
    // ---------------------

    if ( origMat.patternPool_ == NULL ) {

      patternPool_ = NULL;
      patternID_   = NO_PATTERN_ID;
      
      // Allocate memory for internal arrays
      NEWARRAY( colInd_, UInt, numEntries_ );
      NEWARRAY( rowPtr_, UInt, this->nrows_ + 1 );
      NEWARRAY( diagPtr_, UInt, sizeof(origMat.diagPtr_) );
      NEWARRAY( data_, T, numEntries_ );

      // Copy information
      for ( UInt i = 0; i < (UInt)numEntries_; i++ ) {
        data_[i]   = origMat.data_[i];
        colInd_[i] = origMat.colInd_[i];
      }

      for ( UInt i = 0; i < (UInt)this->nrows_ + 1; i++ ) {
        rowPtr_[i] = origMat.rowPtr_[i];
      }
      // only try to allocate if there is diagPtr_ in original matrix
      if (origMat.diagPtr_){
        for ( UInt i = 0; i < sizeof(origMat.diagPtr_); i++ ) {
          diagPtr_[i] = origMat.diagPtr_[i];
        }
      }
    }

    else {
      // -----------------
      // Pattern is shared
      // -----------------
      
      // Copy sparsity pattern from poool
      this->SetSparsityPattern( origMat.patternPool_, origMat.patternID_ );

      // Generate copy of data array
      NEWARRAY( data_, T, numEntries_ );
      for ( UInt i = 0; i < (UInt)numEntries_; i++ ) {
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
    EXCEPTION( "Alternate Constructor for SCRS_Matrix only works for RMatrix1" );
  }


  // Specialisation for the RMatrix1 case
  template<>
  SCRS_Matrix<double>::SCRS_Matrix( CoordFormat<double> &sparseMat ) : diagPtr_( NULL ){


    // Set general matrix pattern and dimension info
    this->nrows_ = sparseMat.GetNumRows();
    this->ncols_ = sparseMat.GetNumCols();
    this->nnz_   = sparseMat.GetNNZ();
    numEntries_ = (this->nnz_ + this->nrows_) / 2;
    UInt storedEntries = sparseMat.GetNumEntries();

    // Test that matrix is square
    if ( this->nrows_ != this->ncols_ ) {
      EXCEPTION( "The input matrix is of dimension " << this->nrows_ << " x "
	       << this->ncols_ << ", but should be square!" );
    }

    // Make a simply consistency check
    if ( sparseMat.GetSymmStorageFlag() == true ) {

      // In this case the number of entries stored in
      // sparseMat should be identical to the number
      // of entries we are going to store
      if ( storedEntries != (UInt)numEntries_ ) {
        EXCEPTION( "Inconsistency detected in input matrix in CoordFormat! "
            << "Matrix claims to be in symmetric storage format, but "
            << "fails test on number of entries! "
            << "storedEntries = " << storedEntries << " , "
            << "numEntries_ = " << numEntries_ );
      }
    }

    else {

      // In this case the number of stored entries should be identical
      // to the number of non-zero matrix entries
      if ( storedEntries != (UInt)this->nnz_ ) {
        EXCEPTION( "Inconsistency detected in input matrix in CoordFormat!\n"
            << "Matrix claims to be in general storage format, but "
            << "fails test on number of entries!" );
      }
    }

    // Allocate memory
    NEWARRAY( rowPtr_, UInt, this->nrows_ + 1  );
    NEWARRAY( colInd_, UInt, numEntries_ );
    NEWARRAY( data_  , Double , numEntries_ );

    // Sort input matrix, such that it conforms to our layout requirements
    sparseMat.Sort( CoordFormat<double>::PURELEX );

    // Generate auxilliary vector
    UInt *auxVec = NULL;
    NEWARRAY( auxVec, UInt, this->nrows_ );

    // Determine the number of non-zero entries in the different rows.
    // In doing this we ignore entries in the lower triangle.
    UInt *entriesPerRow = auxVec;
    for ( UInt k = 0; k < this->nrows_; k++ ) {
      entriesPerRow[k] = 0;
    }
    for ( UInt i = 0; i < storedEntries; i++ ) {
      if ( sparseMat.cidx(i) >= sparseMat.ridx(i) ) {
        entriesPerRow[sparseMat.ridx(i)]++;
      }
    }

    // Now prepare the row pointer array
    rowPtr_[0] = 0;
    for ( UInt k = 0; k < this->nrows_; k++ ) {
      rowPtr_[k+1] = rowPtr_[k] + entriesPerRow[k];
    }

    // Make a copy of values and column indices
    UInt *pos = auxVec;
    UInt curRow;
    for ( UInt i = 0; i < (UInt)this->nrows_; i++ ) {
      pos[i] = 0;
    }
    for ( UInt j = 0; j < storedEntries; j++ ) {
      curRow = sparseMat.ridx(j);
      if ( curRow <= sparseMat.cidx(j) ) {
        data_  [ rowPtr_[curRow] + pos[curRow] ] = sparseMat.val (j);
        colInd_[ rowPtr_[curRow] + pos[curRow] ] = sparseMat.cidx(j);
        pos[ curRow ]++;
      }
    }

    // Destroy auxilliary vector
    delete[] (auxVec);
    NEWARRAY( diagPtr_, UInt, this->nrows_ );
    // Set pattern pool pointer to NULL, since we allocated pattern
    // ourselves
    patternPool_ = NULL;
    patternID_ = NO_PATTERN_ID;
  }


  // Specialisation for the CMatrix1 case
  template<>
  SCRS_Matrix<Complex>::SCRS_Matrix( CoordFormat<Complex> &sparseMat ) : diagPtr_( NULL ) {


    // Set general matrix pattern and dimension info
    this->nrows_ = sparseMat.GetNumRows();
    this->ncols_ = sparseMat.GetNumCols();
    this->nnz_   = sparseMat.GetNNZ();
    numEntries_ = (this->nnz_ + this->nrows_) / 2;
    UInt storedEntries = sparseMat.GetNumEntries();

    // Test that matrix is square
    if ( this->nrows_ != this->ncols_ ) {
      EXCEPTION( "The input matrix is of dimension " << this->nrows_ << " x "
                 << this->ncols_ << ", but should be square!" );
    }

    // Make a simply consistency check
    if ( sparseMat.GetSymmStorageFlag() == true ) {

      // In this case the number of entries stored in
      // sparseMat should be identical to the number
      // of entries we are going to store
      if ( storedEntries != (UInt)numEntries_ ) {
        EXCEPTION( "Inconsistency detected in input matrix in CoordFormat!\n"
            << "Matrix claims to be in symmetric storage format, but "
            << "fails test on number of entries!" );
      }
    }

    else {

      // In this case the number of stored entries should be identical
      // to the number of non-zero matrix entries
      if ( storedEntries != (UInt)this->nnz_ ) {
        EXCEPTION( "Inconsistency detected in input matrix in CoordFormat!\n"
            << "Matrix claims to be in general storage format, but "
            << "fails test on number of entries!" );
      }
    }

    // Allocate memory
    NEWARRAY( rowPtr_, UInt, this->nrows_ + 1  );
    NEWARRAY( colInd_, UInt, numEntries_ );
    NEWARRAY( data_  , Complex, numEntries_ );

    // Sort input matrix, such that it conforms to our layout requirements
    sparseMat.Sort( CoordFormat<Complex>::PURELEX );

    // Generate auxilliary vector
    Integer *auxVec = NULL;
    NEWARRAY( auxVec, Integer, this->nrows_ );

    // Determine the number of non-zero entries in the different rows.
    // In doing this we ignore entries in the lower triangle.
    Integer *entriesPerRow = auxVec;
    for ( UInt k = 0; k < this->nrows_; k++ ) {
      entriesPerRow[k] = 0;
    }
    for ( UInt i = 0; i < storedEntries; i++ ) {
      if ( sparseMat.cidx(i) >= sparseMat.ridx(i) ) {
        entriesPerRow[sparseMat.ridx(i)]++;
      }
    }

    // Now prepare the row pointer array
    rowPtr_[0] = 0;
    for ( UInt k = 0; k < this->nrows_; k++ ) {
      rowPtr_[k+1] = rowPtr_[k] + entriesPerRow[k];
    }

    // Make a copy of values and column indices
    Integer *pos = auxVec;
    UInt curRow;
    for ( UInt i = 0; i < (UInt)this->nrows_; i++ ) {
      pos[i] = 0;
    }
    for ( UInt j = 0; j < storedEntries; j++ ) {
      curRow = sparseMat.ridx(j);
      if ( curRow <= sparseMat.cidx(j) ) {
        data_  [ rowPtr_[curRow] + pos[curRow] ] = sparseMat.val (j);
        colInd_[ rowPtr_[curRow] + pos[curRow] ] = sparseMat.cidx(j);
        pos[ curRow ]++;
      }
    }
    NEWARRAY( diagPtr_, UInt, this->nrows_ );
    // Destroy auxilliary vector
    delete[] (auxVec);

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

    // Obtain block definition from graph object
    this->rowBlocks_ = graph.GetRowBlockDefinition();
    this->colBlocks_= graph.GetColBlockDefinition();

    UInt rs;
    UInt *index;
    UInt i, k, actPos;


    // Check that no pattern was allocated
    if ( rowPtr_ != NULL || colInd_ != NULL || patternPool_ != NULL ) {
      EXCEPTION( "There seems to already be a sparsity pattern!" );
    }

    // Allocate memory for row pointers and initialise first one
    NEWARRAY( rowPtr_, UInt, this->nrows_ + 1 );

    // Allocate memory for column indices
    NEWARRAY( colInd_, UInt, numEntries_ );

    // loop over all rows
    rowPtr_[0] = 0;
    actPos = 0;
    for ( i = 0; i < this->nrows_; i++ ) {

      // get neighbours (indices and numbers) of i-th node
      index = graph.GetGraphRow(i);
      rs = graph.GetRowSize(i);

      // set column indices of non-zero entries in the right
      // part of this row and count their number
      for ( k = 0; k < rs; k++ ) {

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
    if ( numEntries_ != rowPtr_[this->nrows_] ) {
      EXCEPTION( "SCRS_Matrix: Mismatch in number of stored entries!\n"
          << " numEntries_ = " << numEntries_
          << "\n rowPtr_[this->nrows_+1] = " << rowPtr_[this->nrows_] << '\n' );
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
      EXCEPTION( "SCRS_Matrix<T>::SetSparsityPattern: patternPool = NULL! "
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

    SCRS_Pattern *myPattern = NULL;
    myPattern = dynamic_cast<SCRS_Pattern*>( basePattern );
    if ( myPattern == NULL ) {
      EXCEPTION( WRONG_CAST_MSG );
    }

    // Extract pointers
    colInd_ = myPattern->cidx_;
    rowPtr_ = myPattern->rptr_;
  }

  template<typename T>
  void SCRS_Matrix<T>::SetSize( UInt nrows, UInt ncols, UInt nnz ) {
    UInt newNumEntries = (nnz + nrows) / 2;

    if(patternPool_ && ((this->nrows_ != nrows) || (this->numEntries_ != newNumEntries))){ // If we don't own the pattern, deregister and clear pattern
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
      delete[] rowPtr_;
      rowPtr_ = nullptr;
    }
    this->nnz_ = nnz;

    //this is correct iff there are no 0s on the diagonal!
    if ( numEntries_!= newNumEntries) {
      delete[] data_;
      data_ = nullptr;
      delete[] colInd_;
      colInd_ = nullptr;
    }
    numEntries_ = newNumEntries;
    
    if(!rowPtr_){
      rowPtr_  = new UInt[this->nrows_ + 1];
      rowPtr_[0]  = 0;
    }
    if(!diagPtr_){
      diagPtr_ = new UInt[this->nrows_];
      diagPtr_[0] = 0;
    }
    if(!colInd_){
      colInd_ = new UInt[this->numEntries_];
    }
    if(!data_){
      data_ = new T[this->numEntries_];
    }
  }

  // **********************
  //   SetSparsityPatternData
  // **********************
  template<typename T>
  void SCRS_Matrix<T>::SetSparsityPatternData( const StdVector<UInt>& rowP,
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
  SCRS_Matrix<T>::TransferPatternToPool( PatternPool *patternPool ) {


    // Safety checks
    if ( patternPool == NULL ) {
      EXCEPTION( "SCRS_Matrix<T>::TransferPatternToPool: "
               << "patternPool = NULL! This cannot work ;-)" );
    }
    if ( patternPool_ != NULL ) {
      EXCEPTION( "SCRS_Matrix<T>::TransferPatternToPool: "
               << "Do not own pattern, so I refuse to transfer it!" );
    }

    // Generate pattern object for putting into the pool
    SCRS_Pattern *myPattern = new SCRS_Pattern;
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
  void SCRS_Matrix<T>::SetSparsityPattern( const StdMatrix &srcMat ) {

    if ( patternPool_ != NULL ) {
      EXCEPTION("Setting of sparsity pattern not allowed if a pattern pool is present.");
    }
    
    if( srcMat.GetEntryType() == BaseMatrix::DOUBLE ) 
    {
      const SCRS_Matrix<Double>& mat = dynamic_cast< const SCRS_Matrix<Double>& >(srcMat);
      const UInt* srcColInd = mat.GetColPointer();
      const UInt* srcRowPtr = mat.GetRowPointer();

      delete[] colInd_;
      delete[] rowPtr_;
      delete[] data_;

      NEWARRAY( colInd_, UInt, mat.GetNumEntries() );
      NEWARRAY( rowPtr_, UInt, mat.GetNumRows() + 1 );
      NEWARRAY( data_, T, mat.GetNumEntries() );

      // Copy information
      for ( UInt i = 0; i < (UInt)numEntries_; i++ ) {
        colInd_[i] = srcColInd[i];
      }
      
      for ( UInt i = 0; i < (UInt)numEntries_; i++ ) {
        data_[i] = T(0.0);
      }

      for ( UInt i = 0; i < (UInt)this->nrows_ + 1; i++ ) {
        rowPtr_[i] = srcRowPtr[i];
      }
    } else 
    {
    	if(srcMat.GetEntryType() == BaseMatrix::COMPLEX)
    	{
    	      const SCRS_Matrix<Complex>& mat = dynamic_cast< const SCRS_Matrix<Complex>& >(srcMat);
    	      const UInt* srcColInd = mat.GetColPointer();
    	      const UInt* srcRowPtr = mat.GetRowPointer();

    	      delete[] colInd_;
    	      delete[] rowPtr_;
    	      delete[] data_;

    	      NEWARRAY( colInd_, UInt, mat.GetNumEntries() );
    	      NEWARRAY( rowPtr_, UInt, mat.GetNumRows() + 1 );
    	      NEWARRAY( data_, T, mat.GetNumEntries() );

    	      // Copy information
    	      for ( UInt i = 0; i < (UInt)numEntries_; i++ ) {
    	        colInd_[i] = srcColInd[i];
    	      }

    	      for ( UInt i = 0; i < (UInt)numEntries_; i++ ) {
    	        data_[i] = T(0.0);
    	      }

    	      for ( UInt i = 0; i < (UInt)this->nrows_ + 1; i++ ) {
    	        rowPtr_[i] = srcRowPtr[i];
    	      }

    	}	else
    		{
    			EXCEPTION("SCRS_Matrix<T>::SetSparsityPattern not yet implemented for complex matrices.");
    		}
    	}
    }
    
    



  // ***********************************
  //  Matrix-vector multiplication
  //  (specialised for symmetric case)
  // ***********************************
  template<typename T>
  inline void SCRS_Matrix<T>::Mult( const Vector<T> &mvec,
                                    Vector<T> &rvec ) const {

    assert(mvec.GetSize() == rvec.GetSize());
    assert(this->ncols_ == mvec.GetSize());

    UInt k, rs;
    UInt c;
    UInt i, j;

    rvec.Init();

    for ( i = 0; i < this->nrows_; i++ ) {
      k = rowPtr_[i];
      rs = rowPtr_[i+1] - k;
      rvec[i] += this->data_[k] * mvec[i];
      k++;
      for ( j = 1; j < rs; j++ ){
        c = colInd_[k];
        rvec[i] += this->data_[k] * mvec[c];
        rvec[c] += OpType<T>::multT( this->data_[k], mvec[i] );
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

    UInt k,rs;
    UInt c;
    UInt i, j;

    for ( i = 0; i < this->nrows_; i++ ) {
      k = rowPtr_[i];
      rs = rowPtr_[i+1] - k;
      rvec[i] += this->data_[k] * mvec[i];
      k++;
      for ( j = 1; j < rs; j++ ) {
        c = colInd_[k];
        rvec[i] += this->data_[k] * mvec[c];
        rvec[c] += OpType<T>::multT( this->data_[k], mvec[i] );
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


    UInt k,rs;
    UInt c;
    UInt i, j;

    for ( i = 0; i < this->nrows_; i++ ) {
      k = rowPtr_[i];
      rs = rowPtr_[i+1] - k;
      rvec[i] -= this->data_[k] * mvec[i];
      k++;
      for ( j = 1; j < rs; j++ ) {
        c = colInd_[k];
        rvec[i] -= this->data_[k] * mvec[c];
        rvec[c] -= OpType<T>::multT( this->data_[k], mvec[i] );
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


    UInt i;

    // Copy right hand side b to r
    for ( i = 0; i < this->nrows_; i++ ) {
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
  
  // ************
  //   MultTSub
  // ************
  template<typename T>
  inline void SCRS_Matrix<T>::MultTSub( const Vector<T> & mvec,
                                        Vector<T> & rvec ) const {
    MultSub( mvec, rvec );
  }


  // ****************************
  //   Print matrix to a stream
  // ****************************
  template<typename T>
  std::string SCRS_Matrix<T>::ToString( char colSeparator,
                                        char rowSeparator ) const {


    std::stringstream os;
    UInt k,rs;
    UInt c;
    UInt i, j;

    for ( i = 0; i < this->nrows_; i++ ) {
      k = rowPtr_[i];
      rs = rowPtr_[i+1] - k;
      os << std::scientific << i << colSeparator 
                            << i << colSeparator 
                            << data_[rowPtr_[i]] << rowSeparator;
      k++;
      for ( j = 1; j < rs; j++ ) {
        c = colInd_[k];
        os << std::scientific 
           << i << colSeparator 
           << c << colSeparator 
           << data_[k] << rowSeparator;
        os << std::scientific 
           << c << colSeparator 
           << i << colSeparator 
           << data_[k] << rowSeparator;
        k++;
      }
    }

    os << std::endl;
    return os.str();
  }


  template<typename T>
  std::string SCRS_Matrix<T>::Dump() const
  {
    // this is ugly copy&paste from CRS_Matrix :(
    std::stringstream ss;
    // don't use ToString() from this class but the glocal ToString() from tools.hh
    ss << " row=" << ::ToString<unsigned int>(rowPtr_, this->nrows_+1) << std::endl;
    ss << " col=" << ::ToString<unsigned int>(colInd_, this->nnz_)  << std::endl;
    ss << " val=" << ::ToString<T>(data_, this->nnz_);
    return ss.str();
  }


  // *************************
  //   Export matrix to file
  // *************************
  template<typename T>
  void SCRS_Matrix<T>::ExportMatrixMarket(const std::string& fname, const std::string& comment ) const{


    // open output file and check for errors
    FILE *fp = fopen( fname.c_str(), "w" );
    if ( fp == NULL ) {
      EXCEPTION( "Cannot open file " << fname << " for writing!" );
    }

    // ---------------------
    //   Write file header
    // ---------------------

    // Matrix Market Format Specification
    if ( this->GetEntryType() == BaseMatrix::DOUBLE ) {
      fprintf( fp, "%%%%MatrixMarket matrix coordinate real symmetric\n" );
    }
    else if ( this->GetEntryType() == BaseMatrix::COMPLEX ) {
      fprintf( fp, "%%%%MatrixMarket matrix coordinate complex symmetric\n" );
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
    fprintf( fp, "%d\t%d\t%d\n", this->nrows_, this->ncols_,
             this->numEntries_ );


    // ------------------------
    //   Write matrix entries
    // ------------------------

    UInt i, j, nblocks; // TODO: Check if this is still needed
    // Integer i, j, ib, jb, nblocks; // TODO: Check if this is still needed

    // loop over all block rows
    for ( i = 0; i < this->nrows_; i++ ) {

      // get number of blocks in i-th row
      nblocks = rowPtr_[i+1] - rowPtr_[i];

      // loop over all blocks in this row
      for ( j = 0; j < nblocks; j++ ) {

        // store row and column index (note that the MatrixMarket format
        // specifies to store the lower triangle, so we swap the two
        // indices)
        fprintf( fp, "%6d\t%6d\t",
                 ( colInd_[rowPtr_[i]+j]  )  + 1,
                 i  + 1 );

        // store non-zero entry
        OpType<T>::ExportEntry( this->data_[rowPtr_[i]+j], 0, 0,
                                fp );
        fprintf( fp, "\n" );
      }
    }

    // close output file
    if ( fclose( fp ) == EOF ) {
      WARN( "Could not close file " << fname << " after writing!" );
    }
  }



  // *********************************************************************
  //   Copy system matrix to vector containing all upper triangle entries
  // *********************************************************************

  template<typename T>
  void SCRS_Matrix<T>::CopySCRSMatrix2Vec(Complex* &dataVec){

    EXCEPTION( "CopySCRSMatrix2Vec is only implemented for complex valued "
             << "scalar matrices!" );
  }

  // explicit specialization for complex valued matrices
  template<>
  void SCRS_Matrix<Complex>::CopySCRSMatrix2Vec(Complex* &dataVec){

    // UInt k=0;  // TODO: Check if this is still needed

    dataVec = new Complex[nrows_*nrows_];

    // ------------------------
    //   Write matrix entries
    // ------------------------

    UInt i, j, nblocks ;

    // loop over all  rows
    for ( i = 0; i < this->nrows_; i++ ) {

      // get number of blocks in i-th row
      nblocks = rowPtr_[i+1] - rowPtr_[i];

      //       // loop over all blocks in this row
      for ( j = 0; j < nblocks; j++ )
        dataVec[(i-1)*nrows_+colInd_[rowPtr_[i]+j-1]-1] =
          Complex(this->data_[rowPtr_[i]+j-1]);
    }
  }


  // *******************************
  //   Add value to a matrix entry
  // *******************************

  template<typename T>
inline void SCRS_Matrix<T>::AddToMatrixEntry( UInt i, UInt j, const T& v )
  {
    // If the entry lies in the lower triangular part, we simply ignore it
    if ( j < i ) {
      return;
    }

    // Otherwise, we add the value to the existing entry
    else {

      bool found = false;
      // logarithmic search (this has complexity O(log(n)), n=u-l)
      UInt l = rowPtr_[i];
      UInt u = rowPtr_[i+1]-1;
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

      if ( found == false ) {
        std::cerr << " Index pair = (" << i << " , " << j << ")\n";
        EXCEPTION( "Index pair not found!" );
      }
    }
  }


  // *****************************
  //   Get specific matrix entry
  // *****************************
  template<typename T>
  void SCRS_Matrix<T>::GetMatrixEntry( UInt i, UInt j,
                                       T &v ) const {

    // Determine, which of i and j is row and column index for upper
    // triangular part
    UInt row, col;
    if ( j > i ) {
      row = i;
      col = j;
    } else {
      row = j;
      col = i;
    }

    bool found = false;
    // logarithmic search (this has complexity O(log(n)), n=u-l)
    UInt l = rowPtr_[row];
    UInt u = rowPtr_[row+1]-1;
    while(l <= u){
      UInt k = (l+u) >> 1;
      if(colInd_[k] > col){
        u = k-1;
      }else if(colInd_[k] < col){
        l = k+1;
      }else{
        v = data_[k];
        found = true;
        break;
      }
    }
    if ( found == false ) {
      std::cerr << " Index pair = (" << i << " , " << j << ")\n";
      EXCEPTION( "Index pair not found!" );
    }

  }

  // ******************
   //   GetMemoryUsage
   // ******************
   template <typename T>
   Double SCRS_Matrix<T>::GetMemoryUsage() const {
     Double sum = 0.0;
     
     // check if pattern is shared
     if( !patternPool_ ) {
       sum += ( (this->nrows_ + 1)   // rowPtr
           + this->numEntries_ )     // colInd
           * sizeof(UInt);
     }
     sum += this->nnz_ * sizeof(T);  // data
     
     return sum;
   }
  
  // *****************************
  //   Set specific matrix entry
  // *****************************
  template<typename T>
  void SCRS_Matrix<T>::SetMatrixEntry( UInt i, UInt j, const T &v,
                                       bool setCounterPart ) {

    // Check, if a non-diagonal entry should be set without its counterpart
    if ( i!= j && setCounterPart == false ) {
      EXCEPTION( "SCRS_Matrix: Can not set entry at position ("
               << i <<", " << j << ") without its counterPart!" );
    }

    // Determine, which of i and j is row and column index for upper
    // triangular part
    UInt row, col;
    if ( j > i ) {
      row = i;
      col = j;
    } else {
      row = j;
      col = i;
    }

    bool found = false;
    // logarithmic search (this has complexity O(log(n)), n=u-l)
    UInt l = rowPtr_[row];
    UInt u = rowPtr_[row+1]-1;
    while(l <= u){
      UInt k = (l+u) >> 1;
      if(colInd_[k] > col){
        u = k-1;
      }else if(colInd_[k] < col){
        l = k+1;
      }else{
        SyncAccess<SYNC_DATA>::Set(data_[k],v);
        found = true;
        break;
      }
    }
    if ( found == false ) {
      std::cerr << " Index pair = (" << i << " , " << j << ")\n";
      EXCEPTION( "Index pair not found!" );
    }
  }

  // ***********
  //   Insert
  // ***********
  template<typename T>
  void SCRS_Matrix<T>::Insert(UInt row, UInt col,  UInt pos){
    EXCEPTION( "not implemented" );
  }

  // *********
  //   Scale
  // *********
  template<>
  void SCRS_Matrix<Double>::Scale( Double factor ) {
    for ( UInt i = 0; i < numEntries_; i++ ) {
      data_[i] *= factor;
    }
  }
  
  template<>
  void SCRS_Matrix<Complex>::Scale( Double factor ) {
    for ( UInt i = 0; i < numEntries_; i++ ) {
      data_[i] *= factor;
    }
  }

  template<>
  void SCRS_Matrix<Double>::Scale( Complex factor ) {
	  EXCEPTION("CRS: Matrix is Double; you can't multiply be a complex value");
  }


  template<>
  void SCRS_Matrix<Complex>::Scale( Complex factor ) {
    for ( UInt i = 0; i < numEntries_; i++ ) {
      data_[i] *= factor;
    }
  }


  // ************************
  //   Scale on index subset
  // ************************
  template<typename T>
  void SCRS_Matrix<T>::Scale( Double factor,
                              const std::set<UInt>& rowIndices,
                              const std::set<UInt>& colIndices ) {
    EXCEPTION("Implement me");
  }

  // ************************
  //   Add (another matrix)
  // ************************
  template<typename T>
  void SCRS_Matrix<T>::Add( const Double alpha, const StdMatrix& mat )
  {
    // Down-cast input matrix
    const SCRS_Matrix<T>& scrsMat = dynamic_cast<const SCRS_Matrix<T>&>(mat);

    // Obtain pointer to data array of other matrix
    const T *data = scrsMat.GetDataPointer();

    // We now assume that the matrix have matching
    // dimensions and sparsity patterns
    for ( UInt i = 0; i < numEntries_; i++ ) {
      data_[i] += alpha * data[i];
    }
  }

  template<typename T>
  void SCRS_Matrix<T>::Add( const Complex alpha, const StdMatrix& mat )
  {
    EXCEPTION("Complex Add is not implemented for SCRS_Matrix.");
  }



  // Specialization for complex matrices
  template<>
  void SCRS_Matrix<Complex>::Add( const Double alpha, const StdMatrix& mat ) {

    // Check for entry type of mat
    BaseMatrix::EntryType eType = mat.GetEntryType();

    if( eType == BaseMatrix::DOUBLE ) {
      // Down-cast input matrix
      const SCRS_Matrix<Double>& scrsMat = dynamic_cast<const SCRS_Matrix<Double>&>(mat);

      // Obtain pointer to data array of other matrix
      const Double *data = scrsMat.GetDataPointer();

      // We now assume that the matrix have matching
      // dimensions and sparsity patterns
      for ( UInt i = 0; i < numEntries_; i++ ) {
        data_[i] += alpha * Complex(data[i],0.0);
      }

    }
    else
    {
      // Down-cast input matrix
      const SCRS_Matrix<Complex>& scrsMat = dynamic_cast<const SCRS_Matrix<Complex>&>(mat);

      // Obtain pointer to data array of other matrix
      const Complex *data = scrsMat.GetDataPointer();

      // We now assume that the matrix have matching
      // dimensions and sparsity patterns
      for ( UInt i = 0; i < numEntries_; i++ ) {
        data_[i] += alpha * data[i];
      }
    }
  }

  template<>
  void SCRS_Matrix<Complex>::Add( const Complex alpha, const StdMatrix& mat ) {
    EXCEPTION("Complex Add is not implemented for SCRS_Matrix.");
  }
  
  // ******************************************
  //   Add (another matrix, only index subset)
  // ******************************************
  template<typename T>
  void SCRS_Matrix<T>::Add( const Complex alpha, const StdMatrix& mat,
                            const std::set<UInt>& rowIndices,
                            const std::set<UInt>& colIndices) {
    EXCEPTION("Complex Add is not implemented for SCRS_Matrix.");
  }

  template<typename T>
  void SCRS_Matrix<T>::Add( const Double alpha, const StdMatrix& mat,
                            const std::set<UInt>& rowIndices,
                            const std::set<UInt>& colIndices) {
    // Down-cast input matrix
    const SCRS_Matrix<T>& scrsMat = dynamic_cast<const SCRS_Matrix<T>&>(mat);

    // Obtain pointer to data array of other matrix
    const T *data = scrsMat.GetDataPointer();
    
    // Distinguish 4 cases:
    // 1) Neither row- nor col-indices are set (i.e. take all indices)
    //    -> use standard Add methods
    // 2) Row and col-indices are set and contain all rows / cols
    //    -> use standard Add methods
    if( (rowIndices.size() == 0 && colIndices.size() == 0) ||
        (rowIndices.size() == this->nrows_ &&
         colIndices.size() == this->ncols_ ) ) {
      // use standard method
      this->Add(alpha, mat);
      return;
    }
    
    // 3) Both sets are non-empty sub-sets of all indices
    //    --> loop over row / column indices to be set
    std::set<UInt>::const_iterator rowIt, colIt;
    if( rowIndices.size() > 0 && colIndices.size() > 0 ) {
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
      //    -> either loop over all rows and take into account only
      //       selected column
      //    -> or loop over selected rows and take into account

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
        const T *data = dynamic_cast<const SCRS_Matrix<T>&>(mat).GetDataPointer();
        //another optimization check if data in set is continuous
        std::set<UInt>::iterator starting  = colIndices.begin();
        std::set<UInt>::iterator ending = colIndices.end();
        UInt min = *starting;
        UInt max = *colIndices.rbegin();
        UInt size = colIndices.size();
        if( (max - min) == (size-1)){
#pragma omp parallel for
          for ( Integer i = 0; i < (Integer) numEntries_; i++ ) {
            if (colInd_[i] >= min && colInd_[i] <= max)
              data_[i] += alpha * data[i];
          }
        }else{
          // we have to searach the set for every entry (nnz_ * O(log sizeof(colIndices)))
#pragma omp parallel for
          for ( Integer i = 0; i < (Integer) numEntries_; i++ ) {
            if (colIndices.find(colInd_[i]) != ending)
              data_[i] += alpha * data[i];
          }
        }
      } // if column set is nonzero
      //EXCEPTION( "This method is only implemented for the case of "
      //    << "non-empty row and column index sets!" );
    }
  }
  template<typename T>
  void SCRS_Matrix<T>::GetNumBlocks( UInt& nRowBlocks, UInt& nColBlocks, 
                                     UInt& numBlocks ) const {
    nRowBlocks = rowBlocks_.GetSize();
    nColBlocks = colBlocks_.GetSize();
    numBlocks = nRowBlocks * nColBlocks;
  }

  template<typename T>
  void SCRS_Matrix<T>::GetDiagBlock( UInt blockRow, DenseMatrix& diagBlock ) const {
    // obtain row/col indices of given block
    const std::pair<UInt,UInt> & rowSet = rowBlocks_[blockRow];
    const std::pair<UInt,UInt> & colSet = colBlocks_[blockRow];
    
    // copy entries to diagBlock
    UInt sizeRow = rowSet.second - rowSet.first + 1;
    UInt sizeCol = colSet.second - colSet.first + 1;
    
    Matrix<T> & mat = static_cast<Matrix<T>& >(diagBlock);
    mat.Resize( sizeRow, sizeCol );
    mat.Init();
    UInt i, j, k, rs, c, rr, rc;
    for ( i = rowSet.first; i < rowSet.second+1; i++ ) {
      k = rowPtr_[i];
      rs = rowPtr_[i+1] - k;
      rr = i - rowSet.first;
      mat[rr][rr] = this->data_[k]; // diagonal entry
      k++;
      
      for( j = 1; j < rs; ++j ) {
        c = colInd_[k];
        rc = c-colSet.first;
        if( c < colSet.second + 1 && c >= colSet.first ) {
          mat[rr][rc] = this->data_[k]; 
          mat[rc][rr] = this->data_[k]; // transposed part
        }
        k++;
      }
    }
    
  }
// Explicit template instantiation
  template class SCRS_Matrix<Double>;
  template class SCRS_Matrix<Complex>;
#ifdef MSVC
  template void SCRS_Matrix<Double>::Mult(const Vector<Double>&, Vector<Double>&) const;
  template void SCRS_Matrix<Complex>::Mult(const Vector<Complex>&, Vector<Complex>&) const;
  template void SCRS_Matrix<Double>::MultT(const Vector<Double>&, Vector<Double>&) const;
  template void SCRS_Matrix<Complex>::MultT(const Vector<Complex>&, Vector<Complex>&) const;
  template void SCRS_Matrix<Double>::MultAdd(const Vector<Double>&, Vector<Double>&) const;
  template void SCRS_Matrix<Complex>::MultAdd(const Vector<Complex>&, Vector<Complex>&) const;
  template void SCRS_Matrix<Double>::MultTAdd(const Vector<Double>&, Vector<Double>&) const;
  template void SCRS_Matrix<Complex>::MultTAdd(const Vector<Complex>&, Vector<Complex>&) const;
  template void SCRS_Matrix<Double>::MultSub(const Vector<Double>&, Vector<Double>&) const;
  template void SCRS_Matrix<Complex>::MultSub(const Vector<Complex>&, Vector<Complex>&) const;
  template void SCRS_Matrix<Double>::MultTSub(const Vector<Double>&, Vector<Double>&) const;
  template void SCRS_Matrix<Complex>::MultTSub(const Vector<Complex>&, Vector<Complex>&) const;
#endif
}
