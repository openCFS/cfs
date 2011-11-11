#include <cmath>
#include <complex>
#include <string.h>
#include "LapackGBMatrix.hh"

// the following headers are required for Export()
#include <cstdio>

namespace CoupledField {

  // **********************************************
  //   convert a CRS_Matrix into a LapackGBMatrix
  // **********************************************
  template <class entryF, class entryC>
  void LapackGBMatrix<entryF,entryC>::Convert( const CRS_Matrix<entryC>& A ) {

  
    // check the dimensions
    if( A.GetNumRows() != GetNumRows() || A.GetNumCols() != GetNumCols() ) {
      EXCEPTION( "The matrix was created with dimensions different from "
                 << "those of the passed matrix" );
    }

    // set the sparsity pattern analoguously to SetSparsityPattern
    if( data_ != NULL ) {
      EXCEPTION( "Matrix structure has already been defined!" );
    }

    // first obtain bandwidth by hand from the matrix
    const UInt* const pRow = A.GetRowPointer(); // pointer to row indices
    const UInt* const pCol = A.GetColPointer(); // pointer to column indices
    const entryC*  const pDat = A.GetDataPointer();// pointer to the data array
    UInt bw; // variable for temporary storage of the bandwidth of one row
    wlower_ = wupper_ = 0;
    for( UInt i = 1; i <= GetNumRows(); i++ ) {
      // only compute bandwidth, if offdiagonal entries are present
      if( pRow[i] < pRow[i+1] ) {
        // lower bandwidth
        if( i > pCol[pRow[i]+1]   &&
            wlower_ < (bw = i - pCol[pRow[i]+1]  ) )  wlower_ = bw;
        // upper bandwidth
        if( i < pCol[pRow[i+1]-1] &&
            wupper_ < (bw = pCol[pRow[i+1]-1] - i) )  wupper_ = bw;
      }
    }

    // -> followed by the code in SetSparsityPattern
    // we do not have additional storage space
    amEnlarged_ = false;
    offset_     = 0;

    // compute number of columns in packed matrix
    nrowsact_ = ( wlower_ + 1 + wupper_ );

    // compute amount of storage for matrix
    length_ = nrowsact_ * ncols_;

    // allocate memory
    NEWARRAY( data_, entryF, length_ );

    // set matrix entries to zero
    Init();

    // fill the matrix (there might be a more performant way, but
    // I implement it that way first ... in a braver hour differently)
    entryC tempVar = 0;
    for( UInt i = 1; i <= GetNumRows(); i++ ) {
      for( UInt ij = pRow[i]; ij < pRow[i+1]; ij++ ) {
        tempVar = (entryC)pDat[ij];
        SetMatrixEntry( i, (UInt)pCol[ij], tempVar );
      }
    }

    /*
    (*debug) << " LapackGBMatrix: nrows    = " << nrows_      << std::endl
             << "                 ncols    = " << ncols_      << std::endl
             << "                 wlower   = " << wlower_     << std::endl
             << "                 wupper   = " << wupper_     << std::endl
             << "                 nrowsact = " << nrowsact_   << std::endl
             << "                 length   = " << length_     << std::endl
             << "                 enlarged = " << amEnlarged_ << std::endl;
    */
  }


  // **************************
  //   Setup sparsity pattern
  // **************************
  template <class entryF, class entryC>
  void LapackGBMatrix<entryF,entryC>::SetSparsityPattern( BaseGraph &graph ) {

    // Check that matrix has not already been defined
    if ( data_ != NULL ) {
      EXCEPTION( "Matrix structure has already been defined!" );
    }

    // obtain bandwidth info from graph
    UInt avgBw; // not really used
    graph.GetBandwidth( wlower_, wupper_, avgBw );

    // we do not have additional storage space
    amEnlarged_ = false;
    offset_ = 0;

    // compute number of columns in packed matrix
    nrowsact_ = ( wlower_ + 1 + wupper_ );

    // compute amount of storage for matrix
    length_ = nrowsact_ * ncols_;

    /*
    (*debug) << " LapackGBMatrix: nrows    = " << nrows_      << std::endl
             << "                 ncols    = " << ncols_      << std::endl
             << "                 wlower   = " << wlower_     << std::endl
             << "                 wupper   = " << wupper_     << std::endl
             << "                 nrowsact = " << nrowsact_   << std::endl
             << "                 length   = " << length_     << std::endl
             << "                 enlarged = " << amEnlarged_ << std::endl;
    */    
    // allocate memory
    NEWARRAY( data_, entryF, length_ );

    // set matrix entries to zero
    Init();

  }


  // ***************************
  //   Alternative Constructor
  // ***************************
  template <class entryF, class entryC>
  LapackGBMatrix<entryF,entryC>::LapackGBMatrix( UInt nrows,
                                                 UInt ncols,
                                                 UInt wlower,
                                                 UInt wupper,
                                                 BaseMatrix::EntryType etype ) {


    // Ensure that matrix is square
    if ( nrows != ncols ) {
      EXCEPTION( "Request for a non-square LapackGBMatrix" );
    }

    // Set matrix information
    nrows_ = nrows;
    ncols_ = ncols;
    wupper_ = wupper;
    wlower_ = wlower;
    myEntryType_ = etype;
    data_ = NULL;

    // compute amount of storage for matrix and number of rows in packed
    // banded matrix, take work space into account in doing so
    amEnlarged_ = true;
    nrowsact_ = ( 2 * wlower_ + 1 + wupper_ );
    length_ = nrowsact_ * ncols_;
    offset_ = wlower_;

    /*
    (*debug) << " LapackGBMatrix: nrows    = " << nrows_      << std::endl
    << "                 ncols    = " << ncols_      << std::endl
    << "                 wlower   = " << wlower_     << std::endl
    << "                 wupper   = " << wupper_     << std::endl
    << "                 nrowsact = " << nrowsact_   << std::endl
    << "                 length   = " << length_     << std::endl
    << "                 enlarged = " << amEnlarged_ << std::endl;
    */

    // allocate memory
    NEWARRAY( data_, entryF, length_ );

    // set matrix entries to zero
    Init();

  }


  // *****************************
  //   Initialise matrix to zero
  // *****************************
  template <class entryF, class entryC>
  void LapackGBMatrix<entryF,entryC>::Init() {
    entryF fzero;
    entryC czero = 0;
    CC2F77( czero, fzero );
    for ( UInt i = 1; i <= length_; i++ ) {
      data_[i] = fzero;
    }
  }


  // ********************************************
  //   Initialise matrix to zero (F77complex16)
  // ********************************************
  template <>
  void LapackGBMatrix<F77complex16,std::complex<double> >::Init() {
    for ( UInt i = 1; i <= length_; i++ ) {
      data_[i].real = 0.0;
      data_[i].imag = 0.0;
    }
  }


  // *******************************
  //   Set value of a matrix entry
  // *******************************
  template <class entryF, class entryC>
  void LapackGBMatrix<entryF,entryC>::SetMatrixEntry( UInt i, UInt j,
                                                      entryC &v ) {

    /*
    if ( Index(i,j) > length_ || Index(i,j) < 1 ) {
      (*debug) << " FUBAR: Index(i,j) is out of bounds!" << std::endl
      << "         (i,j) = (" << i << ", " << j << ")" << std::endl
      << "    Index(i,j) = " << Index(i,j) << std::endl
      << "        nrows  = " << nrows_  << std::endl
      << "        ncols  = " << ncols_  << std::endl
      << "        wlower = " << wlower_ << std::endl
      << "        wupper = " << wupper_ << std::endl
      << "        length = " << length_ << std::endl;
    }
    */

    entryF val;
    CC2F77( v, val );
    data_[Index(i,j)] = val;
  }


  // *******************************
  //   Add value to a matrix entry
  // *******************************
  template <class entryF, class entryC>
  void LapackGBMatrix<entryF,entryC>::AddToMatrixEntry( UInt i, UInt j,
                                                        entryC &v ) {

    /* 
    if ( Index(i,j) > length_ || Index(i,j) < 1 ) {
      (*debug) << " FUBAR: Index(i,j) is out of bounds!" << std::endl
      << "         (i,j) = (" << i << ", " << j << ")" << std::endl
      << "    Index(i,j) = " << Index(i,j) << std::endl
      << "        nrows  = " << nrows_  << std::endl
      << "        ncols  = " << ncols_  << std::endl
      << "        wlower = " << wlower_ << std::endl
      << "        wupper = " << wupper_ << std::endl
      << "        length = " << length_ << std::endl;
    }
    */

    entryF val;
    CC2F77( v, val );
    data_[Index(i,j)] += val;
  }


  // *************************
  //   Export matrix to file
  // *************************
  template <class entryF, class entryC>
  void LapackGBMatrix<entryF,entryC>::Export( const char *fname,
                                              const char *comment ) const {


    // open output file and check for errors
    FILE *fp = fopen( fname, "w" );
    if ( fp == NULL ) {
      EXCEPTION( "Cannot open file " << fname << " for writing!" );
    }

    // ---------------------
    //   Write file header
    // ---------------------

    // Matrix Market Format Specification
    if ( GetEntryType() == F77REAL8 ) {
      fprintf( fp, "%%%%MatrixMarket matrix coordinate real general\n" );
    }
    else if ( GetEntryType() == F77COMPLEX16 ) {
      fprintf( fp, "%%%%MatrixMarket matrix coordinate complex general\n" );
    }
    else {
      std::cerr << "Unexpected matrix entry type " << GetEntryType()
      << std::endl;
      EXCEPTION( "Unexpected matrix entry type" );
    }

    // User-supplied private comment
    if ( comment != NULL ) {
      fprintf( fp, "%%\n%% %s\n%%\n", comment );
    }
    else {
      fprintf( fp, "%%\n%% Matrix exported by OLAS\n%%\n" );
    }

    // -----------------
    //   Write entries
    // -----------------
    WriteEntries( fp );

    // close output file
    if ( fclose( fp ) == EOF ) {
      WARN( "Could not close file " << fname << " after writing!" );
    }
  }


  // **********************************
  //   Write matrix entries (F77real8)
  // **********************************
  template <class entryF, class entryC>
  void LapackGBMatrix<entryF,entryC>::WriteEntries( FILE *fp ) const {


    UInt i, j;
    UInt colinit, colstop;

    // Count number of non-zero matrix entries
    UInt nnz = 0;
    for ( j = 1; j <= ncols_; j++ ) {

      // compute column bounds
      colinit = (j - wupper_) < 1 ? 1 : (j - wupper_);
      colstop = (j + wlower_) > nrows_ ? nrows_ : (j + wlower_);

      // loop over column entries in band
      for ( i = colinit; i <= colstop; i++ ) {

        if ( data_[Index(i,j)] != static_cast<entryF>(0.0) ) {
          nnz++;
        }
      }
    }

    // Write information on number of rows, columns and entries
    fprintf( fp, "%d\t%d\t%d\n", nrows_, ncols_, nnz );

    // Write non-zero entries
    F77real8 val;
    for ( j = 1; j <= ncols_; j++ ) {

      // compute column bounds
      colinit = (j - wupper_) < 1 ? 1 : (j - wupper_);
      colstop = (j + wlower_) > nrows_ ? nrows_ : (j + wlower_);

      // loop over column entries in band
      for ( i = colinit; i <= colstop; i++ ) {

    /*
        if ( Index(i,j) > length_ || Index(i,j) < 1 ) {
          (*debug) << " FUBAR: Index(i,j) is out of bounds!" << std::endl
          << "         (i,j) = (" << i << ", " << j << ")" << std::endl
          << "    Index(i,j) = " << Index(i,j) << std::endl
          << "        nrows  = " << nrows_  << std::endl
          << "        ncols  = " << ncols_  << std::endl
          << "        wlower = " << wlower_ << std::endl
          << "        wupper = " << wupper_ << std::endl
          << "        length = " << length_ << std::endl;
        }
    */

        val = data_[Index(i,j)];
        if ( val != 0.0 ) {
          fprintf( fp, "%6d\t%6d\t% 22.16e\n", i, j, val );
        }
      }
    }
  }


  // ***************************************
  //   Write matrix entries (F77COMPLEX16)
  // ***************************************
  template <>
  void LapackGBMatrix<F77complex16,std::complex<double> >::
  WriteEntries(FILE *fp) const {


    UInt i, j;
    UInt colinit, colstop;

    // Count number of non-zero matrix entries
    UInt nnz = 0;
    for ( j = 1; j <= ncols_; j++ ) {

      // compute column bounds
      colinit = (j - wupper_) < 1 ? 1 : (j - wupper_);
      colstop = (j + wlower_) > nrows_ ? nrows_ : (j + wlower_);

      // loop over column entries in band
      for ( i = colinit; i <= colstop; i++ ) {

        if ( data_[Index(i,j)].real != 0 || data_[Index(i,j)].imag != 0 ) {
          nnz++;
        }
      }
    }

    // Write information on number of rows, columns and entries
    fprintf( fp, "%d\t%d\t%d\n", nrows_, ncols_, nnz );

    // Write non-zero entries
    F77complex16 val;
    for ( j = 1; j <= ncols_; j++ ) {

      // compute column bounds
      colinit = (j - wupper_) < 1 ? 1 : (j - wupper_);
      colstop = (j + wlower_) > nrows_ ? nrows_ : (j + wlower_);

      // loop over column entries in band
      for ( i = colinit; i <= colstop; i++ ) {
        val = data_[Index(i,j)];
        if ( val.real != 0 || val.imag != 0 ) {
          fprintf( fp, "%6d\t%6d\t% 22.16e\t% 22.16e\n", i, j, val.real,
                   val.imag );
        }
      }
    }
  }


  // ***********************************************
  //   Add a multiple of another matrix (F77real8)
  // ***********************************************
  template <class entryF, class entryC>
  void LapackGBMatrix<entryF,entryC>::Add( const Double factor,
                                           const StdMatrix& mat ) {


    // Down-cast to LapackGBMatrix
      const LapackGBMatrix<entryF,entryC> &lpmat = dynamic_cast<const LapackGBMatrix<entryF,entryC>&>(mat);

#ifdef DEBUG_LAPACKGBMATRIX
      if ( length_ != lpmat.length_ ) {
        (*debug) << " LapackGBMatrix::Add - Inconsistency in matrix structure!"
        << " Matrix A: length of data vector = " << length_
        << std::endl
        << " Matrix B: length of data vector = " << lpmat.length_
        << std::endl;
        std::cerr << "LapackGBMatrix::Add - Inconsistency in matrix structure!"
        << "Matrix A: length of data vector = " << length_
        << std::endl
        << "Matrix B: length of data vector = " << lpmat.length_
        << std::endl;
        exit(1);
      }
#endif

      // Perform the scaled addition
      for ( UInt i = 1; i <= length_; i++ ) {
        data_[i] += factor * lpmat.data_[i];
      }

  }


  // ***************************************************
  //   Add a multiple of another matrix (F77complex16)
  // ***************************************************
  template <>
  void LapackGBMatrix<F77complex16,std::complex<double> >::
  Add( const Double factor, const StdMatrix& mat ) {
    // Down-cast to LapackGBMatrix
    const LapackGBMatrix<F77complex16,std::complex<double> > &lpmat =
      dynamic_cast<const LapackGBMatrix<F77complex16,std::complex<double> >&>(mat);

      // Perform the scaled addition
      for ( UInt i = 1; i <= length_; i++ ) {
        data_[i].real += factor * lpmat.data_[i].real;
        data_[i].imag += factor * lpmat.data_[i].imag;
      }
  }


  // ******************************
  //   Get Maximal Diagonal Entry
  // ******************************
  template <class entryF, class entryC>
  Double LapackGBMatrix<entryF,entryC>::GetMaxDiag() const {


    Double absmax = 0.0;
    entryF fdata;
    entryC cdata;
    for ( UInt i = 1; i <= ncols_; i++ ) {
      fdata = data_[Index(i,i)];
      F772CC( fdata, cdata );
      absmax = std::abs(cdata) > absmax ? std::abs(cdata) : absmax;
    }

    return absmax;
  }


  // ************************
  //   Set a Diagonal Entry
  // ************************
  template <class entryF, class entryC>
  void LapackGBMatrix<entryF,entryC>::SetDiagEntry( UInt i, entryC &v ){
    entryF fdata;
    CC2F77( v, fdata );
    data_[Index(i,i)] = fdata;
  }


  // ************************
  //   Get a Diagonal Entry
  // ************************
  template <class entryF, class entryC>
  void LapackGBMatrix<entryF,entryC>::GetDiagEntry( UInt i,
                                                    entryC &v ) const {
    F772CC( data_[Index(i,i)], v );
  }


  // Explicit template instantiation
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template class LapackGBMatrix< F77real8, Double >;
  template class LapackGBMatrix< F77complex16, Complex >;
#endif
}

