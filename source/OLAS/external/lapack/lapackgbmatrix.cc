// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

/* $Id$ */

#include <cmath>
#include <complex>
#include "external/lapack/lapackgbmatrix.hh"
#include "utils/utils.hh"

// the following headers are required for Export()
#include <cstdio>

namespace OLAS {

  // **********************************************
  //   convert a CRS_Matrix into a LapackGBMatrix
  // **********************************************
  template <class entryF, class entryC>
  void LapackGBMatrix<entryF,entryC>::Convert( const CRS_Matrix<entryC>& A ) {

    ENTER_FCN( "LapackGBMatrix::Convert(CRS_Matrix)" );
  
    // check the dimensions
    if( A.GetNrows() != GetNrows() || A.GetNcols() != GetNcols() ) {
      Error( "The matrix was created with dimensions different from "
             "those of the passed matrix", __FILE__, __LINE__ );
    }

    // set the sparsity pattern analoguously to SetSparsityPattern
    if( data_ != NULL ) {
      Error( "Matrix structure has already been defined!",
             __FILE__, __LINE__ );
    }

    // first obtain bandwidth by hand from the matrix
    const Integer* const pRow = A.GetRowPointer(); // pointer to row indices
    const Integer* const pCol = A.GetColPointer(); // pointer to column indices
    const entryC*  const pDat = A.GetDataPointer();// pointer to the data array
    UInt bw; // variable for temporary storage of the bandwidth of one row
    wlower_ = wupper_ = 0;
    for( UInt i = 1; i <= GetNrows(); i++ ) {
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
    NewArray( data_, entryF, length_ );

    // set matrix entries to zero
    Init();

    // fill the matrix (there might be a more performant way, but
    // I implement it that way first ... in a braver hour differently)
    entryC tempVar = 0;
    for( Integer i = 1; i <= GetNrows(); i++ ) {
      for( Integer ij = pRow[i]; ij < pRow[i+1]; ij++ ) {
	tempVar = (entryC)pDat[ij];
        SetMatrixEntry( i, (Integer)pCol[ij], tempVar );
      }
    }

#ifdef DEBUG_LAPACKGBMATRIX
    (*debug) << " LapackGBMatrix: nrows    = " << nrows_      << std::endl
             << "                 ncols    = " << ncols_      << std::endl
             << "                 wlower   = " << wlower_     << std::endl
             << "                 wupper   = " << wupper_     << std::endl
             << "                 nrowsact = " << nrowsact_   << std::endl
             << "                 length   = " << length_     << std::endl
             << "                 enlarged = " << amEnlarged_ << std::endl;
#endif
  }


  // **************************
  //   Setup sparsity pattern
  // **************************
  template <class entryF, class entryC>
  void LapackGBMatrix<entryF,entryC>::SetSparsityPattern( BaseGraph &graph ) {
    ENTER_FCN( "LapackGBMatrix::SetSparsityPattern" );

    // Check that matrix has not already been defined
    if ( data_ != NULL ) {
      Error( "Matrix structure has already been defined!", __FILE__,
             __LINE__ );
    }

    // obtain bandwidth info from graph
    graph.GetBandwidth( wlower_, wupper_ );

    // we do not have additional storage space
    amEnlarged_ = false;
    offset_ = 0;

    // compute number of columns in packed matrix
    nrowsact_ = ( wlower_ + 1 + wupper_ );

    // compute amount of storage for matrix
    length_ = nrowsact_ * ncols_;

#ifdef DEBUG_LAPACKGBMATRIX
    (*debug) << " LapackGBMatrix: nrows    = " << nrows_      << std::endl
             << "                 ncols    = " << ncols_      << std::endl
             << "                 wlower   = " << wlower_     << std::endl
             << "                 wupper   = " << wupper_     << std::endl
             << "                 nrowsact = " << nrowsact_   << std::endl
             << "                 length   = " << length_     << std::endl
             << "                 enlarged = " << amEnlarged_ << std::endl;
#endif
    
    // allocate memory
    NewArray( data_, entryF, length_ );

    // set matrix entries to zero
    Init();

  }


  // ***************************
  //   Alternative Constructor
  // ***************************
  template <class entryF, class entryC>
  LapackGBMatrix<entryF,entryC>::LapackGBMatrix( Integer nrows,
						 Integer ncols,
						 Integer wlower,
						 Integer wupper,
						 MatrixEntryType etype ) {

    ENTER_FCN( "LapackGBMatrix::LapackGBMatrix" );

    // Ensure that matrix is square
    if ( nrows != ncols ) {
      Error( "Request for a non-square LapackGBMatrix", __FILE__, __LINE__ );
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

#ifdef DEBUG_LAPACKGBMATRIX
    (*debug) << " LapackGBMatrix: nrows    = " << nrows_      << std::endl
	     << "                 ncols    = " << ncols_      << std::endl
	     << "                 wlower   = " << wlower_     << std::endl
	     << "                 wupper   = " << wupper_     << std::endl
	     << "                 nrowsact = " << nrowsact_   << std::endl
	     << "                 length   = " << length_     << std::endl
	     << "                 enlarged = " << amEnlarged_ << std::endl;
#endif

    // allocate memory
    NewArray( data_, entryF, length_ );

    // set matrix entries to zero
    Init();

  }


  // *****************************
  //   Initialise matrix to zero
  // *****************************
  template <class entryF, class entryC>
  void LapackGBMatrix<entryF,entryC>::Init() {
    ENTER_FCN( "LapackGBMatrix::Init" );
    entryF fzero;
    entryC czero = 0;
    CC2F77( czero, fzero );
    for ( Integer i = 1; i <= length_; i++ ) {
      data_[i] = fzero;
    }
  }


  // ********************************************
  //   Initialise matrix to zero (F77complex16)
  // ********************************************
  template <>
  void LapackGBMatrix<F77complex16,std::complex<double> >::Init() {
    ENTER_FCN( "LapackGBMatrix::Init" );
    for ( Integer i = 1; i <= length_; i++ ) {
      data_[i].real = 0.0;
      data_[i].imag = 0.0;
    }
  }


  // *******************************
  //   Set value of a matrix entry
  // *******************************
  template <class entryF, class entryC>
  void LapackGBMatrix<entryF,entryC>::SetMatrixEntry( Integer i, Integer j,
						      entryC &v ) {
    ENTER_IFCN( "LapackGBMatrix::SetMatrixEntry" );

#ifdef DEBUG_LAPACKGBMATRIX
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
#endif

    entryF val;
    CC2F77( v, val );
    data_[Index(i,j)] = val;
  }


  // *******************************
  //   Add value to a matrix entry
  // *******************************
  template <class entryF, class entryC>
  void LapackGBMatrix<entryF,entryC>::AddToMatrixEntry( Integer i, Integer j,
							entryC &v ) {
    ENTER_IFCN( "LapackGBMatrix::AddToMatrixEntry" );

#ifdef DEBUG_LAPACKGBMATRIX
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
#endif
    
    entryF val;
    CC2F77( v, val );
    data_[Index(i,j)] += val;
  }


  // *************************
  //   Export matrix to file
  // *************************
  template <class entryF, class entryC>
  void LapackGBMatrix<entryF,entryC>::Export( const Char *fname,
					      const Char *comment ) const {

    ENTER_FCN( "LapackGBMatrix::Export" );

    // open output file and check for errors
    FILE *fp = fopen( fname, "w" );
    if ( fp == NULL ) {
      Char *errmsg;
      NewArray( errmsg, Char, strlen(fname)+40 );
      sprintf( errmsg, "Cannot open file %s for writing!", fname );
      Error( errmsg, __FILE__, __LINE__ );
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
      Error( "Unexpected matrix entry type", __FILE__, __LINE__ );
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
      Char *errmsg;
      NewArray( errmsg, Char, strlen(fname)+40 );
      sprintf( errmsg, "Could not close file %s after writing!", fname );
      Warning( errmsg, __FILE__, __LINE__ );
    }
  }


  // **********************************
  //   Write matrix entries (F77real8)
  // **********************************
  template <class entryF, class entryC>
  void LapackGBMatrix<entryF,entryC>::WriteEntries( FILE *fp ) const {

    ENTER_FCN( "LapackGBMatrix::WriteEntries" );

    int i, j;
    int colinit, colstop;

    // Count number of non-zero matrix entries
    Integer nnz = 0;
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

#ifdef DEBUG_LAPACKGBMATRIX
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
#endif

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

    ENTER_FCN( "LapackGBMatrix::WriteEntries" );

    int i, j;
    int colinit, colstop;

    // Count number of non-zero matrix entries
    Integer nnz = 0;
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

    ENTER_FCN( "LapackGBMatrix::Add" );

    // Down-cast to LapackGBMatrix
    TRY_CAST {
      const LapackGBMatrix<entryF,entryC> &lpmat = dynamic_cast<const
	LapackGBMatrix<entryF,entryC>&>(mat);

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
      for ( Integer i = 1; i <= length_; i++ ) {
	data_[i] += factor * lpmat.data_[i];
      }

    } CATCH_CAST;

  }


  // ***************************************************
  //   Add a multiple of another matrix (F77complex16)
  // ***************************************************
  template <>
  void LapackGBMatrix<F77complex16,std::complex<double> >::
  Add( const Double factor, const StdMatrix& mat ) {

    ENTER_FCN( "LapackGBMatrix::Add" );

    // Down-cast to LapackGBMatrix
    TRY_CAST {
      const LapackGBMatrix<F77complex16,std::complex<double> > &lpmat =
	dynamic_cast<const LapackGBMatrix<F77complex16,std::complex<double> >&>
	(mat);

      // Perform the scaled addition
      for ( Integer i = 1; i <= length_; i++ ) {
	data_[i].real += factor * lpmat.data_[i].real;
	data_[i].imag += factor * lpmat.data_[i].imag;
      }

    } CATCH_CAST;

  }


  // ******************************
  //   Get Maximal Diagonal Entry
  // ******************************
  template <class entryF, class entryC>
  Double LapackGBMatrix<entryF,entryC>::GetMaxDiag() const {

    ENTER_FCN( "LapackGBMatrix::GetMaxDiag" );

    Double absmax = 0.0;
    entryF fdata;
    entryC cdata;
    for ( Integer i = 1; i <= ncols_; i++ ) {
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
  void LapackGBMatrix<entryF,entryC>::SetDiagEntry( Integer i, entryC &v ){
    ENTER_FCN( "LapackGBMatrix::SetDiagEntry" );
    entryF fdata;
    CC2F77( v, fdata );
    data_[Index(i,i)] = fdata;
  }


  // ************************
  //   Get a Diagonal Entry
  // ************************
  template <class entryF, class entryC>
  void LapackGBMatrix<entryF,entryC>::GetDiagEntry( Integer i,
						    entryC &v ) const {
    ENTER_FCN( "LapackGBMatrix::GetDiagEntry" );
    F772CC( data_[Index(i,i)], v );
  }


  // ***********************
  //  Forced Instantiation
  // ***********************
  template <class entryF, class entryC>
  void LapackGBMatrix<entryF,entryC>::InstantiatePublicMethods() {

    ENTER_FCN( "LapackGBMatrix::InstantiatePublicMethods" );

    Error( "This function should never be called", __FILE__, __LINE__ );

    CRS_Matrix<entryC> dummyCRSMat;
    StdMatrix *dummyMat = NULL;
    BaseGraph dummyGraph( 0, 0, NOREORDERING );
    Integer i = 0, j = 0;
    Double val = 0;
    entryC dummyEntry;
    MatrixStorageType dummyType = NOSTORAGETYPE;
    MatrixEntryType dummyEntryType = NOENTRYTYPE;

    LapackGBMatrix( 1, 1, 2, 2, dummyEntryType );
    Export( "", "" );
    Convert( dummyCRSMat );
    SetSparsityPattern( dummyGraph );
    Init();
    SetMatrixEntry( i, j, dummyEntry );
    AddToMatrixEntry( i, j, dummyEntry );
    Add( val, *dummyMat );
    SetDiagEntry( i, dummyEntry );
    GetDiagEntry( i, dummyEntry );

    dummyType = GetStorageType();
    dummyEntryType = GetEntryType();
    i = GetUpperBandwidth();
    i = GetLowerBandwidth();
    val = GetMaxDiag();
    entryF *dummyArray1 = GetDataPointer0();
    const entryF *dummyArray2 = GetDataPointer0();

  }

#ifdef __sgi
#pragma do_not_instantiate LapackGBMatrix< F77complex16, std::complex<float> >
#endif

}

