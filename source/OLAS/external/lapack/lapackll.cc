#include "utils/utils.hh"
#include "external/lapack/lapackll.hh"
#include "external/lapack/olasf77mapping.hh"

namespace OLAS {

  // ***********************
  //   Default Constructor
  // ***********************
  Lapack_LL::Lapack_LL() {
    ENTER_FCN( "Lapack_LL::Lapack_LL" );
    Error( "Default constructor of LAPACK_LL is forbidden!",
	   __FILE__, __LINE__ );
  }


  // *****************************
  //   Constructor (for factory)
  // *****************************
  Lapack_LL::Lapack_LL( OLAS_Params *myParams, OLAS_Report *myReport ) {

    ENTER_FCN( "Lapack_LL::Lapack_LL" );

    // Set pointers to communication objects
    myParams_ = myParams;
    myReport_ = myReport;

    // Initialise remaining attributes
    facMat_            = NULL;
    lapackRHS_         = NULL;
    lapackRHSCapacity_ = 0;
    sysMatBW_          = 0;
    facMatCapacity_    = 0;
    facMatEntries_     = 0;
    amFactorised_      = false;

  }


  // ********************
  //   Deep Constructor
  // ********************
  Lapack_LL::~Lapack_LL() {
    ENTER_FCN( "Lapack_LL::~Lapack_LL" );
    DeleteArray( lapackRHS_ );
    DeleteArray( facMat_ );
  }


  // ********************************
  //   Setup: Perform Factorisation
  // ********************************
  void Lapack_LL::Setup( BaseMatrix &sysMat ) {

    ENTER_FCN( "Lapack_LL::Setup" );

    // Are we expected to be verbose?
    bool logging = myParams_->GetBoolValue( "LAPACKLL_logging" );

    // Report to logfile
    if ( logging == true ) {
      (*cla) << " --------------------------------------\n"
	     << " LAPACK_LL: Starting factorisation of a "
	     << sysMat.GetNrows() << " x " << sysMat.GetNcols()
	     << " matrix" << std::endl;
    }

    TRY_CAST {

      // Down-cast matrix to standard matrix
      RefCast( sysMat, StdMatrix, stdMat );

      // Check that we have the correct matrix type
      MatrixStorageType mtype = stdMat.GetStorageType();
      if ( mtype != SPARSE_SYM ) {
	(*error) << "Lapack_LL: Expected a sparseSym matrix, but got a "
		 << Enum2String( mtype ) << " matrix";
	Error( __FILE__, __LINE__ );
      }

      // Check the block size of the matrix entries
      Integer blockSize = stdMat.GetBlockSize();
      if ( blockSize != 1 ) {
	(*error) << "Lapack_LL: Block size of matrix entries is " << blockSize
		 << " but should be one";
	Error( __FILE__, __LINE__ );
      }

      // Get the entry type to figure out which Factorisation method to call
      MatrixEntryType etype = stdMat.GetEntryType();

      // Call appropriate factorisation routine
      switch( etype ) {

      case DOUBLE:
	TRY_CAST {
	  FactoriseReal( stdMat );
	} CATCH_CAST;
	break;

      case COMPLEX:
	TRY_CAST {
	  FactoriseComplex( stdMat );
	} CATCH_CAST;
	break;

      default:
	(*error) << "Matrix entry type is neither real nor complex!";
	Error( __FILE__, __LINE__ );
      }

    } CATCH_CAST;

    // now we have a (new) factorisation
    amFactorised_ = true;

    // Report to logfile
    if ( logging == true ) {
      (*cla) << " LAPACK_LL: Finished factorisation\n"
	     << " --------------------------------------" << std::endl;
    }
  }


  // *****************
  //   FactoriseReal
  // *****************
  void Lapack_LL::FactoriseReal( StdMatrix &stdMat ) {

    ENTER_FCN( "Lapack_LL::FactoriseReal" );

    UInt i, k;

    // Are we expected to be verbose?
    bool logging = myParams_->GetBoolValue( "LAPACKLL_logging" );

    // Initialise parameters for DPBTRF
    char lp_uplo = 'L';
    int  lp_n    = stdMat.GetNcols();
    int  lp_kd   = 0;
    int  lp_ldab = 0;
    int  lp_info = 0;

    // Down-cast matrix to special format
    RefCast( stdMat, SCRS_Matrix<Double>, scrsMat );

    // Get hold of column index array
    const Integer *cidx = scrsMat.GetColPointer();

    // Get hold of row pointer index array
    const Integer *rptr = scrsMat.GetRowPointer();

    // Get hold of data array
    const Double *data = scrsMat.GetDataPointer();

    // Compute number of stored matrix entries
    Integer nEntries = ( scrsMat.GetNnz() + lp_n ) / 2;

    // Test, if we must check for a new matrix' bandwidth and
    // maybe also allocate new memory
    if ( amFactorised_ == false || facMat_ == NULL ||
	 myParams_->GetBoolValue( "newMatrixPattern" ) == true ) {

#ifdef DEBUG_LAPACK_LL

      // Are we over-writing an existing factorisation? Check, whether the
      // two indicators agree, if not cry out loud!
      if ( (facmat_+1) != NULL || amFactorised_ == true ) {
	if ( (facmat_+1) == NULL && amFactorised_ == true ) {
	  (*error) << "Lapack_LL: Internal error. facmat_ = NULL, but "
		   << "amFactorised_ = true!";
	  Error( __FILE__, __LINE__ );
	}
	else if ( amFactorised_ == false ) {
	  (*error) << "Lapack_LL: Internal error. facmat_ <> NULL, but "
		   << "amFactorised_ = false!";
	  Error( __FILE__, __LINE__ );
	}
      }

#endif

      // Determine bandwidth of the matrix

      // Since the entries in a row are stored in ascending order
      // we only have to compute the distance between the column index
      // of the final entries in each row and the row index to find
      // the bandwidth. As always we assume that there is a diagonal entry,
      // so no row is empty.
      UInt bwGlobal = 0;
      UInt bwLocal  = 0;
      for ( i = 1; i <= lp_n; i++ ) {

	// Determine bandwidth of this row
	bwLocal = (UInt)(cidx[ rptr[i+1] - 1 ] - i);

	// Compare to global bandwidth
	bwGlobal = bwLocal > bwGlobal ? bwLocal : bwGlobal;

      }

      // Save bandwidth of the system matrix
      sysMatBW_ = bwGlobal;

      // Determine number of entries in band matrix
      facMatEntries_ = ( bwGlobal + 1 ) * lp_n;

      // Report to logfile
      if ( logging == true ) {
	(*cla) << "    Bandwidth = " << bwGlobal << "\n"
	       << "    Sparse matrix has " << scrsMat.GetNnz()
	       << " non-zero entries\n"
	       << "    Banded matrix has " << facMatEntries_ << " real entries"
	       << std::endl;
      }

      // Memory allocation for band matrix: We only perform a
      // re-allocation, if the new size requirements are larger
      // than the old ones
      if ( amFactorised_ == false ) {
	NewArray( facMat_, Double, facMatEntries_ );
      }
      else {
	if ( facMatEntries_ > facMatCapacity_ ) {
	  DeleteArray( facMat_ );
	  NewArray( facMat_, Double, facMatEntries_ );
	  facMatCapacity_ = facMatEntries_;
	}

#ifdef DEBUG_LAPACK_LL
	else {
	  for ( i = facMatEntries_ + 1; i <= facMatCapacity_; i++ ) {
	    facMat_[i] = 0.0;
	  }
	}
#endif

      }
    }

    // Set all entries in band matrix to zero
    for ( i = 1; i <= facMatEntries_; i++ ) {
      facMat_[i] = 0.0;
    }

    // Now we create a symmetric LAPACK band matrix. The SCRS_Matrix stores the
    // upper triangle of the matrix. If we directly insert the rows into the
    // facMat_ array, this gives us in FORTRAN77 format the lower triangle of
    // the matrix, which is fine for LAPACK and is hopefully more cache
    // efficient than converting the rows into FORTRAN format.
    UInt colInit = 0;
    for ( i = 1; i <= lp_n; i++ ) {
      for ( k = rptr[i]; k < rptr[i+1]; k++ ) {
	facMat_[ colInit + cidx[k] - i + 1 ] = data[k];
      }
      colInit += (sysMatBW_ + 1);
    }

    // Prepare parameters for LAPACK
    lp_kd = sysMatBW_;
    lp_ldab = sysMatBW_ + 1;

    // Now call the LAPACK routine to compute Cholesky factorisation
    LP_DPBTRF( &lp_uplo, &lp_n, &lp_kd, facMat_ + 1, &lp_ldab, &lp_info );

    // Process return status of LAPACK (i.e. if it is compiled to return
    // on error)
    if ( lp_info < 0 ) {
      (*error) << "Lapack_LL: DPBTRF reports input value no. " << -lp_info
	       << " is invalid";
      Error( __FILE__, __LINE__ );
    }
    else if ( lp_info > 0 ) {
      (*error) << "Lapack_LL: DPBTRF reports that leading minor of order "
	       << lp_info << " is not positive definite and the factorization"
	       << " could not be completed";
      Error( __FILE__, __LINE__ );
    }
  }


  // ********************
  //   FactoriseComplex
  // ********************
  void Lapack_LL::FactoriseComplex( StdMatrix &stdMat ) {

    ENTER_FCN( "Lapack_LL::FactoriseComplex" );

    UInt i, k;

    // Are we expected to be verbose?
    bool logging = myParams_->GetBoolValue( "LAPACKLL_logging" );

    // Initialise parameters for ZPBTRF
    char lp_uplo = 'L';
    int  lp_n    = stdMat.GetNcols();
    int  lp_kd   = 0;
    int  lp_ldab = 0;
    int  lp_info = 0;

    // Down-cast matrix to special format
    RefCast( stdMat, SCRS_Matrix<Complex>, scrsMat );

    // Get hold of column index array
    const Integer *cidx = scrsMat.GetColPointer();

    // Get hold of row pointer index array
    const Integer *rptr = scrsMat.GetRowPointer();

    // Get hold of data array
    const Complex *data = scrsMat.GetDataPointer();

    // Compute number of stored matrix entries
    Integer nEntries = ( scrsMat.GetNnz() + lp_n ) / 2;

    // Test, if we must check for a new matrix' bandwidth and
    // maybe also allocate new memory
    if ( amFactorised_ == false || facMat_ == NULL ||
	 myParams_->GetBoolValue( "newMatrixPattern" ) == true ) {

#ifdef DEBUG_LAPACK_LL

      // Are we over-writing an existing factorisation? Check, whether the
      // two indicators agree, if not cry out loud!
      if ( (facmat_+1) != NULL || amFactorised_ == true ) {
	if ( (facmat_+1) == NULL && amFactorised_ == true ) {
	  (*error) << "Lapack_LL: Internal error. facmat_ = NULL, but "
		   << "amFactorised_ = true!";
	  Error( __FILE__, __LINE__ );
	}
	else if ( amFactorised_ == false ) {
	  (*error) << "Lapack_LL: Internal error. facmat_ <> NULL, but "
		   << "amFactorised_ = false!";
	  Error( __FILE__, __LINE__ );
	}
      }

#endif

      // Determine bandwidth of the matrix

      // Since the entries in a row are stored in ascending order
      // we only have to compute the distance between the column index
      // of the final entries in each row and the row index to find
      // the bandwidth. As always we assume that there is a diagonal entry,
      // so no row is empty.
      UInt bwGlobal = 0;
      UInt bwLocal  = 0;
      for ( i = 1; i <= lp_n; i++ ) {

	// Determine bandwidth of this row
	bwLocal = (UInt)(cidx[ rptr[i+1] - 1 ] - i);

	// Compare to global bandwidth
	bwGlobal = bwLocal > bwGlobal ? bwLocal : bwGlobal;

      }

      // Save bandwidth of the system matrix
      sysMatBW_ = bwGlobal;

      // Determine number of entries in band matrix
      UInt facMatEntries_ = ( bwGlobal + 1 ) * lp_n;

      // Report to logfile
      if ( logging == true ) {
	(*cla) << "    Bandwidth = " << bwGlobal << "\n"
	       << "    Sparse matrix has " << scrsMat.GetNnz()
	       << " non-zero entries\n"
	       << "    Banded matrix has " << facMatEntries_
	       << " complex entries"
	       << std::endl;
      }

      // Memory allocation for band matrix: We only perform a
      // re-allocation, if the new size requirements are larger
      // than the old ones (real and imag part count separately)
      if ( amFactorised_ == false ) {
	NewArray( facMat_, Double, 2 * facMatEntries_ );
      }
      else {
	if ( 2 * facMatEntries_ > facMatCapacity_ ) {
	  DeleteArray( facMat_ );
	  NewArray( facMat_, Double, 2 * facMatEntries_ );
	  facMatCapacity_ = 2 * facMatEntries_;
	}

#ifdef DEBUG_LAPACK_LL
	else {
	  for ( i = 2 * facMatEntries_ + 1; i <= facMatCapacity_; i++ ) {
	    facMat_[i] = 0.0;
	  }
	}
#endif

      }

    }

    // Set all entries in band matrix to zero
    for ( i = 1; i <= 2 * facMatEntries_; i++ ) {
      facMat_[i] = 0.0;
    }

    // Now we create a symmetric LAPACK band matrix. The SCRS_Matrix stores the
    // upper triangle of the matrix. If we directly insert the rows into the
    // facMat_ array, this gives us in FORTRAN77 format the lower triangle of
    // the matrix, which is fine for LAPACK and is hopefully more cache
    // efficient than converting the rows into FORTRAN format.
    UInt colInit = 0;
    for ( i = 1; i <= lp_n; i++ ) {
      for ( k = rptr[i]; k < rptr[i+1]; k++ ) {
	facMat_[ colInit + 2 * (cidx[k]-i+1) - 1 ] = data[k].real();
	facMat_[ colInit + 2 * (cidx[k]-i+1)     ] = data[k].imag();
      }
      colInit += (facMatEntries_ + 1);
    }

    // Prepare parameters for LAPACK
    lp_kd = facMatEntries_;
    lp_ldab = facMatEntries_ + 1;

    // Now call the LAPACK routine to compute Cholesky factorisation
    LP_ZPBTRF( &lp_uplo, &lp_n, &lp_kd, facMat_ + 1, &lp_ldab, &lp_info );

    // Process return status of LAPACK (i.e. if it is compiled to return
    // on error)
    if ( lp_info < 0 ) {
      (*error) << "Lapack_LL: ZPBTRF reports input value no. " << -lp_info
	       << " is invalid";
      Error( __FILE__, __LINE__ );
    }
    else if ( lp_info > 0 ) {
      (*error) << "Lapack_LL: ZPBTRF reports that leading minor of order "
	       << lp_info << " is not positive definite and the factorization"
	       << " could not be completed";
      Error( __FILE__, __LINE__ );
    }
  }


  // ***********************
  //   Solve linear system
  // ***********************
  void Lapack_LL::Solve( const BaseMatrix &sysMat, const BasePrecond &precond,
			 const BaseVector &rhs, BaseVector &sol ) {

    ENTER_FCN( "Lapack_LL::Solve" );

    // Are we expected to be verbose?
    bool logging = myParams_->GetBoolValue( "LAPACKLL_logging" );

    // Report to logfile
    if ( logging == true ) {
      (*cla) << " --------------------------------------\n"
	     << " LAPACK_LL: Computing solution of linear system"
	     << " with " << sysMat.GetNrows() << " unknowns"
	     << std::endl;
    }

    // Complain, if no factorisation is available
    // or the two indicators disagree
    if ( (facmat_+1) == NULL || amFactorised_ == false ) {
      if ( (facmat_+1) == NULL && amFactorised_ == false ) {
	(*error) << "Lapack_LL: No factorisation is available. Call Setup() "
		 << " first";
	Error( __FILE__, __LINE__ );
      }
      else if ( amFactorised_ == false ) {
	(*error) << "Lapack_LL: Internal error. facmat_ <> NULL, but "
		 << "amFactorised_ = false!";
	Error( __FILE__, __LINE__ );
      }
      else {
	(*error) << "Lapack_LL: Internal error. facmat_ = NULL, but "
		 << "amFactorised_ = true!";
	Error( __FILE__, __LINE__ );
      }
    }

    TRY_CAST {

      // Down-cast matrix to standard matrix
      ConstRefCast( sysMat, StdMatrix, stdMat );

      // Check that we have the correct matrix type
      MatrixStorageType mtype = stdMat.GetStorageType();
      if ( mtype != SPARSE_SYM ) {
	(*error) << "Lapack_LL: Expected a sparseSym matrix, but got a "
		 << Enum2String( mtype ) << " matrix";
	Error( __FILE__, __LINE__ );
      }

      // Get the entry type to figure out which Factorization method to call
      MatrixEntryType etype = stdMat.GetEntryType();

      // Call appropriate solution routine
      switch( etype ) {

      case DOUBLE:
	SolveReal( rhs, sol );
	break;

      case COMPLEX:
	SolveComplex( rhs, sol );
	break;

      default:
	(*error) << "Matrix entry type is neither real nor complex!";
	Error( __FILE__, __LINE__ );
      }

    } CATCH_CAST;

    // Report to logfile
    if ( logging == true ) {
      (*cla) << " --------------------------------------" << std::endl;
    }
  }


  // *************
  //   SolveReal
  // *************
  void Lapack_LL::SolveReal( const BaseVector &rhs, BaseVector &sol ) {

    ENTER_FCN( "Lapack_LL::SolveReal" );

    UInt i;

    // Initialise parameters for ZPBTRS
    char lp_uplo = 'L';
    int  lp_n    = rhs.GetSize();
    int  lp_nrhs = 1;
    int  lp_kd   = (int)sysMatBW_;
    int  lp_ldab = lp_kd + 1;
    int  lp_ldb  = lp_n;
    int  lp_info = 0;

    TRY_CAST {

      // Down-cast the vectors
      ConstRefCast( rhs, Vector<Double>, myRHS );
      RefCast( sol, Vector<Double>, mySol );

      // Get data pointers
      const Double *dataRHS = myRHS.GetPointer();
      Double *dataSol = mySol.GetPointer();

      // See, whether we need to allocate storage for the
      // right-hand side array
      if ( lapackRHSCapacity_ < lp_n ) {
	DeleteArray( lapackRHS_ );
	NewArray( lapackRHS_, Double, lp_n );
	lapackRHSCapacity_ = lp_n;
      }

      // Generate right-hand side for LAPACK
      for ( i = 1; i <= lp_n; i++ ) {
	lapackRHS_[i] = dataRHS[i];
      }

      // Call LAPACK's solution routine
      LP_DPBTRS( &lp_uplo, &lp_n, &lp_kd, &lp_nrhs, facMat_ + 1, &lp_ldab,
		 lapackRHS_ + 1, &lp_ldb, &lp_info );

      // Check return status
      if ( lp_info < 0 ) {
	(*error) << "Lapack_LL: DPBTRS reports input argument number "
		 << -lp_info << " had illegal value";
	Error( __FILE__, __LINE__ );
      }

      // Extract solution from LAPACK array
      for ( i = 1; i <= lp_n; i++ ) {
	dataSol[i] = lapackRHS_[i];
      }

    } CATCH_CAST;

  }


  // ****************
  //   SolveComplex
  // ****************
  void Lapack_LL::SolveComplex( const BaseVector &rhs, BaseVector &sol ) {

    ENTER_FCN( "Lapack_LL::SolveComplex" );

    UInt i;

    // Initialise parameters for ZPBTRS
    char lp_uplo = 'L';
    int  lp_n    = rhs.GetSize();
    int  lp_nrhs = 1;
    int  lp_kd   = (int)sysMatBW_;
    int  lp_ldab = lp_kd + 1;
    int  lp_ldb  = lp_n;
    int  lp_info = 0;

    TRY_CAST {

      // Down-cast the vectors
      ConstRefCast( rhs, Vector<Complex>, myRHS );
      RefCast( sol, Vector<Complex>, mySol );

      // Get data pointers
      const Complex *dataRHS = myRHS.GetPointer();
      Complex *dataSol = mySol.GetPointer();

      // See, whether we need to allocate storage for the
      // right-hand side array
      if ( lapackRHSCapacity_ < 2 * lp_n ) {
	DeleteArray( lapackRHS_ );
	NewArray( lapackRHS_, Double, 2 * lp_n );
	lapackRHSCapacity_ = 2 * lp_n;
      }

      // Generate right-hand side for LAPACK
      for ( i = 1; i <= lp_n; i++ ) {
	lapackRHS_[2*i-1] = dataRHS[i].real();
	lapackRHS_[2*i  ] = dataRHS[i].imag();
      }

      // Call LAPACK's solution routine
      LP_ZPBTRS( &lp_uplo, &lp_n, &lp_kd, &lp_nrhs, facMat_ + 1, &lp_ldab,
		 lapackRHS_ + 1, &lp_ldb, &lp_info );

      // Check return status
      if ( lp_info < 0 ) {
	(*error) << "Lapack_LL: ZPBTRS reports input argument number "
		 << -lp_info << " had illegal value";
	Error( __FILE__, __LINE__ );
      }

      // Extract solution from LAPACK array
      for ( i = 1; i <= lp_n; i++ ) {
	Complex aux( lapackRHS_[2*i-1], lapackRHS_[2*i]  );
	dataSol[i] = aux;
      }

    } CATCH_CAST;

  }

}
