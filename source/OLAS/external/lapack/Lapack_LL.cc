#include "MatVec/BaseMatrix.hh"
#include "MatVec/StdMatrix.hh"
#include "MatVec/SCRS_Matrix.hh"
#include "MatVec/BLASLAPACKInterface.hh"
#include "OLAS/external/lapack/Lapack_LL.hh"

namespace CoupledField {

  // ***********************
  //   Default Constructor
  // ***********************
  Lapack_LL::Lapack_LL() {
    EXCEPTION("Default constructor of LAPACK_LL is forbidden!");
  }


  // *****************************
  //   Constructor (for factory)
  // *****************************
  Lapack_LL::Lapack_LL( PtrParamNode solverNode, PtrParamNode olasInfo ) {
    
    xml_ = solverNode;
    infoNode_ = olasInfo->Get("lapackLL");
    
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
    delete [] ( lapackRHS_ );
    delete [] ( facMat_ );
  }


  // ********************************
  //   Setup: Perform Factorisation
  // ********************************
  void Lapack_LL::Setup( BaseMatrix &sysMat) {

    StdMatrix& stdMat = dynamic_cast<StdMatrix&>(sysMat);

    // Check that we have the correct matrix type
    BaseMatrix::StorageType mtype = stdMat.GetStorageType();
    if ( mtype != BaseMatrix::SPARSE_SYM ) {
      EXCEPTION("Lapack_LL: Expected a sparseSym matrix, but got a "
          << MatrixStorageTypeEnum.ToString( mtype ) << " matrix");
    }

    // Get the entry type to figure out which Factorisation method to call
    BaseMatrix::EntryType etype = stdMat.GetEntryType();

    // Call appropriate factorisation routine
    switch( etype ) {

    case BaseMatrix::DOUBLE:
      FactoriseReal( stdMat );
      break;

    case BaseMatrix::COMPLEX:
      FactoriseComplex( stdMat );
      break;

    default:
      EXCEPTION("Matrix entry type is neither real nor complex!");
    }

    // now we have a (new) factorisation
    amFactorised_ = true;
  }


  // *****************
  //   FactoriseReal
  // *****************
  void Lapack_LL::FactoriseReal( StdMatrix &stdMat ) {


    UInt i, k;

    // Initialise parameters for DPBTRF
    char lp_uplo = 'L';
    int  lp_n    = stdMat.GetNumCols();
    int  lp_kd   = 0;
    int  lp_ldab = 0;
    int  lp_info = 0;

    // Down-cast matrix to special format
    SCRS_Matrix<Double>& scrsMat = dynamic_cast<SCRS_Matrix<Double>&>(stdMat);

    // Get hold of column index array
    const UInt *cidx = scrsMat.GetColPointer();

    // Get hold of row pointer index array
    const UInt *rptr = scrsMat.GetRowPointer();

    // Get hold of data array
    const Double *data = scrsMat.GetDataPointer();

    // Compute number of stored matrix entries
    // COMPWARNING: unused variable Integer nEntries = ( scrsMat.GetNnz() + lp_n ) / 2;

    // TODO: THIS CHECK DOES NOT MAKE SENSE IN MY OPINION SINCE
    //       'newMatrixPattern' is set to false in olasparams.cc
    //       and gets never changed elsewhere. A more intelligent
    //       test would ask the matrix if its pattern did change.

    bool newMatrixPattern = false;
      
    // Test, if we must check for a new matrix' bandwidth and
    // maybe also allocate new memory
    if ( amFactorised_ == false || facMat_ == NULL ||
        newMatrixPattern ) {

#ifdef DEBUG_LAPACK_LL

      // Are we over-writing an existing factorisation? Check, whether the
      // two indicators agree, if not cry out loud!
      if ( (facmat_+1) != NULL || amFactorised_ == true ) {
        if ( (facmat_+1) == NULL && amFactorised_ == true ) {
          EXCEPTION( "Lapack_LL: Internal error. facmat_ = NULL, but "
                     << "amFactorised_ = true!" );
        }
        else if ( amFactorised_ == false ) {
          EXCEPTION( "Lapack_LL: Internal error. facmat_ <> NULL, but "
                     << "amFactorised_ = false!" );
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
      for ( i = 1; i <= static_cast<UInt>(lp_n); i++ ) {

        // Determine bandwidth of this row
        bwLocal = (UInt)(cidx[ rptr[i+1] - 1 ] - i);

        // Compare to global bandwidth
        bwGlobal = bwLocal > bwGlobal ? bwLocal : bwGlobal;

      }

      // Save bandwidth of the system matrix
      sysMatBW_ = bwGlobal;

      // Determine number of entries in band matrix
      facMatEntries_ = ( bwGlobal + 1 ) * lp_n;

      // Memory allocation for band matrix: We only perform a
      // re-allocation, if the new size requirements are larger
      // than the old ones
      if ( amFactorised_ == false ) {
        NEWARRAY( facMat_, Double, facMatEntries_ );
      }
      else {
        if ( facMatEntries_ > facMatCapacity_ ) {
          delete [] ( facMat_ );  facMat_  = NULL;
          NEWARRAY( facMat_, Double, facMatEntries_ );
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
    for ( i = 1; i <= static_cast<UInt>(lp_n); i++ ) {
      for ( k = rptr[i]; k < rptr[i+1]; k++ ) {
        facMat_[ colInit + cidx[k] - i + 1 ] = data[k];
      }
      colInit += (sysMatBW_ + 1);
    }

    // Prepare parameters for LAPACK
    lp_kd = sysMatBW_;
    lp_ldab = sysMatBW_ + 1;

    // Now call the LAPACK routine to compute Cholesky factorisation
    dpbtrf( &lp_uplo, &lp_n, &lp_kd, facMat_ + 1, &lp_ldab, &lp_info );

    // Process return status of LAPACK (i.e. if it is compiled to return
    // on error)
    if ( lp_info < 0 ) {
      EXCEPTION("Lapack_LL: DPBTRF reports input value no. " << -lp_info
                << " is invalid");
    }
    else if ( lp_info > 0 ) {
      EXCEPTION("Lapack_LL: DPBTRF reports that leading minor of order "
                << lp_info << " is not positive definite and the factorization"
                << " could not be completed");
    }
  }


  // ********************
  //   FactoriseComplex
  // ********************
  void Lapack_LL::FactoriseComplex( StdMatrix &stdMat ) {


    UInt i, k;

    // Initialise parameters for ZPBTRF
    char lp_uplo = 'L';
    int  lp_n    = stdMat.GetNumCols();
    int  lp_kd   = 0;
    int  lp_ldab = 0;
    int  lp_info = 0;

    // Down-cast matrix to special format
    SCRS_Matrix<Complex>& scrsMat = dynamic_cast<SCRS_Matrix<Complex>&>(stdMat);

    // Get hold of column index array
    const UInt *cidx = scrsMat.GetColPointer();

    // Get hold of row pointer index array
    const UInt *rptr = scrsMat.GetRowPointer();

    // Get hold of data array
    const Complex *data = scrsMat.GetDataPointer();

    // Compute number of stored matrix entries
    // COMPWARNING: unused variable Integer nEntries = ( scrsMat.GetNnz() + lp_n ) / 2;
    
    // TODO: THIS CHECK DOES NOT MAKE SENSE IN MY OPINION SINCE
    //       'newMatrixPattern' is set to false in olasparams.cc
    //       and gets never changed elsewhere. A more intelligent
    //       test would ask the matrix if its pattern did change.

    bool newMatrixPattern = false;
      
    // Test, if we must check for a new matrix' bandwidth and
    // maybe also allocate new memory
    if ( amFactorised_ == false || facMat_ == NULL ||
        newMatrixPattern ) {

#ifdef DEBUG_LAPACK_LL

      // Are we over-writing an existing factorisation? Check, whether the
      // two indicators agree, if not cry out loud!
      if ( (facmat_+1) != NULL || amFactorised_ == true ) {
        if ( (facmat_+1) == NULL && amFactorised_ == true ) {
          EXCEPTION( "Lapack_LL: Internal error. facmat_ = NULL, but "
                     << "amFactorised_ = true!" );
        }
        else if ( amFactorised_ == false ) {
          EXCEPTION("Lapack_LL: Internal error. facmat_ <> NULL, but "
                    << "amFactorised_ = false!");
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
      for ( i = 1; (int) i <= lp_n; i++ ) {

        // Determine bandwidth of this row
        bwLocal = (UInt)(cidx[ rptr[i+1] - 1 ] - i);

        // Compare to global bandwidth
        bwGlobal = bwLocal > bwGlobal ? bwLocal : bwGlobal;

      }

      // Save bandwidth of the system matrix
      sysMatBW_ = bwGlobal;

      // Determine number of entries in band matrix
      UInt facMatEntries_ = ( bwGlobal + 1 ) * lp_n;

      // Memory allocation for band matrix: We only perform a
      // re-allocation, if the new size requirements are larger
      // than the old ones (real and imag part count separately)
      if ( amFactorised_ == false ) {
        NEWARRAY( facMat_, Double, 2 * facMatEntries_ );
      }
      else {
        if ( 2 * facMatEntries_ > facMatCapacity_ ) {
          delete [] ( facMat_ );  facMat_  = NULL;
          NEWARRAY( facMat_, Double, 2 * facMatEntries_ );
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
    for ( i = 1; (int) i <= lp_n; i++ ) {
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
    zpbtrf( &lp_uplo, &lp_n, &lp_kd, facMat_ + 1, &lp_ldab, &lp_info );

    // Process return status of LAPACK (i.e. if it is compiled to return
    // on error)
    if ( lp_info < 0 ) {
      EXCEPTION("Lapack_LL: ZPBTRF reports input value no. " << -lp_info
                << " is invalid");
    }
    else if ( lp_info > 0 ) {
      EXCEPTION("Lapack_LL: ZPBTRF reports that leading minor of order "
                << lp_info << " is not positive definite and the factorization"
                << " could not be completed");
    }
  }


  // ***********************
  //   Solve linear system
  // ***********************
  void Lapack_LL::Solve( const BaseMatrix &sysMat, const BaseVector &rhs, BaseVector &sol) {

    // Complain, if no factorisation is available
    // or the two indicators disagree
    if ( (facmat_+1) == NULL || amFactorised_ == false ) {
      if ( (facmat_+1) == NULL && amFactorised_ == false ) {
        EXCEPTION("Lapack_LL: No factorisation is available. Call Setup() "
            << " first");
      }
      else if ( amFactorised_ == false ) {
        EXCEPTION("Lapack_LL: Internal error. facmat_ <> NULL, but "
            << "amFactorised_ = false!");
      }
      else {
        EXCEPTION("Lapack_LL: Internal error. facmat_ = NULL, but "
            << "amFactorised_ = true!");
      }
    }

    const StdMatrix& stdMat = dynamic_cast<const StdMatrix&>(sysMat);

    // Check that we have the correct matrix type
    BaseMatrix::StorageType mtype = stdMat.GetStorageType();
    if ( mtype != BaseMatrix::SPARSE_SYM ) {
      EXCEPTION("Lapack_LL: Expected a sparseSym matrix, but got a "
          << MatrixStorageTypeEnum.ToString( mtype ) << " matrix");
    }

    // Get the entry type to figure out which Factorization method to call
    BaseMatrix::EntryType etype = stdMat.GetEntryType();

    // Call appropriate solution routine
    switch( etype ) {

    case BaseMatrix::DOUBLE:
      SolveReal( rhs, sol );
      break;

    case BaseMatrix::COMPLEX:
      SolveComplex( rhs, sol );
      break;

    default:
      EXCEPTION("Matrix entry type is neither real nor complex!");
    }
  }


  // *************
  //   SolveReal
  // *************
  void Lapack_LL::SolveReal( const BaseVector &rhs, BaseVector &sol ) {


    int i;

    // Initialise parameters for ZPBTRS
    char lp_uplo = 'L';
    int  lp_n    = rhs.GetSize();
    int  lp_nrhs = 1;
    int  lp_kd   = (int)sysMatBW_;
    int  lp_ldab = lp_kd + 1;
    int  lp_ldb  = lp_n;
    int  lp_info = 0;

    const Vector<Double>& myRHS = dynamic_cast<const Vector<Double>&>(rhs);
    Vector<Double>& mySol = dynamic_cast<Vector<Double>&>(sol);

      // Get data pointers
      const Double *dataRHS = myRHS.GetPointer();
      Double *dataSol = mySol.GetPointer();

      // See, whether we need to allocate storage for the
      // right-hand side array
      if ( (int) lapackRHSCapacity_ < lp_n ) {
        delete [] ( lapackRHS_ );  lapackRHS_  = NULL;
        NEWARRAY( lapackRHS_, Double, lp_n );
        lapackRHSCapacity_ = lp_n;
      }

      // Generate right-hand side for LAPACK
      for ( i = 1; i <= lp_n; i++ ) {
        lapackRHS_[i] = dataRHS[i];
      }

      // Call LAPACK's solution routine
      dpbtrs( &lp_uplo, &lp_n, &lp_kd, &lp_nrhs, facMat_ + 1, &lp_ldab,
              lapackRHS_ + 1, &lp_ldb, &lp_info );

      // Check return status
      if ( lp_info < 0 ) {
        EXCEPTION("Lapack_LL: DPBTRS reports input argument number "
            << -lp_info << " had illegal value");
      }

      // Extract solution from LAPACK array
      for ( i = 1; i <= lp_n; i++ ) {
        dataSol[i] = lapackRHS_[i];
      }
  }


  // ****************
  //   SolveComplex
  // ****************
  void Lapack_LL::SolveComplex( const BaseVector &rhs, BaseVector &sol ) {


    int i;

    // Initialise parameters for ZPBTRS
    char lp_uplo = 'L';
    int  lp_n    = rhs.GetSize();
    int  lp_nrhs = 1;
    int  lp_kd   = (int)sysMatBW_;
    int  lp_ldab = lp_kd + 1;
    int  lp_ldb  = lp_n;
    int  lp_info = 0;

    const Vector<Complex>& myRHS = dynamic_cast<const Vector<Complex>&>(rhs);
    Vector<Complex>& mySol = dynamic_cast<Vector<Complex>&>(sol);

      // Get data pointers
      const Complex *dataRHS = myRHS.GetPointer();
      Complex *dataSol = mySol.GetPointer();

      // See, whether we need to allocate storage for the
      // right-hand side array
      if ( (int) lapackRHSCapacity_ < 2 * lp_n ) {
        delete [] ( lapackRHS_ );  lapackRHS_  = NULL;
        NEWARRAY( lapackRHS_, Double, 2 * lp_n );
        lapackRHSCapacity_ = 2 * lp_n;
      }

      // Generate right-hand side for LAPACK
      for ( i = 1; i <= lp_n; i++ ) {
        lapackRHS_[2*i-1] = dataRHS[i].real();
        lapackRHS_[2*i  ] = dataRHS[i].imag();
      }

      // Call LAPACK's solution routine
      zpbtrs( &lp_uplo, &lp_n, &lp_kd, &lp_nrhs, facMat_ + 1, &lp_ldab,
              lapackRHS_ + 1, &lp_ldb, &lp_info );

      // Check return status
      if ( lp_info < 0 ) {
        EXCEPTION("Lapack_LL: ZPBTRS reports input argument number "
            << -lp_info << " had illegal value");
      }

      // Extract solution from LAPACK array
      for ( i = 1; i <= lp_n; i++ ) {
        Complex aux( lapackRHS_[2*i-1], lapackRHS_[2*i]  );
        dataSol[i] = aux;
      }
  }
}
