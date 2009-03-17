// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

//#include <omp.h>

#include "MatVec/basematrix.hh"
#include "MatVec/stdmatrix.hh"
#include "MatVec/crs_matrix.hh"
#include "MatVec/scrs_matrix.hh"
#include "OLAS/algsys/olasparams.hh"

#include "pardisosolver.hh"

namespace CoupledField {

#define F77_FUNC(func)   func ## _

  extern "C" {

    int F77_FUNC(pardisoinit) (int *, int *, int *);

    int F77_FUNC(pardiso) (int *, int *, int *, int *, int *, int *,
                           const Double *, const int *, const int *, int *,
                           int *, int *, int *, const Double *,
                           Double *, int *);
  }


  // ***********************
  //   Default Constructor
  // ***********************
  template<typename T>
  PardisoSolver<T>::PardisoSolver() {
    EXCEPTION( "Default constructor of PardisoSolver is forbidden!" );
  }


  template<typename T>
  std::string PardisoSolver<T>::GetErrorString(int err_code) {
    switch (err_code) {
      case NO_ERROR:
        return "No error.";
      case INPUT_INCONSISTENT:
        return "Input inconsistent.";
      case NOT_ENOUGH_MEMORY:
        return "Not enough memory.";
      case REORDERING_PROBLEM:
        return "Reordering problem.";
      case ZERO_PIVOT:
        return "Zero pivot, numerical factorization or iterative refinement problem.";
      case PREORDERING_FAILED:
        return "Preordering failed.";
      case DIAGONAL_MATRIX:
        return "Diagonal matrix problem.";
      case INT_OVERFLOW:
        return "32-bit integer overflow problem.";
      default:
        return "Unclassified (internal) error.";
    }
  }

  // ***************
  //   Constructor
  // ***************
  template<typename T>
  PardisoSolver<T>::PardisoSolver( OLAS_Params *myParams,
                                   OLAS_Report *myReport ) {


    // Set pointers to communication objects
    myParams_  = myParams;
    myReport_  = myReport;

    // Initialise attributes
    firstCall_ = true;
    msgLvl_ = 0;

    // Our private Fortran zeros
    zeroINT_ = 0;
    zeroDBL_ = 0.0;

    // For production runs, we need no identity reordering
    idPerm_     = NULL;
    idPermSize_ = 0;

    // These pointers are used to hold the addresses of the internal
    // (S)CRS matrix structures. The related memory segments must not
    // be altered of deleted. Therefore the pointers are const!
    rowPtr_ = NULL;
    colPtr_ = NULL;
    datPtr_ = NULL;
  }


  // **************
  //   Destructor
  // **************
  template<typename T>
  PardisoSolver<T>::~PardisoSolver() {

    // PARDISO - Last Phase: Cleaning up the parameters
    if ( firstCall_ == false ) {

      int errorFlag = 0;
      int phase = -1;

      F77_FUNC(pardiso) ( pt_, &maxfct, &mnum, &mType_, &phase,
                          &probDim_, theMatrix_, rowPtr_, colPtr_,
                          idPerm_, &nrhs, iparm_+1, &msgLvl_, &zeroDBL_,
                          &zeroDBL_, &errorFlag );

      if ( errorFlag != NO_ERROR) {
        EXCEPTION( "Pardiso: Error occured during cleanup: "
                   << GetErrorString(errorFlag) )
      }
    }

    // Delete identity re-ordering (if exists)
    DELETEARRAY( idPerm_ );
    idPermSize_ = 0;

  }


  // *********
  //   Setup
  // *********
  template<typename T>
  void PardisoSolver<T>::Setup( BaseMatrix &sysMat ) {

    // Flag for check Pardiso's return status
    int errorFlag = 0;

    // Determine, whether we are expected to be verbose
    bool logging = myParams_->GetBoolValue( "PARDISO_logging" );
    if ( logging == true ) {
      (*cla) << " -------------------------------------------------------"
             << "-----------------------\n";
    }

    // ============================================================
    //  Determine which of the two steps: symbolical and numerical
    //  factorisation must be performed
    // ============================================================
    bool facSymbolic = false;
    bool facNumeric = false;

    // No factorisation available, so perform both steps
    if ( firstCall_ == true ) {
      facSymbolic = true;
      facNumeric  = true;
    }

    else {

      // When the matrix pattern has changed, we need to re-do
      // both steps, also the symbolical one
      if( myParams_->GetBoolValue( "newMatrixPattern" ) ) {
        facSymbolic = true;
        facNumeric  = true;
      }

      // If only the values of the matrix entries changed, we
      // can keep the re-ordering and only perform a numerical
      // factorisation
      else {
        facSymbolic = false;
        facNumeric  = true;
      }
    }

    // =====================================
    //  Check matrix entry and storage type
    // =====================================
    TRY_CAST {

      CONSTREFCAST( sysMat, StdMatrix, stdMat );

      etype = stdMat.GetEntryType();
      if ( (etype != BaseMatrix::DOUBLE) && (etype != BaseMatrix::COMPLEX) ) {
        EXCEPTION( "Pardiso: Expected DOUBLE or COMPLEX entries, but got '"
                 << BaseMatrix::entryType.ToString(etype) << "'" );
      }

      stype = stdMat.GetStorageType();
      if ( (stype != BaseMatrix::SPARSE_SYM) &&
          (stype != BaseMatrix::SPARSE_NONSYM) ) {
        EXCEPTION( "Pardiso: Expected a sparseSym or sparseNonSym matrix, "
                 << "but got a '" << BaseMatrix::storageType.ToString( stype ) << "' matrix" );
      }

      // Determine problem size
      probDim_ = stdMat.GetNumRows();

    } CATCH_CAST;


    // ==================================================
    //  Get pointers to arrays containing information on
    //  (S)CRS matrix from the problem matrix object
    // ==================================================
    if ( stype == BaseMatrix::SPARSE_NONSYM ) {
      CONSTREFCAST( sysMat, CRS_Matrix<T>, crsMat );
      rowPtr_ = (Integer*)(crsMat.GetRowPointer());
      colPtr_ = (Integer*)(crsMat.GetColPointer());
      datPtr_ = crsMat.GetDataPointer();
      nnz_    = crsMat.GetNnz();
    }
    else {
      CONSTREFCAST( sysMat, SCRS_Matrix<T>, scrsMat );
      rowPtr_ = (Integer*)(scrsMat.GetRowPointer());
      colPtr_ = (Integer*)(scrsMat.GetColPointer());
      datPtr_ = scrsMat.GetDataPointer();
      nnz_    = scrsMat.GetNumEntries();
    }

    // Increment pointers to make them 0-based
//     rowPtr_++;
//     colPtr_++;
//     datPtr_++;

    // Do C-style casting to fix problem with over-loading vs. extern C
    void *hatred = (void *) datPtr_;
    theMatrix_ = (Double *) hatred;


    // ====================================
    //  Set default parameters for Pardiso
    // ====================================

    // Some flags for determining the type of the matrix
    bool symPard, defPard, herPard, strPard;

    if ( stype == BaseMatrix::SPARSE_SYM ) {
      symPard = true;
    }
    else {
      symPard = false;
    }

    mType_ = 0;

    defPard = myParams_->GetBoolValue( "PARDISO_posDef" );
    herPard = myParams_->GetBoolValue( "PARDISO_hermitian" );
    strPard = myParams_->GetBoolValue( "PARDISO_symStructure" );

    if ( (etype == BaseMatrix::DOUBLE ) && (!symPard) && ( strPard) ) mType_ =  1;
    if ( (etype == BaseMatrix::DOUBLE ) && ( symPard) && ( defPard) ) mType_ =  2;
    if ( (etype == BaseMatrix::DOUBLE ) && ( symPard) && (!defPard) ) mType_ = -2;
    if ( (etype == BaseMatrix::COMPLEX) && (!symPard) && ( strPard) ) mType_ =  3;
    if ( (etype == BaseMatrix::COMPLEX) && ( herPard) && ( defPard) ) mType_ =  4;
    if ( (etype == BaseMatrix::COMPLEX) && ( herPard) && (!defPard) ) mType_ = -4;
    if ( (etype == BaseMatrix::COMPLEX) && ( symPard) ) mType_ = 6;
    if ( (etype == BaseMatrix::DOUBLE ) && (!symPard) && (!strPard) ) mType_ = 11;
    if ( (etype == BaseMatrix::COMPLEX) && (!symPard) && (!strPard) && (!herPard)) {
      mType_ = 13;
    }
    if ( mType_ == 0 ) {
      EXCEPTION( "PardisoSolver: There appears to be an inconsistency in "
               << "the input parameters. I cannot determine correct matrix "
               << "properties for pardiso" );
    }
    else if ( logging == true ) {
      (*cla) << " Pardiso: Classified matrix as mType = "
             << mType_ << std::endl;
    }

    // We do not need to call pardisoinit, if the matrix pattern
    // has not changed
    if ( facSymbolic ) {
      if ( logging == true ) {
        (*cla) << " Pardiso: Calling pardisoinit" << std::endl;
      }
      F77_FUNC(pardisoinit) ( pt_, &mType_, iparm_+1 );
    }


    // =======================================
    //  Alter default parameters to our taste
    // =======================================

    // Avoid that Pardiso over-writes our settings
    iparm_[1] = 1;

    // Specify number of OpenMP threads. OLAS currently has no consistent
    // concept for shared memory parallelism, so we specify one thread here.
    //#if defined (_OPENMP)
    //iparm_[3] = omp_get_num_threads();
    //#else
    iparm_[3] = 1;
    //#endif

    // Determine the re-ordering strategy: We can either fo nested dissection
    // or minimum degree re-ordering or no re-ordering at all (i.e. we use
    // the initial ordering of the linear system, which might already have been
    // re-ordered via the graph)
    ReorderingType ordering;
    myParams_->GetEnumValue( "PARDISO_ordering", ordering );
    switch ( ordering ) {

    case NESTED_DISSECTION:
      iparm_[2] = 2;
      iparm_[5] = 0;
      break;

    case MINIMUM_DEGREE:
      iparm_[2] = 0;
      iparm_[5] = 0;
      break;

    case NOREORDERING:
      // In this case iparm_[2] is irrelevant, we generate an identity
      // permutation and use this one, by setting iparm_[5] and skip the
      // symbolic factorisation
      iparm_[2] = 0;
      iparm_[5] = 1;
      if ( idPermSize_ < probDim_ ) {
        DELETEARRAY( idPerm_ );
        NEWARRAY( idPerm_, int, probDim_ );
        for ( int i = 0; i < probDim_; i++ ) {
          // (i+1), since fortran needs indices starting with 1!!
          idPerm_[i] = i+1;
        }
      }
      break;

    default:
      std::string tmp;
      Enum2String( ordering, tmp );

      EXCEPTION( "Re-ordering of type '" << tmp
               << "' is not available with the PardisoSolver" );
    }

    if ( logging == true && facSymbolic == true ) {
      if ( ordering != NOREORDERING ) {
        std::string tmp;
        Enum2String( ordering, tmp );
        (*cla) << " Pardiso: Analyse phase will determine a '"
               << tmp << "' re-ordering" << std::endl;
      }
      else {
        (*cla) << " Pardiso: Factorisation uses original matrix ordering"
               << std::endl;
      }
    }

    // Do we need to determine MFLOPs for the LU factorisation
    iparm_[19] = 0;
    if ( myParams_->GetBoolValue( "PARDISO_stats" ) == true ) {
      iparm_[19] = -1;
    }

    // Setting pivoting strategy for indefinit problems
    iparm_[21] = myParams_->GetIntValue( "PARDISO_pivoting" );

    // In case we have no positive definite system (especially piezo)
    // we perform additional scaling to enhance the condition for very
    // small off-diagonal entries (iparm_[11]). In addition we enable
    // the method of 'symmetric weighted matchings' (iparam_[13]).
    // For further information, refer to the pardiso user manual.
    if( !defPard ) {
      iparm_[11] = 1;
      iparm_[13] = 1;
    }

    // Pardiso keeps one factorisation in memory (and that is used for
    // the solution phase)
    maxfct = 1;
    mnum = 1;

    // We simultaneously treat a single right hand side
    nrhs = 1;

    // Set the message level (for printing statistics)
    msgLvl_ = 0;
    if ( myParams_->GetBoolValue( "PARDISO_stats" ) == true ) {
      msgLvl_ = 1;
    }


    // we have to icrement the entries of the col- and row-position arrays
    // by one, so that the first col and first row start with index 1 (and
    // not with zero) to be consistent with fortran
    // at the end of the method we will undo it!!
    for (UInt i=0; i < static_cast<UInt>(probDim_+1); i++ )
      rowPtr_[i] += 1;

    for (UInt i=0; i< nnz_; i++ )
       colPtr_[i] += 1;

    // ========================
    //  Symbolic Factorisation
    // ========================
    if ( facSymbolic == true ) {

      // log report
      if ( logging == true ) {
        (*cla) << " Pardiso: Performing analyse phase (symbolic factorisation)"
               << " ... " << std::flush;
      }

      // only analyse
      int phase = 11;

      // let pardiso go for it
      F77_FUNC(pardiso) (pt_, &maxfct, &mnum, &mType_, &phase,
                         &probDim_, theMatrix_, rowPtr_, colPtr_,
                         idPerm_, &nrhs, iparm_+1, &msgLvl_, &zeroDBL_,
                         &zeroDBL_, &errorFlag );

      // Check return status
      if ( errorFlag != NO_ERROR ) {
        EXCEPTION( "Pardiso: Error occured during symbolic factorization: "
                   << GetErrorString(errorFlag) );
      }
      else {
        if ( logging == true ) {
          (*cla) << "done" << std::endl;
        }
      }
    }


    // =========================
    //  Numerical Factorisation
    // =========================
    if ( facNumeric == true ) {

      // log report
      if ( logging == true ) {
        (*cla) << " Pardiso: Performing factorise phase (numerical "
               << "factorisation) ... " << std::flush;
      }

      // only factorise (numerical)
      int phase = 22;

      // let pardiso go for it
      F77_FUNC(pardiso) (pt_, &maxfct, &mnum, &mType_, &phase,
                         &probDim_, theMatrix_, rowPtr_, colPtr_,
                         idPerm_, &nrhs, iparm_+1, &msgLvl_, &zeroDBL_,
                         &zeroDBL_, &errorFlag );

      // Check return status
      if ( errorFlag != NO_ERROR ) {
        EXCEPTION( "Pardiso: Error occured during numerical factorization: "
                   << GetErrorString(errorFlag) );
      }
      else {
        if ( logging == true ) {
          (*cla) << "done" << std::endl;
        }
      }
    }

    // Now we were called once, and a factorisation is available
    firstCall_ = false;

    // Finish log report
    if ( logging == true ) {
    
      (*cla) << " -------------------------------------------------------"
             << "-----------------------\n";
    }

    // now we undo our increment, since on our side the frist col and row
    // has an value of zero!!
    for (UInt i=0; i <  static_cast<UInt>(probDim_+1); i++ )
      rowPtr_[i] -= 1;

    for (UInt i=0; i< nnz_; i++ )
      colPtr_[i] -= 1;

  }



  // *************************
  //   Solve linear system
  // *************************
  template<typename T>
  void PardisoSolver<T>::Solve( const BaseMatrix &sysmat,
                                const BasePrecond &precond,
                                const BaseVector &rhs, BaseVector &sol ) {


    // Determine, whether we are expected to be verbose
    bool logging = myParams_->GetBoolValue( "PARDISO_logging" );
    if ( logging == true ) {
      (*cla) << " -------------------------------------------------------"
             << "-----------------------\n"
             << " Pardiso: Solving linear system ... ";
    }

    if ( firstCall_ == true ) {
      EXCEPTION( "The matrix has not yet been factorised by Pardiso! "
               << "Call Setup() first" );
    }


    // Check that we have the correct vector types and
    // obtain data pointers
    const T *rhsArray;
    T* solArray;
    TRY_CAST {
      CONSTREFCAST( rhs, Vector<T>, myRHS );
      REFCAST( sol, Vector<T>, mySol );
      rhsArray = myRHS.GetPointer();
      solArray = mySol.GetPointer();
    } CATCH_CAST;

    // We must perform a nasty cast in order to be able to interface with
    // a pardiso routine of the same name in both cases, real and complex
    // vector entries.
    void *hatred = NULL;
    Double *theRHS = NULL;
    Double *theSol = NULL;

    hatred = (void *) rhsArray;
    theRHS = (Double *) hatred;

    hatred = (void *) solArray;
    theSol = (Double *) hatred;

    // PARDISO - Third Phase : Solution of the system
    int errorFlag = 0;
    int phase = 33;

    // we have to icrement the entries of the col- and row-position arrays
    // by one, so that the first col and first row start with index 1 (and
    // not with zero) to be consistent with fortran
    // at the end of the method we will undo it!!
    for (UInt i=0; i<static_cast<UInt>(probDim_+1); i++ )
      rowPtr_[i] += 1;
    for (UInt i=0; i< nnz_; i++ )
       colPtr_[i] += 1;

    F77_FUNC(pardiso) (pt_, &maxfct, &mnum, &mType_, &phase,
                       &probDim_, theMatrix_, rowPtr_, colPtr_,
                       idPerm_, &nrhs, iparm_+1, &msgLvl_, theRHS,
                       theSol, &errorFlag );

    // now we undo our increment, since on our side the first col and row
    // has an value of zero!!
    for (UInt i=0; i <  static_cast<UInt>(probDim_+1); i++ )
      rowPtr_[i] -= 1;
    for (UInt i=0; i< nnz_; i++ )
      colPtr_[i] -= 1;

    // Check return status
    if ( errorFlag != NO_ERROR ) {
      EXCEPTION( "Pardiso: Error occured during solution of linear system: "
                 << GetErrorString(errorFlag) );
    }
    else {
      if ( logging == true ) {
        (*cla) << "done" << std::endl;
      }
    }

    // Finish log report
    if ( logging == true ) {
      (*cla) << " number of iterative refinement steps: " << iparm_[7] << std::endl;
      (*cla) << " number of perturbed pivots: " << iparm_[14] << std::endl;
      (*cla) << " number of positive eigenvalues: " << iparm_[22] << std::endl;
      (*cla) << " number of negative eigenvalues: " << iparm_[23] << std::endl;
      (*cla) << " -------------------------------------------------------"
             << "-----------------------\n";
    }

    // Create Report (no sensible things to write for direct solvers yet)
    if ( myReport_ != NULL ) {
      myReport_->SetValue( "numIter", -1 );
      myReport_->SetValue( "finalNorm", -1.0 );
    }
  }

// Explicit template instantiation
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template class PardisoSolver<Double>;
  template class PardisoSolver<Complex>;
#endif

}
