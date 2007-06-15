// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

//#include <omp.h>

#include "utils/utils.hh"
#include "external/pardiso/pardisosolver.hh"

namespace OLAS {

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
    ENTER_FCN( "PardisoSolver::PardisoSolver" );
    Error( "Default constructor of PardisoSolver is forbidden!",
           __FILE__, __LINE__ );
  }


  // ***************
  //   Constructor
  // ***************
  template<typename T>
  PardisoSolver<T>::PardisoSolver( OLAS_Params *myParams,
                                   OLAS_Report *myReport ) {

    ENTER_FCN( "PardisoSolver::PardisoSolver" );

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
    ENTER_FCN( "PardisoSolver::~PardisoSolver" );

    // PARDISO - Last Phase: Cleaning up the parameters
    if ( firstCall_ == false ) {

      int errorFlag = 0;
      int phase = -1;

      F77_FUNC(pardiso) ( pt_, &maxfct, &mnum, &mType_, &phase,
                          &probDim_, theMatrix_, rowPtr_, colPtr_,
                          (idPerm_+1), &nrhs, (iparm_+1), &msgLvl_, &zeroDBL_,
                          &zeroDBL_, &errorFlag );
      
      if ( errorFlag != 0) {
        (*error) << "Pardiso: Error " << errorFlag
                 << " occured during cleanup";
        Error( __FILE__, __LINE__ );
      }
    }

    // Delete identity re-ordering (if exists)
    DeleteArray( idPerm_ );
    idPermSize_ = 0;

  }


  // *********
  //   Setup
  // *********
  template<typename T>
  void PardisoSolver<T>::Setup( BaseMatrix &sysMat ) {

    ENTER_FCN( "PardisoSolver::Setup" );

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

      ConstRefCast( sysMat, StdMatrix, stdMat );

      etype = stdMat.GetEntryType();
      if ( (etype != DOUBLE) && (etype != COMPLEX) ) {
        (*error) << "Pardiso: Expected DOUBLE or COMPLEX entries, but got '"
                 << Enum2String(etype) << "'";
        Error( __FILE__, __LINE__ );
      }

      stype = stdMat.GetStorageType();
      if ( (stype != SPARSE_SYM) && (stype != SPARSE_NONSYM) ) {
        (*error) << "Pardiso: Expected a sparseSym or sparseNonSym matrix, "
                 << "but got a '" << Enum2String(stype) << "' matrix";
        Error( __FILE__, __LINE__ );
      }

      // Determine problem size
      probDim_ = stdMat.GetNrows();

    } CATCH_CAST;


    // ==================================================
    //  Get pointers to arrays containing information on
    //  (S)CRS matrix from the problem matrix object
    // ==================================================
    if ( stype == SPARSE_NONSYM ) {
      ConstRefCast( sysMat, CRS_Matrix<T>, crsMat );
      rowPtr_ = crsMat.GetRowPointer();
      colPtr_ = crsMat.GetColPointer();
      datPtr_ = crsMat.GetDataPointer();
    }
    else {
      ConstRefCast( sysMat, SCRS_Matrix<T>, scrsMat );
      rowPtr_ = scrsMat.GetRowPointer();
      colPtr_ = scrsMat.GetColPointer();
      datPtr_ = scrsMat.GetDataPointer();
    }

    // Increment pointers to make them 0-based
    rowPtr_++;
    colPtr_++;
    datPtr_++;

    // Do C-style casting to fix problem with over-loading vs. extern C
    void *hatred = (void *) datPtr_;
    theMatrix_ = (Double *) hatred;


    // ====================================
    //  Set default parameters for Pardiso
    // ====================================

    // Some flags for determining the type of the matrix
    bool symPard, defPard, herPard, strPard;

    if ( stype == SPARSE_SYM ) {
      symPard = true;
    }
    else {
      symPard = false;
    }

    mType_ = 0;

    defPard = myParams_->GetBoolValue( "PARDISO_posDef" );
    herPard = myParams_->GetBoolValue( "PARDISO_hermitian" );
    strPard = myParams_->GetBoolValue( "PARDISO_symStructure" );

    if ( (etype == DOUBLE ) && (!symPard) && ( strPard) ) mType_ =  1;
    if ( (etype == DOUBLE ) && ( symPard) && ( defPard) ) mType_ =  2;
    if ( (etype == DOUBLE ) && ( symPard) && (!defPard) ) mType_ = -2;
    if ( (etype == COMPLEX) && (!symPard) && ( strPard) ) mType_ =  3;
    if ( (etype == COMPLEX) && ( herPard) && ( defPard) ) mType_ =  4;
    if ( (etype == COMPLEX) && ( herPard) && (!defPard) ) mType_ = -4;
    if ( (etype == COMPLEX) && ( symPard) ) mType_ = 6;
    if ( (etype == DOUBLE ) && (!symPard) && (!strPard) ) mType_ = 11;
    if ( (etype == COMPLEX) && (!symPard) && (!strPard) && (!herPard)) {
      mType_ = 13;
    }
    if ( mType_ == 0 ) {
      (*error) << "PardisoSolver: There appears to be an inconsistency in "
               << "the input parameters. I cannot determine correct matrix "
               << "properties for pardiso";
      Error( __FILE__, __LINE__ );
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
      F77_FUNC(pardisoinit) ( pt_, &mType_, iparm_ + 1 );
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
        DeleteArray( idPerm_ );
        NewArray( idPerm_, int, probDim_ );
        for ( int i = 1; i <= probDim_; i++ ) {
          idPerm_[i] = i;
        }
      }
      break;

    default:
      (*error) << "Re-ordering of type '" << Enum2String(ordering)
               << "' is not available with the PardisoSolver";
      Error( __FILE__, __LINE__ );
    }

    if ( logging == true && facSymbolic == true ) {
      if ( ordering != NOREORDERING ) {
        (*cla) << " Pardiso: Analyse phase will determine a '"
               << Enum2String(ordering) << "' re-ordering" << std::endl;
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
                         (idPerm_+1), &nrhs, (iparm_+1), &msgLvl_, &zeroDBL_,
                         &zeroDBL_, &errorFlag );

      // Check return status
      if ( errorFlag != 0 ) {
        (*error) << "Pardiso: Error " << errorFlag << " occured during "
                 << "symbolic factorization";
        Error( __FILE__, __LINE__ ); 
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
                         (idPerm_+1), &nrhs, (iparm_+1), &msgLvl_, &zeroDBL_,
                         &zeroDBL_, &errorFlag );
        
      // Check return status
      if ( errorFlag != 0 ) {
        (*error) << "Pardiso: Error " << errorFlag << " occured during "
                 << "numerical factorization";
        Error( __FILE__, __LINE__ ); 
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
  }



  // *************************
  //   Solve linear system
  // *************************
  template<typename T>
  void PardisoSolver<T>::Solve( const BaseMatrix &sysmat,
                                const BasePrecond &precond,
                                const BaseVector &rhs, BaseVector &sol ) {

    ENTER_FCN( "PardisoSolver::Solve" );

    // Determine, whether we are expected to be verbose
    bool logging = myParams_->GetBoolValue( "PARDISO_logging" );
    if ( logging == true ) {
      (*cla) << " -------------------------------------------------------"
             << "-----------------------\n"
             << " Pardiso: Solving linear system ... ";
    }

    if ( firstCall_ == true ) {
      (*error) << "The matrix has not yet been factorised by Pardiso! "
               << "Call Setup() first";
      Error( __FILE__, __LINE__ );
    }
    
 
    // Check that we have the correct vector types and
    // obtain data pointers
    const T *rhsArray;
    T* solArray;
    TRY_CAST {
      ConstRefCast( rhs, Vector<T>, myRHS );
      RefCast( sol, Vector<T>, mySol );
      rhsArray = myRHS.GetPointer() + 1;
      solArray = mySol.GetPointer() + 1;
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

    F77_FUNC(pardiso) (pt_, &maxfct, &mnum, &mType_, &phase,
                       &probDim_, theMatrix_, rowPtr_, colPtr_,
                       (idPerm_+1), &nrhs, (iparm_+1), &msgLvl_, theRHS,
                       theSol, &errorFlag );

    // Check return status
    if ( errorFlag != 0 ) {
      (*error) << "Pardiso: Error " << errorFlag << " occured during "
               << "solution of linear system";
      Error( __FILE__, __LINE__ ); 
    }
    else {
      if ( logging == true ) {
        (*cla) << "done" << std::endl;
      }
    }

    // Finish log report
    if ( logging == true ) {
      (*cla) << " -------------------------------------------------------"
             << "-----------------------\n";
    }

    // Create Report (no sensible things to write for direct solvers yet)
    if ( myReport_ != NULL ) {
      myReport_->SetValue( "numIter", -1 );
      myReport_->SetValue( "finalNorm", -1.0 );
    }
  }

}
