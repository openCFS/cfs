#include <fstream>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/exception.hpp>
namespace fs = boost::filesystem;

#include "MatVec/BaseMatrix.hh"
#include "MatVec/StdMatrix.hh"
#include "MatVec/CRS_Matrix.hh"
#include "MatVec/SCRS_Matrix.hh"

#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ProgramOptions.hh"
#include "Domain/Domain.hh"

#include "def_use_pardiso.hh"
#include <def_cfs_fortran_interface.hh>

#include "PardisoSolverPrimitive.hh"

namespace CoupledField {

DEFINE_LOG(pardisoSolverPrim, "olas.solvers.pardiso")


#if PARDISO_API_VER == 3

extern "C" {
    void pardisoinit (void *, int *, int *);

    void pardiso (void *, int *, int *, int *, int *, int *,
                  const Double *, const int *, const int *, int *,
                  int *, int *, int *, const Double *,
                  Double *, int *);
}
#endif

#if PARDISO_API_VER == 4
extern "C" {
  void pardisoinit (void *pt, int *mtype, int *solver,
                              int *iparm, Double *dparm, int *error);

  void pardiso (void *pt, int *maxfct, int *mnum, int *mtype,
                int *phase, int *n, const Double *a, const int *ia, const int *ja,
                int *perm, int *nrhs, int *iparm, int *msglvl,
                const Double *b, Double *x, int *error, double *dparm );
}
#endif

  // ***********************
  //   Default Constructor
  // ***********************
  template<typename T>
  PardisoSolverPrimitive<T>::PardisoSolverPrimitive() {
    EXCEPTION( "Default constructor of PardisoSolverPrimitive is forbidden!" );
  }


  template<typename T>
  std::string PardisoSolverPrimitive<T>::GetErrorString(int err_code) {
    switch (err_code) {
      case PARDISO_NO_ERROR:
        return "No error.";
      case INPUT_INCONSISTENT:
        return "Input inconsistent.";
      case NOT_ENOUGH_MEMORY:
        return "Not enough memory. Try to switch to out-of-core mode if using MKL.";
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
      case NOT_ENOUGH_OOC_MEM:
        return "Not enough memory for out-of-core. Try increasing MKL_PARDISO_OOC_MAX_CORE_SIZE in ./pardiso_ooc.cfg.";
      case NO_LIC_FILE:
        {
          std::string msg;
          msg = "No license file pardiso.lic found.\n";
          msg += "Please get the file at http://www.pardiso-project.org\n";
          msg += "and set the environment variable PARDISO_LIC_PATH to\n";
          msg += "the directory containing pardiso.lic.\n";
          return msg;
        }
      case LIC_EXPIRED:
        return "License is expired.";
      case WRONG_USER_OR_HOSTNAME:
        return "Wrong username or hostname.";
      case MAX_KRYLOV_ITERATIONS:
        return "Reached maximum number of Krylov-subspace iterations in iterative solver.";
      case INSUFF_KRYLOV_CONVERGENCE:
        return "No sufficient convergence in Krylov-subspace iteration within 25 iterations.";
      case KRYLOV_ITERATION_ERROR:
        return "Error in Krylov-subspace iteration.";
      case KRYLOV_BREAKDOWN:
        return "Break-Down in Krylov-subspace iteration.";
      
      default:
        return "Unclassified (internal) error.";
    }
  }

  // ***************
  //   Constructor
  // ***************
  template<typename T>
  PardisoSolverPrimitive<T>::PardisoSolverPrimitive( PtrParamNode solverNode,
                                   PtrParamNode olasInfo ) {


    // Initialise attributes
    firstCall_ = true;
    msgLvl_ = 0;

    // Resize data arrays for Pardiso
    pt_.Resize(64); pt_.Init(NULL);
    iparm_.Resize(64); iparm_.Init(0);
    dparm_.Resize(64); dparm_.Init(0.0);

    // Set default solver type to direct sparse solver
    std::string solverType = "direct";
      
    if(solverType == "direct") {
      mSolver_ = 0;
    } else {
      mSolver_ = 1;
    }

    std::string sRule = "relNormRes0";

    if(mSolver_ && sRule != "relNormRes0") {
      WARN("The iterative solver in PARDISO only supports relative " \
           "residual minimization as stopping rule");
    }

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
  PardisoSolverPrimitive<T>::~PardisoSolverPrimitive() {

    // PARDISO - Last Phase: Cleaning up the parameters
    if ( firstCall_ == false ) {

      int errorFlag = 0;
      int phase = -1;

#if PARDISO_API_VER == 4
      pardiso ( &pt_[0], &maxfct_, &mnum_, &mType_, &phase,
                &probDim_, theMatrix_, rowPtr_, colPtr_,
                idPerm_, &nrhs_, &iparm_[0], &msgLvl_, &zeroDBL_,
                &zeroDBL_, &errorFlag, &dparm_[0] );
#endif

#if PARDISO_API_VER == 3
      pardiso ( &pt_[0], &maxfct_, &mnum_, &mType_, &phase,
                &probDim_, theMatrix_, rowPtr_, colPtr_,
                idPerm_, &nrhs_, &iparm_[0], &msgLvl_, &zeroDBL_,
                &zeroDBL_, &errorFlag );
#endif

      if ( errorFlag != PARDISO_NO_ERROR) {
        std::cout << "Error occured during cleanup: " << GetErrorString(errorFlag) << std::endl;
        domain->GetInfoRoot()->ToFile();
        exit(-1); // no exceptions in desctructors
      }

    // Read iterative solver statistics
    if(mSolver_ && msgLvl_ && fs::exists("pardiso-ml.out")) {
      std::ifstream pardisoMLOut("pardiso-ml.out", std::ifstream::binary);
      
      LOG_TRACE(pardisoSolverPrim) << "Contents of pardiso-ml.out:"
                               << std::endl
                               << pardisoMLOut.rdbuf() << std::endl;
      LOG_TRACE(pardisoSolverPrim) << " -------------------------------------------------------"
                               << "-----------------------";
      pardisoMLOut.close();

      try {
        fs::remove("pardiso-ml.out");
      } catch (std::exception &ex) {
        std::cout << "Error while trying to remove pardiso-ml.out: " << ex.what() << std::endl;
        domain->GetInfoRoot()->ToFile();
        exit(-1); // no exceptions in desctructors
      }      
    }

    }

    // Delete identity re-ordering (if exists)
    delete [] ( idPerm_ );  idPerm_  = NULL;
    idPermSize_ = 0;

  }


  // ******************
  // Set matrix pattern
  // ******************
  template<typename T>
  void PardisoSolverPrimitive<T>::SetNewMatrixPattern(void) {
    newMatrixPattern_ = true;
  }


  // *********
  //   Setup
  // *********
  template<typename T>
  void PardisoSolverPrimitive<T>::Setup( BaseMatrix &sysMat) {

    // Flag for check Pardiso's return status
    int errorFlag = 0;

    // Determine, whether we are expected to be verbose
    LOG_TRACE(pardisoSolverPrim) << " -----------------------------------------"
                             << "-------------------------------------";

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

      // TODO: THIS CHECK DOES NOT MAKE SENSE IN MY OPINION SINCE
      //       'newMatrixPattern' is set to false in olasparams.cc
      //       and gets never changed elsewhere. A more intelligent
      //       test would ask the matrix if its pattern did change.

      //bool newMatrixPattern = false;
      // When the matrix pattern has changed, we need to re-do
      // both steps, also the symbolical one
      if( newMatrixPattern_ ) {
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
    const StdMatrix& stdMat = dynamic_cast<const StdMatrix&>(sysMat);

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

    // ==================================================
    //  Get pointers to arrays containing information on
    //  (S)CRS matrix from the problem matrix object
    // ==================================================
    if ( stype == BaseMatrix::SPARSE_NONSYM ) {
      const CRS_Matrix<T>& crsMat = dynamic_cast<const CRS_Matrix<T>&>(sysMat);
      rowPtr_ = (Integer*)(crsMat.GetRowPointer());
      colPtr_ = (Integer*)(crsMat.GetColPointer());
      datPtr_ = crsMat.GetDataPointer();
      nnz_    = crsMat.GetNnz();
    }
    else {
      const SCRS_Matrix<T>& scrsMat = dynamic_cast<const SCRS_Matrix<T>&>(sysMat);
      rowPtr_ = (Integer*)(scrsMat.GetRowPointer());
      colPtr_ = (Integer*)(scrsMat.GetColPointer());
      datPtr_ = scrsMat.GetDataPointer();
      nnz_    = scrsMat.GetNumEntries();
    }

    // Do C-style casting to fix problem with over-loading vs. extern C
    void *hatred = (void *) datPtr_;
    theMatrix_ = (Double *) hatred;


    // ====================================
    //  Set default parameters for Pardiso
    // ====================================

    // Some flags for determining the type of the matrix
    bool defPard, herPard, strPard;

    bool symPard = stype == BaseMatrix::SPARSE_SYM ? true : false;

    mType_ = 0;

    defPard = false;

    herPard = false;

    strPard = false;

    // see pardiso manual -> MTYPE, page 9
    if(etype == BaseMatrix::DOUBLE)
    {
      if(strPard)
        mType_ = 1;                 // 1: real and structurally symmetric -- wtf is this???
      if(symPard)
        mType_ = defPard ? 2 : -2;  // 2: real and symmetric positive definite, -2: real and symmetric indefinite
      if(!symPard && !strPard)
        mType_ = 11;                // 11: real and nonsymmetric
    }
    else
    {
      if(herPard)
        mType_ = defPard ? 4 : -4; // 4: complex and Hermitian positive definite, -4 : complex and Hermitian indefinite
      if(strPard)
        mType_ = 3;                // 3: complex and structurally symmetric
      if(!strPard && !herPard)
        mType_ = symPard ? 6 : 13; // 6: complex and symmetric, 13: complex  and nonsymmetric
    }

    if(mType_ == 0)
      throw Exception("Pardiso matrix type cannot be determined.");

    if(mSolver_ == 1 && (mType_ == 11 || mType_ == 13))
      throw Exception("The iterative pardiso solver just supports symmetric matrices!" );

    // Set default input values for dparm_;
    dparm_[0] = 300;   // Maximum number of Krylov-subspace iterations

    dparm_[1] = 1e-6;  // Relative residual reduction

    dparm_[2] = 1e-6;  // Coarse Grid Matrix Dimension.

    dparm_[3] = 10;    // Maximum Number of Grid Levels.

    dparm_[4] = 1e-2;  // Dropping value for the incomplete factor.

    dparm_[5] = 5e-5;  // Dropping value for the schurcomplement.

    dparm_[6] = 10;    // Maximum number of ﬁll-in in each column in the factor.

    dparm_[7] = 500;   // Bound for the inverse of the incomplete factor L.

    dparm_[8] = 25;    // Maximum number of non-improvement steps in Krylov-Subspace method

    // Remove output file for iterative solver
    if(mSolver_ && fs::exists("pardiso-ml.out")) {
      try {
        fs::remove("pardiso-ml.out");
      } catch (std::exception &ex) {
        EXCEPTION("Error while trying to remove pardiso-ml.out: " << ex.what());
      }      
    }
    
    // We do not need to call pardisoinit, if the matrix pattern
    // has not changed
    if ( facSymbolic ) {
      LOG_TRACE(pardisoSolverPrim) << " Calling pardisoinit";

#if PARDISO_API_VER == 4
      Integer error = 0;
      pardisoinit ( &pt_[0],  &mType_, &mSolver_, &iparm_[0], &dparm_[0], &error);

      if (error != 0) 
        EXCEPTION(GetErrorString(error));
#endif

#if PARDISO_API_VER == 3
      pardisoinit ( &pt_[0], &mType_, &iparm_[0] );
#endif

    }


    // =======================================
    //  Alter default parameters to our taste
    // =======================================

    // Avoid that Pardiso over-writes our settings
    iparm_[0] = 1;

    // Specify number of OpenMP threads. OLAS currently has no consistent
    // concept for shared memory parallelism, so we specify one thread here.
    //#if defined (_OPENMP)
    //iparm_[2] = omp_get_num_threads();
    //#else
    //iparm_[2] = 1;
    //#endif

    // Determine the re-ordering strategy: We can either fo nested dissection
    // or minimum degree re-ordering or no re-ordering at all (i.e. we use
    // the initial ordering of the linear system, which might already have been
    // re-ordered via the graph)
    BaseOrdering::ReorderingType ordering = BaseOrdering::NESTED_DISSECTION;
    std::string orderStr = "nestedDissection";
    ordering = BaseOrdering::reorderingType.Parse( orderStr );

    switch ( ordering ) {

    case BaseOrdering::NESTED_DISSECTION:
      iparm_[1] = 2;
      iparm_[4] = 0;
      break;

    case BaseOrdering::MINIMUM_DEGREE:
      iparm_[1] = 0;
      iparm_[4] = 0;
      break;

    case BaseOrdering::NOREORDERING:
      // In this case iparm_[1] is irrelevant, we generate an identity
      // permutation and use this one, by setting iparm_[4] and skip the
      // symbolic factorisation
      iparm_[1] = 0;
      iparm_[4] = 1;
      if ( idPermSize_ < probDim_ ) {
        delete [] ( idPerm_ );  idPerm_  = NULL;
        NEWARRAY( idPerm_, int, probDim_ );
        for ( int i = 0; i < probDim_; i++ ) {
          // (i+1), since fortran needs indices starting with 1!!
          idPerm_[i] = i+1;
        }
      }
      break;

    default:
      std::string tmp;
      tmp = BaseOrdering::reorderingType.ToString( ordering );

      EXCEPTION( "Re-ordering of type '" << tmp
               << "' is not available with the PardisoSolverPrimitive" );
    }

    if(!mSolver_) {
      if ( facSymbolic == true ) {
        if ( ordering != BaseOrdering::NOREORDERING ) {
          std::string tmp;
          tmp = BaseOrdering::reorderingType.ToString( ordering );

          LOG_TRACE(pardisoSolverPrim) << " Analyse phase will determine a '"
                                   << tmp << "' re-ordering";
        }
        else {
          LOG_TRACE(pardisoSolverPrim) << " Factorisation uses original matrix ordering";
        }
      }
    }
    
    // Do we need to determine MFLOPs for the LU factorisation
    bool stats = false;
    
    if(stats)
      iparm_[18] = -1;
    else
      iparm_[18] = 0;
  
    // Setting pivoting strategy for indefinite problems
    iparm_[20] = 1;

    // In case we have no positive definite system (especially piezo)
    // we perform additional scaling to enhance the condition for very
    // small off-diagonal entries (iparm_[11]). In addition we enable
    // the method of 'symmetric weighted matchings' (iparam_[13]).
    // For further information, refer to the pardiso user manual.
    if( !defPard ) {
      iparm_[10] = 1;
      iparm_[12] = 1;
    }

    // Pardiso keeps one factorisation in memory (and that is used for
    // the solution phase)
    maxfct_ = 1;
    mnum_ = 1;

    // We simultaneously treat a single right hand side
    nrhs_ = 1;

    // Set the message level (for printing statistics)
    msgLvl_ = 0;
    if ( stats ) {
      msgLvl_ = 1;
    }

    // Switch to iterative solver
    if(mSolver_) {
#if PARDISO_API_VER == 3
      EXCEPTION("This openCFS executable has been linked to a PARDISO 3.x library.\n"
                << "PARDISO implements iterative solvers since version 4.0.\n"
                << "To get iterative solvers you should switch CFS_PARDISO to SCHENK.");
#endif      
      
      iparm_[31] = 1;
    }

    
    // we have to increment the entries of the col- and row-position arrays
    // by one, so that the first col and first row start with index 1 (and
    // not with zero) to be consistent with fortran
    // at the end of the method we will undo it!!
    for (UInt i=0; i < static_cast<UInt>(probDim_+1); ++i )
      rowPtr_[i] += 1;

    for (UInt i=0; i< nnz_; i++ )
       colPtr_[i] += 1;



    // ========================
    //  Symbolic Factorisation
    // ========================
    if ( facSymbolic == true ) {

      tSymfact_.ResetStart();
      // log report
      LOG_TRACE(pardisoSolverPrim) << " Performing analyse phase (symbolic factorisation)"
                               << " ... ";

      // only analyse
      int phase = 11;

      // let pardiso go for it
#if PARDISO_API_VER == 4
      pardiso (&pt_[0], &maxfct_, &mnum_, &mType_, &phase,
               &probDim_, theMatrix_, rowPtr_, colPtr_,
               idPerm_, &nrhs_, &iparm_[0], &msgLvl_, &zeroDBL_,
               &zeroDBL_, &errorFlag, &dparm_[0] );
#endif

#if PARDISO_API_VER == 3
      pardiso (&pt_[0], &maxfct_, &mnum_, &mType_, &phase,
               &probDim_, theMatrix_, rowPtr_, colPtr_,
               idPerm_, &nrhs_, &iparm_[0], &msgLvl_, &zeroDBL_,
               &zeroDBL_, &errorFlag );
#endif

      // Check return status
      if ( errorFlag != PARDISO_NO_ERROR ) {
        EXCEPTION( "Error occured during symbolic factorization:\n"
                   << GetErrorString(errorFlag));
      }
      else {
        LOG_TRACE(pardisoSolverPrim) << "done";
      }

      tSymfact_.Stop();
    }
    // =========================
    //  Numerical Factorisation
    // =========================
    if ( facNumeric == true ) {

      tNumfact_.ResetStart();
      // log report
      LOG_TRACE(pardisoSolverPrim) << " Performing factorise phase (numerical "
                               << "factorisation) ... ";

      // only factorise (numerical)
      int phase = 22;

      // let pardiso go for it
#if PARDISO_API_VER == 4
      pardiso (&pt_[0], &maxfct_, &mnum_, &mType_, &phase,
               &probDim_, theMatrix_, rowPtr_, colPtr_,
               idPerm_, &nrhs_, &iparm_[0], &msgLvl_, &zeroDBL_,
               &zeroDBL_, &errorFlag, &dparm_[0] );
#endif

#if PARDISO_API_VER == 3
      pardiso (&pt_[0], &maxfct_, &mnum_, &mType_, &phase,
               &probDim_, theMatrix_, rowPtr_, colPtr_,
               idPerm_, &nrhs_, &iparm_[0], &msgLvl_, &zeroDBL_,
               &zeroDBL_, &errorFlag );
#endif

      // Check return status
      if ( errorFlag != PARDISO_NO_ERROR ) {
        EXCEPTION( "Error occured during numerical factorization:\n"
                   << GetErrorString(errorFlag) );
      }
      else {
        LOG_TRACE(pardisoSolverPrim) << "done";
        LOG_TRACE(pardisoSolverPrim) << "Memory consumption during numerical factorization and solution: " << iparm_[16] << "kBytes";
      }

      tNumfact_.Stop();;
    }


    // Now we were called once, and a factorisation is available
    firstCall_ = false;

    // now we undo our increment, since on our side the first col and row
    // has an value of zero!!
    for (UInt i=0; i <  static_cast<UInt>(probDim_+1); ++i )
      rowPtr_[i] -= 1;

    for (UInt i=0; i< nnz_; i++ )
      colPtr_[i] -= 1;

  }



  // *************************
  //   Solve linear system
  // *************************
  template<typename T>
  void PardisoSolverPrimitive<T>::Solve( const BaseMatrix &sysmat, const BaseVector &rhs, BaseVector &sol) {

    LOG_TRACE(pardisoSolverPrim) << " -----------------------------------------"
                             << "-------------------------------------";
    LOG_TRACE(pardisoSolverPrim) << " Solving linear system ...";

    if ( firstCall_ == true ) {
      EXCEPTION( "The matrix has not yet been factorised by Pardiso! "
               << "Call Setup() first" );
    }


    // Check that we have the correct vector types and
    // obtain data pointers
    const T *rhsArray;
    T* solArray;
    const Vector<T>& myRHS = dynamic_cast<const Vector<T>&>(rhs);
    Vector<T>& mySol = dynamic_cast<Vector<T>&>(sol);
    rhsArray = myRHS.GetPointer();
    solArray = mySol.GetPointer();

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

    // we have to increment the entries of the col- and row-position arrays
    // by one, so that the first col and first row start with index 1 (and
    // not with zero) to be consistent with fortran
    // at the end of the method we will undo it!!
    for (UInt i=0; i<static_cast<UInt>(probDim_+1); ++i )
      rowPtr_[i] += 1;
    for (UInt i=0; i< nnz_; i++ )
       colPtr_[i] += 1;

#if PARDISO_API_VER == 4
    pardiso (&pt_[0], &maxfct_, &mnum_, &mType_, &phase,
             &probDim_, theMatrix_, rowPtr_, colPtr_,
             idPerm_, &nrhs_, &iparm_[0], &msgLvl_, theRHS,
             theSol, &errorFlag, &dparm_[0] );
#endif

#if PARDISO_API_VER == 3
    pardiso (&pt_[0], &maxfct_, &mnum_, &mType_, &phase,
             &probDim_, theMatrix_, rowPtr_, colPtr_,
             idPerm_, &nrhs_, &iparm_[0], &msgLvl_, theRHS,
             theSol, &errorFlag );
#endif

    // now we undo our increment, since on our side the first col and row
    // has an value of zero!!
    for (UInt i=0; i <  static_cast<UInt>(probDim_+1); ++i )
      rowPtr_[i] -= 1;
    for (UInt i=0; i< nnz_; i++ )
      colPtr_[i] -= 1;

    // Check return status
    if ( errorFlag != PARDISO_NO_ERROR ) {
      EXCEPTION( "Error occured during solution of linear system:\n"
                 << GetErrorString(errorFlag) );
    }
    else {
      LOG_TRACE(pardisoSolverPrim) << "done";
    }

    // Finish log report
    if(!mSolver_) {
      if ( iparm_[7] > 0 && iparm_[6] == iparm_[7] ) {
        WARN( "PARDISO reached the maximum number of iterative refinement steps ("
              << iparm_[7] << ").\n This could indicate that the solution is incorrect!\n"
              << " You may need to increase the iterRefineSteps parameter in"
              << " order to obtain\n a correct solution." );
      }
      LOG_TRACE(pardisoSolverPrim) << " number of iterative refinement steps: " << iparm_[6];
      LOG_TRACE(pardisoSolverPrim) << " number of perturbed pivots: " << iparm_[13];
      LOG_TRACE(pardisoSolverPrim) << " number of positive eigenvalues: " << iparm_[21];
      LOG_TRACE(pardisoSolverPrim) << " number of negative eigenvalues: " << iparm_[22];
    } else {
      LOG_TRACE(pardisoSolverPrim) << " Maximum number of iterations: " << dparm_[0];
      LOG_TRACE(pardisoSolverPrim) << " Number of iterations after solve step: " << dparm_[34];
      LOG_TRACE(pardisoSolverPrim) << " Relative Residual Reduction: " << dparm_[1];
      LOG_TRACE(pardisoSolverPrim) << " Relative residual after Krylov-Subspace convergence: " << dparm_[33];
    }
    LOG_TRACE(pardisoSolverPrim) << " -------------------------------------------------------"
                             << "-----------------------";

  }

  // Explicit template instantiation
  template class PardisoSolverPrimitive<Double>;
  template class PardisoSolverPrimitive<Complex>;

}
