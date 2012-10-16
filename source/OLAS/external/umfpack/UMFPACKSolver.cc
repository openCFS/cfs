#include <fstream>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/exception.hpp>
namespace fs = boost::filesystem;

#include <umfpack.h>

#include "MatVec/BaseMatrix.hh"
#include "MatVec/StdMatrix.hh"
#include "MatVec/CRS_Matrix.hh"
#include "MatVec/SCRS_Matrix.hh"

#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ProgramOptions.hh"

#include "def_use_pardiso.hh"
#include <def_fortran_interface.hh>

#include "UMFPACKSolver.hh"

namespace CoupledField {

  // Declare logging stream and make sure that it is also available in
  // release mode by using BOOST_DECLARE_LOG() instead of DECLARE_LOG()
  BOOST_DECLARE_LOG(umfpackSolver)
  DEFINE_LOG(umfpackSolver, "olas.solvers.pardiso")

#if PARDISO_API_VER == 3
extern "C" {
    void FORTRAN_CALL F77NAME(pardisoinit) (void *, int *, int *);

    void FORTRAN_CALL F77NAME(pardiso) (void *, int *, int *, int *, int *, int *,
                            const Double *, const int *, const int *, int *,
                            int *, int *, int *, const Double *,
                            Double *, int *);
}
#endif

#if PARDISO_API_VER == 4
extern "C" {
  void FORTRAN_CALL F77NAME(pardisoinit) (void *pt, int *mtype, int *solver,
                              int *iparm, Double *dparm, int *error);

  void FORTRAN_CALL F77NAME(pardiso) (void *pt, int *maxfct, int *mnum, int *mtype,
                           int *phase, int *n, const Double *a, const int *ia, const int *ja,
                           int *perm, int *nrhs, int *iparm, int *msglvl,
                           const Double *b, Double *x, int *error, double *dparm );
}
#endif

  // ***********************
  //   Default Constructor
  // ***********************
  template<typename T>
  UMFPACKSolver<T>::UMFPACKSolver() {
    EXCEPTION( "Default constructor of UMFPACKSolver is forbidden!" );
  }


  template<typename T>
  std::string UMFPACKSolver<T>::GetErrorString(int err_code) {
    switch (err_code) {
      case NO_ERROR:
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
  UMFPACKSolver<T>::UMFPACKSolver( PtrParamNode solverNode,
                                   PtrParamNode olasInfo ) {


    // Set pointers to communication objects
    xml_ = solverNode;
    infoNode_ = olasInfo->Get("pardiso",ParamNode::APPEND);

    // Initialise attributes
    firstCall_ = true;
    msgLvl_ = 0;

    // Resize data arrays for Pardiso
    pt_.Resize(64); pt_.Init(NULL);
    iparm_.Resize(64); iparm_.Init(0);
    dparm_.Resize(64); dparm_.Init(0.0);

    // Set default solver type to direct sparse solver
    std::string solverType = "direct";
    xml_->GetValue("type", solverType, ParamNode::INSERT);
      
    if(solverType == "direct") {
      mSolver_ = 0;
    } else {
      mSolver_ = 1;
    }

    PtrParamNode stopRulenode = xml_->Get("stoppingRule", ParamNode::INSERT);
    std::string sRule = "relNormRes0";
    stopRulenode->GetValue("type", sRule, ParamNode::INSERT);

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
  UMFPACKSolver<T>::~UMFPACKSolver() {

    umfpack_di_free_symbolic( &Symbolic );
    umfpack_di_free_numeric( &Numeric );

    // PARDISO - Last Phase: Cleaning up the parameters
    if ( firstCall_ == false ) {

      int errorFlag = 0;
      int phase = -1;

#if PARDISO_API_VER == 4
      F77NAME(pardiso) ( &pt_[0], &maxfct_, &mnum_, &mType_, &phase,
                         &probDim_, theMatrix_, rowPtr_, colPtr_,
                         idPerm_, &nrhs_, &iparm_[0], &msgLvl_, &zeroDBL_,
                         &zeroDBL_, &errorFlag, &dparm_[0] );
#endif

#if PARDISO_API_VER == 3
      F77NAME(pardiso) ( &pt_[0], &maxfct_, &mnum_, &mType_, &phase,
                         &probDim_, theMatrix_, rowPtr_, colPtr_,
                         idPerm_, &nrhs_, &iparm_[0], &msgLvl_, &zeroDBL_,
                         &zeroDBL_, &errorFlag );
#endif

      if ( errorFlag != NO_ERROR) {
        EXCEPTION( "Error occured during cleanup:\n"
                   << GetErrorString(errorFlag) )
      }

    // Read iterative solver statistics
    if(mSolver_ && msgLvl_ && fs::exists("pardiso-ml.out")) {
      std::ifstream pardisoMLOut("pardiso-ml.out", std::ifstream::binary);
      
      LOG_TRACE(umfpackSolver) << "Contents of pardiso-ml.out:"
                               << std::endl
                               << pardisoMLOut.rdbuf() << std::endl;
      LOG_TRACE(umfpackSolver) << " -------------------------------------------------------"
                               << "-----------------------";
      pardisoMLOut.close();

      try {
        fs::remove("pardiso-ml.out");
      } catch (std::exception &ex) {
        EXCEPTION("Error while trying to remove pardiso-ml.out: " << ex.what());
      }      
    }

    }

    // Delete identity re-ordering (if exists)
    delete [] ( idPerm_ );  idPerm_  = NULL;
    idPermSize_ = 0;

  }


  // *********
  //   Setup
  // *********
  template<typename T>
  void UMFPACKSolver<T>::Setup( BaseMatrix &sysMat, PtrParamNode analysis_step ) {

    // Flag for check Pardiso's return status
    //    int errorFlag = 0;

    // Determine, whether we are expected to be verbose
    LOG_TRACE(umfpackSolver) << " -----------------------------------------"
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

      bool newMatrixPattern = false;
      // When the matrix pattern has changed, we need to re-do
      // both steps, also the symbolical one
      if( newMatrixPattern ) {
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
    bool symPard, defPard, herPard, strPard;

    if ( stype == BaseMatrix::SPARSE_SYM ) {
      symPard = true;
    }
    else {
      symPard = false;
    }

    mType_ = 0;

    defPard = false;
    xml_->GetValue("posDef", defPard, ParamNode::INSERT);

    herPard = false;
    xml_->GetValue("hermitean", herPard, ParamNode::INSERT);

    strPard = false;
    xml_->GetValue("symStruct", strPard, ParamNode::INSERT);

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
      EXCEPTION( "UMFPACKSolver: There appears to be an inconsistency in "
               << "the input parameters. I cannot determine correct matrix "
               << "properties for pardiso" );
    }
    else {
      LOG_TRACE(umfpackSolver) << " Classified matrix as mType = " << mType_;
    }

    if(mSolver_ == 1) {
      LOG_TRACE(umfpackSolver) << " Using iterative solver";
      switch(mType_) {
        case 11:
        case 13:
          EXCEPTION( "UMFPACKSolver: The iterative solver just supports symmetric matrices!" );
          break;
      }
    } else {
      LOG_TRACE(umfpackSolver) << " Using direct solver";
    }

    // Set default input values for dparm_;
    dparm_[0] = 300;   // Maximum number of Krylov-subspace iterations
    xml_->GetValue("maxIter", dparm_[0], ParamNode::INSERT);

    dparm_[1] = 1e-6;  // Relative residual reduction
    xml_->GetValue("tol", dparm_[1], ParamNode::INSERT);

    dparm_[2] = 1e-6;  // Coarse Grid Matrix Dimension.
    xml_->GetValue("coarseGridDim", dparm_[2], ParamNode::INSERT);

    dparm_[3] = 10;    // Maximum Number of Grid Levels.
    xml_->GetValue("maxNumGridLevels", dparm_[3], ParamNode::INSERT);

    dparm_[4] = 1e-2;  // Dropping value for the incomplete factor.
    xml_->GetValue("incompFacDropVal", dparm_[4], ParamNode::INSERT);

    dparm_[5] = 5e-5;  // Dropping value for the schurcomplement.
    xml_->GetValue("schurcompDropVal", dparm_[5], ParamNode::INSERT);

    dparm_[6] = 10;    // Maximum number of ﬁll-in in each column in the factor.
    xml_->GetValue("maxNumFillIn", dparm_[6], ParamNode::INSERT);

    dparm_[7] = 500;   // Bound for the inverse of the incomplete factor L.
    xml_->GetValue("invBoundIncompFac", dparm_[0], ParamNode::INSERT);

    dparm_[8] = 25;    // Maximum number of non-improvement steps in Krylov-Subspace method
    xml_->GetValue("maxNumStagnationSteps", dparm_[0], ParamNode::INSERT);

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
      LOG_TRACE(umfpackSolver) << " Calling pardisoinit";

#if 0

#if PARDISO_API_VER == 4
      Integer error = 0;
      F77NAME(pardisoinit) ( &pt_[0],  &mType_, &mSolver_, &iparm_[0], &dparm_[0], &error);

      if (error != 0) 
        EXCEPTION(GetErrorString(error));
#endif

#if PARDISO_API_VER == 3
      F77NAME(pardisoinit) ( &pt_[0], &mType_, &iparm_[0] );
#endif

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
    xml_->GetValue("ordering", orderStr, ParamNode::INSERT);
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
               << "' is not available with the UMFPACKSolver" );
    }

    if(!mSolver_) {
      if ( facSymbolic == true ) {
        if ( ordering != BaseOrdering::NOREORDERING ) {
          std::string tmp;
      tmp = BaseOrdering::reorderingType.ToString( ordering );

          LOG_TRACE(umfpackSolver) << " Analyse phase will determine a '"
                                   << tmp << "' re-ordering";
        }
        else {
          LOG_TRACE(umfpackSolver) << " Factorisation uses original matrix ordering";
        }
      }
    }
    
    // Do we need to determine MFLOPs for the LU factorisation
    bool stats = false;
      xml_->GetValue("stats", stats, ParamNode::INSERT);
    
    if(stats)
      iparm_[18] = -1;
    else
      iparm_[18] = 0;
  
    // Setting pivoting strategy for indefinit problems
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

    // Switch on out-of-core
    if (xml_->Has("outOfCore"))
    {
      if(std::string(CFS_PARDISO) != "MKL") 
      {
        WARN( "The value of outOfCore has no effect for " << CFS_PARDISO << ".\n"
              << "Switch to MKL if want to use out-of-core memory." );
      }
      xml_->GetValue("outOfCore", iparm_[59], ParamNode::INSERT);
    }

    // Switch to iterative solver
    if(mSolver_) {
#if PARDISO_API_VER == 3
      EXCEPTION("This CFS++ executable has been linked to a PARDISO 3.x library.\n"
                << "PARDISO implements iterative solvers since version 4.0.\n"
                << "To get iterative solvers switch CFS_PARDISO to SCHENK.");
#endif      
      
      iparm_[31] = 1;
    }

    /* The facotrisation algorithm in pardiso rows/columns which have the about
     * the same magnitude and structure into one supernode. This may cause
     * numerical errors due to pivotisation, the iterative refine step remedies
     * this. (see O.Schenk et al "Solving unsymmetric sparse systems of linear
     * equations with PARDISO")
     */
    if (xml_->Has("IterRefineSteps"))
    {
      xml_->GetValue("IterRefineSteps", iparm_[7], ParamNode::INSERT);
    } else if (iparm_[7] == 0 && mType_ > 10) {
      WARN( "Iterative refinement steps is not set!\n" 
            << "PARDISO did not set it neither.\n"
            << "It is advisable for unsymmetric to set it in the xml,\n"
            << "since bad pivotisation may result in numerical errors!" );
    }
    
    // we have to icrement the entries of the col- and row-position arrays
    // by one, so that the first col and first row start with index 1 (and
    // not with zero) to be consistent with fortran
    // at the end of the method we will undo it!!
    //    for (UInt i=0; i < static_cast<UInt>(probDim_+1); i++ )
    //    rowPtr_[i] += 1;

    // for (UInt i=0; i< nnz_; i++ )
    //   colPtr_[i] += 1;

    // ========================
    //  Symbolic Factorisation
    // ========================
    if ( facSymbolic == true ) {

      // log report
      LOG_TRACE(umfpackSolver) << " Performing analyse phase (symbolic factorisation)"
                               << " ... ";

      umfpack_di_symbolic (probDim_, probDim_, rowPtr_, colPtr_, theMatrix_, &Symbolic, 0, 0) ;

#if 0
      // only analyse
      int phase = 11;

      // let pardiso go for it
#if PARDISO_API_VER == 4
      F77NAME(pardiso) (&pt_[0], &maxfct_, &mnum_, &mType_, &phase,
                        &probDim_, theMatrix_, rowPtr_, colPtr_,
                        idPerm_, &nrhs_, &iparm_[0], &msgLvl_, &zeroDBL_,
                        &zeroDBL_, &errorFlag, &dparm_[0] );
#endif

#if PARDISO_API_VER == 3
      F77NAME(pardiso) (&pt_[0], &maxfct_, &mnum_, &mType_, &phase,
                        &probDim_, theMatrix_, rowPtr_, colPtr_,
                        idPerm_, &nrhs_, &iparm_[0], &msgLvl_, &zeroDBL_,
                        &zeroDBL_, &errorFlag );
#endif

      // Check return status
      if ( errorFlag != NO_ERROR ) {
        EXCEPTION( "Error occured during symbolic factorization:\n"
                   << GetErrorString(errorFlag));
      }
      else {
        LOG_TRACE(umfpackSolver) << "done";
      }
#endif

    }


    // =========================
    //  Numerical Factorisation
    // =========================
    if ( facNumeric == true ) {

      // log report
      LOG_TRACE(umfpackSolver) << " Performing factorise phase (numerical "
                               << "factorisation) ... ";


      umfpack_di_numeric (rowPtr_, colPtr_, theMatrix_, Symbolic, &Numeric, 0, 0);

      // only factorise (numerical)
#if 0
      int phase = 22;



      // let pardiso go for it
#if PARDISO_API_VER == 4
      F77NAME(pardiso) (&pt_[0], &maxfct_, &mnum_, &mType_, &phase,
                        &probDim_, theMatrix_, rowPtr_, colPtr_,
                        idPerm_, &nrhs_, &iparm_[0], &msgLvl_, &zeroDBL_,
                        &zeroDBL_, &errorFlag, &dparm_[0] );
#endif

#if PARDISO_API_VER == 3
      F77NAME(pardiso) (&pt_[0], &maxfct_, &mnum_, &mType_, &phase,
                        &probDim_, theMatrix_, rowPtr_, colPtr_,
                        idPerm_, &nrhs_, &iparm_[0], &msgLvl_, &zeroDBL_,
                        &zeroDBL_, &errorFlag );
#endif
      
      // Check return status
      if ( errorFlag != NO_ERROR ) {
        EXCEPTION( "Error occured during numerical factorization:\n"
                   << GetErrorString(errorFlag) );
      }
      else {
        LOG_TRACE(umfpackSolver) << "done";
      }

#endif

    }

    // Now we were called once, and a factorisation is available
    firstCall_ = false;

#if 0
    // now we undo our increment, since on our side the frist col and row
    // has an value of zero!!
    for (UInt i=0; i <  static_cast<UInt>(probDim_+1); i++ )
      rowPtr_[i] -= 1;

    for (UInt i=0; i< nnz_; i++ )
      colPtr_[i] -= 1;
#endif

  }



  // *************************
  //   Solve linear system
  // *************************
  template<typename T>
  void UMFPACKSolver<T>::Solve( const BaseMatrix &sysmat,
                                const BaseVector &rhs, BaseVector &sol, PtrParamNode analysis_step ) {

    LOG_TRACE(umfpackSolver) << " -----------------------------------------"
                             << "-------------------------------------";
    LOG_TRACE(umfpackSolver) << " Solving linear system ...";

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


    umfpack_di_solve (UMFPACK_A, rowPtr_, colPtr_, theMatrix_, theSol, theRHS, Numeric, 0, 0) ;

#if 0
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

#if PARDISO_API_VER == 4
    F77NAME(pardiso) (&pt_[0], &maxfct_, &mnum_, &mType_, &phase,
                      &probDim_, theMatrix_, rowPtr_, colPtr_,
                      idPerm_, &nrhs_, &iparm_[0], &msgLvl_, theRHS,
                      theSol, &errorFlag, &dparm_[0] );
#endif

#if PARDISO_API_VER == 3
    F77NAME(pardiso) (&pt_[0], &maxfct_, &mnum_, &mType_, &phase,
                      &probDim_, theMatrix_, rowPtr_, colPtr_,
                      idPerm_, &nrhs_, &iparm_[0], &msgLvl_, theRHS,
                      theSol, &errorFlag );
#endif

    // now we undo our increment, since on our side the first col and row
    // has an value of zero!!
    for (UInt i=0; i <  static_cast<UInt>(probDim_+1); i++ )
      rowPtr_[i] -= 1;
    for (UInt i=0; i< nnz_; i++ )
      colPtr_[i] -= 1;

    // Check return status
    if ( errorFlag != NO_ERROR ) {
      EXCEPTION( "Error occured during solution of linear system:\n"
                 << GetErrorString(errorFlag) );
    }
    else {
      LOG_TRACE(umfpackSolver) << "done";
    }

#endif


    // Finish log report
    if(!mSolver_) {
      LOG_TRACE(umfpackSolver) << " number of iterative refinement steps: " << iparm_[6];
      LOG_TRACE(umfpackSolver) << " number of perturbed pivots: " << iparm_[13];
      LOG_TRACE(umfpackSolver) << " number of positive eigenvalues: " << iparm_[21];
      LOG_TRACE(umfpackSolver) << " number of negative eigenvalues: " << iparm_[22];
    } else {
      LOG_TRACE(umfpackSolver) << " Maximum number of iterations: " << dparm_[0];
      LOG_TRACE(umfpackSolver) << " Number of iterations after solve step: " << dparm_[34];
      LOG_TRACE(umfpackSolver) << " Relative Residual Reduction: " << dparm_[1];
      LOG_TRACE(umfpackSolver) << " Relative residual after Krylov-Subspace convergence: " << dparm_[33];
    }
    LOG_TRACE(umfpackSolver) << " -------------------------------------------------------"
                             << "-----------------------";

    // Create Report (no sensible things to write for direct solvers yet)
    ParamNode::ActionType at = progOpts->DoDetailedInfo() ? ParamNode::APPEND : ParamNode::DEFAULT;
    PtrParamNode out = infoNode_->Get(ParamNode::PROCESS)->Get("solver", at);
    out->Get("numIter")->SetValue(-1);
    out->Get("finalNorm")->SetValue(-1.0);
  }

// Explicit template instantiation
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template class UMFPACKSolver<Double>;
  template class UMFPACKSolver<Complex>;
#endif

}
