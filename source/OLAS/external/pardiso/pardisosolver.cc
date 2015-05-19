// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <stddef.h>
#include <exception>
#include <fstream>

#include "boost/filesystem/exception.hpp"
#include "boost/filesystem/operations.hpp"

namespace fs = boost::filesystem;

#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"
#include "DataInOut/programOptions.hh"
#include "General/Enum.hh"
#include "General/exception.hh"
#include "MatVec/basematrix.hh"
#include "MatVec/crs_matrix.hh" // IWYU pragma: keep
#include "MatVec/scrs_matrix.hh" // IWYU pragma: keep
#include "MatVec/stdmatrix.hh"
#include "OLAS/graph/baseordering.hh"
#include "def_use_pardiso.hh"
#include "pardisosolver.hh"

namespace CoupledField {
class BasePrecond;
class BaseVector;
template <typename T> class Vector;

  // Declare logging stream and make sure that it is also available in
  // release mode by using BOOST_DECLARE_LOG() instead of DECLARE_LOG()
  BOOST_DECLARE_LOG(pardisoSolver)
  DEFINE_LOG(pardisoSolver, "olas.solvers.pardiso")

#define F77_FUNC(func)   func ## _

#if PARDISO_API_VER == 3
extern "C" {
    void F77_FUNC(pardisoinit) (void *, int *, int *);

    void F77_FUNC(pardiso) (void *, int *, int *, int *, int *, int *,
                            const Double *, const int *, const int *, int *,
                            int *, int *, int *, const Double *,
                            Double *, int *);
}
#endif

#if PARDISO_API_VER == 4
extern "C" {
  void F77_FUNC(pardisoinit) (void *pt, int *mtype, int *solver,
                              int *iparm, Double *dparm, int *error);

  void F77_FUNC(pardiso) (void *pt, int *maxfct, int *mnum, int *mtype,
                           int *phase, int *n, const Double *a, const int *ia, const int *ja,
                           int *perm, int *nrhs, int *iparm, int *msglvl,
                           const Double *b, Double *x, int *error, double *dparm );
}
#endif

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
  PardisoSolver<T>::PardisoSolver( PtrParamNode solverNode,
                                   PtrParamNode olasInfo ) {


    // Set pointers to communication objects
    xml_ = solverNode;
    solverInfo_ = olasInfo->Get("pardiso");

    // Initialise attributes
    firstCall_ = true;
    msgLvl_ = 0;

    // Resize data arrays for Pardiso
    pt_.Resize(64);
    pt_.Init(NULL);
    iparm_.Resize(64);
    iparm_.Init(0);
    dparm_.Resize(64);
    dparm_.Init(0.0);

    // Set default solver type to direct sparse solver
    PtrParamNode sNode = xml_->Get("pardiso", ParamNode::INSERT);
    std::string solverType = "direct";
    sNode->GetValue("type", solverType, ParamNode::INSERT);
      
    if(solverType == "direct") {
      mSolver_ = 0;
    } else {
      mSolver_ = 1;
    }

    sNode = xml_->Get("stoppingRule", ParamNode::INSERT);
    std::string sRule = "relNormRes0";
    sNode->GetValue("type", sRule, ParamNode::INSERT);

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
  PardisoSolver<T>::~PardisoSolver() {

    // PARDISO - Last Phase: Cleaning up the parameters
    if ( firstCall_ == false ) {

      int errorFlag = 0;
      int phase = -1;

#if PARDISO_API_VER == 4
      F77_FUNC(pardiso) ( &pt_[0], &maxfct_, &mnum_, &mType_, &phase,
                          &probDim_, theMatrix_, rowPtr_, colPtr_,
                          idPerm_, &nrhs_, &iparm_[0], &msgLvl_, &zeroDBL_,
                          &zeroDBL_, &errorFlag, &dparm_[0] );
#endif

#if PARDISO_API_VER == 3
      F77_FUNC(pardiso) ( &pt_[0], &maxfct_, &mnum_, &mType_, &phase,
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
      
      LOG_TRACE(pardisoSolver) << "Contents of pardiso-ml.out:"
                               << std::endl
                               << pardisoMLOut.rdbuf() << std::endl;
      LOG_TRACE(pardisoSolver) << " -------------------------------------------------------"
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
    delete [] ( idPerm_ );
    idPerm_  = NULL;
    idPermSize_ = 0;

  }


  // *********
  //   Setup
  // *********
  template<typename T>
  void PardisoSolver<T>::Setup( BaseMatrix &sysMat, PtrParamNode analysis_step ) {

    // Flag for check Pardiso's return status
    int errorFlag = 0;

    PtrParamNode sNode;
    sNode = xml_->Get("pardiso", ParamNode::INSERT);
    
    // Determine, whether we are expected to be verbose
    LOG_TRACE(pardisoSolver) << " -----------------------------------------"
                             << "-------------------------------------";

    // ============================================================
    //  Determine which of the two steps: symbolical and numerical
    //  factorisation must be performed
    // ============================================================
    // If only the values of the matrix entries changed, we
    // can keep the re-ordering and only perform a numerical
    // factorisation
    bool facSymbolic = false;
    bool facNumeric = true;
//    bool facNumeric = false;

    // No factorisation available, so perform both steps
    if ( firstCall_ == true ) {
      facSymbolic = true;
      facNumeric  = true;
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
    bool symPard = false;
    if ( stype == BaseMatrix::SPARSE_SYM ) {
      symPard = true;
    }

    bool defPard = false;
    sNode->GetValue("posDef", defPard, ParamNode::INSERT);

    bool herPard = false;
    sNode->GetValue("hermitean", herPard, ParamNode::INSERT);

    bool strPard = false;
    sNode->GetValue("symStruct", strPard, ParamNode::INSERT);

    mType_ = 0;
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
    else {
      LOG_TRACE(pardisoSolver) << " Classified matrix as mType = " << mType_;
    }

    if(mSolver_ == 1) {
      LOG_TRACE(pardisoSolver) << " Using iterative solver";
      switch(mType_) {
        case 11:
        case 13:
          EXCEPTION( "PardisoSolver: The iterative solver just supports symmetric matrices!" );
          break;
      }
    } else {
      LOG_TRACE(pardisoSolver) << " Using direct solver";
    }

    // Set default input values for dparm_;
    dparm_[0] = 300;   // Maximum number of Krylov-subspace iterations
    sNode->GetValue("maxIter", dparm_[0], ParamNode::INSERT);

    dparm_[1] = 1e-6;  // Relative residual reduction
    sNode->GetValue("tol", dparm_[1], ParamNode::INSERT);

    dparm_[2] = 1e-6;  // Coarse Grid Matrix Dimension.
    sNode->GetValue("coarseGridDim", dparm_[2], ParamNode::INSERT);

    dparm_[3] = 10;    // Maximum Number of Grid Levels.
    sNode->GetValue("maxNumGridLevels", dparm_[3], ParamNode::INSERT);

    dparm_[4] = 1e-2;  // Dropping value for the incomplete factor.
    sNode->GetValue("incompFacDropVal", dparm_[4], ParamNode::INSERT);

    dparm_[5] = 5e-5;  // Dropping value for the schurcomplement.
    sNode->GetValue("schurcompDropVal", dparm_[5], ParamNode::INSERT);

    dparm_[6] = 10;    // Maximum number of ﬁll-in in each column in the factor.
    sNode->GetValue("maxNumFillIn", dparm_[6], ParamNode::INSERT);

    dparm_[7] = 500;   // Bound for the inverse of the incomplete factor L.
    sNode->GetValue("invBoundIncompFac", dparm_[0], ParamNode::INSERT);

    dparm_[8] = 25;    // Maximum number of non-improvement steps in Krylov-Subspace method
    sNode->GetValue("maxNumStagnationSteps", dparm_[0], ParamNode::INSERT);

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
      LOG_TRACE(pardisoSolver) << " Calling pardisoinit";

#if PARDISO_API_VER == 4
      Integer error = 0;
      F77_FUNC(pardisoinit) ( &pt_[0],  &mType_, &mSolver_, &iparm_[0], &dparm_[0], &error);

      if (error != 0) 
        EXCEPTION(GetErrorString(error));
#endif

#if PARDISO_API_VER == 3
      F77_FUNC(pardisoinit) ( &pt_[0], &mType_, &iparm_[0] );
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
    iparm_[2] = 1;
    //#endif

    // by default not transposed (IPARM(12))
    std::string transposed = sNode->Has("Transposed") ? sNode->Get("Transposed")->As<std::string>() :  "non_transposed";
    iparm_[11] = 0;
    if(sNode->Has("Transposed")) {
      std::string par = sNode->Get("Transposed")->As<std::string>();
      if(par == "conjugate_transposed")
        iparm_[11] = 1;
      if(par == "transposed")
        iparm_[11] = 2;
    }
    LOG_DBG(pardisoSolver) << "transposed setting IPAR(12) = " << iparm_[11] << " with 0 for non-transposed";


    // Determine the re-ordering strategy: We can either fo nested dissection
    // or minimum degree re-ordering or no re-ordering at all (i.e. we use
    // the initial ordering of the linear system, which might already have been
    // re-ordered via the graph)
    BaseOrdering::ReorderingType ordering = BaseOrdering::NESTED_DISSECTION;
    sNode = xml_->Get("pardiso", ParamNode::INSERT);
    std::string orderStr = "nestedDissection";
    sNode->GetValue("ordering", orderStr, ParamNode::INSERT);
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

      EXCEPTION( "Re-ordering of type '" << tmp  << "' is not available with the PardisoSolver" );
    }

    if(!mSolver_) {
      if ( facSymbolic == true ) {
        if ( ordering != BaseOrdering::NOREORDERING ) {
          std::string tmp;
          tmp = BaseOrdering::reorderingType.ToString( ordering );

          LOG_TRACE(pardisoSolver) << " Analyse phase will determine a '" << tmp << "' re-ordering";
        }
        else {
          LOG_TRACE(pardisoSolver) << " Factorisation uses original matrix ordering";
        }
      }
    }
    
    // Do we need to determine MFLOPs for the LU factorisation
    bool stats = false;
      sNode->GetValue("stats", stats, ParamNode::INSERT);
    
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
    if (sNode->Has("outOfCore"))
    {
      if(std::string(CFS_PARDISO) != "MKL") 
      {
        WARN( "The value of outOfCore has no effect for " << CFS_PARDISO << ".\n"
              << "Switch to MKL if want to use out-of-core memory." );
      }
      sNode->GetValue("outOfCore", iparm_[59], ParamNode::INSERT);
    }

    // Switch to iterative solver
    if(mSolver_) {
#if PARDISO_API_VER == 3
      EXCEPTION("This CFS++ executable has been linked to a PARDISO 3.x library.\n"
                << "PARDISO implements iterative solvers since version 4.0.\n"
                << "To get iterative solvers you should switch CFS_PARDISO to SCHENK.");
#endif      
      
      iparm_[31] = 1;
    }

    /* The facotrisation algorithm in pardiso rows/columns which have the about
     * the same magnitude and structure into one supernode. This may cause
     * numerical errors due to pivotisation, the iterative refine step remedies
     * this. (see O.Schenk et al "Solving unsymmetric sparse systems of linear
     * equations with PARDISO")
     */
    if (sNode->Has("IterRefineSteps"))
    {
      sNode->GetValue("IterRefineSteps", iparm_[7], ParamNode::INSERT);
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
    for (int i = 0; i <= probDim_; ++i )
      rowPtr_[i] += 1;

    for (UInt i=0; i< nnz_; i++ )
       colPtr_[i] += 1;

    // ========================
    //  Symbolic Factorisation
    // ========================
    if ( facSymbolic == true ) {
      tSymfact_.Start();
      // log report
      LOG_TRACE(pardisoSolver) << " Performing analyse phase (symbolic factorisation)"
                               << " ... ";

      // only analyse
      int phase = 11;

      // let pardiso go for it
#if PARDISO_API_VER == 4
      F77_FUNC(pardiso) (&pt_[0], &maxfct_, &mnum_, &mType_, &phase,
                         &probDim_, theMatrix_, rowPtr_, colPtr_,
                         idPerm_, &nrhs_, &iparm_[0], &msgLvl_, &zeroDBL_,
                         &zeroDBL_, &errorFlag, &dparm_[0] );
#endif

#if PARDISO_API_VER == 3
      F77_FUNC(pardiso) (&pt_[0], &maxfct_, &mnum_, &mType_, &phase,
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
        LOG_TRACE(pardisoSolver) << "done";
      }

      tSymfact_.Stop();
      solverInfo_->Get("symbfact/timer/cpu")->SetValue(tSymfact_.GetCPUTime());
      solverInfo_->Get("symbfact/timer/wall")->SetValue(tSymfact_.GetWallTime());
      solverInfo_->Get("symbfact/timer/calls")->SetValue(tSymfact_.GetCalls());
    }

    // =========================
    //  Numerical Factorisation
    // =========================
    if ( facNumeric == true ) {
      tNumfact_.Start();
      // log report
      LOG_TRACE(pardisoSolver) << " Performing factorise phase (numerical "
                               << "factorisation) ... ";

      // only factorise (numerical)
      int phase = 22;

//      std::cout << "pardiso(" << "x," << maxfct_ << "," << mnum_ << "," << mType_ << "," <<
//                phase << "," << probDim_ << "," << theMatrix_ << "," << rowPtr_ << "," <<  colPtr_ << "," <<
//                                         idPerm_ << "," << nrhs_ << "," <<  iparm_[0] << "," <<  msgLvl_ << "," <<  zeroDBL_ << "," <<
//                                         zeroDBL_ << "," << errorFlag << "," <<  dparm_[0] << ")" << "\n";

      // let pardiso go for it
#if PARDISO_API_VER == 4
      F77_FUNC(pardiso) (&pt_[0], &maxfct_, &mnum_, &mType_, &phase,
                         &probDim_, theMatrix_, rowPtr_, colPtr_,
                         idPerm_, &nrhs_, &iparm_[0], &msgLvl_, &zeroDBL_,
                         &zeroDBL_, &errorFlag, &dparm_[0] );
#endif

#if PARDISO_API_VER == 3
      F77_FUNC(pardiso) (&pt_[0], &maxfct_, &mnum_, &mType_, &phase,
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
        LOG_TRACE(pardisoSolver) << "done";
        LOG_TRACE(pardisoSolver) << "Memory consumption during numerical factorization and solution: " << iparm_[16] << "kBytes";
      }

      tNumfact_.Stop();
      solverInfo_->Get("numfact/timer/cpu")->SetValue(tNumfact_.GetCPUTime());
      solverInfo_->Get("numfact/timer/wall")->SetValue(tNumfact_.GetWallTime());
      solverInfo_->Get("numfact/timer/calls")->SetValue(tNumfact_.GetCalls());
    }

    solverInfo_->Get("symbfact/memory/peak")->SetValue(iparm_[14]);
    solverInfo_->Get("symbfact/memory/permanent")->SetValue(iparm_[15]);
    solverInfo_->Get("numfact/memory/peak")->SetValue(iparm_[16]);

    // Now we were called once, and a factorisation is available
    firstCall_ = false;

    // now we undo our increment, since on our side the frist col and row
    // has an value of zero!!
    for (int i = 0; i <= probDim_; ++i )
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
                                const BaseVector &rhs, BaseVector &sol, PtrParamNode analysis_step ) {

    LOG_TRACE(pardisoSolver) << " -----------------------------------------"
                             << "-------------------------------------";
    LOG_TRACE(pardisoSolver) << " Solving linear system ...";

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

    // we have to icrement the entries of the col- and row-position arrays
    // by one, so that the first col and first row start with index 1 (and
    // not with zero) to be consistent with fortran
    // at the end of the method we will undo it!!
    for (int i = 0; i <= probDim_; ++i )
      rowPtr_[i] += 1;
    for (UInt i=0; i< nnz_; i++ )
       colPtr_[i] += 1;

#if PARDISO_API_VER == 4
    F77_FUNC(pardiso) (&pt_[0], &maxfct_, &mnum_, &mType_, &phase,
                       &probDim_, theMatrix_, rowPtr_, colPtr_,
                       idPerm_, &nrhs_, &iparm_[0], &msgLvl_, theRHS,
                       theSol, &errorFlag, &dparm_[0] );
#endif

#if PARDISO_API_VER == 3
    F77_FUNC(pardiso) (&pt_[0], &maxfct_, &mnum_, &mType_, &phase,
                       &probDim_, theMatrix_, rowPtr_, colPtr_,
                       idPerm_, &nrhs_, &iparm_[0], &msgLvl_, theRHS,
                       theSol, &errorFlag );
#endif

    // now we undo our increment, since on our side the first col and row
    // has an value of zero!!
    for (int i = 0; i <= probDim_; ++i )
      rowPtr_[i] -= 1;
    for (UInt i=0; i< nnz_; i++ )
      colPtr_[i] -= 1;

    // Check return status
    if ( errorFlag != NO_ERROR ) {
      EXCEPTION( "Error occured during solution of linear system:\n"
                 << GetErrorString(errorFlag) );
    }
    else {
      LOG_TRACE(pardisoSolver) << "done";
    }

    // Finish log report
    if(!mSolver_) {
      LOG_TRACE(pardisoSolver) << " number of iterative refinement steps: " << iparm_[6];
      LOG_TRACE(pardisoSolver) << " number of perturbed pivots: " << iparm_[13];
      LOG_TRACE(pardisoSolver) << " number of positive eigenvalues: " << iparm_[21];
      LOG_TRACE(pardisoSolver) << " number of negative eigenvalues: " << iparm_[22];
    } else {
      LOG_TRACE(pardisoSolver) << " Maximum number of iterations: " << dparm_[0];
      LOG_TRACE(pardisoSolver) << " Number of iterations after solve step: " << dparm_[34];
      LOG_TRACE(pardisoSolver) << " Relative Residual Reduction: " << dparm_[1];
      LOG_TRACE(pardisoSolver) << " Relative residual after Krylov-Subspace convergence: " << dparm_[33];
    }
    LOG_TRACE(pardisoSolver) << " -------------------------------------------------------"
                             << "-----------------------";

    // Create Report (no sensible things to write for direct solvers yet)
    ParamNode::ActionType at = progOpts->DoDetailedInfo() ? ParamNode::APPEND : ParamNode::DEFAULT;
    PtrParamNode out = solverInfo_->Get(ParamNode::PROCESS)->Get("solver", at);
    out->Get("numIter")->SetValue(-1);
    out->Get("finalNorm")->SetValue(-1.0);
  }

// Explicit template instantiation
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template class PardisoSolver<Double>;
  template class PardisoSolver<Complex>;
#endif

}
