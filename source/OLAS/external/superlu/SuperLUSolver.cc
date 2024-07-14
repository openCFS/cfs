// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

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

#include "def_use_superlu.hh"
#include "SuperLUSolver.hh"

namespace CoupledField {

  DEFINE_LOG(superluSolver, "olas.solvers.superlu")

  // ***********************
  //   Default Constructor
  // ***********************
  template<typename T>
  SuperLUSolver<T>::SuperLUSolver() :
    Astore(NULL),
    Ustore(NULL),
    Lstore(NULL),
    a(NULL),
    asub(NULL),
    xa(NULL),
    perm_c(NULL), /* column permutation vector */
    perm_r(NULL), /* row permutations from partial pivoting */
    etree(NULL),
    work(NULL),
    rhsb(NULL),
    rhsx(NULL),
    xact(NULL),
    R(NULL),
    C(NULL),
    ferr(NULL),
    berr(NULL)
  {
    EXCEPTION( "Default constructor of SuperLUSolver is forbidden!" );
  }


  // ***************
  //   Constructor
  // ***************
  template<typename T>
  SuperLUSolver<T>::SuperLUSolver( PtrParamNode solverNode,
                                   PtrParamNode olasInfo ) {


    // Set pointers to communication objects
    xml_ = solverNode;
    infoNode_ = olasInfo->Get("superlu",progOpts->DoDetailedInfo() ? ParamNode::APPEND : ParamNode::INSERT);

    // Set default solver type to direct sparse solver
    PtrParamNode sNode = xml_->Get("superlu", ParamNode::INSERT);
  }


  // **************
  //   Destructor
  // **************
  template<typename T>
  SuperLUSolver<T>::~SuperLUSolver() {
    SUPERLU_FREE (rhsb);
    SUPERLU_FREE (rhsx);
    SUPERLU_FREE (xact);
    SUPERLU_FREE (etree);
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    SUPERLU_FREE (R);
    SUPERLU_FREE (C);
    SUPERLU_FREE (ferr);
    SUPERLU_FREE (berr);
    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperMatrix_Store(&X);
    if ( lwork == 0 ) {
        Destroy_SuperNode_Matrix(&L);
        Destroy_CompCol_Matrix(&U);
    } else if ( lwork > 0 ) {
        SUPERLU_FREE(work);
    }
  }


  // *********
  //   Setup
  // *********
  template<typename T>
  void SuperLUSolver<T>::Setup( BaseMatrix &sysMat) {

    PtrParamNode sNode;
    sNode = xml_->Get("superlu", ParamNode::INSERT);
    
    // Determine, whether we are expected to be verbose
    LOG_TRACE(superluSolver) << "SuperLUSolver<T>::Setup firstCall=" << firstCall_;

    int *colind = NULL;
    int *rowptr = NULL;
    
    /* Defaults */
    lwork = 0;
    nrhs  = 1;
    equil = YES;	
    u     = 1.0;
    trans = NOTRANS;

    // ============================================================
    //  Determine which of the two steps: symbolical and numerical
    //  factorisation must be performed
    // ============================================================
    bool facSymbolic = false;

    // No factorisation available, so perform both steps
    if ( firstCall_ == true ) {
      facSymbolic = true;
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
      }

      // If only the values of the matrix entries changed, we
      // can keep the re-ordering and only perform a numerical
      // factorisation
      else {
        facSymbolic = false;
      }
    }

    // =====================================
    //  Check matrix entry and storage type
    // =====================================
    const StdMatrix& stdMat = dynamic_cast<const StdMatrix&>(sysMat);

    etype = stdMat.GetEntryType();
    if ( (etype != BaseMatrix::DOUBLE) && (etype != BaseMatrix::COMPLEX) ) {
      EXCEPTION( "SuperLU: Expected DOUBLE or COMPLEX entries, but got '"
                 << BaseMatrix::entryType.ToString(etype) << "'" );
    }

    isComplex_ = ( etype == BaseMatrix::COMPLEX);
    
    stype = stdMat.GetStorageType();
    if( stype != BaseMatrix::SPARSE_NONSYM ) {
      EXCEPTION( "SuperLU: Expected a sparseNonSym matrix, "
                 << "but got a '" << BaseMatrix::storageType.ToString( stype ) << "' matrix" );
    }

    // Determine problem size
    probDim_ = stdMat.GetNumRows();

    // ==================================================
    //  Get pointers to arrays containing information on
    //  (S)CRS matrix from the problem matrix object
    // ==================================================
    void *hatred;
    
    const CRS_Matrix<T>& crsMat = dynamic_cast<const CRS_Matrix<T>&>(sysMat);
    rowptr = (Integer*)(crsMat.GetRowPointer());
    colind = (Integer*)(crsMat.GetColPointer());
    hatred = (void *) crsMat.GetDataPointer();
    nnz_    = crsMat.GetNnz();
    
    // Do C-style casting to fix problem with over-loading vs. extern C
    // void *hatred = (void *) datPtr_;
    a = (Double *) hatred;

    // Do we need to determine MFLOPs for the LU factorisation
    bool stats = false;
    xml_->GetValue("stats", stats, ParamNode::INSERT);
    
    // We simultaneously treat a single right hand side
    nrhs_ = 1;

    /* Set the default values for options argument:
	options.Fact = DOFACT;
        options.Equil = YES;
    	options.ColPerm = COLAMD;
	options.DiagPivotThresh = 1.0;
    	options.Trans = NOTRANS;
    	options.IterRefine = NOREFINE;
    	options.SymmetricMode = NO;
    	options.PivotGrowth = NO;
    	options.ConditionNumber = NO;
    	options.PrintStat = YES;
    */
    set_default_options(&options);

    if ( lwork > 0 ) {
	work = SUPERLU_MALLOC(lwork);
	if ( !work ) {
	    ABORT("DLINSOLX: cannot allocate work[]");
	}
    }

    double* at = NULL;
    int* rowind = NULL;
    int* colptr = NULL;

    // ========================
    //  Symbolic Factorisation
    // ========================
    if ( facSymbolic ) {

      // log report
      LOG_TRACE(superluSolver) << " Performing analyse phase (symbolic factorisation)"
                               << " ... ";

      // First set defaults for SuperLU and generate CSC representation of input matrix
      if(isComplex_) 
      {
#if 0
        Rx.resize(2 * nnz_);

        /* get the default control parameters */
        umfpack_zi_defaults (&Control[0]) ;
        
        /* change the default print level for this demo */
        /* (otherwise, nothing will print) */
        Control[UMFPACK_PRL] = printingLevel_ ;

        status = umfpack_zi_transpose (probDim_, probDim_,
                                       Ap, Ai, Ax, NULL,
                                       (int *) NULL, (int *) NULL,
                                       &Rp[0], &Ri[0], &Rx[0], NULL, 0) ;
        if (status < 0)
        {
          umfpack_zi_report_status (&Control[0], status) ;
          EXCEPTION("umfpack_zi_transpose failed") ;
        }
#endif
      }
      else
      {
        dCompRow_to_CompCol(probDim_, probDim_, nnz_,
                            a, colind, rowptr,
                            &at, &rowind, &colptr);

        Astore = static_cast<NCformat*>(malloc(sizeof(NCformat)));
        Astore->nnz = nnz_;
        Astore->nzval = at;
        Astore->rowind = rowind;
        Astore->colptr = colptr;

        dCreate_CompCol_Matrix(&A, probDim_, probDim_, nnz_, at, asub, xa, SLU_NC, SLU_D, SLU_GE);
        A.Store = Astore;

        printf("Dimension %dx%d; # nonzeros %d\n", A.nrow, A.ncol, Astore->nnz);
      }

      // Finally compute symbolic factorization
      if(isComplex_) 
      {      
#if 0
        umfpack_zi_symbolic (probDim_, probDim_,
                             &Rp[0], &Ri[0], &Rx[0], NULL,
                             &Symbolic, &Control[0], &Info[0]) ;
        if (status < 0)
        {
          umfpack_zi_report_info (&Control[0], &Info[0]) ;
          umfpack_zi_report_status (&Control[0], status) ;
          EXCEPTION("umfpack_zi_symbolic failed") ;
        }

        if(printStats_) 
        { 
          std::cout << "========================================"
                    << "========================================" << std::endl;
          std::cout << "                            "
                    << "Symbolic Factorization " << std::endl;
          std::cout << "========================================"
                    << "========================================" << std::endl << std::endl;
          std::cout << "----------------------------------------"
                    << "----------------------------------------" << std::endl;         
          std::cout << " UMFPACK Control Block: " << std::endl;
          std::cout << "----------------------------------------"
                    << "----------------------------------------" << std::endl;         
          umfpack_zi_report_control(&Control[0]);
          std::cout << std::endl;         

          std::cout << "----------------------------------------"
                    << "----------------------------------------" << std::endl;         
          std::cout << " UMFPACK Info Block: " << std::endl;
          std::cout << "----------------------------------------"
                    << "----------------------------------------" << std::endl;         
          umfpack_zi_report_info (&Control[0], &Info[0]);
          std::cout << std::endl;         

          std::cout << "----------------------------------------"
                    << "----------------------------------------" << std::endl;         
          std::cout << " UMFPACK Symbolic Factorization Report: " << std::endl;
          std::cout << "----------------------------------------"
                    << "----------------------------------------" << std::endl;         
          status = umfpack_zi_report_symbolic (Symbolic, &Control[0]) ;
        }
#endif
      }
      else 
      {    
        if ( !(rhsb = doubleMalloc(probDim_ * 1)) ) ABORT("Malloc fails for rhsb[].");
        if ( !(rhsx = doubleMalloc(probDim_ * 1)) ) ABORT("Malloc fails for rhsx[].");
        dCreate_Dense_Matrix(&B, probDim_, 1, rhsb, probDim_, SLU_DN, SLU_D, SLU_GE);
        dCreate_Dense_Matrix(&X, probDim_, 1, rhsx, probDim_, SLU_DN, SLU_D, SLU_GE);
        xact = doubleMalloc(probDim_ * 1);
        ldx = probDim_;
        dGenXtrue(probDim_, 1, xact, ldx);
        dFillRHS(trans, 1, xact, ldx, &A, &B);
        
        if ( !(etree = intMalloc(probDim_)) ) ABORT("Malloc fails for etree[].");
        if ( !(perm_r = intMalloc(probDim_)) ) ABORT("Malloc fails for perm_r[].");
        if ( !(perm_c = intMalloc(probDim_)) ) ABORT("Malloc fails for perm_c[].");
        if ( !(R = (double *) SUPERLU_MALLOC(A.nrow * sizeof(double))) ) 
          ABORT("SUPERLU_MALLOC fails for R[].");
        if ( !(C = (double *) SUPERLU_MALLOC(A.ncol * sizeof(double))) )
          ABORT("SUPERLU_MALLOC fails for C[].");
        if ( !(ferr = (double *) SUPERLU_MALLOC(1 * sizeof(double))) )
          ABORT("SUPERLU_MALLOC fails for ferr[].");
        if ( !(berr = (double *) SUPERLU_MALLOC(1 * sizeof(double))) ) 
          ABORT("SUPERLU_MALLOC fails for berr[].");
        
        /* Initialize the statistics variables. */
        StatInit(&stat);
        
        /* ONLY PERFORM THE LU DECOMPOSITION */
        B.ncol = 0;  /* Indicate not to solve the system */
        dgssvx(&options, &A, perm_c, perm_r, etree, equed, R, C,
               &L, &U, work, lwork, &B, &X, &rpg, &rcond, ferr, berr,
               &Glu, &mem_usage, &stat, &info);
        
        printf("LU factorization: dgssvx() returns info %d\n", info);
        
        if ( info == 0 || info == probDim_+1 ) {
          
          if ( options.PivotGrowth ) printf("Recip. pivot growth = %e\n", rpg);
          if ( options.ConditionNumber )
	    printf("Recip. condition number = %e\n", rcond);
          Lstore = (SCformat *) L.Store;
          Ustore = (NCformat *) U.Store;
          printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
          printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
          printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz - probDim_);
          printf("FILL ratio = %.1f\n", (float)(Lstore->nnz + Ustore->nnz - probDim_)/nnz_);
          
          printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
                 mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
          fflush(stdout);
          
        } else {
          if ( info > 0 && lwork == -1 ) {
            printf("** Estimated memory: %d bytes\n", info - probDim_);
          }
        }

        if ( options.PrintStat ) {
          StatPrint(&stat);
        }
        StatFree(&stat);
      }
    }

    // Now we were called once, and a factorisation is available
    firstCall_ = false;
       
  }



  // *************************
  //   Solve linear system
  // *************************
  template<typename T>
  void SuperLUSolver<T>::Solve( const BaseMatrix &sysmat,
                                const BaseVector &rhs,
                                BaseVector &sol) {

    LOG_TRACE(superluSolver) << " -----------------------------------------"
                             << "-------------------------------------";
    LOG_TRACE(superluSolver) << " Solving linear system ...";


    if ( firstCall_ == true ) {
      EXCEPTION( "The matrix has not yet been factorised by SuperLU! "
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

    /* ------------------------------------------------------------
       NOW WE SOLVE THE LINEAR SYSTEM USING THE FACTORED FORM OF A.
       ------------------------------------------------------------*/
    options.Fact = FACTORED; /* Indicate the factored form of A is supplied. */
    B.ncol = 1;  /* Set the number of right-hand side */

    /* Initialize the statistics variables. */
    StatInit(&stat);

    if(isComplex_) 
    {      
    }    
    else
    {
      assert(rhsb != nullptr);
      assert(theRHS != nullptr);
      for(Integer idx=0; idx<probDim_; idx++) {
        rhsb[idx] = theRHS[idx];
      }

      dgssvx(&options, &A, perm_c, perm_r, etree, equed, R, C,
             &L, &U, work, lwork, &B, &X, &rpg, &rcond, ferr, berr,
             &Glu, &mem_usage, &stat, &info);

      printf("Triangular solve: dgssvx() returns info %d\n", info);

      if ( info == 0 || info == probDim_+1 ) {

          /* This is how you could access the solution matrix. */
          double *sol = (double*) ((DNformat*) X.Store)->nzval; 

          for(Integer idx=0; idx<probDim_; idx++) {
            theSol[idx] = sol[idx];
          }

	  if ( options.IterRefine ) {
              printf("Iterative Refinement:\n");
	      printf("%8s%8s%16s%16s\n", "rhs", "Steps", "FERR", "BERR");
	      for (i = 0; i < 1; ++i)
	        printf("%8d%8d%16e%16e\n", i+1, stat.RefineSteps, ferr[i], berr[i]);
          }
          fflush(stdout);
      } else {
        if ( info > 0 && lwork == -1 ) {
          printf("** Estimated memory: %d bytes\n", info - probDim_);
        }
      }
    }

    if ( options.PrintStat ) StatPrint(&stat);
    StatFree(&stat);
  }

  // Explicit template instantiation
  template class SuperLUSolver<Double>;
  template class SuperLUSolver<Complex>;

}
