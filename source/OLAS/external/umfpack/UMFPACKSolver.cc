#include <fstream>

#include <filesystem>
namespace fs = std::filesystem;

#include <umfpack.h>

#include "MatVec/BaseMatrix.hh"
#include "MatVec/StdMatrix.hh"
#include "MatVec/CRS_Matrix.hh"
#include "MatVec/SCRS_Matrix.hh"

#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ProgramOptions.hh"

#include "Utils/Timer.hh"

#include "UMFPACKSolver.hh"

namespace CoupledField {

  DEFINE_LOG(umfpackSolver, "olas.solvers.umfpack")

  // ***********************
  //   Default Constructor
  // ***********************
  template<typename T>
  UMFPACKSolver<T>::UMFPACKSolver() {
    EXCEPTION( "Default constructor of UMFPACKSolver is forbidden!" );
  }

  // ***************
  //   Constructor
  // ***************
  template<typename T>
  UMFPACKSolver<T>::UMFPACKSolver( PtrParamNode solverNode,
                                   PtrParamNode olasInfo ) {


    // Set pointers to communication objects
    xml_ = solverNode;
    infoNode_ = olasInfo->Get("umfpack", progOpts->DoDetailedInfo() ? ParamNode::APPEND : ParamNode::INSERT);

    // Initialise attributes
    firstCall_ = true;

    // Resize data arrays for UMFPACK
    Info.resize(UMFPACK_INFO);
    Control.resize(UMFPACK_CONTROL);

    // Set printing level for Control[UMFPACK_PRL]
    printingLevel_ = 1;
    xml_->GetValue("printinglevel", printingLevel_, ParamNode::INSERT);
    if(printingLevel_ < 0 || printingLevel_ > 6)
      printingLevel_ = 1;

    // Set printing of solver statistics
    printStats_ = false;
    xml_->GetValue("stats", printStats_, ParamNode::INSERT);

    // Set kind of ordering and pivoting strategy
    std::string strat = "auto";
    xml_->GetValue("strategy", strat, ParamNode::INSERT);
    if(strat == "auto")        { strategy_ = UMFPACK_STRATEGY_AUTO; }
    if(strat == "symmetric")   { strategy_ = UMFPACK_STRATEGY_SYMMETRIC; }
    if(strat == "unsymmetric") { strategy_ = UMFPACK_STRATEGY_UNSYMMETRIC; }
    
    std::string order = "amd";
    xml_->GetValue("ordering", order, ParamNode::INSERT);
    if(order == "cholmod") { ordering_ = UMFPACK_ORDERING_CHOLMOD; }
    if(order == "amd")     { ordering_ = UMFPACK_ORDERING_AMD; }
    if(order == "none")    { ordering_ = UMFPACK_ORDERING_NONE; }
    if(order == "metis")   { ordering_ = UMFPACK_ORDERING_METIS; }
    if(order == "best")    { ordering_ = UMFPACK_ORDERING_BEST; }

    tol_=0.1;    
    xml_->GetValue("tol", tol_, ParamNode::INSERT);

    symtol_=0.001;
    xml_->GetValue("symtol", symtol_, ParamNode::INSERT);

    irstep_=2;
    xml_->GetValue("irstep", irstep_, ParamNode::INSERT);

    // These pointers are used to hold the addresses of the internal
    // (S)CRS matrix structures. The related memory segments must not
    // be altered of deleted. Therefore the pointers are const!
    Ap = NULL;
    Ai = NULL;
    Ax = NULL;
  }


  // **************
  //   Destructor
  // **************
  template<typename T>
  UMFPACKSolver<T>::~UMFPACKSolver() {

    if(isComplex_) 
    {
      umfpack_zi_free_symbolic( &Symbolic );
      umfpack_zi_free_numeric( &Numeric );
    }
    else
    {  
      umfpack_di_free_symbolic( &Symbolic );
      umfpack_di_free_numeric( &Numeric );
    }
  }


  // *********
  //   Setup
  // *********
  template<typename T>
  void UMFPACKSolver<T>::Setup( BaseMatrix &sysMat)
  {

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
      EXCEPTION( "UMFPACK: Expected DOUBLE or COMPLEX entries, but got '"
                 << BaseMatrix::entryType.ToString(etype) << "'" );
    }

    isComplex_ = ( etype == BaseMatrix::COMPLEX);
    
    stype = stdMat.GetStorageType();
    if ( (stype != BaseMatrix::SPARSE_SYM) &&
         (stype != BaseMatrix::SPARSE_NONSYM) ) {
      EXCEPTION( "UMFPACK: Expected a sparseSym or sparseNonSym matrix, "
                 << "but got a '" << BaseMatrix::storageType.ToString( stype ) << "' matrix" );
    }

    // Determine problem size
    probDim_ = stdMat.GetNumRows();

    // ==================================================
    //  Get pointers to arrays containing information on
    //  (S)CRS matrix from the problem matrix object
    // ==================================================
    void *hatred;
    
    if ( stype == BaseMatrix::SPARSE_NONSYM ) {
      const CRS_Matrix<T>& crsMat = dynamic_cast<const CRS_Matrix<T>&>(sysMat);
      Ap = (Integer*)(crsMat.GetRowPointer());
      Ai = (Integer*)(crsMat.GetColPointer());
      hatred = (void *) crsMat.GetDataPointer();
      nnz_    = crsMat.GetNnz();
    }
    else {
      const SCRS_Matrix<T>& scrsMat = dynamic_cast<const SCRS_Matrix<T>&>(sysMat);
      Ap = (Integer*)(scrsMat.GetRowPointer());
      Ai = (Integer*)(scrsMat.GetColPointer());
      hatred = (void *) scrsMat.GetDataPointer();
      nnz_    = scrsMat.GetNumEntries();
    }

    // Do C-style casting to fix problem with over-loading vs. extern C
    // void *hatred = (void *) datPtr_;
    Ax = (Double *) hatred;

    // Do we need to determine MFLOPs for the LU factorisation
    bool stats = false;
    xml_->GetValue("stats", stats, ParamNode::INSERT);
    
    // We simultaneously treat a single right hand side
    nrhs_ = 1;

    // ========================
    //  Symbolic Factorisation
    // ========================
    if ( facSymbolic ) {

      // log report
      LOG_TRACE(umfpackSolver) << " Performing analyse phase (symbolic factorisation)"
                               << " ... ";
      Rp.resize(probDim_+1);
      Ri.resize(nnz_);

      // First set defaults for UMFPACK and generate CSC representation of input matrix
      if(isComplex_) 
      {      
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
        
#if 0
        int n = 5;
        int job = 1;
        int ipos = 1;
        int Ap[]={0,2,5,9,10,12};
        int Ai[]={0,1,0, 2,4,1,2,3, 4,2,1,4};
        double Ax[] = {2., 3., 3., -1., 4., 4., -3., 1., 2., 2., 6., 1., 
                       3., 3., 3., 3., 3., 3., 3., 3., 3., 3., 3., 3.} ;
        double Az[] = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.} ;
        double b[] = {8., 45., -3., 3., 19.} ;      

        //status = umfpack_zi_transpose (n, n, Ap, Ai, Ax, NULL, (int *) NULL, (int *) NULL, Rp, Ri, Rx, Rz, 0) ;

        for(int i=0; i<=n; i++) Ap[i]++;
        for(int i=0; i<12; i++) Ai[i]++;
        
        csrcsc_double_(&n, &job, &ipos,
                       Ax, Ai, Ap,
                       Rx, Ri, Rp);
        csrcsc_double_(&n, &job, &ipos,
                       Az, Ai, Ap,
                       Rz, Ri, Rp);

        for(int i=0; i<=n; i++) Rp[i]--;
        for(int i=0; i<12; i++) Ri[i]--;
        umfpack_zi_report_matrix (n, n, Rp, Ri, Rx, Rz, 1, &Control[0]) ;
#endif
      }
      else
      {
        Rx.resize(nnz_);

        /* get the default control parameters */
        umfpack_di_defaults (&Control[0]) ;
        
        /* change the default print level for this demo */
        /* (otherwise, nothing will print) */
        Control[UMFPACK_PRL] = printingLevel_ ;      
        
        //    printf ("\nMatrix A: ") ;
        //    (void) umfpack_di_report_matrix (probDim_, probDim_, rowPtr_, colPtr_, theMatrix_, 1, NULL) ;
        
        status = umfpack_di_transpose (probDim_, probDim_,
                                       Ap, Ai, Ax,
                                       (int *) NULL, (int *) NULL,
                                       &Rp[0], &Ri[0], &Rx[0]) ;
        if (status < 0)
        {
          umfpack_di_report_status (&Control[0], status) ;
          EXCEPTION("umfpack_di_transpose failed") ;
        }
      }

      // Determine the pivoting and re-ordering strategy
      if(strategy_ != UMFPACK_STRATEGY_AUTO)
      {
        Control [UMFPACK_STRATEGY] = strategy_;        
      }
      
      // Determine the re-ordering strategy
      if(strategy_ != UMFPACK_ORDERING_AMD)
      {
        Control [UMFPACK_ORDERING] = ordering_;        
      }

      Control[UMFPACK_PIVOT_TOLERANCE] = tol_;
      Control[UMFPACK_SYM_PIVOT_TOLERANCE] = symtol_;
      Control[UMFPACK_IRSTEP] = irstep_;
      
      // Finally compute symbolic factorization
      if(isComplex_) 
      {      
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
      }
      else 
      {    
        umfpack_di_symbolic (probDim_, probDim_,
                             &Rp[0], &Ri[0], &Rx[0],
                             &Symbolic, &Control[0], &Info[0]) ;
        if (status < 0)
        {
          umfpack_di_report_info (&Control[0], &Info[0]) ;
          umfpack_di_report_status (&Control[0], status) ;
          EXCEPTION("umfpack_di_symbolic failed") ;
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
          umfpack_di_report_control(&Control[0]);
          std::cout << std::endl;         

          std::cout << "----------------------------------------"
                    << "----------------------------------------" << std::endl;         
          std::cout << " UMFPACK Info Block: " << std::endl;
          std::cout << "----------------------------------------"
                    << "----------------------------------------" << std::endl;         
          umfpack_di_report_info (&Control[0], &Info[0]);
          std::cout << std::endl;         

          std::cout << "----------------------------------------"
                    << "----------------------------------------" << std::endl;         
          std::cout << " UMFPACK Symbolic Factorization Report: " << std::endl;
          std::cout << "----------------------------------------"
                    << "----------------------------------------" << std::endl;         
          status = umfpack_di_report_symbolic (Symbolic, &Control[0]) ;
        }
      }
    }

    // =========================
    //  Numerical Factorisation
    // =========================
    if ( facNumeric == true ) {

      // log report
      LOG_TRACE(umfpackSolver) << " Performing factorise phase (numerical "
                               << "factorisation) ... ";

      if(isComplex_) 
      {      
        // only factorise (numerical)
        status = umfpack_zi_numeric (&Rp[0], &Ri[0], &Rx[0], NULL,
                                     Symbolic, &Numeric, &Control[0], &Info[0]);
        if (status < 0)
        {
          EXCEPTION ("umfpack_zi_numeric failed with status " << status) ;
        }
        // printf ("\nNumeric factorization of C: ") ;
        // (void) umfpack_zi_report_numeric (Numeric, &Control[0]) ;

        if(printStats_) 
        {          
          std::cout << "========================================"
                    << "========================================" << std::endl;
          std::cout << "                            "
                    << "Numeric Factorization " << std::endl;
          std::cout << "========================================"
                    << "========================================" << std::endl << std::endl;
          std::cout << "----------------------------------------"
                    << "----------------------------------------" << std::endl;         
          std::cout << " UMFPACK Info Block: " << std::endl;
          std::cout << "----------------------------------------"
                    << "----------------------------------------" << std::endl;         
          umfpack_zi_report_info (&Control[0], &Info[0]);
          std::cout << std::endl;         

          std::cout << "----------------------------------------"
                    << "----------------------------------------" << std::endl;         
          std::cout << " UMFPACK Numeric Factorization Report: " << std::endl;
          std::cout << "----------------------------------------"
                    << "----------------------------------------" << std::endl;         
          status = umfpack_zi_report_numeric (Numeric, &Control[0]) ;
        }
      }
      else 
      {      
        status = umfpack_di_numeric (&Rp[0], &Ri[0], &Rx[0],
                                     Symbolic, &Numeric, &Control[0], &Info[0]);
        if (status < 0)
        {
          EXCEPTION ("umfpack_di_numeric failed") ;
        }
        // printf ("\nNumeric factorization of C: ") ;
        // (void) umfpack_di_report_numeric (Numeric, &Control[0]) ;

        if(printStats_) 
        {          
          std::cout << "========================================"
                    << "========================================" << std::endl;
          std::cout << "                            "
                    << "Numeric Factorization " << std::endl;
          std::cout << "========================================"
                    << "========================================" << std::endl << std::endl;
          std::cout << "----------------------------------------"
                    << "----------------------------------------" << std::endl;         
          std::cout << " UMFPACK Info Block: " << std::endl;
          std::cout << "----------------------------------------"
                    << "----------------------------------------" << std::endl;         
          umfpack_di_report_info (&Control[0], &Info[0]);
          std::cout << std::endl;         

          std::cout << "----------------------------------------"
                    << "----------------------------------------" << std::endl;         
          std::cout << " UMFPACK Numeric Factorization Report: " << std::endl;
          std::cout << "----------------------------------------"
                    << "----------------------------------------" << std::endl;         
          status = umfpack_di_report_numeric (Numeric, &Control[0]) ;
        }
      }

    }

    // Now we were called once, and a factorisation is available
    firstCall_ = false;
  }



  // *************************
  //   Solve linear system
  // *************************
  template<typename T>
  void UMFPACKSolver<T>::Solve( const BaseMatrix &sysmat,  const BaseVector &rhs, BaseVector &sol)
  {
    LOG_DBG(umfpackSolver) << "S: Solving linear system ...";
    
    if(firstCall_ == true)
      EXCEPTION( "The matrix has not yet been factorised by UMFPACK! " << "Call Setup() first" );

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


    if(isComplex_) 
    {      
      status = umfpack_zi_solve (UMFPACK_A, &Rp[0], &Ri[0], &Rx[0], NULL,
                                 theSol, NULL,
                                 theRHS, NULL,
                                 Numeric, &Control[0], &Info[0]) ;
      if (status < 0)
      {
        umfpack_zi_report_status (&Control[0], status) ;
        EXCEPTION ("umfpack_zi_solve failed") ;
      }
    }    
    else
    {      
      status = umfpack_di_solve (UMFPACK_A, &Rp[0], &Ri[0], &Rx[0],
                                 theSol,
                                 theRHS,
                                 Numeric, &Control[0], &Info[0]) ;
      if (status < 0)
      {
        umfpack_di_report_status (&Control[0], status) ;
        EXCEPTION ("umfpack_di_solve failed") ;
      }

    }

  }

  // Explicit template instantiation
  template class UMFPACKSolver<Double>;
  template class UMFPACKSolver<Complex>;

}
