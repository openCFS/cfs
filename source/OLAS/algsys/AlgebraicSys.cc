#include <iostream> 
#include <iomanip>
#include <fstream> 
#include <string> 

#include <def_use_metis.hh>
#include <def_use_pardiso.hh>
#include <def_use_lapack.hh>
#include <def_use_ilupack.hh>
#include <def_use_arpack.hh>

#include "OLAS/algsys/AlgebraicSys.hh"
#include "OLAS/algsys/SolStrategy.hh"
#include "MatVec/SBM_Matrix.hh"
#include "MatVec/VBR_Matrix.hh"
#include "MatVec/PatternPool.hh"

#include "OLAS/algsys/generateidbchandler.hh"
#include "OLAS/algsys/BaseIDBC_Handler.hh"
#include "OLAS/graph/GraphManager.hh"

#include "OLAS/solver/BaseSolver.hh"
#include "OLAS/solver/generatesolver.hh"
#include "OLAS/solver/generateEigensolver.hh"
#include "OLAS/solver/BaseEigenSolver.hh"

#include "OLAS/precond/generateprecond.hh"
#include "OLAS/precond/BasePrecond.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ProgramOptions.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

#include "Utils/Timer.hh"
#include "Domain/Domain.hh"
#include "Driver/BaseDriver.hh"
#include "Driver/AnalysisID.hh"

DECLARE_LOG(algSys)
DEFINE_LOG(algSys, "algSys")

namespace CoupledField {

  
 

  // ***********************
  //   Default Constructor
  // ***********************
  AlgebraicSys::AlgebraicSys(PtrParamNode param, PtrParamNode info, bool isSolutionComplex ) 
  {
    myParam_ = param;
    myInfo_ = info;
    graphManager_       = NULL;
    solver_             = NULL;
    eigenSolver_        = NULL;
    precond_            = NULL;
    
    
    numFcts_            = 0;
    numBlocks_          = 0;
    numRegFcts_         = 0;
    
    registrationFinished_ = false;
    systemCreated_        = false;
    distinctMatGraphs_    = false;
    
    rhs_                   =  NULL;
    sol_                   = NULL;
    sbmSymm_               = true;
    sharedPatternPossible_ = true;
    onlyOneMatrixBlock_    = false;
    statCond_              = false;
    isMatrixComplex_       = false;
    isSolutionComplex_     = isSolutionComplex;
    effMat_                = NULL;
    effRhs_                = NULL;
    effSol_                = NULL;
    tmpRHS_                = NULL;
    patternPool_           = NULL; 
    
    eigenValues_           = NULL;
    eigenValError_         = NULL;

    idbcHandler_           = NULL;

    assembleDirichletToSysMat_ = false;
    
    // Default is to always use a system matrix
    matrixTypes_.insert( SYSTEM );
    
    // Setup solution strategy enum
    PtrParamNode stratNode = myParam_->Get("solutionStrategy",ParamNode::INSERT);
    solStrat_ = SolStrategy::Generate(stratNode);
    
    // Set flag for insertion of penalty terms into matrix
    usingPenalty_ = solStrat_->UseDirichletPenalty();
    std::string aux = usingPenalty_ ? "penalty" : "elimination";
    myInfo_->Get("setup")->Get("idbcHandling")->SetValue(aux);
    
    // Set flag for insertion of penalty terms into matrix
    if ( usingPenalty_ ) {
      assembleDirichletToSysMat_ = true;
    }
    else {
      assembleDirichletToSysMat_ = false;
    }
    
    // create timer object
    graphTimer_ = boost::shared_ptr<Timer>(new Timer());
  }


  // **************
  //   Destructor
  // **************
  AlgebraicSys::~AlgebraicSys() {

    delete solver_;
    solver_ = NULL;

    delete precond_;
    precond_ = NULL;
    
    delete eigenSolver_;
    eigenSolver_ = NULL;
    
    delete graphManager_;
    graphManager_ = NULL;

    delete idbcHandler_;
    idbcHandler_ = NULL;
    
    delete eigenValues_; eigenValues_ = NULL;
    delete eigenValError_; eigenValError_ = NULL;

    for( UInt i = 0; i < numBlocks_; ++i )
      delete blockInfo_[i];
    blockInfo_.Clear();
    
    std::map<FEMatrixType, SBM_Matrix*>::iterator it;
    for (it = sysMat_.begin(); it != sysMat_.end(); ++it ) {
      delete it->second;
    }
    sysMat_.clear();
    
    delete rhs_;
    rhs_ = NULL;
    
    delete sol_;
    sol_ = NULL;
    
    // delete also "effective" matrices /vectors
    delete effRhs_;
    effRhs_ = NULL;
    
    delete effSol_;
    effSol_ = NULL;
    
    delete effMat_;
    effMat_ = NULL;

    if(tmpRHS_)
      delete tmpRHS_;
    tmpRHS_ = NULL;
    
    delete patternPool_;
    patternPool_ = NULL;
  }
  
  void AlgebraicSys::UpdateToSolStrategy() {
    LOG_TRACE(algSys) << "Updating parameters due to solution strategy";
    
    // switch according to type of solution strategy
    if( solStrat_->GetType() == SolStrategy::TWO_LEVEL_STRATEGY ) {
      
      // In case of a two-level strategy, we have to adjust the 
      // "effective" matrix / vector size, depending on the step.
      
      if( solStrat_->GetActSolStep() == 1
          && solStrat_->GetNumSolSteps() == 2) {
        // --------------------------
        //  Step 1: Only (1,1) block
        // --------------------------
        LOG_TRACE(algSys) << "\t=> Switching to reduced (1,1)-system";
        
        delete effMat_;
        delete effSol_;
        delete effRhs_;
        effMat_ = new SBM_Matrix( *sysMat_[SYSTEM], 1,1 ); 
        effRhs_ = new SBM_Vector( *rhs_, 1 );
        effSol_ = new SBM_Vector( *sol_, 1 );
        
        // important: also de-activate static condensation
        statCond_ = false;
      } else if( solStrat_->GetActSolStep() == 2 ||
                 solStrat_->GetNumSolSteps() == 1 ) {
        // --------------------------
        //  Step 2: Complete system
        // --------------------------
        LOG_TRACE(algSys) << "\t=> Switching to full system again";
        delete effMat_;
        delete effSol_;
        delete effRhs_;
        
        // re-status of static condensation
        statCond_ = solStrat_->UseStaticCondensation();
        
        if (statCond_) {
          effMat_ = new SBM_Matrix( *sysMat_[SYSTEM], numBlocks_-1, 
                                    numBlocks_-1 );
          effRhs_ = new SBM_Vector( *rhs_, numBlocks_-1 );
          effSol_ = new SBM_Vector( *sol_, numBlocks_-1 );
        } else {
          effMat_ = new SBM_Matrix( *sysMat_[SYSTEM], numBlocks_, 
                                    numBlocks_ ); 
          effRhs_ = new SBM_Vector( *rhs_, numBlocks_ );
          effSol_ = new SBM_Vector( *sol_, numBlocks_ );
        }
      } else {
        EXCEPTION("The two level solution strategy has only two steps");
      }
    } // if TWO_LEVEL strategy
    
  }


  void AlgebraicSys::CreateLinSys() {

    LOG_TRACE(algSys) << "Creating linear system";


    // first check, if registration is finished
    if( !registrationFinished_ ) {
      EXCEPTION("AlgebraicSys::CreateLinSys() may just be called after "
          "registration was finished using "
          "AlgebraicSys::FinishRegistration()" );
    }

    std::set<FEMatrixType>::iterator fIt;
    // ------------------------------
    //  Generation of matrix objects
    // ------------------------------

    // Check initially, if we can share the matrix patterns
    bool sharePattern = false;
    
    // In case of static condensation, we can not share
    // the pattern, as the VBR matrix does not support it
    if (statCond_ ) {
      sharedPatternPossible_ = false;
    }
    
    if( sharedPatternPossible_ ) {
      sharePattern = true;
      patternPool_ = new PatternPool();
      
      // insert default pattern IDs in map
      for( UInt iRow = 0; iRow < numBlocks_; ++iRow ) {
        for( UInt iCol = 0; iCol < numBlocks_; ++iCol ) {
          SubMatrixID id;
          id.rowInd = iRow;
          id.colInd = iCol;
          sbmPatternIds_[id] = NO_PATTERN_ID;
        }
      }
    }
    
    // Obtain some info from parameter file
    BaseMatrix::EntryType entryType = isMatrixComplex_ ? 
        BaseMatrix::COMPLEX :
        BaseMatrix::DOUBLE;
    for ( fIt = matrixTypes_.begin(); fIt != matrixTypes_.end(); fIt++ ) {
      sysMat_[*fIt] = GenerateSBM_Matrix( *fIt, entryType, sharePattern );
    }

    // Log what we will do
    PrintFeMatrixInfo();

    // --------------------------------------------
    //  Treatment of Dirichlet Boundary Conditions
    // --------------------------------------------

    // collect number of Dirichlet values
    if ( usingPenalty_) {
      idbcHandler_ = 
          GenerateIDBC_HandlerObjectPenalty( numDirichletValuesPerBlock_,
                                             entryType );
    }
    else {
      idbcHandler_ = 
          GenerateIDBC_HandlerObject( matrixTypes_, graphManager_,
                                      numDirichletValuesPerBlock_,
                                      entryType );
    }

    // ------------------------------
    //  Generation of vector objects
    // ------------------------------

    BaseMatrix::EntryType solEntryType = isSolutionComplex_ ? 
        BaseMatrix::COMPLEX :
        BaseMatrix::DOUBLE;

    // Generate empty SBM vectors
    rhs_ = dynamic_cast<SBM_Vector*>( GenerateVectorObject( *(sysMat_[SYSTEM]), solEntryType ) );

      sol_ = dynamic_cast<SBM_Vector*>( GenerateVectorObject( *(sysMat_[SYSTEM]), solEntryType ) );
      
      if ( rhs_ == NULL || sol_ == NULL ) {
        EXCEPTION( WRONG_CAST_MSG );
      }
      
      // For the moment we insert a sub-vector for each position.
      // In the case of the right-hand side we might actually be
      // more economic. How do we get the information which sub-vectors
      // are really needed, however?
      StdMatrix *stdMat = NULL;
      BaseVector *bVec = NULL;
      SingleVector *sVec = NULL;
      for ( UInt k =0; k < numBlocks_; k++ ) {

        // Get diag matrix for vector generation
        stdMat = sysMat_[SYSTEM]->GetPointer( k, k );

        if(stdMat == NULL){
          EXCEPTION("SBM-Block was not initialized");
        }
        
        // Insert sub-vector into solution
        bVec = GenerateVectorObject( *stdMat, solEntryType );
        sVec = dynamic_cast<SingleVector*>( bVec );
        sol_->SetSubVector( sVec, k );

        // Insert sub-vector into right-hand side
        bVec = GenerateVectorObject( *stdMat, solEntryType );
        sVec = dynamic_cast<SingleVector*>( bVec );
        rhs_->SetSubVector( sVec, k );
      }

      // ---------------------------------------
      //   Generate effective matrices vectors
      // ---------------------------------------

      // This depends on the status of static condensation
      if (statCond_) {
        effMat_ = new SBM_Matrix( *sysMat_[SYSTEM], numBlocks_-1, 
                                  numBlocks_-1 );
        effRhs_ = new SBM_Vector( *rhs_, numBlocks_-1 );
        effSol_ = new SBM_Vector( *sol_, numBlocks_-1 );
      } else {
        effMat_ = new SBM_Matrix( *sysMat_[SYSTEM], numBlocks_, 
                                  numBlocks_ ); 
        effRhs_ = new SBM_Vector( *rhs_, numBlocks_ );
        effSol_ = new SBM_Vector( *sol_, numBlocks_ );
      }
      
      // -----------------
      //  Memory clean-up
      // -----------------

      // At this point, hopefully, the graph object is no longer
      // required by anyone, so release pointer and delete manager
      // to free memory
      delete graphManager_;
      graphManager_ = NULL;
      
      // set flag
      systemCreated_ = true;
  }

  void AlgebraicSys::CreatePrecond() {
    
    LOG_TRACE(algSys) << "Creating preconditioner";
    
    // check, if system was already created
    if( !systemCreated_ ) {
      EXCEPTION( "Matrices were not created yet. Please call "
                 "AlgebraicSys::CreateLinSys() first!" );
    }
    
    PtrParamNode precondListNode = myParam_->Get("precondList");
    PtrParamNode infoNode = myInfo_->Get("precond");
    
    // if we have just one SBM matrix block, use directly
    // the specialized methods for StdMatrices
    std::string precondId = solStrat_->GetPrecondId();
    if( onlyOneMatrixBlock_ ) {
      precond_ = GenerateStdPrecondObject((*effMat_)(0,0), precondId,
                                          precondListNode, infoNode );
    } else {
      precond_ = GenerateSBMPrecondObject( *(effMat_), precondId,
                                           precondListNode, infoNode );
    }
    
  }

  void AlgebraicSys::CreateSolver() {
    
    LOG_TRACE(algSys) << "Creating solver";
    
    // check, if system was already created
    if( !systemCreated_ ) {
      EXCEPTION( "Matrices were not created yet. Please call "
          "AlgebraicSys::CreateLinSys() first!" );
    }

    PtrParamNode solverListNode = myParam_->Get("solverList");
    PtrParamNode infoNode = myInfo_->Get("solver");
    
    // if we have just one SBM matrix block, use directly
    // the specialized methods for StdMatrices
    if( onlyOneMatrixBlock_ ) {
      solver_ = GenerateSolverObject( (*effMat_)(0,0), solStrat_, 
                                      solverListNode, infoNode);
    } else {
      solver_ = GenerateSolverObject( *(effMat_), solStrat_,
                                      solverListNode, infoNode);
    }
     
  }

  void AlgebraicSys::CreateEigenSolver() {
    
    LOG_TRACE(algSys) << "Creating eigenvalue solver";
    // check, if system was already created
    if( !systemCreated_ ) {
      EXCEPTION( "Matrices were not created yet. Please call "
          "AlgebraicSys::CreateLinSys() first!" );
    }
    
    PtrParamNode esNode = myParam_->Get("eigenSolverList");
    PtrParamNode sNode = myParam_->Get("solverList");
    PtrParamNode pNode = myParam_->Get("precondList");

    // if we have just one SBM matrix block, use directly
    // the specialized methods for StdMatrices
    if( onlyOneMatrixBlock_ ) {
      eigenSolver_ =  
          GenerateEigenSolverObject( (*sysMat_[SYSTEM])(0,0), solStrat_, 
                                     esNode, sNode, pNode,
                                     myInfo_->Get("solve_eigen") ); 
    } else {
      EXCEPTION("Eigenvalue solver can currently only handle SBM "
          << "matrices with 1 block!");
    }
  }

  void AlgebraicSys::SetupPrecond()  {
    
    LOG_TRACE(algSys) << "Setup of preconditioner";
    // check, if system was already created
    if( !systemCreated_ ) {
      EXCEPTION( "Matrices were not created yet. Please call "
          "AlgebraicSys::CreateLinSys() first!" );
    }
    
    // start setup timer of preconditioner
    precond_->GetSetupTimer()->Start();
    
    // if we have just one SBM matrix block, use directly
    // the specialized methods for StdMatrices
    if( onlyOneMatrixBlock_ ) {
      precond_->Setup( (*effMat_)(0,0));
    } else {
      precond_->Setup( *(effMat_));
    }
    
    // stop setup timer of preconditioner
    precond_->GetSetupTimer()->Stop();

  }

  void AlgebraicSys::SetupSolver() {
    
    LOG_TRACE(algSys) << "Setup of solver";
    // check, if system was already created
    if( !systemCreated_ ) {
      EXCEPTION( "Matrices were not created yet. Please call "
          "AlgebraicSys::CreateLinSys() first!" );
    }
       
    // start setup timer of solver
    solver_->GetSetupTimer()->Start();
    
    // if we have just one SBM matrix block, use directly
    // the specialized methods for StdMatrices
    if( onlyOneMatrixBlock_ ) {
      solver_->Setup( (*effMat_)(0,0));
    } else {
      solver_->Setup( *effMat_);
    }
    
   // in any case, pass preconditioner object to solver, if created
    if( precond_ != NULL ) {
      solver_->SetPrecond( precond_ );
    }
    
    ExportLinSys(true, false, false); // setup

    // stop setup timer of solver
    solver_->GetSetupTimer()->Stop();
  }

  void AlgebraicSys::SetupEigenSolver( UInt numFreq, Double shift,
                                       bool isQuadratic, bool bloch ) {
    
    LOG_TRACE(algSys) << "Setup of eigenvalue solver";
    // check, if system was already created
    if( !systemCreated_ ) {
      EXCEPTION( "Matrices were not created yet. Please call "
          "AlgebraicSys::CreateLinSys() first!" );
    }
    
    // Currently we just can solve problems with one SBM block
    if( !onlyOneMatrixBlock_ ) {
      EXCEPTION("Eigenvalue solver can currently only handle SBM "
                << "matrices with 1 block!");
    }
    
    // Determine if a generalized or a quadratic eigenvalue
    // problem has to be solved
    // In the latter case, a damping matrix has to be present,
    // otherwise we issue a warning
    bool dampPresent = (matrixTypes_.find( DAMPING) != matrixTypes_.end());
    bool massPresent = (matrixTypes_.find( MASS) != matrixTypes_.end());

    // bloch only in one configuration
    assert(!bloch || (!isQuadratic && !dampPresent && massPresent));

    if( isQuadratic == true ) {
      if( dampPresent == false ) {
        EXCEPTION("To solve a quadratic eigenvalue problem, a damping " \
                  << "matrix has to be present!");
      }

      // Setup the quadratic eigenvalue solver
      eigenSolver_->Setup( (*sysMat_[STIFFNESS])(0,0), 
                           (*sysMat_[MASS])(0,0),
                           (*sysMat_[DAMPING])(0,0), numFreq, shift );
    } else {
      if( dampPresent == true ) {
        WARN("Although a damping matrix is present, only a generalized "
            << "eigenvalue will be solved (as you have specified)!");
      }

      
      if( massPresent == true ) {
        // Setup the eigenvalue solver for generalized EV problem
        eigenSolver_->Setup( (*sysMat_[STIFFNESS])(0,0), 
                             (*sysMat_[MASS])(0,0),
                             numFreq, shift, bloch);
      } else {
        // Setup the eigenvalue solver for standard EV problem
        eigenSolver_->Setup( (*sysMat_[STIFFNESS])(0,0), 
                             numFreq, shift );
      }
    }

    // Determine some basic properties and create vectors
    // for eigenvalues and related error bounds
    BaseMatrix::EntryType eType;
    if( isQuadratic || bloch ) {
      eType = BaseMatrix::COMPLEX;
    } else {
      eType = BaseMatrix::DOUBLE;
    }

    UInt totalSize = (*sysMat_[SYSTEM])(0,0).GetNumRows();

    if(eigenValues_ == NULL || eigenValError_ == NULL)
    {
      BaseVector *bVec = GenerateSingleVectorObject(eType, totalSize);
      BaseVector *errVec = GenerateSingleVectorObject(BaseMatrix::DOUBLE, totalSize);
      eigenValues_ = dynamic_cast<SingleVector*>(bVec);
      eigenValError_ = dynamic_cast<SingleVector*>(errVec);
    }
    ExportLinSys(true, false, false); // setup

  }


  void AlgebraicSys::AddIDBCToRHS() {
    LOG_TRACE(algSys) << "Add IDBC to RHS ";

    idbcHandler_->AddIDBCToRHS( rhs_ );
  }


  void AlgebraicSys::Solve(bool setIDBC) {
    
    LOG_TRACE(algSys) << "Solving problem";

    // ======================================================================
    //  CHECK FOR CALCULATION OF CONDITION NUMBER
    // ======================================================================
    // Check, if condition number is to be calculated
    bool calcCapa = solStrat_->CalcConditionNumber();

    if( calcCapa ) {
      Double condNumber = 0.0;
      Vector<Double> ev, err;
      PtrParamNode esNode = myParam_->Get("eigenSolverList");
      PtrParamNode sNode = myParam_->Get("solverList");
      PtrParamNode pNode = myParam_->Get("precondList");

      BaseEigenSolver * evs = GenerateEigenSolverObject( *(sysMat_[SYSTEM]), solStrat_,
                                     esNode, sNode, pNode, myInfo_->Get("solve_eigen") );
      PtrParamNode in = myInfo_->Get(ParamNode::PROCESS)->Get("conditionNumber", ParamNode::APPEND);

      try {
        evs->CalcConditionNumber( (*sysMat_[SYSTEM])(0,0), condNumber,
                                  ev, err );
        in->Get("value")->SetValue(condNumber);
        in = in->Get("extremalEigenValues");

        for( UInt i = 0; i < ev.GetSize(); i++ )  {
          PtrParamNode t = in->Get("eigenvalue", ParamNode::APPEND);
          t->Get("value")->SetValue(ev[i]);
          t->Get("tolerance")->SetValue(err[i]);
        }

      } catch (Exception& ex ) {
        WARN("Calculation of condition number for system matrix did not converge");

        in->Get("value")->SetValue(0.0);
      }

      delete evs;
      // in the end prevent-recalculation of evs by 
      //re-setting the value for CalcConditionNumber
      solStrat_->GetSetupNode()->Get("calcConditionNumber")->SetValue( boost::any(false) );
    }
    // ======================================================================

    if(domain->GetBasePDE()->GetName() == "LatticeBoltzmann") // we have to adjust size of solution vector for LBM optimization
    {
    	LOG_DBG(algSys) << "Resized solution vector for LBM optimization to " << (*sysMat_[SYSTEM])(0,0).GetNumRows();
    	sol_->GetPointer(0)->Resize((*sysMat_[SYSTEM])(0,0).GetNumRows());
    }

    myInfo_->GetRoot()->ToFile(); // write current info state

    // start timer of solver
    solver_->GetSolveTimer()->Start();
    
    // If the penalty formulation is used and we have inhomogeneous
    // Dirichlet boundary conditions, then the right-hand side is
    // "contaminated" with penalty terms
    if ( usingPenalty_ ) {
      solver_->SetUsingPenalty( false );
    }

    // Iterative solvers require an initial guess and in the penalty case
    // we should insert the Dirichlet values into it
    if ( dynamic_cast<BaseIterativeSolver*>(solver_) != NULL &&
        usingPenalty_ ) {
      idbcHandler_->SetDofsToIDBC( effSol_ );
      (*cla) << " Inserted Dirichlet values into initial guess"
          << std::endl;
    }

    // Assume that everything will go well
    PtrParamNode out = myInfo_->Get(ParamNode::PROCESS)->Get("solver");
    out->Get("solutionIsOkay")->SetValue(true);

    // Now modify the right-hand side vector.
    // Note: It is mandatory to incorporate the IDBC values to the
    // complete RHS.
    if ( setIDBC ) 
      idbcHandler_->AddIDBCToRHS( rhs_ );



    // -------------------------------------------
    //  Adjust RHS for due to static condensation
    // -------------------------------------------
    if( statCond_) {
      
      // Get inverted block. By definition, this is always
      // the last block on the diagonal
      StdMatrix & S_ii = 
          (*sysMat_[SYSTEM])(numBlocks_-1, numBlocks_-1);
      SingleVector& r_i = (*rhs_)(numBlocks_-1);
      for( UInt r = 0; r < numBlocks_-1; ++r ) {
        // S: system matrix
        // r: row index
        // i: index of inner row
        // calculate rhs_r -= S_ri * S_ii^-1 * rhs_i
        
        // we need a temporary vector
        SingleVector* tmp = CopySingleVectorObject(r_i);
        S_ii.Mult(r_i, *tmp);
        (*sysMat_[SYSTEM])(r,numBlocks_-1).MultSub(*tmp, (*rhs_)(r));

        // delete temporary vector at the end
        delete tmp;
      }
    }
    
    ExportLinSys(false, true, false); // pre_solve

    // Trigger solution
    if( onlyOneMatrixBlock_ )
      solver_->Solve( (*effMat_)(0,0), (*effRhs_)(0), (*effSol_)(0));
    else
      solver_->Solve( *effMat_, *effRhs_, *effSol_);

    // -------------------------------------------
    //  Adjust Sol for due to static condensation
    // -------------------------------------------
    if( statCond_) {
      // Re-insert unknowns related to the inner block back to the solution 
      // vector:
      // sol_i = S_ii^-1 (f_i - sum(S_ic * sol_c) )
      
      // wee need a temporary vector
      SingleVector*  tmp = CopySingleVectorObject((*rhs_)(numBlocks_-1));
      StdMatrix & S_ii = 
              (*sysMat_[SYSTEM])(numBlocks_-1, numBlocks_-1);

      if( sysMat_[SYSTEM]->IsSymmetric() ) {

        // sum up row contributions S_ic * sol_c
        // calculate not directly with S_ic but with S_ci^T
        // (S_ic does not exist in symmetric sbm matrices)
        for(UInt c = 0; c < numBlocks_ -1; ++c ) {
          StdMatrix &stdMat =(*sysMat_[SYSTEM])(c,numBlocks_-1);
          stdMat.MultTSub((*sol_)(c),*tmp);
        }
        S_ii.Mult(*tmp, (*sol_)(numBlocks_-1));
        //(*sol_)(numBlocks_-1).Init();
      } else {
        // sum up row contributions S_ic * sol_c
        for(UInt c = 0; c < numBlocks_ -1; ++c ) {
          StdMatrix &stdMat =(*sysMat_[SYSTEM])(numBlocks_-1,c);
          stdMat.MultSub((*sol_)(c),*tmp);
        }
        S_ii.Mult(*tmp, (*sol_)(numBlocks_-1));
      }
      // delete temporary vector
      delete tmp;
    }

    // Now de-modify the right-hand side vector
    if ( setIDBC ) 
      idbcHandler_->RemoveIDBCFromRHS( rhs_ );

    // Check that solution went fine, if not issue a warning
    if ( out->Get("solutionIsOkay")->As<bool>() == false ) {
      WARN("Solver reports a problem! Consult .las file for "
          << "further diagnostics!");
    }
    

    ExportLinSys(false, false, true); // post_solve

    // stop timer associated with solver
    solver_->GetSolveTimer()->Stop();
  }


  void AlgebraicSys::ExportLinSys(bool setup, bool pre_solve, bool post_solve)
  {
    assert((setup && !pre_solve && !post_solve) || (!setup && pre_solve && !post_solve) || (!setup && !pre_solve && post_solve));

    AnalysisID& id = domain->GetDriver()->GetAnalysisId();

    LOG_DBG(algSys) << "ELS setup=" << setup << " pre=" << pre_solve << " post=" << post_solve << " id=" << id.ToString();

    if(!solStrat_->GetParamNode()->Has("exportLinSys"))
      return;

    PtrParamNode els = solStrat_->GetParamNode()->Get("exportLinSys");

    std::string base = els->Has("baseName") ? els->Get("baseName")->As<std::string>() : progOpts->GetSimName();
    if(id.ToString(true) != "") // filename variant
      base += "_" + id.ToString(true);

    BaseMatrix::OutputFormat mat_format = BaseMatrix::outputFormat.Parse(els->Get("format")->As<std::string>());
    BaseMatrix::OutputFormat vec_format = BaseMatrix::outputFormat.Parse(els->Get("vecFormat")->As<std::string>());

    LOG_DBG(algSys) << "ELS: stiffness=" << (sysMat_.find(STIFFNESS) != sysMat_.end()) << " system=" << (sysMat_.find(SYSTEM) != sysMat_.end());

    if(setup)
    {
      // if(els->Get("format")->As<std::string>()  == "harwell-boeing") ???
      //  EXCEPTION("Harwell-Boeing Format not implemented for SBM-case");

      if(els->Get("precond")->As<bool>())
      {
        assert(precond_ != NULL); // seems to be Id by default !?
        LOG_DBG(algSys) << "ELS precond: " << BasePrecond::precondType.ToString(precond_->GetPrecondType());

        SBM_Matrix* copy = new SBM_Matrix(*(sysMat_[SYSTEM]));
        if(onlyOneMatrixBlock_)
          precond_->GetPrecondSysMat((*copy)(0,0));
        else
          precond_->GetPrecondSysMat(*copy);
        copy->Export(base + "_precond", mat_format);
        delete copy;
      }

      if(els->Get("system")->As<bool>() && sysMat_.find(SYSTEM) != sysMat_.end() && sysMat_[SYSTEM] != NULL)
        sysMat_[SYSTEM]->Export(base + "_sys", mat_format);

      if(els->Get("stiffness")->As<bool>() && sysMat_.find(STIFFNESS) != sysMat_.end() && sysMat_[STIFFNESS] != NULL)
        sysMat_[STIFFNESS]->Export(base + "_stiffness", mat_format);

      if(els->Get("stiffness_update")->As<bool>() && sysMat_.find(STIFFNESS_UPDATE) != sysMat_.end() && sysMat_[STIFFNESS_UPDATE] != NULL)
        sysMat_[STIFFNESS_UPDATE]->Export(base + "_stiffness_update", mat_format);

      if(els->Get("damping")->As<bool>() && sysMat_.find(DAMPING) != sysMat_.end() && sysMat_[DAMPING] != NULL)
        sysMat_[DAMPING]->Export(base + "_damping", mat_format);

      if(els->Get("damping_update")->As<bool>() && sysMat_.find(DAMPING_UPDATE) != sysMat_.end() && sysMat_[DAMPING_UPDATE] != NULL)
        sysMat_[DAMPING_UPDATE]->Export(base + "_damping_update", mat_format);

      if(els->Get("mass")->As<bool>() && sysMat_.find(MASS) != sysMat_.end() && sysMat_[MASS] != NULL)
        sysMat_[MASS]->Export(base + "_mass", mat_format);

      if(els->Get("mass_update")->As<bool>() && sysMat_.find(MASS_UPDATE) != sysMat_.end() && sysMat_[MASS_UPDATE] != NULL)
        sysMat_[MASS_UPDATE]->Export(base + "_mass_update", mat_format);

      if(els->Get("auxiliary")->As<bool>() && sysMat_.find(AUXILIARY) != sysMat_.end() && sysMat_[AUXILIARY] != NULL)
        sysMat_[AUXILIARY]->Export(base + "_aux", mat_format);
    }

    if(pre_solve)
    {
      // rhs is only in harwell-boing included
      if(els->Get("rhs")->As<bool>())
        rhs_->Export(base + "_rhs", vec_format);

      if(els->Get("initialGuess")->As<bool>())
        sol_->Export(base + "_intial_guess", vec_format);
    }

    if(post_solve)
    {
      // Export solution if desired. In the eigenvalue case the solutions are the eigenvectors
      if(els->Get("solution")->As<bool>())
      {
        if(!eigenSolver_)
          sol_->Export(base + "_sol", vec_format);
        else
        {
          // we have not one solution but all the eigenvectors
          assert(eigenValues_ != NULL);
          for(unsigned int i = 0; i < eigenValues_->GetSize(); i++)
          {
            GetEigenMode(i);
            sol_->Export(base + "_mode_" + lexical_cast<std::string>(i+1), vec_format);
          }
        }
      }
    }
  }



  void AlgebraicSys::CalcEigenFrequencies(Vector<Double>& frequencies, Vector<Double>& err)
  {
    LOG_TRACE(algSys) << "Calculating real-valued eigenfrequencies";

    // Trigger calculation of eigenvalues
    eigenSolver_->CalcEigenFrequencies( *eigenValues_, *eigenValError_ );

    // Hard coded cast for double values
    Vector<Double>& valVec = dynamic_cast<Vector<Double>&>(*eigenValues_);
    frequencies = valVec;

    Vector<Double>& errVec = dynamic_cast<Vector<Double>&>(*eigenValError_);
    err = errVec;

    ExportLinSys(false, false, true); // the modes are actually not the solution
  }

  void AlgebraicSys::CalcEigenFrequencies(Vector<Complex>& frequencies, Vector<Double>& err)
  {
    LOG_TRACE(algSys) << "Calculating complex-valued eigenfrequencies";

    // Check, if eigenvalue solver is quadratic, as only in this case
    // this method is well-defined

    if(!eigenSolver_->IsBloch() && !eigenSolver_->IsQuadratic())
      EXCEPTION("When solving a generalized eigenvalue problem, only " \
               << "real-valued results are obtained! Use the second " \
               << "CalcEigenFrequencies()-method!");

    // Trigger calculation of eigenvalues
    eigenSolver_->CalcEigenFrequencies( *eigenValues_, *eigenValError_ );

    // Hard coded cast for double values
    Vector<Complex>& valVec = dynamic_cast<Vector<Complex>&>(*eigenValues_);
    frequencies = valVec;

    Vector<Double>& errVec = dynamic_cast<Vector<Double>&>(*eigenValError_);
    err = errVec;

    ExportLinSys(false, false, true);
  }

  void AlgebraicSys::GetEigenMode( UInt numMode )  {
    
    LOG_DBG(algSys) << "GEM #" << numMode;
    if(eigenSolver_->IsQuadratic() || eigenSolver_->IsBloch()) {
       Vector<Complex>& solHelp = dynamic_cast<Vector<Complex> &>((*sol_)(0));
       eigenSolver_->GetComplexEigenMode(numMode, solHelp);
    } else {
      Vector<Complex>& solHelp = dynamic_cast<Vector<Complex> &>((*sol_)(0));
      eigenSolver_->GetEigenMode(numMode, solHelp);
    }
    LOG_DBG2(algSys) << "GEM -> " << sol_->ToString();
  }

  void AlgebraicSys::GraphSetupInit( UInt numFcts, 
                                     bool useDistinctGraphs ) {
    
    LOG_DBG(algSys) << "Setup matrix graph for " << numFcts << " functions.";
    LOG_DBG(algSys) << "Use distinct graphs:" << useDistinctGraphs << std::endl;
    
    
    // start timer for graph setup
    graphTimer_->Start();
    
    distinctMatGraphs_ = useDistinctGraphs;
    
    // feFunction specific data
    numFcts_ = numFcts;
    numEqnsPerFct_.Resize(numFcts);
    lastFreeEqnPerFct_.Resize(numFcts);
    eqnToSBMBlock_.Resize(numFcts);
    
    // store flag for applying static condensation
    statCond_ = solStrat_->UseStaticCondensation();
    myInfo_->Get("setup")->Get("staticCondensation")->SetValue(statCond_);
    
    // Note: currently static condensation does not work in conjunction
    // with direct coupled / mixed problems.
    if( numFcts > 1 && statCond_ ) {
      WARN("Static condensation is currently just implemented for "
          << "systems with one FeFunction only."
          << "However for mech-acou problems it works as the coupling boundary contains no "
          << "inner degrees of freedom.");
    }
  }

  
  // ***************
  //   ObtainFctId
  // ***************
  FeFctIdType AlgebraicSys::ObtainFctId( const std::string& fctString ) {

   LOG_TRACE(algSys) << "Obtaining FctId for fct '" << fctString << "'";

   // Check, if system was already finalized
   if( registrationFinished_ ) {
     EXCEPTION("Can not register new Functions after "
               "registration is completed" );
   }
   
    // Check if Fct is already registered
    std::map<FeFctIdType,std::string>::iterator it;
    for ( it = fctNames_.begin(); it != fctNames_.end(); it++ ) {
      if ( (*it).second == fctString ) {
        EXCEPTION("A FeFunction with name '" << fctString 
                 << "' was already registered!");
      }
    }
    
    // Create Id
    FeFctIdType id = fctNames_.size();
    fctNames_[id] = fctString;

    LOG_DBG(algSys) << "FctId of '" << fctString << "' is " << id;
    matIsSymm_[id] = true;
    
    return id;
  }


  // ***************
  //   RegisterFct
  // ***************
  void AlgebraicSys::RegisterFct( const FeFctIdType fctId,
                                  UInt const numEqns,
                                  UInt const numLastFreeEqn ) {
     
    LOG_DBG(algSys) << "Registering fct"
             << "\n --> fctId        = '" << fctId << "'"
             << "\n --> numEqns        = " << numEqns
             << "\n --> numLastFreeDof = " << numLastFreeEqn
             << std::endl;

    // Check, if system was already finalized
    if( registrationFinished_ ) {
      EXCEPTION("Can not register new Functions after "
          "registration is completed" );
    }

    // Remember number of equation numbers for each function
    numEqnsPerFct_[fctId] = numEqns;
    lastFreeEqnPerFct_[fctId] = numLastFreeEqn;
    eqnToSBMBlock_[fctId].Resize( numEqns );
    eqnToSBMBlock_[fctId].Init(0);
  }
  
  
  Integer AlgebraicSys::
  DefineSBMMatrixBlock( const std::map<FeFctIdType,std::set<Integer> >& eqns,
                        bool isInnerBlock ) {
    
    LOG_TRACE(algSys) << "Defining new SBM block #" << numBlocks_;

    // Check, if system was already finalized
    if( registrationFinished_ ) {
      EXCEPTION("Can not define new SBM matrix blocks after "
                "registration is completed" );
    }
    
    // Just logging output
    if (IS_LOG_ENABLED(algSys, dbg3)) {
      LOG_DBG3(algSys) << "Mapping is as follows:";
      LOG_DBG3(algSys) << "\tfctId\teqn";
      std::map<FeFctIdType,std::set<Integer> >::const_iterator it = eqns.begin();
      
      // loop over all fctIds
      for( ; it != eqns.end(); ++it ) {
        const std::set<Integer> & fctEqns = it->second;
        // loop over all eqns
        std::set<Integer>::const_iterator setIt = fctEqns.begin();
        for( ; setIt != fctEqns.end(); ++setIt) {
          LOG_DBG3(algSys) << "\t" << it->first << "\t" << *setIt;
        }
      }
    } // if logging enabled
    
    // return value: valid sbmIndex, if block is non-empty
    UInt sbmIndex = 0;
    
    // check, if map contains any entries at all
    if (eqns.size() == 0 || eqns.begin()->second.size() == 0) {
      LOG_TRACE(algSys) << "\tBlock is empty, leaving";
      // in addition, if this block is supposed to be the static condensation block,
      // we deactivate it
      
      // Notify also solution strategy about empty block ....
      WARN("Implement dropping of empty block ...");
      
      if( isInnerBlock && statCond_) {
        statCond_ = false;
        LOG_TRACE(algSys) << "\tDeactivating static condensation";
      }
      return -1;
    }
    
    sbmIndex = numBlocks_;
    numBlocks_++;
    isDiagBlockSymm_.Push_back(true);
    //if( !isInnerBlock ) {
      numDirichletValuesPerBlock_.Push_back(0);
    //}
    blockInfo_.Push_back( new GraphManager::SBMBlockInfo() );
    
    // counters for indices  
    UInt index = 0;
    UInt numLastFreeIndex = 0;
    
    // temporary map for storing fixedEqns
    std::map<FeFctIdType,std::set<UInt> > fixedEqnsToIndex;
    
    // get current SBM-block info
    GraphManager::SBMBlockInfo & bi = *(blockInfo_[sbmIndex]);
    bi.eqnToIndex.Resize(eqns.size());
    
    // loop over all entries
    std::map<FeFctIdType,std::set<Integer> >::const_iterator it = eqns.begin();
    
    // loop over all fctIds
    for( ; it != eqns.end(); ++it ) {
      FeFctIdType fctId = it->first;
      const std::set<Integer> & fctEqns = it->second;
      
      
     // loop over all equations of one fct
      std::set<Integer>::const_iterator eqnIt = fctEqns.begin();
      for( ; eqnIt != fctEqns.end(); ++eqnIt) {
        UInt actEqn = std::abs(*eqnIt);
        
        // omit homogeneous BCs
        if( actEqn == 0)
          continue;
        
        // check, if current equation is larger than the last free one
        if( actEqn > lastFreeEqnPerFct_[fctId]) {
          fixedEqnsToIndex[fctId].insert(actEqn);
        } else {
          // free equation
          bi.eqnToIndex[fctId][actEqn] = ++index;
        }        
        // remember sbm-block of this (fctId,eqn)-combination
        eqnToSBMBlock_[fctId][actEqn-1] = sbmIndex;
        
        // also remember the block, in which this functionId occurs
        fctIdsInBlocks_[fctId].insert(sbmIndex);
      } // eqns
    } // fctIds
      
    // now we know the number of unknown and fixed entries in this block
    // and we can renumber the fixed equations to start with numLastFreeEqn+1 
    numLastFreeIndex = index;

    std::map<FeFctIdType,std::set<UInt> >::const_iterator fixedFctIt = 
        fixedEqnsToIndex.begin();

    // loop over all fctIds
    for( ; fixedFctIt != fixedEqnsToIndex.end(); ++fixedFctIt ) {
      FeFctIdType fctId = fixedFctIt->first;
      const std::set<UInt> & fixedEqns = fixedFctIt->second;

      // loop over all fixed equations of this fctId
      std::set<UInt>::const_iterator fixedEqnIt = fixedEqns.begin();
      for( ; fixedEqnIt != fixedEqns.end(); ++fixedEqnIt ) {
        bi.eqnToIndex[fctId][*fixedEqnIt] = ++index; 
      } // fixed eqns
    } // fctIds

    // Store number of Dirichlet values per SBM block
    // Note: if we perform static condensation, we have one block
    // less for the IDBC-handler
    if( sbmIndex < numDirichletValuesPerBlock_.GetSize() ) {
      numDirichletValuesPerBlock_[sbmIndex] = index - numLastFreeIndex;
    }

    // set total number of unknowns for this matrix block
    if (usingPenalty_ ) {
      // Case a): If we use penalty approach, we do not create an IDBC-graph
      //          and thus we do not have to split equation numbers
      bi.size = index;
      bi.numLastFreeIndex = index;
    } else {
      // Case b): If we use elimination approach, we 
      //          create an additional IDBC-graph and sort the indices /
      //          equations according to the numLastFreeIndex value
      bi.size = index;
      bi.numLastFreeIndex = numLastFreeIndex;
    }

    // Initially we assume, there are no sub-blocks defined on the sparse matrix
    bi.hasSubBlocks = false;

    // Print final definition of SBM-block i.e. print
    // (index) -> (fctId,eqnNr) mapping
    if (IS_LOG_ENABLED(algSys, dbg3)) {
      StdVector<std::pair<UInt,UInt> > indexToIdEqn(bi.size+1);

      // loop over all functions Ids
      for( UInt iFct = 0; iFct < bi.eqnToIndex.GetSize(); ++iFct ) {
        FeFctIdType fctId = iFct;
        boost::unordered_map<UInt, UInt> & eqnToIndex = bi.eqnToIndex[iFct];
        boost::unordered_map<UInt, UInt>::iterator eqnToIndexIt = eqnToIndex.begin(); 
        for( ; eqnToIndexIt != eqnToIndex.end(); ++eqnToIndexIt ) {
          UInt eqnNr = eqnToIndexIt->first;
          UInt index = eqnToIndexIt->second;
          indexToIdEqn[index] = std::pair<UInt,UInt>(fctId,eqnNr);
        }
      }

      LOG_DBG3(algSys) << "Final block SBM #" << sbmIndex << " has "
          << bi.size << " rows/cols and " << bi.numLastFreeIndex 
          << " as highest free entry.";
      LOG_DBG3(algSys) << "\tindex\t(fctId,eqnNr)";
      for( UInt i = 1; i <= bi.size; ++i ) {
        LOG_DBG3(algSys) << "\t" << i << "\t(" << indexToIdEqn[i].first
            << ", " << indexToIdEqn[i].second << ")";
      }
    } // if logging enabled
    
    // return newly created block
    return sbmIndex;
  }

  void AlgebraicSys::RegisterSubMatrixBlocks( UInt sbmIndex, UInt numMinorBlocks ) {
    
    LOG_TRACE(algSys) << "Registering " << numMinorBlocks 
        << " sub-matrix blocks for SBM block #" << sbmIndex;

    if( registrationFinished_ ) {
      EXCEPTION("Can not register new submatrix matrix blocks after "
          "registration is completed" );
    }

    blockInfo_[sbmIndex]->indexBlocks.Resize( numMinorBlocks );
    blockInfo_[sbmIndex]->hasSubBlocks = true;
  } 


  void AlgebraicSys::DefineSubMatrixBlocks( UInt sbmIndex, UInt blockIndex,
                                            const StdVector<FeFctIdType>& fctIds,
                                            const StdVector<Integer>& eqns ) {

    LOG_DBG2(algSys) << "Defining sub-matrix block " << blockIndex 
                    << " for SBM block #" << sbmIndex;
    
    if( registrationFinished_ ) {
      EXCEPTION("Can not register new submatrix matrix blocks after "
          "registration is completed" );
    }
    
    // check for sensible block size
    if( fctIds.GetSize() == 0 ) {
      EXCEPTION("SubMatrixBlock definition " << blockIndex << " within sbmBlock #" 
                << sbmIndex << " has 0 size!");
    }
    
    // map (fctIds,eqnNrs) to (sbmIndex, row/colIndex)
    StdVector<UInt> blockNums, indices;
    MapFctIdEqnToIndex(fctIds, eqns, blockNums, indices);
    UInt subBlockSize = fctIds.GetSize();
    
    for( UInt i = 0; i < subBlockSize; ++i ) {
      GraphManager::SBMBlockInfo & bi = *blockInfo_[blockNums[i]];
      
      if( indices[i] != 0  && indices[i] <= bi.numLastFreeIndex ) {
        bi.indexBlocks[blockIndex].Push_back(indices[i]-1);
        LOG_DBG3(algSys) << "\tsbmBlock #" << blockNums[i] 
                                           << ", block " << blockIndex
                                           << ": adding " << indices[i];
      }
      if( bi.indexBlocks.GetSize() > 0 ) {
        bi.indexBlocks[blockIndex].Trim();
      }
    }
  }


  // *******************
  //   GetFEMatrixType
  // *******************
  void AlgebraicSys::GetFEMatrixTypes( std::set<FEMatrixType> &matrixTypes ) const{

    matrixTypes = matrixTypes_;

  }

  void AlgebraicSys::SetFEMatrixType( const FEMatrixType matrixType,
                                      const bool isSymmetric,
                                      const bool isComplex,
                                      const FeFctIdType fctId1,
                                      const FeFctIdType fctId2 ) {

    LOG_TRACE(algSys) << "Setting matrix type '" << feMatrixType.ToString(matrixType)
                           << "' for fct-Ids (" << fctId1 << ", " << fctId2 << ")";
                           
    LOG_TRACE(algSys) << "\tsymmetry: " << isSymmetric;
    LOG_TRACE(algSys) << "\tcomplex values: " << isComplex;

    
    // Note: The "isSymmetric" attribute is a bit misleading, as its meaning differs,
    // depending if we we have a diagonal block (fctId1 == fctId2) or an off-diagonal
    // block.
    // For a diagonal block, it really denotes, if the (quadratic) matrix is symmetric.
    // For off-diagonal blocks (which are in general rectangular shaped and not
    // "symmetric"), it denots, if the transposed (rectangular) matrix is also set
    // and thus the overall system is symmetric again.

    // Do nothing for nothing :)
    if( matrixType != NOTYPE ) {

      // Insert FE-matrix type
      matrixTypes_.insert( matrixType );

      // Construct sub-matrix identifier
      SubMatrixID sID;
      sID.rowInd = fctId1;
      sID.colInd = fctId2;

      // Insert sub-matrix identifier into
      // corresponding FE-Matrix set
      feSubMatricesByFctId_[matrixType].insert( sID );

      // Treat symmetry type  
      if( fctId1 == fctId2 ) {

        // === DIAGONAL BLOCK ===
        
        // if we have a diagonal block, we can check, if we need
        // a symmetric sparse matrix by looking at the symmetry type. 
        // Note: The sbmSymmetry is not affected in this case.
        if( !isSymmetric) {
          this->matIsSymm_[fctId1] = false;
          LOG_TRACE(algSys) << "\t=> matrix will begic";
        }
      
      } else {
        // === OFF-DIAGONAL BLOCK ===
        
        // For an off-diagonal entry the symmetric-flag denotes, if the 
        // transposed matrix is set as well. Thus, we can determine the overall
        // symmetry of the SBM-Matrix. In case at least one integrator
        // is non-symmetric, so will be the SBM matrix.
        this->sbmSymm_ &= isSymmetric;
        LOG_TRACE(algSys) << "\t=> SBM-symmetry: " << this->sbmSymm_;

        // If matrix is symmetric
        if( isSymmetric) {
          sID.colInd = fctId1;
          sID.rowInd = fctId2;

          // Insert sub-matrix identifier into
          // corresponding FE-Matrix set
          feSubMatricesByFctId_[matrixType].insert( sID );
        }
      }

      // Determine entry type of matrices. Currently we proceed as follows:
      // If at least one matrix contains complex entries, we assume that
      // the whole system will have complex entries.
      if( isComplex )
        this->isMatrixComplex_ = true;
    }
  }

  // =======================================================================
  //   Methods for creating and assembling the matrix graph
  // =======================================================================


  // ******************
  //   GraphSetupDone
  // ******************
  void AlgebraicSys::GraphSetupDone() {

    LOG_TRACE(algSys) << "Finished setup of graph";

    std::set<FEMatrixType>::iterator fIt;
    std::set<SubMatrixID,SortSubMatrixID>::iterator sIt;

    // -----------------------------------------------------------
    // Up to now, the feSubMatrices are defined on a fctId level.
    // However we need the mapping on the block level, so we have to map
    // for each matrixType (SYSTEM, STIFFNESS etc.), how the fctIds are 
    // spread to sbmBlocks and initialize the structure 
    // feSubMatricesByBlocks_.
    // -----------------------------------------------------------

    // loop over all matrix Types in feSubMatricesByFctId
    std::map<FEMatrixType, SubMatrixSet>::const_iterator itMatByFct =
        feSubMatricesByFctId_.begin();
    for( ; itMatByFct != feSubMatricesByFctId_.end(); ++itMatByFct ) {

      const FEMatrixType & matrixType = itMatByFct->first;
      const SubMatrixSet & sbmSet = itMatByFct->second;
      SubMatrixSet::const_iterator sbmIt = sbmSet.begin();

      // loop over all blocks
      for( ; sbmIt != sbmSet.end(); ++sbmIt ) {

        const FeFctIdType rowFctId = sbmIt->rowInd;
        const FeFctIdType colFctId = sbmIt->colInd;

        // get SBM blocks, in which the current (rowFctId, colFctId)
        // occurs
        std::set<UInt>& rowBlocks = fctIdsInBlocks_[rowFctId];
        std::set<UInt>& colBlocks = fctIdsInBlocks_[colFctId];

        // insert all combinations (rowBlock,colBlock) feSubMatricesByBlock_
        std::set<UInt>::const_iterator rowIt = rowBlocks.begin();
        for( ; rowIt != rowBlocks.end(); ++rowIt ) {
          std::set<UInt>::const_iterator colIt = colBlocks.begin();
          for( ; colIt != colBlocks.end(); ++colIt ) {
            SubMatrixID sID;
            sID.rowInd = *rowIt;
            sID.colInd = *colIt;

            // Insert sub-matrix identifier into corresponding FE-Matrix set
            feSubMatricesByBlocks_[matrixType].insert( sID );
          } // loop over cols
        } // loop over rows
      } // loop over subblocks (byFctId)
    } // loop over matrixType

    // --------------------------------------------------------------------
    // Determine the set of sub-matrices that we need for the SYSTEM matrix.
    // There are two different cases:
    //
    // 1) We have other FE matrices, then we generate the set by creating
    //    the union of all other sets
    //
    // 2) There are no other sets, so we generate a sub-matrix for each
    //    sub-graph that the graph manager stores
    // --------------------------------------------------------------------
    if ( matrixTypes_.size() > 1 ) {
      for ( fIt = matrixTypes_.begin(); fIt != matrixTypes_.end(); fIt++ ) {
        if ( *fIt != SYSTEM ) {
          feSubMatricesByBlocks_[SYSTEM].insert(
              feSubMatricesByBlocks_[*fIt].begin(),
              feSubMatricesByBlocks_[*fIt].end() );
        }
      }
    }
    else {
      SubMatrixID sID;
      for ( UInt i = 0; i < numBlocks_; i++ ) {
        for ( UInt j = 0; j < numBlocks_; j++ ) {
          if ( graphManager_->SubGraphExists( i, j ) == true ) {
            sID.rowInd = i;
            sID.colInd = j;
            feSubMatricesByBlocks_[SYSTEM].insert( sID );
          }
        }
      }
    }

    // determine overall number of entries (including fixed equations)
    size_ = 0;
    for( UInt i = 0; i < numBlocks_; i++ ) {
      size_ += blockInfo_[i]->size;
    }
    StdVector< StdVector<UInt> > toBeCopied(numBlocks_);
    rowIndList1_.Set(toBeCopied);
    rowList1_.Set(toBeCopied);
    rowIndList2_.Set(toBeCopied);
    rowList2_.Set(toBeCopied);
    colIndList1_.Set(toBeCopied);
    colList1_.Set(toBeCopied);
    colIndList2_.Set(toBeCopied);
    colList2_.Set(toBeCopied);


    // Determine symmetry type of diagonal SBM-blocks.
    // Up to now, we know only the symmetry of the matrices w.r.t. the
    // FeFunctions (see matIsSymm_). 
    std::map<FeFctIdType,bool>::iterator fctIt = matIsSymm_.begin();
    for( ; fctIt != matIsSymm_.end(); ++fctIt ){

      FeFctIdType fctId = fctIt->first;
      bool isFctSymm = fctIt->second;

      // loop over all diagonal SBMBlocks, where this fctId occurs
      // and perform logical AND operation regarding symmetry (i.e.
      // if at least one function in this block is non-symmetric, so
      // will be the complete block)
      std::set<UInt>& affectedBlocks = fctIdsInBlocks_[fctId];
      std::set<UInt>::iterator blockIt = affectedBlocks.begin();
      for( ; blockIt != affectedBlocks.end(); ++blockIt) {
        isDiagBlockSymm_[*blockIt] &= isFctSymm; 
      }

    }

    // --------------------------------------------------------
    //  Perform Consistency Check Before Finalizing The System
    // --------------------------------------------------------
    // Check for:
    // - symmetry of system
    // - compatible solver type
    // - compatible preconditioner type
    // - compatible reordering type
    CheckConsistency();

    // Print information about registered functions
    PrintRegistrationInfo( );

    // Collect reordering of different matrices and assemble vector
    StdVector<BaseOrdering::ReorderingType> reorder(numBlocks_);
    for( UInt i = 0; i < numBlocks_; ++i ) {
      std::string orderString = solStrat_->GetMatrixNode(i)->
          Get("reordering")->As<std::string>();
      reorder[i] = BaseOrdering::reorderingType.Parse(orderString);
    }
    
    // Finalize graph manager setup
    graphManager_->SetupDone(reorder);
    
    // Now we have all graphs and IDBC in their re-ordered state,
    // so we have to fetch the reordering array from the GraphManager and 
    // update information in the blockInfo array for all SBM-Blocks 
    
    // Loop over all blocks
    for( UInt iBlock = 0; iBlock < numBlocks_; ++iBlock ) {
      GraphManager::SBMBlockInfo &bi = *blockInfo_[iBlock];
      UInt numFcts = bi.eqnToIndex.GetSize();
      
      // Obtain reordering vector
      StdVector<UInt> newOrder;
      graphManager_->GetReordering(iBlock, newOrder);

      // Loop over all functions
      for( UInt iFct = 0; iFct < numFcts; ++iFct ) {
        boost::unordered_map<UInt, UInt> & eqnToIndex = bi.eqnToIndex[iFct];
        boost::unordered_map<UInt, UInt>::iterator it = eqnToIndex.begin();

        // Here we have to distinguish two cases:
        // a) PENALTY: We have to reorder all equations, as also the IDBC
        //             eqns are within the normal matrix. So we do not
        //             have to maintain the splitting w.r.t. lastFreeEqn
        // b) Elimination: We are just allowed to reorder equations <
        //                 numLastFreeEqnPerFct, as this is the size of
        //                 the underlying matrix. All fixed equations 
        //                 are handled by the IDBC graph.
        if( usingPenalty_) {
          for( ; it != eqnToIndex.end(); ++it ){
            it->second = newOrder[it->second-1];
          }
        } else {
          for( ; it != eqnToIndex.end(); ++it ){
            // Loop over all (eqn)->(index) entries
            if( it->first <= lastFreeEqnPerFct_[iFct]) 
              it->second = newOrder[it->second-1];
          } // loop eqns
        } // if clause
      } // loop functions
    } // loop blocks
    
    // stop timer for graph setup
    graphTimer_->Stop();
  }


  // **********************
  //   FinishRegistration
  // **********************
  void AlgebraicSys::FinishRegistration( ) {

    LOG_DBG(algSys) << "FinishRegistration";
    
    // create new graph manager object and initialize it
    graphManager_ = new GraphManager();
    graphManager_->SetupInit( numBlocks_, distinctMatGraphs_ );

    // loop over all blocks and register them with the graph manager
    for( UInt sbmIndex = 0; sbmIndex < numBlocks_; ++sbmIndex ) {
      graphManager_->RegisterBlock( sbmIndex, blockInfo_[sbmIndex]  );
    }

    // determine, if we have a "real" SBM-system with more than 1
    // block in the sbm-matrix
    if( numBlocks_ == 1  || 
        ( numBlocks_ == 2 && solStrat_->UseStaticCondensation() ) ) {
      onlyOneMatrixBlock_ = true;
    }
    
    // set flag for registration
    registrationFinished_ = true;
  }


  // *****************
  //   SetElementPos
  // *****************
  void AlgebraicSys::SetElementPos( const FeFctIdType fctId1,
                                    const StdVector<Integer>& eqnNrs1,
                                    const FeFctIdType fctId2,
                                    const StdVector<Integer>& eqnNrs2,
                                    FEMatrixType matrixType,
                                    bool setCounterPart ) {
    
    LOG_DBG(algSys) << "Setting element position for fctIds ("
                     << fctId1 << ", " << fctId2 << ")";
    LOG_DBG2(algSys) << "matrixType: " << feMatrixType.ToString(matrixType);
    LOG_DBG2(algSys) << "counterPart: " << (setCounterPart ? "yes" : "no");
    LOG_DBG2(algSys) << "EqnVec1: " << eqnNrs1.ToString();
    LOG_DBG2(algSys) << "EqnVec2: " << eqnNrs2.ToString();

    // check, if registration was already finished
#ifndef NDEBUG
    if( !registrationFinished_ ) {
      EXCEPTION("Element connectivity can only be set after "
                "AlgebraicSys::FinishRegistration() was called" );
    }
#endif
    StdVector<UInt>& rowBlocks    = rowBlocks_.Mine();
    StdVector<UInt>& colBlocks    = colBlocks_.Mine();
    StdVector<UInt>& rowNums      = rowNums_.Mine();
    StdVector<UInt>& colNums      = colNums_.Mine();
      
    // Re-map entries from (fctId,eqnNr) -> (blockNum,index)
    MapFctIdEqnToIndex(fctId1, eqnNrs1, rowBlocks, rowNums);
    MapFctIdEqnToIndex(fctId2, eqnNrs2, colBlocks, colNums);

    // Quirk: If fctId1 == fctId2, we normally
    // need not set the counterPart. If however, they are now in
    // off diagonal matrices, we have to ensure, that originally
    // symmetric matrix structures remain symmetric also in the 
    // blockSystem
    if( fctId1 == fctId2 ) 
      setCounterPart = true;

    graphManager_->SetElementPos( rowBlocks, rowNums,
                                  colBlocks, colNums,
                                  matrixType,
                                  setCounterPart );

  }

  void AlgebraicSys::MapFctIdEqnToIndex( const FeFctIdType fctId,
                                         const StdVector<Integer>& eqns,
                                         StdVector<UInt>& blockNums,
                                         StdVector<UInt>& indices ) {
    LOG_DBG(algSys) << "Mapping fctId,eqnNr to blockNum,indices";
    
    blockNums.Resize(eqns.GetSize());
    indices.Resize(eqns.GetSize());
    
    // get hold of fct-specific map
    StdVector<UInt>& eqnToBlock = eqnToSBMBlock_[fctId];
    
    UInt numEqns = eqns.GetSize();
    for( UInt iEqn = 0; iEqn < numEqns; ++iEqn ) {
      const UInt & eqnNr = std::abs(eqns[iEqn]);
      
      // take care of homogeneous BCs
      if( eqnNr == 0) {
        blockNums[iEqn] = 0;
        indices[iEqn] = 0;
      } else {
        const UInt & blockNum = eqnToBlock[eqnNr-1];
        blockNums[iEqn] = blockNum;
        indices[iEqn] = blockInfo_[blockNum]->eqnToIndex[fctId][eqnNr];
      }
    }
    
    if( IS_LOG_ENABLED(algSys,dbg3)) {
      LOG_DBG3(algSys) << "\t(fctId,eqnNr) -> (blockNum,index)";
      for( UInt i = 0; i < numEqns; ++i ) {
        LOG_DBG3(algSys) << "\t(" << fctId << ", " << eqns[i]
                         << ") -> (" << blockNums[i] << ", "
                         << indices[i] << ")";
      }
    }
  }
  void AlgebraicSys::MapFctIdEqnToIndex( const StdVector<FeFctIdType>& fctIds,
                                         const StdVector<Integer>& eqns,
                                         StdVector<UInt>& blockNums,
                                         StdVector<UInt>& indices ) {
    LOG_DBG(algSys) << "Mapping fctId,eqnNr to blockNum,indices";

    blockNums.Resize(eqns.GetSize());
    indices.Resize(eqns.GetSize());

    UInt numEqns = eqns.GetSize();
    for( UInt iEqn = 0; iEqn < numEqns; ++iEqn ) {
      const UInt & eqnNr = std::abs(eqns[iEqn]);
      const FeFctIdType & fctId = fctIds[iEqn];
      
      // get hold of fct-specific map
      StdVector<UInt>& eqnToBlock = eqnToSBMBlock_[fctId];

      // take care of homogeneous BCs
      if( eqnNr == 0) {
        blockNums[iEqn] = 0;
        indices[iEqn] = 0;
      } else {
        const UInt & blockNum = eqnToBlock[eqnNr-1];
        blockNums[iEqn] = blockNum;
        indices[iEqn] = blockInfo_[blockNum]->eqnToIndex[fctId][eqnNr];
      }
    }

    if( IS_LOG_ENABLED(algSys,dbg3)) {
      LOG_DBG3(algSys) << "\t(fctId,eqnNr) -> (blockNum,index)";
      for( UInt i = 0; i < numEqns; ++i ) {
        LOG_DBG3(algSys) << "\t(" << fctIds[i] << ", " << eqns[i]
                         << ") -> (" << blockNums[i] << ", "
                         << indices[i] << ")";
      }
    }
  }
  
  void AlgebraicSys::MapCompleteFctIdToIndex( const FeFctIdType fctId,
                                              StdVector<UInt>& blockNums,
                                              StdVector<UInt>& indices ) {

    std::set<FeFctIdType> fctIds;
    UInt size = 0;
    if ( fctId == NO_FCT_ID ) {
      std::map<FeFctIdType, std::set<UInt> >::iterator it;
      it = fctIdsInBlocks_.begin();
      for( ; it != fctIdsInBlocks_.end(); ++it ){
        fctIds.insert(it->first);
        size += numEqnsPerFct_[it->first];
      }
    } else {
      fctIds.insert(fctId);
      size = numEqnsPerFct_[fctId];
    }

    blockNums.Resize(size);
    indices.Resize(size);
    
    // loop over all functionIds
    std::set<FeFctIdType>::const_iterator fctIt = fctIds.begin();
    for( ; fctIt != fctIds.end(); ++fctIt ) {
      
      const FeFctIdType actFctId = *fctIt;
      
      // loop over all blocks
      for( UInt iBlock = 0; iBlock < numBlocks_; ++iBlock ) {
        const boost::unordered_map<UInt, UInt> & eqnToIndexSet = 
            blockInfo_[iBlock]->eqnToIndex[actFctId];
        boost::unordered_map<UInt, UInt>::const_iterator eqnIt = eqnToIndexSet.begin();
        // loop over equations
        for( ; eqnIt != eqnToIndexSet.end(); ++eqnIt ) {
          indices[eqnIt->first-1] = eqnIt->second;
          blockNums[eqnIt->first-1] = iBlock;
        } // loop over equations
      } // loop over blocks
    } // loop over functions
  }
  
  
  void AlgebraicSys::MapFctIdEqnToIndex( const FeFctIdType fctId,
                                         const Integer eqnNr,
                                         UInt& blockNum,
                                         UInt& index ) {
    if( eqnNr == 0 ) {
      blockNum = 0;
      index = 0;
    } else {
      StdVector<UInt>& eqnToBlock = eqnToSBMBlock_[fctId];
      blockNum = eqnToBlock[std::abs(eqnNr)-1];
      index = blockInfo_[blockNum]->eqnToIndex[fctId][std::abs(eqnNr)];
    }
  }


  void AlgebraicSys::
  MapCompleteFctIdToIndex(const FeFctIdType fctId,
                          std::map<UInt, std::set<UInt> >& freeIndPerBlock,
                          std::map<UInt, std::set<UInt> >& fixedIndPerBlock,
                          bool indexZeroBased ) {
    std::set<FeFctIdType> fctIds;
    UInt size = 0;
    if ( fctId == NO_FCT_ID ) {
      std::map<FeFctIdType, std::set<UInt> >::iterator it;
      it = fctIdsInBlocks_.begin();
      for( ; it != fctIdsInBlocks_.end(); ++it ){
        fctIds.insert(it->first);
        size += numEqnsPerFct_[it->first];
      }
    } else {
      fctIds.insert(fctId);
      size = numEqnsPerFct_[fctId];
    }

    //blockNums.Resize(size);
    //indices.Resize(size);

    Integer offset = 0;
    if( indexZeroBased )
      offset = -1;
    
    
    // loop over all functionIds
    std::set<FeFctIdType>::const_iterator fctIt = fctIds.begin();
    for( ; fctIt != fctIds.end(); ++fctIt ) {

      const FeFctIdType actFctId = *fctIt;

      // loop over all blocks
      for( UInt iBlock = 0; iBlock < numBlocks_; ++iBlock ) {
        const boost::unordered_map<UInt, UInt> & eqnToIndexSet = 
            blockInfo_[iBlock]->eqnToIndex[actFctId];
        boost::unordered_map<UInt, UInt>::const_iterator eqnIt = eqnToIndexSet.begin();
        // loop over equations
        for( ; eqnIt != eqnToIndexSet.end(); ++eqnIt ) {
          const UInt index = eqnIt->second;
          if( index > blockInfo_[iBlock]->numLastFreeIndex ) {
            // fixed index
            fixedIndPerBlock[iBlock].insert( 
                index - blockInfo_[iBlock]->numLastFreeIndex + offset );
          } else {
            // free index
            freeIndPerBlock[iBlock].insert( index + offset );
          }
        } // loop over equations
      } // loop over blocks
    } // loop over functions
  }


  // all methods regarding assembly 
  void AlgebraicSys::InitMatrix( FEMatrixType matrixType,
                                 const FeFctIdType fctId ) {
    
    LOG_TRACE(algSys) << "Initializing matrix " << feMatrixType.ToString(matrixType)
                      << " for fctId " << fctId;
    
    // If matrix specified init this one
    if ( matrixType != NOTYPE )
    {
      //TODO: ok, that the matrix is zero is not uncommon in case 
      //      of transient analysis. Two ways to avoid this: 
      //      first: asselmble calls this function only with 
      //             needed matrix types.
      //      second: we check here if matrix type is contained 
      //              in matrixTypes_ map
      if(matrixTypes_.find(matrixType) != matrixTypes_.end()){
      //assert(sysMat_[matrixType] != NULL);
        sysMat_[matrixType]->Init();
        idbcHandler_->InitMatrix(matrixType);
      }
    }

    // Otherwise init all matrices
    else {
      std::set<FEMatrixType>::iterator it;
      for ( it = matrixTypes_.begin(); it != matrixTypes_.end(); it++ ) {
        sysMat_[*it]->Init();
        idbcHandler_->InitMatrix(*it);
      }
    }

    // Do we need to re-insert penalty values into matrix?
    if ( ( matrixType == NOTYPE || matrixType == SYSTEM ) &&
        usingPenalty_ ) {
      assembleDirichletToSysMat_ = true;
    }
  }

  void AlgebraicSys::InitRHS( const FeFctIdType fctId ) {
    
    LOG_TRACE(algSys) << "Initializing RHS for fctId " << fctId;
    
    if ( fctId == NO_FCT_ID ) {
      // in this case initialize complete RHS   
      rhs_->Init();
    }
    else {
      // find out affected blocks
      std::set<UInt> & blockNums = fctIdsInBlocks_[fctId];
      std::set<UInt>::iterator it = blockNums.begin();
      for( ; it != blockNums.end(); ++it ) {
        rhs_->GetPointer( *it )->Init();
      }
    }
  }
  
  void AlgebraicSys::InitRHS( const SBM_Vector& newRHS ) {
    
    LOG_TRACE(algSys) << "Initializing RHS with new vector";

    
    // ensure that the RHS vector to set consists of as many
    // sub-vectors as the RHS of the system
    if( newRHS.GetSize() != numFcts_ && domain->GetBasePDE()->GetName() != "LatticeBoltzmann") {
      EXCEPTION( "New rhs consists of " << newRHS.GetSize() << " sub-vectors, the RHS of the algebraic system of " << rhs_->GetSize() << " entries." )
    }
    
    if (domain->GetBasePDE()->GetName() == "LatticeBoltzmann"){
    	*(rhs_) = newRHS;
    	LOG_DBG(algSys) << "InitRHS: Initilized rhs with vector v=" << rhs_->GetPointer(0)->ToString(0,',');
    	return;
    }

    // loop over all feFctIDs
    for(UInt i = 0; i < numFcts_; ++i ) {

      // get all (blockId,index)-combinations for the current fctId
      StdVector<UInt> blockNums, indices;
      MapCompleteFctIdToIndex( i, blockNums, indices);
      UInt size = blockNums.GetSize();
  
      // security check: ensure that sub-vector has the same size
      // as the block indices
      if( newRHS(i).GetSize() != indices.GetSize() && (domain->GetBasePDE()->GetName() != "LatticeBoltzmann")) {
        EXCEPTION( "Number of entries of " << i << "-th sub-vector and number of indices do not match!");
      }
      
      if( newRHS.GetEntryType() == BaseMatrix::DOUBLE ) {
        Vector<Double> & nRHS = 
            dynamic_cast<Vector<Double>&>( newRHS(i) );
        for( UInt j = 0; j < size; ++j ) {
          // omit entries for Dirichlet values
          if( indices[j] <= blockInfo_[blockNums[j]]->numLastFreeIndex) {
            rhs_->GetPointer(blockNums[j])
                ->SetEntry(indices[j]-1, nRHS[j] );
          }
        }
        
      } else {
        EXCEPTION("Implement me. Dont worry: mostly C&P code");
      }
    }
  }
  
  void AlgebraicSys::InitSol( const FeFctIdType fctId ) {
    
    LOG_TRACE(algSys) << "Initializing solution of fctId " << fctId;
    
    if ( fctId == NO_FCT_ID ) {
      // in this case initialize complete RHS   
      sol_->Init();
    }
    else {
      // find out affected blocks
      std::set<UInt> & blockNums = fctIdsInBlocks_[fctId];
      std::set<UInt>::iterator it = blockNums.begin();
      for( ; it != blockNums.end(); ++it ) {
        sol_->GetPointer( *it )->Init();
      }
    }
  }
  
  void AlgebraicSys::InitSol( const SBM_Vector& newSol ) {
    
    LOG_TRACE(algSys) << "Initializing solution with new vector";
    REFACTOR;
  }
  
  
  template<typename T>
  void AlgebraicSys::SetElementMatrix( FEMatrixType matrixType, 
                                       Matrix<T>& elemMat,
                                       FeFctIdType fctId1,
                                       const StdVector<Integer>& eqnNrs1,
                                       FeFctIdType fctId2,
                                       const StdVector<Integer>& eqnNrs2,
                                       bool setCounterPart,
                                       bool noStaticCond,
                                       bool isDiagonal) {
    
    LOG_DBG(algSys) << "Setting element matrix for fctIds ("
                     << fctId1 << ", " << fctId2 << ")";
    LOG_DBG2(algSys) << "Matrix: " << feMatrixType.ToString(matrixType);
    LOG_DBG2(algSys) << "EqnVec1: (" << eqnNrs1.GetSize() << "): " << eqnNrs1.ToString();
    LOG_DBG2(algSys) << "EqnVec2: (" << eqnNrs2.GetSize() << "): " << eqnNrs2.ToString();
    LOG_DBG3(algSys) << "elemMat (" << elemMat.GetNumRows() << ", " << elemMat.GetNumCols() << "):\n " << elemMat;
    LOG_DBG3(algSys) << "noStaticCond is:\n " << noStaticCond;
    LOG_DBG2(algSys) << "isDiagonal: " << (isDiagonal ? "yes" : "no");
    
    // Security check: check if we have as many equations as numRows/Cols
    // of the matrix
    if(eqnNrs1.GetSize() != elemMat.GetNumRows()){
      EXCEPTION("dummy1 " << eqnNrs1.GetSize() << " : eMat " << elemMat.GetNumRows() )
    }
    if(eqnNrs2.GetSize() != elemMat.GetNumCols()){
      EXCEPTION("dummy2 " << eqnNrs2.GetSize() << " : eMat " << elemMat.GetNumCols() )
    }
    assert( eqnNrs1.GetSize() == elemMat.GetNumRows());
    assert( eqnNrs2.GetSize() == elemMat.GetNumCols());
    
    //obtain thread local cache lists
    UInt tNum = 0;
#ifdef USE_OPENMP
    tNum = omp_get_thread_num();
#endif
    StdVector< StdVector<UInt> >& rowIndList1  = rowIndList1_.Mine(tNum);
    StdVector< StdVector<UInt> >& rowList1     = rowList1_.Mine(tNum);
    StdVector< StdVector<UInt> >& rowIndList2  = rowIndList2_.Mine(tNum);
    StdVector< StdVector<UInt> >& rowList2     = rowList2_.Mine(tNum);
    StdVector< StdVector<UInt> >& colIndList1  = colIndList1_.Mine(tNum);
    StdVector< StdVector<UInt> >& colList1     = colList1_.Mine(tNum);
    StdVector< StdVector<UInt> >& colIndList2  = colIndList2_.Mine(tNum);
    StdVector< StdVector<UInt> >& colList2     = colList2_.Mine(tNum);
    StdVector<UInt>& rowBlocks                 = rowBlocks_.Mine(tNum);
    StdVector<UInt>& colBlocks                 = colBlocks_.Mine(tNum);
    StdVector<UInt>& rowNums                   = rowNums_.Mine(tNum);
    StdVector<UInt>& colNums                   = colNums_.Mine(tNum);

    // Re-map entries from (fctId,eqnNr) -> (blockNum,index)
    MapFctIdEqnToIndex(fctId1, eqnNrs1, rowBlocks, rowNums);
    if( isDiagonal ) {
      colBlocks = rowBlocks;
      colNums = rowNums;
    } else {
      MapFctIdEqnToIndex(fctId2, eqnNrs2, colBlocks, colNums);
    }
    // initialize empty vectors
    for(UInt i = 0; i < numBlocks_; ++i ) {
      rowIndList1[i].Clear(true);
      rowList1[i].Clear(true);
      rowIndList2[i].Clear(true);
      rowList2[i].Clear(true);
      colIndList1[i].Clear(true);
      colList1[i].Clear(true);
      colIndList2[i].Clear(true);
      colList2[i].Clear(true);
    }
    UInt numRows = rowBlocks.GetSize();
    UInt numCols = colBlocks.GetSize();
    
    // Loop over all rows
    for( UInt iRow = 0; iRow < numRows; ++iRow ) {
      // get hold of block numbers and indices
      const UInt & rowBlock = rowBlocks[iRow];
      const UInt & rowNum = rowNums[iRow];

      // Compute index of graph in graph pointer matrix
      // get hold of vertex and edgelists
      StdVector<UInt> & rList1 = rowList1[rowBlock];
      StdVector<UInt> & rIndList1 = rowIndList1[rowBlock];
      StdVector<UInt> & rList2 = rowList2[rowBlock];
      StdVector<UInt> & rIndList2 = rowIndList2[rowBlock];
      // get limits of free indices
      const UInt & lastFreeRowIndex = blockInfo_[rowBlock]->numLastFreeIndex;

      // STEP 1: Generate row index list from first connect array, dropping
      //         equation numbers for dofs fixed by (in)homogeneous Dirichlet
      //         boundary conditions and changing the sign of those fixed by
      //         constraints.
      if ( rowNum > 0 ) {
        if ( rowNum > lastFreeRowIndex ) {
          rList2.Push_back( rowNum - lastFreeRowIndex - 1 );
          rIndList2.Push_back( iRow );
        } else {
          rList1.Push_back( rowNum - 1);
          rIndList1.Push_back( iRow );
        }
      }
    }

    // Loop over all columns
    for( UInt iCol = 0; iCol < numCols; ++iCol ) {

      // get hold of block numbers and indices
      const UInt & colBlock = colBlocks[iCol];
      const UInt & colNum = colNums[iCol];

      // get hold of vertex and edgelists
      StdVector<UInt> & cList1 = colList1[colBlock];
      StdVector<UInt> & cIndList1 = colIndList1[colBlock];
      StdVector<UInt> & cList2 = colList2[colBlock];
      StdVector<UInt> & cIndList2 = colIndList2[colBlock];

      // get limits of free indices
      const UInt & lastFreeColIndex = blockInfo_[colBlock]->numLastFreeIndex;

      // STEP 2: Split the second connect array into two edge lists, one for
      //         the graph and one for the IDBCgraph (which handles the indices
      //         fixed by inhomogeneous Dirichlet boundary conditions)
      if( colNum > 0 ) {
        if ( colNum > lastFreeColIndex ) {
          cList2.Push_back( colNum - lastFreeColIndex - 1);
          cIndList2.Push_back( iCol );
        }
        else {
          cList1.Push_back( colNum - 1);
          cIndList1.Push_back( iCol );
        }
      }
    } // loop over cols
    // ======================================================================
    //  End of C&P-Code
    // ======================================================================

    SBM_Matrix * actMat = sysMat_[matrixType];
    UInt rowInd, colInd;
    
    // ======================================================================
    //  S T A T I C   C O N D E N S A T I O N
    // ======================================================================
    // Perform check, if static condensation has to be performed. If so:
    // 1) Split element matrix in contributions to different sbm blocks 
    // 2) invert matrix K_II
    // 3) loop over all entries and calculate K_RR - K_RC^T * K_II^-1 * K_RC
    // 4) store back the matrices
    
    // check for static condensation and if inner block has non-zero size
    // additionally added a switch which allows to disable static cond
    // (needed for transient case where the system matrix can be build with static condensation
    //  but the matrix parts, which are needed for the calculation of the rhs, will not be condensed)
    if( statCond_ && rowIndList1[numBlocks_-1].GetSize() &&
        (noStaticCond == false) ) {
      
      LOG_DBG(algSys) << "Performing static condensation";

      // generate logging output
      LOG_DBG3(algSys) << "\tComplete matrix is:\n" << elemMat << std::endl;
      
      
      // For static condensation, we need all matrices explicitly.
      // Loop over all blocks and extract matrices
      // --------------------------------
      //  Extract all free-free matrices
      // --------------------------------
      StdVector<Matrix<T> > matrices(numBlocks_*numBlocks_);
      for( UInt iRow = 0; iRow < numBlocks_; iRow++ ) {
        for( UInt iCol = 0; iCol < numBlocks_; iCol++ ) {
          Matrix<T>& mat = matrices[numBlocks_ * iRow + iCol];
          elemMat.GetSubMatrixByInd( mat, rowIndList1[iRow],
                                     colIndList1[iCol]);
          if( IS_LOG_ENABLED(algSys, dbg3) ){
            StdVector<UInt> ind1(rowIndList1[iRow]);
            StdVector<UInt> ind2(colIndList1[iCol]);

            LOG_DBG3(algSys) << "\tMatrix (" << iRow << ", " << iCol << "), "
                << "consisting of rowIndices: "  
                << ind1.ToString() 
                << "\n\tand colIndices: " 
                << ind2.ToString()
                << " is:\n"
                << mat << std::endl;
          }
        } // col blocks
      } // row blocks
      
//      // --------------------------------
//      //  Extract all free-idbc matrices
//      // --------------------------------
//      // Extract also Free-Dirichlet block
//      StdVector<Matrix<T> > idbcMat(numBlocks_*numBlocks_);
//      for( UInt iRow = 0; iRow < numBlocks_; iRow++ ) {
//        for( UInt iCol = 0; iCol < numBlocks_; iCol++ ) {
//          Matrix<T>& mat = idbcMat[numBlocks_ * iRow + iCol];
//          elemMat.GetSubMatrixByInd( mat, rowIndList1[iRow], 
//                                     colIndList2[iCol]);
//          LOG_DBG3(algSys) << "\tMatrix (" << iRow << ", " << iCol << ") is \n"
//              << mat << std::endl;
//        } // col blocks
//      } // row blocks

      // Invert inner block and store it back into element matrix
      matrices.Last().Invert_Lapack();
      LOG_DBG3(algSys) << "\tInverted block  (" << numBlocks_-1 << ", " 
                       << numBlocks_-1 << ") is \n"
                       << matrices.Last() << std::endl;
      elemMat.SetSubMatrixByInd( matrices.Last(),
                                 rowIndList1[numBlocks_-1],
                                 colIndList1[numBlocks_-1] );

      Matrix<T>& invMat = matrices.Last();
      Matrix<T> temp;
      
      // Loop over all blocks again, multiply matrices and re-insert afterwards
      for( UInt iRow = 0; iRow < numBlocks_-1; iRow++ ) {
        for( UInt iCol = 0; iCol < numBlocks_-1; iCol++ ) {
          // calc:  K_rc = K_rc - K_ri * K_ii^-1 * K_ic
          Matrix<T>& k_rc = matrices[numBlocks_ * iRow + iCol];
          Matrix<T>& k_ri = matrices[numBlocks_ * iRow + numBlocks_-1];
          Matrix<T>& k_ic = matrices[numBlocks_ * (numBlocks_-1) + iCol];
          temp.Resize(invMat.GetNumRows(),k_ic.GetNumCols());
          
          // Fast method: use BLAS 
          invMat.Mult_Blas(k_ic, temp, false, false, 1.0, 0.0);
          k_ri.Mult_Blas(temp, k_rc, false, false, -1.0, 1.0);
          // Alternative solution without BLAS
//            temp = invMat * k_ic;
//            k_rc -= k_ri*temp;
          elemMat.SetSubMatrixByInd( k_rc, rowIndList1[iRow],
                                     colIndList1[iCol]);
        } // col blocks
      } // row blocks
      
      
      // Note: Honestly I do not understand, why this block is NOT 
      // necessary, as we also need to orthogonalize the free-idbc block w.r.t.
      // to the inner one
      
//      // Loop over all blocks again, multiply matrices and re-insert afterwards
//      for( UInt iRow = 0; iRow < numBlocks_-1; iRow++ ) {
//        for( UInt iCol = 0; iCol < numBlocks_-1; iCol++ ) {
//          // just continue, if Dirichlet values are present
//          
//          if( colIndList2[iCol].size() == 0) continue;
//          // calc:  K_rc = K_rc - K_ri * K_ii^-1 * K_ic
//          Matrix<T>& k_rc = idbcMat[numBlocks_ * iRow + iCol];
//          Matrix<T>& k_ri = matrices[numBlocks_ * iRow + numBlocks_-1];
//          Matrix<T>& k_ic = idbcMat[numBlocks_ * (numBlocks_-1) + iCol];
//          temp.Resize(invMat.GetNumRows(),k_ic.GetNumCols());
//
//          // Fast method: use BLAS 
//          //invMat.Mult_Blas(k_ic, temp,false,false,1.0,0.0);
//          //k_ri.Mult_Blas(temp,k_rc,false,false,-1.0,1.0);
//          // Alternative solution without BLAS
//          LOG_DBG3(algSys) << "k_rc before was:\n" << k_rc;
//          temp = invMat * k_ic;
//          k_rc -= k_ri*temp;
//          LOG_DBG3(algSys) << "k_rc afterwards is:\n" << k_rc;
//          elemMat.SetSubMatrixByInd( k_rc,rowIndList1[iRow],
//                                     colIndList2[iCol]);
//        } // col blocks
//      } // row blocks

      
      LOG_DBG3(algSys) << "matrix AFTER static condensation is:\n " << elemMat;
    } // if static cond
    
    
    // ======================================================================
    
    // loop over all blocks and pass for every block the information to
    // the corresponding graph / IDBC graph
    LOG_DBG3(algSys) << "setting matrix entries";
    for( UInt sbmRow = 0; sbmRow < numBlocks_; ++sbmRow ) {
      
      // Maybe we have to switch here depending on transpose
      for( UInt sbmCol = 0; sbmCol < numBlocks_; ++sbmCol ) {
        LOG_DBG3(algSys) << "\tsetting SBM block (" << sbmRow 
                          << "," << sbmCol << ")";
        
        StdMatrix * stdMat = actMat->GetPointer(sbmRow, sbmCol);
        const StdVector<UInt> & rList1 = rowList1[sbmRow];
        const StdVector<UInt> & rList2 = rowList2[sbmRow];
        const StdVector<UInt> & cList1 = colList1[sbmCol];
        const StdVector<UInt> & cList2 = colList2[sbmCol];
        const StdVector<UInt> & rIndList1 = rowIndList1[sbmRow];
        const StdVector<UInt> & rIndList2 = rowIndList2[sbmRow];
        const StdVector<UInt> & cIndList1 = colIndList1[sbmCol];
        const StdVector<UInt> & cIndList2 = colIndList2[sbmCol];

        // Attention: This check is not really implemented in a clean way!
        if( stdMat != NULL ) {
          LOG_DBG3(algSys) << "\t1) free-free entries:";
          LOG_DBG3(algSys) << "\t\trowIndices: " << rList1.ToString();
          LOG_DBG3(algSys) << "\t\tcolIndices: " << cList1.ToString();
          LOG_DBG3(algSys) << "\t\tmat: " << stdMat->ToInfoString();
          // 2) Assemble all free <-> free entries
          // loop over all rows/col

          // Note: The following statement is experimental and not
          // thorougly tested!

          for ( UInt i = 0; i < rList1.GetSize(); i++ ) {
            rowInd = rIndList1[i];
            for ( UInt j = 0; j < cList1.GetSize(); j++ ) {
              colInd = cIndList1[j];
              stdMat->AddToMatrixEntry( rList1[i], cList1[j],
                                        elemMat[rowInd][colInd] );
            } //j
          } //i


          // 2) if sbmRow == sbmCol and transposed should be set,
          // we have to assemble the transposed by hand
          // loop over all rows/col
          if( sbmRow == sbmCol && setCounterPart ) {
            LOG_DBG3(algSys) << "\t2) free-free entries (transposed):";
            LOG_DBG3(algSys) << "\t\trowIndices: " << cList1.ToString();
            LOG_DBG3(algSys) << "\t\tcolIndices: " << rList1.ToString();
            for ( UInt i = 0; i < rList1.GetSize(); i++ ) {
              rowInd = rIndList1[i];
              for ( UInt j = 0; j < cList1.GetSize(); j++ ) {
                colInd = cIndList1[j];
                stdMat->AddToMatrixEntry( cList1[j], rList1[i],
                                          elemMat[rowInd][colInd] );
              } //j
            } //i
          } // sbmRow == sbmCol
        } // stdMat != NULL

        // 3) Assemble all free <-> fixed entries
        if( cList2.GetSize() ) {
          LOG_DBG3(algSys) << "\t3) free-fixed entries:";
          LOG_DBG3(algSys) << "\t\trowIndices: " << rList1.ToString();
          LOG_DBG3(algSys) << "\t\tcolIndices: " << cList2.ToString();

          for ( UInt i = 0; i < rList1.GetSize(); i++ ) {
            rowInd = rIndList1[i];
            for ( UInt j = 0; j < cList2.GetSize(); j++ ) {
              colInd = cIndList2[j];
              idbcHandler_->AddWeightFixedToFree( matrixType, sbmRow, sbmCol,
                                                  rList1[i], cList2[j],
                                                  elemMat[rowInd][colInd]);
            } // j
          } // i
        } // if cList2.GetSize()

        // 4) Assemble all free <-> fixed entries ( TRANSPOSED )
        if( rList2.GetSize() ) {
          if( sbmRow == sbmCol && setCounterPart == true) {
            LOG_DBG3(algSys) << "\t4) free-fixed entries (transposed):";
            LOG_DBG3(algSys) << "\t\trowIndices: " << cList1.ToString();
            LOG_DBG3(algSys) << "\t\tcolIndices: " << rList2.ToString();

            for ( UInt i = 0; i < rList2.GetSize(); i++ ) {
              rowInd = rIndList2[i];
              for ( UInt j = 0; j < cList1.GetSize(); j++ ) {
                colInd = cIndList1[j];
                idbcHandler_->AddWeightFixedToFree( matrixType, sbmCol, sbmRow,
                                                    cList1[j], rList2[i], 
                                                    elemMat[rowInd][colInd]);
              } // j
            } // i
          } // sbmCol == sbmRow
        } // rList2.GetSize()

      } //smbCol
    }// sbmRow
  }

  template<typename T>
  void AlgebraicSys::SetElementRHS( const Vector<T>& elemRHS, 
                                    const FeFctIdType fctId,
                                    StdVector<Integer>& eqnNrs ) {

    LOG_DBG(algSys) << "Setting element RHS for fctId ("<< fctId << ")";
    LOG_DBG2(algSys) << "EqnVec: " << eqnNrs.ToString();
    LOG_DBG3(algSys) << "vector is:\n " << elemRHS.ToString();
    
    // Ensure that there are as many equations as vector entries
    assert( eqnNrs.GetSize() == elemRHS.GetSize());
    
    // Re-map entries from (fctId,eqnNr) -> (blockNum,index)
    StdVector<UInt>& rowBlocks    = rowBlocks_.Mine();
    StdVector<UInt>& rowNums      = rowNums_.Mine();
    MapFctIdEqnToIndex(fctId, eqnNrs, rowBlocks, rowNums);
    
    // Now, dismantle equations
     UInt numRows = rowBlocks.GetSize();
     
     // Loop over all rows
     for( UInt iRow = 0; iRow < numRows; ++iRow ) {
       // get hold of block numbers and indices
       const UInt & rowBlock = rowBlocks[iRow];
       const UInt & rowNum = rowNums[iRow];

       // get limits of free indices
       const UInt & lastFreeRowIndex = blockInfo_[rowBlock]->numLastFreeIndex;
       
       // get vector
       SingleVector &vec = (*rhs_)(rowBlock);
       
       if ( rowNum > 0 && rowNum <= lastFreeRowIndex ) {
         if ( rowNum <= lastFreeRowIndex ) {
           vec.AddToEntry( rowNum-1, elemRHS[iRow]);
         }
       } // loop over rows
     } // loop over blocks 
  } 


  template<typename T>
  void AlgebraicSys::SetNodeRHS(T val, FeFctIdType fctId,
                                Integer eqnNr) {
    
    LOG_DBG(algSys) << "Setting node RHS of " << eqnNr << " for fct " 
                    << fctId << " to " << val;
    UInt block,idx;
    this->MapFctIdEqnToIndex(fctId,eqnNr,block,idx);
    rhs_->GetPointer(block)->AddToEntry(idx-1,val);
  }

  template<typename T>
  void AlgebraicSys::SetFncRHS(  const Vector<T>& fncRHS, FeFctIdType fctId ) {

    LOG_DBG(algSys) << "Setting Function RHS for fctId ("<< fctId << ")";

    // Re-map entries from (fctId,eqnNr) -> (blockNum,index)
    StdVector<UInt> blockNums, indices;
    MapCompleteFctIdToIndex( fctId, blockNums, indices);

    // Now, dismantle equations
     UInt numRows = blockNums.GetSize();

     // Loop over all rows
     for( UInt iRow = 0; iRow < numRows; ++iRow ) {

       // get hold of block numbers and indices
       const UInt & rowBlock = blockNums[iRow];
       const UInt & rowNum = indices[iRow];

       // get limits of free indices
       const UInt & lastFreeRowIndex = blockInfo_[rowBlock]->numLastFreeIndex;

       // get vector
       //SingleVector &vec = (*rhs_)(rowBlock);

       if ( rowNum > 0 && rowNum <= lastFreeRowIndex ) {
         rhs_->GetPointer(rowBlock)->AddToEntry(rowNum-1, fncRHS[iRow]);
       } // loop over rows
     } // loop over blocks
  }


  void AlgebraicSys::UpdateRHS(FEMatrixType matrixType, 
                               const SBM_Vector& fup,bool SysMatUpdated) {
    
    LOG_TRACE(algSys) << "Updating RHS of matrix " 
                      << feMatrixType.ToString(matrixType);

    if(matrixTypes_.find(matrixType) == matrixTypes_.end())
      return;

    // ensure that the RHS vector to set consists of as many
    // sub-vectors as the RHS of the system
    if( fup.GetSize() != numFcts_ ) {
      EXCEPTION( "New rhs consists of " << fup.GetSize()
                 << " sub-vectors, the RHS of the algebraic system of "
                 << rhs_->GetSize() << " entries." )
    }

    // loop over all feFctIDs and create a converted rhs
    if(tmpRHS_==NULL)
      tmpRHS_ = dynamic_cast<SBM_Vector*> ( GenerateVectorObject( *(sysMat_[SYSTEM]) ) );

    tmpRHS_->Init();

    for(UInt i = 0; i < numFcts_; ++i ) {

      // get all (blockId,index)-combinations for the current fctId
      StdVector<UInt> blockNums, indices;
      MapCompleteFctIdToIndex( i, blockNums, indices);
      UInt size = blockNums.GetSize();

      // security check: ensure that sub-vector has the same size
      // as the block indices
      if( fup(i).GetSize() != indices.GetSize() ) {
        EXCEPTION( "Number of entries of " << i << "-th sub-vector and number "
                   "of indices do not match!");
      }

      if( fup.GetEntryType() == BaseMatrix::DOUBLE ) {
        Vector<Double> & nRHS =
            dynamic_cast<Vector<Double>&>( fup(i) );

        for( UInt j = 0; j < size; ++j ) {
          // omit entries for Dirichlet values
          if( indices[j] <= blockInfo_[blockNums[j]]->numLastFreeIndex) {
            tmpRHS_->GetPointer(blockNums[j])
                ->AddToEntry(indices[j]-1, nRHS[j] );
          }else if(!usingPenalty_){
            idbcHandler_->AddFixedToFreeRHS(matrixType,blockNums[j],
			    		indices[j],rhs_,nRHS[j]);
          }
        }

      } else {
        EXCEPTION("Implement me. Dont worry: mostly C&P code");
      }

    }
    //now just perform multiplication
    sysMat_[matrixType]->MultAdd(*tmpRHS_,*rhs_);

  }
  
  template<typename T>
  void AlgebraicSys::AddToDiagMatrixEntry( FEMatrixType matrixType,
                                           const FeFctIdType fctId,
                                           Integer eqnNum,
                                           T val ) {
    
    LOG_DBG(algSys) << "Adding " << val << " to diagonal eq " << eqnNum
                    << " of fctId " << fctId << " for matrix " 
                    << feMatrixType.ToString(matrixType);
    REFACTOR;
  }

  template<typename T>
  void AlgebraicSys::GetMatrixEntry( FEMatrixType matrixType,
                                     const FeFctIdType rowFctId,
                                     Integer eqnNum1,
                                     const FeFctIdType colFctId,
                                     Integer eqnNum2,
                                     T & val ) {
    
    LOG_DBG(algSys) << "Getting real-valued matrix entry";
    REFACTOR;
  }

  template<typename T>
  void AlgebraicSys::SetMatrixEntry( FEMatrixType matrixId,
                                     const FeFctIdType rowFctId,
                                     Integer rowEqnNum,
                                     const FeFctIdType colFctId,
                                     Integer colEqnNum,
                                     T val,
                                     bool setCounterPart ) {
   
    LOG_DBG(algSys) << "Setting real-valued matrix entry";
    REFACTOR;
  }

  void AlgebraicSys::
  ConstructEffectiveMatrix( const FeFctIdType fctId, 
                            const std::map<FEMatrixType,Double> &matFactors ) {

    LOG_TRACE(algSys) << "Constructing effective system matrix for feFunction "
        << "with id " << fctId;
    if (IS_LOG_ENABLED(algSys, dbg)) {
      LOG_DBG(algSys) << "Factors are:";
      std::map<FEMatrixType,Double>::const_iterator it = matFactors.begin();
      for( ; it != matFactors.end(); ++it ) {
        LOG_DBG(algSys) << feMatrixType.ToString(it->first) << ": " << it->second;
      }
    }

    factorMap::const_iterator it;
    SBM_Matrix *sys = sysMat_[SYSTEM];
    
    // As one functionId can be spread over many SBM blocks, we
    // have to map the fctId to (sbmBlocks,indices)
    std::map<UInt, std::set<UInt> > freeIndPerBlock, fixedIndPerBlock;
    std::map<UInt, std::set<UInt> > dummyFreeSet;
    MapCompleteFctIdToIndex(fctId, freeIndPerBlock, fixedIndPerBlock, true);

    // If there are no affected free dofs, leave immediately
    if( freeIndPerBlock.size() == 0 ) {
      return;
    }
    
    // It's okay, if there are no factors, if there is only a system
    // matrix and no other ones
    if ( matFactors.empty() == true ) {
      if ( matrixTypes_.size() == 1 && sys != NULL ) {
        // Also assemble the effective auxilliary system matrix for moving
        // IDBCs to the right-hand side
        idbcHandler_->BuildSystemMatrix( matFactors, freeIndPerBlock );

        // Now we are done
        return;
      }
      else {
        WARN("SBM_System::ConstructEffectiveMatrix: "
            << "Map with factors is empty, but there are "
            << matrixTypes_.size() << " FE matrices in the game!");
      }
    }
    
    for ( it = matFactors.begin(); it != matFactors.end(); it++ ) {
      if ( sysMat_[(*it).first] != NULL  && (*it).second != 0.0 ) {
        std::map<UInt, std::set<UInt> > dummyFreeSet;
        sys->Add( (*it).second, *sysMat_[(*it).first], 
                  dummyFreeSet, freeIndPerBlock );
      }
    }

    // Also assemble the effective auxilliary system matrix for moving
    // IDBCs to the right-hand side
    idbcHandler_->BuildSystemMatrix( matFactors, freeIndPerBlock );

  }

  template<typename T>
  void AlgebraicSys::SetDirichlet( const FeFctIdType fctId,
                                   Integer eqnNr, const T& val ) {
    
    LOG_DBG(algSys) << "Setting Dirichlet value " << val << " for eqn " << eqnNr
                    << " of fctId " << fctId;
    
    if(eqnNr == 0)
      return;

    UInt blockNr, index;
    MapFctIdEqnToIndex( fctId, eqnNr, blockNr, index );
    
    // Delegate work to IDBC handler
    idbcHandler_->SetIDBC( blockNr, index,  val );
  }

  void AlgebraicSys::BuildInDirichlet() {

    LOG_TRACE(algSys) << "Incorporating Dirichlet values into system matrix";

    // If necessary modify matrix diagonal for penalty approach
    if ( assembleDirichletToSysMat_ == true ) {
      idbcHandler_->AdaptSystemMatrix( *(sysMat_[SYSTEM]) );
      assembleDirichletToSysMat_ = false;
    }
  }

  void AlgebraicSys::GetSolutionVal( SingleVector& ptSol,
                                     const FeFctIdType fctId,
                                     bool setIDBC) {

    LOG_TRACE(algSys) << "Getting solution values of fct " << fctId;

    // get all (blockId,index)-combinations for the current fctId
    StdVector<UInt> blockNums, indices;
    MapCompleteFctIdToIndex( fctId, blockNums, indices);
    UInt size = blockNums.GetSize();
    ptSol.Resize(size);
    ptSol.Init();

    if (domain->GetBasePDE()->GetName() == "LatticeBoltzmann")
    {
    	//FIXME Dirty code, can be improved!
    	ptSol.Resize(size);
    	Vector<Double> & retVec = dynamic_cast<Vector<Double>&>( ptSol );
    	retVec = *(sol_->GetPointer(0));

    	return;
  	}

    if( ptSol.GetEntryType() == BaseMatrix::DOUBLE ) {
      Vector<Double> & retVec = dynamic_cast<Vector<Double>&>( ptSol );
      Double entry = 0.0;

      // #pragma omp parallel for private (entry)
      for( UInt i = 0; i < size; ++i )
      {
        // if index number is larger the lastFree dof, insert Dirichlet value
        // (just if setIDBC is true!)
        if( indices[i] > blockInfo_[blockNums[i]]->numLastFreeIndex )
        {
          if ( setIDBC ) 
            idbcHandler_->GetIDBC(blockNums[i],  indices[i], entry);
          else 
            entry = 0.0;
        }
        else
        {
          sol_->GetPointer(blockNums[i])->GetEntry(indices[i]-1,entry);
        }
        retVec[i] = entry;
      }

    }
    else
    {
      Vector<Complex> & retVec = dynamic_cast<Vector<Complex>&>( ptSol );
      Complex entry = 0.0;
      // #pragma omp parallel for private (entry)
      for( UInt i = 0; i < size; ++i ) {

        // if index number is larger the lastFree dof, insert Dirichlet value
        // (just if setIDBC is true!)
        if( indices[i] > blockInfo_[blockNums[i]]->numLastFreeIndex ) {
          if ( setIDBC) 
            idbcHandler_->GetIDBC(blockNums[i],  indices[i], entry);
          else 
            entry = 0.0;
        } else {
          sol_->GetPointer(blockNums[i])->GetEntry(indices[i]-1,entry);
        }
        retVec[i] = entry;
      }
    }
  }

  
  void AlgebraicSys::GetSolutionVal( SBM_Vector& solVec, bool setIDBC ) {
    
    // resize solVec to match number of functions
    solVec.Resize( numFcts_);
    
    // loop over all feFctIDs
    for(UInt i = 0; i < numFcts_; ++i ) {
    
      // call specialized GetSolutionVal method
      GetSolutionVal(solVec(i), i, setIDBC);
    }
  }
  
  void AlgebraicSys::GetRHSVal( SingleVector &ptRhs,
                                const FeFctIdType fctId  ) {
    
    LOG_TRACE(algSys) << "Getting RHSvalue of fct " << fctId;
    
    // get all (blockId,index)-combinations for the current fctId
    StdVector<UInt> blockNums, indices;
    MapCompleteFctIdToIndex( fctId, blockNums, indices);
    UInt size = blockNums.GetSize();
    ptRhs.Resize(size);
    ptRhs.Init();

    if( ptRhs.GetEntryType() == BaseMatrix::DOUBLE ) {
      Vector<Double> & retVec = dynamic_cast<Vector<Double>&>( ptRhs );
      Double entry = 0.0;

      for( UInt i = 0; i < size; ++i ) {

        // if index number is larger the lastFree dof, insert 0
        if( indices[i] > blockInfo_[blockNums[i]]->numLastFreeIndex) {
          entry = 0.0;
        } else {
          rhs_->GetPointer(blockNums[i])->GetEntry(indices[i]-1,entry);
        }
        retVec[i] = entry;
      }

    } else {

      Vector<Complex> & retVec = dynamic_cast<Vector<Complex>&>( ptRhs );
      Complex entry = 0.0;

      for( UInt i = 0; i < size; ++i ) {

        // if index number is larger the lastFree dof, insert 0
        if( indices[i] > blockInfo_[blockNums[i]]->numLastFreeIndex) {
          entry = 0.0;
        } else {
          rhs_->GetPointer(blockNums[i])->GetEntry(indices[i]-1,entry);
        }
        retVec[i] = entry;
      }
    }
  }
  
  void AlgebraicSys::GetRHSVal( SBM_Vector& rhsVec ) {
    
    // resize rhsVec to match number of functions
    rhsVec.Resize( numFcts_);

    // loop over all feFctIDs
    for(UInt i = 0; i < numFcts_; ++i ) {

      // call specialized GetSolutionVal method
      GetRHSVal(rhsVec(i), i);
    }
      
  }
  
  SBM_Matrix* AlgebraicSys::GenerateSBM_Matrix( FEMatrixType matType,
                                                BaseMatrix::EntryType entryType,
                                                bool sharePattern ) {
    
    LOG_TRACE(algSys) << "Generating SBMMatrix of size " << numBlocks_ 
        << " for matrix type " << feMatrixType.ToString(matType);

    // STEP 1: Generate empty SBM_Matrix
    SBM_Matrix *retMat = NULL;
    retMat = new SBM_Matrix( numBlocks_, numBlocks_, sbmSymm_ );
    if ( retMat == NULL ) {
      EXCEPTION( "SBM_System::GenerateSBM_Matrix: "
          << "This is the end my friend!\n"
          << "Generation of empty SBM_Matrix failed!" );
    }
    
    // STEP 2: Populate with sub-matrices
    std::set<SubMatrixID,SortSubMatrixID>::iterator sIt;
    BaseGraph *graph = NULL;
    for ( sIt = feSubMatricesByBlocks_[matType].begin();
        sIt != feSubMatricesByBlocks_[matType].end(); sIt++ ) {

      // Determine row / col
      UInt sbmRow = (*sIt).rowInd;
      UInt sbmCol = (*sIt).colInd;

      // Determine number of matrix rows and columns
      UInt nrows = blockInfo_[sbmRow]->numLastFreeIndex;
      UInt ncols = blockInfo_[sbmCol]->numLastFreeIndex;

      // Check for necessity of generation
      if ( sbmRow <= sbmCol || sbmSymm_ == false ) {

        graph = graphManager_->GetGraph( sbmRow, sbmCol );

        // Trigger generation of sub-matrix
        if ( sbmRow == sbmCol && sbmSymm_ == true ) {
          // for diagonal blocks we allow a variable
          // matrix layout which we query at the
          // sol-strategy object
          
          BaseMatrix::StorageType sT = solStrat_->GetStorageType(sbmRow);
          LOG_DBG(algSys) << "storage Type of matrix (" << sbmRow +1
              << ", " << sbmCol+1 << ") is " 
              << BaseMatrix::storageType.ToString(sT);
          
          // If we perform static condensation and this is the 
          // inner-inner block, we use the variable block row
          // format to increase performance.
          if( statCond_ && sbmRow == numBlocks_-1 ) {
            sT = BaseMatrix::VAR_BLOCK_ROW;
          }
          
          retMat->SetSubMatrix ( sbmRow, sbmCol, entryType, 
                                 sT,
                                 nrows, ncols, graph->GetNNE() );
        } else {
          // Off-diagonal entries are by nature rectangular and thus, only
          // two possibilities remain: either sparse_nonsym crs format
          // of the variable block row format. 
          // CRS is the more general case and will be preferred.
          retMat->SetSubMatrix ( sbmRow, sbmCol, 
                                 entryType, BaseMatrix::SPARSE_NONSYM,
                                 nrows, ncols, graph->GetNNE() );
          LOG_DBG(algSys) << "storage Type of matrix (" << sbmRow +1
                        << ", " << sbmCol+1 << ") is " 
                        << BaseMatrix::storageType.ToString(BaseMatrix::SPARSE_NONSYM);
        }

        // check, if matrix pattern can be shared and 
        // obtain matrix graph
        if( sharePattern ) {
          LOG_DBG(algSys) << "\tSharing pattern";
          PatternIdType patternID = sbmPatternIds_[(*sIt)]; 
          if( sbmPatternIds_[(*sIt)] != NO_PATTERN_ID ) {
            LOG_DBG(algSys) << "\tObtaining pattern '" << patternID << "' from pool";
            (*retMat)( sbmRow, sbmCol ).SetSparsityPattern( patternPool_, patternID );
          } else {
            (*retMat)( sbmRow, sbmCol ).SetSparsityPattern( *graph );
            sbmPatternIds_[(*sIt)] = 
                (*retMat)( sbmRow, sbmCol ).TransferPatternToPool( patternPool_ );
            LOG_DBG(algSys) << "\tPutting pattern '" << sbmPatternIds_[(*sIt)] << "' to pool";
          }
        } else {
          LOG_DBG(algSys) << "\tUsing no shared sparsity pattern";
          // Set sparsity pattern of sub-matrix
          (*retMat)( sbmRow, sbmCol ).SetSparsityPattern( *graph );  
        }
        
        
      }
    }

    return retMat;
  }
  
  
  void AlgebraicSys::CheckConsistency() {
    
    // First check, if we have a true SBM system.
    // consisting of more than one SBM-Block
    if( numBlocks_ == 1 ||
        (numBlocks_ == 2 && statCond_) ) {
      
      // ========================
      //  Only one block present 
      // ========================
      
      // --------------------------
      //  Check Symmetry of Matrix
      // --------------------------
      // ... if not defined -> choose depending on symmetry of system
      // ... if defined -> check with symmetry of matrix
      PtrParamNode matNode = solStrat_->GetMatrixNode(0);
      BaseMatrix::StorageType storType = BaseMatrix::NOSTORAGETYPE;
      
      // generate preferred storage format, based on symmetry type
      std::string storageString = "sparseSym";
      if(!isDiagBlockSymm_[0] ) 
              storageString = "sparseNonSym";
      
      // check, if symmetry type was set or if we can change it
      bool canChangeMatFormat = true;
      if( matNode->Has("storage") && 
          matNode->Get("storage")->As<std::string>() != "noStorageType" ) {
        canChangeMatFormat = false; 
        matNode->GetValue("storage",storageString);
      } else {
        matNode->Get("storage",ParamNode::INSERT)->SetValue(storageString);
      }
            
      matNode->GetValue("storage",storageString, ParamNode::INSERT);
      storType = BaseMatrix::storageType.Parse(storageString);
      
      // check, if unphysical setting was set
      if( !isDiagBlockSymm_[0] && 
          storType == BaseMatrix::SPARSE_SYM) {
        EXCEPTION( "Can not use storage format 'sparseSym' for an "
                   "unsymmetric system" );
      }
      
      // -------------------------
      //  Check Eigenvalue Solver 
      // -------------------------
      std::string eSolverId = solStrat_->GetEigenSolverId();
      PtrParamNode eSolverList = myParam_->Get("eigenSolverList", 
                                               ParamNode::INSERT);
      ParamNodeList esNodes =  eSolverList->GetChildren();
      PtrParamNode eSolverNode;
      for( UInt i = 0; i < esNodes.GetSize(); ++i ) {
        if( esNodes[i]->Get("id")->As<std::string>() == eSolverId ) {
          eSolverNode = esNodes[i]; 
        }
      }

      // BaseEigenSolver::EigenSolverType est;
      if( !eSolverNode ) {
        // -------------------------------------
        //  no eigensolver set -> use default direct 
        // -------------------------------------
        // est = BaseEigenSolver::ARPACK;
        eSolverList->Get("arpack",ParamNode::INSERT)->
            Get("id",ParamNode::INSERT)->SetValue(eSolverId);
      } else {
        // .. nothing to be done yet
      }

      // --------------
      //  Check Solver
      // --------------
      // .. if not defined -> LDL / ILDL (depending on symmetry)
      // .. if defined -> Check if solver suites symmetry type of matrix
      std::string solverId = solStrat_->GetSolverId();
      PtrParamNode solverList = myParam_->Get("solverList", 
                                              ParamNode::INSERT);
      ParamNodeList sNodes =  solverList->GetChildren();
      PtrParamNode solverNode;
      for( UInt i = 0; i < sNodes.GetSize(); ++i ) {
        if( sNodes[i]->Get("id")->As<std::string>() == solverId ) {
          solverNode = sNodes[i]; 
        }
      }
      BaseSolver::SolverType st;
      // set for allowed matrix types of the solver
      std::set<BaseMatrix::StorageType> solverStorTypes;
      if( !solverNode ) {
        // -------------------------------------
        //  no solver set -> use default direct 
        // -------------------------------------
        st = BaseSolver::PARDISO_SOLVER;
        solverList->Get("pardiso",ParamNode::INSERT)->
          Get("id",ParamNode::INSERT)->SetValue(solverId);
      } else {
        // ---------------------------------------------------
        //  solver set -> check for compatibility with matrix
        // ---------------------------------------------------

        // convert solver string to enum
       st = BaseSolver::solverType.Parse(solverNode->GetName());

        // obtain list of allowed matrix format
        solverStorTypes = GetSolverCompatMatrixFormats(st);

        // check, if current matrix format is in allowed list
        if( solverStorTypes.find(storType) == solverStorTypes.end() &&
            solverStorTypes.size() != 0 ) {
          //  matrix format is not allowed

          //  a) we can change matrix -> do so
          storType = *solverStorTypes.begin();
          if( canChangeMatFormat) {
            storageString = BaseMatrix::storageType.ToString(storType);
            matNode->Get("storage")->SetValue(storageString);
          } else {
            EXCEPTION("Solver '" << solverNode->GetName() 
                      << "' can not operate on matrix with storage type '" 
                      << storageString << "'. \nChange format to '"
                      << BaseMatrix::storageType.ToString(storType)
                      << "'.");
            // b) we can not change matrix -> EXCEPTION
          } // canChangeFormat
        } // find storageType
      } // no solver set

      // ---------------
      //  Check Precond
      // ---------------
      std::string precondId = solStrat_->GetPrecondId();
      PtrParamNode precondList = myParam_->Get("precondList", 
                                               ParamNode::INSERT);
      ParamNodeList pNodes =  precondList->GetChildren();
      PtrParamNode precondNode;
      for( UInt i = 0; i < pNodes.GetSize(); ++i ) {
        if( pNodes[i]->Get("id")->As<std::string>() == precondId ) {
          precondNode = pNodes[i]; 
        }
      }
      
      BaseSolver::PrecondType pt;
      if( !precondNode ) {
        // ----------------------------------
        //  no precond set -> use default ID 
        // ----------------------------------
        pt = BasePrecond::ID;
        precondList->Get("Id",ParamNode::INSERT)->
            Get("id",ParamNode::INSERT)->SetValue(precondId);
      } else {
        // ---------------------------------------------------
        //  precond set -> check for compatibility with matrix
        // ---------------------------------------------------

        // convert precond string to enum
        pt = BaseSolver::precondType.Parse(precondNode->GetName());

        // obtain list of allowed matrix format
        std::set<BaseMatrix::StorageType> mf =
            GetPrecondCompatMatrixFormats(pt);

        // check, if current matrix format is in allowed list
        if( mf.find(storType) == mf.end() &&
            mf.size() != 0 ) {
          //  matrix format is not allowed

          //  a) we can change matrix AND (!!!) the requested
          //     matrix layout is compatible with the solver -> change it

          storType = *mf.begin();
          bool isCompatibleWithSolver =
              (solverStorTypes.size() == 0 ||
                  solverStorTypes.find(storType) != solverStorTypes.end() );

          if( canChangeMatFormat && isCompatibleWithSolver) {
            storageString = BaseMatrix::storageType.ToString(storType);
            matNode->Get("storage")->SetValue(storageString);
          } else {
            EXCEPTION("Precond '" << precondNode->GetName() 
                      << "' can not operate on matrix with storage type '" 
                      << storageString << "'. \nChange format to '"
                      << BaseMatrix::storageType.ToString(storType)
                      << "'.");
            // b) we can not change matrix -> EXCEPTION
          } // canChangeFormat
        } // find storageType
      } // no precond set

      // ---------------------------------------------------
      //  sensibility test for preconditioner
      // ---------------------------------------------------
      // a) all direct solver do not need any preconditioner
      if(( st == BaseSolver::LDL_SOLVER ||
          st == BaseSolver::LU_SOLVER  || 
          st == BaseSolver::LAPACK_LU  ||
          st == BaseSolver::LAPACK_LL  ||
          st == BaseSolver::PARDISO_SOLVER )
          && !(pt == BasePrecond::ID ||
              pt == BasePrecond::NOPRECOND) ) {
        EXCEPTION( "A direct solver only works with the Identity (ID) "
                   "preconditioner." );
      }
      
      // --------------------------
      //  Check for shared pattern
      // --------------------------
      if( st == BaseSolver::DIAGSOLVER ) {
        sharedPatternPossible_ = false;
      }
      
      // ---------------
      //  Check Reordering
      // ---------------
      BaseOrdering::ReorderingType ot = BaseOrdering::SLOAN;
#ifdef USE_METIS
      ot = BaseOrdering::METIS;
#endif
      bool canChangeReordering = true;
      if (matNode->Has("reordering") && 
          matNode->Get("reordering")->As<std::string>() != "_default_" ) {
        ot = BaseOrdering::reorderingType.Parse(
            matNode->Get("reordering")->As<std::string>());
        canChangeReordering = false;
      }
      
      // a) for our own direct solvers we activate re-ordering
      if( (st == BaseSolver::LU_SOLVER ||
           st == BaseSolver::LDL_SOLVER || 
           st == BaseSolver::LAPACK_LL || 
           st == BaseSolver::LAPACK_LU ) &&
           ot == BaseOrdering::NOREORDERING && 
          canChangeReordering == true ) {
#ifdef USE_METIS
        ot = BaseOrdering::METIS;
#else
        ot = BaseOrdering::SLOAN;
#endif
      }
      
      // b) pardiso and most external solvers need no reordering or have their own
      if( st == BaseSolver::PARDISO_SOLVER &&
          st == BaseSolver::UMFPACK &&
          st == BaseSolver::ILUPACK &&
//          st == BaseSolver::LIS &&
          st == BaseSolver::SUPERLU &&
          st == BaseSolver::SPOOLES &&
          ot != BaseOrdering::NOREORDERING &&
          canChangeReordering == true ) {
        ot = BaseOrdering::NOREORDERING;
      }
      
      // c) ilu-based preconditioners prefer reordering
      if( ( pt == BasePrecond::ILUK ||
            pt == BasePrecond::ILUTP ||
            pt == BasePrecond::ILDLK ||
            pt == BasePrecond::ILDLTP ||
            pt == BasePrecond::ILDLCN ) &&
            ot == BaseOrdering::NOREORDERING && 
            canChangeReordering == true ) {
#ifdef USE_METIS
        ot = BaseOrdering::METIS;
#else
        ot = BaseOrdering::SLOAN;
#endif
      }
      
      // in the end store back the reordering type
      matNode->Get("reordering",ParamNode::INSERT)->
          SetValue(BaseOrdering::reorderingType.ToString(ot));
      
    } else {
      // ======================
      //  True SBM Case 
      // ======================
      WARN("This section is not yet implemented");
    }
    
    // Dump tree in the end
//    std::cerr << "Dump of parameter tree at the end of "
//        << "AlgebraicSys::CheckConsistency():\n";
//    myParam_->Dump();
  }
  
  void AlgebraicSys::PrintFeMatrixInfo( ) {

    LOG_TRACE(algSys) << "Print matrix information";
    
    PtrParamNode setupNode = myInfo_->Get("setup");
    
    // add timer
    setupNode->Get("setupTime")->SetValue(graphTimer_);
    
    // Print overview of defined matrices
    setupNode->SetComment("List of defined matrices");
    PtrParamNode matrixListNode = setupNode->Get("matrices");
    matrixListNode->SetComment("Memory is in MByte");
    
    // Loop over all FeMatrixTypes
    std::map<FEMatrixType, SubMatrixSet>::iterator feMatIt
    = feSubMatricesByBlocks_.begin();
    for( ; feMatIt != feSubMatricesByBlocks_.end(); ++feMatIt ) {
    
      // generate <STIFFNES> .. node
      PtrParamNode matTypeNode = 
          matrixListNode->Get(feMatrixType.ToString(feMatIt->first));
      
      SubMatrixSet & smSet = feMatIt->second;
      
      UInt totalNumRows, totalNumCols, totalNumNonZeros;
      totalNumRows = totalNumCols = totalNumNonZeros = 0;
      Double totalFillLevelPerCent = 0.0;
      Double totalMemory = 0.0;
      // Loop overall matrix blocks
      SubMatrixSet::iterator smIt = smSet.begin();
      for( ; smIt != smSet.end(); ++smIt ) {

        // obtain sub-matrix
        SBM_Matrix & sbmMat = *sysMat_[feMatIt->first];
        StdMatrix * stdMat = NULL;
        SubMatrixID smId = *smIt;
        stdMat = sbmMat.GetPointer(smId.rowInd,smId.colInd);

        // if std-matrix is not defined, leave
        if( !stdMat)
          continue;

        PtrParamNode mNode = matTypeNode->Get("matrix",ParamNode::APPEND);
        

        mNode->Get("blockRow")->SetValue(smId.rowInd);
        mNode->Get("blockCol")->SetValue(smId.colInd);
              
        std::string storageType = 
            BaseMatrix::storageType.ToString(stdMat->GetStorageType());
        mNode->Get("storageType")->SetValue(storageType);
        mNode->Get("numRows")->SetValue(stdMat->GetNumRows());
        mNode->Get("numCols")->SetValue(stdMat->GetNumCols());
        mNode->Get("numNonZeros")->SetValue(stdMat->GetNnz());
        Double fillLevel = stdMat->GetNnz() * 100 / 
            (Double(stdMat->GetNumRows()) * Double(stdMat->GetNumCols()) ); 
        mNode->Get("fillLevelPerCent")->SetValue(fillLevel);
        
        Double actMemory = stdMat->GetMemoryUsage();
        mNode->Get("memory")->SetValue( actMemory/1024./1024. );
        totalMemory += actMemory;
        // get reordering strategy from solution strategy object
        
        
        // special check for VBR-block matrix type
        if( stdMat->GetStorageType() == BaseMatrix::VAR_BLOCK_ROW ) {
          UInt nbRows, nbCols, nBlocks, numOffDiagEntries, effNNZ;
          Double avgBlockSize;
          
          if( stdMat->GetEntryType() == BaseMatrix::DOUBLE ) {
            VBR_Matrix<Double> & vMat = 
                dynamic_cast<VBR_Matrix<Double> &>(*stdMat);
            vMat.GetNumBlocks( nbRows, nbCols, nBlocks );
            effNNZ = vMat.GetEffNnz();
            avgBlockSize = vMat.GetAvgRowBlockSize();
            numOffDiagEntries = vMat.GetNumOffDiagonalEntries();
          } else {
            VBR_Matrix<Complex> & vMat = 
                dynamic_cast<VBR_Matrix<Complex> &>(*stdMat);
            vMat.GetNumBlocks( nbRows, nbCols, nBlocks );
            effNNZ = vMat.GetEffNnz();
            avgBlockSize = vMat.GetAvgRowBlockSize();
            numOffDiagEntries = vMat.GetNumOffDiagonalEntries();
          }
          mNode->Get("numBlockRows")->SetValue(nbRows);
          mNode->Get("numBlockCols")->SetValue(nbCols);
          mNode->Get("numBlocks")->SetValue(nBlocks);
          mNode->Get("avgRowBlockSize")->SetValue(avgBlockSize);
          mNode->Get("numNonZerosEff")->SetValue(effNNZ);
          mNode->Get("numOffDiagEntries")->SetValue(numOffDiagEntries);
        }
        
       
        totalNumNonZeros += stdMat->GetNnz();
        BaseGraph * graph = graphManager_->GetGraph(smId.rowInd,smId.colInd); 
        
        // bandwidth and reordering gets just written for diagonal blocks 
        if( smId.rowInd == smId.colInd ) {
          UInt bwLow = 0, bwUp = 0, bwAvg = 0;
          graph->GetBandwidth(bwLow, bwUp, bwAvg);
          mNode->Get("upperBandWidth")->SetValue(bwUp);
          mNode->Get("lowerBandWidth")->SetValue(bwLow);
          mNode->Get("avgBandWidth")->SetValue(bwAvg);
          mNode->Get("symCounterPart")->SetValue("no");
          totalNumRows += stdMat->GetNumRows();
          totalNumCols += stdMat->GetNumCols();
        } else {
          mNode->Get("symCounterPart")->SetValue(sbmMat.IsSymmetric());
        }
        
        std::string orderStr = 
            BaseOrdering::reorderingType.ToString(graph->GetReorderType());
        mNode->Get("orderingType")->SetValue(orderStr);
        
      } // matrix blocks
      
      totalFillLevelPerCent = totalNumNonZeros * 100 / (Double( totalNumRows) * 
                              Double(totalNumCols));
      matTypeNode->Get("totalNumRows")->SetValue(totalNumRows);
      matTypeNode->Get("totalNumCols")->SetValue(totalNumCols);
      matTypeNode->Get("totalNumNonZeros")->SetValue(totalNumNonZeros);
      matTypeNode->Get("totalFillLevelPerCent")->SetValue(totalFillLevelPerCent);
      matTypeNode->Get("totalMemory")->SetValue(totalMemory/1024./1024.);
    } // feMatrixType
  }

  // *************************
  //   PrintRegistrationInfo
  // *************************
  void AlgebraicSys::PrintRegistrationInfo( ) const {
    
    LOG_TRACE(algSys) << "Print registration info";
    
    PtrParamNode setupNode = myInfo_->Get("setup");
    
    // Print overview of feFunctions
    setupNode->SetComment("List of registered FeFunctions");
    PtrParamNode fctListNode = setupNode->Get("feFunctions");
    std::map<FeFctIdType,std::string>::const_iterator it = fctNames_.begin();
    
    UInt totalSize = 0, totalNumFreeEqns = 0, totalNumDirichlet = 0;
    
    for( ; it != fctNames_.end(); ++it ) {
      PtrParamNode fctNode = fctListNode->Get("function",ParamNode::APPEND);
      FeFctIdType fctId = it->first;
      fctNode->Get("id")->SetValue(fctId);
      fctNode->Get("name")->SetValue(it->second);
      
      fctNode->Get("numEqns")->SetValue(numEqnsPerFct_[fctId]);
      fctNode->Get("lastFreeEqn")->SetValue(lastFreeEqnPerFct_[fctId]);
      
      UInt numDirichlet = numEqnsPerFct_[fctId] - lastFreeEqnPerFct_[fctId];
      fctNode->Get("numDirichlet")->SetValue(numDirichlet);
      
      totalNumFreeEqns += lastFreeEqnPerFct_[fctId];
      totalSize += numEqnsPerFct_[fctId];
      totalNumDirichlet += numDirichlet;
      
      // List also feBlocks, in which this functions occurs
      std::map<FeFctIdType, std::set<UInt> >::const_iterator mapIt;
      mapIt = fctIdsInBlocks_.find(it->first);
      const std::set<UInt> & blockSet = (mapIt->second);
      std::set<UInt>::const_iterator blockIt = blockSet.begin();
      PtrParamNode blockListNode = fctNode->Get("usedInBlocks");
      for( ; blockIt != blockSet.end(); ++blockIt ) {
        blockListNode->Get("block",ParamNode::APPEND)
            ->Get("id")->SetValue(*blockIt);
      }
    }
    
    // add totalSize, totalNumDirichlet
    fctListNode->Get("totalNumEqns")->SetValue(totalSize);
    fctListNode->Get("totalNumFreeEqns")->SetValue(totalNumFreeEqns);
    fctListNode->Get("totalNumDirichlet")->SetValue(totalNumDirichlet);
    
    // Print overview of blocks
    setupNode->SetComment("List of SBM-blocks");
    PtrParamNode blockListNode = setupNode->Get("sbmBlocks");
    
    for( UInt i = 0; i < numBlocks_; ++i ) {
      PtrParamNode blockNode = blockListNode->Get("block",ParamNode::APPEND);
      blockNode->Get("id")->SetValue(i);
      blockNode->Get("size")->SetValue(blockInfo_[i]->size);
      blockNode->Get("lastFreeIndex")->SetValue(blockInfo_[i]->numLastFreeIndex);
    }
  }
  
  // ========================================================================
  // Explicit template instantiation
  // ========================================================================
  template void AlgebraicSys::
  SetElementMatrix( FEMatrixType, Matrix<Double>&, 
                    FeFctIdType, const StdVector<Integer>& ,
                    FeFctIdType, const StdVector<Integer>& , bool, bool, bool);
  template void AlgebraicSys::
  SetElementMatrix( FEMatrixType, Matrix<Complex>&, 
                    FeFctIdType, const StdVector<Integer>& ,
                    FeFctIdType, const StdVector<Integer>& , bool, bool, bool);
  
  template void AlgebraicSys::
  SetElementRHS( const Vector<Double>&, const FeFctIdType, 
                 StdVector<Integer>&);
  template void AlgebraicSys::
    SetElementRHS( const Vector<Complex>&, const FeFctIdType, 
                   StdVector<Integer>&);
  
  template void AlgebraicSys::SetNodeRHS(Double, FeFctIdType, Integer );
  template void AlgebraicSys::SetNodeRHS(Complex, FeFctIdType, Integer );
  
  template void AlgebraicSys::SetFncRHS(  const Vector<Double>& fncRHS, FeFctIdType fctId );
  template void AlgebraicSys::SetFncRHS(  const Vector<Complex>& fncRHS, FeFctIdType fctId );

  template void AlgebraicSys::
  AddToDiagMatrixEntry( FEMatrixType, const FeFctIdType, Integer, Double );
  template void AlgebraicSys::
  AddToDiagMatrixEntry( FEMatrixType, const FeFctIdType, Integer, Complex );
  
  template void AlgebraicSys::
  GetMatrixEntry( FEMatrixType, const FeFctIdType, Integer,
                  const FeFctIdType, Integer, Double& );
  template void AlgebraicSys::
  GetMatrixEntry( FEMatrixType, const FeFctIdType, Integer,
                  const FeFctIdType, Integer, Complex& );

   template void AlgebraicSys::
   SetMatrixEntry( FEMatrixType, const FeFctIdType, Integer,
                   const FeFctIdType, Integer, Double, bool );
   template void AlgebraicSys::
   SetMatrixEntry( FEMatrixType, const FeFctIdType, Integer,
                   const FeFctIdType, Integer, Complex, bool );
   
   template void AlgebraicSys::
   SetDirichlet( const FeFctIdType, Integer, const Double& );
   template void AlgebraicSys::
   SetDirichlet( const FeFctIdType, Integer, const Complex& );
  
}// end of Namespace
