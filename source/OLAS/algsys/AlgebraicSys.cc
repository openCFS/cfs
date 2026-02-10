#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include <def_use_metis.hh>
#include <def_use_pardiso.hh>
#include <def_use_arpack.hh>
#include <def_use_phist_cg.hh>
#include <def_use_phist_ev.hh>
#include <def_reordering.hh>

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

#include "MatVec/CRS_Matrix.hh"
#include "MatVec/SCRS_Matrix.hh"

DEFINE_LOG(algSys, "algSys")

namespace CoupledField {




  // ***********************
  //   Default Constructor
  // ***********************
  AlgebraicSys::AlgebraicSys(PtrParamNode param, PtrParamNode info, bool isSolutionComplex, bool isMultHarm)
  {
    size_ = 0;
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

    isMultHarm_ = false;
    isMultHarm_ = isMultHarm;

    // Setup solution strategy enum
    PtrParamNode stratNode = myParam_->Get("solutionStrategy",ParamNode::INSERT);
    solStrat_ = SolStrategy::Generate(stratNode);

    // Inform solution strategy if we use multiharmonic approach
    solStrat_->SetMultHarm(isMultHarm_);

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

    //delete eigenValues_;
    eigenValues_ = NULL;
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

  StdVector<std::string> AlgebraicSys::GetFeFunctionsNames()
  {
    // std::map<FeFctIdType,std::string> fctNames_;
    StdVector<std::string> res(fctNames_.size());
    for(int i = 0; i < (int) res.GetSize(); i++) // we assume the fe functions to be ordered from 0
    {
      assert(fctNames_.find(i) != fctNames_.end());
      res[i] = fctNames_[(FeFctIdType)i];
    }
    return res;
  }

  void AlgebraicSys::UpdateToSolStrategy() {
    LOG_DBG(algSys) << "Updating parameters due to solution strategy";

    if(isMultHarm_) EXCEPTION("AlgebraicSys::UpdateToSolStrategy() cannot handle multiharmonic case yet!")

    // switch according to type of solution strategy
    if( solStrat_->GetType() == SolStrategy::TWO_LEVEL_STRATEGY ) {

      // In case of a two-level strategy, we have to adjust the
      // "effective" matrix / vector size, depending on the step.

      if( solStrat_->GetActSolStep() == 1
          && solStrat_->GetNumSolSteps() == 2) {
        // --------------------------
        //  Step 1: Only (1,1) block
        // --------------------------
        LOG_DBG(algSys) << "\t=> Switching to reduced (1,1)-system";

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
        LOG_DBG(algSys) << "\t=> Switching to full system again";
        delete effMat_;
        delete effSol_;
        delete effRhs_;

        // re-status of static condensation
        statCond_ = solStrat_->UseStaticCondensation();

        if (statCond_) {
          effMat_ = new SBM_Matrix( *sysMat_[SYSTEM], numBlocks_-1, numBlocks_-1 );
          effRhs_ = new SBM_Vector( *rhs_, numBlocks_-1 );
          effSol_ = new SBM_Vector( *sol_, numBlocks_-1 );
        } else {
          effMat_ = new SBM_Matrix( *sysMat_[SYSTEM], numBlocks_, numBlocks_ );
          effRhs_ = new SBM_Vector( *rhs_, numBlocks_ );
          effSol_ = new SBM_Vector( *sol_, numBlocks_ );
        }
      } else {
        EXCEPTION("The two level solution strategy has only two steps");
      }
    } // if TWO_LEVEL strategy

  }


  void AlgebraicSys::CreateLinSys() {

    LOG_DBG(algSys) << "Creating linear system";

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
      // for multiharmonic analysis, we have to hack a bit
      if( isMultHarm_ ){
        // ------------------------------
        //  Multiharmonic Case
        // ------------------------------
        UInt M = solStrat_->GetNumHarmM();

        SubMatrixID id;
        UInt Nmax = domain->GetDriver()->GetNumFreq();
        UInt a = (domain->GetDriver()->IsFullSystem())? (M+1) : ((M-1)/2 + 1);
        for( UInt iRow = 0; iRow < Nmax; ++iRow ) {
          id.rowInd = iRow;
          id.colInd = iRow;
          sbmPatternIds_[id] = NO_PATTERN_ID;
          for( UInt iCol = iRow + 1; iCol < iRow + a; ++iCol ) {
            if( iCol <  Nmax){
              id.rowInd = iRow;
              id.colInd = iCol;
              LOG_DBG2(algSys) << "(row,col)=("<<iRow<<","<<iCol<<")";
              sbmPatternIds_[id] = NO_PATTERN_ID;
              // transposed
              id.rowInd = iCol;
              id.colInd = iRow;
              LOG_DBG2(algSys) << "(row,col)=("<<iCol<<","<<iRow<<")";
              sbmPatternIds_[id] = NO_PATTERN_ID;
            }
          }
        }
      }else{
        // ------------------------------
        //  Classic Case
        // ------------------------------
        for( UInt iRow = 0; iRow < numBlocks_; ++iRow ) {
          for( UInt iCol = 0; iCol < numBlocks_; ++iCol ) {
            SubMatrixID id;
            id.rowInd = iRow;
            id.colInd = iCol;
            sbmPatternIds_[id] = NO_PATTERN_ID;
          }
        }
      }
    } // endif sharedPatterPossible_

    // Obtain some info from parameter file
    BaseMatrix::EntryType entryType = isMatrixComplex_ ?
        BaseMatrix::COMPLEX :
        BaseMatrix::DOUBLE;
    for ( fIt = matrixTypes_.begin(); fIt != matrixTypes_.end(); fIt++ ) {
      sysMat_[*fIt] = GenerateSBM_Matrix( *fIt, entryType, sharePattern );
    }

    // Log what we will do.
    // This should work always, but with commit f487c43e Florian mentioned
    // "commented out because it produces a segfault - only for small problems ???"
    // if this is fixed, one could print the information always
    if(progOpts->DoDetailedInfo())
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
      LOG_DBG3(algSys) << "GenerateIDBC_HandlerObject"<< std::endl;
      LOG_DBG3(algSys) << "numDirichletValuesPerBlock_= "<< numDirichletValuesPerBlock_ <<std::endl;

      if (isMultHarm_){
        //For each entry in the multiharmonic matrix a IDBC-graph gets generated
        UInt M = solStrat_->GetNumHarmM();
        UInt a = (domain->GetDriver()->IsFullSystem())? (M+1) : ((M-1)/2 + 1);

        idbcHandler_ = GenerateIDBC_HandlerObjectMH( matrixTypes_, graphManager_,
                                              numDirichletValuesPerBlock_,
                                              entryType, M, a );
      }else {
      idbcHandler_ =
          GenerateIDBC_HandlerObject( matrixTypes_, graphManager_,
                                      numDirichletValuesPerBlock_,
                                      entryType );
      }
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
    UInt nB = (isMultHarm_)? domain->GetDriver()->GetNumFreq() : numBlocks_;
    for ( UInt k = 0; k < nB; k++ ) {
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
      effMat_ = new SBM_Matrix( *sysMat_[SYSTEM], nB, nB );
      effRhs_ = new SBM_Vector( *rhs_, nB );
      effSol_ = new SBM_Vector( *sol_, nB );
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

    LOG_DBG(algSys) << "Creating preconditioner";

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

    LOG_DBG(algSys) << "Creating solver";

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

    LOG_DBG(algSys) << "Creating eigenvalue solver";
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
          GenerateEigenSolverObject(  solStrat_,
                                     esNode, sNode, pNode,
                                     myInfo_->Get("solve_eigen") );
    } else {
      EXCEPTION("Eigenvalue solver can currently only handle SBM "
          << "matrices with 1 block!");
    }
  }

  void AlgebraicSys::SetupPrecond()  {

    LOG_DBG(algSys) << "Setup of preconditioner";
    // check, if system was already created
    if( !systemCreated_ ) {
      EXCEPTION( "Matrices were not created yet. Please call "
          "AlgebraicSys::CreateLinSys() first!" );
    }

    // start setup timer of preconditioner
    precond_->GetSetupTimer()->Start();

    // if we have just one SBM matrix block, use directly
    // the specialized methods for StdMatrices
    if( useAMG_ ){
        precond_->SetupMG( (*effMat_)(0,0), *auxMatAMG_, amgType_, edgeIndNode_, nodeNumIndex_);
    }else{
      if( onlyOneMatrixBlock_ ) {
        precond_->Setup( (*effMat_)(0,0));
      } else {
        precond_->Setup( *(effMat_));
      }
    }

    // stop setup timer of preconditioner
    precond_->GetSetupTimer()->Stop();

  }

  bool AlgebraicSys::HasPrecond() {
    if( precond_->GetPrecondType() == BasePrecond::NOPRECOND || precond_->GetPrecondType() == BasePrecond::ID ) {
      return false;
    } else {
      return true;
    }
  }

  void AlgebraicSys::SetupSolver() {

    LOG_DBG(algSys) << "Setup of solver";
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

  void AlgebraicSys::SetupEigenSolver(UInt numFreq, Double shift, bool isQuadratic, bool sort, bool bloch) {

    LOG_DBG(algSys) << "Setup of eigenvalue solver";
    // check, if system was already created
    if( !systemCreated_ ) {
      EXCEPTION( "Matrices were not created yet. Please call "
          "AlgebraicSys::CreateLinSys() first!" );
    }

    // start setup timer of solver
    eigenSolver_->GetSetupTimer()->Start();

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
      eigenSolver_->Setup((*sysMat_[STIFFNESS])(0,0), (*sysMat_[MASS])(0,0), (*sysMat_[DAMPING])(0,0), numFreq, shift, sort);
    } else {
      if( dampPresent == true ) {
        WARN("Although a damping matrix is present, only a generalized "
            << "eigenvalue will be solved (as you have specified)!");
      }


      if( massPresent == true ) {
        // Setup the eigenvalue solver for generalized EV problem
        eigenSolver_->Setup((*sysMat_[STIFFNESS])(0,0), (*sysMat_[MASS])(0,0), numFreq, shift, sort, bloch);
      } else {
        // Setup the eigenvalue solver for standard EV problem
        eigenSolver_->Setup((*sysMat_[STIFFNESS])(0,0), numFreq, shift, sort);
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

    // stop setup timer of solver
    eigenSolver_->GetSetupTimer()->Stop();
  }

  void AlgebraicSys::SetOldDirichletValues() {
   // idbcHandler_->ToString();
    idbcHandler_->SetOldDirichletValues();
  }

  void AlgebraicSys::AddIDBCToRHS(bool deltaIDBC) {
    LOG_DBG(algSys) << "Add IDBC to RHS ";

    idbcHandler_->AddIDBCToRHS( rhs_, deltaIDBC );
  }

  void AlgebraicSys::Solve(bool setIDBC, bool deltaIDBC) {

    LOG_DBG(algSys) << "Solving problem";

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

      BaseEigenSolver * evs = GenerateEigenSolverObject( solStrat_,
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
      solStrat_->GetSetupNode()->Get("calcConditionNumber")->SetValue( std::any(false) );
    }
    // ======================================================================

    if(domain->GetBasePDE()->GetName() == "LatticeBoltzmann") // we have to adjust size of solution vector for LBM optimization
    {
      LOG_DBG(algSys) << "Resized solution vector for LBM optimization to " << (*sysMat_[SYSTEM])(0,0).GetNumRows();
      sol_->GetPointer(0)->Resize((*sysMat_[SYSTEM])(0,0).GetNumRows());
    }

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
      idbcHandler_->SetDofsToIDBC( effSol_, deltaIDBC );
    }

    // Assume that everything will go well
    PtrParamNode out = myInfo_->Get(ParamNode::PROCESS)->Get("solver");
    out->Get("solutionIsOkay")->SetValue(true);

    // Now modify the right-hand side vector.
    // Note: It is mandatory to incorporate the IDBC values to the
    // complete RHS.
    if ( setIDBC ){
      idbcHandler_->AddIDBCToRHS( rhs_, deltaIDBC );
    }


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
      idbcHandler_->RemoveIDBCFromRHS( rhs_, deltaIDBC );

    // Check that solution went fine, if not issue a warning
    if ( out->Get("solutionIsOkay")->As<bool>() == false ) {
      WARN("Solver reports a problem! Consult .las file for "
          << "further diagnostics!");
    }


    ExportLinSys(false, false, true); // post_solve

    // stop timer associated with solver
    solver_->GetSolveTimer()->Stop();
  }

  void AlgebraicSys::ExportMHSys(int step){
    std::string logString = "sys_step_";
    logString.append( boost::lexical_cast<std::string>(step));
    sysMat_[SYSTEM]->Export(logString, BaseMatrix::OutputFormat::MATRIX_MARKET, "sys");

    //logString = "mass_step_";
    //logString.append( boost::lexical_cast<std::string>(step));
    //sysMat_[DAMPING]->Export(logString, BaseMatrix::OutputFormat::MATRIX_MARKET, "mass");
  }

  void AlgebraicSys::ExportLinSys(bool setup, bool pre_solve, bool post_solve)
  {
    assert((setup && !pre_solve && !post_solve) || (!setup && pre_solve && !post_solve) || (!setup && !pre_solve && post_solve));

    AnalysisID& id = domain->GetDriver()->GetAnalysisId();

    LOG_DBG(algSys) << "ELS setup=" << setup << " pre=" << pre_solve << " post=" << post_solve << " id=" << id.ToString();

    PtrParamNode els = GetExportLinSysParam();
    if(els == NULL)
      return;

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

      if(els->Get("geometric_stiffness")->As<bool>() && sysMat_.find(GEOMETRIC_STIFFNESS) != sysMat_.end() && sysMat_[GEOMETRIC_STIFFNESS] != NULL)
        sysMat_[GEOMETRIC_STIFFNESS]->Export(base + "_geometric_stiffness", mat_format);

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

  PtrParamNode AlgebraicSys::GetExportLinSysParam(){
    if(!solStrat_->GetParamNode()->Has("exportLinSys")){
      return NULL;
    }
    else {
      return solStrat_->GetParamNode()->Get("exportLinSys");
    }
  }



  void AlgebraicSys::CalcEigenFrequencies(Vector<Double>& frequencies, Vector<Double>& err)
  {
    LOG_DBG(algSys) << "Calculating real-valued eigenfrequencies";

    // start timer of solver
    eigenSolver_->GetSolveTimer()->Start();

    // Trigger calculation of eigenvalues
    eigenSolver_->CalcEigenFrequencies( *eigenValues_, *eigenValError_ );

    // Hard coded cast for double values
    Vector<Double>& valVec = dynamic_cast<Vector<Double>&>(*eigenValues_);
    frequencies = valVec;

    Vector<Double>& errVec = dynamic_cast<Vector<Double>&>(*eigenValError_);
    err = errVec;

    ExportLinSys(false, false, true); // the modes are actually not the solution

    // stop timer associated with solver
    eigenSolver_->GetSolveTimer()->Stop();
  }

  void AlgebraicSys::CalcEigenFrequencies(Vector<Complex>& frequencies, Vector<Double>& err)
  {
    LOG_DBG(algSys) << "Calculating complex-valued eigenfrequencies: bloch=" << eigenSolver_->IsBloch() << " quadratic=" << eigenSolver_->IsQuadratic();

    // start timer of solver
    eigenSolver_->GetSolveTimer()->Start();

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

    // stop timer associated with solver
    eigenSolver_->GetSolveTimer()->Stop();
  }

  void AlgebraicSys::CalcEigenValues(BaseVector &sol, BaseVector &err, Double minVal, Double maxVal){
    // start timer of solver
    eigenSolver_->GetSolveTimer()->Start();

    eigenSolver_->CalcEigenValues(sol,err,minVal,maxVal);
    //SingleVector ev = dynamic_cast< SingleVector & > sol;
    //SingleVector * sv = dynamic_cast<BaseVector &>(sol);
    //eigenValues_ = dynamic_cast<SingleVector &>(sol);
    //eigenValues_ = dynamic_cast<SingleVector>(*sol);
    //eigenValues_ = dynamic_cast<SingleVector>(&sol);// target is not pinter or ref
    eigenValues_ = dynamic_cast<SingleVector *>(&sol);
    ExportLinSys(false, false, true);

    // stop timer associated with solver
    eigenSolver_->GetSolveTimer()->Stop();
  }

  void AlgebraicSys::CalcEigenValues(BaseVector &sol, BaseVector &err){
    // start timer of solver
    eigenSolver_->GetSolveTimer()->Start();

    eigenSolver_->CalcEigenValues(sol,err);
    //SingleVector ev = dynamic_cast< SingleVector & > sol;
    //SingleVector * sv = dynamic_cast<BaseVector &>(sol);
    //eigenValues_ = dynamic_cast<SingleVector &>(sol);
    //eigenValues_ = dynamic_cast<SingleVector>(*sol);
    //eigenValues_ = dynamic_cast<SingleVector>(&sol);// target is not pinter or ref
    eigenValues_ = dynamic_cast<SingleVector *>(&sol);
    ExportLinSys(false, false, true);

    // stop timer associated with solver
    eigenSolver_->GetSolveTimer()->Stop();
  }

  void AlgebraicSys::GetEigenMode( UInt numMode, bool right )  {

    LOG_DBG(algSys) << "GEM #" << numMode;
    Vector<Complex>& solHelp = dynamic_cast<Vector<Complex> &>((*sol_)(0));

    eigenSolver_->GetNormalizedEigenMode(numMode, solHelp, right);

    LOG_DBG2(algSys) << "GEM -> " << sol_->ToString();
  }

  void AlgebraicSys::GraphSetupInit( UInt numFcts,
                                     bool useDistinctGraphs ) {

    LOG_DBG(algSys) << "Setup matrix graph for " << numFcts << " functions.";
    LOG_DBG(algSys) << "Use distinct graphs:" << useDistinctGraphs << std::endl;

    distinctMatGraphs_ = useDistinctGraphs;

    // feFunction specific data
    numFcts_ = numFcts;
    numEqnsPerFct_.Resize(numFcts);
    lastFreeEqnPerFct_.Resize(numFcts);
    eqnToSBMBlock_.Resize(numFcts);

    // store flag for applying algebraic multigrid precond/solver
    useAMG_ = false;
    if(myParam_->Has("precondList")){
      useAMG_ = myParam_->Get("precondList")->Has("MG");
    }

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

    // Neither does it work with multiharmonic analysis
    if(isMultHarm_ && statCond_){
      EXCEPTION("Static condensation cannot be applied in a multiharmonic analysis!");
    }

  }


  // ***************
  //   ObtainFctId
  // ***************
  FeFctIdType AlgebraicSys::ObtainFctId( const std::string& fctString ) {

   LOG_DBG(algSys) << "Obtaining FctId for fct '" << fctString << "'";

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


  Integer AlgebraicSys::DefineSBMMatrixBlock( const std::map<FeFctIdType,std::set<Integer> >& eqns, bool isInnerBlock) {
    LOG_DBG(algSys) << "Defining new SBM block #" << numBlocks_;

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
      LOG_DBG(algSys) << "\tBlock is empty, leaving";
      // in addition, if this block is supposed to be the static condensation block,
      // we deactivate it

      // Notify also solution strategy about empty block ....
      WARN("Implement dropping of empty block ...");

      if( isInnerBlock && statCond_) {
        statCond_ = false;
        LOG_DBG(algSys) << "\tDeactivating static condensation";
      }
      return -1;
    }

    sbmIndex = numBlocks_;
    numBlocks_++;
    isDiagBlockSymm_.Push_back(true);
    numDirichletValuesPerBlock_.Push_back(0);

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
        std::unordered_map<UInt, UInt> & eqnToIndex = bi.eqnToIndex[iFct];
        std::unordered_map<UInt, UInt>::iterator eqnToIndexIt = eqnToIndex.begin();
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

    // return index of newly created block
    return sbmIndex;
  }

  void AlgebraicSys::RegisterSubMatrixBlocks( UInt sbmIndex, UInt numMinorBlocks ) {

    LOG_DBG(algSys) << "Registering " << numMinorBlocks << " sub-matrix blocks for SBM block #" << sbmIndex;

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

    LOG_DBG(algSys) << "Setting matrix type '" << feMatrixType.ToString(matrixType)
                           << "' for fct-Ids (" << fctId1 << ", " << fctId2 << ")";

    LOG_DBG(algSys) << "\tsymmetry: " << isSymmetric;
    LOG_DBG(algSys) << "\tcomplex values: " << isComplex;


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
          LOG_DBG(algSys) << "\t=> matrix will begic";
        }

      } else {
        // === OFF-DIAGONAL BLOCK ===

        // For an off-diagonal entry the symmetric-flag denotes, if the
        // transposed matrix is set as well. Thus, we can determine the overall
        // symmetry of the SBM-Matrix. In case at least one integrator
        // is non-symmetric, so will be the SBM matrix.
        this->sbmSymm_ &= isSymmetric;
        LOG_DBG(algSys) << "\t=> SBM-symmetry: " << this->sbmSymm_;

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

    LOG_DBG(algSys) << "Finished setup of graph";

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
    std::map<FEMatrixType, SubMatrixSet>::const_iterator itMatByFct;
    for( itMatByFct = feSubMatricesByFctId_.begin(); itMatByFct != feSubMatricesByFctId_.end(); ++itMatByFct ) {

      const FEMatrixType & matrixType = itMatByFct->first;
      const SubMatrixSet & sbmSet = itMatByFct->second;
      SubMatrixSet::const_iterator sbmIt;

      // loop over all blocks
      for( sbmIt = sbmSet.begin() ; sbmIt != sbmSet.end(); ++sbmIt ) {

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
        std::unordered_map<UInt, UInt> & eqnToIndex = bi.eqnToIndex[iFct];
        std::unordered_map<UInt, UInt>::iterator it = eqnToIndex.begin();

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
  }




  // ******************
  //   GraphSetupDoneMH
  // ******************
  void AlgebraicSys::GraphSetupDoneMH() {

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

    // number of harmonics
    UInt M = solStrat_->GetNumHarmM();
    UInt a = (domain->GetDriver()->IsFullSystem())? (M+1) : ((M-1)/2 + 1);

    // loop over all matrix Types in feSubMatricesByFctId
    std::map<FEMatrixType, SubMatrixSet>::const_iterator itMatByFct;
    for( auto & itMatByFct : feSubMatricesByFctId_ ) {
      const FEMatrixType & matrixType = itMatByFct.first;

      // Fetch all SBM blocks, in which the current (rowFctId, colFctId) occurs.
      // In the multiharmonic case, give it all the non-zero (sbmRow, sbmCol) combinations
      for( UInt sbmRow = 0; sbmRow < domain->GetDriver()->GetNumFreq(); ++sbmRow ) {
        for( UInt sbmCol = sbmRow; sbmCol < sbmRow + a; ++sbmCol ) {
          if( sbmCol < domain->GetDriver()->GetNumFreq() ){
            SubMatrixID sID;
            sID.rowInd = sbmRow;
            sID.colInd = sbmCol;
            // Insert sub-matrix identifier into corresponding FE-Matrix set
            feSubMatricesByBlocks_[matrixType].insert( sID );

            //also set the transposed
            if(sbmRow != sbmCol){
              sID.rowInd = sbmCol;
              sID.colInd = sbmRow;
              // Insert sub-matrix identifier into corresponding FE-Matrix set
              feSubMatricesByBlocks_[matrixType].insert( sID );
            }
          }
        }
      }

    } // loop over matrixType

    // --------------------------------------------------------------------
    // Determine the set of sub-matrices that we need for the SYSTEM matrix.
    // There are two different cases:
    //
    // 1) We have other FE matrices: NOT IMPLEMENTED IN MULTIHARMONIC CASE
    //
    // 2) There are no other sets, so we generate a sub-matrix for each
    //    sub-graph that the graph manager stores
    // --------------------------------------------------------------------
    UInt nnzBlocks = 0; //number of nonzero blocks
    if ( matrixTypes_.size() > 3 ) {
      EXCEPTION("AlgebraicSys::GraphSetupDoneMH() more than two matrix type...not implemented yet \n"
                "for the multiharmonic case!" );
    }
    else {
      SubMatrixID sID;
      for (  UInt sbmRow = 0; sbmRow < domain->GetDriver()->GetNumFreq(); ++sbmRow ) {
        for ( UInt sbmCol = sbmRow; sbmCol < sbmRow + a; ++sbmCol ) {
          if( sbmCol < domain->GetDriver()->GetNumFreq()){
            ++nnzBlocks;
            if ( graphManager_->SubGraphExists( sbmRow, sbmCol ) == true ) {
              sID.rowInd = sbmRow;
              sID.colInd = sbmCol;
              feSubMatricesByBlocks_[SYSTEM].insert( sID );
              feSubMatricesByBlocks_[DAMPING].insert( sID );
            }

            // also handle the transpose
            if(sbmRow != sbmCol){
              ++nnzBlocks;
              if ( graphManager_->SubGraphExists( sbmCol, sbmRow ) == true ) {
                sID.rowInd = sbmCol;
                sID.colInd = sbmRow;
                feSubMatricesByBlocks_[SYSTEM].insert( sID );
                feSubMatricesByBlocks_[DAMPING].insert( sID );
              }
            }

          }
        }
      }
    }

    // determine overall number of entries (including fixed equations)
    size_ = 0;
    for( UInt i = 0; i < nnzBlocks; i++ ) {
      size_ += blockInfo_[0]->size;
    }
    StdVector< StdVector<UInt> > toBeCopied( domain->GetDriver()->GetNumFreq()  * domain->GetDriver()->GetNumFreq() );
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
    for( auto & fctIt : matIsSymm_  ){

      FeFctIdType fctId = fctIt.first;
      bool isFctSymm = fctIt.second;

      // loop over all diagonal SBMBlocks, where this fctId occurs
      // and perform logical AND operation regarding symmetry (i.e.
      // if at least one function in this block is non-symmetric, so
      // will be the complete block)
      std::set<UInt>& affectedBlocks = fctIdsInBlocks_[fctId];
      for( auto & blockIt : affectedBlocks) {
        isDiagBlockSymm_[blockIt] &= isFctSymm;
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
    graphManager_->SetupDoneMH(reorder, solStrat_->GetNumHarmN(), solStrat_->GetNumHarmM() );


    // Now we have all graphs and IDBC in their re-ordered state,
    // so we have to fetch the reordering array from the GraphManager and
    // update information in the blockInfo array for all SBM-Blocks
    GraphManager::SBMBlockInfo &bi = *blockInfo_[0];
    UInt numFcts = bi.eqnToIndex.GetSize();
    StdVector<UInt> newOrder;

    // Since all diagonal blocks have the same reordering and only the
    // blockInfo_ with index 0 is valid, we only have to perform the following
    // loop once. It is left as a loop for further implementations, where we
    // might apply different reorderings to the blocks
    for (  UInt sbmRow = 0; sbmRow < 1; ++sbmRow ) {
      // ----------------------------------
      //   D I A G O N A L    B L O C K S
      // ----------------------------------
      // Obtain reordering vector
      // take care, the following method needs the row-number
      graphManager_->GetReordering(sbmRow, newOrder);


      // Loop over all functions
      for( UInt iFct = 0; iFct < numFcts; ++iFct ) {
        std::unordered_map<UInt, UInt> & eqnToIndex = bi.eqnToIndex[iFct];
        std::unordered_map<UInt, UInt>::iterator it = eqnToIndex.begin();

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

    }

  }

  // **********************
  //   FinishRegistration
  // **********************
  void AlgebraicSys::FinishRegistration( ) {

    LOG_DBG(algSys) << "FinishRegistration";

    // create new graph manager object and initialize it
    graphManager_ = new GraphManager();

    // Different setup for multiharmonic analysis
    if( isMultHarm_ ){
      UInt numSBMRows = domain->GetDriver()->GetNumFreq();
      graphManager_->SetupInit( numSBMRows, distinctMatGraphs_, true,
                                solStrat_->GetNumHarmN(), solStrat_->GetNumHarmM(),
                                domain->GetDriver()->GetNumFreq(), solStrat_->IsFullSystem()  );

      // In the multiharmonic case, we have only one set of equations
      // but they are present several times in the final system matrix
      // (for the different frequencies)
      graphManager_->RegisterBlockMultHarm( blockInfo_[0]);

      // "real" SBM case
      onlyOneMatrixBlock_ = false;

    }else{
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

//    std::cout << "Setting element position for fctIds ("
//                     << fctId1 << ", " << fctId2 << ")" << std::endl;
//    std::cout << "matrixType: " << feMatrixType.ToString(matrixType) << std::endl;
//    std::cout << "counterPart: " << (setCounterPart ? "yes" : "no") << std::endl;
//    std::cout << "EqnVec1: " << eqnNrs1.ToString() << std::endl;
//    std::cout << "EqnVec2: " << eqnNrs2.ToString() << std::endl;

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
    if( isMultHarm_ ){
      MapFctIdEqnToIndex_MultHarm(fctId1, eqnNrs1, rowBlocks, rowNums, nnzSBMInd_);
      MapFctIdEqnToIndex_MultHarm(fctId2, eqnNrs2, colBlocks, colNums, nnzSBMInd_);
    }else{
//      std::cout << "MapFctIdEqnToIndex" << std::endl;
      MapFctIdEqnToIndex(fctId1, eqnNrs1, rowBlocks, rowNums);
      MapFctIdEqnToIndex(fctId2, eqnNrs2, colBlocks, colNums);
    }

    // Quirk: If fctId1 == fctId2, we normally
    // need not set the counterPart. If however, they are now in
    // off diagonal matrices, we have to ensure, that originally
    // symmetric matrix structures remain symmetric also in the
    // blockSystem
    if( fctId1 == fctId2 )
      setCounterPart = true;


//    std::cout << "setCounterPart? " << setCounterPart << std::endl;
    graphManager_->SetElementPos( rowBlocks, rowNums,
                                  colBlocks, colNums,
                                  matrixType,
                                  setCounterPart);

  }


  void AlgebraicSys::MapFctIdEqnToIndex_MultHarm( const FeFctIdType fctId,
                                                   const StdVector<Integer>& eqns,
                                                   StdVector<UInt>& blockNums,
                                                   StdVector<UInt>& indices,
                                                   const StdVector<UInt>& sbmIndices) {
    LOG_DBG(algSys) << "MFIETI Mapping fctId,eqnNr to blockNum,indices";

    blockNums.Resize(eqns.GetSize() * sbmIndices.GetSize() );
    indices.Resize(eqns.GetSize() * sbmIndices.GetSize() );

    UInt numEqns = eqns.GetSize();
    UInt blockCnt = 0;
    for( UInt iEqn = 0; iEqn < numEqns; ++iEqn ) {
      const UInt & eqnNr = std::abs(eqns[iEqn]);

      // take care of homogeneous BCs
      if( eqnNr == 0) {
        // TODO check if this is correct
        //WARN("Homogeneous BC's not yet tested for multiharmonic analysis!!!");
        for(UInt i = 0; i < sbmIndices.GetSize(); ++i){
          blockNums[blockCnt] = 0;
          // multiharmonic only one blockInfo
          indices[blockCnt] = 0;
          ++blockCnt;
        }
      } else {
        for(auto blockInd : sbmIndices){
          blockNums[blockCnt] = blockInd;
          // multiharmonic only one blockInfo
          indices[blockCnt] = blockInfo_[0]->eqnToIndex[fctId][eqnNr];
          ++blockCnt;
        }

      }
    }
  }


  void AlgebraicSys::MapFctIdEqnToIndex_MultHarm( const FeFctIdType fctId,
                                                  const StdVector<Integer>& eqns,
                                                  StdVector<UInt>& blockNums,
                                                  StdVector<UInt>& indices ) {
    LOG_DBG(algSys) << "MFIETI Mapping fctId,eqnNr to blockNum,indices";

    blockNums.Resize(eqns.GetSize());
    indices.Resize(eqns.GetSize());


    UInt numEqns = eqns.GetSize();
    for( UInt iEqn = 0; iEqn < numEqns; ++iEqn ) {
      const UInt & eqnNr = std::abs(eqns[iEqn]);

      // take care of homogeneous BCs
      if( eqnNr == 0) {
        // TODO check if this is correct
        //WARN("Homogeneous BC's not yet tested for multiharmonic analysis!!!");
        blockNums[iEqn] = 0;
        indices[iEqn] = 0;
      } else {
        // note: in multiharmonic analysis only one blockInfo
        const UInt & blockNum = 0;
        blockNums[iEqn] = blockNum;
        indices[iEqn] = blockInfo_[blockNum]->eqnToIndex[fctId][eqnNr];
      }
    }

  }

  void AlgebraicSys::MapFctIdEqnToIndex( const FeFctIdType fctId,
                                         const StdVector<Integer>& eqns,
                                         StdVector<UInt>& blockNums,
                                         StdVector<UInt>& indices ) {
    LOG_DBG(algSys) << "MFIETI Mapping fctId,eqnNr to blockNum,indices";

    blockNums.Resize(eqns.GetSize());
    indices.Resize(eqns.GetSize());
//
//    std::cout << "eqns.GetSize() = " << eqns.GetSize() << std::endl;
//
    // get hold of fct-specific map
    StdVector<UInt>& eqnToBlock = eqnToSBMBlock_[fctId];

//    std::cout << "FeFctIdType = " << fctId << std::endl;
//    std::cout << "eqnToBlock (size) = " << eqnToBlock.GetSize() << std::endl;
//    std::cout << "blockInfo_ (size) = " << blockInfo_.GetSize() << std::endl;

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
        const std::unordered_map<UInt, UInt> & eqnToIndexSet =
            blockInfo_[iBlock]->eqnToIndex[actFctId];
        std::unordered_map<UInt, UInt>::const_iterator eqnIt = eqnToIndexSet.begin();
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

    if(fctId == NO_FCT_ID)
    {
      std::map<FeFctIdType, std::set<UInt> >::iterator it;
      it = fctIdsInBlocks_.begin();
      for(; it != fctIdsInBlocks_.end(); ++it)
        fctIds.insert(it->first);
    }
    else
      fctIds.insert(fctId);


    //blockNums.Resize(size);
    //indices.Resize(size);

    Integer offset = 0;
    if(indexZeroBased)
      offset = -1;


    // loop over all functionIds
    std::set<FeFctIdType>::const_iterator fctIt = fctIds.begin();
    for( ; fctIt != fctIds.end(); ++fctIt ) {

      const FeFctIdType actFctId = *fctIt;

      // loop over all blocks
      for( UInt iBlock = 0; iBlock < numBlocks_; ++iBlock ) {
        const std::unordered_map<UInt, UInt> & eqnToIndexSet =
            blockInfo_[iBlock]->eqnToIndex[actFctId];
        std::unordered_map<UInt, UInt>::const_iterator eqnIt = eqnToIndexSet.begin();
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

    LOG_DBG(algSys) << "Initializing matrix " << feMatrixType.ToString(matrixType)
                      << " for fctId " << fctId;

    // If matrix specified init this one
    if ( matrixType != NOTYPE )
    {
      //TODO: ok, that the matrix is zero is not uncommon in case
      //      of transient analysis. Two ways to avoid this:
      //      first: assemble calls this function only with
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

    LOG_DBG(algSys) << "Initializing RHS for fctId " << fctId;

    if ( fctId == NO_FCT_ID ) {
      // in this case initialize complete RHS
      rhs_->Init();
    }
    else {
      if(solStrat_->IsMultHarm()){
        EXCEPTION("This branch of AlgebraicSys::InitSol is not"
                  "meant to be reached in a multiharmonic analysis!");
      }
      // find out affected blocks
      std::set<UInt> & blockNums = fctIdsInBlocks_[fctId];
      std::set<UInt>::iterator it = blockNums.begin();
      for( ; it != blockNums.end(); ++it ) {
        rhs_->GetPointer( *it )->Init();
      }
    }
  }




  void AlgebraicSys::InitRHS( const SBM_Vector& newRHS ) {

    LOG_DBG(algSys) << "Initializing RHS with new vector";

    // ensure that the RHS vector to set consists of as many
    // sub-vectors as the RHS of the system
    if( (newRHS.GetSize() != numFcts_ && domain->GetBasePDE()->GetName() != "LatticeBoltzmann") && !solStrat_->IsMultHarm() ) {
      EXCEPTION( "New rhs consists of " << newRHS.GetSize() << " sub-vectors, the RHS of the algebraic system of " << rhs_->GetSize() << " entries." )
    }

    if (domain->GetBasePDE()->GetName() == "LatticeBoltzmann"){
      *(rhs_) = newRHS;
      LOG_DBG(algSys) << "InitRHS: Initilized rhs with vector v=" << rhs_->GetPointer(0)->ToString();
      return;
    }

    if(solStrat_->IsMultHarm() ){

      // get all (blockId,index)-combinations for the current fctId
      StdVector<UInt> blockNums, indices;
      MapCompleteFctIdToIndex( 0, blockNums, indices);
      UInt size = blockNums.GetSize();

      for(UInt i = 0; i < newRHS.GetSize(); ++i ) {
        // security check: ensure that sub-vector has the same size
        // as the block indices
        if( newRHS(i).GetSize() != indices.GetSize() && (domain->GetBasePDE()->GetName() != "LatticeBoltzmann")) {
          EXCEPTION( "Number of entries of " << i << "-th sub-vector and number of indices do not match!");
        }

        Vector<Complex> & nRHS = dynamic_cast<Vector<Complex>&>( newRHS(i) );
        for( UInt j = 0; j < size; ++j ) {
          // omit entries for Dirichlet values
          if( indices[j] <= blockInfo_[blockNums[j]]->numLastFreeIndex) {
            rhs_->GetPointer(i)->SetEntry(indices[j]-1, nRHS[j] );
          }
        }
      }


      LOG_DBG(algSys) << "InitRHS: Initilized rhs with vector v=" << rhs_->ToString();
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
        Vector<Double> & nRHS = dynamic_cast<Vector<Double>&>( newRHS(i) );
        for( UInt j = 0; j < size; ++j ) {
          // omit entries for Dirichlet values
          if( indices[j] <= blockInfo_[blockNums[j]]->numLastFreeIndex) {
            rhs_->GetPointer(blockNums[j])->SetEntry(indices[j]-1, nRHS[j] );
          }
        }
      } else if( newRHS.GetEntryType() == BaseMatrix::COMPLEX ) {
        Vector<Complex>& nRHS = dynamic_cast<Vector<Complex>&>(newRHS(i));
        for( UInt j = 0; j < size; ++j )
        {
          // omit entries for Dirichlet values
          if( indices[j] <= blockInfo_[blockNums[j]]->numLastFreeIndex) {
            rhs_->GetPointer(blockNums[j])->SetEntry(indices[j]-1, nRHS[j]);
          }
        }
      } else {
        EXCEPTION("Data type must be either double or complex");
      }
    }
  }

  void AlgebraicSys::InitSol( const FeFctIdType fctId ) {

    LOG_DBG(algSys) << "Initializing solution of fctId " << fctId;

    if ( fctId == NO_FCT_ID ) {
      // in this case initialize complete RHS
      sol_->Init();
    }
    else {
      if(solStrat_->IsMultHarm()){
        EXCEPTION("This branch of AlgebraicSys::InitSol is not"
                  "meant to be reached in a multiharmonic analysis!");
      }
      // find out affected blocks
      std::set<UInt> & blockNums = fctIdsInBlocks_[fctId];
      std::set<UInt>::iterator it = blockNums.begin();
      for( ; it != blockNums.end(); ++it ) {
        sol_->GetPointer( *it )->Init();
      }
    }
  }

  void AlgebraicSys::InitSol( const SBM_Vector& newSol ) {

    LOG_DBG(algSys) << "Initializing solution with new vector";

    // ensure that the solution vector to set consists of as many
    // sub-vectors as the solution of the system
    if( (newSol.GetSize() != numFcts_ && domain->GetBasePDE()->GetName() != "LatticeBoltzmann") && !solStrat_->IsMultHarm() ) {
      EXCEPTION( "New rhs consists of " << newSol.GetSize() << " sub-vectors, the RHS of the algebraic system of " << rhs_->GetSize() << " entries." )
    }

    if(solStrat_->IsMultHarm() ){

      // get all (blockId,index)-combinations for the current fctId
      StdVector<UInt> blockNums, indices;
      MapCompleteFctIdToIndex( 0, blockNums, indices);
      UInt size = blockNums.GetSize();

      for(UInt i = 0; i < newSol.GetSize(); ++i ) {
        // security check: ensure that sub-vector has the same size
        // as the block indices
        if( newSol(i).GetSize() != indices.GetSize()) {
          EXCEPTION( "Number of entries of " << i << "-th sub-vector and number of indices do not match!");
        }

        Vector<Complex> & nSol = dynamic_cast<Vector<Complex>&>( newSol(i) );
        for( UInt j = 0; j < size; ++j ) {
          // omit entries for Dirichlet values
          if( indices[j] <= blockInfo_[blockNums[j]]->numLastFreeIndex) {
            sol_->GetPointer(i)->SetEntry(indices[j]-1, nSol[j] );
          }
        }
      }
    }else{
      EXCEPTION("AlgebraicSys::InitSol currently only allowed in multiharmonic analysis!!!")
    }

    LOG_DBG(algSys) << "InitSol: Initilized solution with vector v=" << sol_->ToString();
    return;
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

    LOG_DBG(algSys) << "Setting element matrix for fctIds (" << fctId1 << ", " << fctId2 << ")";
    LOG_DBG2(algSys) << "Matrix: " << feMatrixType.ToString(matrixType);
    LOG_DBG2(algSys) << "EqnVec1: (" << eqnNrs1.GetSize() << "): " << eqnNrs1.ToString();
    LOG_DBG2(algSys) << "EqnVec2: (" << eqnNrs2.GetSize() << "): " << eqnNrs2.ToString();
    LOG_DBG3(algSys) << "elemMat (" << elemMat.GetNumRows() << ", " << elemMat.GetNumCols() << ")"; //:\n " << elemMat;
    LOG_DBG3(algSys) << "noStaticCond is: " << noStaticCond;
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
     LOG_DBG3(algSys) << "sbmInd/sbmIndices : " << sbmRow +1 << "/[" << numBlocks_ << "]" ;

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
  void AlgebraicSys::SetElementMatrix_MultHarm( FEMatrixType matrixType,
                                                Matrix<T>& elemMat,
                                                FeFctIdType fctId1,
                                                const StdVector<Integer>& eqnNrs1,
                                                FeFctIdType fctId2,
                                                const StdVector<Integer>& eqnNrs2,
                                                bool setCounterPart,
                                                const StdVector<UInt>& sbmIndices) {

    if(fctId1 != fctId2) EXCEPTION("AlgebraicSys::SetElementMatrix_MultHarm function Id's don't match!");

    // lambda for converting flattened sbm-index to (row, col) tuple
    UInt size = domain->GetDriver()->GetNumFreq();
    auto DeflattenIndex = [size](UInt ind) { std::vector<UInt> a = {ind / (size), ind % (size)}; return a;};
    //auto FlattenIndex = [N](UInt row, UInt col) { return N * row + col;};

    std::string t;
    if (IS_LOG_ENABLED(algSys, dbg3)) {
      // construct logging output
      for(auto s : sbmIndices){
        t.append("(");
        t.append( boost::lexical_cast<std::string>(DeflattenIndex(s)[0]) );
        t.append( ", ");
        t.append( boost::lexical_cast<std::string>(DeflattenIndex(s)[1]) );
        t.append("), ");
      };
    }

    LOG_DBG(algSys) << "Setting element matrix for fctIds ("
                     << fctId1 << ", " << fctId2 << ")";
    LOG_DBG2(algSys) << "Matrix: " << feMatrixType.ToString(matrixType);
    LOG_DBG2(algSys) << "EqnVec1: (" << eqnNrs1.GetSize() << "): " << eqnNrs1.ToString();
    LOG_DBG2(algSys) << "EqnVec2: (" << eqnNrs2.GetSize() << "): " << eqnNrs2.ToString();
    LOG_DBG3(algSys) << "elemMat (" << elemMat.GetNumRows() << ", " << elemMat.GetNumCols() << "):\n " << elemMat;
    LOG_DBG3(algSys) << "in SBM blocks: " << t <<"\n";

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


    // TODO still open problem (performance) :
    /* Problem in multiharmonic analysis
     * In the classic SetElementMatrix method, every (fctId,eqnNr) occurs
     * exactly once. But here it occurs in every block,
     * according to sbmIndices.
     *
     * Therefore we currently call the special method MapFctIdEqnToIndex_MultHarm,
     * which simply appends the CfsTLS vectors and then we continue as in the classic case.
     *
     * Maybe there is a more performant way but we would need to introduce a different concept.
     *
     */

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

    // Re-map entries from (fctId,eqnNr) -> (index)
    // Re-map entries from (fctId,eqnNr) -> (index)
    MapFctIdEqnToIndex_MultHarm(fctId1, eqnNrs1, rowBlocks, rowNums);
    MapFctIdEqnToIndex_MultHarm(fctId2, eqnNrs2, colBlocks, colNums);

    rowIndList1[0].Clear(true);
    rowList1[0].Clear(true);
    rowIndList2[0].Clear(true);
    rowList2[0].Clear(true);
    colIndList1[0].Clear(true);
    colList1[0].Clear(true);
    colIndList2[0].Clear(true);
    colList2[0].Clear(true);

    UInt numRows = rowBlocks.GetSize();
    UInt numCols = colBlocks.GetSize();

    // Compute index of graph in graph pointer matrix
    // get hold of vertex and edgelists
    StdVector<UInt> & rList1 = rowList1[0];
    StdVector<UInt> & rIndList1 = rowIndList1[0];
    StdVector<UInt> & rList2 = rowList2[0];
    StdVector<UInt> & rIndList2 = rowIndList2[0];
    // get hold of vertex and edgelists
    StdVector<UInt> & cList1 = colList1[0];
    StdVector<UInt> & cIndList1 = colIndList1[0];
    StdVector<UInt> & cList2 = colList2[0];
    StdVector<UInt> & cIndList2 = colIndList2[0];

    // Loop over all indices
    for( UInt i = 0; i < numRows; ++i ) {
      // get hold of block numbers and indices
      const UInt & rowNum = rowNums[i];
      // get limits of free indices
      // remember: in multiharmonic analysis we only have one blockInfo
      const UInt & lastFreeRowIndex = blockInfo_[0]->numLastFreeIndex;

      // STEP 1: Generate row index list from first connect array, dropping
      //         equation numbers for dofs fixed by (in)homogeneous Dirichlet
      //         boundary conditions and changing the sign of those fixed by
      //         constraints.
      if ( rowNum > 0 ) {
        if ( rowNum > lastFreeRowIndex ) {
          rList2.Push_back( rowNum - lastFreeRowIndex - 1 );
          rIndList2.Push_back( i );
        } else {
          rList1.Push_back( rowNum - 1);
          rIndList1.Push_back( i );
        }
      }
    }

    // Loop over all columns
    for( UInt i = 0; i < numCols; ++i ) {
      // get hold of block numbers and indices
      const UInt & colNum = colNums[i];

      // get limits of free indices
      // remember: in multiharmonic analysis we only have one blockInfo
      const UInt & lastFreeColIndex = blockInfo_[0]->numLastFreeIndex;

      // STEP 2: Split the second connect array into two edge lists, one for
      //         the graph and one for the IDBCgraph (which handles the indices
      //         fixed by inhomogeneous Dirichlet boundary conditions)
      if( colNum > 0 ) {
        if ( colNum > lastFreeColIndex ) {
          cList2.Push_back( colNum - lastFreeColIndex - 1);
          cIndList2.Push_back( i );
        }
        else {
          cList1.Push_back( colNum - 1);
          cIndList1.Push_back( i );
        }
      }
    } // loop over cols


    SBM_Matrix * actMat = sysMat_[matrixType];
    UInt rowInd, colInd;


    // ======================================================================

    // loop over all blocks and pass for every block the information to
    // the corresponding graph / IDBC graph
    LOG_DBG3(algSys) << "setting matrix entries";

    // now loop over every sbm-block specified in sbmIndices parameter
    for(auto sbmInd : sbmIndices){
      LOG_DBG3(algSys) << "sbmInd/sbmIndices : " << sbmInd << "/[" << sbmIndices.ToString()<< "]" ;
      UInt sbmRow = DeflattenIndex(sbmInd)[0];
      UInt sbmCol = DeflattenIndex(sbmInd)[1];
      LOG_DBG3(algSys) << "\tsetting SBM block (" << sbmRow
                << "," << sbmCol << ") with sbm-index " << sbmInd;

      StdMatrix * stdMat = actMat->GetPointer(sbmRow, sbmCol);

      LOG_DBG3(algSys) << "\t1) free-free entries:";
      LOG_DBG3(algSys) << "\t\trowIndices: " << rList1.ToString();
      LOG_DBG3(algSys) << "\t\tcolIndices: " << cList1.ToString();
      LOG_DBG3(algSys) << "\t\tmat: " << stdMat->ToInfoString();

      // Attention: This check is not really implemented in a clean way!
      if( stdMat != NULL ) {
        LOG_DBG3(algSys) << "\t1) free-free entries:";
        LOG_DBG3(algSys) << "\t\trowIndices: " << rList1.ToString();
        LOG_DBG3(algSys) << "\t\tcolIndices: " << cList1.ToString();
        LOG_DBG3(algSys) << "\t\tmat: " << stdMat->ToInfoString();

        // 2) Assemble all free <-> free entries

        for ( UInt i = 0; i < rList1.GetSize(); i++ ) {
          rowInd = rIndList1[i];
          for ( UInt j = 0; j < cList1.GetSize(); j++ ) {
            colInd = cIndList1[j];
            stdMat->AddToMatrixEntry( rList1[i], cList1[j], elemMat[rowInd][colInd] );
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
              stdMat->AddToMatrixEntry( cList1[j], rList1[i], elemMat[rowInd][colInd] );
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
            idbcHandler_->AddWeightFixedToFree( matrixType, sbmRow, sbmCol, rList1[i], cList2[j], elemMat[rowInd][colInd]);
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
              idbcHandler_->AddWeightFixedToFree( matrixType, sbmCol, sbmRow, cList1[j], rList2[i], elemMat[rowInd][colInd]);
            } // j
          } // i
        } // sbmCol == sbmRow
      } // rList2.GetSize()

    }

  }


  template<typename T>
  void AlgebraicSys::SetElementRHS( const Vector<T>& elemRHS,
      const FeFctIdType fctId,
      StdVector<Integer>& eqnNrs, UInt& harm ) {

    LOG_DBG(algSys) << "SER: Setting element RHS for fctId ("<< fctId << ")";
    LOG_DBG2(algSys) << "SER: vector " << elemRHS.GetSize() << " -> " << elemRHS.ToString();
    LOG_DBG2(algSys) << "SER: eqnNrs " << eqnNrs.GetSize() << " -> " << eqnNrs.ToString();
    // Ensure that there are as many equations as vector entries
    //assert(eqnNrs.GetSize() == elemRHS.GetSize());

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

      // Get vector, differentiate between normal and multiharmonic case
      // index:     [  0     1     2  ... N-1  N    N+1   N+2 ...  2N ]
      // harmonic:  [ -N   -N+1  -N+2 ... -1   0     1     2  ...   N ]
      // NOTE: If the performance optimized version is used,
      // only odd harmonics are considered and therefore the above mapping looks like
      // index:     [  0     1     2  ... (N-1)/2   (N+1)/2   (N-1)/2+2   ...  N+1 ]
      // harmonic:  [ -N   -N+2  -N+4 ...    -1        1           3      ...   N ]
      if( isMultHarm_ ){
        SingleVector &vecP = (*rhs_)( domain->GetDriver()->IndexOfHarmonic(harm));
        SingleVector &vecN = (*rhs_)( domain->GetDriver()->IndexOfHarmonic(-harm));
        if(vecP.GetEntryType() == BaseMatrix::COMPLEX && vecN.GetEntryType() == BaseMatrix::COMPLEX){
          if ( rowNum > 0 && rowNum <= lastFreeRowIndex ) {
            if ( rowNum <= lastFreeRowIndex &&  elemRHS[iRow] != (Complex)0.0 ) {
              // If we want excitation in the real part (corresponds to cosine excitation)
              vecP.AddToEntry( rowNum-1, elemRHS[iRow]/2.0);
              //this entry must be conjugate complex
              vecN.AddToEntry( rowNum-1, std::conj(elemRHS[iRow])/2.0);
            }
          } // loop over rows
        }else EXCEPTION("This error when filling the multiharm rhs should not happen");

        LOG_DBG3(algSys) << "SER: rhs is:\n " << (*rhs_).ToString();
        LOG_DBG3(algSys) << "vecP: \n " << vecP.ToString();
        LOG_DBG3(algSys) << "vecN: \n " << vecP.ToString();

      } else{
        SingleVector &vec = (*rhs_)(rowBlock);
        if ( rowNum > 0 && rowNum <= lastFreeRowIndex ) {
          if ( rowNum <= lastFreeRowIndex ) {
            vec.AddToEntry( rowNum-1, elemRHS[iRow]);
          }
        } // loop over rows
      }


    } // loop over blocks
  }


  template<typename T>
  void AlgebraicSys::SetNodeRHS(T val, FeFctIdType fctId,
                                Integer eqnNr) {

    LOG_DBG(algSys) << "Setting node RHS of " << eqnNr << " for fct "
                    << fctId << " to " << val;

    if(isMultHarm_) EXCEPTION("AlgebraicSys::SetNodeRHS cannot handle multiharmonic case yet!")

    UInt block,idx;
    this->MapFctIdEqnToIndex(fctId,eqnNr,block,idx);
    rhs_->GetPointer(block)->AddToEntry(idx-1,val);
  }

  template<typename T>
  void AlgebraicSys::SetFncRHS(  const Vector<T>& fncRHS, FeFctIdType fctId ) {

    LOG_DBG(algSys) << "Setting Function RHS for fctId ("<< fctId << ")";

    if(isMultHarm_) EXCEPTION("AlgebraicSys::SetFncRHS cannot handle multiharmonic case yet!")

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

  void AlgebraicSys::BackupCurrentSystemMatrix(FEMatrixType storageLocation) {
    CopyMatrixToOther(FEMatrixType::SYSTEM, storageLocation);
  }

  void AlgebraicSys::LoadBackupToCurrentSystemMatrix(FEMatrixType storageLocation){

    if(matrixTypes_.find(storageLocation) == matrixTypes_.end()) {
      WARN("AlgebraicSystem::LoadBackupToCurrentSystemMatrix --- Cannot find backup at storageLocation "<<FEMatrixType(storageLocation));
      return;
    }
    assert(matrixTypes_.find(FEMatrixType::SYSTEM) != matrixTypes_.end());
//    std::cout << "Restore deep copy of current system (" << storageLocation << ")" << std::endl;
    CopyMatrixToOther(storageLocation, FEMatrixType::SYSTEM);
  }

  void AlgebraicSys::CopyMatrixToOther(FEMatrixType matrix, FEMatrixType other, bool add) {
    if(matrixTypes_.find(other) == matrixTypes_.end()){
      matrixTypes_.insert(other);
      // create deep copy of matrix
      SBM_Matrix* copy = new SBM_Matrix(*(sysMat_[matrix]));
      sysMat_[other] = copy;
    } else {
      if (!add)
        // clear storage
        sysMat_[other]->Init();
      // create deep copy is not a good idea as sysmat will point to a new object afterwards!
      sysMat_[other]->Add(1.0,*(sysMat_[matrix]));
    }
  }

  void AlgebraicSys::ComputeSysMatTimesVector(FEMatrixType matrixType, SBM_Vector& inputVec, SBM_Vector& outputVec, bool transpose){
    /*
     * -> use already existing functionality to perform multiplication of
     *    inputVec with system matrix, i.e. perform multiplication
     *    by taking the route over UpdateRHS and GetRHS
     * -> Reason:
     *      the orderings of inputVec (sorted by equation number) and system matrix
     *      (ordered/sorted by inner/outer dof or some other special way) do
     *      not match
     *      > we have to resort inputVec first, multiply and then sort back
     *
     *      additionally, we have to take care about inhomogeneours dirichlet
     *      boundary conditions as those influence the system matrix, too
     *      (e.g. by elimination of rows/columns)
     *      > correct treatment and addtion of dof to rhs has to be considered
     *
     *    > this is already done in updateRHS
     */

    /*
     * a) backup current rhs and load outputVec into rhs
     */
    SBM_Vector rhsBAK(BaseMatrix::DOUBLE);
    GetRHSVal(rhsBAK);
    rhs_->Init();

    /*
     * b) compute matrix*inputVec
     * (the first flag SysMatUpdated is not checked in updateRHS at all)
     */
    UpdateRHS(matrixType,inputVec,true,transpose);

    /*
     * c) load rhs into outputVec, then restore rhs
     * (has to be done via GetRHSVal as we need reordering of dofs again)
     */

    GetRHSVal(outputVec);
    InitRHS(rhsBAK);
  }

  void AlgebraicSys::UpdateRHS(FEMatrixType matrixType,
                               const SBM_Vector& fup,bool SysMatUpdated,bool useTransposed) {

    LOG_DBG(algSys) << "Updating RHS of matrix "
                      << feMatrixType.ToString(matrixType);

    if(isMultHarm_) EXCEPTION("AlgebraicSys::UpdateRHS cannot handle multiharmonic case yet!")

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

    /*
     * Important note:
     * we cannot just simply multiply the fup with the matrix as they use
     * different order of the unknowns
     * > so in a first step, go over fup and rearrange entries!
     * > at the same time, we can also deal with idbc value which have to be
     * sorted out and put to the rhs (except for penalty case)
     */
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
          LOG_DBG2(algSys) << "func i= " << i << ", j= " << j << ", nRHS[j]= " << nRHS[j];
        }

      }
      else if( fup.GetEntryType() == BaseMatrix::COMPLEX ) {
          Vector<Complex> & nRHS =
              dynamic_cast<Vector<Complex>&>( fup(i) );

          for( UInt j = 0; j < size; ++j ) {
            // omit entries for Dirichlet values
            if( indices[j] <= blockInfo_[blockNums[j]]->numLastFreeIndex) {
              tmpRHS_->GetPointer(blockNums[j])
                            ->AddToEntry(indices[j]-1, nRHS[j] );
            }else if(!usingPenalty_){
              idbcHandler_->AddFixedToFreeRHS(matrixType,blockNums[j],
                      indices[j],rhs_,nRHS[j]);
            }
            LOG_DBG2(algSys) << "func i= " << i << ", j= " << j << ", nRHS[j]= " << nRHS[j];
          }

      }
      else {
        EXCEPTION("Implement me. Dont worry: mostly C&P code");
      }

    }

    //now just perform multiplication using the REORDERED vector tmpRHS
    if(useTransposed == false){
      sysMat_[matrixType]->MultAdd(*tmpRHS_,*rhs_);
    } else {
      // for Linesearch algorithms based on Wolfe conditions we need to compute
      // Jac^T * residual
      // where Jac is approximated by the systemmatrix if e.g. Newton method is
      // used or in case of Hysteresisis using a deltaMat formulation
      sysMat_[matrixType]->MultTAdd(*tmpRHS_,*rhs_);
    }

  }



  void AlgebraicSys::UpdateRHS_MultHarm(FEMatrixType matrixType,
                                        const SBM_Vector& fup,bool SysMatUpdated) {

    LOG_TRACE(algSys) << "Updating multiharmonic RHS of matrix "<< feMatrixType.ToString(matrixType);

    if(matrixTypes_.find(matrixType) == matrixTypes_.end())
      return;

    // ensure that the RHS vector to set consists of as many
    // sub-vectors as the RHS of the system
    if( fup.GetSize() != rhs_->GetSize() ) {
      EXCEPTION( "New rhs consists of " << fup.GetSize()
                 << " sub-vectors, the RHS of the algebraic system of "
                 << rhs_->GetSize() << " entries." )
    }

    if( fup.GetEntryType() == BaseMatrix::DOUBLE ) {
      EXCEPTION(" AlgebraicSys::UpdateRHS_MultHarm called with a Double vector\n"
                " this should not happen!");
    }

    // loop over all harmonics and create a converted rhs
    if(tmpRHS_== NULL)
      tmpRHS_ = dynamic_cast<SBM_Vector*> ( GenerateVectorObject( *(sysMat_[SYSTEM]) ) );

    tmpRHS_->Init();

    if(numFcts_ != 1){ EXCEPTION("AlgebraicSys::UpdateRHS_MultHarm currently "
                                 "only implemented for ONE FunctionID");}

    // get all (blockId,index)-combinations for the current fctId
    StdVector<UInt> blockNums, indices;
    MapCompleteFctIdToIndex( 0, blockNums, indices);
    UInt size = blockNums.GetSize();

    for(UInt i = 0; i < domain->GetDriver()->GetNumFreq() ; ++i ) {
      // security check: ensure that sub-vector has the same size
      // as the block indices
      if( fup(i).GetSize() != indices.GetSize() ) {
        EXCEPTION( "Number of entries of " << i << "-th sub-vector and number "
                   "of indices do not match!");
      }

      Vector<Complex> & nRHS = dynamic_cast<Vector<Complex>&>( fup(i) );
      for( UInt j = 0; j < size; ++j ) {
        // omit entries for Dirichlet values
        if( indices[j] <= blockInfo_[blockNums[j]]->numLastFreeIndex) {
          tmpRHS_->GetPointer(i)->AddToEntry(indices[j]-1, nRHS[j] );
        }else if(!usingPenalty_){
          idbcHandler_->AddFixedToFreeRHS(matrixType, i, indices[j],rhs_,nRHS[j]);
        }
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

  template<typename T>
  void AlgebraicSys::ConstructEffectiveMatrix( const FeFctIdType fctId,
                            const std::map<FEMatrixType,T> &matFactors,
                            const bool isMultHarm) {

    LOG_DBG(algSys) << "Constructing effective system matrix for feFunction "
        << "with id " << fctId;
    if (IS_LOG_ENABLED(algSys, dbg)) {
      LOG_DBG(algSys) << "Factors are:";
      typename std::map<FEMatrixType,T>::const_iterator it = matFactors.begin();
      for( ; it != matFactors.end(); ++it ) {
        LOG_DBG(algSys) << feMatrixType.ToString(it->first) << ": " << it->second;
      }
    }

    typename std::map<FEMatrixType,T>::const_iterator it;
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
    }

    // In multiharmonic analysis, we have to adapt the mass part
    if(isMultHarm){
      SBM_Matrix *mass = sysMat_[DAMPING];
      for(UInt iRow = 0; iRow < domain->GetDriver()->GetNumFreq(); ++iRow){
        StdMatrix* sysSub = sys->GetPointer(iRow, iRow);
        StdMatrix* massSub = mass->GetPointer(iRow, iRow);
        // multiply massSub with harmonic number and add it to system matrix
        //Integer fac = (iRow - solStrat_->GetNumHarmN() < 0)? -1 : 1;
        Double fac = 1.0;
        sysSub->Add(fac, *massSub);
      }
    }else{
      for ( it = matFactors.begin(); it != matFactors.end(); it++ ) {
        if ( sysMat_[(*it).first] != NULL  && (*it).second != 0.0 ) {
          std::map<UInt, std::set<UInt> > dummyFreeSet;
          sys->Add( (*it).second, *sysMat_[(*it).first],
              dummyFreeSet, freeIndPerBlock );
        }
      }
    }
    // Also assemble the effective auxilliary system matrix for moving
    // IDBCs to the right-hand side
    idbcHandler_->BuildSystemMatrix( matFactors, freeIndPerBlock );
  }

  void AlgebraicSys::PrintKeff(){
    if(effMat_ == NULL){
      std::cout << "effMat == NULL" << std::endl;
    } else {
      //std::cout << effMat_->ToString() << std::endl;
      std::cout << effMat_->ToString() << std::endl;
    }
  }

  void AlgebraicSys::PrintMatrixPart(FEMatrixType matrixPart){
    if(sysMat_[matrixPart] == NULL){
      std::cout << "sysMat_["<<matrixPart<<"] == NULL" << std::endl;
    } else {
      //std::cout << effMat_->ToString() << std::endl;
      std::cout << sysMat_[matrixPart]->ToString() << std::endl;
    }
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

  template<typename T>
  void AlgebraicSys::SetDirichletMH( const FeFctIdType fctId,
                                   Integer eqnNr, const T& val , UInt& harmInt) {

    LOG_DBG(algSys) << "Setting Dirichlet value " << val << " for eqn " << eqnNr
                    << " of fctId " << fctId;

    if(eqnNr == 0)
      return;

    UInt blockNr, index;
    MapFctIdEqnToIndex( fctId, eqnNr, blockNr, index );

    Complex halfVal;

    if (harmInt ==0){
      halfVal = val;
    }else{
      halfVal = val / 2.0;
    }
    UInt N = domain->GetDriver()->GetNumFreq()/2;
    // Delegate work to IDBC handler
    idbcHandler_->SetIDBC( N+harmInt, index,  halfVal );
    idbcHandler_->SetIDBC( N-harmInt, index,  std::conj(halfVal) );
  }

  void AlgebraicSys::BuildInDirichlet() {

    LOG_DBG(algSys) << "Incorporating Dirichlet values into system matrix";

    // If necessary modify matrix diagonal for penalty approach
    if ( assembleDirichletToSysMat_ == true ) {
      idbcHandler_->AdaptSystemMatrix( *(sysMat_[SYSTEM]) );
      assembleDirichletToSysMat_ = false;
    }
  }

// Remove zero non-zero elements from system matrix
  template <typename T>
  void AlgebraicSys::GetRidOfZeros(Double tol)
  { 
    PtrParamNode infoGetRidOfZeros = myInfo_->Get("GetRidOfZeros", ParamNode::APPEND);
    Timer timeGetRidOfZeros;
    timeGetRidOfZeros.Start();
    
    // TODO
    // This routine works as it is, but it might be better to introduce to combine it with the function setting IDBCs.
    // Since assembling IDBCs might change some entries, we should check beforehand which entries should not be touched.
    // Otherwise, we have to reiterate over all entries which might take some time. By redefining a pattern compliant 
    // with the IDBC handler, we could store and reuse the location of the non-zero entries. This could potentially save
    // us the call to the solver where we have to set a new matrix pattern. This might not be restricted to IDBCs and 
    // should be thoroughly tested.
    // effMat_ eventually gets passed to the solver, but this is just a weak copy of sysMat_[SYSTEM], hence we modify the system matrix itself
    
    // create a deep copy and work with the copy
    SBM_Matrix * actMat = new SBM_Matrix(*(sysMat_[SYSTEM]));
    UInt nRows = actMat->GetNumRows();
    UInt nCols = actMat->GetNumCols();

    // store the old sysMat so that we can access it later on
    SetSysMatBackup(actMat);

    LOG_DBG(algSys) << "\tGet rid of Zeros in sys matrix. Number of Rows of SBMMatrix: " << nRows << " Number of Coloumns: " << nCols;
    for ( UInt row = 0; row < nRows; row++ ){
      for ( UInt col = 0; col < nRows; col++ ){
        // get the current sub matrix
        LOG_DBG(algSys) << "\tEditing SBM sub matrix (" << row << "," << col << ")";
        StdMatrix * stdMat = actMat->GetPointer(row,col);

        UInt* RowVec;
        UInt* ColVec;
        UInt NumRows;
        UInt newNnz = 0;
        T* DataVec;
        // check for SCRS/CRS
        switch(stdMat->GetStorageType())
        {
        case BaseMatrix::SPARSE_SYM:
        {
          SCRS_Matrix<T> &tmpSym = dynamic_cast<SCRS_Matrix<T>&>(*stdMat);
          RowVec = tmpSym.GetRowPointer();
          ColVec = tmpSym.GetColPointer();
          DataVec = tmpSym.GetDataPointer();
          NumRows = tmpSym.GetNumRows();
        }
          break;
        case BaseMatrix::SPARSE_NONSYM:
        {
          CRS_Matrix<T> &tmpNonSym = dynamic_cast<CRS_Matrix<T>&>(*stdMat);
          RowVec = tmpNonSym.GetRowPointer();
          ColVec = tmpNonSym.GetColPointer();
          DataVec = tmpNonSym.GetDataPointer();
          NumRows = tmpNonSym.GetNumRows();
        }
          break;
        
        default:
          EXCEPTION("AlgebraicSys::GetRidOfZeros: Could not find suitable matrix conversion to get rid of zeros, skipping this matrix");
          break;
        }
        // now loop over the entries of the sub matrix and check for zero entries --> compute new NNZ
        /* for (UInt i = 0; i < size; i++) {
          //if (DataVec[i] != 0.0) newNnz++;
          if (abs(DataVec[i]) >= tol) newNnz++;
        } */
        // above you can find the old implementation which can cause probmes because sometimes - for
        // whatever reason - the dataVec contains quite a bit more entries than those which are
        // actually accessed via the row- and column vector. Hence, we also use those vectors to
        // compute our NNZ
        for (UInt i = 0; i < NumRows; i++) {
          for (UInt j = RowVec[i]; j < RowVec[i+1]; j++) {
            if (abs(DataVec[j]) >= tol) newNnz++;
          }
        }
        UInt oldNnz = stdMat->GetNnz();
        // debug info
        //std::cout << "Old NNZ:" << oldNnz << std::endl << "New NNZ: " << newNnz << std::endl;
        // log intermediate results
        LOG_DBG(algSys) << "\tNNZ according to initial guess "<< oldNnz;
        LOG_DBG(algSys) << "\tNNZ actually counted in the overall vector" << newNnz;
        // check if we would actually gain anything from this procedure
        if ( newNnz < oldNnz ) {
          LOG_DBG(algSys) << "\tProceeding to create new matrix";
          // we have less non-zero entries, proceed to construct the new matrix
          StdVector<T> DataVecNew;
          StdVector<UInt> ColVecNew;
          StdVector<UInt> RowVecNew;
          StdVector<UInt> DiagVecNew;
          UInt k{0};
          RowVecNew.push_back(RowVec[0]);
          UInt matCount = 0;
          UInt nnzMatCount = 0;
          for (UInt i = 0; i < NumRows; i++) {
            // debug
            //std::cout << "RowVec[i]: " << RowVec[i] << "; RowVec[i+1]: " << RowVec[i+1] << std::endl;
            // Assume that there is no diagonal entry (and correct
            // it below, if we find one)
            DiagVecNew.push_back(0);
            for (UInt j = RowVec[i]; j < RowVec[i+1]; j++) {
              matCount++;
              //if (DataVec[j] != 0.0) {
              // add non-zero entries as well as diagonal entries (regardless of their actual value)
              if (abs(DataVec[j]) >= tol || i==ColVec[j]) {
                DataVecNew.push_back(DataVec[j]);
                ColVecNew.push_back(ColVec[j]);
                k++;
                nnzMatCount++;
              }
              // replace zero with actual diagonal value
              if (i==ColVec[j]) {
                DiagVecNew[i] = j;
              }
            }
            RowVecNew.push_back(k);
          }
          // debug
          /* for (UInt i = 0; i < RowVecNew.size(); i++) {
            std::cout << "RowVecNew[" << i << "] = " << RowVecNew[i] << std::endl;
          } */
          //std::cout << "Looped over " << matCount << " elements" << std::endl;
          //std::cout << "Found " << nnzMatCount << " non-zero elements (including potential zeros on the diagonal)" << std::endl;
          
          // log
          LOG_DBG(algSys) << "\tLooped over " << matCount << " elements";
          LOG_DBG(algSys) << "\tNNZ actually counted in the used matrix (including potential zeros on the diagonal)" << nnzMatCount;
          // debug
          //stdMat->Export("oldMat.mtx", BaseMatrix::MATRIX_MARKET);
          // delete the old effMat_ since otherwise we'll have some not so nice heap corruption
          delete effMat_;
          effMat_ = nullptr;
          // overwrite old matrix with new (resized) empty matrix (the last bool forces an overwrite)
          LOG_DBG(algSys) << "\tOverwriting old sub matrix";
          // generate new empty SBM-Matrix
          actMat->SetSubMatrix(row, col, stdMat->GetEntryType(), stdMat->GetStorageType(),
                                stdMat->GetNumRows(), stdMat->GetNumCols(), nnzMatCount, true, true);
          stdMat = actMat->GetPointer(row,col);
          
          // now set the sparsity pattern data for the newly created submatrix
          switch(stdMat->GetStorageType())
          {
          case BaseMatrix::SPARSE_SYM:
          {
            SCRS_Matrix<T> &tmpSym = dynamic_cast<SCRS_Matrix<T>&>(*stdMat);
            tmpSym.SetSparsityPatternData(RowVecNew, ColVecNew, DataVecNew);
            // here we also have to correct the wrongly initialized numEntries variable
            tmpSym.SetNumEntries(nnzMatCount);
          }
            break;
          case BaseMatrix::SPARSE_NONSYM:
          {
            CRS_Matrix<T> &tmpNonSym = dynamic_cast<CRS_Matrix<T>&>(*stdMat);
            tmpNonSym.SetSparsityPatternData(RowVecNew, ColVecNew, DataVecNew);
          }
            break;
          default:
            EXCEPTION("This should not be possible")
            break;
          }
          // delete the old sysMat_[SYSTEM]
          delete sysMat_[SYSTEM];
          // now we have a full and cleaned copy of sysMat_[SYSTEM] (actMat), we just need to pass it back
          sysMat_[SYSTEM] = actMat;
          LOG_DBG(algSys) << "\tCreated new sysMat_ entry";
          // debug
          //actMat->Export("actMat", BaseMatrix::OutputFormat::MATRIX_MARKET, "actMat");
          //sysMat_[SYSTEM]->Export("sysMat", BaseMatrix::OutputFormat::MATRIX_MARKET, "sysMat");
          // rebuild effMat_
          // TODO: we might have to do this outside of the loop, should be only important for harmonic analysis which is not working anyways at the moment
          UInt nB = (isMultHarm_)? domain->GetDriver()->GetNumFreq() : numBlocks_;
          for ( UInt k = 0; k < nB; k++ ) {
            if (statCond_) {
              delete effMat_;
              effMat_ = new SBM_Matrix( *sysMat_[SYSTEM], numBlocks_-1,
                                        numBlocks_-1 );
            } else {
              delete effMat_;
              effMat_ = new SBM_Matrix( *sysMat_[SYSTEM], nB, nB );
            }
          }
          LOG_DBG(algSys) << "\tRebuilt effMat_";
          LOG_DBG(algSys) << "\tSuccessfully reduced complexity of the system by eliminating unnecessary zeros";
          infoGetRidOfZeros->Get("zeros_removed", ParamNode::APPEND)->SetValue(oldNnz - newNnz);
        } else {
          // no zeros found, we skip everything
          WARN("AlgebraicSys::GetRidOfZeros: No zeros to remove with defined tolerance were found")
          UInt intZero = 0;
          infoGetRidOfZeros->Get("zeros_removed", ParamNode::APPEND)->SetValue(intZero);
          LOG_DBG(algSys) << "No zero entries found, skipping this matrix";

          // free the matrix, otherwise the PatternPool destructor will throw an error
          delete actMat;
        }
      }
    }
    timeGetRidOfZeros.Stop();
    infoGetRidOfZeros->Get("cpu", ParamNode::APPEND)->SetValue(timeGetRidOfZeros.GetCPUTime());
    infoGetRidOfZeros->Get("wall", ParamNode::APPEND)->SetValue(timeGetRidOfZeros.GetWallTime());
  }

  template void AlgebraicSys::GetRidOfZeros<Double>(Double tol);
  template void AlgebraicSys::GetRidOfZeros<Complex>(Double tol);

  void AlgebraicSys::SetSysMatBackup( SBM_Matrix* actMat ){
    delete sysMat_[BACKUP];
    sysMat_[BACKUP] = new SBM_Matrix(*actMat);
  }

  void AlgebraicSys::RestoreSystemMatrixFromBackup(){
    delete effMat_;
    effMat_ = nullptr;
    delete sysMat_[SYSTEM];
    sysMat_[SYSTEM] = new SBM_Matrix(*sysMat_[BACKUP]);

    UInt nB = (isMultHarm_)? domain->GetDriver()->GetNumFreq() : numBlocks_;
    for ( UInt k = 0; k < nB; k++ ) {
      if (statCond_) {
        delete effMat_;
        effMat_ = new SBM_Matrix( *sysMat_[SYSTEM], numBlocks_-1,
                                  numBlocks_-1 );
      } else {
        delete effMat_;
        effMat_ = new SBM_Matrix( *sysMat_[SYSTEM], nB, nB );
      }
    }
    LOG_DBG(algSys) << "\tRestored sysMat_[SYSTEM] and effMat_";
  }

  void AlgebraicSys::ClearIDBCFromSolutionVal( SBM_Vector& solVec ){
    // resize solVec to match number of functions
    solVec.Resize( numFcts_);

    // loop over all feFctIDs
    for(UInt i = 0; i < numFcts_; ++i ) {

      // call specialized GetSolutionVal method
      ClearIDBCFromSolutionVal(solVec(i), i);
    }
  }

  void AlgebraicSys::ClearIDBCFromSolutionVal( SingleVector& ptSol,const FeFctIdType fctId){
    LOG_DBG(algSys) << "Clearing IDBC nodes from solution values of fct " << fctId;

    // get all (blockId,index)-combinations for the current fctId
    StdVector<UInt> blockNums, indices;
    MapCompleteFctIdToIndex( fctId, blockNums, indices);
    UInt size = blockNums.GetSize();

    if( ptSol.GetEntryType() == BaseMatrix::DOUBLE ) {
      Vector<Double> & retVec = dynamic_cast<Vector<Double>&>( ptSol );
      Double entry = 0.0;

      // #pragma omp parallel for private (entry)
      for( UInt i = 0; i < size; ++i )
      {
        // if index number is larger the lastFree dof, set value in vector to 0
        // i.e. overwrite IDBC values; otherwise, do nothing (normal nodes shall be preserved!)
        if( indices[i] > blockInfo_[blockNums[i]]->numLastFreeIndex )
        {
          retVec[i] = entry;
        }
      }
    }
    else
    {
      Vector<Complex> & retVec = dynamic_cast<Vector<Complex>&>( ptSol );
      Complex entry = 0.0;

      // #pragma omp parallel for private (entry)
      for( UInt i = 0; i < size; ++i )
      {
        // if index number is larger the lastFree dof, set value in vector to 0
        // i.e. overwrite IDBC values; otherwise, do nothing (normal nodes shall be preserved!)
        if( indices[i] > blockInfo_[blockNums[i]]->numLastFreeIndex )
        {
          retVec[i] = entry;
        }
      }
    }
  }

  void AlgebraicSys::GetSolutionVal( SingleVector& ptSol,
                                     const FeFctIdType fctId,
                                     bool setIDBC, bool deltaIDBC) {

    LOG_DBG(algSys) << "Getting solution values of fct " << fctId;

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
            idbcHandler_->GetIDBC(blockNums[i],  indices[i], entry, deltaIDBC);
          else
            entry = 0.0;
        }
        else
        {
          sol_->GetPointer(blockNums[i])->GetEntry(indices[i]-1,entry);
        }
        retVec[i] = entry;
        LOG_DBG2(algSys) << "retVec[" << i << "]= " << entry;
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
            idbcHandler_->GetIDBC(blockNums[i],  indices[i], entry, deltaIDBC);
          else
            entry = 0.0;
        } else {
          sol_->GetPointer(blockNums[i])->GetEntry(indices[i]-1,entry);
        }
        retVec[i] = entry;
        LOG_DBG2(algSys) << "retVec[" << i << "]= " << entry;
      }
    }
  }


  void AlgebraicSys::GetSolutionVal( SingleVector& ptSol,
                                     const UInt& block,
                                     bool setIDBC,
                                     bool deltaIDBC,
                                     const bool ident) {

    LOG_TRACE(algSys) << "Getting multiharmonic solution values ";

    //idbcHandler_->PrintIDBCvec();

    StdVector<UInt> blockNums, indices;

    MapCompleteFctIdToIndex( 0, blockNums, indices);
    UInt size = blockNums.GetSize();
    ptSol.Resize(size);
    ptSol.Init();

    if( ptSol.GetEntryType() == BaseMatrix::DOUBLE ) {
      EXCEPTION("AlgebraicSys::GetSolutionVal Double-Version for multiharmonic not implemented");
    }else{
      Vector<Complex> & retVec = dynamic_cast<Vector<Complex>&>( ptSol );
      Complex entry = 0.0;
      for( UInt i = 0; i < size; ++i ) {
        // if index number is larger the lastFree dof, insert Dirichlet value
        // (just if setIDBC is true!)
        if( indices[i] > blockInfo_[0]->numLastFreeIndex ) {
          if ( setIDBC){
            idbcHandler_->GetIDBC(block,  indices[i], entry, deltaIDBC);
          }else{
            entry = 0.0;
          }
        } else {
          sol_->GetPointer(block)->GetEntry(indices[i]-1,entry);
        }
        retVec[i] = entry;
      }
    }
  }

  void AlgebraicSys::GetSolutionVal( SBM_Vector& solVec, bool setIDBC, bool deltaIDBC ) {
    solVec.Resize( numFcts_);
    // loop over all feFctIDs
    for(UInt i = 0; i < numFcts_; ++i )
      GetSolutionVal(solVec(i), i, setIDBC, deltaIDBC);
  }


  void AlgebraicSys::GetSolutionVal( const UInt& h,  SBM_Vector& solVec, bool setIDBC, bool deltaIDBC ) {
    if(numFcts_ != 1){
      EXCEPTION("AlgebraicSys::GetSolutionVal This shouldn't happen");
    }

    solVec.Resize( 1 );
    // Get the correct block-vector
    // call specialized GetSolutionVal method, the boolean
    // has no effect, it's just an identifier to call the correct method
    GetSolutionVal(solVec(0), h, setIDBC, deltaIDBC, true);
  }


  void AlgebraicSys::GetFullMultiHarmSolutionVal(SBM_Vector& solVec, bool setIDBC, bool deltaIDBC ) {
    solVec.Resize( domain->GetDriver()->GetNumFreq() );
    // solVec gets initialized in GetSolutionVal method

    // loop over all block vector and call specialized GetRHSVal method, the boolean
    // has no effect, it's just an identifier to call the correct method
    for( UInt i = 0; i < solVec.GetSize(); ++i){
      GetSolutionVal(solVec(i), i, setIDBC, deltaIDBC, true);
    }
  }


  void AlgebraicSys::GetRHSVal( SingleVector &ptRhs,
                                const FeFctIdType fctId  ) {

    LOG_DBG(algSys) << "Getting RHSvalue of fct " << fctId;

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
        LOG_DBG2(algSys) << "i= " << i << " retVec[i]= " << entry;
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
      LOG_DBG(algSys) << "retVec " << retVec.ToString();
    }
  }



  void AlgebraicSys::GetRHSVal( SingleVector &ptRhs,
                                const UInt& blockVec,
                                const bool ident) {

    LOG_TRACE(algSys) << "Getting multiharmonic RHSvalue";

    StdVector<UInt> blockNums, indices;
    MapCompleteFctIdToIndex( NO_FCT_ID, blockNums, indices);
    UInt size = blockNums.GetSize();
    ptRhs.Resize(size);
    ptRhs.Init();

    if( ptRhs.GetEntryType() == BaseMatrix::DOUBLE ) {
      EXCEPTION("AlgebraicSys::GetRHSVal This method shall only be called in multiharmonic analysis!!");
    } else {
      Vector<Complex> & retVec = dynamic_cast<Vector<Complex>&>( ptRhs );
      Complex entry = 0.0;

      for( UInt i = 0; i < size; ++i ) {
        // if index number is larger the lastFree dof, insert 0
        // Remember: in multiharmonic analysis only one blockInfo_ !
        if( indices[i] > blockInfo_[0]->numLastFreeIndex) {
          entry = 0.0;
        } else {
          rhs_->GetPointer(blockVec)->GetEntry(indices[i]-1,entry);
        }
        retVec[i] = entry;
      }
    }
  }

  void AlgebraicSys::GetRHSVal( SBM_Vector& rhsVec ) {
    // resize rhsVec to match number of functions
    rhsVec.Resize( numFcts_);
    // loop over all feFctIDs
    for(UInt i = 0; i < numFcts_; ++i ) GetRHSVal(rhsVec(i), i);
  }



  void AlgebraicSys::GetRHSVal( const UInt& h, SBM_Vector& rhsVec ) {
    rhsVec.Resize( numFcts_ );

    // Get the correct block-vector
    // call specialized GetRHSVal method, the boolean
    // has no effect, it's just an identifier to call the correct method
    GetRHSVal(rhsVec(0), h, true);
  }


  void AlgebraicSys::GetFullMultiHarmRHSVal(SBM_Vector& rhsVec ) {
    rhsVec.Resize( domain->GetDriver()->GetNumFreq() );
    rhsVec.Init();

    // loop over all block vector and call specialized GetRHSVal method, the boolean
    // has no effect, it's just an identifier to call the correct method
    for( UInt i = 0; i < rhsVec.GetSize(); ++i){
      GetRHSVal(rhsVec(i), i, true);
    }
  }


  SBM_Matrix* AlgebraicSys::GenerateSBM_Matrix( FEMatrixType matType,
                                                BaseMatrix::EntryType entryType,
                                                bool sharePattern ) {

    LOG_DBG(algSys) << "Generating SBMMatrix of size " << numBlocks_
        << " for matrix type " << feMatrixType.ToString(matType);

    // STEP 1: Generate empty SBM_Matrix
    SBM_Matrix *retMat = NULL;


    /* Strategy for multiharmonic system:
     * Generate (N+1 x N+1) or (N+2 x N+2) SBM matrix...driver privides the correct size
     * in the variable numFreq_  but only populate the
     * nonzero blocks, the sbm class should recognize that ...
     */
    if( isMultHarm_ ){
      UInt s =  domain->GetDriver()->GetNumFreq(); //size

      // Multiharmonic system matrix has size [numFreq, numFreq]
      bool sysMatSym = false; //system matrix is not symmetric
      retMat = new SBM_Matrix( s, s, sysMatSym );
    }else{
      retMat = new SBM_Matrix( numBlocks_, numBlocks_, sbmSymm_ );
    }


    if ( retMat == NULL ) {
      EXCEPTION( "SBM_System::GenerateSBM_Matrix: "
          << "This is the end my friend!\n"
          << "Generation of empty SBM_Matrix failed!" );
    }


    // STEP 2: Populate with sub-matrices
    std::set<SubMatrixID,SortSubMatrixID>::iterator sIt;
    BaseGraph *graph = NULL;

    // again differentiate between multiharmonic and classic case
    if( isMultHarm_ ){
      //if( feSubMatricesByBlocks_[matType].size() > 1 ){
      //  EXCEPTION("AlgebraicSys::GenerateSBM_Matrix, only one block allowed")
      //}

      // we have only one graph, since all blocks consist of the same function
      // Determine number of matrix rows and columns
      UInt nrows = blockInfo_[0]->numLastFreeIndex;
      UInt ncols = blockInfo_[0]->numLastFreeIndex;

      //UInt N = solStrat_->GetNumHarmN();
      UInt M = solStrat_->GetNumHarmM();
      UInt a = (domain->GetDriver()->IsFullSystem())? (M+1) : ((M-1)/2 + 1);
      for( UInt sbmRow = 0; sbmRow < domain->GetDriver()->GetNumFreq(); ++sbmRow ) {
        for( UInt sbmCol = sbmRow ; sbmCol < sbmRow + a; ++sbmCol ) {
          if( sbmCol < domain->GetDriver()->GetNumFreq()){
            graph = graphManager_->GetGraph( sbmRow, sbmCol );
            // we only allow nonsymmetric storage scheme
            BaseMatrix::StorageType sT = BaseMatrix::SPARSE_NONSYM;
            LOG_DBG(algSys) << "storage Type of matrix (" << sbmRow +1
                << ", " << sbmCol+1 << ") is "
                << BaseMatrix::storageType.ToString(sT);
            retMat->SetSubMatrix ( sbmRow, sbmCol, entryType, sT, nrows, ncols, graph->GetNNE() );
            LOG_DBG(algSys)<< "Inserting a SubMatrix at: " << "[" << sbmRow << "," << sbmCol << "]";
            LOG_DBG(algSys)<< "SubMatrix has size of: " << "[" <<nrows<< "," << ncols << "]";

            // check, if matrix pattern can be shared and
            // obtain matrix graph
            if( sharePattern ) {
              SubMatrixID id;
              id.rowInd = sbmRow;
              id.colInd = sbmCol;
              LOG_DBG(algSys) << "\tSharing pattern";
              PatternIdType patternID = sbmPatternIds_[id];
              if( sbmPatternIds_[id] != NO_PATTERN_ID ) {
                LOG_DBG(algSys) << "\tObtaining pattern '" << patternID << "' from pool";
                (*retMat)( sbmRow, sbmCol ).SetSparsityPattern( patternPool_, patternID );
              } else {
                (*retMat)( sbmRow, sbmCol ).SetSparsityPattern( *graph );
                sbmPatternIds_[id] = (*retMat)( sbmRow, sbmCol ).TransferPatternToPool( patternPool_ );
                LOG_DBG(algSys) << "\tPutting pattern '" << sbmPatternIds_[id] << "' to pool";
              }
            } else {
              LOG_DBG(algSys) << "\tUsing no shared sparsity pattern";
              // Set sparsity pattern of sub-matrix
              (*retMat)( sbmRow, sbmCol ).SetSparsityPattern( *graph );
            }

            // also set the transposed...if it's not the diagonal block
            if(sbmRow != sbmCol){
              LOG_DBG(algSys) << "storage Type of matrix (" << sbmCol +1
                  << ", " << sbmRow+1 << ") is "
                  << BaseMatrix::storageType.ToString(sT);
              retMat->SetSubMatrix ( sbmCol, sbmRow, entryType, sT, nrows, ncols, graph->GetNNE() );

              if( sharePattern ) {
                SubMatrixID id;
                id.rowInd = sbmCol;
                id.colInd = sbmRow;
                LOG_DBG(algSys) << "\tSharing pattern";
                PatternIdType patternID = sbmPatternIds_[id];
                if( sbmPatternIds_[id] != NO_PATTERN_ID ) {
                  LOG_DBG(algSys) << "\tObtaining pattern '" << patternID << "' from pool";
                  (*retMat)( sbmCol, sbmRow ).SetSparsityPattern( patternPool_, patternID );
                } else {
                  (*retMat)( sbmCol, sbmRow ).SetSparsityPattern( *graph );
                  sbmPatternIds_[id] = (*retMat)( sbmCol, sbmRow ).TransferPatternToPool( patternPool_ );
                  LOG_DBG(algSys) << "\tPutting pattern '" << sbmPatternIds_[id] << "' to pool";
                }
              } else {
                LOG_DBG(algSys) << "\tUsing no shared sparsity pattern";
                // Set sparsity pattern of sub-matrix
                (*retMat)( sbmCol, sbmRow ).SetSparsityPattern( *graph );
              }
            }

          } // endif sbmCol < size
        } // loop over sbmCols
      } // loop over sbmRows



    }else{

      for ( sIt = feSubMatricesByBlocks_[matType].begin();
          sIt != feSubMatricesByBlocks_[matType].end(); sIt++ ) {

        // Determine row / col
        UInt sbmRow = (*sIt).rowInd;
        UInt sbmCol = (*sIt).colInd;

        // Determine number of matrix rows and columns
        UInt nrows = blockInfo_[sbmRow]->numLastFreeIndex;
        UInt ncols = blockInfo_[sbmCol]->numLastFreeIndex;
        LOG_DBG(algSys) << "The shape of following matrix is: " << "[" << nrows << "," << ncols <<"]";

        // Check for necessity of generation
        if ( sbmRow <= sbmCol || sbmSymm_ == false ) {

          graph = graphManager_->GetGraph( sbmRow, sbmCol );
          //sbmSymm_ = false;
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
    }
    return retMat;
  }


  BaseOrdering::ReorderingType AlgebraicSys::ReorderingDefault(PtrParamNode pn, bool& can_change)
  {
    LOG_DBG(algSys) << "RD: xml reordering: " << (pn->Has("reordering") ? pn->Get("reordering")->As<string>() : "none") << " RD=" << REORDERING_DEFAULT;

    // we have only no default when we set im xml other than "default"
    if(pn->Has("reordering") && pn->Get("reordering")->As<string>() != "default")
    {
      can_change = false; // no, we know what we want
      return BaseOrdering::reorderingType.Parse(pn->Get("reordering")->As<string>());
    }

    // we have default, select from compile settings but allow for change
    can_change = true; // holds for all comming returns

    // do we use CFS_REORDERING with a manual setting?
    if(string(REORDERING_DEFAULT).size() > 0 && string(REORDERING_DEFAULT) != "default")
      return  BaseOrdering::reorderingType.Parse(REORDERING_DEFAULT);

    // use METIS if present - might be changed anyway depending on solver
    // what we set by compile option
    #ifdef USE_METIS
      return BaseOrdering::METIS;
    #else
      return BaseOrdering::SLOAN;
    #endif
  }

  void AlgebraicSys::CheckConsistency()
  {
    // First check, if we have a true SBM system.
    // consisting of more than one SBM-Block
    if( (numBlocks_ == 1 || (numBlocks_ == 2 && statCond_) ) && !isMultHarm_ ) {

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
      PtrParamNode eSolverList = myParam_->Get("eigenSolverList", ParamNode::INSERT);
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
      PtrParamNode solverList = myParam_->Get("solverList", ParamNode::INSERT);
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

      // check if a solver is specified
      if(!solverNode)
      {
        // no solver set -> use default direct solver in order of availability: pardiso, cholmod, directLDL
#ifdef USE_PARDISO
        st = BaseSolver::PARDISO_SOLVER;
#elif USE_CHOLMOD
        st = BaseSolver::CHOLMOD_SOLVER;
#else
        if(storType == BaseMatrix::SPARSE_SYM)
        {
          st = BaseSolver::LDL_SOLVER; // symmetric case
        } else {
          st = BaseSolver::LU_SOLVER; // unsymmetric case
        }
#endif
        solverList->Get(BaseSolver::solverType.ToString(st),ParamNode::INSERT)->Get("id",ParamNode::INSERT)->SetValue(solverId);
      }
      else
      {
        //  solver set -> check for compatibility with matrix

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
      PtrParamNode precondList = myParam_->Get("precondList", ParamNode::INSERT);
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
      // this is the standard reordering from compile settings. We generally use it
      // - if not the user selects something differen
      // - if not we keeep "default" in xml and we have solver which does it by itself (or does not profit from it)
      bool canChangeReordering = false; // overwritten in ReorderingDefault()
      BaseOrdering::ReorderingType ot = ReorderingDefault(matNode, canChangeReordering);

      LOG_DBG(algSys) << "CS: pre handling: " << BaseOrdering::reorderingType.ToString(ot) << " change=" << canChangeReordering;

      // ot can now be be sload, metis from default or noreodering from matNode. No idea how minimumDegree or nestedDissection can be selected?!

      // a) for our own direct solvers ot stays metis or sloan in the default case.
      //    I the non-default case, it can also be noreordering when selected by the user

      // b) pardiso and most external solvers need no reordering or have their own
      // if we are a solver which does not need reordering and if we are default, we switch reordering off
      if((st == BaseSolver::PARDISO_SOLVER ||
          st == BaseSolver::CHOLMOD ||
          st == BaseSolver::UMFPACK ||
          st == BaseSolver::LIS ||
          st == BaseSolver::SUPERLU ||
          st == BaseSolver::PETSC ||
          st == BaseSolver::GINKGO) && (canChangeReordering == true) )
      {
        ot = BaseOrdering::NOREORDERING;
      }
      // c) ilu-based preconditioners prefer reordering: stay with the default ot

      // in the end store back the reordering type
      matNode->Get("reordering",ParamNode::INSERT)->SetValue(BaseOrdering::reorderingType.ToString(ot));
    }

    if( isMultHarm_ ){
        // =========================================================================
        //  True SBM Case, e.g. for multiharmonic analysis
        // =========================================================================
        WARN("The implementation of this section is not yet finished");

        // --------------------------
        //  Check Symmetry of Matrix
        // --------------------------
        PtrParamNode matNode = solStrat_->GetMatrixNode(0);
        BaseMatrix::StorageType storType = BaseMatrix::NOSTORAGETYPE;

        // we only allow nonsymmetric storage format
        std::string storageString = "sparseNonSym";
        if( matNode->Has("storage")) {
          if( matNode->Get("storage")->As<std::string>() != "sparseNonSym" ){
           EXCEPTION(" You probably perform a multiharmonic analysis, therefore please set nonsymmetric storage type! ");
          }
        }

        bool canChangeMatFormat = false;

        matNode->GetValue("storage",storageString, ParamNode::INSERT);
        storType = BaseMatrix::storageType.Parse(storageString);


        // -----------------------------------------------
        //  Check of Eigenvalue Solver not yet implemented
        // -----------------------------------------------

        // ------------------------------------------------
        //  Check Solver
        // ------------------------------------------------
        std::string solverId = solStrat_->GetSolverId();
        PtrParamNode solverList = myParam_->Get("solverList", ParamNode::INSERT);
        ParamNodeList sNodes =  solverList->GetChildren();
        PtrParamNode solverNode;
        for( UInt i = 0; i < sNodes.GetSize(); ++i ){
          if( sNodes[i]->Get("id")->As<std::string>() == solverId ){
            solverNode = sNodes[i];
          }
        }
        BaseSolver::SolverType st;
        // set for allowed matrix types of the solver
        std::set<BaseMatrix::StorageType> solverStorTypes;
        if( !solverNode ) {
          // ---------------------------------------------------
          //  no solver set -> use default direct
          // ---------------------------------------------------
          st = BaseSolver::PARDISO_SOLVER;
          solverList->Get("pardiso",ParamNode::INSERT)->
            Get("id",ParamNode::INSERT)->SetValue(solverId);
        }else{
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
            EXCEPTION("Solver '" << solverNode->GetName()
                      << "' can not operate on matrix with storage type '"
                      << storageString << "'. \nChange format to '"
                      << BaseMatrix::storageType.ToString(storType)
                      << "'.");
          }
        }

        // -------------------------------------------------------
        //  Check Precond
        // -------------------------------------------------------
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
          // -------------------------------------------------------
          //  no precond set -> use default ID
          // -------------------------------------------------------
          pt = BasePrecond::ID;
          precondList->Get("Id",ParamNode::INSERT)->
              Get("id",ParamNode::INSERT)->SetValue(precondId);
        }else{
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
            bool isCompatibleWithSolver = (solverStorTypes.size() == 0 || solverStorTypes.find(storType) != solverStorTypes.end() );

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
          //EXCEPTION("Preconditioning for true SBM case not yet implemented!");
        }

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

        // --------------------------------------------------------
        //  Check for shared pattern
        // --------------------------------------------------------
        if( st == BaseSolver::DIAGSOLVER ) {
          sharedPatternPossible_ = false;
        }

        // ---------------
        //  Check Reordering
        // ---------------
        bool canChangeReordering = false; // overwritten in ReorderingDefault()
        BaseOrdering::ReorderingType ot = ReorderingDefault(matNode, canChangeReordering);

        // a) for our own direct solvers we use the default or selected reordering

        // b) pardiso and most external solvers need no reordering or have their own - switch it off when we have default
        if(canChangeReordering == true &&
           (st == BaseSolver::PARDISO_SOLVER || st == BaseSolver::UMFPACK || st == BaseSolver::SUPERLU /* || st == BaseSolver::LIS */))
        {
          ot = BaseOrdering::NOREORDERING;
        }

        // c) ilu-based preconditioners prefer reordering we stick to the default selection

        // in the end store back the reordering type
        matNode->Get("reordering",ParamNode::INSERT)->SetValue(BaseOrdering::reorderingType.ToString(ot));
    }

  }

  void AlgebraicSys::PrintFeMatrixInfo( ) {

    LOG_DBG(algSys) << "Print matrix information";

    PtrParamNode setupNode = myInfo_->Get("setup");

    // Print overview of defined matrices
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
        if( smId.rowInd == smId.colInd && !isMultHarm_) {
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

    LOG_DBG(algSys) << "Print registration info";

    PtrParamNode setupNode = myInfo_->Get("setup");

    // Print overview of feFunctions
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
        blockListNode->Get("block",ParamNode::APPEND)->Get("id")->SetValue(*blockIt);
      }
    }

    // add totalSize, totalNumDirichlet
    fctListNode->Get("totalNumEqns")->SetValue(totalSize);
    fctListNode->Get("totalNumFreeEqns")->SetValue(totalNumFreeEqns);
    fctListNode->Get("totalNumDirichlet")->SetValue(totalNumDirichlet);

    // Print overview of blocks
    PtrParamNode blockListNode = setupNode->Get("sbmBlocks");

    for( UInt i = 0; i < numBlocks_; ++i ) {
      PtrParamNode blockNode = blockListNode->Get("block",ParamNode::APPEND);
      blockNode->Get("id")->SetValue(i);
      blockNode->Get("size")->SetValue(blockInfo_[i]->size);
      blockNode->Get("lastFreeIndex")->SetValue(blockInfo_[i]->numLastFreeIndex);
    }
  }


  void AlgebraicSys::SetGeomIndexMap(const StdVector< Vector<Double> >& indGeomMap,
                                     const UInt& dim){
    dim_ = dim;
    edge_ = false;
    //loop over all indices
    geomInd_.Resize(indGeomMap.GetSize());
    for( UInt ind = 0; ind < indGeomMap.GetSize(); ++ind ) {
      geomInd_[ind] = indGeomMap[ind];
    }
  }

  void AlgebraicSys::SetEdgeIndexMap(std::unordered_map<Integer, Double>& lengths,
                                     std::unordered_map< Integer, StdVector<Integer> >& eNodes){
    dim_ = 1;
    edge_ = true;

    //loop over all indices
    std::unordered_map< Integer , Double >::const_iterator lIt = lengths.begin();
    std::unordered_map< Integer , StdVector<Integer> >::const_iterator eIt = eNodes.begin();
    while(lIt != lengths.end() ){
      StdVector<Integer> e = eIt->second;
      Double l = lIt->second;
      geomIndEdge_[lIt->first].eNodes = e;
      geomIndEdge_[lIt->first].length = l;
      lIt++;
      eIt++;
    }
  }

  void AlgebraicSys::BuildAMGAuxMatrix(){
    std::cout << "++ Assembling auxiliary matrix for AMG-preconditioner "<<std::endl;
    if(!edge_){
      BuildAMGLagrangeAuxMatrix();
    }else{
      BuildAMGEdgeAuxMatrix();
    }
  }

  void AlgebraicSys::BuildAMGEdgeAuxMatrix(){
    amgType_ = EDGE;

    /************ Perform some checks on the system matrix ***********/
    // only singlefield...only one submatrix
    BaseMatrix& b = (*effMat_)(0,0);

    // Check if we have a StdMatrix
    if ( b.GetStructureType() != BaseMatrix::SPARSE_MATRIX ) {
      EXCEPTION( "AlgebraicSys::BuildAMGEdgeAuxMatrix(): "
          << " AMG cannot deal with matrices other than StdMatrix" );
    }

    StdMatrix& stdMat = dynamic_cast<StdMatrix&>(b);

    // Now test the storage layout
    BaseMatrix::StorageType sType = stdMat.GetStorageType();
    if ( sType != BaseMatrix::SPARSE_NONSYM ) {
      EXCEPTION( "AlgebraicSys::BuildAMGEdgeAuxMatrix(): AMG requires the system matrix"
          << " to be a CRS_Matrix i.e. sparseNonSym. The system matrix you supplied is a "
          << " matrix in " << BaseMatrix::storageType.ToString( sType )
          << " format." );
    }


    /**************** Build the auxiliary graph ********************/
    /* since we don't know the non-zero entries in the auxiliary matrix
     * a priori, we have to introduce a temporary graph-like structure:
     * For every node i, we store its neighbour-nodes j (connected via edges).
     * Every neighbour-entry corresponds also to a length (edge-length between
     * node i and j).
     * We can't use geomIndEdge_ because it's edge-based, we need
     * a node-based graph
     */

    // fill the indexNodeNum_ map... nodeNum is the key for entry in auxiliary matrix
    // also fill edgeIndNode_, basically the same as geomIndEdge_ only without the
    // unnecessary stuff
    UInt nRows = 0;
    std::unordered_map<Integer, StdVector<UInt> > tempNodesOfEdge;
    std::unordered_map<Integer, StdVector<Double> > tempLengthOfEdge;
    edgeIndNode_.Resize(geomIndEdge_.size());
    for(UInt i = 0; i < geomIndEdge_.size(); ++i ){
      StdVector<Integer> nodes = geomIndEdge_[i].eNodes;
      edgeIndNode_[i] = nodes;
      StdVector<UInt>& t1 = tempNodesOfEdge[nodes[0]];
      StdVector<UInt>& t2 = tempNodesOfEdge[nodes[1]];
      StdVector<Double>& et1 = tempLengthOfEdge[nodes[0]];
      StdVector<Double>& et2 = tempLengthOfEdge[nodes[1]];
      t1.Push_back( nodes[1] );
      t2.Push_back( nodes[0] );
      et1.Push_back( geomIndEdge_[i].length );
      et2.Push_back( geomIndEdge_[i].length );
      //tempTest[nodes[0]] = 2 * i;
      //tempTest[2 * i + 1] = nodes[1];
      if( !nodeNumIndex_.Contains( nodes[0]) ){
        nodeNumIndex_.Push_back( nodes[0] );
        indexNodeNum_[ nodes[0] ] = nRows;
        nRows++;
      }
      if( !nodeNumIndex_.Contains(nodes[1]) ){
        nodeNumIndex_.Push_back( nodes[1] );
        indexNodeNum_[ nodes[1] ] = nRows;
        nRows++;
      }
    }

    StdVector< StdVector<UInt> > nodeList;
    StdVector< StdVector<Double> > edgeList;
    nodeList.Resize(nRows);
    edgeList.Resize(nRows);
    // row pointer of auxiliary matrix
    StdVector<UInt> rP;
    rP.Push_back(0);
    // column index of auxiliary matrix
    StdVector<UInt> cI;
    // data pointer of auxiliary matrix
    StdVector<Double> dataAux;
    UInt rSize = 0;
    //std::unordered_map< Integer, EdgeGeom>::const_iterator eIt = geomIndEdge_.begin(); //old version
    // loop over every node and collect connected neighbours
    for(UInt nIt = 0; nIt < nRows; ++nIt){
      Integer nNum = nodeNumIndex_[nIt]; //real nodeNumber
      StdVector<UInt>& nL = nodeList[nIt];
      StdVector<Double>& eL = edgeList[nIt];
      //eIt = geomIndEdge_.begin(); //old version
      // insert template for diagonal element in auxiliary matrix
      cI.Push_back(nIt);
      dataAux.Push_back(0.0);

      // connected edges of node nNum
      StdVector<UInt>& e = tempNodesOfEdge[nNum];
      StdVector<Double>& l = tempLengthOfEdge[nNum];
      for( UInt eI = 0; eI < e.GetSize(); ++eI){
        if(!nL.Contains(indexNodeNum_[ e[eI] ])){
          nL.Push_back(indexNodeNum_[ e[eI] ]);
          eL.Push_back( l[eI] );
          cI.Push_back(indexNodeNum_[ e[eI] ]);
          dataAux.Push_back(0.0);
        }
      }

      rSize += nL.GetSize() +1;//+1 for diagonal
      rP.Push_back(rSize);
    }


    /**************** Build the auxiliary matrix ********************/
    // row- and column pointer already filled, now we have to assemble dataAux
    for(UInt i = 0; i < nRows; ++i){
      const StdVector<Double>& eL = edgeList[i];
      Double diag = 0.0;
      UInt d = 0;
      for(UInt j = rP[i] +1; j < rP[i+1]; ++j){
        Double invL = 1.0 / eL[d];
        dataAux[j] -= invL;
        diag += invL;
        ++d;
        }
      // insert diagonal element LEX_DIAG_FIRST
      dataAux[rP[i]] = diag;
      }

    UInt nnz = dataAux.GetSize();
    auxMatAMG_ = GenerateStdMatrixObject(BaseMatrix::EntryType::DOUBLE,
                                         BaseMatrix::SPARSE_NONSYM,
                                         nRows,
                                         nRows,
                                         nnz);

    // Generate auxiliary matrix
    LOG_DBG(algSys) << " Generating auxiliary matrix for edge-AMG solver\n"
                       " with " << nRows << " rows, " << nRows << " columns and "
                       << nnz << " non-zero entries";

    CRS_Matrix<Double>& auxMat = dynamic_cast<CRS_Matrix<Double>&>(*auxMatAMG_);
    auxMat.SetSparsityPatternData(rP, cI, dataAux);
    if(auxMat.GetCurrentLayout() != CRS_Matrix<Double>::LEX_DIAG_FIRST){
      auxMat.ChangeLayout(CRS_Matrix<Double>::LEX_DIAG_FIRST);
    }

  }




  void AlgebraicSys::BuildAMGLagrangeAuxMatrix(){

    // only singlefield...only one submatrix
    BaseMatrix& b = (*effMat_)(0,0);

    // Check that we have a StdMatrix
    if ( b.GetStructureType() != BaseMatrix::SPARSE_MATRIX ) {
      EXCEPTION( "AlgebraicSys::BuildAMGLagrangeAuxMatrix(): "
          << " AMG cannot deal with matrices other than StdMatrix" );
    }

    StdMatrix& stdMat = dynamic_cast<StdMatrix&>(b);

    // Now test the storage layout
    BaseMatrix::StorageType sType = stdMat.GetStorageType();
    if ( sType != BaseMatrix::SPARSE_NONSYM ) {
      EXCEPTION( "AlgebraicSys::BuildAMGLagrangeAuxMatrix(): AMG requires the"
          << " system matrix to be a CRS_Matrix i.e. sparseNonSym. "
          << " The system matrix you supplied is a "
          << " matrix in " << BaseMatrix::storageType.ToString( sType ) << " format." );
    }

    // Down-cast to CRS_Matrix
    CRS_Matrix<Double>& crsMat = dynamic_cast<CRS_Matrix<Double>&>(b);
    if( b.GetEntryType() == BaseMatrix::EntryType::COMPLEX){
      EXCEPTION("AlgebraicSys::BuildAMGLagrangeAuxMatrix() Complex version not yet implemented!");
    }


    // Get handles and info from system-matrix

    // Query number of matrix rows and columns
    UInt nRows = crsMat.GetNumRows();
    // Get hold of column index array
    const UInt *colP = crsMat.GetColPointer();
    // Get hold of row pointer index array
    const UInt *rowInd = crsMat.GetRowPointer();
    // Get hold of data array
    const Double* dataA = crsMat.GetDataPointer();

    UInt nnzSys = crsMat.GetNnz();

    // dimension of the problem (number of entries per node)
    UInt dim = dim_;

    // Declare some variables
    StdVector<UInt> rP, cI;
    StdVector<Double> dataAux;

    UInt sizeB;
    StdVector<Double> diagData;
    UInt rowSize = 0;
    StdVector<UInt> t;
    UInt diagPos = 0;
    Double rSum;
    Vector<Double> dist;

    switch(dim){
      case 1:
        amgType_ = SCALAR;
        /************ case of scalar nodal H1 FE ***************/
        if( crsMat.GetCurrentLayout() != CRS_Matrix<Double>::LEX_DIAG_FIRST ){
          crsMat.ChangeLayout(CRS_Matrix<Double>::LEX_DIAG_FIRST);
        }
        rP.Resize(nRows + 1, 0);
        UInt rIt;
        rIt = 0;
        for (UInt i = 0; i < nRows; i++ ) {
          for(UInt nzIt = rowInd[i]; nzIt < rowInd[i + 1]; ++nzIt){
            cI.Push_back(colP[nzIt]);
            dataAux.Push_back(dataA[nzIt]);
            rIt += 1;
          }
          rP[i + 1] = rIt;
        }

        auxMatAMG_ = GenerateStdMatrixObject(crsMat.GetEntryType(),
                                             BaseMatrix::SPARSE_NONSYM,
                                             nRows,
                                             nRows,
                                             nnzSys);
        break;

      case 2:
        amgType_ = VECTORIAL;
        /************ 2D nodal H1 FE ***************/
        /************ the "geometric way" **********/
        // sanity check
        if( nRows % 2 != 0) EXCEPTION("AlgebraicSys::BuildAMGAuxMatrix(): matrix is not 2D!!")
        // want to have the matrix in LEX_DIAG_FIRST layout
        if( crsMat.GetCurrentLayout() != CRS_Matrix<Double>::LEX_DIAG_FIRST ){
          crsMat.ChangeLayout(CRS_Matrix<Double>::LEX_DIAG_FIRST);
        }
        sizeB = nRows / 2;

        diagData.Resize(sizeB, 0.0);


        /* nz_entry          *     *     *     *     *     *
         * col-ind           0     1     2     3     4     5
         *                  |_______ |  |________|  |________|
         * auxMatAMG-col i      0           1            2
         *
         * conversion between nz-entry and auxMatAMG_ column (node): 2*i<= nz_entry<=2*(i+1)-1
         *                                               or inverse: (k-1)/2<=i<=k/2
         */
        rP.Resize(sizeB + 1, 0);
        rP[0] = 0;

        t.Resize(sizeB);
        t.Init(0);
        for(UInt i = 0; i < nRows; ++i){
          for(UInt j = rowInd[i]; j < rowInd[i + 1]; ++j){
            // find out to which node this entry belongs
            if( colP[j] % dim == 0 ){
              // node = colP[j] / 2
              if( t[colP[j]/2] == 0){
                t[colP[j]/2] = 1;
                cI.Push_back( colP[j] / 2 );
                rowSize++;
              }
            }else{
              // node = (colP[j]-1) / 2
              if(t[(colP[j]-1)/2] == 0){
                t[(colP[j]-1)/2] = 1;
                cI.Push_back( (colP[j]-1) / 2 );
                rowSize++;
              }
            }

          }// loop over rowsize
          if(i==1){
            rP[1] = rowSize;
            t.Init(0);
          }
          if( ((i-1) % 2 == 0 && i > 1)){
            rP[i/2 + 1] = rowSize;
            t.Init(0);
          }
          if(i == nRows-1){
            rP[sizeB] = rowSize;
            t.Init(0);
          }
        }// loop over every row


        /* Now we can start to build the actual auxiliary-matrix:
         * for every off-diagonal entry (i,j) with i != j, we calculate
         * the geometric distance between nodes i and j via geomInd_ .
         * The entry auxMatAMG_[i,j] = -1.0 / distance(i,j)^2
         * The diagonal entry i is -sum(auxMatAMG_(i,:)) ... negative
         * row-sum of off-diagonals
         * Defining the auxiliary-matrix in this way, we obtain a
         * strong diagonally dominant SPD matrix
         */
        dataAux.Resize(cI.GetSize());
        dataAux.Init(0.0);
        for( UInt i = 0; i < sizeB; ++i ){
          rSum = 0;
          Vector<Double> ci = geomInd_[dim * i];
          for( UInt j = rP[i]; j < rP[i+1]; ++j ){
            Vector<Double> cj = geomInd_[dim * cI[j] ];
            dist.Resize(cj.GetSize(), 0.0);
            dist.Add(1.0, ci, -1.0, cj);
            // if diagonal, store the position
            if( i == cI[j]){
              diagPos = j;
            }else{
              dataAux[j] =  -1.0 / dist.NormL2();
              rSum += -1.0 / dist.NormL2();
            }
          }
          // fill diagonal
          dataAux[diagPos] = -rSum;
        }

        auxMatAMG_ = GenerateStdMatrixObject(crsMat.GetEntryType(),
            BaseMatrix::SPARSE_NONSYM,
            sizeB,
            sizeB,
            cI.GetSize());
        break;

      case 3:
        amgType_ = VECTORIAL;
        /************ 3D nodal H1 FE ***************/
        /************ the "geometric way" **********/
        // sanity check
        if( nRows % 3 != 0) EXCEPTION("AlgebraicSys::BuildAMGAuxMatrix(): matrix is not 3D!!")
        // want to have the matrix in LEX_DIAG_FIRST layout
        if( crsMat.GetCurrentLayout() != CRS_Matrix<Double>::LEX_DIAG_FIRST ){
          crsMat.ChangeLayout(CRS_Matrix<Double>::LEX_DIAG_FIRST);
        }
        sizeB = nRows / 3;

        diagData.Resize(sizeB, 0.0);

        //TODO I simply can not access the basegraph, this would make
        // things way easier...
        // every nz-entry in auxMatAMG_ corresponds to a
        // dim x dim entry in the system-matrix

        // ok the plan now is to mark every dim x dim entry
        // if it is non zero, this means if the dim x dim entry
        // in the system matrix is nz, also the entry of dim 1
        // in the auxMatAMG_ is nz
        // therefore we build our own sparsity-pattern by introducing
        // a CRS-style col- and row-pointer but no data-pointer
        // since we only need the information that it's a nz entry

        /* nz_entry          *     *     *     *     *     *     *     *     *     *
         * col-ind           0     1     2     3     4     5     6     7     8     9
         *                  |______________|  |______________|  |______________|  |___
         * auxMatAMG-col i         0                  1                2
         *
         * conversion between nz-entry and auxMatAMG_ column (node): 3*i<= nz_entry<=3*(i+1)-1
         *                                               or inverse: (k-2)/3<=i<=k/3
         */
        rP.Resize(sizeB + 1, 0);
        rP[0] = 0;

        t.Resize(sizeB);
        t.Init(0);
        for(UInt i = 0; i < nRows; ++i){
          for(UInt j = rowInd[i]; j < rowInd[i + 1]; ++j){
            // find out to which node this entry belongs
            if( colP[j] % dim == 0 ){
              // node = colP[j] / dim
              if( t[colP[j]/3] == 0){
                t[colP[j ]/3] = 1;
                cI.Push_back( colP[j] / 3 );
                rowSize++;
              }
            }else{
              if( (colP[j]-2) % dim == 0 && colP[j] > 1 ){
                // node = (colP[j]-2) / dim
                if(t[(colP[j]-2)/3] == 0){
                  t[(colP[j]-2)/3] = 1;
                  cI.Push_back( (colP[j]-2) / 3 );
                  rowSize++;
                }
              }else{
                // node = (colP[j]-1) / dim
                // same as node = (colP[j]+1) / dim
                // because it's a middle-entry
                if(t[(colP[j]-1)/3] == 0){
                  t[(colP[j]-1)/3] = 1;
                  cI.Push_back( (colP[j]-1) / 3 );
                  rowSize++;
                }
              }
            }
          }// loop over rowsize
          if(i==2){
            rP[1] = rowSize;
            t.Init(0);
          }
          if( ((i-2) % 3 == 0 && i > 2)){
            rP[i/3 + 1] = rowSize;
            t.Init(0);
          }
          if(i == nRows-1){
            rP[sizeB] = rowSize;
            t.Init(0);
          }
        }// loop over every row

        /* Now we can start to build the actual auxiliary-matrix:
         * for every off-diagonal entry (i,j) with i != j, we calculate
         * the geometric distance between nodes i and j via geomInd_ .
         * The entry auxMatAMG_[i,j] = -1.0 / distance(i,j)^2
         * The diagonal entry i is -sum(auxMatAMG_(i,:)) ... negative
         * row-sum of off-diagonals
         * Defining the auxiliary-matrix in this way, we obtain a
         * strong diagonally dominant SPD matrix
         */
        dataAux.Resize(cI.GetSize());
        dataAux.Init(0.0);
        for( UInt i = 0; i < sizeB; ++i ){
          rSum = 0;
          Vector<Double> ci = geomInd_[dim * i];
          for( UInt j = rP[i]; j < rP[i+1]; ++j ){
            Vector<Double> cj = geomInd_[dim * cI[j] ];
            dist.Resize(cj.GetSize(), 0.0);
            dist.Add(1.0, ci, -1.0, cj);
            // if diagonal, store the position
            if( i == cI[j]){
              diagPos = j;
            }else{
              dataAux[j] =  -1.0 / dist.NormL2();
              rSum += -1.0 / dist.NormL2();
            }
          }
          // fill diagonal
          dataAux[diagPos] = -rSum;
        }

        auxMatAMG_ = GenerateStdMatrixObject(crsMat.GetEntryType(),
                                             BaseMatrix::SPARSE_NONSYM,
                                             sizeB,
                                             sizeB,
                                             cI.GetSize());



        crsMat.ChangeLayout(CRS_Matrix<Double>::LEX);

        break;

    }// switch dimension

    // Generate auxiliary matrix
    LOG_DBG(algSys) << " Generating auxiliary matrix for AMG solver\n"
                       " with " << nRows << " rows, " << nRows << " columns and "
                       << nnzSys << " non-zero entries";

    CRS_Matrix<Double>& auxMat = dynamic_cast<CRS_Matrix<Double>&>(*auxMatAMG_);
    auxMat.SetSparsityPatternData(rP, cI, dataAux);
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
  SetElementMatrix_MultHarm( FEMatrixType, Matrix<Double>&,
                    FeFctIdType, const StdVector<Integer>& ,
                    FeFctIdType, const StdVector<Integer>& , bool,
                    const StdVector<UInt>&);
  template void AlgebraicSys::
  SetElementMatrix_MultHarm( FEMatrixType, Matrix<Complex>&,
                    FeFctIdType, const StdVector<Integer>& ,
                    FeFctIdType, const StdVector<Integer>& , bool,
                    const StdVector<UInt>&);

  template void AlgebraicSys::
  SetElementRHS( const Vector<Double>&, const FeFctIdType,
                 StdVector<Integer>&, UInt&);
  template void AlgebraicSys::
    SetElementRHS( const Vector<Complex>&, const FeFctIdType,
                   StdVector<Integer>&, UInt&);

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

   template void AlgebraicSys::
   SetDirichletMH( const FeFctIdType, Integer, const Double& , UInt&);
   template void AlgebraicSys::
   SetDirichletMH( const FeFctIdType, Integer, const Complex&, UInt&);

  template void AlgebraicSys::ConstructEffectiveMatrix( const FeFctIdType fctId, const std::map<FEMatrixType,Complex> &matFactors, const bool isMultHarm);
  template void AlgebraicSys::ConstructEffectiveMatrix( const FeFctIdType fctId, const std::map<FEMatrixType,Double> &matFactors, const bool isMultHarm);
}// end of Namespace
