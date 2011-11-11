#include <iostream> 
#include <iomanip>
#include <fstream> 
#include <string> 

#include "OLAS/algsys/AlgebraicSys.hh"
#include "MatVec/SBM_Matrix.hh"
#include "MatVec/VBR_Matrix.hh"

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


DECLARE_LOG(algSys)
DEFINE_LOG(algSys, "algSys")

namespace CoupledField {

  
 

  // ***********************
  //   Default Constructor
  // ***********************
  AlgebraicSys::AlgebraicSys(PtrParamNode param, PtrParamNode info ) 
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
    
    rhs_ = NULL;
    sol_ = NULL;
    sbmSymm_            = true;
    statCond_           = false;
    effMat_             = NULL;
    effRhs_             = NULL;
    effSol_             = NULL;
    
    idbcHandler_    = NULL;
    assembleDirichletToSysMat_ = false;
    
    // Default is to always use a system matrix
    matrixTypes_.insert( SYSTEM );
    
    // Set flag for insertion of penalty terms into matrix
    PtrParamNode setupNode;
    setupNode = myParam_->Get("setup", ParamNode::INSERT );

    
    std::string aux = "penalty";
    setupNode->GetValue("idbcHandling", aux, ParamNode::INSERT );
    usingPenalty_ = aux == "penalty" ? true : false;
    
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
  }


  void AlgebraicSys::CreateLinSys() {
    
    LOG_TRACE(algSys) << "Creating linear system";
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
          std::set<UInt>::const_iterator colIt = colBlocks.begin();
          for( ; rowIt != rowBlocks.end(); ++rowIt ) {
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
      // Determine the set of sub-matrices that we need for the SYSTEM matrix
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

      // ------------------------------
      //  Generation of matrix objects
      // ------------------------------

      // Obtain some info from parameter file
      BaseMatrix::EntryType entryType;
      std::string entryStr = "double";
      entryType = BaseMatrix::entryType.Parse( entryStr );

      for ( fIt = matrixTypes_.begin(); fIt != matrixTypes_.end(); fIt++ ) {
        sysMat_[*fIt] = GenerateSBM_Matrix( *fIt, entryType );
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

      // Generate empty SBM vectors
      rhs_ = dynamic_cast<SBM_Vector*>
        ( GenerateVectorObject( *(sysMat_[SYSTEM]) ) );

      sol_ = dynamic_cast<SBM_Vector*>
        ( GenerateVectorObject( *(sysMat_[SYSTEM]) ) );
      
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

        // Insert sub-vector into solution
        bVec = GenerateVectorObject( *stdMat );
        sVec = dynamic_cast<SingleVector*>( bVec );
        sol_->SetSubVector( sVec, k );

        // Insert sub-vector into right-hand side
        bVec = GenerateVectorObject( *stdMat );
        sVec = dynamic_cast<SingleVector*>( bVec );
        rhs_->SetSubVector( sVec, k );
      }

      // ---------------------------------------
      //   Generate efficient matrices vectors
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
  }

  void AlgebraicSys::CreatePrecond() {
    
    LOG_TRACE(algSys) << "Creating preconditioner";
    
    PtrParamNode precondNode = myInfo_->Get("precond");
    
    // if we have just one SBM matrix block, use directly
    // the specialized methods for StdMatrices
    if( numBlocks_ == 1) {
      precond_ = GenerateStdPrecondObject((*effMat_)(0,0),
                                          myParam_, precondNode );
    } else {
      precond_ = GenerateSBMPrecondObject( *(effMat_), 
                                           myParam_, precondNode );
    }
    
  }

  void AlgebraicSys::CreateSolver() {
    
    LOG_TRACE(algSys) << "Creating solver";
    PtrParamNode solverNode = myInfo_->Get("solver");
    
    // if we have just one SBM matrix block, use directly
    // the specialized methods for StdMatrices
    if( numBlocks_ == 1) {
      solver_ = GenerateSolverObject( (*effMat_)(0,0), 
                                      myParam_, solverNode);
    } else {
      solver_ = GenerateSolverObject( *(effMat_), 
                                      myParam_, solverNode);
    }
     
  }

  void AlgebraicSys::CreateEigenSolver() {
    
    LOG_TRACE(algSys) << "Creating eigenvalue solver";
    REFACTOR;
  }

  void AlgebraicSys::SetupPrecond(PtrParamNode analysis_id)  {
    
    LOG_TRACE(algSys) << "Setup of preconditioner";
    
    // if we have just one SBM matrix block, use directly
    // the specialized methods for StdMatrices
    if( numBlocks_ == 1) {
      precond_->Setup( (*effMat_)(0,0), analysis_id);
    } else {
      precond_->Setup( *(effMat_), analysis_id);
    }

  }

  void AlgebraicSys::SetupSolver(PtrParamNode analysis_id) {
    
    LOG_TRACE(algSys) << "Setup of solver";
    
    // if we have just one SBM matrix block, use directly
    // the specialized methods for StdMatrices
    if( numBlocks_ == 1 ) {
      solver_->Setup( (*effMat_)(0,0), analysis_id );
    } else {
      solver_->Setup( *effMat_, analysis_id);
    }
    
   // in any case, pass preconditioner object to solver, if created
    if( precond_ != NULL ) {
      solver_->SetPrecond( precond_ );
    }
  }

  void AlgebraicSys::SetupEigenSolver( UInt numFreq, Double shift,
                                       bool quadratic ) {
    
    LOG_TRACE(algSys) << "Setup of eigenvalue solver";
    REFACTOR;
  }

  void AlgebraicSys::Solve(PtrParamNode analysis_id) {
    
    LOG_TRACE(algSys) << "Solving problem";
    
    // ======================================================================
    //  CHECK FOR CALCULATION OF CONDITION NUMBER
    // ======================================================================
    // Check, if condition number is to be calculated
      bool calcCapa = false;
      myParam_->Get("matrix", ParamNode::INSERT)
      ->GetValue("calcConditionNumber", calcCapa, ParamNode::INSERT );
      
      if( calcCapa ) {
        Double condNumber = 0.0;
        Vector<Double> ev, err;
        BaseEigenSolver * evs = 
          GenerateEigenSolverObject( *(sysMat_[SYSTEM]), myParam_, 
                                     myInfo_->Get("solve_eigen") );
        PtrParamNode in = myInfo_->
            Get(ParamNode::PROCESS)->Get("conditionNumber", ParamNode::APPEND);
        in->Get("analysis_id")->SetValue(analysis_id);
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
        // in the end prevent-relcaulation of evs by re-setting the value for CalcConditionNumber
        myParam_->Get("matrix")->Get("calcConditionNumber")->SetValue( boost::any(false) );
      }
    // ======================================================================
    
    // If the penalty formulation is used and we have inhomogeneous
    // Dirichlet boundary conditions, then the righ-hand side is
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

    // Now modifiy the right-hand side vector
    idbcHandler_->AddIDBCToRHS( effRhs_ );

    // Remove the export linear system stuff, it has changed in standardsys.cc an as below
    // the solve part is commentet out, I see no reason to export linsys also here, it would
    // require a generalization anyway. Fabian 16.11.07
    // check if we do export stuff
    PtrParamNode els = myParam_->Get("exportLinSys", ParamNode::PASS );
    std::string file;
    std::string base;

    // TODO: This is most ugly copy & paste from standardsys.cc -> Generalize common parts!!
    // need it common even when exclusive solution
    if(els) {
      std::ostringstream os;
      std::string name = els->Has("baseName") ? els->Get("baseName")->As<std::string>() : progOpts->GetSimName();
      os << name;
      //          std::string id = analysis_id->Get("analysis_id")->As<std::string>();
      //          boost::replace_all(id, ":", "_");
      //          os << "_" << id;
      base = os.str();
    }

    // check if we do not only want the solution
    if(els && els->Get("solution")->As<std::string>() != "exclusive")
    {
      // two formats. The harwell-boing format includes the rhs!
      if(els->Get("format")->As<std::string>()  == "harwell-boeing")
      {
        EXCEPTION( "Harwell-Boeing Format not implemented for SBM-case" );
      }
      else // classical (default) matrix-market
      {
        sysMat_[SYSTEM]->Export((base+".mtx").c_str() );

        // HARD-CODED: Export also preconditioner
        SBM_Matrix * copy = new SBM_Matrix(*(sysMat_[SYSTEM]));
        if( numBlocks_ == 1 ) {
          precond_->GetPrecondSysMat((*copy)(0,0));
        } else {
          precond_->GetPrecondSysMat(*copy);
        }
        copy->Export((base+"_precond.mtx").c_str());
        
        
        if(els->HasByVal("damping", true) && sysMat_[DAMPING] != NULL)
          sysMat_[DAMPING]->Export((base+"_damping.mtx").c_str() );

        if(els->HasByVal("auxiliary", true) && sysMat_[AUXILIARY] != NULL)
          sysMat_[AUXILIARY]->Export((base+"_aux.mtx").c_str() );

        // rhs is only in harwell-boing included
        rhs_->Export((base+".vec").c_str() );
      }
    }
    if(els && els->HasByVal("initialGuess", true))
      sol_->Export((base+"_intial_guess.vec").c_str());

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
    
    // Trigger solution
    if( numBlocks_ == 1 ) { 
      solver_->Solve( (*effMat_)(0,0), 
                      (*effRhs_)(0), (*effSol_)(0),
                      analysis_id);
    } else {
      solver_->Solve( *effMat_, *effRhs_, 
                      *effSol_, analysis_id );
    }
    
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
        for(UInt c = 0; c < numBlocks_ -1; ++c ) {
          StdMatrix &stdMat =(*sysMat_[SYSTEM])(c,numBlocks_-1);
          stdMat.MultTSub((*sol_)(c),*tmp);
        }
        S_ii.Mult(*tmp, (*sol_)(numBlocks_-1));
      } else {
        EXCEPTION("Non-symmetric case not yet implemented");
      }
      // delete temporary vector
      delete tmp;
    }

    // Export solution if desired
    if(els && els->Get("solution")->As<std::string>() != "no")
      sol_->Export((base+".sol.vec").c_str());

    // Now de-modifiy the right-hand side vector
    idbcHandler_->RemoveIDBCFromRHS( effRhs_ );

    // Check that solution went fine, if not issue a warning
    if ( out->Get("solutionIsOkay")->As<bool>() == false ) {
      WARN("Solver reports a problem! Consult .las file for "
          << "further diagnostics!");
    }
  }

  void AlgebraicSys::CalcEigenFrequencies( Vector<Double>& frequencies,
                                           Vector<Double>& err ) {
    
    LOG_TRACE(algSys) << "Calculating real-valued eigenfrequencies";
    REFACTOR;
  }

  void AlgebraicSys::CalcEigenFrequencies( Vector<Complex>& frequencies,
                                           Vector<Double>& err ) {
    
    LOG_TRACE(algSys) << "Calculating complex-valued eigenfrequencies";
    REFACTOR;
  }

  void AlgebraicSys::CalcEigenMode( UInt numMode )  {
    
    LOG_TRACE(algSys) << "Calculating eigenmode #" << numMode;
    REFACTOR;

  }

  void AlgebraicSys::GraphSetupInit( UInt numFcts, UInt numBlocks,
                                     bool useDistinctGraphs,
                                     bool staticCondensation ) {
    
    LOG_DBG(algSys) << "Setup matrix graph for " << numFcts 
                    << " function(s) with " << numBlocks << " block(s)";
    LOG_DBG(algSys) << "Use distinct graphs:" << useDistinctGraphs << std::endl;
    LOG_DBG(algSys) << "Static Condensation: " 
                    << (staticCondensation ? "true" : "false"); 
    
    
    // feFunction specific data
    numFcts_ = numFcts;
    numEqnsPerFct_.Resize(numFcts);
    lastFreeEqnPerFct_.Resize(numFcts);
    eqnToSBMBlock_.Resize(numFcts);
    
    // store flag for applying static condensation
    statCond_ = staticCondensation;
    
    // block specific data
    numBlocks_ = numBlocks;
    
    // If we have static condensation, the last diagonal block is the 
    // inner-inner one, which will not be connected to a Dirichlet block. 
    if( staticCondensation ) {
      numDirichletValuesPerBlock_.Resize(numBlocks-1);
    } else {
      numDirichletValuesPerBlock_.Resize(numBlocks);
    }
    blockInfo_.Resize(numBlocks_);
    for( UInt i = 0; i < numBlocks_; ++i ) {
      blockInfo_[i] = new GraphManager::SBMBlockInfo();
    }
    
    // create new graph manager object and initialize it
    graphManager_ = new GraphManager();
    graphManager_->SetupInit( numBlocks_, useDistinctGraphs );
  }

  
  // ***************
  //   ObtainFctId
  // ***************
  FeFctIdType AlgebraicSys::ObtainFctId( const std::string& fctString ) {

   LOG_TRACE(algSys) << "Obtaining FctId for fct '" << fctString << "'";

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


    // Remember number of equation numbers for each function
    numEqnsPerFct_[fctId] = numEqns;
    lastFreeEqnPerFct_[fctId] = numLastFreeEqn;
  }
  
  
  void AlgebraicSys::
  DefineSBMMatrixBlock( UInt sbmIndex, 
                        const std::map<FeFctIdType,std::set<Integer> >& eqns ) {
    
    // Just logging output
    LOG_TRACE(algSys) << "Defining SBM block #" << sbmIndex;
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
        eqnToSBMBlock_[fctId][actEqn] = sbmIndex;
        
        // also remember the block, in which this functionId occurs
        fctIdsInBlocks_[fctId].insert(sbmIndex);
      }
      
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
        }
      }
      
      // Store number of Dirichlet values per SBM block
      // Note: if we perform static condensation, we have one block
      // less for the IDBC-hanlder
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
        // Case b): If we use elimination approach, we create the graphmanager
        //          creates an additional IDBC-graph and sorts the indices /
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
          std::map<UInt, UInt> & eqnToIndex = bi.eqnToIndex[iFct];
          std::map<UInt, UInt>::iterator eqnToIndexIt = eqnToIndex.begin(); 
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
    } // loop over fctIds
    
    // Determine reordering 
    // Note: At the moment we just allow for one reordering type, which is  the one
    // taken from the <matrix> element.
    BaseOrdering::ReorderingType reorder = BaseOrdering::NOREORDERING;
    std::string reorderStr = "noReordering";
    myParam_->Get( "matrix",ParamNode::INSERT)
                    ->GetValue( "reordering", reorderStr, ParamNode::INSERT );
    reorder = BaseOrdering::reorderingType.Parse( reorderStr ); 

    // Now register block with graph manager
    graphManager_->RegisterBlock( sbmIndex, blockInfo_[sbmIndex], reorder );
    
  }

  void AlgebraicSys::RegisterSubMatrixBlocks( UInt sbmIndex, UInt numMinorBlocks ) {
    
    LOG_TRACE(algSys) << "Registering " << numMinorBlocks 
        << " sub-matrix blocks for SBM block #" << sbmIndex;

    blockInfo_[sbmIndex]->indexBlocks.Resize( numMinorBlocks );
    blockInfo_[sbmIndex]->hasSubBlocks = true;
  } 


  void AlgebraicSys::DefineSubMatrixBlocks( UInt sbmIndex, UInt blockIndex,
                                            const StdVector<FeFctIdType>& fctIds,
                                            const StdVector<Integer>& eqns ) {

    LOG_DBG2(algSys) << "Defining sub-matrix block " << blockIndex 
                    << " for SBM block #" << sbmIndex;
    
    
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
    }
  }


  // *******************
  //   GetFEMatrixType
  // *******************
  void AlgebraicSys::GetFEMatrixTypes( std::set<FEMatrixType> &matrixTypes ) const{

    matrixTypes = matrixTypes_;

  }

  void AlgebraicSys::SetFEMatrixType( const FEMatrixType matrixType,
                                      const FeFctIdType fctId1,
                                      const FeFctIdType fctId2) {

    LOG_TRACE(algSys) << "Setting matrix type '" << feMatrixType.ToString(matrixType)
                           << "' for fct-Ids (" << fctId1 << ", " << fctId2 << ")";


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

      // -----------------------------------------------------------
      // NOTE: We temporarily insert also the transposed sub-matrix
      if ( fctId1 != fctId2 ) {
        sID.colInd = fctId1;
        sID.rowInd = fctId2;

        // Insert sub-matrix identifier into
        // corresponding FE-Matrix set
        feSubMatricesByFctId_[matrixType].insert( sID );
      }
      // -----------------------------------------------------------

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

    // Print information about registered functions
    PrintRegistrationInfo( );
    
    // Finalize graph manager setup
    graphManager_->SetupDone();
    
    // Now we have all graphs and IDBC in their re-order state,
    // so we have to fetch the reordering array from the graphmanager and 
    // update information in the blockInfo array. 
    
    // Loop over all blocks
    for( UInt iBlock = 0; iBlock < numBlocks_; ++iBlock ) {
      GraphManager::SBMBlockInfo &bi = *blockInfo_[iBlock];
      UInt numFcts = bi.eqnToIndex.GetSize();
      
      // Obtain reordering vector
      StdVector<UInt> newOrder;
      graphManager_->GetReordering(iBlock, newOrder);

      // Loop over all functions
      for( UInt iFct = 0; iFct < numFcts; ++iFct ) {
        std::map<UInt, UInt> & eqnToIndex = bi.eqnToIndex[iFct];
        std::map<UInt, UInt>::iterator it = eqnToIndex.begin();

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


  // ****************
  //   AssembleInit
  // ****************
  void AlgebraicSys::AssembleInit( const FeFctIdType fctId1,
                                   const FeFctIdType fctId2,
                                   bool isSymmetric,
                                   bool assemblingTranspose ) {
    
    LOG_DBG(algSys) << "Starting assembly of block ("
        << fctId1 << ", " << fctId2 << ")";
    
    // Basic idea: if at least one combination (fctId1, fctId2)
    // is not symmetric, so will be also the final SBM-system
    if( !isSymmetric ) {
      sbmSymm_ = false;
    }
  }


  // ****************
  //   AssembleDone
  // ****************
  void AlgebraicSys::AssembleDone( const FeFctIdType fctId1,
                                   const FeFctIdType fctId2,
                                   bool assemblingTranspose ) {

    LOG_DBG(algSys) << "Finishing assembly of block ("
                    << fctId1 << ", " << fctId2 << ")";
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
    LOG_DBG2(algSys) << "EqnVec2: " << eqnNrs1.ToString();
    
    // Re-map entries from (fctId,eqnNr) -> (blockNum,index)
    StdVector<UInt> blockNums1, blockNums2, indices1, indices2;
    MapFctIdEqnToIndex(fctId1, eqnNrs1, blockNums1, indices1);
    MapFctIdEqnToIndex(fctId2, eqnNrs2, blockNums2, indices2);

    // Quirk: If fctId1 == fctId2, we normally
    // need not set the counterPart. If however, the are now in
    // off diagonal matrices, we have to ensure, that originally
    // symmetric matrix structures remain symmetric also in the 
    // blockSystem
    if( fctId1 == fctId2 ) 
      setCounterPart = true;
    graphManager_->SetElementPos( blockNums1, indices1, 
                                  blockNums2, indices2,
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
    std::map<Integer,UInt>& eqnToBlock = eqnToSBMBlock_[fctId];
    
    UInt numEqns = eqns.GetSize();
    for( UInt iEqn = 0; iEqn < numEqns; ++iEqn ) {
      const UInt & eqnNr = eqns[iEqn];
      
      // take care of homogeneous BCs
      if( eqnNr == 0) {
        blockNums[iEqn] = 0;
        indices[iEqn] = 0;
      } else {
        const UInt & blockNum = eqnToBlock[eqnNr];
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
      const UInt & eqnNr = eqns[iEqn];
      const FeFctIdType & fctId = fctIds[iEqn];
      
      // get hold of fct-specific map
      std::map<Integer,UInt>& eqnToBlock = eqnToSBMBlock_[fctId];

      // take care of homogeneous BCs
      if( eqnNr == 0) {
        blockNums[iEqn] = 0;
        indices[iEqn] = 0;
      } else {
        const UInt & blockNum = eqnToBlock[eqnNr];
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
        const std::map<UInt, UInt> & eqnToIndexSet = 
            blockInfo_[iBlock]->eqnToIndex[actFctId];
        std::map<UInt, UInt>::const_iterator eqnIt = eqnToIndexSet.begin();
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
      std::map<Integer,UInt>& eqnToBlock = eqnToSBMBlock_[fctId];
      blockNum = eqnToBlock[eqnNr];
      index = blockInfo_[blockNum]->eqnToIndex[fctId][eqnNr];
    }
  }



  // all methods regarding assembly 
  void AlgebraicSys::InitMatrix( FEMatrixType matrixType,
                                 const FeFctIdType fctId ) {
    
    LOG_TRACE(algSys) << "Initializing matrix " << feMatrixType.ToString(matrixType)
                      << " for fctId " << fctId;
    
    // If matrix specified init this one
    if ( matrixType != NOTYPE )
    {
      assert(sysMat_[matrixType] != NULL);
      sysMat_[matrixType]->Init();
      idbcHandler_->InitMatrix(matrixType);
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
    if( newRHS.GetSize() != numFcts_ ) {
      EXCEPTION( "New rhs consists of " << newRHS.GetSize()
                 << " sub-vectors, the RHS of the algebraic system of "
                 << rhs_->GetSize() << " entries." )
    }
    
    // loop over all feFctIDs
    for(UInt i = 0; i < numFcts_; ++i ) {

      // get all (blockId,index)-combinations for the current fctId
      StdVector<UInt> blockNums, indices;
      MapCompleteFctIdToIndex( i, blockNums, indices);
      UInt size = blockNums.GetSize();
  
      // security check: ensure that sub-vector has the same size
      // as the block indices
      if( newRHS(i).GetSize() != indices.GetSize() ) {
        EXCEPTION( "Number of entries of " << i << "-th sub-vector and number "
                   "of indices do not match!");
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
                                       bool setCounterPart ) {
    
    LOG_DBG(algSys) << "Setting element matrix for fctIds ("
                     << fctId1 << ", " << fctId2 << ")";
    LOG_DBG2(algSys) << "Matrix: " << feMatrixType.ToString(matrixType);
    LOG_DBG2(algSys) << "EqnVec1: " << eqnNrs1.ToString();
    LOG_DBG2(algSys) << "EqnVec2: " << eqnNrs1.ToString();
    LOG_DBG3(algSys) << "matrix is:\n " << elemMat.ToString();
    
    // Re-map entries from (fctId,eqnNr) -> (blockNum,index)
    StdVector<UInt> rowBlocks, colBlocks, rowNums, colNums;
    MapFctIdEqnToIndex(fctId1, eqnNrs1, rowBlocks, rowNums);
    MapFctIdEqnToIndex(fctId2, eqnNrs2, colBlocks, colNums);
    
    
    // Now, dismantle equations
    std::vector<std::vector<UInt> > rowList1(numBlocks_);
    std::vector<std::vector<UInt> > rowIndList1(numBlocks_);
    std::vector<std::vector<UInt> > rowList2(numBlocks_);
    std::vector<std::vector<UInt> > rowIndList2(numBlocks_);
    std::vector<std::vector<UInt> > colList1(numBlocks_);
    std::vector<std::vector<UInt> > colIndList1(numBlocks_);
    std::vector<std::vector<UInt> > colList2(numBlocks_);
    std::vector<std::vector<UInt> > colIndList2(numBlocks_);
    UInt numRows = rowBlocks.GetSize();
    UInt numCols = colBlocks.GetSize();
    
    // Loop over all rows
    for( UInt iRow = 0; iRow < numRows; ++iRow ) {
      // get hold of block numbers and indices
      const UInt & rowBlock = rowBlocks[iRow];
      const UInt & rowNum = rowNums[iRow];

      // Compute index of graph in graph pointer matrix
      // get hold of vertex and edgelists
      std::vector<UInt> & rList1 = rowList1[rowBlock];
      std::vector<UInt> & rIndList1 = rowIndList1[rowBlock];
      std::vector<UInt> & rList2 = rowList2[rowBlock];
      std::vector<UInt> & rIndList2 = rowIndList2[rowBlock];
      // get limits of free indices
      const UInt & lastFreeRowIndex = blockInfo_[rowBlock]->numLastFreeIndex;

      // STEP 1: Generate row index list from first connect array, dropping
      //         equation numbers for dofs fixed by (in)homogeneous Dirichlet
      //         boundary conditions and changing the sign of those fixed by
      //         constraints.
      if ( rowNum > 0 ) {
        if ( rowNum > lastFreeRowIndex ) {
          rList2.push_back( rowNum - lastFreeRowIndex - 1 );
          rIndList2.push_back( iRow );
        } else {
          rList1.push_back( rowNum - 1);
          rIndList1.push_back( iRow );
        }
      }
    }

    // Loop over all columns
    for( UInt iCol = 0; iCol < numCols; ++iCol ) {

      // get hold of block numbers and indices
      const UInt & colBlock = colBlocks[iCol];
      const UInt & colNum = colNums[iCol];

      // get hold of vertex and edgelists
      std::vector<UInt> & cList1 = colList1[colBlock];
      std::vector<UInt> & cIndList1 = colIndList1[colBlock];
      std::vector<UInt> & cList2 = colList2[colBlock];
      std::vector<UInt> & cIndList2 = colIndList2[colBlock];

      // get limits of free indices
      const UInt & lastFreeColIndex = blockInfo_[colBlock]->numLastFreeIndex;

      // STEP 2: Split the second connect array into two edge lists, one for
      //         the graph and one for the IDBCgraph (which handles the indices
      //         fixed by inhomogeneous Dirichlet boundary conditions)
      if( colNum > 0 ) {
        if ( colNum > lastFreeColIndex ) {
          cList2.push_back( colNum - lastFreeColIndex - 1);
          cIndList2.push_back( iCol );
        }
        else {
          cList1.push_back( colNum - 1);
          cIndList1.push_back( iCol );
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
    if( statCond_ ) {
      LOG_DBG(algSys) << "Performing static condensation";
      
      // For static condensation, we need all matrices explicitly.
      // Loop over all blocks and extract matrices
      StdVector<Matrix<T> > matrices(numBlocks_*numBlocks_);
      for( UInt iRow = 0; iRow < numBlocks_; iRow++ ) {
        for( UInt iCol = 0; iCol < numBlocks_; iCol++ ) {
          Matrix<T>& mat = matrices[numBlocks_ * iRow + iCol];
          elemMat.GetSubMatrixByInd( mat, rowIndList1[iRow], 
                                     colIndList1[iCol]);
        } // col blocks
      } // row blocks
      
      // Invert inner block and store it back into element matrix
      matrices.Last().Invert_Lapack();
      elemMat.SetSubMatrixByInd( matrices.Last(),
                                 rowIndList1[numBlocks_-1],
                                 colIndList1[numBlocks_-1] );

      Matrix<T>& invMat = matrices.Last();
      Matrix<T> temp;
      
      // Loop over all blocks again, multiply matrices and re-insert afterwards
      for( UInt iRow = 0; iRow < numBlocks_; iRow++ ) {
        for( UInt iCol = 0; iCol < numBlocks_; iCol++ ) {
          // calc:  K_rc = K_rc - K_ri * K_ii^-1 * K_ic
          Matrix<T>& k_rc = matrices[numBlocks_ * iRow + iCol];
          Matrix<T>& k_ri = matrices[numBlocks_ * iRow + numBlocks_-1];
          Matrix<T>& k_ic = matrices[numBlocks_ * (numBlocks_-1) + iCol];
          temp.Resize(invMat.GetNumRows(),k_ic.GetNumCols());
          
          // Fast method: use BLAS 
          invMat.Mult_Blas(k_ic, temp,false,false,1.0,0.0);
          k_ri.Mult_Blas(temp,k_rc,false,false,-1.0,1.0);
          // Alternative solution without BLAS
//            temp = invMat * k_ic;
//            k_rc -= k_ri*temp;
          elemMat.SetSubMatrixByInd( k_rc,rowIndList1[iRow],
                                     colIndList1[iCol]);
        } // col blocks
      } // row blocks
    } // if static cond
    // ======================================================================
    
    // loop over all blocks and pass for every block the information to
    // the corresponding graph / IDBC graph
    for( UInt sbmRow = 0; sbmRow < numBlocks_; ++sbmRow ) {
      
      // Maybe we have to switch here depending on transpose
      for( UInt sbmCol = 0; sbmCol < numBlocks_; ++sbmCol ) {
        
        StdMatrix * stdMat = actMat->GetPointer(sbmRow, sbmCol);
        std::vector<UInt> & rList1 = rowList1[sbmRow];
        std::vector<UInt> & rList2 = rowList2[sbmRow];
        std::vector<UInt> & cList1 = colList1[sbmCol];
        std::vector<UInt> & cList2 = colList2[sbmCol];
        std::vector<UInt> & rIndList1 = rowIndList1[sbmRow];
        std::vector<UInt> & rIndList2 = rowIndList2[sbmRow];
        std::vector<UInt> & cIndList1 = colIndList1[sbmCol];
        std::vector<UInt> & cIndList2 = colIndList2[sbmCol];

        // Attention: This check is not really implemented in a clean way!
        if( stdMat != NULL ) {
          // 1) Assemble all free <-> free entries
          // loop over all rows/col
          for ( UInt i = 0; i < rList1.size(); i++ ) {
            rowInd = rIndList1[i];
            for ( UInt j = 0; j < cList1.size(); j++ ) {
              colInd = cIndList1[j];
              stdMat->AddToMatrixEntry( rList1[i], cList1[j],
                                        elemMat[rowInd][colInd] );
            } //j
          } //i
        }
        
        // 2) Assemble all free <-> fixed entries
        for ( UInt i = 0; i < rList1.size(); i++ ) {
           rowInd = rIndList1[i];
           for ( UInt j = 0; j < cList2.size(); j++ ) {
             colInd = cIndList2[j];
             idbcHandler_->AddWeightFixedToFree( matrixType, sbmRow, sbmCol,
                                                rList1[i], cList2[j],
                                                elemMat[rowInd][colInd]);
           } // j
        } // i
        
        // 3) Assemble all free <-> fixed entries ( TRANSPOSED )
        if( sbmRow != sbmCol && setCounterPart == true) {
          for ( UInt i = 0; i < rList2.size(); i++ ) {
            rowInd = rIndList2[i];
            for ( UInt j = 0; j < cList1.size(); j++ ) {
              colInd = cIndList1[j];
              idbcHandler_->AddWeightFixedToFree( matrixType, sbmCol, sbmRow,
                                                  rList1[sbmCol], cList2[sbmRow],
                                                  elemMat[rowInd][colInd]);
            } // j
          } // i
        }
        
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
    
    // Re-map entries from (fctId,eqnNr) -> (blockNum,index)
    StdVector<UInt> rowBlocks, rowNums;
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
    
    REFACTOR;
  }


  void AlgebraicSys::UpdateRHS(FEMatrixType matrixType, 
                               const SBM_Vector& fup) {
    
    LOG_TRACE(algSys) << "Updating RHS of matrix " 
                      << feMatrixType.ToString(matrixType);
    REFACTOR;
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

  void AlgebraicSys::ConstructEffectiveMatrix( const std::map<FEMatrixType,Double> &matFactors ) {

    LOG_TRACE(algSys) << "Constructing efffective system matrix";
    if (IS_LOG_ENABLED(algSys, dbg)) {
      LOG_DBG(algSys) << "Factors are:";
      std::map<FEMatrixType,Double>::const_iterator it = matFactors.begin();
      for( ; it != matFactors.end(); ++it ) {
        LOG_DBG(algSys) << feMatrixType.ToString(it->first) << ": " << it->second;
      }
    }

    factorMap::const_iterator it;
    SBM_Matrix *sys = sysMat_[SYSTEM];

    // It's okay, if there are no factors, if there is only a system
    // matrix and no other ones
    if ( matFactors.empty() == true ) {
      if ( matrixTypes_.size() == 1 && sys != NULL ) {
        // Also assemble the effective auxilliary system matrix for moving
        // IDBCs to the right-hand side
        idbcHandler_->BuiltSystemMatrix( matFactors );

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
        sys->Add( (*it).second, *sysMat_[(*it).first] );
      }
    }

    // Also assemble the effective auxilliary system matrix for moving
    // IDBCs to the right-hand side
    idbcHandler_->BuiltSystemMatrix( matFactors );

  }

  template<typename T>
  void AlgebraicSys::SetDirichlet( const FeFctIdType fctId,
                                   Integer eqnNr, const T& val ) {
    
    LOG_DBG(algSys) << "Setting Dirichlet value " << val << " for eqn " << eqnNr
                    << " of fctId " << fctId;
    
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
                                     const FeFctIdType fctId) {
    
    LOG_TRACE(algSys) << "Getting solution values of fct " << fctId;

    // get all (blockId,index)-combinations for the current fctId
    StdVector<UInt> blockNums, indices;
    MapCompleteFctIdToIndex( fctId, blockNums, indices);
    UInt size = blockNums.GetSize();
    ptSol.Resize(size);

    if( ptSol.GetEntryType() == BaseMatrix::DOUBLE ) {
      Vector<Double> & retVec = dynamic_cast<Vector<Double>&>( ptSol );
      Double entry = 0.0;
      for( UInt i = 0; i < size; ++i ) {

        // if index number is larger the lastFree dof, insert Dirichlet value
        if( indices[i] > blockInfo_[blockNums[i]]->numLastFreeIndex) {
          idbcHandler_->GetIDBC(blockNums[i],  indices[i], entry);
        } else {
          sol_->GetPointer(blockNums[i])->GetEntry(indices[i]-1,entry);
        }
        retVec[i] = entry;
      }

    } else {
      Vector<Complex> & retVec = dynamic_cast<Vector<Complex>&>( ptSol );
      Complex entry = 0.0;
      for( UInt i = 0; i < size; ++i ) {

        // if index number is larger the lastFree dof, insert Dirichlet value
        if( indices[i] > blockInfo_[blockNums[i]]->numLastFreeIndex) {
          idbcHandler_->GetIDBC(blockNums[i],  indices[i], entry);
        } else {
          sol_->GetPointer(blockNums[i])->GetEntry(indices[i]-1,entry);

        }
        retVec[i] = entry;
      }
    }
  }

  
  void AlgebraicSys::GetSolutionVal( SBM_Vector& solVec ) {
    
    // resize solVec to match number of functions
    solVec.Resize( numFcts_);
    
    // loop over all feFctIDs
    for(UInt i = 0; i < numFcts_; ++i ) {
    
      // call specialized GetSolutionVal method
      GetSolutionVal(solVec(i), i);
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
                                              BaseMatrix::EntryType entryType ) {

    // STEP 1: Generate empty SBM_Matrix
    SBM_Matrix *retMat = NULL;
    retMat = new SBM_Matrix( numBlocks_, numBlocks_, sbmSymm_ );
    if ( retMat == NULL ) {
      EXCEPTION( "SBM_System::GenerateSBM_Matrix: "
          << "This is the end my friend!\n"
          << "Generation of empty SBM_Matrix failed!" );
    }

    
    // TEMPORARY: find out matrix storage type by looking at the matrix element
    // determine matrix storage type
    BaseMatrix::StorageType sT = BaseMatrix::SPARSE_SYM;
    PtrParamNode matrixNode = myParam_->Get("matrix");
    if( matrixNode ) {
      std::string storageString;
      matrixNode->GetValue("storage", storageString );
      sT = BaseMatrix::storageType.Parse(storageString);
    }
    
    // STEP 2: Populate with sub-matrices
    std::set<SubMatrixID,SortSubMatrixID>::iterator sIt;
    BaseGraph *graph = NULL;
    for ( sIt = feSubMatricesByBlocks_[matType].begin();
        sIt != feSubMatricesByBlocks_[matType].end(); sIt++ ) {

      // Determine row / col
      UInt sbmRow = (*sIt).rowInd;
      UInt sbmCol= (*sIt).colInd;

      // Determine number of matrix rows and columns
      UInt nrows = blockInfo_[sbmRow]->numLastFreeIndex;
      UInt ncols = blockInfo_[sbmCol]->numLastFreeIndex;

      // Check for necessity of generation
      if ( sbmRow <= sbmCol || sbmSymm_ == false ) {

        // Trigger generation of sub-matrix
        graph = graphManager_->GetGraph( sbmRow, sbmCol );
        if ( sbmRow == sbmCol && sbmSymm_ == true ) {
          // for diagonal blocks we allow a variable
          // matrix layout
          
          // ====================
          //  HARD-CODED SECTION
          // ====================
          if( sbmRow == 0) {
            sT = BaseMatrix::SPARSE_SYM;
            
            
          } else {
            sT = BaseMatrix::VAR_BLOCK_ROW;
          }
          
          
          retMat->SetSubMatrix ( sbmRow, sbmCol, entryType, 
                                 sT,
                                 nrows, ncols, graph->GetNNE() );
        } else {
           //Test: Check, if performance increases by use of VBR-matrix
//          if( sbmRow > 0 && sbmCol > 0) {
//          retMat->SetSubMatrix ( sbmRow, sbmCol, 
//                                 entryType, BaseMatrix::VAR_BLOCK_ROW,
//                                 nrows, ncols, graph->GetNNE() );
//          } else {
          retMat->SetSubMatrix ( sbmRow, sbmCol, 
                                 entryType, BaseMatrix::SPARSE_NONSYM,
                                 nrows, ncols, graph->GetNNE() );
//          }
        }

        // Set sparsity pattern of sub-matrix
        (*retMat)( sbmRow, sbmCol ).SetSparsityPattern( *graph );
      }
    }

    return retMat;
  }
  
  void AlgebraicSys::PrintFeMatrixInfo( ) {

    LOG_TRACE(algSys) << "Print matrix information";
    
    PtrParamNode setupNode = myInfo_->Get("setup");

    // Print overview of defined matrices
    setupNode->SetComment("List of defined matrices");
    PtrParamNode matrixListNode = setupNode->Get("matrices");
    
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
            Double(stdMat->GetNumRows() * stdMat->GetNumCols()); 
        mNode->Get("fillLevelPerCent")->SetValue(fillLevel);
        
        
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
        
        // bandwidth gets just written for diagonal blocks 
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
      
      totalFillLevelPerCent = totalNumNonZeros * 100 / Double( totalNumRows * 
                              totalNumCols);
      matTypeNode->Get("totalNumRows")->SetValue(totalNumRows);
      matTypeNode->Get("totalNumCols")->SetValue(totalNumCols);
      matTypeNode->Get("totalNumNonZeros")->SetValue(totalNumNonZeros);
      matTypeNode->Get("totalFillLevelPerCent")->SetValue(totalFillLevelPerCent);
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
                    FeFctIdType, const StdVector<Integer>& , bool);
  template void AlgebraicSys::
  SetElementMatrix( FEMatrixType, Matrix<Complex>&, 
                    FeFctIdType, const StdVector<Integer>& ,
                    FeFctIdType, const StdVector<Integer>& , bool);
  
  template void AlgebraicSys::
  SetElementRHS( const Vector<Double>&, const FeFctIdType, 
                 StdVector<Integer>&);
  template void AlgebraicSys::
    SetElementRHS( const Vector<Complex>&, const FeFctIdType, 
                   StdVector<Integer>&);
  
  template void AlgebraicSys::SetNodeRHS(Double, FeFctIdType, Integer );
  template void AlgebraicSys::SetNodeRHS(Complex, FeFctIdType, Integer );
  
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
