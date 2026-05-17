#include "MatVec/opdefs.hh"
#include "MatVec/Vector.hh"
#include "MatVec/generatematvec.hh"
#include "OLAS/precond/BasePrecond.hh"
#include "OLAS/solver/COCRSolver.hh"

#include "DataInOut/Logging/LogConfigurator.hh"

namespace CoupledField {

  DEFINE_LOG(cocrsolver, "cocrsolver")

  // ***************
  //   Constructor
  // ***************
  template<typename T>
  COCRSolver<T>::COCRSolver( PtrParamNode solverNode,
                             PtrParamNode olasInfo ) {
    xml_      = solverNode;
    infoNode_ = olasInfo->Get( "cocr" );
  }

  // **********************
  //   Default destructor
  // **********************
  template<typename T>
  COCRSolver<T>::~COCRSolver() {
    delete r_;
    delete p_;
    delete q_;
    delete z_;
    delete qtld_;
    delete az_;
  }

  // **************************
  //   Setup (public version)
  // **************************
  template<typename T>
  void COCRSolver<T>::Setup( BaseMatrix &sysMat ) {
    PrivateSetup( sysMat );
  }

  // ***************************
  //   Setup (private version)
  // ***************************
  template<typename T>
  void COCRSolver<T>::PrivateSetup( const BaseMatrix &sysMat ) {

    // SBM matrices are not yet supported. Fail loudly rather than silently
    // miscompute via a dynamic_cast that would later return NULL.
    if ( sysMat.GetStructureType() == BaseMatrix::SBM_MATRIX ) {
      EXCEPTION( "COCRSolver currently supports StdMatrix (single-field) "
                 "problems only. SBM_Matrix (multi-physics) is not yet "
                 "supported." );
    }

    UInt newDim = sysMat.GetNumCols();

    if ( newDim > vectorLength_ ) {

      bool allocate = false;

      if ( r_ != NULL ) {

        // Try to resize the existing vectors first
        Vector<T> *rr  = dynamic_cast<Vector<T>*>( r_    );
        Vector<T> *pp  = dynamic_cast<Vector<T>*>( p_    );
        Vector<T> *qq  = dynamic_cast<Vector<T>*>( q_    );
        Vector<T> *zz  = dynamic_cast<Vector<T>*>( z_    );
        Vector<T> *qqt = dynamic_cast<Vector<T>*>( qtld_ );
        Vector<T> *aaz = dynamic_cast<Vector<T>*>( az_   );

        if ( rr != NULL ) {
          rr->Resize( newDim );
          pp->Resize( newDim );
          qq->Resize( newDim );
          zz->Resize( newDim );
          qqt->Resize( newDim );
          aaz->Resize( newDim );
        }
        else {
          // Wrong concrete vector type; throw the old ones away
          delete r_;   delete p_;    delete q_;
          delete z_;   delete qtld_; delete az_;
          allocate = true;
        }
      }
      else {
        allocate = true;
      }

      if ( allocate ) {
        r_    = GenerateVectorObject( sysMat );
        p_    = GenerateVectorObject( sysMat );
        q_    = GenerateVectorObject( sysMat );
        z_    = GenerateVectorObject( sysMat );
        qtld_ = GenerateVectorObject( sysMat );
        az_   = GenerateVectorObject( sysMat );
      }

      vectorLength_ = newDim;
    }

    // Reset all auxiliary vectors to zero
    r_->Init();
    p_->Init();
    q_->Init();
    z_->Init();
    qtld_->Init();
    az_->Init();
  }

  // *********
  //   Solve
  // *********
  template<typename T>
  void COCRSolver<T>::Solve( const BaseMatrix &sysMat,
                             const BaseVector &rhs, BaseVector &sol ) {

    PrivateSetup( sysMat );

    // ------------------------
    //   Read XML parameters
    // ------------------------
    Integer maxIter            = 1000;
    Double  eps                = 1e-6;
    Integer recomputeTrueResEvery = 0;
    bool    consoleConvergence = false;

    if ( xml_ != NULL ) {
      xml_->GetValue( "maxIter", maxIter, ParamNode::INSERT );
      xml_->GetValue( "tol", eps, ParamNode::INSERT );
      xml_->GetValue( "recomputeTrueResEvery", recomputeTrueResEvery, ParamNode::INSERT );
      xml_->GetValue( "consoleConvergence", consoleConvergence, ParamNode::INSERT );
    }

    // ------------------------
    //   Initial residual & threshold (computed once, reused across retries)
    // ------------------------
    sysMat.CompRes( *r_, sol, rhs );
    Double resNorm   = 0;
    Double threshold = ComputeThreshold( eps, rhs, *r_, resNorm, false );

    LOG_DBG(cocrsolver) << "COCR initial ||r||_2 = " << resNorm
                        << ", threshold = " << threshold;

    Integer niter       = 0;
    bool    converged   = ( resNorm <= threshold );
    bool    refreshTriggered = false;
    Double  trueResNorm = resNorm;

    // ----------------------------------------------------------
    //   Outer retry loop. Runs once normally; up to twice when
    //   the preconditioner requests an abort-and-refresh.
    // ----------------------------------------------------------
    bool restartRequested = false;
    do {
      restartRequested = false;
      if ( converged ) break;

      // On retry, refresh r from the current iterate (sol carries
      // partial progress from the aborted attempt).
      if ( refreshTriggered ) {
        sysMat.CompRes( *r_, sol, rhs );
        resNorm = r_->NormL2();
        if ( resNorm <= threshold ) {
          converged = true;
          break;
        }
      }

      // ------------------------
      //   Re-initialise p, q, z with the *current* preconditioner
      // ------------------------
      // p_0 = M^{-1} r_0   (uses freshly refactorised M on retry)
      ptPrecond_->Apply( sysMat, *r_, *p_ );
      // q_0 = A p_0
      sysMat.Mult( *p_, *q_ );
      // z_0 = p_0
      z_->Init();
      z_->Add( *p_ );

      const Double breakdownEps = 1e-30;
      T alpha, beta, rho, dot_rq, dot_zq;

      Vector<T> &rVec    = dynamic_cast<Vector<T>&>( *r_    );
      Vector<T> &qVec    = dynamic_cast<Vector<T>&>( *q_    );
      Vector<T> &qtldVec = dynamic_cast<Vector<T>&>( *qtld_ );
      Vector<T> &azVec   = dynamic_cast<Vector<T>&>( *az_   );

      // ------------------------
      //   Main loop
      // ------------------------
      while ( niter < maxIter ) {
        niter++;

        // qtld = M^{-1} q
        ptPrecond_->Apply( sysMat, *q_, *qtld_ );

        // rho = (qtld, q)_bilinear
        rho = qtldVec * qVec;
        if ( Abs( rho ) < breakdownEps ) {
          WARN( "COCRSolver: breakdown at iteration " << niter
                << " (|rho| = " << Abs( rho )
                << "). Returning best iterate so far." );
          break;
        }

        dot_rq = rVec * qtldVec;
        alpha  = dot_rq / rho;

        // x_k = x_{k-1} + alpha p_{k-1}; r_k = r_{k-1} - alpha q_{k-1}
        sol.Add( alpha, *p_ );
        r_->Add( -alpha, *q_ );

        if ( recomputeTrueResEvery > 0 && ( niter % recomputeTrueResEvery ) == 0 ) {
          sysMat.CompRes( *r_, sol, rhs );
        }

        resNorm = r_->NormL2();
        if ( consoleConvergence ) {
          std::cout << "COCR iter " << niter
                    << "  ||r|| = " << resNorm << std::endl;
        }
        LOG_DBG(cocrsolver) << "COCR iter " << niter
                            << " ||r|| = " << resNorm;

        // Convergence check
        if ( resNorm <= threshold ) {
          converged = true;
          break;
        }

        // ------------------------------------------------
        //   Adaptive abort-and-refresh check
        //   (returns true at most once per Solve)
        // ------------------------------------------------
        if ( !refreshTriggered && ptPrecond_ != nullptr
             && ptPrecond_->ShouldAbortAndRefresh( static_cast<UInt>(niter) ) ) {
          refreshTriggered = true;
          restartRequested = true;
          // Refactorise the preconditioner on the current system matrix.
          // PardisoSolver::Setup will see forceRefresh_ and do a numeric refac.
          ptPrecond_->Setup( const_cast<BaseMatrix&>( sysMat ) );
          break;
        }

        // ------------------------
        //   Rest of CR-style recurrences
        // ------------------------
        z_->Add( -alpha, *qtld_ );
        sysMat.Mult( *z_, *az_ );
        dot_zq = azVec * qtldVec;
        beta   = -dot_zq / rho;

        // p_k = z_k + beta p_{k-1}
        p_->Add( T(1), *z_, beta, *p_ );
        // q_k = az_k + beta q_{k-1}
        q_->Add( T(1), *az_, beta, *q_ );
      } // inner while
    } while ( restartRequested );

    // ------------------------------------------------
    //   False-convergence check at termination
    // ------------------------------------------------
    sysMat.CompRes( *r_, sol, rhs );
    trueResNorm = r_->NormL2();
    if ( converged && trueResNorm > threshold ) {
      WARN( "COCRSolver: false convergence detected.\n"
            << " Recurrence residual norm = " << resNorm << '\n'
            << " True residual norm       = " << trueResNorm );
      converged = false;
    }

    // ------------------------
    //   Solution report
    // ------------------------
    PtrParamNode out = infoNode_->Get( ParamNode::PROCESS )->Get( "solver", ParamNode::APPEND );
    out->Get( "numIter"          )->SetValue( niter );
    out->Get( "finalNorm"        )->SetValue( trueResNorm );
    out->Get( "solutionIsOkay"   )->SetValue( converged );
    out->Get( "refreshTriggered" )->SetValue( refreshTriggered );

    // ------------------------
    //   Running statistics
    // ------------------------
    numCalls_++;
    accIters_ += niter;
    if ( scalFac_ > 0 ) {
      accReduction_ += trueResNorm / scalFac_;
    }
    PtrParamNode stat = infoNode_->Get( ParamNode::SUMMARY )->Get( "statistics" );
    stat->Get( "avgIterations"   )->SetValue( accIters_ / numCalls_ );
    stat->Get( "avgResReduction" )->SetValue( accReduction_ / numCalls_ );
  }

  // Explicit template instantiation
  template class COCRSolver<Double>;
  template class COCRSolver<Complex>;

}