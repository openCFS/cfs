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
    if ( ptPrecond_ == NULL ) {
      EXCEPTION( "COCRSolver requires a preconditioner. Set <Precond><Id/></Precond> for unpreconditioned operation.");
    }

    // ------------------------
    //   Read XML parameters
    // ------------------------

    // Default values for parameters; will be overwritten if specified in XML
    Integer maxIter            = 2000;
    Double  eps                = 1e-4;
    bool    consoleConvergence = false;

    if ( xml_ != NULL ) {
      xml_->GetValue( "maxIter", maxIter, ParamNode::INSERT );
      xml_->GetValue( "tol", eps, ParamNode::INSERT );
      xml_->GetValue( "consoleConvergence", consoleConvergence, ParamNode::INSERT );
    }

    // ------------------------
    //   Initialisation
    // ------------------------
    // r_0 = b - A x_0
    sysMat.CompRes( *r_, sol, rhs );

    // Stopping threshold via inherited ComputeThreshold(); also fills resNorm
    Double resNorm   = 0;
    Double threshold = ComputeThreshold( eps, rhs, *r_, resNorm, false );
    // Hard floor for breakdown detection on rho = (qtld, q)
    const Double breakdownEps = 1e-30;

    LOG_DBG(cocrsolver) << "COCR initial ||r||_2 = " << resNorm
                        << ", threshold = " << threshold;

    Integer niter       = 0;
    bool    converged   = false;
    Double  trueResNorm = resNorm;

    // Early exit if the initial guess already satisfies the stopping test
    if ( resNorm <= threshold ) {
      converged = true;
    }
    else {

      // p_0 = M^-1 r_0
      ptPrecond_->Apply( sysMat, *r_, *p_ );
      // q_0 = A p_0
      sysMat.Mult( *p_, *q_ );
      // z_0 = p_0   (z is zero from PrivateSetup; one in-place add does it)
      z_->Add( *p_ );

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

        // qtld = M^-1 q
        ptPrecond_->Apply( sysMat, *q_, *qtld_ );

        // rho = (qtld, q)_bilinear
        rho = qtldVec * qVec;
        
        // Breakdown check: if rho is zero (or very close to it), the algorithm breaks down.
        // A breekdown free version of COCR exists and could be implemented in the future
        if ( Abs( rho ) < breakdownEps ) {
          WARN( "COCRSolver: break down at iteration " << niter
                << " (|rho| = " << Abs( rho )
                << "). Returning best solution so far." );
          break;
        }

        // dot_rq = (r, qtld)_bilinear
        dot_rq = rVec * qtldVec;

        // alpha = dot_rq / rho
        alpha = dot_rq / rho;

        // x_k = x_{k-1} + alpha * p_{k-1}
        sol.Add( alpha, *p_ );
        // r_k = r_{k-1} - alpha * q_{k-1}
        r_->Add( -alpha, *q_ );

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

        // z_k = z_{k-1} - alpha * qtld
        z_->Add( -alpha, *qtld_ );

        // az = A z_k
        sysMat.Mult( *z_, *az_ );

        // dot_zq = (az, qtld)_bilinear
        dot_zq = azVec * qtldVec;

        // beta = -dot_zq / rho
        beta = -dot_zq / rho;

        // p_k = 1*z_k + beta * p_{k-1}.  p_ appears on both sides of Add,
        // which is safe here: Vector<T>::Add writes data_[i] only after
        // reading both source operands at index i within the same iteration
        // (see Vector.cc).
        p_->Add( T(1), *z_, beta, *p_ );

        // q_k = az + beta * q_{k-1}
        q_->Add( T(1), *az_, beta, *q_ );
      }

      // ------------------------------------------------
      //   False-convergence check at termination.
      //   The recurrence residual drifts in finite
      //   precision; the only authoritative norm is
      //   the explicitly computed one.
      // ------------------------------------------------
      sysMat.CompRes( *r_, sol, rhs );
      trueResNorm = r_->NormL2();

      if ( converged && trueResNorm > threshold ) {
        WARN( "COCRSolver: false convergence detected.\n"
              << " Recurrence residual norm = " << resNorm << '\n'
              << " True residual norm       = " << trueResNorm );
        converged = false;
      }
    }

    // ------------------------
    //   Solution report
    // ------------------------
    PtrParamNode out = infoNode_->Get( ParamNode::PROCESS )->Get( "solver", ParamNode::APPEND );
    out->Get( "numIter"        )->SetValue( niter );
    out->Get( "finalNorm"      )->SetValue( trueResNorm );
    out->Get( "solutionIsOkay" )->SetValue( converged );

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