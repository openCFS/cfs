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

    r_    = NULL;
    p_    = NULL;
    q_    = NULL;
    z_    = NULL;
    qtld_ = NULL;
    az_   = NULL;
    vectorLength_ = 0;
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


  // ******************************************************
  //   Bilinear dot product (NOT the Hermitian inner prod)
  // ******************************************************
  template<typename T>
  T COCRSolver<T>::BilinearInner( const BaseVector &x,
                                  const BaseVector &y ) const {

    // Vector<T>::operator* performs the bilinear product x^T y without
    // complex conjugation (see Vector.hh: "Vector product, but not the
    // scalar product, which would use the complex conjugate as Inner()
    // does!"). This is exactly the dot product COCR requires.
    const Vector<T> &xt = dynamic_cast<const Vector<T>&>( x );
    const Vector<T> &yt = dynamic_cast<const Vector<T>&>( y );
    return xt * yt;
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
    Integer resDirectly        = 0;       // 0 = never refresh true residual
    bool    consoleConvergence = false;

    if ( xml_ != NULL ) {
      xml_->GetValue( "maxIter",            maxIter,            ParamNode::INSERT );
      xml_->GetValue( "tol",                eps,                ParamNode::INSERT );
      xml_->GetValue( "resDirectly",        resDirectly,        ParamNode::INSERT );
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

      // Hard floor for breakdown detection on rho = (qtld, q)
      const Double breakdownEps = 1e-30;

      T alpha, beta, rho, dot_rq, dot_zq;

      // ------------------------
      //   Main loop
      // ------------------------
      while ( niter < maxIter ) {

        niter++;

        // qtld = M^-1 q
        ptPrecond_->Apply( sysMat, *q_, *qtld_ );

        // rho = (qtld, q)_bilinear
        rho = BilinearInner( *qtld_, *q_ );

        if ( Abs( rho ) < breakdownEps ) {
          WARN( "COCRSolver: breakdown at iteration " << niter
                << " (|rho| = " << Abs( rho )
                << "). Returning best iterate so far." );
          break;
        }

        // dot_rq = (r, qtld)_bilinear
        dot_rq = BilinearInner( *r_, *qtld_ );

        // alpha = dot_rq / rho
        alpha = dot_rq / rho;

        // x_k = x_{k-1} + alpha * p_{k-1}
        sol.Add( alpha, *p_ );
        // r_k = r_{k-1} - alpha * q_{k-1}
        r_->Add( -alpha, *q_ );

        // Periodic true-residual refresh for monitoring (CG-style).
        // We replace the recurrence r by b - A x; the z/p/q recurrence
        // continues unchanged, which matches the established CG pattern.
        if ( resDirectly > 0 && ( niter % resDirectly ) == 0 ) {
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

        // z_k = z_{k-1} - alpha * qtld
        z_->Add( -alpha, *qtld_ );

        // az = A z_k
        sysMat.Mult( *z_, *az_ );

        // dot_zq = (az, qtld)_bilinear
        dot_zq = BilinearInner( *az_, *qtld_ );

        // beta = -dot_zq / rho
        beta = -dot_zq / rho;

        // p_k = z_k + beta * p_{k-1}
        // Self-aliasing is safe for the 4-arg Add: data_[i] is written
        // AFTER both reads idvec1[i] (= z) and idvec2[i] (= p) within the
        // same loop iteration (see Vector.cc Vector<T>::Add).
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
    PtrParamNode out = infoNode_->Get( ParamNode::PROCESS )
                                ->Get( "solver", ParamNode::APPEND );
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

    PtrParamNode stat = infoNode_->Get( ParamNode::SUMMARY )
                                 ->Get( "statistics" );
    stat->Get( "avgIterations"   )->SetValue( accIters_ / numCalls_ );
    stat->Get( "avgResReduction" )->SetValue( accReduction_ / numCalls_ );
  }


  // Explicit template instantiation
  template class COCRSolver<Double>;
  template class COCRSolver<Complex>;

}
