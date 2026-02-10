// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "MatVec/opdefs.hh"
#include "MatVec/generatematvec.hh"

#include "OLAS/precond/BasePrecond.hh"
#include "OLAS/solver/CGSolver.hh"

namespace CoupledField {


  // **************
  //   Destructor
  // **************
  template<typename T>
  CGSolver<T>::~CGSolver() {
    delete r_;
    delete s_;
    delete d_;
    delete q_;
  }

  // ****************
  //   Setup method
  // ****************
  template<typename T>
    void CGSolver<T>::Setup( BaseMatrix &sysmat ) {
    
    if( r_ != NULL ) {
      delete r_;
      r_ = NULL;
      delete d_;
      d_ = NULL;
      delete s_;
      s_ = NULL;
      delete q_;
      q_ = NULL;
    }
    
    // If not yet done, create auxilliary vectors
    if ( r_ == NULL ) {
      r_ = GenerateVectorObject(sysmat);
      d_ = GenerateVectorObject(sysmat);
      s_ = GenerateVectorObject(sysmat);
      q_ = GenerateVectorObject(sysmat);
    }    
    
  }
  
  // ****************
  //   Solve method
  // ****************
  template<typename T>
  void CGSolver<T>::Solve( const BaseMatrix &sysmat,
                           const BaseVector &rhs, BaseVector &sol ) {

    // Tracing information

    // Set auxilliary vectors to zero
    r_->Init();
    d_->Init();
    s_->Init();
    q_->Init();

    // We need some scalar variables
    T delta_new = 0;
    T delta_old = 0;
    T alpha = 0;
    T beta = 0;
    T aux = 0;

    // Variables for loop control and related stuff
    bool loop = true;
    Integer niter = 0;
    Double resNorm  = 0;
    Double tol      = 0;

    // set defaults:
    int    maxiter = 50;
    double eps     = 1e-6;
    int    tmp     = 50; // resDirectly
    bool consoleConvergence = false;

    // overwrite if set in xml
    if(xml_ != NULL)
    {
      xml_->GetValue("maxIter", maxiter, ParamNode::INSERT);
      xml_->GetValue("tol", eps, ParamNode::INSERT);
      xml_->GetValue("resDirectly", tmp, ParamNode::INSERT);
      xml_->GetValue("consoleConvergence", consoleConvergence, ParamNode::INSERT);
    } 
    if ( tmp <= 0 ) {
      EXCEPTION( "CGSolver::CGSolver: The current value of "
               << "CG_resDirectly = " << aux
               << "! Please choose a positive value!" );
    }
    else {
      resDirectly_ = (UInt)tmp;
    }

    // =====================
    //   Setup phase of CG
    // =====================

    // Compute residual of initial guess
    sysmat.CompRes( *r_, sol, rhs );

    // Compute threshold for stopping test
    tol = ComputeThreshold( eps, rhs, *r_, resNorm, false );

    // Compute d by applying preconditioner
    ptPrecond_->Apply( sysmat, *r_, *d_ );

    // Compute new delta as inner product of r and d
    r_->Inner( *d_, delta_new );


    // Compute Euclidean norm of preconditioned initial residual


    // If Euclidean norm of initial preconditioned residual is too small
    // do not start CG loop
    if ( resNorm < tol || resNorm == 0 ) {
      loop = false;
    }


    // ====================
    //   Loop Phase of CG
    // ====================
    while( (niter < maxiter) && (loop == true) ) {

      // We start a new iteration
      niter++;

      // Determine new q <- A*d
      
      sysmat.Mult( *d_, *q_ );

      // Compute the new parameter alpha <- delta_new / (d^T * q)
      d_->Inner( *q_, aux );
      alpha = delta_new / aux;

      // Determine new approximation x <- x + alpha * d
      sol.Add( alpha, *d_ );

      // Compute residual of new approximation, either by recursion or
      // directly from the definition r = b - Ax
      if ( (niter % resDirectly_) != 0 ) {
        // We use recursion r <- r - alpha * q
        r_->Add( -alpha, *q_ );
      }
      else {
        // We use the definition
        sysmat.CompRes( *r_, sol, rhs );
      }

      // Determine norm of new residual
      // and log progress, if required
      resNorm = r_->NormL2();

      if(consoleConvergence == true){
        std::cout<<"Residual L2 norm of iteration "<<niter<<" = "<<resNorm<<std::endl;
      }

      // Compute s = M^-1*r by applying preconditioner
      ptPrecond_->Apply( sysmat, *r_, *s_ );

      // Save old delta and compute new one
      delta_old = delta_new;
      r_->Inner( *s_, delta_new );

      // Compute new beta and then new search direction d <- s + beta * d
      beta = delta_new / delta_old;
      d_->Axpy( beta, *s_ );

      // Log progress

      // Check stopping criterion
      if ( resNorm < tol ) {
        loop = false;
      }
    }

    // ============================
    //   Generate solution report
    // ============================
    Double reduction = resNorm / scalFac_;
    PtrParamNode out = infoNode_->Get(ParamNode::PROCESS)->Get("solver", ParamNode::APPEND);
    out->Get("numIter")->SetValue(niter);
    out->Get("finalNorm")->SetValue(resNorm);
    if ( loop == false ) {
      out->Get("solutionIsOkay")->SetValue(true);
    }
    else {
      out->Get("solutionIsOkay")->SetValue(std::any(false));
    }
    
    // Calculate average number of iterations and residual error reduction
    numCalls_++;
    accIters_ += niter;
    accReduction_ += reduction;
    
    PtrParamNode stat = infoNode_->Get(ParamNode::SUMMARY)->Get("statistics");
    stat->Get("avgIterations")->SetValue(accIters_ / numCalls_);
    stat->Get("avgResReduction")->SetValue( accReduction_ / numCalls_);
  }

  // Explicit template instantiation
  template class CGSolver<Double>;
  template class CGSolver<Complex>;
  
}
