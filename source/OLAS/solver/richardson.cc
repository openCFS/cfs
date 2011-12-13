// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "General/environment.hh"
#include "General/exception.hh"
#include "MatVec/basematrix.hh"
#include "MatVec/basevector.hh"
#include "MatVec/generatematvec.hh"
#include "OLAS/precond/baseprecond.hh"
#include "OLAS/solver/richardson.hh"

namespace CoupledField {


  // **************
  //   Destructor
  // **************
  template<typename T>
  RichardsonSolver<T>::~RichardsonSolver() {
    delete r_;
    delete w_;
  }


  // ****************
  //   Solve method
  // ****************
  template<typename T>
  void RichardsonSolver<T>::Solve( const BaseMatrix &sysmat,
				   const BasePrecond &precond,
				   const BaseVector &rhs, BaseVector &sol, PtrParamNode analysis_step ) {

    EXCEPTION("The Richardson solver has not been in use for a very long time."
              << "Please check if it is still working for you!");
    
    // Tracing information
    (*cla) << "### preconditioned Richardson Solver" << std::endl;

    // If not yet done, create auxilliary vectors
    if ( r_ == NULL ) {
      r_ = GenerateVectorObject(sysmat);
      w_ = GenerateVectorObject(sysmat);
    }

    // Set auxilliary vectors to zero
    r_->Init();
    w_->Init();


    // Variables for loop control
    bool loop = true;
    Integer niter = 0;
    Double norm_new;
    Double norm_old;

    // Query parameter object for values
    PtrParamNode pn = xml_->Get("solver", ParamNode::INSERT); 
    pn = pn->Get("richardson", ParamNode::INSERT); 
    Integer maxiter = 1; 
    pn->GetValue("maxIter", maxiter, ParamNode::INSERT);
    Double eps      = 1e-6;
    pn->GetValue("tol", eps, ParamNode::INSERT);
    Double epsmach  = 1e-20;
    pn->GetValue("epsmach", epsmach, ParamNode::INSERT);
    Double omega    = 1.0;
    pn->GetValue("omega", omega, ParamNode::INSERT);

#ifdef DEBUG_RICHARDSON
    (*debug) << " ------- START RICHARDSON ITERATION -------- " << std::endl;
#endif

    // precond.Apply(sysmat,rhs,sol);

    // Compute residual of initial guess
    sysmat.CompRes( *r_, sol, rhs );

#ifdef DEBUG_RICHARDSON
    (*debug) << " initial residual: " << std::endl;
    r_->Print(*debug);
#endif

    // compute the euclidean norm of the residual
    norm_new = r_->NormL2();

#ifdef DEBUG_RICHARDSON
    (*debug) << " initial residual norm: " << norm_new << std::endl;
#endif

    // Compute preconditioned residual
    precond.Apply( sysmat, *r_, *w_ );


    // TEST: Due to the use of the penalty method we currently
    // follow the same approach as in LAS and use a relative
    // stopping criterion
    Double tol = norm_new * eps;

    Double relativeNorm = 1.0;
    Double initialNorm = norm_new;

    norm_old = norm_new; //just for output

    // Log progress
    (*cla) << " Iteration " << niter << ": res norm = " << norm_new << " ,";
    (*cla)   << " rel. to last: " << norm_new / norm_old << 
      "rel. to first: " << relativeNorm << std::endl;


    //check if we have to start at all
    if ( norm_new < eps || norm_new < epsmach ) {
      loop = false;
      (*cla) << "### Norm is small enough, we do not start PRichardson"
	     << std::endl;
    }


    // ====================
    //   Loop Phase of Richardson
    // ====================
    while( (niter < maxiter) && (loop == true) ) {

      // We start a new iteration
      niter++;

      // Determine new approximation x <- x + omg * (P^-1 r)
      sol.Add( omega, *w_ );

      // Compute residual of new approximation
      sysmat.CompRes( *r_, sol, rhs );

      // Compute w = P^-1*r by applying preconditioner
      precond.Apply( sysmat, *r_, *w_ );

      norm_old = norm_new;

      // compute the euclidean norm of the residual
      norm_new = r_->NormL2();
	  

      // Check stopping criterion
      if ( norm_new < tol || norm_new < epsmach ) {
	if ( norm_new < tol ) {
	  (*cla) << std::endl
		 << "### Terminating iterations since norm < eps = " << eps
		 << std::endl;
	}
	else {
	  (*cla) << std::endl
		 << "### Terminating iterations since norm < epsmach = "
		 << epsmach << std::endl;
	}
	loop = false;
      }
      
      relativeNorm = norm_new / initialNorm;
      
      // Log progress
      (*cla) << " Iteration " << niter << ": res norm = " << norm_new << " ,";
      (*cla)   << " rel. to last: " << norm_new / norm_old << 
	"rel. to first: " << relativeNorm << std::endl;

    }//iter


#ifdef DEBUG_RICHARDSON
    (*debug) << " ------- END RICHARDSON ITERATION -------- " << std::endl;
#endif

    // ****************************
    //   Generate solution report
    // ****************************
    PtrParamNode out = solverInfo_->Get(ParamNode::PROCESS)->Get("solver", ParamNode::APPEND);
    out->Get("numIter")->SetValue(niter);
    out->Get("finalPrecondResNorm")->SetValue(norm_new);
    
  }

// Explicit template instantiation
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template class RichardsonSolver<Double>;
  template class RichardsonSolver<Complex>;
#endif
  
}
