#include "MatVec/generatematvec.hh"

#include "OLAS/precond/BasePrecond.hh"
#include "OLAS/solver/RichardsonSolver.hh"

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
				   const BaseVector &rhs, BaseVector &sol ) {

    EXCEPTION("The Richardson solver has not been in use for a very long time."
              << "Please check if it is still working for you!");
    
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

    // Query parameter object for values
    Integer maxiter = 1; 
    xml_->GetValue("maxIter", maxiter, ParamNode::INSERT);
    Double eps      = 1e-6;
    xml_->GetValue("tol", eps, ParamNode::INSERT);
    Double epsmach  = 1e-20;
    xml_->GetValue("epsmach", epsmach, ParamNode::INSERT);
    Double omega    = 1.0;
    xml_->GetValue("omega", omega, ParamNode::INSERT);

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
    ptPrecond_->Apply( sysmat, *r_, *w_ );


    // TEST: Due to the use of the penalty method we currently
    // follow the same approach as in LAS and use a relative
    // stopping criterion
    Double tol = norm_new * eps;

    //check if we have to start at all
    if ( norm_new < eps || norm_new < epsmach ) {
      loop = false;
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
      ptPrecond_->Apply( sysmat, *r_, *w_ );

      // compute the euclidean norm of the residual
      norm_new = r_->NormL2();
	  

      // Check stopping criterion
      if ( norm_new < tol || norm_new < epsmach ) {
        loop = false;
      }
      
    }//iter


#ifdef DEBUG_RICHARDSON
    (*debug) << " ------- END RICHARDSON ITERATION -------- " << std::endl;
#endif

    // ****************************
    //   Generate solution report
    // ****************************
    PtrParamNode out = infoNode_->Get(ParamNode::PROCESS)->Get("solver", ParamNode::APPEND);
    out->Get("numIter")->SetValue(niter);
    out->Get("finalPrecondResNorm")->SetValue(norm_new);
    
  }

// Explicit template instantiation
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template class RichardsonSolver<Double>;
  template class RichardsonSolver<Complex>;
#endif
  
}
