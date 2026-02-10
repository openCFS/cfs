#include "MatVec/generatematvec.hh"

#include "OLAS/precond/BasePrecond.hh"
#include "OLAS/solver/RichardsonSolver.hh"

#include "DataInOut/Logging/LogConfigurator.hh"

namespace CoupledField
{
DEFINE_LOG(richardsonsolver, "richardsonsolver")


  // **************
  //   Destructor
  // **************
  template<typename T>
  RichardsonSolver<T>::~RichardsonSolver() {
    delete r_;
    delete w_;
  }

// ****************
//   Setup method
// ****************
template<typename T>
  void RichardsonSolver<T>::Setup( BaseMatrix &sysmat ) {

  if( r_ != NULL ) {
    delete r_;
    r_ = NULL;
    delete w_;
    w_ = NULL;
  }

  // If not yet done, create auxilliary vectors
  if ( r_ == NULL ) {
    r_ = GenerateVectorObject(sysmat);
    w_ = GenerateVectorObject(sysmat);
  }

}


  // ****************
  //   Solve method
  // ****************
  template<typename T>
  void RichardsonSolver<T>::Solve( const BaseMatrix &sysmat,
				   const BaseVector &rhs, BaseVector &sol ) {

    // Set auxilliary vectors to zero
    r_->Init();
    w_->Init();


    // Variables for loop control
    bool loop = true;
    Integer niter = 0;
    Double resNorm = 0.0;

    // Query parameter object for values
    Integer maxiter = 1; 
    Double tol      = 1e-6;
    Double omegaRe    = 1.0;
    bool consoleConvergence = false;

    // overwrite if set in xml
    if(xml_ != NULL)
    {
      xml_->GetValue("maxIter", maxiter, ParamNode::INSERT);
      xml_->GetValue("tol", tol, ParamNode::INSERT);
      xml_->GetValue("consoleConvergence", consoleConvergence, ParamNode::INSERT);
      xml_->GetValue("omega", omegaRe, ParamNode::INSERT);
    }

    // To call the correct vector methods
    T omega = (T) omegaRe;


    // precond.Apply(sysmat,rhs,sol);

    // Compute residual of initial guess
    sysmat.CompRes( *r_, sol, rhs );

    // compute the euclidean norm of the residual
    resNorm = r_->NormL2();

    LOG_DBG(richardsonsolver) << "Initial residual L2 norm: " << resNorm;

    // Compute preconditioned residual
    ptPrecond_->Apply( sysmat, *r_, *w_ );


    // If Euclidean norm of initial preconditioned residual is too small
    // do not start the loop
    if ( resNorm < tol || resNorm == 0 ) {
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
      resNorm = r_->NormL2();
	  
      if(consoleConvergence == true){
        std::cout<<"Residual L2 norm of iteration "<<niter<<" = "<<resNorm<<std::endl;
      }

      // Check stopping criterion
      if ( resNorm < tol ) {
        loop = false;
      }
      
    } // loop


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
  template class RichardsonSolver<Double>;
  template class RichardsonSolver<Complex>;
  
}
