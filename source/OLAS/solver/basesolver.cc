// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "utils/utils.hh"
#include "basesolver.hh"

namespace OLAS {

  // ******************************
  //   Instantiate Public Methods
  // ******************************
  void BaseSolver::InstantiatePublicMethods( BaseMatrix &sysMat ) {

    ENTER_FCN( "BaseSolver::InstantiatePublicMethods" );

    Error( "This method should _never_ be called!", __FILE__, __LINE__ );

    BasePrecond *dummyPrecond = new IdPrecond;
    BaseVector *dummyVec = GenerateVectorObject( sysMat );

    Setup( sysMat );
    Solve( sysMat, *dummyPrecond, *dummyVec, *dummyVec );
    SolverType myType;
    myType = GetSolverType();

    // Allow derived classes to offer additional public methods
    InstantiateAdditionalPublicMethods( sysMat );
  }


  // *****************************************
  //   Instantiate Additional Public Methods
  // *****************************************
  void BaseSolver::InstantiateAdditionalPublicMethods( BaseMatrix &sysMat ) {

    ENTER_FCN( "BaseSolver::InstantiateAdditionalPublicMethods" );

    Error( "This method should _never_ be called!", __FILE__, __LINE__ );
  }


  // *******************
  //   Log Convergence
  // *******************
  void BaseIterativeSolver::LogConvergence( Double rk, UInt step,
					    bool firstCall ) {

    static double lastNorm = 0;

    ENTER_FCN( "BaseIterativeSolver::LogConvergence" );

    // Write header
    if ( firstCall == true ) {
      (*cla) << "\n "
	     << Enum2String( GetSolverType() )
	     << ": Starting iterations\n\n"
             << "# iter \t resid.norm \trel.res.norm\treduction\n"
             << "#------\t------------\t------------\t---------\n"
             << std::scientific << "  " << step << "\t" << rk << std::endl;
      lastNorm = rk;

      // Test that scalFac_ was initialised
      if ( scalFac_ < 0.0 ) {
	(*error) << "Class attribute scalFac_ was not properly initialised!"
		 << " scalFac_ = " << scalFac_
		 << ".\n Did you call ComputeThreshold?";
	Error( __FILE__, __LINE__ );
      }
    }

    // Write results for current iteration step
    else {
      (*cla) << std::scientific << "  "
             << step << "\t"
             << rk << "\t"
             << rk/scalFac_ << "\t"
             << std::fixed << rk/lastNorm << std::endl;
      cla->unsetf( std::ios::fixed );
      lastNorm = rk;
    }
  }


  // ********************
  //   ComputeThreshold
  // ********************
  Double BaseIterativeSolver::ComputeThreshold( Double eps,
						const BaseVector &rhs,
						const BaseVector &res,
						Double &resNorm,
						bool beVerbose ) {

    ENTER_FCN( "BaseIterativeSolver::ComputeThreshold" );

    // Compute norm of initial residual
    resNorm = res.NormEuclid();

    // Test for the unlikely event, that the inital
    // guess already satisfies the linear system
    if ( resNorm == 0 ) {
      (*warning) << "I like zeros! You too?";
      Warning( __FILE__, __LINE__ );
    }

    // Query user's wish for the stopping criterion
    StopCritType stopCrit = NOSTOPCRITTYPE;
    myParams_->GetEnumValue( "StoppingCriterion", stopCrit );

    // Report this to log file, if required
    if ( beVerbose == true ) {
      (*cla) << "\n Checking stopping rule:\n"
             << " ----------------------\n"
             << " User specified '" << Enum2String( stopCrit ) << "'\n"
	     << std::scientific << " User supplied epsilon = " << eps << '\n';
    }

    switch( stopCrit ) {

      // Now, if the user desires to use an absolute threshold on the
      // Euclidean norm of the residual, we do not modify the tolerance
      // he/she supplied.
    case ABSNORM:
      scalFac_ = 1.0;
      break;

      // User wants to use the norm of the residual relative to that of
      // the right hand side. This can only work, if the latter is non-zero
      // and we are not using the penalty formulation.
      // In case we cannot use RELNORM_RHS we go for RELNORM_RES0 instead.
    case RELNORM_RHS:

      if ( myParams_->GetBoolValue( "RHSwithPenalty" ) == true ) {
	(*cla) << " --> Detected Penalty Formulation\n"
	       << " --> Changing from RELNORM_RHS to RELNORM_RES0\n"
	       << std::endl;
	scalFac_ = resNorm;
	myParams_->SetValue( "StoppingCriterion", RELNORM_RES0 );
	if ( beVerbose == true ) {
	  (*cla) << " Using || r_0 ||_2 = " << scalFac_ << " for scaling\n";
	}
      }
      else {
	scalFac_ = rhs.NormEuclid();
	if ( scalFac_ == 0 ) {
	  (*cla) << " --> Norm of right-hand side is zero\n"
		 << " --> Changing from RELNORM_RHS to RELNORM_RES0\n"
		 << std::endl;
	  scalFac_ = resNorm;
	  myParams_->SetValue( "StoppingCriterion", RELNORM_RES0 );
	  if ( beVerbose == true ) {
	    (*cla) << " Using || r_0 ||_2 = " << scalFac_ << " for scaling\n";
	  }
	}
	else {
	  if ( beVerbose == true ) {
	    (*cla) << " Using || b ||_2 = " << scalFac_ << " for scaling\n";
	  }
	}
      }
      break;

      // User wants to use the norm of the residual relative to that of
      // the initial residual. The latter cannot be zero, since we already
      // tested this. So we can simply go ahead
    case RELNORM_RES0:
      scalFac_ = resNorm;
      if ( beVerbose == true ) {
	(*cla) << " Using || r_0 ||_2 = " << scalFac_ << '\n';
      }
      break;

    default:
      (*error) << "No valid stopping criterion supplied";
      Error( __FILE__, __LINE__ );

    }

    // Now finally we can compute the threshold
    Double threshold = eps * scalFac_;

    // Report threshold to log file, if required
    if ( beVerbose == true ) {
      (*cla) << " Stopping test uses tau = " << threshold << std::endl;
    }

    // That's all
    return threshold;
  }

}
