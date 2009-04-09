// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "MatVec/basevector.hh"
#include "OLAS/algsys/olasparams.hh"

#include "Utils/tools.hh"

#include "OLAS/solver/basesolver.hh"

namespace CoupledField {

  // *******************
  //   Log Convergence
  // *******************
  void BaseIterativeSolver::LogConvergence( Double rk, UInt step,
                                            bool firstCall ) {

    static double lastNorm = 0;


    // Write header
    if ( firstCall == true ) {
      std::string tmp;
      Enum2String( GetSolverType(), tmp );

      (*cla) << "\n "
      << tmp
      << ": Starting iterations\n\n"
             << "# iter \t resid.norm \trel.res.norm\treduction\n"
             << "#------\t------------\t------------\t---------\n"
             << std::scientific << "  " << step << "\t" << rk << std::endl;
      lastNorm = rk;

      // Test that scalFac_ was initialised
      if ( scalFac_ < 0.0 ) {
        EXCEPTION( "Class attribute scalFac_ was not properly initialised!"
            << " scalFac_ = " << scalFac_
            << ".\n Did you call ComputeThreshold?" );
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


    // Compute norm of initial residual
    resNorm = res.NormL2();

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
      std::string tmp;
      Enum2String( stopCrit, tmp );

      (*cla) << "\n Checking stopping rule:\n"
             << " ----------------------\n"
             << " User specified '" << tmp << "'\n"
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
        scalFac_ = rhs.NormL2();
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
      EXCEPTION( "No valid stopping criterion supplied" );

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
