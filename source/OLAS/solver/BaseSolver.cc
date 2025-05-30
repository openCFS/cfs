// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "MatVec/BaseVector.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Utils/tools.hh"
#include "Utils/Timer.hh"
#include "Domain/Domain.hh"
#include "OLAS/solver/BaseSolver.hh"

namespace CoupledField {

  static EnumTuple solverTypeTuples[] = 
  {
    EnumTuple( BaseSolver::NOSOLVER, "noSolver" ),
    EnumTuple( BaseSolver::RICHARDSON, "richardson" ),
    EnumTuple( BaseSolver::DIAGSOLVER, "diagsolver"),
    EnumTuple( BaseSolver::CG, "cg"),
    EnumTuple( BaseSolver::GMRES, "gmres" ),
    EnumTuple( BaseSolver::MINRES, "minres" ),
    EnumTuple( BaseSolver::SYMMLQ, "symmlq"),
    EnumTuple( BaseSolver::LAPACK_LU, "lapackLU"),
    EnumTuple( BaseSolver::LAPACK_LL, "lapackLL" ),
    EnumTuple( BaseSolver::LU_SOLVER, "directLU" ),
    EnumTuple( BaseSolver::LDL_SOLVER, "directLDL"),
    EnumTuple( BaseSolver::LDL_SOLVER2, "directLDL2"),
    EnumTuple( BaseSolver::PARDISO_SOLVER, "pardiso" ),
    EnumTuple( BaseSolver::UMFPACK, "umfpack" ),
    EnumTuple( BaseSolver::CHOLMOD, "cholmod"),
    EnumTuple( BaseSolver::LIS, "lis"),
    EnumTuple( BaseSolver::GINKGO, "ginkgo"),
    EnumTuple( BaseSolver::PETSC, "petsc"),
    EnumTuple( BaseSolver::SUPERLU, "superlu" ),
    EnumTuple( BaseSolver::PHIST, "phist_linSolv"),
    EnumTuple( BaseSolver::EXTERNAL_SOLVER, "externalSolver")

  };

  Enum<BaseSolver::SolverType> BaseSolver::solverType = \
  Enum<BaseSolver::SolverType>("Solver Types",
      sizeof(solverTypeTuples) / sizeof(EnumTuple),
      solverTypeTuples); 

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
      WARN("I like zeros! You too?");
    }

    // Query user's wish for the stopping criterion
    StopCritType stopCrit = NOSTOPCRITTYPE;
    std::string stopCritStr = "relNormRes0";
    PtrParamNode stopRuleNode = xml_->Get("stoppingRule", ParamNode::INSERT );
    stopRuleNode->GetValue("type", stopCritStr, ParamNode::INSERT);
    String2Enum( stopCritStr, stopCrit );
    
    // Report this to log file, if required
    if ( beVerbose == true ) {
      std::string tmp;
      Enum2String( stopCrit, tmp );
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

      if ( usingPenalty_ ) {
        scalFac_ = resNorm;

        stopRuleNode->Get("type")->SetValue("relNormRes0");
        if ( beVerbose == true ) { } // removed logging
      }
      else {
        scalFac_ = rhs.NormL2();
        if ( scalFac_ == 0 ) {
          scalFac_ = resNorm;
        stopRuleNode->Get("type")->SetValue("relNormRes0");
        }
      }
      break;

      // User wants to use the norm of the residual relative to that of
      // the initial residual. The latter cannot be zero, since we already
      // tested this. So we can simply go ahead
    case RELNORM_RES0:
      scalFac_ = resNorm;
      break;

    default:
      EXCEPTION( "No valid stopping criterion supplied" );
    }
    // Now finally we can compute the threshold
    return eps * scalFac_;
  }
  
  BaseSolver::~BaseSolver()
  { }
  
  void BaseSolver::PostInit()
  {
    // Assert that info Node is set
    assert( infoNode_ );

    setupTimer_ = boost::shared_ptr<Timer>(new Timer("setup_" + solverType.ToString(GetSolverType())));
    infoNode_->Get(ParamNode::SUMMARY)->Get("setup/timer")->SetValue(setupTimer_);

    solveTimer_ = boost::shared_ptr<Timer>(new Timer("solve_" + solverType.ToString(GetSolverType())));
    infoNode_->Get(ParamNode::SUMMARY)->Get("solve/timer")->SetValue(solveTimer_);
  }
  
  void BaseSolver::SetPrecond( BasePrecond* precond ) {
    ptPrecond_ = precond;
  }
  
   void BaseSolver::Apply(const BaseMatrix& sysmat, const BaseVector& r, 
                          BaseVector& z) {
     
     this->Solve(sysmat, r, z);
   }
   
  void BaseSolver::CheckParameter(PtrParamNode out, double* val, const char* param_name)
  {
    PtrParamNode tmp = out->Get(param_name);
    tmp->Get("default")->SetValue(*val);
    if (xml_ != NULL && xml_->Has(param_name))
    {
      *val = xml_->Get(param_name)->As<Double>();
      tmp->Get("set")->SetValue(*val);
    }
  }

  void BaseSolver::CheckParameter(PtrParamNode out, int* val, const char* param_name)
  {
    PtrParamNode tmp = out->Get(param_name);
    tmp->Get("default")->SetValue(*val);
    if (xml_ != NULL && xml_->Has(param_name))
    {
      *val = xml_->Get(param_name)->As<Integer>();
      tmp->Get("set")->SetValue(*val);
    }
  }

  void BaseSolver::CheckParameter(PtrParamNode out, size_t* val, const char* param_name)
  {
    PtrParamNode tmp = out->Get(param_name);
    tmp->Get("default")->SetValue(*val);
    if (xml_ != NULL && xml_->Has(param_name))
    {
      *val = xml_->Get(param_name)->As<Integer>();
      tmp->Get("set")->SetValue(*val);
    }
  }

  void BaseSolver::CheckParameter(PtrParamNode out, bool* val, const char* param_name)
  {
    // by convention we interpret this as "integer"
    int* int_ptr = reinterpret_cast<int*>(val);
    
    PtrParamNode tmp = out->Get(param_name);
    tmp->Get("default")->SetValue(*val);
    if (xml_ != NULL && xml_->Has(param_name))
    {
      *int_ptr = xml_->Get(param_name)->As<bool>() == false ? 0 : 1;
      tmp->Get("set")->SetValue(*int_ptr == 0 ? false : true);
    }
  }
}
