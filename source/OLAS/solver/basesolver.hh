// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef OLAS_BASESOLVER_HH
#define OLAS_BASESOLVER_HH

#include "General/environment.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField {

  class BasePrecond;
  class BaseMatrix;
  class BaseVector;
  class OLAS_Report;
  class ParamNode;
  class Timer;

  // =========================================================================
  // BASE SOLVER
  // =========================================================================


  //! Base class for algebraic system solver
  class BaseSolver {
    
  public:
    //! Type of solver

    //! This enumeration data type describes the type of solver which
    //! is applied to solve the system. Note that not all solvers can handle
    //! all types of matrices. The enumeration contains the following values:
    //! - NOSOLVER
    //! - DIRECT
    //! - RICHARDSON
    //! - CG
    //! - LANCZOS
    //! - QMR
    //! - GMRES
    //! - MINRES
    //! - LAPACK_LU
    //! - LAPACK_LL
    //! - PARDISO
    //! - ILUPACK_SOLVER
    //! - LU_SOLVER
    //! - LDL_SOLVER
    //! - LDL_SOLVER2
    typedef enum {NOSOLVER, DIRECT, RICHARDSON, CG, LANCZOS, QMR, GMRES,
                  MINRES, SYMMLQ, LAPACK_LU, LAPACK_LL, PARDISO,  ILUPACK, LU_SOLVER, CHOLMOD, 
                  LDL_SOLVER, LDL_SOLVER2, DIAGSOLVER } SolverType;    
    static Enum<SolverType> solverType;

  public:

    //! Default Constructor
    BaseSolver() {
      setupTimer_ = NULL;
      solveTimer_ = NULL;
      usingPenalty_ = false;
    }

    //! Default Destructor
    virtual ~BaseSolver() {
    }

    /** Does constructor stuff only possible after child constructors are called */
    void PostInit();

    //! General setup routine
  
    /** For direct solvers this might involve a factorization, for Krylov
     * solvers just construction of some vectors etc.
     * @param analysis_id references to the info/analysis/progress/step(/substep) 
     *                   element with the "analysis_id" attribute */
    virtual void Setup( BaseMatrix &sysmat, PtrParamNode analysis_id = PtrParamNode()) = 0;

    /** Solve the linear system sysmat*sol=rhs for sol
     * @param analysis_id @see Setup() */
    virtual void Solve( const BaseMatrix &sysmat, const BasePrecond &precond,
			const BaseVector &rhs, BaseVector &sol, PtrParamNode analysis_id = PtrParamNode()) = 0;

    //! Query type of the solver

    //! This method can be used to query the type of this solver. The answer
    //! is encoded as a value of the enumeration data type SolverType.
    //! \return type of the solver
    virtual SolverType GetSolverType() = 0;

    /** Gives the timer located within PtrParamNode */
    Timer* GetSetupTimer() { return setupTimer_; }
    Timer* GetSolveTimer() { return solveTimer_; }

    void SetUsingPenalty(bool usingPenalty) {
      usingPenalty_ = usingPenalty;
    }
    
  protected:
    
    /** Helper method for InitParameters. Takes the default-value, prints it to the ParamNode (via
     * param_name).
     * If there is a parameter with param_name in the xml file, it is used and printed.
     * @param param_name simple xpath chain via slash */
    void CheckParameter(PtrParamNode out, double* val, const char* param_name);
    void CheckParameter(PtrParamNode out, char** val, const char* param_name);
    void CheckParameter(PtrParamNode out, int* val, const char* param_name);  
    void CheckParameter(PtrParamNode out, size_t* val, const char* param_name);  
    void CheckParameter(PtrParamNode out, bool* val, const char* param_name);
    
    template<typename E> // needed in header because of instantiation
    void CheckParameter(PtrParamNode out, Enum<E>* en, int* val, const char* param_name)
    {
      PtrParamNode tmp = out->Get(param_name);
      tmp->Get("default")->SetValue(en->ToString(static_cast<E>(*val)));
      if (xml_ && xml_->Has(param_name))
      {
        *val = en->Parse(xml_->Get(param_name)->As<std::string>());
        tmp->Get("set")->SetValue(en->ToString(static_cast<E>(*val)));
      }
    }
    
    /** This is the description of the solver part in XML
     * This may not be NULL since proper default values have to be set
     * before the solver gets constructed. */
    PtrParamNode xml_;
    
    /** This stores the general solver Information (HEADER, SUMMARY) ->
     * For the current solve steps the pointer is given */
    PtrParamNode solverInfo_;

    /** This is a pointer to the setup timer. Located within PtrParamNode*/
    Timer* setupTimer_;
    Timer* solveTimer_;

    bool usingPenalty_;
  };


  // =========================================================================
  // BASE ITERATIVE SOLVER
  // =========================================================================


  //! Base class for iterative linear system solvers in OLAS
  class BaseIterativeSolver : public BaseSolver {

  public:

    //! Default Constructor
    BaseIterativeSolver() {
      scalFac_ = -1.0;
    }

    //! Default Destructor
    ~BaseIterativeSolver() {
    }

  protected:

    //! Auxilliary routine for logging convergence info

    //! This routine will log the Euclidean norm of the residual of the
    //! current approximation and the relative residual to the standard
    //! %OLAS report stream (*cla). Relative residual here means
    //! \f$\|r\|_2/\|b\|_2\f$, where \f$b\f$ is the right hand side of
    //! the linear system or the residual of the initial guess, if the right
    //! hand side is zero. If the solver allows for left-preconditioning,
    //! then the residual might be the preconditioned residual.
    //! \param rk        norm of residual of the current approximate,
    //!                  i.e. \f$\|r^{(k)}\|_2\f$
    //! \param step      number \f$k\f$ of current iteration step1
    //! \param firstCall on first call we also write a header
    void LogConvergence( Double rk, UInt step, bool firstCall = false );

    //! Compute threshold for stopping criterion

    //! This method will compute the threshold for the stopping criterion.
    //! We use the following test to terminate the iteration. An approximate
    //! solution \f$x^{(k)}\f$ is accepted as good enough, if
    //! \f[
    //! \| r^{(k)} \|_2 \leq \mbox{tol} := \varepsilon \cdot \tau
    //! \f]
    //! where the parameter \f$\tau\f$ depends on the chosen stopping rule.
    //! The latter is determined from the <b>StoppingCriterion</b> setting
    //! in the myParams_ object. We use
    //! \f[
    //! \tau = \left\{ \begin{array}{cl}
    //! 1 & \mbox{ for ABSNORM}\\[2ex]
    //! \|b\|_2 & \mbox{ for RELNORM\_RHS (and $b\neq 0$)}\\[2ex]
    //! \|r^{(0)}\|_2 & \mbox{ for RELNORM\_RES0 or for RELNORM\_RHS if $b=0$
    //!  or penalty formulation is used}
    //! \end{array}\right.
    //! \f]
    //! Here \f$\varepsilon\f$ is the user supplied value for the stopping
    //! test. The method queries the <b>RHSwithPenalty</b> steering
    //! parameter to determine whether the penalty formulation is active.
    //! \param eps       the user supplied stopping value
    //! \param rhs       right-hand side \f$b\f$ of linear system
    //! \param res       residual \f$r^{(0)}\f$ of initial guess
    //! \param resNorm   the method returns in this parameter the norm of
    //!                  the initial residual
    //! \param beVerbose this flag controls the verbosity of the method
    //! \return The threshold tol to be used in the stopping test
    virtual Double ComputeThreshold( Double eps, const BaseVector &rhs,
				     const BaseVector &res, Double &resNorm,
				     bool beVerbose );

    //! Scaling factor used to compute the threshold

    //! This attribute stores the scaling factor \f$\| v \|_2\f$ computed
    //! by ComputeThreshold in order to obtain the threshold for the stopping
    //! criterion. We store this, so that we can use it e.g. in the
    //! LogConvergence method.
    Double scalFac_;
  };


  // =========================================================================
  // BASE DIRECT SOLVER
  // =========================================================================


  //! Base class for direct linear system solvers in OLAS
  class BaseDirectSolver : public BaseSolver {

  public:

    //! Default Constructor
    BaseDirectSolver() {
    }

    //! Default Destructor
    ~BaseDirectSolver() {
    }

  };


} // namespace



#endif // OLAS_BASESOLVER_HH
