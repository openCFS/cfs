#ifndef OLAS_BASESOLVER_HH
#define OLAS_BASESOLVER_HH

#include <set>
#include "General/Environment.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "OLAS/precond/BasePrecond.hh"
#include "MatVec/BaseMatrix.hh"

namespace CoupledField {

  class BaseMatrix;
  class BaseVector;
  class OLAS_Report;
  
  class Timer;

  // =========================================================================
  // BASE SOLVER
  // =========================================================================


  //! Base class for algebraic system solver
  class BaseSolver : public BasePrecond {
    
  public:
    //! Type of solver

    //! This enumeration data type describes the type of solver which
    //! is applied to solve the system. Note that not all solvers can handle
    //! all types of matrices.
    typedef enum {NOSOLVER, RICHARDSON, CG, LANCZOS, QMR, GMRES,
                  MINRES, SYMMLQ, LAPACK_LU, LAPACK_LL, PARDISO_SOLVER,
                  UMFPACK, CHOLMOD, LIS,PETSC, SUPERLU, SPOOLES,
                  LDL_SOLVER, LU_SOLVER, DIAGSOLVER ,PHIST, EXTERNAL_SOLVER} SolverType;
    static Enum<SolverType> solverType;

  public:

    //! Default Constructor
    BaseSolver() {
      usingPenalty_ = false;
      ptPrecond_ = NULL;
    }

    //! Default Destructor
    virtual ~BaseSolver();

    /** Does constructor stuff only possible after child constructors are called */
    void PostInit();

    //! Set preconditioner object 
    void SetPrecond( BasePrecond* precond);
    
    //! Get preconditioner object
    BasePrecond* GetPrecond() {
      return ptPrecond_;
    }
    
    //! Notify the solver that a new matrix pattern has been set
    // If we use the GetRidOfZeros functionality, the solver will be set up before we can check 
    // the actual matrix for zero entries. In this case (and if we actually find zeros that 
    // shall be removed), we have to notify the solver that the matrix and therefore the pattern 
    // has changed, which re-triggers e.g. symbolic factorization.
    virtual void SetNewMatrixPattern() = 0;
    
    //! \see BasePrecond::Setup
    virtual void Setup(BaseMatrix& sysmat) = 0;
    
  
    //! Solve the linear system sysmat*sol=rhs for sol
     //! \param analysis_id @see Setup() */
    virtual void Solve( const BaseMatrix &sysmat, const BaseVector &rhs, BaseVector &sol) = 0;

    //! Apply the solver as preconditioner
    //! \see BasePrecond::Apply
    void Apply(const BaseMatrix& sysmat, const BaseVector& r,  BaseVector& z);
    
    //! Query type of the solver

    //! This method can be used to query the type of this solver. The answer
    //! is encoded as a value of the enumeration data type SolverType.
    //! \return type of the solver
    virtual SolverType GetSolverType() = 0;

    /** Gives the timer located within PtrParamNode */
    boost::shared_ptr<Timer> GetSetupTimer() { return setupTimer_; }
    boost::shared_ptr<Timer> GetSolveTimer() { return solveTimer_; }

    void SetUsingPenalty(bool usingPenalty) {
      usingPenalty_ = usingPenalty;
    }
    
  protected:
    
    /** Helper method for InitParameters. Takes the default-value, prints it to the ParamNode (via
     * param_name).
     * If there is a parameter with param_name in the xml file, it is used and printed.
     * @param param_name simple xpath chain via slash */
    void CheckParameter(PtrParamNode out, double* val, const char* param_name);
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
        PtrParamNode tmp2 = xml_->Get(param_name);
        *val = en->Parse( tmp2->As<std::string>() );
        tmp->Get("set")->SetValue(en->ToString(static_cast<E>(*val)));
      }
    }
    
    //! Pointer to preconditioner object
    BasePrecond* ptPrecond_;
    
    /** This is a pointer to the setup timer. Located within PtrParamNode*/
    boost::shared_ptr<Timer> setupTimer_;
    boost::shared_ptr<Timer> solveTimer_;

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
      numCalls_ = 0;
      accIters_ = 0;
      accReduction_ = 0.0;
    }

    //! Default Destructor
    ~BaseIterativeSolver() {
    }

  protected:

    //! Auxilliary routine for logging convergence info

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
    
    //! Number of times the solver was called
    UInt numCalls_;
    
    //! Accumulated number of solver iterations
    UInt accIters_;
    
    //! Accumulated residual error reduction
    Double accReduction_;
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

    //! Default Destructor. 
    ~BaseDirectSolver()  {
    }

  };


} // namespace



#endif // OLAS_BASESOLVER_HH
