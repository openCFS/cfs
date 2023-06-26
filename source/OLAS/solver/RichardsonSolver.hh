#ifndef OLAS_RICHARDSON_HH
#define OLAS_RICHARDSON_HH

#include <iostream>
#include <fstream>



#include "BaseSolver.hh"

namespace CoupledField {

  class BasePrecond;
  

  //! Class for a Preconditioned Richardson scheme

  /*! This class represents a preconditioned Richardson iteration to solve
   * \f$P^{-1}Au=P^{-1}f \f$. In addition to the preconditioning step we use a 
   * relaxation parameter \f$\omega\f$.
   * In every iteration the residual is calculated as \f$r^k = f-Au^k\f$ and
   * the current approximation is updated as
   * \f$u^{k+1}=u^k + \omega P^{-1}r^k\f$.
   * The template parameter T relates to the scalar entry type of the system 
   * matrix in the linear system. The solver is applied to.
   * The RichardsonSolver object reads the following parameters from the 
   * xml node:
   * - stoppingRule = absNorm, relNormRHS or relNormRes0
   * - tol = specified tolerance for the stopping rule
   * - maxIter = maximal number of iterations
   * - logging = yes or no to write detailed output to info node
   * - consoleConvergence = yes or no to write the current tolerance for every iteration into terminal
   * - omega = relaxation parameter
   *
   * and writes the following keys into the report object myReport_:
   * - numIter = number of iterations performed
   * - finalPrecondResNorm = the Euclidean norm of the final preconditoned
   *   residual
   */
  template<typename T>
  class RichardsonSolver : public BaseIterativeSolver {

  public:
    //! Constructor

    //! This constructor does nothing but initialise the internal array
    //! pointers to NULL and set the pointers to the communication objects.
    //! \param myParams pointer to parameter object with steering parameters
    //!                 for this solver
    //! \param myReport pointer to report object for storing general
    //!                 information on solution process
    RichardsonSolver( PtrParamNode solverNode, PtrParamNode olasInfo )
      : r_(NULL),
        w_(NULL) {

      xml_ = solverNode;
      infoNode_ = olasInfo->Get("richardson");
    };

    //! Default destructor

    //! The default destructor needs to be deep to free dynamically allocated
    //! memory.
    ~RichardsonSolver();

    //! Solve the linear system sysmat*sol=rhs for sol

    //! \param sysmat System matrix of the linear system to be solved
    //! \param rhs Right-hand side vector of the linear system
    //! \param sol Solution vector of linear system
    void Solve(const BaseMatrix& sysmat,
	       const BaseVector& rhs, BaseVector& sol);

    //! Dummy method: Notify the solver that a new matrix pattern has been set
    void SetNewMatrixPattern() {EXCEPTION("SetNewMatrixPattern not implemented!");};
    
    //! Dummy setup method

    //! This method implements the pure virtual setup function defined in the
    //! BaseSolver class. In the case of the Richardson solver there is nothing
    //! to be done.
    void Setup( BaseMatrix &sysmat );

    //! Query type of the solver

    //! This method can be used to query the type of this solver. The answer
    //! is encoded as a value of the enumeration data type SolverType.
    //! \return RICHARDSON
    virtual SolverType GetSolverType() {
      return RICHARDSON;
    };

  private:

    //!\name auxiliary vectors
    //@{
    BaseVector *r_;  //!< current residual vector
    BaseVector *w_;  //!< the preconditioned residual vector
    //@}

    //! Default Constructor

    //! The default constructor is private in order to disallow its use. The
    //! reason is that we need pointers to the parameter and report objects to
    //! perform correct setup of an object of this class.
    RichardsonSolver(): r_(NULL), w_(NULL) {};

  };

} //namespace

#endif // OLAS_RichardsonSOLVER_HH
