#ifndef OLAS_DIAGSOLVER_HH
#define OLAS_DIAGSOLVER_HH

#include <iostream>
#include <fstream>



#include "BaseSolver.hh"

namespace CoupledField {

  //! Class for a solver, where the system matrix is a diagonal matrix

  /*! This class represents a a simple solver, where the system
   * matrix is a diagonal matrix. Therewith, we just have to divide
   * the RHS-vector by the diagonal entries
   */
  template<typename T>
  class DiagSolver : public BaseIterativeSolver {

  public:
    //! Constructor

    //! This constructor does nothing but initialise the internal array
    //! pointers to NULL and set the pointers to the communication objects.
    //! \param myParams pointer to parameter object with steering parameters
    //!                 for this solver
    //! \param myReport pointer to report object for storing general
    //!                 information on solution process
    DiagSolver( PtrParamNode xml, PtrParamNode olasInfo ) {
      xml_ = xml;
      infoNode_ = olasInfo->Get("diagsolver");
    };

    //! Default destructor

    //! The default destructor needs to be deep to free dynamically allocated
    //! memory.
    ~DiagSolver();

    //! Solve the linear system sysmat*sol=rhs for sol

    //! \param sysmat System matrix of the linear system to be solved
    //! \param rhs Right-hand side vector of the linear system
    //! \param sol Solution vector of linear system
    void Solve(const BaseMatrix& sysmat,
	       const BaseVector& rhs, BaseVector& sol );

    //! Dummy method: Notify the solver that a new matrix pattern has been set
    void SetNewMatrixPattern() {EXCEPTION("SetNewMatrixPattern not implemented!");};

    //! Dummy setup method

    //! This method implements the pure virtual setup function defined in the
    //! BaseSolver class. In our case there is nothing to be done.
    void Setup( BaseMatrix &sysmat ) {};

    //! Query type of the solver

    //! This method can be used to query the type of this solver. The answer
    //! is encoded as a value of the enumeration data type SolverType.
    //! \return RICHARDSON
    virtual SolverType GetSolverType() {
      return DIAGSOLVER;
    };

  private:

    //! Default Constructor

    //! The default constructor is private in order to disallow its use. The
    //! reason is that we need pointers to the parameter and report objects to
    //! perform correct setup of an object of this class.
    DiagSolver(){};

  };

} //namespace

#endif // OLAS_RichardsonSOLVER_HH
