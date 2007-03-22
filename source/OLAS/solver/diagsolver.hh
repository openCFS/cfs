// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef OLAS_DIAGSOLVER_HH
#define OLAS_DIAGSOLVER_HH

#include <iostream>
#include <fstream>

#include "utils/utils.hh"
#include "algsys/olasparams.hh"
#include "solver/basesolver.hh"
#include "matvec/matvec.hh"
#include "precond/precond.hh"

namespace OLAS {

  //! Class for a solver, where the system matrix is a diagonal matrix

  /*! This class represents a a simple solver, where the system
   * matrix is a diagonal matrix. Therewith, we just have to divide
   * the RHS-vector by the diagonal entries
   */
  template<typename T>
  class DiagSolver : public BaseIterativeSolver {

  public:

    //!typename of matrix entries (=T)
    typedef typename assocType<T>::T_Mtype T_Mtype;	

    //!tiny vector of the same dimension as matrix block
    typedef typename assocType<T>::T_Vtype T_Vtype;	

    //scalar of the same primitive data type as matrix
    typedef typename assocType<T>::T_Stype T_Stype;	

    //! Constructor

    //! This constructor does nothing but initialise the internal array
    //! pointers to NULL and set the pointers to the communication objects.
    //! \param myParams pointer to parameter object with steering parameters
    //!                 for this solver
    //! \param myReport pointer to report object for storing general
    //!                 information on solution process
    DiagSolver( OLAS_Params *myParams, OLAS_Report *myReport ) {
      ENTER_FCN( "DiagSolver::DiagSolver" );
      myParams_ = myParams;
      myReport_ = myReport;
    };

    //! Default destructor

    //! The default destructor needs to be deep to free dynamically allocated
    //! memory.
    ~DiagSolver();

    //! Solve the linear system sysmat*sol=rhs for sol

    //! \param sysmat System matrix of the linear system to be solved
    //! \param premat The preconditioning object
    //! \param rhs Right-hand side vector of the linear system
    //! \param sol Solution vector of linear system
    void Solve(const BaseMatrix& sysmat, const BasePrecond& premat, 
	       const BaseVector& rhs, BaseVector& sol );

    //! Dummy setup method

    //! This method implements the pure virtual setup function defined in the
    //! BaseSolver class. In our case there is nothing to be done.
    void Setup( BaseMatrix &sysmat ) {};

    //! Query type of the solver

    //! This method can be used to query the type of this solver. The answer
    //! is encoded as a value of the enumeration data type SolverType.
    //! \return RICHARDSON
    virtual SolverType GetSolverType() {
      ENTER_IFCN( "DiagSolver::GetSolverType" );
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
