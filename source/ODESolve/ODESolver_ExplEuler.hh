// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef ODE_SOLVER_EXPL_EULER_HH
#define ODE_SOLVER_EXPL_EULER_HH

#include "BaseODESolver.hh"
#include "General/defs.hh"

namespace CoupledField {

  //! Explizit Euler method to solve ODEs
class BaseODEProblem;
template <class TYPE> class StdVector;

  class ODESolver_ExplEuler : public BaseODESolver {

  public:

    //! Default Constructor
    ODESolver_ExplEuler() {
    }

    //! Default Destructor
    virtual ~ODESolver_ExplEuler() {
    }

    //! Compute the solution of the initial value problem
    //! \param tInit     initial time
    //! \param tStop     final time
    //! \param y         containing on input the initial values and on output
    //!                  the solution
    //! \param myODE     object containing information on the right hand side
    //!                  function of the ODE
    //! \param hInit     Suggestion for size of first time step
    //! \param hMin      Minimal allowed size for time step
    //! \param hMax      Maximal allowed size for time step
    void Solve( const Double tInit,
                const Double tStop,
                StdVector<Double> &y,
                BaseODEProblem &myODE,
                Double hInit = -1.0,
                Double hMin = -1.0,
                Double hMax = -1.0);

  };

}


#endif
