#ifndef ODE_SOLVER_EXPL_EULER_HH
#define ODE_SOLVER_EXPL_EULER_HH

#include "General/environment.hh"
#include "Utils/StdVector.hh"
#include "ODEDescr/BaseODEProblem.hh"
#include "BaseODESolver.hh"

namespace CoupledField {

  //! Explizit Euler method to solve ODEs
  class ODESolver_ExplEuler : public BaseODESolver {

  public:

    //! Default Constructor
    ODESolver_ExplEuler() {
      ENTER_FCN( "ODESolver_ExplEuler::ODESolver_ExplEuler" );
    }

    //! Default Destructor
    virtual ~ODESolver_ExplEuler() {
      ENTER_FCN( "ODESolver_ExplEuler::~ODESolver_ExplEuler" );
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
