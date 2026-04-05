#ifndef ODE_SOLVER_RKF45_HH
#define ODE_SOLVER_RKF45_HH

#include "General/Environment.hh"
#include "Utils/StdVector.hh"
#include "Utils/ToolsFull.hh"
#include "BaseODESolver.hh"  

namespace CoupledField {


  //! Runge-Kutta solver for systems of first order ode's.
  class ODESolver_RKF45 : public BaseODESolver {

  public:

    //! Default Constructor
    ODESolver_RKF45()  {
    }

    //! Default Destructor
    ~ODESolver_RKF45() {
    }

    //! Compute the solution of the initial value problem
    //! \param tInit     initial time
    //! \param tStop     final time
    //! \param yIintOut  containing on input the initial values and on output
    //!                  the solution
    //! \param myODE     object containing information on the right hand side
    //!                  function of the ODE
    //! \param hInit     Suggestion for size of first time step
    //! \param hMin      Minimal allowed size for time step
    //! \param hMax      Maximal allowed size for time step
    void Solve( const Double tInit,
		const Double tStop,
		StdVector<Double> &yInitOut,
		BaseODEProblem &myODE,
		Double &hInit,
		Double hMin = 0.0,
		Double hMax = -1.0 );




    //   void SetNumEl (Integer numEl){
    // numEl_ = numEl;
    //}

  private:

    // Integer numEl_ ;
    //! Controlls RKF45 step by monitoring local truncation error
    //! and adjusting stepsize
    //! \param y         containing on input the initial values and on output
    //!                  the solution
    //! \param dydt      containing the known derivatives of y at time t
    //! \param t         current time
    //! \param htry      Suggested step size
    //! \param vector    yscal: vector against which the error is scaled 
    //! \param hdid      Stepsize that was actually accomplished
    //! \param hnext     Estimation for next time step
    //! \param myODE     object containing information on the right hand side
    //!                  function of the ODE 
    void RKAdaptiveStepsize (StdVector<Double> &y,
                             StdVector<Double> &dydt,
                             Double &t,
                             const Double hTry,
                             const StdVector<Double> &yScal,
                             Double &hDid,
                             Double &hNext,
                             BaseODEProblem &myODE);


    //! Compute one Cash-Karp Runge-Kutta step
    //! \param y         containing the given value at time t
    //! \param dydt      containing the known derivatives of y at time t
    //! \param t         current time
    //! \param h         step size
    //! \param vector    yout contains the incremented output values
    //! \param vector    yerror Estimation of local truncation error in yerror
    //! \param myODE     object containing information on the right hand side
    //!                  function of the ODE 
    void RKCashKarp (const StdVector<Double> &y,
                     const StdVector<Double> &dydt,
                     const Double t,
                     const Double h,
                     StdVector<Double> &yOut,
                     StdVector<Double> &yError,
                     BaseODEProblem &myODE);
                     
  };


#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class ODESolver_RKF45
  //! 
  //! \purpose
  //! Integrator for solution of an initial value problem using a
  //! Runge-Kutta-Fehlberg method of 4th order with a step size control of
  //! 5th order using the Runge-Kutta parameters by Cash-Karp. See Numerical
  //! Recipes.
  //! 
  //! \collab
  //! Needs a user provided routine to compute the right hand side
  //! of the ode dy/dt = f(t,y); here this is done via the pointer
  //! to the bubble-class ( Gilmoredimlos, Gilmore, Keller-Miksis).
  //! 
  //! \implement 
  //! 
  //! \status In use
  //! 
  //! \unused 
  //! 
  //! \improve


#endif

}


#endif
