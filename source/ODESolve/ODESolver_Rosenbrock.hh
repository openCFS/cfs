#ifndef ODE_SOLVER_ROSENBROCK_HH
#define ODE_SOLVER_ROSENBROCK_HH

#include "General/Environment.hh"
#include "Utils/StdVector.hh"
#include "Utils/ToolsFull.hh"
#include "BaseODESolver.hh"  
#include "MatVec/Matrix.hh"

namespace CoupledField {

  //! Implicit Rosenbrock solver for systems of first order ode's.
  class ODESolver_Rosenbrock : public BaseODESolver {

  public:

    //! Default Constructor
    ODESolver_Rosenbrock()  {


    }

    //! Default Destructor
    ~ODESolver_Rosenbrock() {
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
		Double hMin = -1.0,
		Double hMax = -1.0 );





  private:


    //! Computes and controlls Rosenbrock step by monitoring local 
    //! truncation error and adjusting stepsize
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
    void RosenbAdaptiveStepsize (StdVector<Double> &y,
				 StdVector<Double> &dydt,
				 Double &t,
				 const Double hTry,
				 const StdVector<Double> &yScal,
				 Double &hDid,
				 Double &hNext,
				 BaseODEProblem &myODE);


    // LU-Decomposition
    void LUDecomposition(Matrix<Double> &a,
			 UInt n,
			 StdVector<UInt> &indx,
			 Double &d);
                     


    // LU-BAcksubstitution
    void LUBackSub(Matrix<Double> &a,
		   UInt n,
		   StdVector<UInt> &indx,
		   StdVector<Double> &b);
    
  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class ODESolver_Rosenbrock
  //! 
  //! \purpose
  //! Integrator for solution of an initial value problem using a
  //! semi-implicit Rosenbrock method of 4th order with a step size control of
  //! 5th order. See Numerical Recipes.
  //! 
  //! \collab
  //! Needs a user provided routine to compute the right hand side
  //! of the ode dy/dt = f(t,y) and computation of the Jacobian Matrix; 
  //! here this is done via the pointer
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
