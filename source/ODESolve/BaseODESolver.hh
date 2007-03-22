// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef BASE_ODE_SOLVER_HH
#define BASE_ODE_SOLVER_HH

#include "General/environment.hh"
#include "Utils/StdVector.hh"
#include "ODEDescr/BaseODEProblem.hh"

namespace CoupledField {

  //! Base class from which all ODE solvers are derived
  class BaseODESolver {

  public:

    //! Default Constructor
    BaseODESolver() {
      ENTER_FCN( "BaseODESolver::BaseODESolver" );
      successLastSolve_     = false;
      numStepsLastSolve_    = 0;
      numBadStepsLastSolve_ = 0;
      eps_                  = 1e-6;
      maxSteps_             = 40000000;
      safetyFac_            = 0.5;
    }

    //! Default Destructor
    virtual ~BaseODESolver() {
      ENTER_FCN( "BaseODESolver::~BaseODESolver" );
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
    virtual void Solve( const Double tInit,
                        const Double tStop,
                        StdVector<Double> &y,
                        BaseODEProblem &myODE,
                        Double &hInit,
                        Double hMin = -1.0,
                        Double hMax = -1.0) = 0;

    //! Query status information on last solve
    //! \param success     Was the last solve successful?
    //! \param numSteps    Number of time steps for last solve
    //! \param numBadSteps Number of rejected time steps in last solve
    void GetStatus( bool &success, Integer &numSteps, Integer &numBadSteps ) {
      ENTER_FCN( "BaseODESolver::GetStatus" );
      success     = successLastSolve_;
      numSteps    = numStepsLastSolve_;
      numBadSteps = numBadStepsLastSolve_;
    }

    Double GetEps (){ return eps_;}

    void SetEps (Double epsNew){eps_ = epsNew;}

    Integer GetMaxSteps (){ return maxSteps_;}

    void SetMaxSteps(Integer maxStepsNew) {
      ENTER_IFCN( "BaseODESolver::SetMaxSteps" );
      maxSteps_ = maxStepsNew;
    }

    Double GetSafetyFac (){ return safetyFac_;}

    void SetSafetyFac (Double safetyFacNew){safetyFac_ = safetyFacNew;}

    void SetNumEl (Integer numEl){
      ENTER_IFCN( "ODESolver_RKF45::SetNumEl" );
      numEl_ = numEl;
    }


  protected:

    Integer numEl_ ;

    //! Was the last solve attempt successful?
    bool successLastSolve_;

    //! Number of time steps in last solve
    Integer numStepsLastSolve_;

    //! Number of rejected time steps in last solve
    Integer numBadStepsLastSolve_;

    //! Threshold for stopping test
    Double eps_;

    //! \param maxSteps  Maximal allowed number of time steps
    Integer maxSteps_;

    //! \param safetyFac Safety factor for step size computation
    Double safetyFac_;


    bool RadiusGroesserNull_;

  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class BaseODESolver
  //! 
  //! \purpose
  //! 
  //! \collab
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
