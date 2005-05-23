#ifndef BASE_ODE_Problem_HH
#define BASE_ODE_Problem_HH

#include "General/environment.hh"
#include "Utils/StdVector.hh"

namespace CoupledField {

  //! Base class from which all ODE problems are derived
  class BaseODEProblem {

  public:

    //! Default Constructor
    BaseODEProblem() {
      ENTER_FCN( "BaseODEProblem::BaseODEProblem" );
    }

    //! Default Destructor
    virtual ~BaseODEProblem() {
      ENTER_FCN( "BaseODEProblem::~BaseODEProblem" );
    }

    //! Compute the right hand side for dy/dt=f(t,y)
    //! \param t      Current time step (can also be used as current position
    //! \param vector y cointains starting values
    //! \param vector dydt contains on return the resulting rhs
    virtual void CompDeriv(const Double &t,
                           const StdVector<Double> &y,
                           StdVector<Double> &dydt) = 0;


  };

}


#endif //end of namespace
