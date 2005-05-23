#ifndef BUBBLE_ODE_HH
#define BUBBLE_ODE_HH

#include "General/environment.hh"
#include "Utils/StdVector.hh"
#include "BaseODEProblem.hh"

namespace CoupledField{

  //! Class derived from BaseODEProblem
  //! Base class from which all bubble problems are derived
  class BubbleODE : public BaseODEProblem {

  public:

    //! Default Constructor
    BubbleODE() {
      ENTER_FCN( "BubbleODE::BubbleODE" );
    }

    //! Default Destructor
    virtual ~BubbleODE() {
      ENTER_FCN( "BubbleODE::~BubbleODE" );
    }

    //! Compute the right hand side for dy/dt=f(t,y)
    //! \param t      Current time step (can also be used as current position
    //! \param vector y cointains starting values
    //! \param vector dydt contains on return the resulting rhs
    virtual void CompDeriv(const Double &t,
                           const StdVector<Double> &y,
                           StdVector<Double> &dydt) = 0;

    virtual Double GetP () = 0;

    virtual void SetP (Double p) = 0;

    virtual  Double GetDpdt () = 0;

    virtual void SetDpdt (Double dpdt) = 0;


  };

}


#endif //end of namespace
