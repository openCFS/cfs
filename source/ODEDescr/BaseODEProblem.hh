// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef BASE_ODE_Problem_HH
#define BASE_ODE_Problem_HH

#include "General/environment.hh"
#include "Utils/StdVector.hh"
#include "Matrix/matrix.hh"

namespace CoupledField {

  //! Base class from which all ODE problems are derived
  class BaseODEProblem {

  public:

    //! Default Constructor
    BaseODEProblem() {
    }

    //! Default Destructor
    virtual ~BaseODEProblem() {
    }

    //! Compute the right hand side for dy/dt=f(t,y)
    //! \param t      Current time step (can also be used as current position)
    //! \param vector y cointains starting values
    //! \param vector dydt contains on return the resulting rhs
    virtual void CompDeriv(const Double &t,
                           const StdVector<Double> &y,
                           StdVector<Double> &dydt) = 0;

    //! Compute the Jacobian Matrix of the ode
    virtual void Jacobi(StdVector<Double> &y,
			Matrix<Double> &dfdy,
			Double &t) //= 0;
    {
      std::cerr<< "Jacobi method is not yet implemented for this problem" << std::endl;
    }


  };



#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class BaseODEProblem
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
  //! 

#endif

}


#endif //end of namespace
