#ifndef FILE_PIEZOSWITCHODE_HH
#define FILE_PIEZOSWITCHODE_HH

#include "BaseODEProblem.hh"

#include "General/Environment.hh"
#include "Utils/StdVector.hh"
#include "MatVec/Matrix.hh"


namespace CoupledField {

  //! Base class for piezoelectric switching models
  class PiezoSwitchODE : public BaseODEProblem {

  public:

    //! Default Constructor
    PiezoSwitchODE(UInt numberOfODEs, Double initT, Double alpha);

    //! Default Destructor
    virtual ~PiezoSwitchODE() {
    }

    void setDynamicCoefficients(Matrix<Double>& coeff ) {
      coeffs_ = coeff;
    };

    //! Compute the right hand side for dy/dt=f(t,y)
    //! \param t      Current time step (can also be used as current position)
    //! \param vector y cointains starting values
    //! \param vector dydt contains on return the resulting rhs
    void CompDeriv(const Double &t,
                           const StdVector<Double> &y,
                           StdVector<Double> &dydt);

    //! Compute the Jacobian Matrix of the ode
    void Jacobi(StdVector<Double> &y,
			Matrix<Double> &dfdy,
			Double &t);
    //!
    UInt numODEs_;
    //!
    Double dtPDE_;
    //!
    Double alpha_;
    //!
    Matrix<Double> coeffs_;
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
