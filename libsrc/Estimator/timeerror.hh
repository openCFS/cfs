#ifndef FILE_TIMEERROR_2001
#define FILE_TIMEERROR_2001

#include "Utils/tools.hh"
#include "PDE/basepde.hh"

namespace CoupledField
{

//! Class where we implement time error estimator
class TimeErrorEstimator
{
public:

  //! constructor with pointer to PDE
  TimeErrorEstimator(BasePDE* aptpde);

  //! Deconstructor
  virtual ~TimeErrorEstimator(){;}

  //!  return true, if error is more than tolerance
  /*!
    \param dt time step
  */
  virtual Boolean TestError(const Double dt)=0;

  //! change time step
  virtual void ChangeStep(Double &)=0;

  //! calculation of third derivative
  virtual void CalcThirdDer()=0;

  //! calculation of error
  virtual void CalcError(const Double dt)=0;  

protected:
  
   //! pointer to PDE
   BasePDE* ptPDE_;

   //! relative error
   Double relativeerror; 
};

} // end of namespace

#endif
