#ifndef FILE_TIMESTEPPING_2001
#define FILE_TIMESTEPPING_2001

#include <General/environment.hh>
#include <multigrid.hh>
#include <Utils/array.hh>

namespace CoupledField
{

//! a base class for time stepping

class TimeStepping
{
public:
  //! constructor
  /*!
    \param apdename name of PDE
    \param algebraicsystem pointer to algebraic system used by PDE
  */
  TimeStepping(std::string apdename, BaseSystem * algebraicsystem);

   //! deconstructor
  virtual ~TimeStepping();
  
  //! initilization
  virtual void Init(Double * matrix_factors, Double dt)=0;

  //! perform predictor step
  virtual void Predictor(Array<Double>& solold)=0;

  //! perform corrector step
  virtual void Corrector(Array<Double>& solnew)=0;

  //! perform an update to RHS
  virtual void UpdateRHS()=0;

 //!  return pointer to vector with first derivative of solution
  virtual const Array<Double>& GetDeriv1() const { return solderiv1_;}

  //! return pointer to vector with second derivative of solution
  virtual const Array<Double>& GetDeriv2() const { return solderiv2_;}

protected:

  std::string pdename_;  //<! name of PDE
  BaseSystem * algsys_;  //<! pointer to algebraic system

  Array<Double> solderiv1_, solderiv2_;

private:
   


};

}

#endif // FILE_TIMESTEPPING
