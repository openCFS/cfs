#ifndef FILE_NEWMARKDAMP_2001
#define FILE_NEWMARKDAMP_2001

#include <General/environment.hh>
#include <DataInOut/conffile.hh>
#include "timestepping.hh"
#include <Utils/array.hh>

#ifdef USE_OLAS
#include <olas.hh>
#else
#include <multigrid.hh>
#endif

namespace CoupledField
{

//! class for time stepping of hyperbolic PDE: method is NewmarkDamp

class NewmarkDamp: public TimeStepping
{
public:
  //! constructor
  /*!
    \param apdename name of PDE
    \param algebraicsystem pointer to algebraic system used by PDE
  */
  NewmarkDamp(std::string apdename, BaseSystem * algebraicsystem, Integer dofsprenode, 
	  Integer numnode, Integer adamping, Integer afrac_memory, Double damp);

   //! deconstructor
  virtual ~NewmarkDamp();
  
  //! initilization
  virtual void Init(Double * matrix_factors, Double dt);

  //! perform predictor step
  virtual void Predictor(Array<Double>& solold);

  //! perform corrector step
  virtual void Corrector(Array<Double>& solnew);

  //! perform an update to RHS
  virtual void UpdateRHS();

  //! compute parameters for multiplication
  void CalcParameters(Double dt);

  //! set the time step
  void SetTimeStep(Double dt) 
  { dt_ = dt;};

private:
   
  Double dt_; //<! time step
  Double alpha_, gamma_, beta_;  //<! integration parameters
  Double a0_,a1_,a2_,a3_,a4_,a5_,a6_,a7_; //<! coefficients from NewmarkDamp method

  Double y_;  //!< y-parameter
  Integer frac_memory_;

  Integer damping_;

  Array<Double> solpred_, solderiv1pred_;
  Array<Double> *solfrac_;
};

}

#endif // FILE_NEWMARKDAMP
