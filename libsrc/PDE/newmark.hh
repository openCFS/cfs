#ifndef FILE_NEWMARK_2001
#define FILE_NEWMARK_2001

#include <General/environment.hh>
#include <DataInOut/conffile.hh>
#include "timestepping.hh"
#include <multigrid.hh>
#include <Utils/array.hh>

namespace CoupledField
{

//! class for time stepping of hyperbolic PDE: method is Newmark

class Newmark: public TimeStepping
{
public:
  //! constructor
  /*!
    \param apdename name of PDE
    \param algebraicsystem pointer to algebraic system used by PDE
  */
  Newmark(std::string apdename, BaseSystem * algebraicsystem, Integer dofsprenode, 
	  Integer numnode, Integer adamping);

   //! deconstructor
  virtual ~Newmark();
  
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
  Double a0_,a1_,a2_,a3_,a4_,a5_,a6_,a7_; //<! coefficients from Newmark method

  Integer damping_;

  Array<Double> solpred_, solderiv1pred_;
};

}

#endif // FILE_NEWMARK
