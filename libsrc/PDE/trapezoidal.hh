#ifndef FILE_TRAPEZOIDAL_2004
#define FILE_TRAPEZOIDAL_2004

#include <General/environment.hh>
#include <DataInOut/conffile.hh>
#include "timestepping.hh"
#include <multigrid.hh>
#include <Utils/array.hh>
#include <General/environment.hh>

namespace CoupledField
{

//! class for time stepping of parabolic PDE: method is Trapezoidal 

class Trapezoidal: public TimeStepping
{
public:
  //! constructor
  /*!
    \param apdename name of PDE
    \param algebraicsystem pointer to algebraic system used by PDE
    \param dofspernode number of degree of freedom per node
    \param numnode number of nodes in PDE
    \param damp type of damping
  */
  Trapezoidal(std::string apdename, BaseSystem * algebraicsystem, Integer dofsprenode, 
	  Integer numnode);

   //! deconstructor
  virtual ~Trapezoidal();
  
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
  Double gamma_;  //<! integration parameter
  Double a0_,a1_,a2_,a3_,a4_,a5_,a6_,a7_; //<! coefficients from Trapezoidal method

  Array<Double> solpred_;
};

}

#endif // FILE_NEWMARK
