#ifndef FILE_NEWMARK_2001
#define FILE_NEWMARK_2001

#include <General/environment.hh>
#include <DataInOut/conffile.hh>
#include "timestepping.hh"
#include <multigrid.hh>
#include <Utils/array.hh>
#include <General/environment.hh>

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
    \param dofspernode number of degree of freedom per node
    \param numnode number of nodes in PDE
    \param damp type of damping
  */
  Newmark(std::string apdename, BaseSystem * algebraicsystem, Integer dofsprenode, 
	  Integer numnode, DampingType damp);

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

  //! perform an update to RHS with actual solution (for nonlin calculation)
  virtual void UpdateRHS(std::vector<Double>& actSol);  

  //! compute parameters for multiplication
  void CalcParameters(Double dt);

  //! set the time step
  void SetTimeStep(Double dt) 
  { dt_ = dt;};

private:
   
  Double dt_; //<! time step
  Double alpha_, gamma_, beta_;  //<! integration parameters
  Double a0_,a1_,a2_,a3_,a4_,a5_,a6_,a7_; //<! coefficients from Newmark method

  DampingType damping_;

  Array<Double> solpred_, solderiv1pred_;
};

}

#endif // FILE_NEWMARK
