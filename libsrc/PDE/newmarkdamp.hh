#ifndef FILE_NEWMARKDAMP_2001
#define FILE_NEWMARKDAMP_2001

#include <General/environment.hh>
#include <DataInOut/ParamHandling/ConfFile.hh>
#include <Domain/grid.hh>

#include "timestepping.hh"

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
  NewmarkDamp (std::string apdename, BaseSystem * algebraicsystem, NodeEQN * ptEQN, 
	       Grid * aptgrid, StdVector<std::string> adampingList, 
	       StdVector<std::string> subdomainList, Integer adamping, 
	       Integer afrac_memory, Double damp, Boolean isaxi);

  //! deconstructor
  virtual ~NewmarkDamp();
  
  //! initilization
  virtual void Init(Double * matrix_factors, Double dt);

  //! perform predictor step
  virtual void Predictor(Vector<Double>& solold);

  //! perform corrector step
  virtual void Corrector(Vector<Double>& solnew);

  //! perform an update to RHS
  virtual void UpdateRHS();

  //! compute parameters for multiplication
  void CalcParameters(Double dt);

  //! set the time step
  void SetTimeStep(Double dt) 
  { dt_ = dt;};

private:

  //
  void GetElemSolution (const Vector<Double>& sol, Vector<Double>& elemsol, const StdVector<Integer> & connectPDE);
  
  Double dt_; //<! time step
  Double alpha_, gamma_, beta_;  //<! integration parameters
  Double a0_,a1_,a2_,a3_,a4_,a5_,a6_,a7_; //<! coefficients from NewmarkDamp method

  Double y_;  //!< y-parameter
  Integer frac_memory_;

  Integer damping_;

  Vector<Double> solpred_, solderiv1pred_;
  Vector<Double> *solfrac_;

  //
  Grid * ptgrid_;

  //
  StdVector<std::string> dampingList_;

  //
  StdVector<std::string> subdoms_;

  //
  Boolean isaxi_;
  
  };

}

#endif // FILE_NEWMARKDAMP
