#ifndef FILE_TIMESTEPPING_2001
#define FILE_TIMESTEPPING_2001

#include "General/environment.hh"
#include "Utils/nodestoresol.hh"
#include "PDE/nodeEQN.hh"

#ifdef USE_OLAS
#include <olas.hh>
#else
#include <multigrid.hh>
#endif

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
  TimeStepping(std::string apdename, BaseSystem * algebraicsystem, NodeEQN * ptEQN);

   //! deconstructor
  virtual ~TimeStepping();
  
  //! initilization
  virtual void Init(Double * matrix_factors, Double dt)=0;

  //! perform predictor step
  virtual void Predictor(Vector<Double>& solold)=0;
  
  //! perform corrector step
  virtual void Corrector(Vector<Double>& solnew)=0;
  
  //! perform an update to RHS
  virtual void UpdateRHS()=0;

  //! perform an update to RHS with actual solution (for nonlin calculation)
  virtual void UpdateRHS(Vector<Double>& actSol)
  {Error("Error not implemented!",__FILE__,__LINE__);};

  //! set vector with first derivative
  virtual void SetDeriv1(const Vector<Double> & deriv1)
  {solderiv1_ = deriv1;}
  
  //! set vector with second derivative
  virtual void SetDeriv2(const Vector<Double> & deriv2)
  {solderiv2_ = deriv2;}

  //!  return pointer to vector with first derivative of solution
  virtual const Vector<Double>& GetDeriv1() const { return solderiv1_;}
  
  //! return pointer to vector with second derivative of solution
  virtual const Vector<Double>& GetDeriv2() const { return solderiv2_;}

  //! store solution to solution array (especially for effective mass formulation)
  
  virtual void StoreSolution(NodeStoreSol<Double> & solArr) const
  {Error("Not implemented in base class!", __FILE__, __LINE__);};

  NodeEQN * getNodeEQN(){return ptEQN_;};


protected:

  std::string pdename_;  //<! name of PDE
  BaseSystem * algsys_;  //<! pointer to algebraic system
  NodeEQN * ptEQN_;      //<! pointer to eqn-object

  Vector<Double> solderiv1_, solderiv2_;

private:
   


};

}

#endif // FILE_TIMESTEPPING
