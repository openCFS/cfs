#ifndef FILE_TIMESTEPPING_2001
#define FILE_TIMESTEPPING_2001

#include <General/environment.hh>
#include <Utils/array.hh>
#include <Utils/storesol.hh>

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
  TimeStepping(std::string apdename, BaseSystem * algebraicsystem);

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

  

  //!  return pointer to vector with first derivative of solution
  virtual const Vector<Double>& GetDeriv1() const { return solderiv1_;}
  
  //! return pointer to vector with second derivative of solution
  
  virtual const Vector<Double>& GetDeriv2() const { return solderiv2_;}

  //! store solution to solution array (especially for effective mass formulation)
  
  virtual const void StoreSolution(StoreSol<Double> & solArr) const 
  {Error("Not implemented in base class!", __FILE__, __LINE__);};


protected:

  std::string pdename_;  //<! name of PDE
  BaseSystem * algsys_;  //<! pointer to algebraic system

  Vector<Double> solderiv1_, solderiv2_;

private:
   


};

}

#endif // FILE_TIMESTEPPING
