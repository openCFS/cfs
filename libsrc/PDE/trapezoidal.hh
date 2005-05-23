#ifndef FILE_TRAPEZOIDAL_2004
#define FILE_TRAPEZOIDAL_2004

#include <General/environment.hh>
#include "timestepping.hh"


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
    Trapezoidal(std::string apdename, BaseSystem * algebraicsystem, NodeEQN * ptEQN);

    //! deconstructor
    virtual ~Trapezoidal();
  
    //! initilization
    virtual void Init(Double * matrix_factors, Double dt);

    //! perform predictor step
    virtual void Predictor(Vector<Double>& solold);

    //! perform corrector step
    virtual void Corrector(Vector<Double>& solnew);

    //! perform an update to RHS
    virtual void UpdateRHS();

    //! perform an update to RHS with actual solution (for nonlin calculation)
    virtual void UpdateRHS(Vector<Double>& actSol);

    //! compute parameters for multiplication
    void CalcParameters(Double dt);

  private:
   
    Double gamma_;  //<! integration parameter
    Double a0_,a1_,a2_,a3_,a4_,a5_,a6_,a7_; //<! coefficients from Trapezoidal method

    Vector<Double> solpred_;
  };

}

#endif // FILE_NEWMARK
