// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

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
    //! \param algebraicsystem pointer to algebraic system 
    Trapezoidal(  BaseSystem * algebraicsystem );

    //! destructor
    virtual ~Trapezoidal();
  
    //! initilization
    //! \param rhsSIze total number of entries in the rhs vector
    void Init( Double dt, UInt rhsSize );

    //! perform predictor step
    void Predictor(Vector<Double>& solold);

    //! perform corrector step
    void Corrector(Vector<Double>& solnew);

    //! perform an update to RHS
    void UpdateRHS();

    //! perform an update to RHS with actual solution (for nonlin calculation)
    void UpdateRHS(Vector<Double>& actSol);

    //! Substract Stiffness from RHS
    virtual void SubstractStiffnessFromRHS(Vector<Double>& actSol)
    {;};

  private:
    
    //! compute parameters for multiplication
    void CalcParameters(Double dt);
   
    //! integration parameter
    Double gamma_;  

    //@{
    //! coefficients from Trapezoidal method
    Double a0_,a1_,a2_,a3_,a4_,a5_,a6_,a7_;
    //@}

    //! predictor of solution
    Vector<Double> solpred_;
  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class Trapezoidal
  //! 
  //! \purpose 
  //! 
  //! \collab 
  //! 
  //! \implement 
  //! 
  //! \status In use.
  //! 
  //! \unused 
  //! 
  //! \improve
  //! 

#endif

}

#endif // FILE_NEWMARK
