#ifndef FILE_BDF2_2007
#define FILE_BDF2_2004

#include <General/environment.hh>
#include "timestepping.hh"


namespace CoupledField
{

  //! class for time stepping of parabolic PDE: method is Trapezoidal 

  class Bdf2: public TimeStepping
  {
  public:
    //! constructor
    //! \param algebraicsystem pointer to algebraic system 
    Bdf2(  BaseSystem * algebraicsystem );

    //! destructor
    virtual ~Bdf2();
  
    //! initilization
    //! \param rhsSIze total number of entries in the rhs vector
    void Init( Double dt, UInt rhsSize );

    //! perform predictor step
    void Predictor(Vector<Double>& solold);

    //! perform predictor step
    void RelaxedPredictor(Vector<Double>& solold, 
                          Vector<Double>& solpredRelaxed, 
                          Vector<Double>& solderiv1predRelaxed, 
                          Double omega){;};

    //! perform corrector step
    void Corrector(Vector<Double>& solnew);

    //! perform an update to RHS
    void UpdateRHS();

    //! Substract Stiffness from RHS
    virtual void SubstractStiffnessFromRHS(Vector<Double>& actSol)
    {;};

  private:
    
    //! compute parameters for multiplication
    void CalcParameters(Double dt);
   
    //! integration parameter
    Double gamma_;  

    bool firstTime_;

    //@{
    //! coefficients from Bdf2 method
    Double a0_,a1_,a2_,a3_,a4_,a5_,a6_,a7_;
    //@}

    //! predictor of solution
    Vector<Double> solpred_;
  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class Bdf2
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

#endif // FILE_BDF2
