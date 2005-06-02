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
    //! \param rhsSIze total number of entries in the rhs vector
    Trapezoidal(  BaseSystem * algebraicsystem, UInt rhsSize );

    //! destructor
    virtual ~Trapezoidal();
  
    //! initilization
    void Init( std::map<FEMatrixType,Double> & matrix_factors, 
               Double dt );

    //! perform predictor step
    void Predictor(Vector<Double>& solold);

    //! perform corrector step
    void Corrector(Vector<Double>& solnew);

    //! perform an update to RHS
    void UpdateRHS();

    //! perform an update to RHS with actual solution (for nonlin calculation)
    void UpdateRHS(Vector<Double>& actSol);

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
