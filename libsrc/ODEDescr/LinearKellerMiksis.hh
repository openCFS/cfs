#ifndef FILE_LINEARKELLERMIKSIS
#define FILE_LINEARKELLERMIKSIS

#include "General/environment.hh"
#include "Utils/StdVector.hh"
#include "BubbleODE.hh"
#include "Matrix/matrix.hh"

namespace CoupledField{

  //! Class to determine the bubble dynamics using the LinearKellerMiksis-equations. 


  class LinearKellerMiksis : public BubbleODE {
 
  public:

    //! Constructor
    //! \param RadiusInit Equilibrium radius of bubble
    //! \param density    Density of the undisturbed fluid
    //! \param sonocVel   Sonic velocity of the undisturbed fluid
    //! \param pStatic    Static pressure
    //! \param pVapour    Vapour pressure of fluid
    //! \param surfacTen  Surfacetension of fluid
    //! \param polytrop   Polytropic index
    //! \param viscosity  Viscosity of fluid
    LinearKellerMiksis(Double  RadiusInit, 
            Double  density,
            Double  sonicVel,
            Double  pStatic,
            Double  pVapour,
            Double  surfacTen,
            Double  polytrop,
            Double  viscosity); 
     
    //! Default Destructor
    virtual ~LinearKellerMiksis() {
      ENTER_FCN( "LinearKellerMiksis::~LinearKellerMiksis" );
    }

    //! Compute the right hand side for dy/dt=f(t,y)
    //! \param t      Current time step (can also be used as current position
    //! \param y      cointains starting values
    //! \param dydt   contains on return the resulting rhs
    void CompDeriv(const Double &t,
                   const StdVector<Double> &y,
                   StdVector<Double> &dydt);

    void  Jacobi(StdVector<Double> &y,
		 Matrix<Double> &dfdy,
		 Double &t);
      

    //! Get the pressure used in this class
    Double GetP (){
      ENTER_IFCN( "LinearKellerMiksis::GetP" );      
      return p_;
    }

    //! Set the pressure 
    void SetP (Double p){
      ENTER_IFCN( "LinearKellerMiksis::SetP" );
      p_ = p;
    }

    //! Get the derivative of the pressure used in this class
    Double GetDpdt (){ 
      ENTER_IFCN( "LinearKellerMiksis::GetDpdt" );
      return dpdt_;
    }

    //! Set the derivative of the pressure 
    void SetDpdt (Double dpdt){
      ENTER_IFCN( "LinearKellerMiksis::SetDpdt" );
      dpdt_ = dpdt;
    }
    

   
  private:

    //  Constants for LinearKellerMiksis-equations

    //! \param RadiusInit_ Equilibrium radius of bubble
    Double   RadiusInit_; 
    //! \param density_    Density of fluid
    Double   density_;
    //! \param sonocVel_   Sonic velocity of fluid
    Double   sonicVel_;
    //! \param pStatic_    Static pressure
    Double   pStatic_;
    //! \param pVapour_    Vapour pressure of fluid
    Double   pVapour_;
    //! \param surfacTen_  Surfacetension of fluid
    Double   surfacTen_;
    //! \param polytrop_   Polytropic index
    Double   polytrop_;
    //! \param viscosity_  Viscosity of fluid    
    Double   viscosity_;

    //! \param p_          Pressure computed by acoustics
    Double   p_;
    //! \param dpdt_       Time derivative of pressure computed by acoustics
    Double   dpdt_;
    //! Later p and dpdt will become vector to be able to contain the data for
    //! all elements in one time step



    
  };


#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class LinearKellerMiksis
  //! 
  //! \purpose
  //! The equation is known as the  linearised KellerMiksis model, it
  //! uses the Tait-equation.
  //! 
  //! \collab
  //! Class is derived from BubbleODE. 
  //! Class ist called from ODESolver_RKF45,
  //! ODESolver_Rosenbrock and solveStepAcousticBubble 
  //! 
  //! \implement 
  //! 
  //! \status In use
  //! 
  //! \unused 
  //! 
  //! \improve


#endif


} // end of namespace

#endif

