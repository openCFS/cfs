#ifndef FILE_GILMOREDIMLOS
#define FILE_GILMOREDIMLOS

#include "General/environment.hh"
#include "Utils/StdVector.hh"
#include "BubbleODE.hh"
#include "Matrix/matrix.hh"

namespace CoupledField{

  //! Class to determine the bubble dynamics using the dimensionless Gilmore-equations. 



  class Gilmoredimlos : public BubbleODE {
 
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
    Gilmoredimlos(Double  RadiusInit, 
		  Double  density,
		  Double  sonicVel,
		  Double  pStatic,
		  Double  pVapour,
		  Double  surfacTen,
		  Double  polytrop,
		  Double  viscosity); 

     
    //! Default Destructor
    virtual ~Gilmoredimlos() {
      ENTER_FCN( "Gilmoredimlos::~Gilmoredimlos" );
    }

    //! Compute the right hand side for dy/dt=f(t,y)
    //! \param t      Current time step (can also be used as current position
    //! \param y      cointains starting values
    //! \param dydt   contains on return the resulting rhs
    void CompDeriv(const Double &t,
                   const StdVector<Double> &y,
                   StdVector<Double> &dydt);

    //! Get the pressure used in this class
    Double GetP (){
      ENTER_IFCN( "Gilmoredimlos::GetP" );      
      return p_;
    }

    //! Set the pressure 
    void SetP (Double p){
      ENTER_IFCN( "Gilmoredimlos::SetP" );
      p_ = p;
    }

    //! Get the derivative of the pressure used in this class
    Double GetDpdt (){ 
      ENTER_IFCN( "Gilmoredimlos::GetDpdt" );
      return dpdt_;
    }

    //! Set the derivative of the pressure 
    void SetDpdt (Double dpdt){
      ENTER_IFCN( "Gilmoredimlos::SetDpdt" );
      dpdt_ = dpdt;
    }
    
    //! Compute or approximate the Jacobian Matrix
    //! of the Gilmore-equation
    void  Jacobi(StdVector<Double> &y,
		 Matrix<Double> &dfdy,
		 Double &t);

    //! Help function to compute the Jacobian
    Double df2dy1(StdVector<Double> &y,
		  Double &t);

    //! Help function to compute the Jacobian
    Double df2dy2(StdVector<Double> &y,
		  Double &t);
   
   
  private:

    //  Constants for Gilmore-equations

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



    
    //! Constants needed for computation of enthalpy and new sonicvelocity
    Double A_;
    Double B_;
    Double n_;

    //! \param H_ Enthapydifference of the fluid close to the bubble
    //! and far away
    Double H_;

    //! \param C_ Sonicvelocity of the bubbly-liquid 
    Double C_;

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Dimensionless constants
    Double sigma_;
    Double c0_;
    Double mu_;
    Double a_;
    Double b_;
    Double pV_;

  };



#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class Gilmoredimlos
  //! 
  //! \purpose
  //! The equation is known as the  Gilmore-Akulichev model, it
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
  //! So far the choice of computing or approximating
  //! the Jacobian Matrix has to done by 
  //! directly commenting in and out

#endif

} // end of namespace

#endif

