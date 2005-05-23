#ifndef FILE_GILMORE
#define FILE_GILMORE

#include "General/environment.hh"
#include "Utils/StdVector.hh"
#include "BubbleODE.hh"

namespace CoupledField{

  //! Class for bubble dynamics using the equation of Gilmore together 
  //! with the Tait-equation; also known as Gilmore-Akulichev model.
  //! Class is derived from BubbleODE 

  class Gilmore : public BubbleODE {
 
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
    Gilmore(Double  RadiusInit, 
            Double  density,
            Double  sonicVel,
            Double  pStatic,
            Double  pVapour,
            Double  surfacTen,
            Double  polytrop,
            Double  viscosity); 
     
    //! Default Destructor
    virtual ~Gilmore() {
      ENTER_FCN( "Gilmore::~Gilmore" );
    }

    //! Compute the right hand side for dy/dt=f(t,y)
    //! \param t      Current time step (can also be used as current position
    //! \param y      cointains starting values
    //! \param dydt   contains on return the resulting rhs
    void CompDeriv(const Double &t,
                   const StdVector<Double> &y,
                   StdVector<Double> &dydt);

    Double GetP (){
      ENTER_IFCN( "Gilmore::GetP" );      
      return p_;
    }

    void SetP (Double p){
      ENTER_IFCN( "Gilmore::SetP" );
      p_ = p;
    }

    Double GetDpdt (){ 
      ENTER_IFCN( "Gilmore::GetDpdt" );
      return dpdt_;
    }

    void SetDpdt (Double dpdt){
      ENTER_IFCN( "Gilmore::SetDpdt" );
      dpdt_ = dpdt;
    }


   
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
    //! Later p and dpdt will become vector to be able to contain the data for
    //! all elements in one time step


    
    //! Constants needed for computation of enthalpy and new sonicvelocity
    //! Later these constants should be read from a file, either xml or mat 
    Double A_;
    Double B_;
    Double n_;

    //! \param H_     Enthapy of the fluid
    Double H_;

    //! \param sonicVelMix_ Sonicvelocity of the bubbly-liquid 
    Double sonicVelMix_;

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Dimensionless constants
    Double D_;
    Double v_;
    Double mu_;
    Double a_;
    Double b_;
    
  };

} // end of namespace

#endif

