#include <fstream>
#include <iostream>
#include <string>

#include <cmath>     


#include "LinearKellerMiksis.hh"


namespace CoupledField
{

  LinearKellerMiksis::LinearKellerMiksis(Double  RadiusInit, 
                             Double  density,
                             Double  sonicVel,
                             Double  pStatic,
                             Double  pVapour,
                             Double  surfacTen,
                             Double  polytrop,
                             Double  viscosity) {

    ENTER_FCN( "LinearKellerMiksis::LinearKellerMiksis" );

    RadiusInit_ = RadiusInit;
    density_    = density;
    sonicVel_   = sonicVel;
    pStatic_    = pStatic;
    pVapour_    = pVapour;
    surfacTen_  = surfacTen;
    polytrop_   = polytrop;
    viscosity_  = viscosity;

//     //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//     //Compute dimensionless constants needed in CompDeriv

//     sigma_ = 2 * surfacTen_ / (RadiusInit_ * pStatic_);
//     c0_    = sonicVel_ / (sqrt( pStatic_ / density_ ));
//     mu_    = viscosity_ / ( pStatic_ * RadiusInit_) * sqrt( pStatic_ / density_);
//     pV_    = pVapour_ / pStatic_;


  }

  void  LinearKellerMiksis::CompDeriv(const Double &t,
                                const StdVector<Double> &y,
                                StdVector<Double> &dydt){
    ENTER_FCN( "LinearKellerMiksis::CompDeriv" );

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // "Dimesionsbehaftet" computation of the right hand sight of the 
    // LinearKellerMiksis equation y=dydt(y);
    // Standard input/output values used in bubbledriver.cc and 
    // solveStepAcousticBubble.cc are by now dimensionless.
    // Therefore a transformation to the real values is nessecary.
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    // For testing only:-------------------------------------------

//      Double tdim;
//      StdVector<Double> ydim(2);
//      Double ampl, freq;


//      ydim[0] = y[0] * RadiusInit_ ;
//      ydim[1] = y[1] * sqrt(pStatic_ / density_); 
//      tdim = t * RadiusInit_ / (sqrt( pStatic_ / density_));
//      printf("%10.6e \t %10.6e \t %10.6e \n",tdim, ydim[0],ydim[1]);



//      ampl  = 24.0e5;            //Anregungsamplitude
//      freq  = 27000.0;           //Anregungsfrequenz
//     // ampl  = ampl / pStatic_;
//     // freq  = freq *  RadiusInit_ / (sqrt( pStatic_ / density_));
//      p_    = -ampl *sin(2 * (PI) *freq * (t));                        
//      dpdt_ = -ampl *(2 *(PI) * freq)* cos(2 * (PI) *freq * (t));
    // -------------------------------------------------------------

    // Backtransformation +++++++++++++++++
    Double tdim;
    StdVector<Double> ydim(2);

    ydim[0] = y[0] * RadiusInit_ ;
    ydim[1] = y[1] * sqrt(pStatic_ / density_); 

    Double p,dpdt;
    p = p_ * pStatic_;
    dpdt = dpdt_ * pStatic_/RadiusInit_  * sqrt(pStatic_ / density_);
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    StdVector<Double> yklein(2);
    yklein[0]=ydim[0]- RadiusInit_;
    yklein[1]=ydim[1];
    //NOCHMAL UEBERLEGEN; OB ICH HIER RICHTIG MACHE

    Double omega02, alpha;

    dydt[0] = ydim[1];
    
    omega02 = ( 3.0 * polytrop_ * pStatic_ +  2.0 * surfacTen_ / RadiusInit_ * ( 3.0 * polytrop_ - 1.0)) 
      / ( density_ * RadiusInit_ * RadiusInit_);

    alpha = 4.0 * viscosity_ / ( density_ * RadiusInit_ * RadiusInit_ ) - omega02 * RadiusInit_ / sonicVel_;

    dydt[1] = - p / density_ / RadiusInit_ - alpha * yklein[1] - omega02 * yklein[0]; 




    // Forwardtransformation +++++++++++++++++
    dydt[0] = dydt[0]/ sqrt(pStatic_ / density_); 
    dydt[1] = dydt[1] / ( pStatic_ / ( density_ * RadiusInit_ )); 
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  }

  void  LinearKellerMiksis::Jacobi(StdVector<Double> &y,
   			      Matrix<Double> &dfdy,
   			      Double &t){
  
    ENTER_FCN( "LinearKellerMiksis::Jacobi" );

    Double omega02, alpha;

    omega02 = ( 3.0 * polytrop_ * pStatic_ +  2.0 * surfacTen_ / RadiusInit_ * ( 3.0 * polytrop_ - 1.0)) 
      / ( density_ * RadiusInit_ * RadiusInit_);

    alpha = 4.0 * viscosity_ / ( density_ * RadiusInit_ * RadiusInit_ ) - omega02 * RadiusInit_ / sonicVel_;

    dfdy(0,0) = 0;
    dfdy(0,1) = 1;
    dfdy(1,0) = -omega02;
    dfdy(1,1) = -alpha;
 
  } 



 
} // end of namespace

