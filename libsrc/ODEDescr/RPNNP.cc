#include <fstream>
#include <iostream>
#include <string>

#include <cmath>     


#include "RPNNP.hh"


namespace CoupledField
{

  RPNNP::RPNNP(Double  RadiusInit, 
                             Double  density,
                             Double  sonicVel,
                             Double  pStatic,
                             Double  pVapour,
                             Double  surfacTen,
                             Double  polytrop,
                             Double  viscosity) {

    ENTER_FCN( "RPNNP::RPNNP" );

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

  void  RPNNP::CompDeriv(const Double &t,
                                const StdVector<Double> &y,
                                StdVector<Double> &dydt){
    ENTER_FCN( "RPNNP::CompDeriv" );

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // "Dimesionsbehaftet" computation of the right hand sight of the 
    // RPNNP equation y=dydt(y);
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

    double temp1, temp2, temp3, temp4, temp5, temp6, temp7;

    dydt[0] = ydim[1];
    //  std::cout<<dydt[0]<<std::endl;
    
    temp1  = (pStatic_ - pVapour_ + 2.0 *  surfacTen_ / RadiusInit_) * std::pow((RadiusInit_/ydim[0]),(3.0*polytrop_));
    temp1 += (-2.0) * surfacTen_ / ydim[0] - 4.0 * viscosity_ * ydim[1] / ydim[0] + pVapour_;
    
    dydt[1] = - 3.0 / 2.0 * ydim[1] * ydim[1] / ydim[0] + 1.0 / density_ / ydim[0] * ( temp1 - pStatic_ - p);



    // Forwardtransformation +++++++++++++++++
    dydt[0] = dydt[0]/ sqrt(pStatic_ / density_); 
    dydt[1] = dydt [1] / ( pStatic_ / ( density_ * RadiusInit_ )); 
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  }

   void  RPNNP::Jacobi(StdVector<Double> &y,
   			      Matrix<Double> &dfdy,
   			      Double &t){
  
    ENTER_FCN( "RPNNP::Jacobi" );
    std::cerr<< "Jacobi method is not yet implemented for this problem" << std::endl;

 
  } 



 
} // end of namespace

