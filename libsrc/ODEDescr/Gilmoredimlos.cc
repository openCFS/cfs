#include <fstream>
#include <iostream>
#include <string>

#include <cmath>     


#include "Gilmore.hh"


namespace CoupledField
{

  Gilmore::Gilmore(Double  RadiusInit, 
		   Double  density,
		   Double  sonicVel,
		   Double  pStatic,
		   Double  pVapour,
		   Double  surfacTen,
		   Double  polytrop,
		   Double  viscosity) {
    ENTER_FCN( "Gilmore::Gilmore" );
    RadiusInit_ = RadiusInit;
    density_    = density;
    sonicVel_   = sonicVel;
    pStatic_    = pStatic;
    pVapour_    = pVapour;
    surfacTen_  = surfacTen;
    polytrop_   = polytrop;
    viscosity_  = viscosity;

    n_ = 7.0;

    A_ = 3.001e8;                                    // Philipp and Dähnke
    //A_ = sonicVel_ * sonicVel_ * density_ / n_;          // meine Vermutung
 
    B_ = 3e8;                                        // Philipp and Dähnke
    //B_ = A_ - pStatic_;                                         



    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //Compute dimensionless constants need in CompDeriv
    D_	= 2 * surfacTen_ / (RadiusInit_ * pStatic_);
    v_  = sonicVel_ / (sqrt( pStatic_ / density ));
    mu_ = viscosity_ / ( pStatic_ * RadiusInit_) * sqrt( pStatic_ / density_); 
    a_  = A_ / pStatic_;
    b_  = B_ / pStatic_;

  }



  void  Gilmore::CompDeriv(const Double &t,
			   const StdVector<Double> &y,
			   StdVector<Double> &dydt){
    ENTER_FCN( "Gilmore::CompDeriv" );

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Dimensionless variables;
    Double x;
    Double dxdtau;
    Double p;
    Double dpdtau;
    // tau = t / RadiusInit_ * sqrt(pStatic_ / density_);

    //  x       = y[0] / RadiusInit_;
    //dxdtau  = y[1] / ( sqrt( pStatic_ / density_));
    //p       = p_ / pStatic_;
    //dpdtau  = dpdt_ * RadiusInit_ / pStatic_ / (sqrt( pStatic_ / density_));

    x      = y[0];
    dxdtau  = y[1];
    p      = p_;
    dpdtau = dpdt_;

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    // Dimensionless Gilmore model without pVapour


    // dimensionles pressure at bubble wall
    Double pR;
    // Dimensionles pressure far away from bubble
    Double pInf;
    // Dimensionless enthalpy
    Double H;
    // Dimensionless velocity of sound of the mixture
    Double sonicVelMix;



    Double temp1, temp2, temp3, temp4, temp5;

    pR = ( 1.0 + D_) /  (std::pow((x ), ( 3.0 * polytrop_ ))) - D_ / x 
      - 4 * mu_ / x * dxdtau;

    pInf = 1 + p;





    temp1  = std::pow( (pR + b_)  , (( n_ - 1.0 ) / n_ ));
    temp2  = std::pow( (pInf + b_)  , (( n_ - 1.0 ) / n_ )); 


    H = std::pow( a_ , ( 1.0 / n_ )) *  n_ / ( n_ - 1.0 ) * (temp1 - temp2);

    sonicVelMix = sqrt( v_ * v_ + (n_ - 1.0) * H);

    dydt[0]  = dxdtau;

    temp3   =  ( 1.0 + D_) /  (std::pow((x ), ( 3.0 * polytrop_ ))) 
      * (- 3.0) * polytrop_ * dxdtau / x  ;

    temp3 += D_ / ( x * x ) * dxdtau;
 
    temp3   += 4 * mu_ / ( x * x)  * dxdtau * dxdtau;
 
    temp3   *= std::pow( (pR + b_) , ((- 1.0) / n_ ));
   
    temp4   = std::pow( ( pInf + b_ ) ,( (- 1.0) / n_ )) * dpdtau;
 
    dydt[1] = ( temp3 - temp4 ) * std::pow( a_ , ( 1.0 / n_ ));
 
    dydt[1] *= x / sonicVelMix * ( 1.0 - (dxdtau / sonicVelMix));

    dydt[1] -= 3.0 / 2.0 * dxdtau * dxdtau 
      * ( 1.0 - (dxdtau  / ( 3.0 * sonicVelMix )));

    dydt[1] += ( 1.0 + dxdtau / sonicVelMix ) * H;

    temp5   = ( 1.0 - (dxdtau / sonicVelMix)) * 
      ( x + std::pow( a_ , ( 1.0 / n_ )) * 4 * mu_ / sonicVelMix 
	* std::pow( (pR + b_) , ((- 1.0) / n_ )));
					       
    dydt[1] = dydt[1] / temp5;


    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    // Restoring dimension 
    //    dydt[0] = dydt[0] * (sqrt( pStatic_ / density_));

    //dydt[1] = dydt[1] * pStatic_ / RadiusInit_ / density_;



  }


 
} // end of namespace

