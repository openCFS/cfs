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

    n_ = 7;

    // There are different method to obtain the constants A and B;
    // it is not yet cleared which one is to be used 
    //A_ = 3.001e8;                                    // Philipp and Dähnke
    //A_ = sonicVel_ * sonicVel_ * density_ / ( n * pStatic_); // Sapozhnikov
    A_ = sonicVel_ * sonicVel_ * density_ / n_;          // meine Vermutung
 
    // B_ = 3e8;                                        // Philipp and Dähnke
    B_ = A_ - pStatic_;                                // Sapozhnikov   
    B_ = A_ -1;                                         // Church


    // log to screen 
    // std::cerr << " constants\n";
    // std::cerr << "RadiusInit_" << RadiusInit_ << " \n";
    // std::cerr << "density_"    << density_    << " \n";
    // std::cerr << "SonicVel_"   << sonicVel_   << " \n";
    // std::cerr << "pStatic_"    << pStatic_    << " \n";
    // std::cerr << "pVapour_"    << pVapour_    << " \n";
    // std::cerr << "surfacTen_"  << surfacTen_  << " \n";
    // std::cerr << "polytrop_"   << polytrop_   << " \n";
    // std::cerr << "viscosity_"  << viscosity_  << " \n";
    // std::cerr << "A_"          << A        _  << " \n";
    // std::cerr << "B_"          << B_          << " \n";
    // std::cerr << "n_"          << n_          << " \n";


  }

  void  Gilmore::CompDeriv(const Double &t,
			   const StdVector<Double> &y,
			   StdVector<Double> &dydt){
    ENTER_FCN( "Gilmore::CompDeriv" );

    Double temp1, temp2, temp3, temp4;



    // Computation of enthalpy with pressure at wall 
    // (pStatic + 2 * surfacTen_ / RadiusInit ) 
    //  * std::pow(( RadiusInit_ / y[0] ), ( 3.0 * polytrop_ )) 
    // pVapour neglegted, so is diffusion;
    temp1   = ( pStatic_ + 2.0 * surfacTen_ / RadiusInit_ ) 
      * std::pow(( RadiusInit_ / y[0] ), ( 3.0 * polytrop_ ));

    temp1  -=  2.0 * surfacTen_ / y[0];

    temp1  -=  4.0 * viscosity_ * y[1] / y[0] + B_;

    H_      = std::pow( temp1 , (( n_ - 1 ) / n_ ));

    temp2   =  std::pow( ( pStatic_ + p_ + B_ ) ,(( n_ - 1 ) / n_ ));

    H_     -=  temp2;
  
    H_  *=  n_ / ( n_ - 1 ) / density_ * std::pow( A_ , ( 1 / n_ ));


    sonicVelMix_ = sqrt ( sonicVel_ * sonicVel_ + ( n_ - 1 ) * H_ );


    dydt[0]  = y[1];

    temp3    = ( pStatic_  + 2.0 * surfacTen_ / RadiusInit_ )
      * std::pow(( RadiusInit_ / y[0] ), ( 3.0 * polytrop_ )) 
      * ( -3.0) * polytrop_ * y[1] / y[0];

    temp3   += 2.0 * surfacTen_ / ( y[0] * y[0] ) * y[1];

    temp3   += 4.0 * viscosity_ *  y[1] * y[1]  / y[0];

    temp3   *= std::pow( temp1 , (- 1 / n_ ));

    temp4    = std::pow( ( pStatic_ + p_ + B_ ) ,( - 1 / n_ )) + dpdt_;

    dydt[1]  = ( temp3 - temp4 ) * std::pow( A_ , ( 1 / n_ )) / density_;

    dydt[1] *= ( 1.0 - y[1] / sonicVelMix_ ) * y[0] /sonicVelMix_;

    dydt[1] -= 3.0 / 2.0 * y[1] * y[1] * ( 1.0 - y[1] 
					   / ( 3.0 * sonicVelMix_ ));

    dydt[1] += ( 1.0 + y[1] / sonicVelMix_ ) * H_;

    dydt[1] = dydt[1] / (( 1.0 - ( y[1] / sonicVel_ )) * y[0] 
			 * ( 1.0 + 4.0 * viscosity_ / ( sonicVelMix_ * y[1])
    			     * std::pow( temp1 , (- 1 / n_ )) 
    			     * std::pow( A_ , ( 1 / n_ )) / density_));



  }


 
} // end of namespace

