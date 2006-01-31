#include <fstream>
#include <iostream>
#include <string>
#include <stdio.h>

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
    // A_ = 3.001e8;
    // B_ = 3e8;
    // The following constants can be determined 
    // as in the paper of C. Church
    // or one uses the above values, which are
    // approximations made of experience
    A_ = sonicVel_ * sonicVel_ * density_ / n_;
    B_ = A_ - pStatic_;
  }

  void  Gilmore::CompDeriv(const Double &t,
                           const StdVector<Double> &y,
                           StdVector<Double> &dydt){
    ENTER_FCN( "Gilmore::CompDeriv" );

    Double temp1, temp2, temp3, temp4, temp5, temp6, temp7;

    // For testing only:-------------------------------------------

//         Double ampl, freq;

//         printf("%10.6e \t %10.6e \t %10.6e \n",t, y[0],y[1]);

//         ampl  = 24.0e5;            //Anregungsamplitude
//         freq  = 27000.0;           //Anregungsfrequenz
//         p_    = -ampl *sin(2 * (PI) *freq * (t));                        
//         dpdt_ = -ampl *(2 *(PI) * freq)* cos(2 * (PI) *freq * (t));
    // -------------------------------------------------------------

    dydt[0] = y[1];
    

    temp1  = (pStatic_ - pVapour_ + 2.0 * surfacTen_ / RadiusInit_) 
      * std::pow((RadiusInit_/y[0]),(3.0*polytrop_));
    temp1 += (-2.0) * surfacTen_ / y[0] - 4.0* viscosity_ * y[1] / y[0] + pVapour_ + B_;
    temp2  = std::pow(temp1,((n_-1.0)/n_));
    temp3  = std::pow((pStatic_ + p_ + B_),((n_-1.0)/n_));
    H_      = (temp2 - temp3)*std::pow(A_,(1.0/n_)) / density_ * n_/(n_-1.0);
    C_      = sqrt(sonicVel_ * sonicVel_ + (n_-1.0) * H_);


    temp4  = (pStatic_ - pVapour_ + 2.0 * surfacTen_ / RadiusInit_) 
      *std::pow((RadiusInit_/y[0]),(3.0*polytrop_))*(-3.0) * polytrop_ * y[1] / y[0];
    temp4 += 2.0 * surfacTen_ / (y[0] * y[0]) *y[1] 
      + 4.0* viscosity_ *(y[1] * y[1]) / (y[0] * y[0]);
    temp5  = temp4 * std::pow(temp1,(-1.0/n_));
    temp6  = std::pow((pStatic_ + p_ + B_),(-1.0/n_)) * dpdt_;

    dydt[1]= (temp5 - temp6) * std::pow(A_,(1.0/n_)) 
      / density_ * (1.0-(y[1]/C_)) * y[0]/ C_;
    dydt[1]+= (1.0+ (y[1]/C_)) * H_ - 3.0 / 2.0 * y[1] * y[1]*(1.0 - (y[1]/(3.0*C_)));

    

    temp7  = (1.0 - (y[1] / C_)) * y[0]* 
      (1.0+( 4.0 * viscosity_ / (C_*y[0])
	     *std::pow(temp1,(-1.0/n_)) * std::pow(A_,(1.0/n_))/density_));

    dydt[1]= dydt[1] / temp7;


  }

 
} // end of namespace

