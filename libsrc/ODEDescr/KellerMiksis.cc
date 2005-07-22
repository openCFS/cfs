#include <fstream>
#include <iostream>
#include <string>

#include <cmath>     


#include "KellerMiksis.hh"


namespace CoupledField
{

  KellerMiksis::KellerMiksis(Double  RadiusInit, 
                             Double  density,
                             Double  sonicVel,
                             Double  pStatic,
                             Double  pVapour,
                             Double  surfacTen,
                             Double  polytrop,
                             Double  viscosity) {
    ENTER_FCN( "KellerMiksis::KellerMiksis" );
    RadiusInit_ = RadiusInit;
    density_    = density;
    sonicVel_   = sonicVel;
    pStatic_    = pStatic;
    pVapour_    = pVapour;
    surfacTen_  = surfacTen;
    polytrop_   = polytrop;
    viscosity_  = viscosity;



  }

  void  KellerMiksis::CompDeriv(const Double &t,
                                const StdVector<Double> &y,
                                StdVector<Double> &dydt){
    ENTER_FCN( "KellerMiksis::CompDeriv" );


    dydt[0]  = y[1];

    dydt[1]  = -3.0 / 2.0 * y[1] * y[1] * ( 1.0 - y[1] / ( 3.0 * sonicVel_ ));


    dydt[1] += 1.0 / density_ * ( pStatic_ - pVapour_ 
                                  + 2.0 * surfacTen_ / RadiusInit_ )
      * std::pow(( RadiusInit_ / y[0] ), ( 3.0 * polytrop_ )) 
      * ( 1.0 + ( 1.0 - 3.0 * polytrop_ ) * y[1] / sonicVel_);

    dydt[1] -= ( 1.0 + y[1] / sonicVel_ ) / density_ 
      *( pStatic_ + p_ - pVapour_); 

    dydt[1] -= 2.0 * surfacTen_ / (density_ * y[0] );

    dydt[1] -= 4.0 * viscosity_ * y[1] / ( density_ * y[0] );

    dydt[1] -= y[0] / ( density_ * sonicVel_ ) * dpdt_;    

    dydt[1] = dydt[1] / (( 1.0 - ( y[1] / sonicVel_ )) * y[0]
                         + 4.0 * viscosity_ / ( density_ * sonicVel_));

 

    //-------------------------------------------------------------
    //-------------------------------------------------------------
    // Test problem for ODE-Solver and rough structure of class

    //------- testproblem 1 (for Euler and RKF45)
    // dydt[1] = - y[0] * PI * PI;
        
    // Necessary set initial conditions in bubbledriver to
    // RADIUSINIT = 0; VELINIT = PI;
    // Watch out for print commando
    // Watch out for stepsize in Solvecall in bubbledriver
    // t0=0 tStop=2
    // solution y(0) = sin (pi*t) 
    //------------------------------------------

    //------- testproblem 2 ( only RKF45)
    // Tree body problem, 
    // startpoint and endpoint need to be the same in plot y(0):y(1)

    // Double mu_str, help1, help2;
    // Double MU;

    // MU      =  0.0122774711;
    // mu_str = 1. - MU;

    // help1 = ( y[0] + MU )*( y[0] + MU ) + y[1]*y[1];
    // help1 = help1 * help1 * help1;
    // help1 = sqrt( help1 );

    // help2 = ( y[0] - mu_str )*( y[0] - mu_str ) + y[1]*y[1];
    // help2 = help2 * help2 * help2;
    // help2 = sqrt( help2 );

    // dydt[0] = y[2];
    // dydt[1] = y[3];
    // dydt[2] = y[0]+2.*y[3]-mu_str*(y[0]+MU)/help1-MU*(y[0]-mu_str)/help2;
    // dydt[3] = y[1]-2.*y[2]-mu_str*y[1]/help1-MU*y[1]/help2;
 
    // Set initial condition in bubbledriver:
    // Testanfangwerte für Dreikörperproblem
    // y_.Push_back(0.994);
    // y_.Push_back( 0.0);
    // y_.Push_back(0.0);
    // y_.Push_back(-2.001585106);
    // Watch out for print commando in bubbledriver
    // t0 = 0; tStop   = 17.0652166; used 6000 steps
    //-------------------------------------------------------------


  }


 
} // end of namespace

