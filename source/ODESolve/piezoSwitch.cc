// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <cmath>     

#include "piezoSwitch.hh"
#include "stdio.h"


namespace CoupledField
{


  PiezoSwitchODE::PiezoSwitchODE(UInt numberOfODEs, Double initT, Double alpha) {

    numODEs_ = numberOfODEs; 
    dtPDE_   = initT;
    alpha_   = alpha;
  }

  void  PiezoSwitchODE::CompDeriv(const Double &t,
				 const StdVector<Double> &y,
				 StdVector<Double> &dydt) {
    Double val;
    UInt len1 = y.GetSize();

    dydt.Init();
    for (UInt i=0; i<len1; i++ ) {
      val = std::pow( y[i], alpha_ );
      for (UInt j=0; j<len1; j++ ) {
        if ( i != j) {
          dydt[i] +=  coeffs_[j][i] * std::pow( y[j], alpha_ ) 
            - coeffs_[i][j] * val;
        }
      }
    }

    //    std::cout << "RHS in Rosen:\n" << dydt << std::endl; 

  }


  void  PiezoSwitchODE::Jacobi(StdVector<Double> &y,
   			      Matrix<Double> &dfdy,
   			      Double &t){

  
    UInt len1 = y.GetSize();
    
    //std::cout << "In Rosen: y: \n " << y << std::endl;
    for (UInt i=0; i<len1; i++ ) {
      for (UInt j=0; j<len1; j++ ) {
        if ( i != j) {
          dfdy[i][j] = coeffs_[j][i] * alpha_ * std::pow( y[j], (alpha_-1) );
        }
      }
      //sum = 0.0;
      for (UInt j=0; j<len1; j++ ) {
        if ( i != j) {
          //sum -= coeffs_[i][j];
          dfdy[i][i] -= coeffs_[i][j] * alpha_ * std::pow( y[i], (alpha_-1) );
        }
      }
      //      dfdy[i][i] = sum * alpha_ * std::pow( y[i], (alpha_-1) );
      //coeffs_[i][j] * alpha_ * std::pow( y[i], (alpha_-1) );
    }
     
    // std::cout << "Jacobi in Rosen:\n" << dfdy << std::endl; 

  }
  
} // end of namespace

