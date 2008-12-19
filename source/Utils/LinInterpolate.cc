// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>

#include "LinInterpolate.hh"

namespace CoupledField
{ 
  LinInterpolate::LinInterpolate(std::string nlFileName, ApproxCurveType curveType )
    : ApproxData(nlFileName,curveType)
  {


  }

  LinInterpolate::~LinInterpolate()
  {

  }


  Double LinInterpolate::EvaluateFuncInv(double inVal)
  {
   
    if ( inVal < y_[0] )
      Error("Wrong evaluation: input is smaller as defined in nonlinear file",__FILE__,__LINE__);

    //if inVal is larger as defined, return the last value
    if ( inVal > y_[numMeas_-1] ) 
      return x_[numMeas_-1];

    //loop over array 
		// initialize or receive compiler warning
    Double yPrev = 0.0, yAfter = 0.0, xPrev = 0.0, xAfter = 0.0;

    yAfter = y_[0];
    for ( UInt k=1; k<numMeas_; k++ ) {
      yPrev = yAfter;
      yAfter = y_[k];

      xPrev   = x_[k-1];
      xAfter  = x_[k];

      if (inVal == yAfter) return x_[k];

      if ( inVal < yAfter) break;
     
    }
 
    // linear interpolation
    Double deltay, deltax, erg;
    deltay = yAfter - yPrev;
    deltax = xAfter - xPrev;
 
    if ( abs(deltay) < 1e-12 ) {
      erg= xPrev;
    }
    else {
      erg = xPrev + ( deltax / deltay ) * inVal;
    }
 
    return erg;
 
  }


  Double LinInterpolate::EvaluatePrimeInv(double inVal)
  {
  
    if ( inVal < y_[0] )
      Error("Wrong evaluation: input is smaller as defined in nonlinear file",__FILE__,__LINE__);

    //if inVal is larger as defined, return the last value
    if ( inVal > y_[numMeas_-1] ) 
      return x_[numMeas_-1];

    //loop over array 
		// initialize or receive compiler warning
    Double yPrev = 0.0, yAfter = 0.0, xPrev = 0.0, xAfter = 0.0;

    yAfter = y_[0];
    for ( UInt k=1; k<numMeas_; k++ ) {
      yPrev = yAfter;
      yAfter = y_[k];

      xPrev   = x_[k-1];
      xAfter  = x_[k];

      if (inVal == yAfter) return x_[k];

      if ( inVal < yAfter) break;
     
    }
 
    // linear interpolation
    Double deltay, deltax, erg;
    deltay = yAfter - yPrev;
    deltax = xAfter - xPrev;
 
    if ( abs(deltay) < 1e-12 ) {
      erg= 0;
    }
    else {
      erg =  deltax / deltay;
    }

    std::cout << "hp=" << erg << std::endl; 
    return erg;
 
  }

}
