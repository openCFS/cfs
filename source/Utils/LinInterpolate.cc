#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>

#include "LinInterpolate.hh"

namespace CoupledField
{ 
  LinInterpolate::LinInterpolate(std::string nlFileName, MaterialType matType )
    : ApproxData(nlFileName,matType)
  {


  }

  LinInterpolate::~LinInterpolate()
  {

  }

  Double LinInterpolate::EvaluateFunc( Double xEntry ) {

    Double yValue = 0.0;

    // get index of last element
    const UInt kend = x_.GetSize() - 1;

    // if coordinate is out of bounds, return boundary value 
    // (i.e.first or last) 
    if ( xEntry > x_[kend] ) {
      yValue = y_[kend];
    }
    else if ( xEntry < x_[0] ) {
      yValue = y_[0];
    }
    else {
      UInt klo,khi,k;
      klo=0;
      khi=kend;
      // We will find the right place in the table by means of bisection.
      //  klo and khi bracket the input value of xEntry
      while (khi-klo > 1) {
        k=(khi+klo) >> 1; // binary right shift
        if (x_[k] > xEntry)
          khi=k;
        else
          klo=k;
      }

      // size of x interval
      Double dxVal = x_[khi] - x_[klo];

      // The x-values must be distinct!
      if (dxVal == 0.0) {
        EXCEPTION("You cannot have two equal x values!" );
      }

      // relative distance of xEntry to x-Value bounds
      Double a = ( x_[khi] - xEntry )/dxVal;
      Double b = ( xEntry - x_[klo] )/dxVal;

      //linear interpolation
      yValue = a * y_[klo] + b * y_[khi];
    }

    return yValue;
  }


  Double LinInterpolate::EvaluateFuncInv(double inVal)
  {
   
    if ( inVal < y_[0] )
      EXCEPTION("Wrong evaluation: input is smaller as defined in nonlinear file");

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
      EXCEPTION("Wrong evaluation: input is smaller as defined in nonlinear file");

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
