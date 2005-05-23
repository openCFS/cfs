#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>

#include "LinInterpolate.hh"

namespace CoupledField
{ 
  LinInterpolate :: LinInterpolate(std::string nlFileName)
    : ApproxData(nlFileName)
  {


  }

  LinInterpolate :: ~LinInterpolate()
  {

  }


  Double LinInterpolate::EvaluateFuncInv(double inVal)
  {
    ENTER_FCN( "LinInterpolate::EvaluateFuncInv" );
  
    Integer     i;
 
    if ( inVal < y[0] )
      Error("Wrong evaluation: input is smaller as defined in nonlinear file",__FILE__,__LINE__);

    //if inVal is larger as defined, return the last value
    if ( inVal > y[nummeas-1] ) 
      return x[nummeas-1];

    //loop over array 
    Double yPrev, yAfter, xPrev, xAfter;

    yAfter = y[0];
    for ( Integer k=1; k<nummeas; k++ ) {
      yPrev = yAfter;
      yAfter = y[k];

      xPrev   = x[k-1];
      xAfter  = x[k];

      if (inVal == yAfter) return x[k];

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
    ENTER_FCN( "LinInterpolate::EvaluateFuncInv" );
  
    Integer     i;
 
    if ( inVal < y[0] )
      Error("Wrong evaluation: input is smaller as defined in nonlinear file",__FILE__,__LINE__);

    //if inVal is larger as defined, return the last value
    if ( inVal > y[nummeas-1] ) 
      return x[nummeas-1];

    //loop over array 
    Double yPrev, yAfter, xPrev, xAfter;

    yAfter = y[0];
    for ( Integer k=1; k<nummeas; k++ ) {
      yPrev = yAfter;
      yAfter = y[k];

      xPrev   = x[k-1];
      xAfter  = x[k];

      if (inVal == yAfter) return x[k];

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
