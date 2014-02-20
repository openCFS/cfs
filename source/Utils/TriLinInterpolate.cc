#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>

#include "TriLinInterpolate.hh"
#include "DataInOut/Logging/LogConfigurator.hh"  

namespace CoupledField
{ 

DECLARE_LOG(trilinappx)  
DEFINE_LOG(trilinappx, "trilinappx")

  TriLinInterpolate::TriLinInterpolate(std::string nlFileName, MaterialType matType )
    : ApproxData(nlFileName,matType, 3)
  {


  }

  TriLinInterpolate::~TriLinInterpolate()
  {

  }

  Double TriLinInterpolate::EvaluateFunc( Double x0Entry , Double x1Entry, Double x2Entry) {
    Double zValue = 0.0;
    EXCEPTION("not implemented yet");
    return zValue;
  }

  Double TriLinInterpolate::EvaluateFunc( Double x0Entry , Double x1Entry) {
    Double zValue = 0.0;
    EXCEPTION("you need to provide x, y and z for TriLinear interpolation");

    // get index of last element
    const UInt kx0end = x_.GetSize() - 1;
    const UInt kx1end = x1_.GetSize() - 1;

    // if both coordinates are out of bounds, return boundary value 
    // (i.e.first or last) 
    if ( x0Entry > x_[kx0end]  && x1Entry > x1_[kx1end]) {
      zValue = z_[kx0end][kx1end]; 
    }
    else if ( x0Entry > x_[kx0end]  && x1Entry < x1_[0]) {
      zValue = z_[kx0end][0]; 
    }
    else if ( x0Entry < x_[0]  && x1Entry < x1_[0]) {
      zValue = z_[0][0]; 
    }
    else if ( x0Entry < x_[0]  && x1Entry > x1_[kx1end]) {
      zValue = z_[0][kx1end]; 
    } 
    // if one coordinate is out of bound, do linear interpolation
    else if ( x0Entry > x_[kx0end] || x0Entry < x_[0] ) {
      // x is out of bound, interpolate x1
      UInt k1lo, k1hi;
      Double diff1;
      findBracketIndices(x1Entry, x1_, k1lo, k1hi, diff1);
      // relative distance of xEntry to x-Value bounds
      Double a = ( x1_[k1hi] - x1Entry )/diff1;
      Double b = ( x1Entry - x1_[k1lo] )/diff1;

      UInt x0ind = (x0Entry < x_[0]) ? 0 : kx0end;
      //linear interpolation
      zValue = a * z_[x0ind][k1lo] + b * z_[x0ind][k1hi];
    }
    else if ( x1Entry > x1_[kx1end] || x1Entry < x1_[0] ) {
      // x is out of bound, interpolate x1
      UInt k0lo, k0hi;
      Double diff0;
      findBracketIndices(x0Entry, x_, k0lo, k0hi, diff0);
      // relative distance of xEntry to x-Value bounds
      Double a = ( x_[k0hi] - x0Entry )/diff0;
      Double b = ( x0Entry - x_[k0lo] )/diff0;

      UInt x1ind = (x1Entry < x1_[0]) ? 0 : kx1end;
      //linear interpolation
      zValue = a * z_[k0lo][x1ind] + b * z_[k0hi][x1ind];
    }
    else {

      UInt k0lo, k0hi, k1lo, k1hi;
      Double diff0, diff1;
      findBracketIndices(x0Entry, x_, k0lo, k0hi, diff0);
      findBracketIndices(x1Entry, x1_, k1lo, k1hi, diff1);

      // bilinear interpolation
      zValue = 1/(diff0*diff1) * (
               z_[k0lo][k1lo] * (x_[k0hi] - x0Entry) * (x1_[k1hi] - x1Entry) +
               z_[k0hi][k1lo] * (x0Entry - x_[k0lo]) * (x1_[k1hi] - x1Entry) +
               z_[k0lo][k1hi] * (x_[k0hi] - x0Entry) * (x1Entry - x1_[k1lo]) +
               z_[k0hi][k1hi] * (x0Entry - x_[k0lo]) * (x1Entry - x1_[k1lo])
	       );

    }

    LOG_DBG(trilinappx) << "Eval trilinear interpolator at points (" << x0Entry << ", " << x1Entry << "): " << zValue;
    return zValue;
  }

  void TriLinInterpolate::findBracketIndices(const double &x, const Vector<Double> & axis, UInt & klo, UInt & khi, double &diff) {
    const UInt kend = axis.GetSize() - 1;
    UInt k;
    klo=0;
    khi=kend;
    // We will find the right place in the table by means of bisection.
    //  klo and khi bracket the input value of xEntry
    while (khi-klo > 1) {
      k=(khi+klo) >> 1; // binary right shift
      if (axis[k] > x)
        khi=k;
      else
        klo=k;
    }
    // size of x interval
    diff = axis[khi] - axis[klo];

    // The x-values must be distinct!
    if (diff == 0.0) {
      EXCEPTION("You cannot have two equal x values!" );
    }

  }

  Double TriLinInterpolate::EvaluateFunc( Double xEntry ) {

    EXCEPTION("you need to provide x, y and z for TriLinear interpolation");
    Double yValue = 0.0;
    return yValue;
  }


  Double TriLinInterpolate::EvaluateFuncInv(double inVal)
  {
   
     EXCEPTION("not implemented");
     Double erg = -1.;
    return erg;
 
  }


  Double TriLinInterpolate::EvaluatePrimeInv(double inVal)
  {
  
     EXCEPTION("not implemented");
     Double erg = -1.;
    return erg;
 
  }

}
