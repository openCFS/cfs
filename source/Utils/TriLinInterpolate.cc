#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>

#include "TriLinInterpolate.hh"
#include "DataInOut/Logging/LogConfigurator.hh"  

namespace CoupledField
{ 
DEFINE_LOG(trilinappx, "trilinappx")

  TriLinInterpolate::TriLinInterpolate(std::string nlFileName, MaterialType matType )
    : ApproxData(nlFileName,matType, 3)
  {
    // structure slicesFiles_ contains list of filenames for 2d slices 
    slices_.Resize(numMeas_);

    for (UInt i=0; i<numMeas_; i++) {
      slices_[i] = new BiLinInterpolate(slicesFiles_[i], matType);
    }

  }

  TriLinInterpolate::~TriLinInterpolate()
  {

  }

  Double TriLinInterpolate::EvaluateFunc( Double x0Entry , Double x1Entry, Double x2Entry) {
    Double zValue = 0.0;
    const UInt kx0end = x_.GetSize() - 1;
    UInt k0lo, k0hi;
    double diff0;
    // e.g x2=vds, x1jk
    if ( x0Entry > x_[kx0end]) {
      k0lo = kx0end-1;
      k0hi = kx0end;
      diff0 = x_[kx0end] - x_[kx0end-1];
    }
    else if ( x0Entry < x_[0]) {
      k0lo = 0;
      k0hi = 1;
      diff0 = x_[1] - x_[0];
    }
    else {
      findBracketIndices(x0Entry, x_, k0lo, k0hi, diff0);
    }
    // now we know which slices to ask for their bilinear interpolation
    
    Double c0, c1;
    c0 = slices_[k0lo]->EvaluateFunc(x1Entry, x2Entry);
    c1 = slices_[k0hi]->EvaluateFunc(x1Entry, x2Entry);

    // linear interpolation on  x0
    zValue = c0 * (1 - (x0Entry - x_[k0lo])/diff0) + c1 * (x0Entry - x_[k0lo])/diff0;
      
    LOG_DBG(trilinappx) << "Eval trilinear interpolator at points (" << x0Entry << ", " << x1Entry << ", " << x2Entry << "): " << zValue*factor_;
    return zValue*factor_;
  }

  Double TriLinInterpolate::EvaluateFunc( Double x0Entry , Double x1Entry) {
    Double zValue = 0.0;
    EXCEPTION("you need to provide x, y and z for TriLinear interpolation");
    return zValue;
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
