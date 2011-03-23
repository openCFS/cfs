// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "HCurlElems.hh"


// include calculation of Legendre polynomials
#include "polynomials.hh"

namespace CoupledField {


  // ========================================================================
//  FeCurl
  // ========================================================================
  
  FeHCurl::FeHCurl() {
  }
  
  FeHCurl::~FeHCurl() {
  }
  
  void FeHCurl::GetCurlShFcn( Matrix<Double>& s, const LocPoint& lp, 
                              const Elem* ptElem,  UInt comp) {
    Warning("Implement me");
  }
      
  //! Return global derivative of shape functions
  void FeHCurl::GetGlobDerivShFnc( Matrix<Double>& deriv, LocPointMapped& lp,
                                   const Elem* elem, UInt comp ) {

  }

  
  
} // namespace CoupledField
