// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "HCurlElems.hh"


// include calculation of Legendre polynomials
#include "FeBasis/Polynomials.hh"

namespace CoupledField {


  // ========================================================================
//  FeCurl
  // ========================================================================
  
  FeHCurl::FeHCurl() {
  }
  
  FeHCurl::~FeHCurl() {
  }
  
  void FeHCurl::GetShFnc( Matrix<Double>& shape, LocPointMapped& lpm,
                          const Elem* elem, UInt comp ) {
    
    // Perform local->global gradient transformation
    Matrix<Double> locShape;
    this->CalcLocShFnc(locShape, lpm, elem, comp);
    shape =  Transpose(lpm.jacInv) * locShape;
  }
  
  void FeHCurl::GetCurlShFnc( Matrix<Double>& curl, LocPointMapped& lpm,
                         const Elem* elem, UInt comp ) {
    // Perform local->global curl transformation
    Matrix<Double> locCurl;    
    this->CalcLocCurlShFnc( locCurl, lpm, elem, comp );
    curl = lpm.jac * locCurl;
    curl *= ( 1.0 / std::abs(lpm.jacDet) );
  }
      

  
  
} // namespace CoupledField
