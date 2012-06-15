#include "H1Elems.hh"

#include <algorithm>

namespace CoupledField {


  // ========================================================================
  //  FeH1
  // ========================================================================
  
  FeH1::FeH1() {
  }
  
  FeH1::~FeH1() {
    // clear map with pre-computed shape function
    shapeFncsAtIp_.clear();
    shapeFncsAtIp_.clear();
  }
  
  void FeH1::GetShFnc( Vector<Double>& S, const LocPoint& lp,
                 const Elem* ptElem,  UInt comp ){
    
    //check if the shape function is already computed
    if(shapeFncsAtIp_.find(lp.number) == shapeFncsAtIp_.end() || comp !=1 ){
      CalcShFnc( S, lp.coord, ptElem, comp);
      shapeFncsAtIp_[lp.number] = S;
    }else{
      S = shapeFncsAtIp_[lp.number];
    }
  }
  
  void FeH1::GetGlobDerivShFnc( Matrix<Double>& deriv, const LocPointMapped& lpm,
                                const Elem* elem, UInt comp ){
    // Get local derivative
    Matrix<Double> locDeriv;
    
    //check if the shfunction is already computed
    if(shapeFncDerivsAtIp_.find(lpm.lp.number) == 
        shapeFncDerivsAtIp_.end() || comp !=1 ){
      CalcLocDerivShFnc( locDeriv, lpm.lp.coord, elem, comp);
      
      //add them to the map only, if we are allowed to!
      if( preComputShFnc_ && lpm.lp.number != LocPoint::NOT_SET ) {
        shapeFncDerivsAtIp_[lpm.lp.number] = locDeriv;
      }
    }else{
      locDeriv = shapeFncDerivsAtIp_[lpm.lp.number];
    }
    deriv = locDeriv * lpm.jacInv;
  }
  
  void FeH1::GetLocDerivShFnc( Matrix<Double>& deriv, const LocPoint& lp,
                               const Elem* elem, UInt comp  ) {
    if(shapeFncDerivsAtIp_.find(lp.number) == shapeFncDerivsAtIp_.end()){
      CalcLocDerivShFnc( deriv, lp.coord, elem, comp);
      //add them to the map only, if we are allowed to!
      if( preComputShFnc_ && lp.number != LocPoint::NOT_SET) {
        shapeFncDerivsAtIp_[lp.number] = deriv;
      }
    }else{
      deriv = shapeFncDerivsAtIp_[lp.number];
    }
  }

  } // namespace CoupledField
