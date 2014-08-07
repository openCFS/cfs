#include "H1Elems.hh"

#include <algorithm>

namespace CoupledField {


  // ========================================================================
  //  FeH1
  // ========================================================================
  
  FeH1::FeH1() {
    hasICModes_ = false;
    completeType_ = SERENDIPITY_TYPE;
  }
  
  FeH1::~FeH1() {
    // clear map with pre-computed shape function
    shapeFncsAtIp_.clear();
    shapeFncsAtIp_.clear();
    icModesAtIp_.clear();
    icModesDerivsAtIp_.clear();
  }
  
  void FeH1::GetShFnc( Vector<Double>& S, const LocPoint& lp,
                       const Elem* ptElem,  UInt comp ){
    //check if the shape function is already computed
    if(shapeFncsAtIp_.find(lp.number) == shapeFncsAtIp_.end() || comp !=1 ){

      CalcShFnc( S, lp.coord, ptElem, comp);
      
      //add them to the map only, if we are allowed to!
      if( preComputShFnc_ && lp.number != LocPoint::NOT_SET ) {
        shapeFncsAtIp_[lp.number] = S;    
      }
      
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
  
  // ========================================================================
  //  Incompatible Modes
  // ========================================================================
  void FeH1::GetShFncICModes( Vector<Double>& S, const LocPoint& lp,
                              const Elem* ptElem,  UInt comp ){
    //check if the shape function is already computed
    if(icModesAtIp_.find(lp.number) == icModesAtIp_.end() || comp !=1 ){

      CalcShFncICModes( S, lp.coord, ptElem, comp);

      //add them to the map only, if we are allowed to!
      if( preComputShFnc_ && lp.number != LocPoint::NOT_SET ) {
        icModesAtIp_[lp.number] = S;    
      }

    }else{
      S = icModesAtIp_[lp.number];
    }
  }

  void FeH1::GetGlobDerivShFncICModes( Matrix<Double>& deriv, const LocPointMapped& lpm,
                                       const Elem* elem, UInt comp ){
    // Get local derivative
    Matrix<Double> locDeriv;

    //check if the shfunction is already computed
    if(icModesDerivsAtIp_.find(lpm.lp.number) == 
        icModesDerivsAtIp_.end() || comp !=1 ){
      CalcLocDerivShFncICModes( locDeriv, lpm.lp.coord, elem, comp);

      //add them to the map only, if we are allowed to!
      if( preComputShFnc_ && lpm.lp.number != LocPoint::NOT_SET ) {
        icModesDerivsAtIp_[lpm.lp.number] = locDeriv;
      }
    }else{
      locDeriv = icModesDerivsAtIp_[lpm.lp.number];
    }
    deriv = locDeriv * lpm.jacInv;
  }

  void FeH1::GetLocDerivShFncICModes( Matrix<Double>& deriv, const LocPoint& lp,
                                      const Elem* elem, UInt comp  ) {
    if(icModesDerivsAtIp_.find(lp.number) == icModesDerivsAtIp_.end()){
      CalcLocDerivShFnc( deriv, lp.coord, elem, comp);
      //add them to the map only, if we are allowed to!
      if( preComputShFnc_ && lp.number != LocPoint::NOT_SET) {
        icModesDerivsAtIp_[lp.number] = deriv;
      }
    }else{
      deriv = icModesDerivsAtIp_[lp.number];
    }
  }

  } // namespace CoupledField
