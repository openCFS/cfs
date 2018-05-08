#include <iostream>
#include <fstream>


#include "Hysteresis.hh"

#include "Utils/tools.hh"

namespace CoupledField
{ 
  Hysteresis::Hysteresis(Integer numElem)
  {

    //numElements_ = numElem;
    Hysteresis(numElem, 1.0,1.0, 0.0,0.0,0.0, false);
  }

  Hysteresis::Hysteresis(Integer numElem, Double XSaturated, Double YSaturated, Double anhystA, Double anhystB, Double anhystC, bool anhystOnly)
  {
    numElements_ = numElem;
    anhyst_A_ = anhystA;
    anhyst_B_ = anhystB;
    anhyst_C_ = anhystC;
    anhystOnly_ = anhystOnly;
    XSaturated_ = XSaturated;
    YSaturated_ = YSaturated;
  }
  
  
  Hysteresis::~Hysteresis()
  {
  }
  
	Double Hysteresis::bisectForAnhyst(Double Ytarget, 
		Double Xdown, Double Xup, Double Poffset, Double eps_mu, Double tol){
		
		return XSaturated_*bisectForAnhyst_normalized(Ytarget/YSaturated_, 
			Xdown/XSaturated_, Xup/XSaturated_, Poffset/YSaturated_, XSaturated_*eps_mu/YSaturated_, tol/XSaturated_);
		
	}
	
	
	Double Hysteresis::bisectForAnhyst_normalized(Double Ytarget_normalized, 
		Double Xdown_normalized, Double Xup_normalized, Double Poffset_normalized, Double eps_mu_normalized, Double tol){
		
		/*
		 *	Solve
		 *		YTarget = eps_mu*X + anhystPart(X) + Poffset
		 *	via bisection with 
		 *		X in range [Xstart, Xend]
		 * 
		 *	Poffset is value of Hyst operator
		 *	
		 *	Usage:
		 *	case1: pure anhysteretic case > no polarization, i.e. Poffset = 0
		 *	case2: saturation case > polarization maximal/minimal, i.e. Poffset = +/-1
		 * 
		 *	Note1: all input and output values have to be saturated accordingly!
		 */
		Double Xout_normalized, Xmid_normalized;
		Double resUp, resDown, resMid;
		
		resUp = Poffset_normalized +
			(eps_mu_normalized*Xup_normalized + evalAnhystPart_normalized(Xup_normalized)) - 
			Ytarget_normalized;
		
		resDown = Poffset_normalized +			
			(eps_mu_normalized*Xdown_normalized + evalAnhystPart_normalized(Xdown_normalized)) -
			Ytarget_normalized;
		
    if(resUp*resDown > 0){
      EXCEPTION("Solution not in expected interval!");
    }
    if(abs(resUp) < tol){
      Xout_normalized = Xup_normalized;
    } else if(abs(resDown) < tol){
      Xout_normalized = Xdown_normalized;
    } else {
      UInt maxIter = 50;
      Xmid_normalized = 0.0;
      for(UInt i = 0; i < maxIter; i++){
        Xmid_normalized = (Xup_normalized + Xdown_normalized)/2.0;
        resMid = Poffset_normalized +
					(eps_mu_normalized*Xmid_normalized + evalAnhystPart_normalized(Xmid_normalized)) -
					Ytarget_normalized;
        
        if(abs(resMid) < tol){
          break;
        }
        if(resMid*resDown > 0){
          Xdown_normalized = Xmid_normalized;
          resDown = resMid;
        } else {
          Xup_normalized = Xmid_normalized;
          resUp = resMid;
        }
      }
      Xout_normalized = Xmid_normalized;
//      LOG_DBG(scalpreisachInversion) << "bisection failed; reamining residual: " << resMid;
    }
    return Xout_normalized;

	}
//	
//	
//	
//  Double Hysteresis::bisectForSaturation(Double Yin, Double eps_mu, Double tol, bool negSaturation){
////    if(negSaturation){
////      LOG_DBG(scalpreisachInversion) << "perform bisection for negative saturation case";
////    } else {
////      LOG_DBG(scalpreisachInversion) << "perform bisection for postive saturation case";
////    }
//    // solve
//    // yIn = Ysaturated_ + eps_mu*xOut + YSaturated_*(anhyst_A_*std::atan(anhyst_B_*xOut/Xsaturated_) + anhyst_C_*xOut/Xsaturated_)
//    // for xOut
//    // > solution should be between XSaturated (lower bound) and Xout = (Yin - YSaturated_)/eps_mu (upper bound)
//    // > use bisection
//    Double resUp, resDown, resMid;
//    Double Xup,Xdown,Xmid,Xout;
//    Double factor;
//    if(negSaturation){
//      Xup = (Yin + YSaturated_)/eps_mu;
//      Xdown = -XSaturated_;
//      factor = -1.0;
//    } else {
//      Xup = (Yin - YSaturated_)/eps_mu;
//      Xdown = XSaturated_;
//      factor = 1.0;
//    }
//
//    resUp = factor*YSaturated_ + eps_mu*Xup + YSaturated_*evalAnhystPart_normalized(Xup/XSaturated_) - Yin;
//    // YSaturated_*(anhyst_A_*std::atan(anhyst_B_*Xup/XSaturated_) + anhyst_C_*Xup/XSaturated_) - Yin;
//    resDown = factor*YSaturated_ + eps_mu*Xdown + YSaturated_*evalAnhystPart_normalized(Xdown/XSaturated_) - Yin;
//    // YSaturated_*(anhyst_A_*std::atan(anhyst_B_*Xdown/XSaturated_) + anhyst_C_*Xdown/XSaturated_) - Yin;
//    if(resUp*resDown > 0){
//      EXCEPTION("Solution not in expected interval!");
//    }
//    if(abs(resUp) < tol){
//      Xout = Xup;
//    } else if(abs(resDown) < tol){
//      Xout = Xdown;
//    } else {
//      UInt maxIter = 50;
//      Xmid = 0;
//      for(UInt i = 0; i < maxIter; i++){
//        Xmid = (Xup + Xdown)/2.0;
//        resMid = factor*YSaturated_ + eps_mu*Xmid + YSaturated_*evalAnhystPart_normalized(Xmid/XSaturated_) - Yin;
//        
//        if(abs(resMid) < tol){
//          break;
//        }
//        if(resMid*resDown > 0){
//          Xdown = Xmid;
//          resDown = resMid;
//        } else {
//          Xup = Xmid;
//          resUp = resMid;
//        }
//      }
//      Xout = Xmid;
////      LOG_DBG(scalpreisachInversion) << "bisection failed; reamining residual: " << resMid;
//    }
//    return Xout;
//  }
  
}
