#include <iostream>
#include <fstream>


#include "Hysteresis.hh"

#include "Utils/tools.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

namespace CoupledField
{ 
  DECLARE_LOG(vecpreisachInversion)
  DEFINE_LOG(vecpreisachInversion, "vecpreisachInversion")
  DECLARE_LOG(vecpreisachlinesearch)
  DEFINE_LOG(vecpreisachlinesearch, "vecpreisachlinesearch")

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
    PSaturated_ = YSaturated;
    dim_ = 0;
  }
    
  Hysteresis::~Hysteresis()
  {
  }
    
	Double Hysteresis::bisectForAnhyst(Double Ytarget, 
		Double Xdown, Double Xup, Double Poffset, Double eps_mu, Double tol, Vector<Double> dir, UInt idx){
		
		return XSaturated_*bisectForAnhyst_normalized(Ytarget/PSaturated_, 
			Xdown/XSaturated_, Xup/XSaturated_, Poffset/PSaturated_, XSaturated_*eps_mu/PSaturated_, tol/XSaturated_, dir, idx);
		
	}

	Double Hysteresis::bisectForAnhyst_normalized(Double Ytarget_normalized, 
          Double Xdown_normalized, Double Xup_normalized, Double Poffset_normalized, 
          Double eps_mu_normalized, Double tol, Vector<Double> dir, UInt idx){
		
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
    /*
     * 14.5.2018
     *  Extension needed for Mayergoyz vector model 
     *  -> if dir != zeroVector, Poffset is not assumed to be constant; instead, the 
     *      hyst operator is evaluated 
     *  -> this is required as the Mayergoyz vector model does not return YSst*dir
     *      in general
     */
    
		Double Xout_normalized, Xmid_normalized;
		Double resUp, resDown, resMid;
		
    bool evalRequired = false;
    Vector<Double> Pout = Vector<Double>(dim_);
    Vector<Double> Pin = Vector<Double>(dim_);
    if(dir.NormL2() != 0){
      evalRequired = true;
    }
    int successFlag = 0;
    bool debugOut = false;
    if(evalRequired){
      Pin.Init();
      Pin.Add(Xup_normalized*XSaturated_,dir);
      Pout = computeValue_vec(Pin, idx, false, true, debugOut, successFlag);
      Pout.Inner(dir,Poffset_normalized);
      Poffset_normalized /= PSaturated_;
      // as computeValue_vec already contains anhyst part, we do not have to add it
      resUp = Poffset_normalized +
              (eps_mu_normalized*Xup_normalized) - 
              Ytarget_normalized;
    } else {
      resUp = Poffset_normalized +
              (eps_mu_normalized*Xup_normalized + evalAnhystPart_normalized(Xup_normalized)) - 
              Ytarget_normalized;
		}
    
    if(evalRequired){
      Pin.Init();
      Pin.Add(Xdown_normalized*XSaturated_,dir);
      Pout = computeValue_vec(Pin, idx, false, true, debugOut, successFlag);
      Pout.Inner(dir,Poffset_normalized);
      Poffset_normalized /= PSaturated_;
      
      resDown = Poffset_normalized +			
              (eps_mu_normalized*Xdown_normalized) -
              Ytarget_normalized;
    } else {
      resDown = Poffset_normalized +			
              (eps_mu_normalized*Xdown_normalized + evalAnhystPart_normalized(Xdown_normalized)) -
              Ytarget_normalized;
		}
    
    if(resUp*resDown > 0){
      std::cout << "resUp: " << resUp << std::endl;
      std::cout << "resDown: " << resDown << std::endl;
      
      std::cout << "Xup_normalized: " << Xup_normalized << std::endl;
      std::cout << "Xdown_normalized: " << Xdown_normalized << std::endl;
      
      std::cout << "Ytarget_normalized: " << Ytarget_normalized << std::endl;
      
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
        
        if(evalRequired){
          Pin.Init();
          Pin.Add(Xmid_normalized*XSaturated_,dir);
          Pout = computeValue_vec(Pin, idx, false, true, debugOut, successFlag);
          Pout.Inner(dir,Poffset_normalized);
          Poffset_normalized /= PSaturated_;
          
          resMid = Poffset_normalized +
                  (eps_mu_normalized*Xmid_normalized) -
                  Ytarget_normalized;
        } else {
          
          resMid = Poffset_normalized +
                  (eps_mu_normalized*Xmid_normalized + evalAnhystPart_normalized(Xmid_normalized)) -
                  Ytarget_normalized;
        }
        
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
//    // yIn = Ysaturated_ + eps_mu*xOut + PSaturated_*(anhyst_A_*std::atan(anhyst_B_*xOut/Xsaturated_) + anhyst_C_*xOut/Xsaturated_)
//    // for xOut
//    // > solution should be between XSaturated (lower bound) and Xout = (Yin - PSaturated_)/eps_mu (upper bound)
//    // > use bisection
//    Double resUp, resDown, resMid;
//    Double Xup,Xdown,Xmid,Xout;
//    Double factor;
//    if(negSaturation){
//      Xup = (Yin + PSaturated_)/eps_mu;
//      Xdown = -XSaturated_;
//      factor = -1.0;
//    } else {
//      Xup = (Yin - PSaturated_)/eps_mu;
//      Xdown = XSaturated_;
//      factor = 1.0;
//    }
//
//    resUp = factor*PSaturated_ + eps_mu*Xup + PSaturated_*evalAnhystPart_normalized(Xup/XSaturated_) - Yin;
//    // PSaturated_*(anhyst_A_*std::atan(anhyst_B_*Xup/XSaturated_) + anhyst_C_*Xup/XSaturated_) - Yin;
//    resDown = factor*PSaturated_ + eps_mu*Xdown + PSaturated_*evalAnhystPart_normalized(Xdown/XSaturated_) - Yin;
//    // PSaturated_*(anhyst_A_*std::atan(anhyst_B_*Xdown/XSaturated_) + anhyst_C_*Xdown/XSaturated_) - Yin;
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
//        resMid = factor*PSaturated_ + eps_mu*Xmid + PSaturated_*evalAnhystPart_normalized(Xmid/XSaturated_) - Yin;
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
  
  bool Hysteresis::checkConvergence(Vector<Double>& res, Matrix<Double>& jacT, Double& errorNorm, Double tol){
    // According to Dahmen&Reusken - Numerik partieller DFG
    // the residual is a non-sufficient condition
    // instead we have to check the norm of jacT*res
    Vector<Double> errorVec = Vector<Double>(dim_);
    jacT.Mult(res,errorVec); 
    errorNorm = errorVec.NormL2();
    
    //std::cout << "CheckConvergence" << std::endl;
    //std::cout << "resIn: " << res.ToString() << std::endl;
    //std::cout << "jacT: " << jacT.ToString() << std::endl;
    //std::cout << "ErrorVec: " << errorVec.ToString() << std::endl;
    
    if(errorNorm <= tol){
      return true;
    } else {
      return false;
    }
    
  }
  	
	Double Hysteresis::computeRho(Vector<Double>& xNew, Vector<Double>& xUpdate, 
				Vector<Double>& res, Vector<Double>& resShifted, Matrix<Double>& jac){
		
    std::stringstream computeRhoTrace;
		// According to Dahmen & Reusken, we have to check
		// our increment by computing
		// 
		// rho = (||F(x)||^2 - ||F(x + increment)||^2) /
		//       (||F(x||^2 - ||F(x) + F'(x)*increment||^2)
		//
		// F = residual
		// F' = Jacobian
		Double resNorm = res.NormL2();
		Double resShiftedNorm = resShifted.NormL2();
		Vector<Double> tmp = Vector<Double>(dim_);
		jac.Mult(xUpdate,tmp);
		
		tmp = tmp + res;
		Double tmpNorm = tmp.NormL2();
		Double a = resNorm*resNorm - resShiftedNorm*resShiftedNorm;
		Double b = resNorm*resNorm - tmpNorm*tmpNorm;
		
    LOG_TRACE(vecpreisachInversion) << "Compute rho: a = " << a;
    LOG_TRACE(vecpreisachInversion) << "Compute rho: b = " << b;
    
    computeRhoTrace << "Nominator: ||F(x)||^2 - ||F(x + increment)||^2 = " << a << std::endl;
    computeRhoTrace << "Denominator: ||F(x||^2 - ||F(x) + F'(x)*increment||^2 = " << b << std::endl;
    
		/*
		 * Note: during various tests, rho reached values between
		 *       -200 to 2
		 * > if stepping is too small: rho << 0
		 * > if stepping is too large: rho ~ 1-2
		 * > consider this for creating an appropriate adation of alpha
     * 
     * > EDIT: according to Dahmen Reusken, rho < 0 is never acceptable
     *          and b should always be >= 0
		 */
		Double rho;
		if(b != 0){
			rho = a/b;
		} else {
			rho = a;
      if(a == 0){
        //LOG_DBG(vecpreisach) << "Compute rho: a = 0, too!";
        // rho is actually undefined here (as increment is basically 0)
        // > alpha seems to be way too large (otherwise increment would not be 0)
        // > return rho >= 1 in order to trigger a much smaller alpha
        // > rho >= 10 > reset alpha
        rho = 10.0;
      }
		}
    bool output = false;
//    /*
//     * NEW: check if increment does lead to a usable change in residual at all
//     *      (a/b might be ok, but the actual change in residual (b) might be so small
//     *      that we do not step forward)
//     * 
//     */
//    if(b < 1e-25){
//      b = 0;
//    }
//    
//    if((b/(resNorm*resNorm) < 1e-12)&&(b != 0)){
//      computeRhoTrace << "b/(resNorm*resNorm) = " << b  <<"/"<< (resNorm*resNorm) << " = " << b/(resNorm*resNorm) << " < 1e-12" << std::endl; 
//      computeRhoTrace << " > residual did not really decrease > reset alpha to 1.0" << std::endl; 
//      rho = 10.0;
//      output = true;
//    }
//    
//    
    
    if(b < 0){
      // due to the model itsel, b should always be > 0 (which is the case except when it is 0 an numerical noise comes into play (< 1e-34)
      // however, during linesearch we force the solution to be below saturation; this forcing of the update might actually lead to b < 0
      computeRhoTrace << "b < 0 > SHOULD not happen" << std::endl; 
      output = !true;
    }
    if( (b == 0)&&(a != 0) ){
      computeRhoTrace << "b == 0 > a is taken directly, correct?" << std::endl; 
      output = true;
    }
    if(output == true){
      LOG_TRACE(vecpreisachInversion) << computeRhoTrace.str();
    }
  
		//LOG_DBG(vecpreisach) << "Computed rho: " << rho;
		return rho;
	}
		  
  Integer Hysteresis::checkIncrement(Vector<Double>& xNew, Vector<Double>& xUpdate, 
		Vector<Double>& res, Vector<Double>& resShifted, Matrix<Double>& jac, Double& alpha, int stayBelowSat){
    /*
     * Inspired by "Levenberg-Marquart iteative regularization for th epulse-type impact-force reconstruction"
     * 
     * Trust region method with four regions for rho
     * 
     * rho = (||F(x)||^2 - ||F(x + increment)||^2) /
     *       (||F(x||^2 - ||F(x) + F'(x)*increment||^2)
     * 
     * rho < beta0 : discard update; increase alpha
     * beta0 <= rho < betaLow : accept update; increase alpha
     * betaLow <= rho < betaUp : accept update; keep alpha
     * betaUp <= rho : accept update; decrease alpha
     * 
     * Return value:
     * -1 : dicard update
     *  0 : keep update
     *  1 : keep update
     */
    Double beta0, betaLow, betaUp;
    beta0 = 0.15;
//    betaLow = 0.25;
//    betaUp = 0.75;
    // best for mayergoyz so far
    betaLow = 0.35;
    betaUp = 0.85;
    Double factorUp = 2;
    Double factorDown = 0.5;
    
    //LOG_DBG(vecpreisach) << "Check increment";
    Integer success = 0;
    if((xNew.NormL2() > XSaturated_)&&(stayBelowSat==1)){
      // input yVal is not in saturation (otherwise we would not be here but
      // use the simple case)
      // > new value cannot be in saturation either
      // > current update will definitely not match
      // > problem here: if xNew is in saturation, further comvergence using
      //    the jacobian will not work as the slope will be such that we 
      //    cannot leave this (wrongly obtained) saturation
      // > discard update and try with different alpha again
      //LOG_DBG(vecpreisach) << "Intermediate solution (" << xNew.ToString() << ") has gone into saturation > discard!";
      alpha = alpha*factorUp;
      success = -3;
    } else if((xNew.NormL2() < XSaturated_)&&(stayBelowSat==-1)){
      alpha = alpha/factorUp;
      success = -2;
    } else {
      //std::cout << "Safe case: " << std::endl;
      //std::cout << "Check increment according to residual" << std::endl;
      Double rho = computeRho(xNew, xUpdate, res, resShifted, jac);
      LOG_TRACE(vecpreisachInversion) << "rho: " << rho;
      if(rho < beta0){
        success = -1;
        alpha = alpha*factorUp;
      } else if(rho < betaLow){
        success = 0;
        alpha = alpha*factorUp;
      } else if(rho < betaUp){
        success = 0;
      } else {
        success = 1;
        alpha = alpha*factorDown;
      }
    }
   
    return success;
  }

  Vector<Double> Hysteresis::computeAbsResidualX(Vector<Double>& xVal, Vector<Double>& yVal, Vector<Double>& hystVal, Matrix<Double> mu_inv){
    /*
     *    yVal = mu*xVal + hystVal(xVal)
     * 
     *  > res = xVal + hystVal(xVal)/mu - yVal/mu
     */
    Vector<Double> res = Vector<Double>(dim_);
    Vector<Double> tmp = Vector<Double>(dim_);
    tmp = hystVal;
    tmp.Add(-1.0,yVal);
    mu_inv.Mult(tmp,res);
    res.Add(1.0,xVal);
    
    //    for(UInt i = 0; i < dim_; i++){
    //      // assume that mu is diagonal and has entries != 0
    //      if(mu[i][i] != 0){
    //        res[i] += 1.0/mu[i][i]*hystVal[i];
    //        res[i] += -1.0/mu[i][i]*yVal[i];
    //      }
    //    }
    
    return res;
  }
  
  Vector<Double> Hysteresis::computeResidual(Vector<Double>& xVal, Vector<Double>& yVal, Vector<Double>& hystVal, 
          Matrix<Double> mu, Matrix<Double> mu_inv, bool wrtX, bool relative){
    
    Vector<Double> ret = computeAbsResidualX(xVal, yVal, hystVal, mu_inv);
    if(wrtX){
      /*
       * wrtX = true
       * 
       * > resAbsX = xVal - yVal/mu + hystVal/mu
       * > resRelX = resAbsX/xVal.NormL2()
       */
      if(relative){
        Double xValNorm = xVal.NormL2();
        if(xValNorm != 0){
          ret.ScalarDiv(xValNorm);
        } else {
          WARN("relative residual wrt x requested, but xValNorm == 0; return abs value instead");
        }
      }
    } else {
      /*
       *  wrtX = false
       * 
       * > resAbsY = yVal - hystVal - mu*xVal = -mu*resAbsX
       * > resRelY = reaAbsY/yVal.NormL2()
       * 
       */
      Vector<Double> tmp = Vector<Double>(dim_);
      //ret = ret*(-mu);
      mu.Mult(ret,tmp);
      tmp.ScalarMult(-1.0);
      ret = tmp;
      //      for(UInt i = 0; i < dim_; i++){
      //        ret[i] = -1.0*ret[i]*mu[i][i];
      //      }
      
      if(relative){
        Double yValNorm = yVal.NormL2();
        if(yValNorm != 0){
          ret.ScalarDiv(yValNorm);
        } else {
          WARN("relative residual wrt y requested, but yValNorm == 0; return abs value instead");
        }
      }
    }
    
    return ret;
  }
  
  Matrix<Double> Hysteresis::computeJacobianOfAbsResidualX(Vector<Double>& xVal, Vector<Double>& hystVal, 
          Matrix<Double> mu_inv, Integer operatorIdx, Double sign, UInt implementation, 
					bool overwriteMemory, bool overwriteDirection, int stayBelowSat) {
    
    LOG_DBG(vecpreisachInversion) << " --------- computeJacobianOfAbsResidualX --------- ";
//    LOG_DBG(vecpreisach) << "VecPreisach::computeJacobianOfAbsResidualX";
    if((xVal.NormL2() > XSaturated_)&&(stayBelowSat==1)){
      EXCEPTION("xVal.NormL2() > XSaturated_");
    }
    if((xVal.NormL2() < XSaturated_)&&(stayBelowSat==-1)){
      EXCEPTION("xVal.NormL2() < XSaturated_");
    }
    
    Double deltaX = 0.0;
		Double scal = INV_jacobiResolution_;
    Double deltaXmin = scal*XSaturated_; 
   
//		Double deltaXmin = 1e-8*XSaturated_; 
//    Double scal = 1e-8;//1e-5;
//    
//		
//    bool overwriteMemory = false;
    Vector<Double> xShifted;
    Vector<Double> hystShifted;
    Vector<Double> deltaHyst;
    Matrix<Double> jac = Matrix<Double>(dim_,dim_);
    jac.Init();

    if(implementation == 0){
//      LOG_DBG(vecpreisach) << "Use forward/backward differences for approximation of Jacobian";
      /*
       * Full Jacobian using forward/backward differences
       */
      for(UInt i = 0; i < dim_; i++){
        xShifted = xVal;
        
        //        if( xVal[i] < 0 ){
        //          deltaX = sign*std::min( -scal*XSaturated_, -deltaXmin );
        //        } else {
        //          deltaX = sign*std::max( scal*XSaturated_, deltaXmin );
        //        }
        if( xVal[i] < 0 ){
          deltaX = sign*std::min( scal*xVal[i], -deltaXmin );
        } else {
          deltaX = sign*std::max( scal*xVal[i], deltaXmin );
        }
        
        xShifted[i] += deltaX;
//        if((xShifted[i] > XSaturated_)&&(stayBelowSat==1)){
//          xShifted[i] -= 2*abs(deltaX);
//        } 
//        if((xShifted[i] < -XSaturated_)&&(stayBelowSat==1)){
//          xShifted[i] += 2*abs(deltaX);
//        } 
//        
//        if((xShifted[i] < XSaturated_)&&(stayBelowSat==-1)){
//          xShifted[i] += 2*abs(deltaX);
//        } 
//        if((xShifted[i] > -XSaturated_)&&(stayBelowSat==-1)){
//          xShifted[i] = -XSaturated_-*abs(deltaX);
//        } 
          
        Vector<Double> dXvec = xShifted;
        dXvec.Add(-1.0,xVal);
        LOG_DBG(vecpreisachInversion) << " xVal " << xVal.ToString();
        LOG_DBG(vecpreisachInversion) << " xShifted " << xShifted.ToString(); 
        LOG_TRACE(vecpreisachInversion) << " dXvec " << dXvec.ToString(); 
        int successFlag = 0;
        bool debugOut = false;
        hystShifted = computeValue_vec(xShifted, operatorIdx, overwriteMemory, overwriteDirection, debugOut, successFlag);
        /*
         * Compute Jacobian for residual wrt x
         * 
         * jac_ji = + delta_ji + mu_inv[j][:]*dhystVal/dxVal_i
         */   
        deltaHyst = hystShifted;
        deltaHyst.Add(-1.0,hystVal);
        
        LOG_DBG(vecpreisachInversion) << " hystVal " << hystVal.ToString();
        LOG_DBG(vecpreisachInversion) << " hystShifted " << hystShifted.ToString();
        LOG_TRACE(vecpreisachInversion) << " deltaHyst " << deltaHyst.ToString();
//        Vector<Double> hystValRecomputed = computeValue_vec(xVal, operatorIdx, overwriteMemory, overwriteDirection);
//        LOG_TRACE(vecpreisachInversion) << " hystValRecomputed " << hystValRecomputed.ToString();
        
//        Vector<Double> xOppShift = xShifted;
//        xOppShift[i] -= 2*deltaX;
//        dXvec = xOppShift;
//        dXvec.Add(-1.0,xVal);
//        LOG_TRACE(vecpreisachInversion) << " xOppShift " << xOppShift.ToString(); 
//        LOG_TRACE(vecpreisachInversion) << " dXvec " << dXvec.ToString(); 
//        Vector<Double> hystOppShifted = computeValue_vec(xOppShift, operatorIdx, overwriteMemory, overwriteDirection);
//        LOG_TRACE(vecpreisachInversion) << " hystOppShifted " << hystOppShifted.ToString();
        
        Vector<Double> curCol = Vector<Double>(dim_);
        mu_inv.Mult(deltaHyst,curCol);
        curCol.ScalarDiv(deltaX);
        
        jac[i][i] = 1.0;
        for(UInt j = 0; j < dim_; j++){
          jac[j][i] += curCol[j];
          
          //          //jac[j][i] += (hystShifted[j]-hyst[j])/deltaX/mu; 
          //          // assume matrix mu to be diagonal
          //          // > here we had an error by taking mu[i][i] instead of mu[j][j]
          //          if(mu[j][j]!=0){
          //            jac[j][i] += (hystShifted[j]-hystVal[j])/(xShifted[i]-xShifted_opp[i])/mu[j][j]; 
          //          }
        }
      }
      
      LOG_DBG(vecpreisachInversion) << " Jacobian " << jac.ToString();
      Matrix<Double> jacT = Matrix<Double>(dim_,dim_);
      jac.Transpose(jacT);
      Matrix<Double> jacTjac = Matrix<Double>(dim_,dim_);
      jacT.Mult(jac,jacTjac);
      LOG_DBG(vecpreisachInversion) << " jacTjac " << jacTjac.ToString();
      
      Double det;
      jacTjac.Determinant(det);
      LOG_DBG(vecpreisachInversion) << " det(jacTjac) " << det;
//      if(abs(det) < 1e-14){
//        // repeat with debug output
//        for(UInt i = 0; i < dim_; i++){
//          xShifted = xVal;
//          
//          //        if( xVal[i] < 0 ){
//          //          deltaX = sign*std::min( -scal*XSaturated_, -deltaXmin );
//          //        } else {
//          //          deltaX = sign*std::max( scal*XSaturated_, deltaXmin );
//          //        }
//          if( xVal[i] < 0 ){
//            deltaX = sign*std::min( scal*xVal[i], -deltaXmin );
//          } else {
//            deltaX = sign*std::max( scal*xVal[i], deltaXmin );
//          }
//          
//          xShifted[i] += deltaX;
//          
//          Vector<Double> dXvec = xShifted;
//          dXvec.Add(-1.0,xVal);
//          LOG_TRACE(vecpreisachInversion) << " xVal " << xVal.ToString();
//          LOG_TRACE(vecpreisachInversion) << " xShifted " << xShifted.ToString(); 
//          LOG_TRACE(vecpreisachInversion) << " dXvec " << dXvec.ToString(); 
//          hystShifted = computeValue_vec(xShifted, operatorIdx, overwriteMemory, overwriteDirection,true);
//          /*
//           * Compute Jacobian for residual wrt x
//           * 
//           * jac_ji = + delta_ji + mu_inv[j][:]*dhystVal/dxVal_i
//           */   
//          deltaHyst = hystShifted;
//          deltaHyst.Add(-1.0,hystVal);
//          
//          LOG_TRACE(vecpreisachInversion) << " hystVal " << hystVal.ToString();
//          LOG_TRACE(vecpreisachInversion) << " hystShifted " << hystShifted.ToString();
//          LOG_TRACE(vecpreisachInversion) << " deltaHyst " << deltaHyst.ToString();
//          Vector<Double> hystValRecomputed = computeValue_vec(xVal, operatorIdx, overwriteMemory, overwriteDirection,true);
//          LOG_TRACE(vecpreisachInversion) << " hystValRecomputed " << hystValRecomputed.ToString();
//          
//          Vector<Double> xOppShift = xShifted;
//          xOppShift[i] -= 2*deltaX;
//          dXvec = xOppShift;
//          dXvec.Add(-1.0,xVal);
//          LOG_TRACE(vecpreisachInversion) << " xOppShift " << xOppShift.ToString(); 
//          LOG_TRACE(vecpreisachInversion) << " dXvec " << dXvec.ToString(); 
//          Vector<Double> hystOppShifted = computeValue_vec(xOppShift, operatorIdx, overwriteMemory, overwriteDirection,true);
//          LOG_TRACE(vecpreisachInversion) << " hystOppShifted " << hystOppShifted.ToString();
//          
//          Vector<Double> curCol = Vector<Double>(dim_);
//          mu_inv.Mult(deltaHyst,curCol);
//          curCol.ScalarDiv(deltaX);
//          
//          jac[i][i] = 1.0;
//          for(UInt j = 0; j < dim_; j++){
//            jac[j][i] += curCol[j];
//            
//            //          //jac[j][i] += (hystShifted[j]-hyst[j])/deltaX/mu; 
//            //          // assume matrix mu to be diagonal
//            //          // > here we had an error by taking mu[i][i] instead of mu[j][j]
//            //          if(mu[j][j]!=0){
//            //            jac[j][i] += (hystShifted[j]-hystVal[j])/(xShifted[i]-xShifted_opp[i])/mu[j][j]; 
//            //          }
//          }
//        }
//        
//      }
      
    } else if(implementation == 1){
//      LOG_DBG(vecpreisach) << "Use central differences for approximation of Jacobian";
      /*
       * Full Jacobian using central differences
       */
      Vector<Double> xShifted_opp;
      Vector<Double> hystShifted_opp;
      
      for(UInt i = 0; i < dim_; i++){
        xShifted = xVal;
        
        //        if( xVal[i] < 0 ){
        //          deltaX = sign*std::min( -scal*XSaturated_, -deltaXmin );
        //        } else {
        //          deltaX = sign*std::max( scal*XSaturated_, deltaXmin );
        //        }
        if( xVal[i] < 0 ){
          deltaX = sign*std::min( scal*xVal[i], -deltaXmin );
        } else {
          deltaX = sign*std::max( scal*xVal[i], deltaXmin );
        }
        
        int successFlag = 0;
        bool debugOut = false;
        xShifted[i] += deltaX;
        hystShifted = computeValue_vec(xShifted, operatorIdx, overwriteMemory, overwriteDirection, debugOut, successFlag); 
        xShifted_opp = xVal;
        xShifted_opp[i] -= deltaX;
        hystShifted_opp = computeValue_vec(xShifted_opp, operatorIdx, overwriteMemory, overwriteDirection, debugOut, successFlag);
        
        /*
         * Compute Jacobian for residual wrt x
         * 
         * jac_ji = + delta_ji + mu_inv[j][:]*dhystVal/dxVal_i
         */   
        deltaHyst = hystShifted;
        deltaHyst.Add(-1.0,hystShifted_opp);
        
        Vector<Double> curCol = Vector<Double>(dim_);
        mu_inv.Mult(deltaHyst,curCol);
        curCol.ScalarDiv(2*deltaX);
        
        jac[i][i] = 1.0;
        for(UInt j = 0; j < dim_; j++){
          jac[j][i] += curCol[j];
          
          //          //jac[j][i] += (hystShifted[j]-hyst[j])/deltaX/mu; 
          //          // assume matrix mu to be diagonal
          //          // > here we had an error by taking mu[i][i] instead of mu[j][j]
          //          if(mu[j][j]!=0){
          //            jac[j][i] += (hystShifted[j]-hystVal[j])/(xShifted[i]-xShifted_opp[i])/mu[j][j]; 
          //          }
        }
      }
    } 
    
//    LOG_DBG(vecpreisach) << "Retrieved Jacobian: " << jac.ToString();
    return jac;
  }
  
  Matrix<Double> Hysteresis::computeJacobian(Vector<Double>& xVal, Vector<Double>& yVal, Vector<Double>& hyst, Vector<Double>& resX,
          Matrix<Double> mu, Matrix<Double> mu_inv, Integer operatorIdx, Double sign, bool wrtX, bool relative, 
					UInt implementation, bool overwriteMemory, bool overwriteDirection, int stayBelowSat){

    Matrix<Double> jac = computeJacobianOfAbsResidualX(xVal, hyst, mu_inv, 
												operatorIdx, sign, implementation, overwriteMemory, overwriteDirection, stayBelowSat);
    
    if(wrtX){
      /*
       * wrtX = true
       * 
       * jacAbsX_ij = d resAbsX_i / d x_j
       *            = d ( x_i - y_i/mu + hyst_i/mu ) / d x_j
       *            ~ delta_ij + 1/mu (hyst(x+Dx_j*e_j) - hyst(x))_i / Dx_j 
       *            > computeJacobianOfAbsResidualX(xVal, hyst, mu, idElem, sign);
       * 
       * jacRelX_ij = d resRelX_i / d x_j
       *            
       *            = d (reaAbsX_i / xVal.NormL2()) / d x_j
       * 
       *            = ( xVal.NormL2()* d resAbsX_i / d x_j - resAbsX_i* d xVal.NormL2() / d x_j ) * (1/xVal.NormL2()^2)
       * 
       *            = jacAbsX_ij / xVal.NormL2() - resAbsX_i * xVal_j / (xVal.NormL2()^3) 
       * 
       *            = jacAbsX_ij / xVal.NormL2() - resRelX_i * xVal_j / (xVal.NormL2()^2) 
       * 
       */
      if(relative){
        Double xValNorm = xVal.NormL2();
        if(xValNorm != 0){
          jac = jac*(1.0/xValNorm);
          
          Double xValNormQuad = xValNorm*xValNorm;
          for(UInt i = 0; i < dim_; i++){
            for(UInt j = 0; j < dim_; j++){
              // resX has to be relative wrt x, too!
              jac[i][j] -= (resX[i]*xVal[j])/xValNormQuad;
            }
          }
        } else {
          WARN("Jacobian of relative residual wrt x requested, but xValNorm == 0; return Jacobian of abs value instead");
        }
      }
    } else {
      /*
       * wrtX = false
       * 
       * jacAbsY_ij = d resAbsY_i / d x_j
       *            = d (-mu*resAbsX_i) / d x_j
       *            = -mu * jacAbsX_ij
       * 
       * jacRelY_ij = d resRelY_i / d x_j
       *            = 1/yVal.normL2() * jacAbsY_ij
       *            
       */
      //jac = jac*(-mu);
      Matrix<Double> tmp = Matrix<Double>(dim_,dim_);
      tmp = jac;
      mu.Mult(tmp,jac);
      jac = jac*(-1.0);
      
      if(relative){
        Double yValNorm = yVal.NormL2();
        if(yValNorm != 0){
          jac = jac*(1.0/yValNorm);
        } else {
          WARN("relative residual wrt y requested, but yValNorm == 0; return abs value instead");
        }
      }
    }
    
    return jac;
  }
  
  bool Hysteresis::performLinesearch(Vector<Double>& xVal, Vector<Double>& yVal, Vector<Double>& res, 
		Vector<Double>& xUpdate, Matrix<Double>& jac, Matrix<Double>& jacT, Matrix<Double> mu, Matrix<Double> mu_inv, 
		Integer operatorIdx, bool overwriteMemory, bool overwriteDirection,
		Double& alpha, Double alphaMin, Double alphaMax, bool wrtX, bool relative, 
    UInt& numberOfIterations,Vector<Double>& xStart, Double factorToSat, int stayBelowSat, Vector<Double> sol){
    
    LOG_TRACE(vecpreisachlinesearch) << " --------- START LINESEARCH --------- ";
    LOG_DBG(vecpreisachlinesearch) << "Starting xVal = " << xVal.ToString();
    LOG_DBG(vecpreisachlinesearch) << "Target yVal = " << yVal.ToString();
    LOG_DBG(vecpreisachlinesearch) << "alpha, alphaMin, alphaMax = " << alpha << ", " << alphaMin << ", " << alphaMax << std::endl;
    
    if((xVal.NormL2() > XSaturated_)&&(stayBelowSat==1)){
      EXCEPTION("xInput to Linesearch already above saturation > must not be the case!");
    }
    if((xVal.NormL2() < XSaturated_)&&(stayBelowSat==-1)){
      EXCEPTION("xInput to Linesearch already below saturation > must not be the case!");
    }

		UInt maxIter = 50;
    UInt itCnt = 0;
    
    Matrix<Double> matToInvert = Matrix<Double>(dim_,dim_);
    Matrix<Double> matInverted = Matrix<Double>(dim_,dim_);
    Matrix<Double> jacTjac = Matrix<Double>(dim_,dim_);
    jacT.Mult(jac,jacTjac);
    Vector<Double> jacTres_neg = Vector<Double>(dim_);
    Vector<Double> resNew = Vector<Double>(dim_);
    Vector<Double> xNew = Vector<Double>(dim_);
    Vector<Double> hystNew = Vector<Double>(dim_);
    
    jacT.Mult(res,jacTres_neg);
    jacTres_neg = jacTres_neg*(-1.0);
       
		/*
		 * new: success is now an integer
		 * <0 no success
		 * >0 success but alpha was too large
		 * =0 success alpha good
		 */
    Integer success = -1;
    bool discard = false;
    
    xNew.Init();
    xNew.Add(xVal);
    Double alphaMinLocal = alphaMin;
    Double detMatToInvert;
    UInt cnt = 0;
            
    LOG_DBG(vecpreisachlinesearch) << "Determine minimal alpha that allows inversion of jacTjac";
    LOG_DBG(vecpreisachlinesearch) << "Given alphaMin: " << alphaMin;        
            
    while(true){
      matToInvert = jacTjac;
      jacTjac.Determinant(detMatToInvert);
      LOG_DBG(vecpreisachlinesearch) << " det(jacTjac) =  " << detMatToInvert;
      
      for(UInt i = 0; i < dim_; i++){
        matToInvert[i][i] += alphaMinLocal*alphaMinLocal;
      }
      matToInvert.Determinant(detMatToInvert);
      if(detMatToInvert != 0){
        break;
      } else {
        alphaMinLocal = alphaMinLocal*2;
      }
      cnt++;
      if(cnt > 200){
        EXCEPTION("LM, Linesearch: Cannot find alpha, such that jacTjac becomes invertible!");
      }
    }
    LOG_DBG(vecpreisachlinesearch) << "Found alphaMinLocal: " << alphaMinLocal;
    
    if(alpha < alphaMinLocal){
      LOG_DBG(vecpreisachlinesearch) << "Starting alpha < alphaMinLocal";
      alpha = alphaMinLocal;
    }
    
//    // now reset alpha to alphaMax (note: alphaMax is passed as copy, i.e. if we increase it further
//    // below, the actual value of alphaMax is restored during next call
//    if(alpha > alphaMax){
//      LOG_DBG(vecpreisachlinesearch) << "Current alpha > alphaMax";
//      alpha = alphaMax;
//    }

    // start at alphamin
    alpha = alphaMinLocal;
    
    while(true){
      itCnt++;
			// for statistics
			numberOfIterations=itCnt;
      
      LOG_DBG(vecpreisachlinesearch) << "Start Iteration " << itCnt << " with: ";
      LOG_DBG(vecpreisachlinesearch) << "xVal = " << xVal.ToString();
      LOG_DBG(vecpreisachlinesearch) << "last found x = " << xNew.ToString();
      LOG_DBG(vecpreisachlinesearch) << "alpha = " << alpha;
      
      for(UInt i = 0; i < dim_; i++){
        matToInvert[i][i] += alpha*alpha;
      }
      matToInvert.Determinant(detMatToInvert);
      assert(detMatToInvert != 0);
      matToInvert.Invert(matInverted);
      
      if(INV_useTikhonov_){
        Vector<Double> tikhonovReg = Vector<Double>(dim_);
        tikhonovReg.Init();
        tikhonovReg.Add(alpha*alpha,xStart,-alpha*alpha,xNew);
        tikhonovReg.Add(1.0,jacTres_neg);
        matInverted.Mult(tikhonovReg,xUpdate);
      } else {
        matInverted.Mult(jacTres_neg,xUpdate);
      }
      assert(!xUpdate.ContainsNaN() && !xUpdate.ContainsInf());
      LOG_DBG(vecpreisachlinesearch) << "Computed xUpdate = " << xUpdate.ToString();
      
      xNew.Init();
      xNew.Add(1.0,xVal,1.0,xUpdate);
      LOG_DBG(vecpreisachlinesearch) << "Computed xNew = " << xNew.ToString();
      LOG_DBG(vecpreisachlinesearch) << "Actual solution = " << sol.ToString();
      
      /*
       * New treatments
       * a) check if xUpdate.Norm > XSaturated
       *    > this will probably be an unusable update
       *    > increase alpha, go to next iteration
       *    > increase alphaMax if necessary
       * b) check if xNew = xOld + xUpdate stays above or below saturation
       *    > jumping between these two regimes leads to issues as slope
       *      of hyst operator is jumping here
       *    > if we move to another regime, xUpdate was too large
       *    > increase alpha, go to next iteration
       *    > increase alphaMax if necessary
       * c) check increment according to trust region method
       *    > if returned alpha > alphaMax or < alphaMin > stop
       *    > if check says increment is ok > stop
       *    > else go to next iteration
       * 
       * Problems:
       *  1. the larger alpha gets, the smaller the updates will be up to the
       *      point where we do not improve at all
       *  2. if alpha is very large, we will not come back in a reasonable amount
       *      of iterations
       *  3. the larger alpha gets, the more the update will have the same direction
       *      as the residual (as matrix gets more and more diagonal); this point
       *      is especially critical as rotations of the solution vector (e.g. from
       *      ex to ey) are highly depenedent on the non-diagonal entries
       * 
       * Conclusion:
       *  1. check trust region first; if trust regions says update is ok use it!
       *  2. if updated vector is out of bounds (i.e. going into sat althoug we
       *      want to stay below), we cut down the NEW vector, not the update;
       *      by doing so, we hopefully can retain the change in direction that
       *      is caused by the update
       *  
       */
            
      int successFlag = 0;
      bool debugOut = false;
      Vector<Double> hystSol = computeValue_vec(sol, operatorIdx, overwriteMemory, overwriteDirection, debugOut, successFlag);
      Vector<Double> hystOld = computeValue_vec(xVal, operatorIdx, overwriteMemory, overwriteDirection, debugOut, successFlag);
      Vector<Double> resSol = computeResidual(sol,yVal,hystSol,mu,mu_inv,wrtX,relative);
      
      hystNew = computeValue_vec(xNew, operatorIdx, overwriteMemory, overwriteDirection, debugOut, successFlag);
      resNew = computeResidual(xNew,yVal,hystNew,mu,mu_inv,wrtX,relative);
      
      LOG_DBG(vecpreisachlinesearch) << "hyst vector for sol: " << hystSol.ToString();
      LOG_DBG(vecpreisachlinesearch) << "hystNew: " << hystNew.ToString();
      LOG_DBG(vecpreisachlinesearch) << "hystOld: " << hystOld.ToString();
      LOG_DBG(vecpreisachlinesearch) << "Old res vector: " << res.ToString();
      LOG_DBG(vecpreisachlinesearch) << "New res vector: " << resNew.ToString();
      LOG_DBG(vecpreisachlinesearch) << "res vector for sol: " << resSol.ToString();
      // set stayBelowSat flag to 0 to disable checking
      success = checkIncrement(xNew, xUpdate, res, resNew, jac, alpha,0);
      LOG_DBG(vecpreisachlinesearch) << "Check trust region - return code: " << success; 
      
      if(success >= 0){
        LOG_DBG(vecpreisachlinesearch) << "Update ok according to trust region check; now check if resulting vector may be out of bounds";
        bool wasCut = false;
        if((xNew.NormL2() > XSaturated_)&&(stayBelowSat==1)){
          // too large update
          // > set xNew back to saturation
          LOG_DBG(vecpreisachlinesearch) << "2: xNew.NormL2() > Saturation but stayBelowSat = +1 > cut down xNew";
          
          // restrict; take amplitude of Y into account
          Double factor = 0.99999/PSaturated_*yVal.NormL2();
          //0.99999
          xNew.ScalarMult(factor*XSaturated_/xNew.NormL2());
          wasCut = true;
        }
        
        if((xNew.NormL2() < XSaturated_)&&(stayBelowSat==-1)){
          // too large update
          // > set xNew back to saturation
          LOG_DBG(vecpreisachlinesearch) << "3: xNew.NormL2() < Saturation but stayBelowSat = -1 (i.e. stayAboveSat = true) > cut down xNew";
          
          Double factor = 1.0/0.99999/PSaturated_*yVal.NormL2();
          //1.0/0.99999
          xNew.ScalarMult(factor*XSaturated_/xNew.NormL2());
          wasCut = true;
        }
        
        if(wasCut){
          
          // retrieve updated
          xUpdate = xNew;
          xUpdate -= xVal;
          LOG_DBG(vecpreisachlinesearch) << "Cut down xUpdate = " << xUpdate.ToString();
          // check again
          hystNew = computeValue_vec(xNew, operatorIdx, overwriteMemory, overwriteDirection, debugOut, successFlag);
          resNew = computeResidual(xNew,yVal,hystNew,mu,mu_inv,wrtX,relative);
          success = checkIncrement(xNew, xUpdate, res, resNew, jac, alpha,stayBelowSat);
          LOG_DBG(vecpreisachlinesearch) << "xNew (after cut): " << xNew.ToString();
          LOG_DBG(vecpreisachlinesearch) << "hystNew (after cut): " << hystNew.ToString();
          LOG_DBG(vecpreisachlinesearch) << "New res vector (after cut): " << resNew.ToString();
          LOG_DBG(vecpreisachlinesearch) << "Check trust region - return code (after cut): " << success; 
          if(success >= 0){
            LOG_DBG(vecpreisachlinesearch) << "Update still ok > take it";
            break;
          } else {
            LOG_DBG(vecpreisachlinesearch) << "Update no longer ok > go to next iteration except alpha out of limits";
            // any further
            if(alpha > alphaMax){
              LOG_TRACE(vecpreisachlinesearch) << "Alpha max reached; take update > stop";
              alpha = alphaMax;
              discard = false;
              break;
            }
            if(alpha < alphaMinLocal){
              LOG_TRACE(vecpreisachlinesearch) << "Alpha min reached; discard update > stop";
              alpha = alphaMinLocal;
              discard = true;
              break;
            }
            if(itCnt >= maxIter){
              LOG_TRACE(vecpreisachlinesearch) << "Linesearch was not successful; Discard update.";
              discard = true;
              break;
            }
          }
        } else {
          LOG_DBG(vecpreisachlinesearch) << "Update not out of bounds > take it";
          break;
        }
      } else {
				if(alpha > alphaMax){
					LOG_TRACE(vecpreisachlinesearch) << "Alpha max reached; take update > stop";
					alpha = alphaMax;
					discard = true;
					break;
				}
				if(alpha < alphaMinLocal){
					LOG_TRACE(vecpreisachlinesearch) << "Alpha min reached; discard update > stop";
					alpha = alphaMinLocal;
					discard = true;
					break;
				}
				
        if(itCnt >= maxIter){
          LOG_TRACE(vecpreisachlinesearch) << "Linesearch was not successful; Discard update.";
          discard = true;
          break;
        }
      }

      // might be needed again
//      if((xUpdate.NormL2() >= XSaturated_)){//&&(stayBelowSat==true)){
//        // too large update
//        // > set xNew back to saturation
//        // use only 1/2 of it
//        // > else oscilation between two near saturation states possible
//        LOG_DBG(vecpreisachlinesearch) << "1: xUpdate.NormL2() > Saturation > no sucess; increase alpha";
//        alpha = alpha*2;
//        
//        if(alpha > alphaMax){
//          LOG_DBG(vecpreisachlinesearch) << "1: alpha > alphaMax; increase alphaMax, too";
//          alphaMax = alpha;
//        }
//        continue;
//      } 
    }

    return discard;
  }
//  Initial; 17.05.2018
//  bool Hysteresis::performLinesearch(Vector<Double>& xVal, Vector<Double>& yVal, Vector<Double>& res, 
//		Vector<Double>& xUpdate, Matrix<Double>& jac, Matrix<Double>& jacT, Matrix<Double> mu, Matrix<Double> mu_inv, 
//		Integer operatorIdx, bool overwriteMemory, bool overwriteDirection,
//		Double& alpha, Double alphaMin, Double alphaMax, bool wrtX, bool relative, 
//    UInt& numberOfIterations,Vector<Double>& xStart, Double factorToSat, int stayBelowSat, Vector<Double> sol){
//    
//    LOG_TRACE(vecpreisachlinesearch) << " --------- START LINESEARCH --------- ";
//    LOG_DBG(vecpreisachlinesearch) << "Starting xVal = " << xVal.ToString();
//    LOG_DBG(vecpreisachlinesearch) << "Target yVal = " << yVal.ToString();
//    LOG_DBG(vecpreisachlinesearch) << "alpha, alphaMin, alphaMax = " << alpha << ", " << alphaMin << ", " << alphaMax << std::endl;
//    
//    if((xVal.NormL2() > XSaturated_)&&(stayBelowSat==1)){
//      EXCEPTION("xInput to Linesearch already above saturation > must not be the case!");
//    }
//    if((xVal.NormL2() < XSaturated_)&&(stayBelowSat==-1)){
//      EXCEPTION("xInput to Linesearch already below saturation > must not be the case!");
//    }
//
//		UInt maxIter = 50;
//    UInt itCnt = 0;
//    
//    Matrix<Double> matToInvert = Matrix<Double>(dim_,dim_);
//    Matrix<Double> matInverted = Matrix<Double>(dim_,dim_);
//    Matrix<Double> jacTjac = Matrix<Double>(dim_,dim_);
//    jacT.Mult(jac,jacTjac);
//    Vector<Double> jacTres_neg = Vector<Double>(dim_);
//    Vector<Double> resNew = Vector<Double>(dim_);
//    Vector<Double> xNew = Vector<Double>(dim_);
//    Vector<Double> hystNew = Vector<Double>(dim_);
//    
//    jacT.Mult(res,jacTres_neg);
//    jacTres_neg = jacTres_neg*(-1.0);
//       
//		/*
//		 * new: success is now an integer
//		 * <0 no success
//		 * >0 success but alpha was too large
//		 * =0 success alpha good
//		 */
//    Integer success = -1;
//    bool discard = false;
//    
//    xNew.Init();
//    xNew.Add(xVal);
//    Double alphaMinLocal = alphaMin;
//    Double detMatToInvert;
//    UInt cnt = 0;
//            
//    LOG_DBG(vecpreisachlinesearch) << "Determine minimal alpha that allows inversion of jacTjac";
//    LOG_DBG(vecpreisachlinesearch) << "Given alphaMin: " << alphaMin;        
//            
//    while(true){
//      matToInvert = jacTjac;
//      jacTjac.Determinant(detMatToInvert);
//      LOG_DBG(vecpreisachlinesearch) << " det(jacTjac) =  " << detMatToInvert;
//      
//      for(UInt i = 0; i < dim_; i++){
//        matToInvert[i][i] += alphaMinLocal*alphaMinLocal;
//      }
//      matToInvert.Determinant(detMatToInvert);
//      if(detMatToInvert != 0){
//        break;
//      } else {
//        alphaMinLocal = alphaMinLocal*2;
//      }
//      cnt++;
//      if(cnt > 200){
//        EXCEPTION("LM, Linesearch: Cannot find alpha, such that jacTjac becomes invertible!");
//      }
//    }
//    LOG_DBG(vecpreisachlinesearch) << "Found alphaMinLocal: " << alphaMinLocal;
//    
//    if(alpha < alphaMinLocal){
//      LOG_DBG(vecpreisachlinesearch) << "Starting alpha < alphaMinLocal";
//      alpha = alphaMinLocal;
//    }
//    
//    // now reset alpha to alphaMax (note: alphaMax is passed as copy, i.e. if we increase it further
//    // below, the actual value of alphaMax is restored during next call
//    if(alpha > alphaMax){
//      LOG_DBG(vecpreisachlinesearch) << "Current alpha > alphaMax";
//      alpha = alphaMax;
//    }
//
//    while(true){
//      itCnt++;
//			// for statistics
//			numberOfIterations=itCnt;
//      
//      LOG_DBG(vecpreisachlinesearch) << "Start Iteration " << itCnt << " with: ";
//      LOG_DBG(vecpreisachlinesearch) << "xVal = " << xVal.ToString();
//      LOG_DBG(vecpreisachlinesearch) << "last found x = " << xNew.ToString();
//      LOG_DBG(vecpreisachlinesearch) << "alpha = " << alpha;
//      
//      for(UInt i = 0; i < dim_; i++){
//        matToInvert[i][i] += alpha*alpha;
//      }
//      matToInvert.Determinant(detMatToInvert);
//      assert(detMatToInvert != 0);
//      matToInvert.Invert(matInverted);
//      
//      if(INV_useTikhonov_){
//        Vector<Double> tikhonovReg = Vector<Double>(dim_);
//        tikhonovReg.Init();
//        tikhonovReg.Add(alpha*alpha,xStart,-alpha*alpha,xNew);
//        tikhonovReg.Add(1.0,jacTres_neg);
//        matInverted.Mult(tikhonovReg,xUpdate);
//      } else {
//        matInverted.Mult(jacTres_neg,xUpdate);
//      }
//      assert(!xUpdate.ContainsNaN() && !xUpdate.ContainsInf());
//      LOG_DBG(vecpreisachlinesearch) << "Computed xUpdate = " << xUpdate.ToString();
//      
//      /*
//       * New treatments
//       * a) check if xUpdate.Norm > XSaturated
//       *    > this will probably be an unusable update
//       *    > increase alpha, go to next iteration
//       *    > increase alphaMax if necessary
//       * b) check if xNew = xOld + xUpdate stays above or below saturation
//       *    > jumping between these two regimes leads to issues as slope
//       *      of hyst operator is jumping here
//       *    > if we move to another regime, xUpdate was too large
//       *    > increase alpha, go to next iteration
//       *    > increase alphaMax if necessary
//       * c) check increment according to trust region method
//       *    > if returned alpha > alphaMax or < alphaMin > stop
//       *    > if check says increment is ok > stop
//       *    > else go to next iteration
//       */
//       
//      if((xUpdate.NormL2() >= XSaturated_)){//&&(stayBelowSat==true)){
//        // too large update
//        // > set xNew back to saturation
//        // use only 1/2 of it
//        // > else oscilation between two near saturation states possible
//        LOG_DBG(vecpreisachlinesearch) << "1: xUpdate.NormL2() > Saturation > no sucess; increase alpha";
//        alpha = alpha*2;
//        
//        if(alpha > alphaMax){
//          LOG_DBG(vecpreisachlinesearch) << "1: alpha > alphaMax; increase alphaMax, too";
//          alphaMax = alpha;
//        }
//        continue;
//      } 
//
//      xNew.Init();
//      xNew.Add(1.0,xVal,1.0,xUpdate);
//      LOG_DBG(vecpreisachlinesearch) << "Computed xNew = " << xNew.ToString();
//      LOG_DBG(vecpreisachlinesearch) << "Actual solution = " << sol.ToString();
//      
//      if((xNew.NormL2() > XSaturated_)&&(stayBelowSat==1)){
//        // too large update
//        // > set xNew back to saturation
//        LOG_DBG(vecpreisachlinesearch) << "2: xNew.NormL2() > Saturation but stayBelowSat = true > no sucess; increase alpha";
//        alpha = alpha*2;
//        
//        if(alpha > alphaMax){
//          LOG_DBG(vecpreisachlinesearch) << "2: alpha > alphaMax; increase alphaMax, too";
//          alphaMax = alpha;
//        }
//        continue;
//      }
//      
//      if((xNew.NormL2() < XSaturated_)&&(stayBelowSat==-1)){
//        // too large update
//        // > set xNew back to saturation
//        LOG_DBG(vecpreisachlinesearch) << "3: xNew.NormL2() < Saturation but stayBelowSat = false > no sucess; increase alpha";
//        alpha = alpha*2;
//        
//        if(alpha > alphaMax){
//          LOG_DBG(vecpreisachlinesearch) << "3: alpha > alphaMax; increase alphaMax, too";
//          alphaMax = alpha;
//        }
//        continue;
//      }
//      
//      Vector<Double> hystSol = computeValue_vec(sol, operatorIdx, overwriteMemory, overwriteDirection);
//      Vector<Double> hystOld = computeValue_vec(xVal, operatorIdx, overwriteMemory, overwriteDirection);
//      Vector<Double> resSol = computeResidual(sol,yVal,hystSol,mu,mu_inv,wrtX,relative);
//      
//      hystNew = computeValue_vec(xNew, operatorIdx, overwriteMemory, overwriteDirection);
//      resNew = computeResidual(xNew,yVal,hystNew,mu,mu_inv,wrtX,relative);
//      
//      LOG_DBG(vecpreisachlinesearch) << "hyst vector for sol: " << hystSol.ToString();
//      LOG_DBG(vecpreisachlinesearch) << "hystNew: " << hystNew.ToString();
//      LOG_DBG(vecpreisachlinesearch) << "hystOld: " << hystOld.ToString();
//      LOG_DBG(vecpreisachlinesearch) << "Old res vector: " << res.ToString();
//      LOG_DBG(vecpreisachlinesearch) << "New res vector: " << resNew.ToString();
//      LOG_DBG(vecpreisachlinesearch) << "res vector for sol: " << resSol.ToString();
//      success = checkIncrement(xNew, xUpdate, res, resNew, jac, alpha,stayBelowSat);
//      LOG_DBG(vecpreisachlinesearch) << "Check trust region - return code: " << success;
//      
//      // after trust region checking, we do not increase alphaMax or diecrease alphaMin
//      // any further
//      if(alpha > alphaMax){
//        alpha = alphaMax;
//        break;
//      }
//      if(alpha < alphaMinLocal){
//        alpha = alphaMinLocal;
//        break;
//      }
//	
//      if(success >= 0){
//        LOG_TRACE(vecpreisachlinesearch) << "Linesearch was successful after " << itCnt << " iterations";
//        break;
//      } else {
//        if(itCnt >= maxIter){
//          LOG_TRACE(vecpreisachlinesearch) << "Linesearch was not successful; Discard update.";
//          discard = true;
//          break;
//        }
//      }
//    }
//
//    return discard;
//  }
  
  
  
  
//  
//  bool Hysteresis::performLinesearch(Vector<Double>& xVal, Vector<Double>& yVal, Vector<Double>& res, 
//		Vector<Double>& xUpdate, Matrix<Double>& jac, Matrix<Double>& jacT, Matrix<Double> mu, Matrix<Double> mu_inv, 
//		Integer operatorIdx, bool overwriteMemory, bool overwriteDirection,
//		Double& alpha, Double alphaMin, Double alphaMax, bool wrtX, bool relative, 
//    UInt& numberOfIterations,Vector<Double>& xStart, Double factorToSat, bool stayBelowSat){
//    
//    LOG_DBG(vecpreisachInversion) << " --------- START LINESEARCH --------- ";
//    if((xVal.NormL2() >= XSaturated_)&&(stayBelowSat==true)){
//      EXCEPTION("xInput to Linesearch already above saturation > must not be the case!");
//    }
//    if((xVal.NormL2() < XSaturated_)&&(stayBelowSat==false)){
//      EXCEPTION("xInput to Linesearch already below saturation > must not be the case!");
//    }
//    
//    LOG_TRACE(vecpreisachInversion) << "Start Linesearch with xVal = " << xVal.ToString();
//    
//		//LOG_DBG(vecpreisach) << "Old alpha " << alpha;
//		
//		UInt maxIter = 25;
//    UInt itCnt = 0;
//    
//    Matrix<Double> matToInvert = Matrix<Double>(dim_,dim_);
//    Matrix<Double> matInverted = Matrix<Double>(dim_,dim_);
//    Matrix<Double> jacTjac = Matrix<Double>(dim_,dim_);
//    jacT.Mult(jac,jacTjac);
//    Vector<Double> jacTres_neg = Vector<Double>(dim_);
//    Vector<Double> resNew = Vector<Double>(dim_);
//    Vector<Double> xNew = Vector<Double>(dim_);
//    Vector<Double> hystNew = Vector<Double>(dim_);
//    
//    jacT.Mult(res,jacTres_neg);
//    jacTres_neg = jacTres_neg*(-1.0);
//    
////    Double alphaMax = 256.0; //1e0;//1e1;//e10;
////    // too small not working, too large not working either e-18,e-12,e-4 > no; e-14,e-15 ok
////    Double alphaMin = 1.0/1024; //1e-10;//15; 
//    
//		/*
//		 * new: success is now an integer
//		 * <0 no success
//		 * >0 success but alpha was too large
//		 * =0 success alpha good
//		 */
//    Integer success = -1;
//    bool discard = false;
//    
//    while(true){
//      LOG_TRACE(vecpreisachInversion) << "Start Iteration with xVal = " << xVal.ToString();
//      itCnt++;
//			// for statistics
//			numberOfIterations=itCnt;
//
//      Double detMatToInvert;
//      UInt cnt = 0;
//      while(true){
//        matToInvert = jacTjac;
//        LOG_DBG(vecpreisachInversion) << " jacTjac =  " << jacTjac.ToString();
//        jacTjac.Determinant(detMatToInvert);
//        LOG_DBG(vecpreisachInversion) << " det(jacTjac) =  " << detMatToInvert;
//        
//        for(UInt i = 0; i < dim_; i++){
//          matToInvert[i][i] += alpha*alpha;
//        }
//        matToInvert.Determinant(detMatToInvert);
//        if(detMatToInvert != 0){
//          break;
//        } else {
//          alpha = alpha*2;
//        }
//        cnt++;
//        if(cnt > 200){
//          EXCEPTION("Cannot get matrix for inversion that is invertible!");
//        }
//      }
//      
//      LOG_TRACE(vecpreisachInversion) << " alpha =  " << alpha;
//      LOG_DBG(vecpreisachInversion) << " matToInvert =  " << matToInvert.ToString();
//      matToInvert.Determinant(detMatToInvert);
//      LOG_DBG(vecpreisachInversion) << " det(matToInvert) =  " << detMatToInvert;
//      matToInvert.Invert(matInverted);
//      LOG_DBG(vecpreisachInversion) << " matInverted =  " << matInverted.ToString();
//      
//      if(INV_useTikhonov_){
//        Vector<Double> tikhonovReg = Vector<Double>(dim_);
//        tikhonovReg.Init();
//        tikhonovReg.Add(alpha*alpha,xStart,-alpha*alpha,xVal);
//        tikhonovReg.Add(1.0,jacTres_neg);
//        matInverted.Mult(tikhonovReg,xUpdate);
//      } else {
//        matInverted.Mult(jacTres_neg,xUpdate);
//      }
//      assert(!xUpdate.ContainsNaN() && !xUpdate.ContainsInf());
//      xNew.Init();
//      
//      if((xUpdate.NormL2() >= XSaturated_)){//&&(stayBelowSat==true)){
//        // too large update
//        // > set xNew back to saturation
//        // use only 1/2 of it
//        // > else oscilation between two near saturation states possible
//        LOG_TRACE(vecpreisachInversion) << "xUpdate > Saturation > scale down";
//        LOG_DBG(vecpreisachInversion) <<  "xUpdate = " << xUpdate.ToString();
//        xUpdate.ScalarDiv(2.0);
//        LOG_DBG(vecpreisachInversion) <<  "xUpdate (after scaling down) = " << xUpdate.ToString();
//        alpha = alpha*16;
//      } 
//      LOG_TRACE(vecpreisachInversion) << "xUpdate = " << xUpdate.ToString();
//
//      xNew.Add(1.0,xVal,1.0,xUpdate);
//      LOG_TRACE(vecpreisachInversion) << "xNew = " << xNew.ToString();
//      if((xNew.NormL2() >= XSaturated_)&&(stayBelowSat==true)){
//        // too large update
//        // > set xNew back to saturation
//        LOG_TRACE(vecpreisachInversion) << "Reset xUpdate as x above saturation";
//        // go closer and closer to saturation
//        xNew.ScalarMult(factorToSat*XSaturated_/xNew.NormL2());
//        LOG_TRACE(vecpreisachInversion) << "Set solution to " << factorToSat << " times Saturation = " << xNew.ToString();
//        xUpdate = xNew;
//        xUpdate -= xVal;
//        alpha = alpha*2.0;
//        LOG_TRACE(vecpreisachInversion) << "new alpha " << alpha;
////        if(alpha >= alphaMax){
////          LOG_TRACE(vecpreisachInversion) << "alpha > alphaMax";
////          //        LOG_DBG(vecpreisach) << "Maximal alpha reached > stop";
////          // maximal alpha used; stop here (regardless of success)
////          alphaMax*=2;
//////          break;
////        }
//      }
//      if((xNew.NormL2() < XSaturated_)&&(stayBelowSat==false)){
//        LOG_TRACE(vecpreisachInversion) << "Reset xUpdate as x below saturation";
//        // go closer and closer to saturation
//        xNew.ScalarMult(1.0/factorToSat*XSaturated_/xNew.NormL2());
//        LOG_TRACE(vecpreisachInversion) << "Set solution to " << 1.0/factorToSat << " times Saturation = " << xNew.ToString();
//        xUpdate = xNew;
//        xUpdate -= xVal;
//        alpha = alpha*2;
//        LOG_TRACE(vecpreisachInversion) << "new alpha " << alpha;
//      }
//      
//  
//      hystNew = computeValue_vec(xNew, operatorIdx, overwriteMemory, overwriteDirection);
//      resNew = computeResidual(xNew,yVal,hystNew,mu,mu_inv,wrtX,relative);
//        
//      //std::cout << "Current iteration: " << itCnt << std::endl;
//      //std::cout << "Alpha pre: " << alpha << std::endl;
//      success = checkIncrement(xNew, xUpdate, res, resNew, jac, alpha,stayBelowSat);
//      //std::cout << "Alpha post: " << alpha << std::endl;
////      LOG_DBG(vecpreisach) << "New alpha " << alpha;
////      LOG_DBG(vecpreisach) << "Success? " << success;
//      if(alpha > alphaMax){
////        LOG_DBG(vecpreisach) << "Maximal alpha reached > stop";
//        // maximal alpha used; stop here (regardless of success)
//        alpha = alphaMax;
//        break;
//      }
//      if(alpha < alphaMin){
////        LOG_DBG(vecpreisach) << "Minimal alpha reached > stop";
//        // minimal alpha used; stop here
//        alpha = alphaMin;
//        break;
//      }
//	
//      if(success >= 0){
////        LOG_DBG(vecpreisach) << "Linesearch was successful after " << itCnt << " iterations";
//        break;
//      } else {
//        if(itCnt >= maxIter){
////          LOG_DBG(vecpreisach) << "Linesearch was not successful; Discard update.";
//          discard = true;
//          break;
//        }
//      }
//    }
//
//    return discard;
//  }
  
	Vector<Double> Hysteresis::computeInput_vec_withPrevStates(Vector<Double> yVal, Vector<Double> prevYval,
      Vector<Double> prevXval, Vector<Double> prevHystval, Integer operatorIndex, 
      Matrix<Double> mu, bool overwriteDirection, bool fieldsAlignedAboveSat, 
      bool hystOutputRestrictedToSat, int& successFlag){
		
		UInt totalNumberOfLMIterations=0;
		UInt totalNumberOfLinesearchIterations=0;
		UInt maximalNumberOfLinesearchIterations=0;
    Double minAlphaStatistics,maxAlphaStatistics,avgAlphaStatistics;
		Vector<Double>sol = Vector<Double>(dim_);
    
		return computeInput_vec_withStatistics(yVal, prevYval, prevXval, prevHystval, 
			operatorIndex, mu, overwriteDirection, fieldsAlignedAboveSat, hystOutputRestrictedToSat,
      totalNumberOfLMIterations, totalNumberOfLinesearchIterations, 
      maximalNumberOfLinesearchIterations, successFlag, minAlphaStatistics,maxAlphaStatistics,avgAlphaStatistics,sol);
	}
    
  Vector<Double> Hysteresis::computeInput_vec_withStatistics(Vector<Double> yVal, Vector<Double> prevYval,
          Vector<Double> prevXval, Vector<Double> prevHystval, Integer operatorIndex, 
          Matrix<Double> mu, bool overwriteDirection, bool fieldsAlignedAboveSat, bool hystOutputRestrictedToSat,
          UInt& totalNumberOfLMIterations, UInt& totalNumberOfLinesearchIterations, 
          UInt& maximalNumberOfLinesearchIterations, int& successFlag, Double& minAlphaStatistics, 
					Double& maxAlphaStatistics, Double& avgAlphaStatistics, Vector<Double> sol){
    
    assert(!yVal.ContainsNaN() && !yVal.ContainsInf());
		/*
     * IMPORTANT: do not pass yVal, prevYval, prevXval nor prevHystval as reference
     *						> we do not want to overwrite these values by accident!
     */
    LOG_TRACE(vecpreisachInversion) << " --------- START IVERSION --------- ";
    LOG_TRACE(vecpreisachInversion) << " yVal = " << yVal.ToString(); 
    bool debug = !true;

    std::stringstream traceMsg;
    if(debug){
      traceMsg << " --- TRACE MSG --- " << std::endl;
      traceMsg << "VecPreisach::computeInput_vec for operator index: " << operatorIndex << std::endl;
      traceMsg << "Input value yVal: " << yVal.ToString() << std::endl;
      traceMsg << "Previous input value yVal_old: " << prevYval.ToString() << std::endl;
      traceMsg << "Previous hyst value hVal_old: " << prevHystval.ToString() << std::endl;
      traceMsg << "Previous retrieved input value xVal_old: " << prevXval.ToString() << std::endl;
      traceMsg << "Previous solution: Percentage of saturation (|xVal|/XSaturated): " << prevXval.NormL2()/XSaturated_ << std::endl;
      traceMsg << "Material Tensor: " << mu.ToString() << std::endl;
      traceMsg << "ySat: " << PSaturated_ << std::endl;
      traceMsg << "xSat: " << XSaturated_ << std::endl;
      traceMsg << "INV_maxIter_: " << INV_maxIter_ << std::endl;
      traceMsg << "INV_resTolH_: " << INV_resTolH_ << std::endl;
      traceMsg << "INV_resTolB_: " << INV_resTolB_ << std::endl;
      traceMsg << "INV_jacobiResolution_: " << INV_jacobiResolution_ << std::endl;
      traceMsg << "INV_useTikhonov_: " << INV_useTikhonov_ << std::endl;
      traceMsg << "INV_alphaLSStart_: " << INV_alphaLSStart_ << std::endl;     
      traceMsg << "INV_alphaLSMin_: " << INV_alphaLSMin_ << std::endl;  
      traceMsg << "INV_alphaLSMax_: " << INV_alphaLSMax_ << std::endl;  
    }
    
    // for statistics
    int successFlagForward = 0;
    successFlag = -1; 
    // -1 = fail
    //  0 = reuse value
    //  1 = anhyst only
    //  2 = bisection
    //  3-6 only for vector implementation using Levenberg Marquardt
    //  3 = reamnence
    //  4 = passed dut to error tolerance
    //  5 = passed due to tolerance wrt x
    //  6 = passed due to tolerance wrt y
    
    totalNumberOfLMIterations = 0;
		totalNumberOfLinesearchIterations = 0;
		maximalNumberOfLinesearchIterations = 0;
    
    minAlphaStatistics = 1e16;
    maxAlphaStatistics = -1e16;
    avgAlphaStatistics = 0;
    
		// during inversion, never overwrite hyst-memory
		bool overwriteMemory = false;
    
    // return value
    Vector<Double> xVal = Vector<Double>(dim_);
    
    /*
     * Check if yVal is significantly different from previous value
     */
    Vector<Double> diff;
		diff = yVal;
    diff -= prevYval;
    
    if(diff.NormL2() < INV_resTolB_){
      traceMsg << "--A-- Inversion: Reuse old value" << std::endl;
      xVal = prevXval;
      LOG_TRACE(vecpreisachInversion) << "Reused value xVal: " << xVal.ToString();
      successFlag = 0;
      return xVal;
    }
		
    /*
     * Invert eps/mu for later usage
     */
    Matrix<Double> mu_inv = Matrix<Double>(dim_,dim_);
    Vector<Double> xTMP = Vector<Double>(dim_);
    mu.Invert(mu_inv);
    Double yNorm = yVal.NormL2();
    
    Vector<Double> yDir = Vector<Double>(dim_);
    yDir.Init();
    if(yNorm > 0){
      yDir = yVal;
      yDir.ScalarDiv(yNorm);
    } else {
      yDir[0] = 1.0; // default > take x-direction
      // only required to compute eps_mu
    }
    // get mu in current direction
    Vector<Double> tmp = Vector<Double>(dim_);
    Double eps_mu = 0.0;
    mu.Mult(yDir,tmp);
    yDir.Inner(tmp,eps_mu);

    if(anhystOnly_ == true){
      traceMsg << "--S-- Inversion: Special case, only anhysteretic part > solve by bisection" << std::endl;
      successFlag = 1;
      xVal.Init();
      if(yNorm == 0){
        // anhyst part has no remanence
        return xVal;
      }
      // 
      Double tol = 1e-12;
      Double Xup, Xdown, Poffset, xScal;
      Xup = yNorm/eps_mu;
      Xdown = 0;
      Poffset = 0;
      xScal = bisectForAnhyst(yNorm, Xdown, Xup, Poffset, eps_mu, tol); 
      xVal.Add(xScal,yDir);
      
      return xVal;
    }    
    
    /*
     * Check if yVal is beyond saturation > easy case
     */
    /*
     * NOTE: when we are in saturation (or if we assume to be for checking)
     * we have y and x aligned, i.e. we can check for this case the same way
     * as we do for scalar problems
     * However, we only have to check for positive saturation as we check the
     * amplitude which always should be >= 0
     */
    /*
     * 14.05.2018
     * NEW addition due to Mayergoyz hyst model - evalAboveSaturation
     * Background: 
     *  Mayergoyz Vector Hyst model does not always return YSat*inputDir
     *  if X = XSat*inputDir; this is especially the case, if remament components
     *  perpendicular to inputDir are present in the history of the hyst operator;
     *  unlike Sutors Vector model, these remanent parts vanish only if X is
     *  going towards infinity (Sutor model: remanent parts get 0 for input = sat)
     * Consequence:
     *  If y-input is larger than saturation, we cannot simply state that
     *  output of hyst operator (excluding anhysteretic terms) is Ysat*inputDir;
     *  nevertheless, we can say that X will be aligned with y-input, i.e.
     *  X = Xampl * inputDir; this still allows the reduction to a 1d linesearch
     *  along inputDir; the only difference is, that we have to evaluate the 
     *  hyst operator instead of setting its amplitude to Ysat 
     */
        
    // stayBelowSat
    //  1: stay below 
    //  0: undefined
    // -1: stay above
    int stayBelowSat = 0;
    
    LOG_TRACE(vecpreisachInversion) << "yNorm: " << yNorm;
    if(yNorm == 0){
      stayBelowSat = 1;
    } else {
         
      Double anhystPartPosSat = PSaturated_*evalAnhystPart_normalized(1.0);
      LOG_TRACE(vecpreisachInversion) << "anhystPartPosSat: " << anhystPartPosSat;
      LOG_TRACE(vecpreisachInversion) << "(PSaturated_ + eps_mu*XSaturated_ + anhystPartPosSat): " << (PSaturated_ + eps_mu*XSaturated_ + anhystPartPosSat);
      LOG_TRACE(vecpreisachInversion) << "yNorm - (PSaturated_ + eps_mu*XSaturated_ + anhystPartPosSat): " << yNorm - (PSaturated_ + eps_mu*XSaturated_ + anhystPartPosSat);
      
      
      // THESE CHECKS are not valid for Mayergoyz model as XDir != YDir
      // > this is also the reason why bisection is not possible > we simply do not know the actual direction of the output!
      
      // TODO:
      // make use of triangle inequality
      //  if |y| > |P(Xsat)| + |eps*(Xsat)| + |anhyst(XSat)| >= |P(Xsat) + eps*(Xsat) + anhyst(XSat)|, then y is in saturation
      //    no matter what direction Xsat has; too bad, that we do not know the direction of X
      // > therefore we cannot compute |P(Xsat)| + |eps*(Xsat)| + |anhyst(XSat)|
      // but: we can compute |P(Xsat) + eps*(Xsat) + anhyst(XSat)| and check
      //  if |y| < |P(Xsat) + eps*(Xsat) + anhyst(XSat)|; if this is the case, we are definitely below saturation
      //
      // in case of sutor mode, P(Xsat) will have direction of Xsat which is the direction of y
      // in that case |P(Xsat)| + |eps*(Xsat)| + |anhyst(XSat)| = |P(Xsat) + eps*(Xsat) + anhyst(XSat)|
      // ( assuming anhystPart and eps > 0 )
      //
      // if |y| > |ySat + eps_mu*xSat + anhystPart(xSat)| (i.e. anhystPart and hystPart are aligned)
      //    x has to be above saturation  > stayBelowSat = false = -1
      // if |y| < |ySat - eps_mu*xSat - anhystPart(xSat)| (i.e. anhystPart and hystPart are antiparallel)
      //    x cannot be in saturation > stayBelowSat = true = 1
      // if |ySat + eps_mu*xSat + anhystPart(xSat)| > |y| > |ySat - eps_mu*xSat - anhystPart(xSat)|
      //    x could be above saturation or below > stayBelowSat = unknown = 0
      //
      // note: for sutor model, hystPart and anhystPart are always aligned above saturation
      //        > only case 1 relevant for checking
      //       for mayergoyz model, we have to check case 2 and 3, too
      //
      // note 22.05.2018: in case of sutor model, PSaturated might be unreachable if rotRes < 0 (for revised model)
      //                  > compute actual value at XSaturated instead of using PSaturated;
      //                  > by this, we also do not need to add anhystPartPosSat anymore, as this is done by evaluating the hyst operator
      int tmp = 0; 
      Vector<Double> satInput = Vector<Double>(dim_);
      satInput.Init();
      satInput.Add(XSaturated_,yDir);
      Vector<Double> hystValAtXSat = computeValue_vec(satInput, operatorIndex, false, true, false, tmp);
      
      // check saturation in direction of yIn might solve system
      //Double diffSat = yNorm - (PSaturated_ + eps_mu*XSaturated_ + anhystPartPosSat);
      Double diffSat = yNorm - (hystValAtXSat.NormL2() + eps_mu*XSaturated_);
      if(abs(diffSat) < INV_resTolB_){
        traceMsg << "--B Special-- Inversion: Exact Saturation found" << std::endl;
//        xVal.Init();
//        xVal.Add(XSaturated_,yDir);
        successFlag = 2;
        return satInput;
      }

      bool useBisection = false;
      //if(yNorm >= (PSaturated_ + eps_mu*XSaturated_ + anhystPartPosSat)){
      if(yNorm >= (hystValAtXSat.NormL2() + eps_mu*XSaturated_)){
        // |y| > |ySat + eps_mu*xSat + anhystPart(xSat)|
        traceMsg << "--I 1-- Inversion: Input above Saturation found" << std::endl;
        stayBelowSat = -1;
        if(fieldsAlignedAboveSat == true){
          useBisection = true;
        }
				if(hystOutputRestrictedToSat == false){
					// here we cannot actually say if input really has to be above saturation
					// as yNorm might not suffice anymore
					stayBelowSat = 0;
				}
      //} else if (yNorm <= (PSaturated_ - eps_mu*XSaturated_ - anhystPartPosSat)){
      // here we have to subtract anhystPart twice as the one inside hystValAtXSat is aligned with yDir but we want it to be antiparallel
      } else if (yNorm <= (hystValAtXSat.NormL2() - eps_mu*XSaturated_ - 2*anhystPartPosSat)){
        traceMsg << "--I 2-- Inversion: Input below Saturation found" << std::endl;
        // |y| < |ySat - eps_mu*xSat - anhystPart(xSat)|
        stayBelowSat = 1;
      } else {
        traceMsg << "--I 3-- Inversion: Unclear if input leads to X above or belwo Saturation" << std::endl;
        // |ySat + eps_mu*xSat + anhystPart(xSat)| > |y| > |ySat - eps_mu*xSat - anhystPart(xSat)|
        if(fieldsAlignedAboveSat == true){
          traceMsg << "--I 3a-- Inversion: Fields aligned; X must stay below sat" << std::endl;
          // if fields are aligned, the lower bound |ySat - eps_mu*xSat - anhystPart(xSat)| is not
          // possible, i.e. the previous state and this state merge 
          stayBelowSat = 1;
        } else {
          traceMsg << "--I 3a-- Inversion: Fields may be misaligned; X needs not stay below sat" << std::endl;
          stayBelowSat = 0;
        }
      }
      
      if(useBisection){
//        std::cout << "UseBisection" << std::endl;
        Double xScal = 0.0;
        
        if(anhystPartPosSat == 0){
          traceMsg << "--B1-- Inversion: Anhysteretic part zero > solve by simple division" << std::endl;
          //xScal = (yNorm - PSaturated_)/eps_mu;
          xScal = (yNorm - hystValAtXSat.NormL2())/eps_mu;
          xVal.Init();
          xVal.Add(xScal,yDir);
          successFlag = 2;
          return xVal;
        } else {
          traceMsg << "--B2-- Inversion: Anhysteretic non-zero > solve by bisection" << std::endl;
          Double tol = 1e-12;
          Double Xup, Xdown, Poffset;
//          Xup = (yNorm - PSaturated_)/eps_mu;
//          Xdown = XSaturated_;
//          Poffset = PSaturated_;
          Xup = (yNorm - hystValAtXSat.NormL2())/eps_mu;
          Xdown = hystValAtXSat.NormL2();
          Poffset = hystValAtXSat.NormL2();
          xScal = bisectForAnhyst(yNorm, Xdown, Xup, Poffset, eps_mu, tol);
          xVal.Init();
          xVal.Add(xScal,yDir);
          successFlag = 2;
          return xVal;
        }
      }
    } // yNorm == 0?
        
    /*
     * Use Levenberg-Marquart algorithm as presented in Dahmen&Reusken - Numerik partialler DFG
     */   
     // LM
    Vector<Double> hystVal;
    Vector<Double> hystVal_rem;
    /*
     * Check for remanence
     * > apparently, LM seems to work fine except for the remanence case
     * i.e. yval = hystval, xval = 0
     * > check if xVal = 0 would be a proper solution
     */
    xTMP.Init();
    bool debugOut = false;
    hystVal_rem = computeValue_vec(xTMP, operatorIndex, overwriteMemory, overwriteDirection, debugOut, successFlagForward);
    
    diff.Init();
    diff.Add(1.0,yVal,-1.0,hystVal_rem);
    
    LOG_DBG(vecpreisachInversion) << "Check difference between yVal and hystVal(0); remanence?";
    LOG_DBG(vecpreisachInversion) << "Diff Vector: " << diff.ToString();
    LOG_DBG(vecpreisachInversion) << "Norm of diff: " << diff.NormL2();
    
    if(diff.NormL2() < INV_resTolB_){
      traceMsg << "--C-- Inversion: Remanence detected" << std::endl;
      xVal = xTMP;
      LOG_DBG(vecpreisachInversion) << "Set xVal to 0: " << xVal.ToString();
      successFlag = 3;
      //	compCase = 2;
    } else {
      traceMsg << "--D-- Inversion: Use Levenberg Marquart" << std::endl;
      
      // tolerance wrt y > 1e-10 or 1e-12 seems good > takes 2-3 its
      // only problem: y-x-loops look ugly as x can be quite off!
      //Double tolError = 1e-11;  
      // tolerance for reevalution
      
      UInt itCnt = 0;
      Double alpha = INV_alphaLSStart_;
      Double alphaMin = INV_alphaLSMin_;//1.0/256.0;
      Double alphaMax = INV_alphaLSMax_;//8192;//512.0;
      			
      /*
       * only residual wrt X in abs form tested
       */
      bool wrtX = true;
      bool relError = !true;
      
      if (diff.NormL2() < INV_resTolB_){
        // we are quite close to remanence, but have not yet reached it
        // start from 0
        traceMsg << "--D-- Inversion: Start at 0" << std::endl;
        xVal = xTMP;
        //  startAtZero=true;
        //	compCase = 3;
      } else {
        //	compCase = 4;
        xVal = prevXval;
        traceMsg << "--D-- Inversion: Start at previous X value = " << xVal.ToString() << std::endl;
        // different approach:
        // to this point we know:
        // old and new Y value
        // old X value
        // old H value
        // H value at 0 (H_rem)
        // > by this we can check:
        // a) newY > oldY > X should increase, too
        // b) newY < oldY < X should decrease, too
        //      b1) H_rem > newY > X = 0 is still too large but maybe good starting point
        //      b2) H_rem < newY > X between old X value and 0
        // > note: the steps above (sort of bisect) would work well for scalar model; for
        //         vector model we might have a change in direction but at similiar amplitude
        //         such that we do not match aboves cases
        //         still, this ssems more reasonable than fp iteration that nearly always drives
        //         x way above saturation due to large mu_inv
        
        
        traceMsg << "Try to obtain better starting value" << std::endl;
        if(yVal.NormL2() > prevYval.NormL2()){
          traceMsg << "Y increased in norm > increase X ,too" << std::endl;
          // increase X a bit (does not work if X is zero, though)
          xVal = prevXval;
          
          if(xVal.NormL2() != 0){
            Double factor;
            if(prevYval.NormL2() != 0){
              factor = yVal.NormL2()/prevYval.NormL2();
            } else {
              factor = 1.05;
            }
            traceMsg << "Scale X by " << factor << std::endl;
            xVal.ScalarMult(factor);
          } else {
            traceMsg << "X is zero; scaling will not help; try a scaled version of yVal instead" << std::endl;
            xVal = yVal;
            xVal.ScalarMult(xVal.NormL2()*XSaturated_/PSaturated_);
          }
        } else {
          traceMsg << "Y decreased in norm > decrease X ,too" << std::endl;
          //            if(hystVal_rem.NormL2() > yVal.NormL2()){
          //              LOG_DBG(vecpreisach) << "Y smaller than remanence (norm-wise) > start at remanence";
          //              // start at 0
          //              xVal.Init();
          //            } else {
          //              LOG_DBG(vecpreisach) << "Y larger than remanence but smaller than current hyst-value (norm-wise) > start in between";
          // take midpoint between remanence and currenc xval
          xVal = prevXval;
          
          if(xVal.NormL2() != 0){
            Double factor;
            if(prevYval.NormL2() != 0){
              factor = yVal.NormL2()/prevYval.NormL2();
            } else {
              factor = 1.0/1.05;
            }
            traceMsg << "Scale X by " << factor << std::endl;
            xVal.ScalarMult(factor);
          } else {
            traceMsg << "X is zero; scaling will not help; try a scaled version of yVal instead" << std::endl;
            xVal = yVal;
            xVal.ScalarMult(xVal.NormL2()*XSaturated_/PSaturated_);
          }

          //xVal.ScalarMult(0.85);
          //xVal.Init();
          //            }
        }
        
      }
      
      //xVal.Init();
      
      traceMsg << "starting Xval: " << xVal.ToString() << std::endl;
      assert(!xVal.ContainsNaN() && !xVal.ContainsInf());
      /*
       * In the followng, prev*val_ = value of previous iteration
       *                      *val_ = current value
       */
      
      // perform some initial steps with fp iteation
      /*
       * Note: one iteration seems fine; two might already move solution so
       *       far off that it does not recover
       * > 1 works ok; 
       * (mostly due to bringing x into sat and then resetting it
       */
      //        ClipDirection(xVal);
      hystVal = computeValue_vec(xVal, operatorIndex, overwriteMemory, overwriteDirection, debugOut, successFlagForward);
      LOG_DBG(vecpreisachInversion) << "starting Hval: " << hystVal.ToString();
      //hystVal = prevHystval_[idElem];
      
      /*
       * As y is not in saturation (that case was already ruled out, 
       * the actual xVal cannot be in saturation neither
       * > It seems that initial fp iteration only works due to the following capping;
       * FP with the original mu nearly always put the system over saturation; the following
       * reset will cap it below saturation and from this point on the LM method seems
       * to converge much better than using the previous value of x
       * TODOL
       * > check if fp works better if xsat/ysat is used instead of mu
       * > check if lm works better if reset is done to a smaller percentage of xsat (currently 0.9)
       */
      /*
       * IDEA: take direction of new Y and scale it
       *  > only if y != 0!
       *  if y == 0, take xDirecttion
       * 
       */
      Vector<Double> dirX = Vector<Double>(dim_);
      if(xVal.NormL2() != 0){
        dirX = xVal;
        dirX.ScalarDiv(xVal.NormL2());
      } else if(yNorm > 0){
        dirX = yVal;
        dirX.ScalarDiv(yNorm);
      } else {
        dirX.Init();
      }
            
      if((xVal.NormL2() >= XSaturated_)&&(stayBelowSat==1)){
        traceMsg << "Reset xVal as its value is above Xsaturation but yVal is not" << std::endl;

        xVal = dirX;
        xVal.ScalarMult(1.0/1.2*XSaturated_);
        
        //xVal.ScalarMult(percent*XSaturated_/xVal.NormL2());
        // or reset to 0? > 0.9 works good; 
        //xVal.Init();
        //resetAfterSat=true;
        //          ClipDirection(xVal);
        hystVal = computeValue_vec(xVal, operatorIndex, overwriteMemory, overwriteDirection, debugOut, successFlagForward);
      } 
      
      if((xVal.NormL2() < XSaturated_)&&(stayBelowSat==-1)){
        traceMsg << "Reset xVal as its value is below Xsaturation but yVal is above" << std::endl;

        xVal = dirX;
        xVal.ScalarMult(1.2*XSaturated_);
        
        //xVal.ScalarMult(percent*XSaturated_/xVal.NormL2());
        // or reset to 0? > 0.9 works good; 
        //xVal.Init();
        //resetAfterSat=true;
        //          ClipDirection(xVal);
        hystVal = computeValue_vec(xVal, operatorIndex, overwriteMemory, overwriteDirection, debugOut, successFlagForward);
      } 

      Vector<Double> xStart = xVal;
      
      Double sign = 1.0;
      bool successError = false;
      bool successX = false;
      bool successY = false;
      bool discardUpdate = false;
      
      Vector<Double> xUpdate = Vector<Double>(dim_);
      Vector<Double> res = Vector<Double>(dim_);
      Matrix<Double> jac = Matrix<Double>(dim_,dim_);
      Matrix<Double> jacT = Matrix<Double>(dim_,dim_);
      
      Double errorNorm;
      Double errorNormResX;
      Double errorNormResY;
      UInt implementation = 0;
      Vector<Double> bestSol = Vector<Double>(dim_);
      Double bestErrorNorm = 1e16;
      
      Vector<Double> bestSolRes = Vector<Double>(dim_);
      Double bestErrorNormRes = 1e16;
      
      /*
       * Start actual LM iteration
       */
      while(true){ 
        itCnt++;
        totalNumberOfLMIterations++;
        
        if(debug){
          traceMsg << "--OUTER ITERATION-- (" << itCnt << ")" << std::endl;
          traceMsg << "Current solution: " << xVal.ToString() << std::endl;
        }
        // do not override hyst memory here
        //          ClipDirection(xVal);
        hystVal = computeValue_vec(xVal, operatorIndex, overwriteMemory, overwriteDirection, debugOut, successFlagForward);
        
        res = computeResidual(xVal,yVal,hystVal,mu,mu_inv,wrtX,relError);
        jac = computeJacobian(xVal,yVal,hystVal,res,mu,mu_inv,operatorIndex,
                sign,wrtX,relError,implementation,overwriteMemory,overwriteDirection,stayBelowSat);
        jac.Transpose(jacT);
        
        if(debug){
          traceMsg << "Current hystVal: " << hystVal.ToString() << std::endl;
          traceMsg << "Current Residual (wrt x): " << res.ToString() << std::endl;
          traceMsg << "Current startvalue for alpha_LS: " << alpha << std::endl;
          traceMsg << "Current JacobiMatrix: " << jac.ToString() << std::endl;
          
          Vector<Double> jacTres = Vector<Double>(dim_);
          jacT.Mult(res,jacTres);
          traceMsg << "jacTres: " << jacTres.ToString() << std::endl;
          traceMsg << "tolForErrorVector: " << INV_resTolH_ << std::endl;
          traceMsg << "tolForResidual X: " << INV_resTolH_ << std::endl;
          traceMsg << "tolForResidual Y: " << INV_resTolB_ << std::endl;
        }
        
        successError = checkConvergence(res,jacT,errorNorm,INV_resTolH_);
        errorNormResX = res.NormL2();
        successX = (errorNormResX <= INV_resTolH_);
        successY = checkInversionOutput(xVal, yVal, mu, INV_resTolB_, errorNormResY,
                operatorIndex, overwriteMemory, overwriteDirection);
        
        if(successError){
          // main criterion > would be best if this one could be satisfied
          traceMsg << "Success! Error estimate |jacT*ResX| = " << errorNorm << " < " << INV_resTolH_ << std::endl;
          successFlag = 4;
          LOG_TRACE(vecpreisachInversion) << "Success! Error estimate |jacT*ResX| = " << errorNorm << " < " << INV_resTolH_ << std::endl;
          break;
        } else if(successX){
          // failback; still top if this works
          traceMsg << "Success! Residual norm wrt X = |ResX| = " << errorNormResX << " < " << INV_resTolH_ << std::endl;
          successFlag = 5;
          LOG_TRACE(vecpreisachInversion) << "Success! Residual norm wrt X = |ResX| = " << errorNormResX << " < " << INV_resTolH_ << std::endl;
          break;
        } else {
          if( (totalNumberOfLMIterations%10 == 0)||(itCnt >= INV_maxIter_) ){
            if(successY){
              // failback; might still have large res error in x buz might be the best we can find
              // check only each 10th iteration
              LOG_TRACE(vecpreisachInversion) << "Success! Residual norm wrt Y = |ResY| = " << errorNormResY << " < " << INV_resTolB_ << std::endl;
              traceMsg << "Success! Residual norm wrt Y = |ResY| = " << errorNormResY << " < " << INV_resTolB_ << std::endl;
              successFlag = 6;
              break;
            } else if (itCnt >= INV_maxIter_) {
              LOG_TRACE(vecpreisachInversion) << "NO Success! Remaining residual norm wrt Y = |ResY| = " << errorNormResY << " < " << INV_resTolB_ << std::endl;
              if(debug){
                traceMsg << "Max number of iterations reached." << std::endl;
                traceMsg << "Last found solution xFound = " << xVal.ToString() << std::endl;
                traceMsg << "p(xFound) = " << hystVal.ToString() << std::endl;
                traceMsg << "Remaining residual wrt x-value (xFound - nu*(yIn - p(xFound)): " << res.ToString() << std::endl;
                traceMsg << "Remaining error norm |jacT*ResX|: " << errorNorm << std::endl;
                traceMsg << "Remaining residual-norm wrt x: " << errorNormResX << std::endl;
                traceMsg << "Remaining residual-norm wrt y: " << errorNormResY << std::endl;
              }
              //                  std::cout << "Remaining error norm |jacT*ResX|: " << errorNorm << std::endl;
              //                  std::cout << "Remaining residual-norm wrt x: " << errorNormResX << std::endl;
              //                  std::cout << "Remaining residual-norm wrt y: " << errorNormResY << std::endl;
              successFlag = 0;
              
//              std::cout << "previousStateAboveSat? " << previousStateAboveSat << std::endl;
//              std::cout << "currentStateAboveSat? " << currentStateAboveSat << std::endl;
//              
//              if(previousStateAboveSat != currentStateAboveSat){
//                successCode = 10;
//              }
              
              break;
            }
          }
        }
        
        if(errorNormResX < bestErrorNormRes){
          bestErrorNormRes = errorNormResX;
          bestSolRes = xVal;                
        }
        
        if(errorNorm < bestErrorNorm){
          bestErrorNorm = errorNorm;
          bestSol = xVal;                
        }
        
        LOG_DBG(vecpreisachInversion) << "Solution not appropriate yet > Perform linesearch";
        UInt numberOfIterations = 0;
        
        // we have troubles if we come close to saturation
        // in that case we often shoot over saturation and have to step back
        // if we step back to the same value each time, we will never come forward, however
        Double factorToSat = 0.1*Double(itCnt-1)/Double(INV_maxIter_) + 0.9;
        
        discardUpdate = performLinesearch(xVal, yVal, res, xUpdate, jac, jacT, mu, mu_inv, 
                operatorIndex, overwriteMemory, overwriteDirection, alpha, alphaMin, alphaMax, wrtX, 
                relError, numberOfIterations,xStart,factorToSat,stayBelowSat,sol);
        
        if(alpha < minAlphaStatistics){
          minAlphaStatistics = alpha;
        }
        if(alpha > maxAlphaStatistics){
          maxAlphaStatistics = alpha;
        }
        avgAlphaStatistics += alpha;
	        
        traceMsg << "Computed update: " << xUpdate.ToString() << std::endl;
        traceMsg << "Discard update? " << discardUpdate << std::endl;
        
        totalNumberOfLinesearchIterations += numberOfIterations;
        
        if(numberOfIterations > maximalNumberOfLinesearchIterations){
          maximalNumberOfLinesearchIterations = numberOfIterations;
        }
        
        /*
         * New idea: start with tighter bounds for alpha but shift these limits
         * after each 10 iterations
         */
//        if(totalNumberOfLMIterations%10 == 0){
//          alphaMin = alphaMin/2.0;
//          alphaMax = alphaMax*2.0;
//          // for Tikhonov > set xStart to xVal
//          xStart = bestSol; // own test
//        } 
        
        LOG_DBG(vecpreisachInversion) << "Computed update: " << xUpdate.ToString();
        if(!discardUpdate){
          xVal = xVal+xUpdate;
        } else {
          LOG_DBG(vecpreisachInversion) << "Discard update";
					//std::cout << "Discard update; reset to best solution so far" << std::endl;
					alphaMin = alphaMin/2.0;
          alphaMax = alphaMax*2.0;
          // for Tikhonov > set xStart to xVal
          //xVal = bestSol; // own test
//					// reset solution
//					if(stayBelowSat==1){
//						std::cout << "Reset to 0" << std::endl;
//						xVal.Init();
//					} else if(stayBelowSat==-1){
//						std::cout << "Reset to random direction*2*Xsaturation" << std::endl;
//						xVal.Init();
//						for(UInt i = 0; i < dim_; i++){
//							xVal[i] = (rand()%200-100)/100.0;
//						}
//						if(xVal.NormL2() == 0){
//							xVal[0] = 1;
//						}
//						xVal.ScalarMult(2.0*XSaturated_/xVal.NormL2());
//						
//						std::cout << "xVal: " << xVal.ToString() << std::endl;
//					} else {
//						std::cout << "Reset to random direction*Xsaturation" << std::endl;
//						xVal.Init();
//						for(UInt i = 0; i < dim_; i++){
//							xVal[i] = (rand()%200-100)/100.0;
//						}
//						if(xVal.NormL2() == 0){
//							xVal[0] = 1;
//						}
//						xVal.ScalarMult(1.0*XSaturated_/xVal.NormL2());
//						
//						std::cout << "xVal: " << xVal.ToString() << std::endl;
//					}	
        }
        sign = sign*(-1.0);
      }
      
      avgAlphaStatistics /= totalNumberOfLMIterations;
      
      if((xVal.NormL2() > XSaturated_)&&(stayBelowSat==1)){
        EXCEPTION("LM lead xVal into saturation although input is below > must not be the case!");
      }
      
      if((xVal.NormL2() < XSaturated_)&&(stayBelowSat==-1)){
        EXCEPTION("LM lead xVal below saturation although input is above> must not be the case!");
      }
    }
    // end LM
    
    /*
     * set values for next time
     */
    if(overwriteMemory == true){
      EXCEPTION("Memory should not be overridden here");
    }
    
    if(debug && successFlag==0){
      Double resYNorm = 0.0;
      checkInversionOutput(xVal, yVal, mu, INV_resTolB_, resYNorm, operatorIndex, overwriteMemory, overwriteDirection,true);
      LOG_TRACE(vecpreisachInversion) << traceMsg.str();
    }
    
    LOG_TRACE(vecpreisachInversion) << " --------- END IVERSION --------- ";
    
    return xVal;
  } 
  
	// returns true, if yTarget can be reached by mu*xComputed + hyst(xComputed)
	bool Hysteresis::checkInversionOutput(Vector<Double>& xComputed, Vector<Double>& yTarget, 
          Matrix<Double>& mu, Double tol, Double& resYNorm, Integer operatorIdx, bool overwriteMemory, bool overwriteDirection, bool output){
		
		Vector<Double> yCheck = Vector<Double>(dim_);
    bool debugOut = false;
    int successCode = 0;
    Vector<Double> hCheck = computeValue_vec(xComputed, operatorIdx, overwriteMemory, overwriteDirection, debugOut, successCode);
    
		mu.Mult(xComputed,yCheck);
		yCheck.Add(hCheck);
		
		Vector<Double> diff;
		diff = yCheck;
		diff -= yTarget;
    resYNorm = diff.NormL2();
    if(output){
      LOG_TRACE(vecpreisachInversion) << "CheckInversionOutput - norm of residual wrt y-vector: " << diff.NormL2();
    }
		if(diff.NormL2() < tol){
      if(output){
        LOG_TRACE(vecpreisachInversion) << "Solution OK; |diff| = " << diff.NormL2() << " < " << tol;
      }
			return true;
		} else {
      if(output){
        LOG_TRACE(vecpreisachInversion) << "Solution NOT OK; |diff| = " << diff.NormL2() << " > " << tol;
        LOG_TRACE(vecpreisachInversion) << "yRetrieved = " << yCheck.ToString();
        LOG_TRACE(vecpreisachInversion) << "yIn = " << yTarget.ToString();
      }
			return false;
		}
	}
  
}
