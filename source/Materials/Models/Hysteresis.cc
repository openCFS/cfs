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
      Pout = computeValue_vec(Pin, idx, false, debugOut, successFlag);
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
      Pout = computeValue_vec(Pin, idx, false, debugOut, successFlag);
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
          Pout = computeValue_vec(Pin, idx, false, debugOut, successFlag);
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
  
  Integer Hysteresis::checkIncrementTrustRegion(Vector<Double>& x_new, 
          Vector<Double>& res_cur, Vector<Double>& res_new,
          Vector<Double>& jac_dx, Double& alpha, int stayBelowSat){
    /*
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
     * <0 : dicard update
     * =0 : keep update
     * =1 : keep update
     * 
     * Input:
     *  x_new = current x + increment dx
     *  res_cur = F(x) > residual evaluated at current x
     *  res_new = F(x + increment) > residual evaluated at new x
     *  jac_dx = jacobian evaluated at current x times increment
     *  alpha = stepping value that was used in increment
     *          > Levenberg Marquardt: alpha from regularization
     *          > Newton, Krylov-Newton: alpha from linesearch
     *  stayBelowSat = extra flag for narrowing down acceptable 
     *                ´+1 -> discard update if x_new > sat
     *                 -1 -> discard update if x_new < sat
     *                  0 -> do not check x_new
     */
    Double beta0, betaLow, betaUp;
    beta0 = 0.15;
    betaLow = 0.35;
    betaUp = 0.85;
    
//    beta0 = 0.2;
//    betaLow = 0.4;
//    betaUp = 0.8;
    
    Double factorUp = 2;
    Double factorDown = 0.5;
    
    Integer success = 0;
    if((x_new.NormL2() > XSaturated_)&&(stayBelowSat==1)){
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
    } else if((x_new.NormL2() < XSaturated_)&&(stayBelowSat==-1)){
      alpha = alpha/factorUp;
      success = -2;
    } else {
      /*
       * rho = (||F(x)||^2 - ||F(x + increment)||^2) /
       *       (||F(x||^2 - ||F(x) + F'(x)*increment||^2)
       */
      Vector<Double> res_jac_dx = Vector<Double>(dim_);
      res_jac_dx.Add(1.0,res_cur,1.0,jac_dx);
      Double res_cur_square = res_cur.NormL2_squared();
      
      Double nominator = res_cur_square - res_new.NormL2_squared();
      Double denominator = res_cur_square - res_jac_dx.NormL2_squared();
      
      LOG_TRACE(vecpreisachInversion) << "Nominator: ||F(x)||^2 - ||F(x + increment)||^2 = " << nominator;
      LOG_TRACE(vecpreisachInversion) << "Denominator: ||F(x||^2 - ||F(x) + F'(x)*increment||^2 = " << denominator;
      
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
      if(denominator != 0){
        rho = nominator/denominator;
      } else {
        rho = nominator;
        if(nominator == 0){
          //LOG_DBG(vecpreisach) << "Compute rho: a = 0, too!";
          // rho is actually undefined here (as increment is basically 0)
          // > alpha seems to be way too large (otherwise increment would not be 0)
          // > return rho >= 1 in order to trigger a much smaller alpha
          // > rho >= 10 > reset alpha
          rho = 10.0;
        }
      }
      
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
		Double resNorm = res.NormL2_squared();
		Double resShiftedNorm = resShifted.NormL2_squared();
		Vector<Double> tmp = Vector<Double>(dim_);
    //    tmp = computeJacobianTimesVector(x, Vector<Double>& v, 
    //          Vector<Double>& y, Vector<Double>& hyst_x, Matrix<Double> mu_inv, Integer operatorIdx)
		jac.Mult(xUpdate,tmp);
		
		tmp = tmp + res;
		Double tmpNorm = tmp.NormL2_squared();
		Double a = resNorm - resShiftedNorm;
		Double b = resNorm - tmpNorm;
		
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
  
  Vector<Double> Hysteresis::computeResidual(Vector<Double>& xVal, Vector<Double>& yVal, Vector<Double>& hystVal, Matrix<Double> mu_inv){
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
    
    return res;
  }
  
  Vector<Double> Hysteresis::obtainStartVector(Vector<Double> previousXVector, Vector<Double> previousYVector,
          Vector<Double> currentYVector, Double YdiffToRemancence, Double YdiffToSaturation, int stayBelowSat){
    
    UInt size = previousXVector.GetSize();
    Vector<Double> startXVector = Vector<Double>(size);
    
    if (YdiffToRemancence < INV_resTolB_){
      // we are quite close to remanence, but have not yet reached it
      // start from 0
      startXVector.Init();
      //  startAtZero=true;
      //	compCase = 3;
      //    } else if (YdiffToSaturation < 100*INV_resTolB_){
      //      startXVector.Init();
      //      
      //      Double yNorm = currentYVector.NormL2();
      //      Vector<Double> directionYVector = Vector<Double>(size);
      //      directionYVector.Init();
      //      
      //      if(yNorm == 0){
      //        yNorm = previousYVector.NormL2();
      //        if(yNorm == 0){
      //          yNorm = 1.0;
      //        }
      //        directionYVector.Add(1.0/yNorm,previousYVector); 
      //      } else {
      //        directionYVector.Add(1.0/yNorm,currentYVector); 
      //      }
      //      
      //      startXVector.Add( (1.0-YdiffToSaturation)*XSaturated_, directionYVector);
    } else {
      //	compCase = 4;
      Double factor;
      if(startXVector.NormL2() != 0){
        startXVector = previousXVector;
        
        if(previousYVector.NormL2() != 0){
          factor = currentYVector.NormL2()/previousYVector.NormL2();
        } else {
          if(currentYVector.NormL2() > previousYVector.NormL2()){
            factor = 1.05;
          } else {
            factor = 1.0/1.05;
          }
        }
        startXVector.ScalarMult(factor);
      } else {
        //            traceMsg << "X is zero; scaling will not help; try a scaled version of yVal instead" << std::endl;
        startXVector = currentYVector;
        startXVector.ScalarMult(startXVector.NormL2()*XSaturated_/PSaturated_);
      }
    }
    //      traceMsg << "starting Xval: " << xVal.ToString() << std::endl;
    assert(!startXVector.ContainsNaN() && !startXVector.ContainsInf());
    
    Vector<Double> dir = Vector<Double>(dim_);
    dir.Init();
    if(startXVector.NormL2() != 0){
      dir.Add(1.0/startXVector.NormL2(),startXVector);
    } else if(currentYVector.NormL2() != 0){
      dir.Add(1.0/currentYVector.NormL2(),currentYVector);
    } 
    
    //    std::cout << "startXVector.NormL2(): " << startXVector.NormL2() << std::endl;
    //    std::cout << "currentYVector.NormL2(): " << currentYVector.NormL2() << std::endl;
    //    std::cout << "dir: " << dir.ToString() << std::endl;
    //    std::cout << "XSaturated_: " << XSaturated_ << std::endl;
    if((startXVector.NormL2() >= XSaturated_)&&(stayBelowSat==1)){
      //        traceMsg << "Reset xVal as its value is above Xsaturation but yVal is not" << std::endl;
      //      std::cout << "scale down" << std::endl;
      startXVector = dir;
      startXVector.ScalarMult(1.0/1.2*XSaturated_);
    } 
    
    if((startXVector.NormL2() < XSaturated_)&&(stayBelowSat==-1)){
      //        traceMsg << "Reset xVal as its value is below Xsaturation but yVal is above" << std::endl;
      //      std::cout << "scale up" << std::endl;
      startXVector = dir;
      startXVector.ScalarMult(1.2*XSaturated_);        
    }
    //    std::cout << "startXVector.NormL2(): " << startXVector.NormL2() << std::endl;
    return startXVector;
  }
  
//  Hysteresis::performLinesearch(){
//    
//    /*
//     * Two cases:
//     *  A) Newton/Krylov-Newton:
//     *    1. compute update dx
//     *    2. while()
//     *       2.1. get scaled update dx_scal = alpha*dx
//     *       2.2. check dx_scal using trustregion
//     *          2.2.1. update ok > take dx_scal and leave
//     *          2.2.2. update not ok > take new alpha and try again
//     *          2.2.3. update has to be discarded > problem!
//     * 
//     *  B) Levenberg-Marquardt
//     *    1. while()
//     *       1.1. compute update dx_alpha via regualirzed system depending on alpha
//     *       1.2. check dx_alpha using trustregion
//     *          1.2.1. update ok > take dx_alpha and leave
//     *          1.2.2. update not ok > take new alpha and try again
//     *          1.2.3. update has to be discarded > problem!
//     */       
//    
//  }
//  
  
  Vector<Double> Hysteresis::computeUpdate_Newton(Vector<Double>& x, Vector<Double>& y, 
          Vector<Double>& hyst_x, Matrix<Double> mu_inv, Integer operatorIdx){
    
    /*
     * Jac*dX = -res
     * > solve directly
     */
    Matrix<Double> jac = computeJacobian(x, hyst_x, 
            mu_inv, operatorIdx, 1.0, 0, false, 0); 
    UInt vecSize = x.GetSize();
    Matrix<Double> jacInv = Matrix<Double>(vecSize,vecSize);
    jac.Invert(jacInv);
    
    Vector<Double> res = computeResidual(x,y,hyst_x,mu_inv);
    Vector<Double> dx = Vector<Double>(vecSize);
    jacInv.Mult(res,dx);
    
    dx.ScalarMult(-1.0);
    return dx;
  }
  
  Vector<Double> Hysteresis::computeUpdate_Krylov(Vector<Double>& x, Vector<Double>& y, 
          Vector<Double>& hyst_x, Matrix<Double> mu_inv, Integer operatorIdx){
    
    /*
     * Jac*dX = -res
     * > solve iteratively by Krylov space method
     * 
     * Compute update for Newton-type iteration using a Krylov-space approach
     * > source: 
     *  "Jacobian-free Newton-Krylov methods: a survey of approaches and applications"
     * > i.e. use GMRES
     * x = current solution
     * deltaX_0 = initial guess 
     * r0 = -res(x) - Jac(x) deltaX_0
     * deltaX_j = deltaX_0 + sum_i=0^j-1 beta_i Jac^i r_0
     *  > continue till deltaX_j is ok
     *  > beta_i might be found via linesearch
     * xNew = x + deltaX_j
     */
    
    /*
     * GMRES solver taken from Wikipedia
     * https://en.wikipedia.org/wiki/Generalized_minimal_residual_method
     */
    
    Double tol = 1e-15;
    // Solve 
    //  Jac*dX = -res 
    //  A*x = b
    // r0 = b - A*x = -res - Jac*dX
    Vector<Double> dx = Vector<Double>(dim_);
    dx.Init();
    Vector<Double> res = computeResidual(x,y,hyst_x,mu_inv);
    Vector<Double> jacdx = computeJacobianTimesVector(x, dx, y, hyst_x, mu_inv, operatorIdx);
    
    Vector<Double> r0 = Vector<Double>(dim_);
    r0.Add(-1.0,res,-1.0,jacdx);
    
    if(r0.NormL2() < tol){
      return dx;
    }
    UInt n = dim_;
    
    Vector<Double>* vArray = new Vector<Double>[n+1];
    Vector<Double>* wArray = new Vector<Double>[n];
    Matrix<Double> h = Matrix<Double>(n+1,n);
    Vector<Double> gamma = Vector<Double>(n+1);
    Vector<Double> coefs = Vector<Double>(n);
    
    vArray[0] = Vector<Double>(n);
    vArray[0].Init();
    vArray[0].Add(1.0/r0.NormL2(),r0);
    
    gamma[0] = r0.NormL2();
    
    Vector<Double> q = Vector<Double>(n);
    Vector<Double> c = Vector<Double>(n+1);
    Vector<Double> s = Vector<Double>(n+1);
    Double tmp1,tmp2,beta;
    UInt j;
    for(j = 0; j < n; j++ ){
      
      // q = Av_j = Jac*v_j
      q = computeJacobianTimesVector(x, vArray[j], y, hyst_x, mu_inv, operatorIdx);
      
      /*
       * Setup Arnoldi space
       */
      for(UInt i = 0; i <= j; i++){
        h[i][j] = vArray[i].Inner(q);
        q.Add(-h[i][j],vArray[i]);
      }
      
      wArray[j] = Vector<Double>(n);
      wArray[j] = q;
      
      h[j+1][j] = wArray[j].NormL2();
      
      for(UInt i = 0; i < j; i++){
        tmp1 = c[i+1]*h[i][j] + s[i+1]*h[i+1][j];
        tmp2 = -s[i+1]*h[i][j] + c[i+1]*h[i+1][j];
        h[i][j] = tmp1;
        h[i+1][j] = tmp2;
      }
      
      /*
       * Givens rotation
       */
      beta = std::sqrt(h[j][j]*h[j][j] + h[j+1][j]*h[j+1][j]);
      s[j+1] = h[j+1][j]/beta;
      c[j+1] = h[j][j]/beta;
      h[j][j] = beta;
      
      gamma[j+1] = -s[j+1]*gamma[j];
      gamma[j] = c[j+1]*gamma[j];
      
      if(abs(gamma[j+1]) >= tol){
        vArray[j+1] = Vector<Double>(n);
        vArray[j+1].Init();
        vArray[j+1].Add(1.0/h[j+1][j],wArray[j]);
      } 
      else {
        break;
      }
    }
    
    /*
     * construct update
     */
    for(int i = (int) j; i >= 0; i--){
      coefs[i] = gamma[i];
      for(UInt k = i+1; k <= j; k++){
        coefs[i] -= h[i][k]*coefs[k];
      }
      if(abs(h[i][i])<1e-19){
        coefs[i] = 0.0;
      } else {
        coefs[i] /= h[i][i];
      }
    }
    for(UInt i = 0; i <= j; i++){
      dx.Add(coefs[i],vArray[i]);
    }
    return dx;
  }
  
  Vector<Double> Hysteresis::computeJacobianTimesVector(Vector<Double>& x, Vector<Double>& v, 
          Vector<Double>& y, Vector<Double>& hyst_x, Matrix<Double> mu_inv, Integer operatorIdx){
    
    /*
     * Function for Jacobian-Free-Newton-Krylov
     * Instead of computing/approximating the Jacobian and then compute
     * Jacobian*vecForMultiplication
     * we instead approximate
     * Jacobian*vecForMultiplication directly via
     * 
     *  Jac*v \approx [F(x + eps*v) - F(x)]/eps
     * 
     * with 
     *  F(x) = residual wrt x = x + mu_inv( hyst_x/mu - y )
     *  x = current solution
     *  v = vector for multiplication
     *  hyst_x = Hystoperator evaluated at x
     *  y = target vector
     */
    Double eps;
    if(v.NormL2() != 0){
      eps = std::sqrt((1.0 + x.NormL2())*1e-15)/v.NormL2();
    } else {
      eps = 1e-16;
    }
    //    std::cout << "--- eps: " << eps << std::endl;
    UInt vecSize = x.GetSize();
    Vector<Double> res_x = computeResidual(x,y,hyst_x,mu_inv);
    Vector<Double> xv = Vector<Double>(vecSize);
    xv.Init();
    xv.Add(1.0,x,eps,v);
    
    bool overwriteMem = false;
    bool debugOut = false;
    int successFlag = 0;
    Vector<Double> hyst_xv = computeValue_vec(xv, operatorIdx, overwriteMem, debugOut, successFlag);
    Vector<Double> res_xv = computeResidual(xv,y,hyst_xv,mu_inv);
    Vector<Double> jacv = Vector<Double>(vecSize);
    jacv.Init();
    jacv.Add(1.0/eps,res_xv,-1.0/eps,res_x);
    
    /*
     * for testing:
     * compute Jacobian and multiply by v
     */
    bool testing = false;
    if(testing){
      Matrix<Double> jac = computeJacobian(x, hyst_x, 
              mu_inv, operatorIdx, 1.0, 0, false, 0); 
      
      Vector<Double>tmp = Vector<Double>(vecSize);
      jac.Mult(v,tmp);
      Vector<Double>jacv2 = Vector<Double>(vecSize);
      jacv2 = tmp;
      
      tmp.Add(-1.0,jacv);
      //      std::cout << "--- jac*v - jacv = " << tmp.ToString() << std::endl;
      
      return jacv2;
    }
    
    return jacv;
  }
  
  Matrix<Double> Hysteresis::computeJacobian(Vector<Double>& xVal, Vector<Double>& hystVal, 
          Matrix<Double> mu_inv, Integer operatorIdx, Double sign, UInt implementation, 
					bool overwriteMemory, int stayBelowSat) {
    
    LOG_DBG(vecpreisachInversion) << " --------- computeJacobian --------- ";
    //    LOG_DBG(vecpreisach) << "VecPreisach::computeJacobian";
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
        hystShifted = computeValue_vec(xShifted, operatorIdx, overwriteMemory, debugOut, successFlag);
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
        
        Vector<Double> curCol = Vector<Double>(dim_);
        mu_inv.Mult(deltaHyst,curCol);
        curCol.ScalarDiv(deltaX);
        
        jac[i][i] = 1.0;
        for(UInt j = 0; j < dim_; j++){
          jac[j][i] += curCol[j];
        }
      }      
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
        hystShifted = computeValue_vec(xShifted, operatorIdx, overwriteMemory, debugOut, successFlag); 
        xShifted_opp = xVal;
        xShifted_opp[i] -= deltaX;
        hystShifted_opp = computeValue_vec(xShifted_opp, operatorIdx, overwriteMemory, debugOut, successFlag);
        
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
        }
      }
    } 
    //    LOG_DBG(vecpreisach) << "Retrieved Jacobian: " << jac.ToString();
    return jac;
  }
  
  Vector<Double> Hysteresis::computeUpdate_LM(Vector<Double> jacTres_neg, 
          Matrix<Double>& jacTjac, Double& alpha, Double& alphaAcc, Double& alphaMinReq, Double alphaMax){
    
    Matrix<Double> matToInvert = Matrix<Double>(dim_,dim_);
    Matrix<Double> matInverted = Matrix<Double>(dim_,dim_);
    Vector<Double> dx = Vector<Double>(dim_);
    
    Double minDeterminant = 1e-12;  
    Double incrFactor = std::sqrt(2.0);
    Double detMatToInvert;
    UInt cnt = 0;
    
    matToInvert = jacTjac;
    
    /*
     * Add regularization
     * > note: it works better if alpha is accumulated for some reason
     */
    alphaAcc += alpha*alpha;
    for(UInt i = 0; i < dim_; i++){
      matToInvert[i][i] += alphaAcc;
    }
    matToInvert.Determinant(detMatToInvert);
    
    if(abs(detMatToInvert)<minDeterminant){
      /*
       * current alpha will not work
       * > reset alphaAcc and find out minimal alpha
       */
      while(true){
        matToInvert = jacTjac;
        // no += here; here we actually reset alphaAcc
        alphaAcc = alpha*alpha;
        
        for(UInt i = 0; i < dim_; i++){
          matToInvert[i][i] += alphaAcc;
        }
        matToInvert.Determinant(detMatToInvert);
        
        if(abs(detMatToInvert)>=minDeterminant){
          break;
        } else {
          alpha = alpha*incrFactor;
        }
        cnt++;
        if( (cnt > 200)||(alpha > alphaMax) ){
          EXCEPTION("LM: Cannot find alpha, such that jacTjac becomes invertible!");
        }
      }
    }
    
    alphaMinReq = alpha;
    matToInvert.Invert(matInverted);
    matInverted.Mult(jacTres_neg,dx);
    
    assert(!dx.ContainsNaN() && !dx.ContainsInf());
    
    return dx;
  }
  
  bool Hysteresis::computeUpdate_LM_full(Vector<Double>& xVal, Vector<Double>& yVal, Vector<Double>& res, 
          Vector<Double>& xUpdate, Matrix<Double>& jac, Matrix<Double>& jacT, Matrix<Double> mu, Matrix<Double> mu_inv, 
          Integer operatorIdx, bool overwriteMemory,
          Double& alpha, Double alphaMin, Double alphaMax,
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
    
		UInt maxIter = 150;
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
 
//    
//    
//    LOG_DBG(vecpreisachlinesearch) << "Determine minimal alpha that allows inversion of jacTjac";
//    LOG_DBG(vecpreisachlinesearch) << "Given alphaMin: " << alphaMin;        
//    Double minDeterminant = 1e-12;  
//    Double incrFactor = std::sqrt(2.0);
//    while(true){
//      matToInvert = jacTjac;
//      jacTjac.Determinant(detMatToInvert);
//      LOG_DBG(vecpreisachlinesearch) << " det(jacTjac) =  " << detMatToInvert;
//      
//      for(UInt i = 0; i < dim_; i++){
//        matToInvert[i][i] += alphaMinLocal*alphaMinLocal;
//      }
//      matToInvert.Determinant(detMatToInvert);
//      if(abs(detMatToInvert)>=minDeterminant){
//        //      if(detMatToInvert != 0){
//        break;
//      } else {
//        alphaMinLocal = alphaMinLocal*incrFactor;
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
    
    
    
    
    //    // now reset alpha to alphaMax (note: alphaMax is passed as copy, i.e. if we increase it further
    //    // below, the actual value of alphaMax is restored during next call
    //    if(alpha > alphaMax){
    //      LOG_DBG(vecpreisachlinesearch) << "Current alpha > alphaMax";
    //      alpha = alphaMax;
    //    }
    
    // start at alphamin
    alpha = alphaMinLocal;
    Double alphaAcc = alpha*alpha;
    while(true){
      itCnt++;
			// for statistics
			numberOfIterations=itCnt;
      
      LOG_DBG(vecpreisachlinesearch) << "Start Iteration " << itCnt << " with: ";
      LOG_DBG(vecpreisachlinesearch) << "xVal = " << xVal.ToString();
      LOG_DBG(vecpreisachlinesearch) << "last found x = " << xNew.ToString();
      LOG_DBG(vecpreisachlinesearch) << "alpha = " << alpha;
      
      
//      
//      matToInvert = jacTjac;
//      // do not reset alphaAcc here > with accumulated alpha it works better
//      //alphaAcc = 0;
//      alphaAcc += alpha*alpha;
//      for(UInt i = 0; i < dim_; i++){
//        matToInvert[i][i] += alphaAcc;
//      }
//      matToInvert.Determinant(detMatToInvert);
//      
//      if(abs(detMatToInvert)<minDeterminant){
//        //      if(detMatToInvert == 0){
//        //        std::cout << "Matrix not invertible! " << std::endl;
//        //        std::cout << "matToInvert: " << matToInvert.ToString() << std::endl;
//        //        std::cout << "Min alpha for inversion: " << alphaMinLocal << std::endl;
//        //        std::cout << "Current alpha: " << alpha << std::endl;
//        //        std::cout << "Try to find alpha that works again" << std::endl;
//        while(true){
//          matToInvert = jacTjac;
//          alphaAcc = 0;
//          alphaAcc += alpha*alpha;
//          jacTjac.Determinant(detMatToInvert);
//          LOG_DBG(vecpreisachlinesearch) << " det(jacTjac) =  " << detMatToInvert;
//          
//          for(UInt i = 0; i < dim_; i++){
//            matToInvert[i][i] += alphaAcc;
//          }
//          matToInvert.Determinant(detMatToInvert);
//          if(abs(detMatToInvert)>=minDeterminant){
//            //      if(detMatToInvert != 0){
//            break;
//          } else {
//            alpha = alpha*incrFactor;
//          }
//          cnt++;
//          if(cnt > 200){
//            EXCEPTION("LM, Linesearch: Cannot find alpha, such that jacTjac becomes invertible!");
//          }
//        }
//        alphaMinLocal = alpha;
//      }
//      assert(detMatToInvert != 0);
//      matToInvert.Invert(matInverted);
//      
//      matInverted.Mult(jacTres_neg,xUpdate);
//      
//      assert(!xUpdate.ContainsNaN() && !xUpdate.ContainsInf());
//      LOG_DBG(vecpreisachlinesearch) << "Computed xUpdate = " << xUpdate.ToString();
//      
      
      
      
      
      xUpdate = computeUpdate_LM(jacTres_neg, 
          jacTjac, alpha, alphaAcc, alphaMinLocal, 1e16);
//      
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
//      Vector<Double> hystSol = computeValue_vec(sol, operatorIdx, overwriteMemory, debugOut, successFlag);
//      Vector<Double> hystOld = computeValue_vec(xVal, operatorIdx, overwriteMemory, debugOut, successFlag);
//      Vector<Double> resSol = computeResidual(sol,yVal,hystSol, mu_inv);
      //      Vector<Double> resSol = computeResidual(sol,yVal,hystSol,mu,mu_inv,wrtX,relative);
      
      hystNew = computeValue_vec(xNew, operatorIdx, overwriteMemory, debugOut, successFlag);
      resNew = computeResidual(xNew,yVal,hystNew,mu_inv);
      //      resNew = computeResidual(xNew,yVal,hystNew,mu,mu_inv,wrtX,relative);
      
//      LOG_DBG(vecpreisachlinesearch) << "hyst vector for sol: " << hystSol.ToString();
//      LOG_DBG(vecpreisachlinesearch) << "hystNew: " << hystNew.ToString();
//      LOG_DBG(vecpreisachlinesearch) << "hystOld: " << hystOld.ToString();
//      LOG_DBG(vecpreisachlinesearch) << "Old res vector: " << res.ToString();
//      LOG_DBG(vecpreisachlinesearch) << "New res vector: " << resNew.ToString();
//      LOG_DBG(vecpreisachlinesearch) << "res vector for sol: " << resSol.ToString();
      // set stayBelowSat flag to 0 to disable checking
      
      Vector<Double> jac_dx = Vector<Double>(dim_);
      jac.Mult(xUpdate,jac_dx);
      success = checkIncrementTrustRegion(xNew, res, resNew, jac_dx, alpha, stayBelowSat); 
      
      //      success = checkIncrement(xNew, xUpdate, res, resNew, jac, alpha,0);
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
          hystNew = computeValue_vec(xNew, operatorIdx, overwriteMemory, debugOut, successFlag);
          resNew = computeResidual(xVal,yVal,hystNew,mu_inv);
          //          resNew = computeResidual(xNew,yVal,hystNew,mu,mu_inv,wrtX,relative);
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
//          std::cout << "Alpha max reached" << std::endl;
					LOG_TRACE(vecpreisachlinesearch) << "Alpha max reached; take update > stop";
					alpha = alphaMax;
					discard = true;
					break;
				}
				if(alpha < alphaMinLocal){
//          std::cout << "Alpha min reached" << std::endl;
					LOG_TRACE(vecpreisachlinesearch) << "Alpha min reached; discard update > stop";
					alpha = alphaMinLocal;
					discard = true;
					break;
				}
				
        if(itCnt >= maxIter){
//          std::cout << "Max number of iterations used" << std::endl;
          LOG_TRACE(vecpreisachlinesearch) << "Linesearch was not successful; Discard update.";
          discard = true;
          break;
        }
      }
    }
    
//    std::cout << "alpha: " << alpha << std::endl;
//    std::cout << "alphaMax: " << alphaMax << std::endl;
//    std::cout << "alphaAcc: " << std::sqrt(alphaAcc) << std::endl;
    
    return discard;
  }
  
	Vector<Double> Hysteresis::computeInput_vec_withPrevStates(Vector<Double> yVal, Vector<Double> prevYval,
          Vector<Double> prevXval, Vector<Double> prevHystval, Integer operatorIndex, 
          Matrix<Double> mu, bool fieldsAlignedAboveSat, 
          bool hystOutputRestrictedToSat, int& successFlag){
		
		UInt totalNumberOfLMIterations=0;
		UInt totalNumberOfLinesearchIterations=0;
		UInt maximalNumberOfLinesearchIterations=0;
    Double minAlphaStatistics,maxAlphaStatistics,avgAlphaStatistics;
		Vector<Double>sol = Vector<Double>(dim_);
    
		return computeInput_vec_withStatistics(yVal, prevYval, prevXval, prevHystval, 
            operatorIndex, mu, fieldsAlignedAboveSat, hystOutputRestrictedToSat,
            totalNumberOfLMIterations, totalNumberOfLinesearchIterations, 
            maximalNumberOfLinesearchIterations, successFlag, minAlphaStatistics,maxAlphaStatistics,avgAlphaStatistics,sol);
	}
  
  Vector<Double> Hysteresis::computeInput_vec_withStatistics(Vector<Double> yVal, Vector<Double> prevYval,
          Vector<Double> prevXval, Vector<Double> prevHystval, Integer operatorIndex, 
          Matrix<Double> mu, bool fieldsAlignedAboveSat, bool hystOutputRestrictedToSat,
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
    // note 24.05.2018: the treatment as described in note from 22.05. helps a lot; nevertheless some inversion fails
    //                  most probably as we further down (and during subfunctions) also assume that the maximal output is PSaturated
    //                  > two options possible
    //                    a) search all these spots and replace PSaturated+anhystPartPosSat with hystValAtXSat.NormL2()
    //                    b) add scaling to vector model, such that PSaturated is actually reached at saturation
    //                  > currently option b) is used!
    // note 9.7.2018: if option b) is used, the parameter rotRes has a much smaller influence as expected; 
    //                  the actual determination of rotRes from measurement data is harder compared to allowing
    //                  the rotation state to change AFTER saturation by first computing rotRes*inputNorm and then
    //                  cutting it down to 1;
    //                this treatment has the severe consequence, that input and output might no longer be aligned at
    //                saturation because parts of the rotational state get overwritten by values beyond saturation
    //                > sutor model must be treated similar to Mayergoyz model, i.e. fieldsAlignedAboveSat = false
    //                > but: sutor model will never exceed saturtion, i.e. hystOutputRestrictedToSat = false
    Vector<Double> satInput = Vector<Double>(dim_);
    satInput.Init();
    satInput.Add(XSaturated_,yDir);
    int successCodeTmp = 0;
    Vector<Double> hystValAtXSat = computeValue_vec(satInput, operatorIndex, false, false, successCodeTmp);
    Double anhystPartPosSat = PSaturated_*evalAnhystPart_normalized(1.0);
    //      std::cout << "yNorm: " << yNorm << std::endl;
    //      std::cout << "hystValAtXSat: " << hystValAtXSat << std::endl;
    //      std::cout << "XSaturated_: " << XSaturated_ << std::endl;
    //      std::cout << "PSaturated_: " << PSaturated_ << std::endl;
    //      std::cout << "eps_mu: " << eps_mu << std::endl;
    //      
    // check saturation in direction of yIn might solve system
    //Double diffSat = yNorm - (PSaturated_ + eps_mu*XSaturated_ + anhystPartPosSat);
    Double diffSat = yNorm - (hystValAtXSat.NormL2() + eps_mu*XSaturated_);
    
    LOG_TRACE(vecpreisachInversion) << "yNorm: " << yNorm;
    if(yNorm == 0){
      stayBelowSat = 1;
    } else {   
      if(abs(diffSat) < INV_resTolB_){
        traceMsg << "--B Special-- Inversion: Exact Saturation found" << std::endl;
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
          //Xup = (yNorm - PSaturated_)/eps_mu;
          //Poffset = PSaturated_;
          Xup = (yNorm - hystValAtXSat.NormL2())/eps_mu;
          Xdown = XSaturated_;
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
     * Check for remanence
     * > apparently, LM seems to work fine except for the remanence case
     * i.e. yval = hystval, xval = 0
     * > check if xVal = 0 would be a proper solution
     */
    Vector<Double> hystVal_rem;
    xTMP.Init();
    bool debugOut = false;
    hystVal_rem = computeValue_vec(xTMP, operatorIndex, overwriteMemory, debugOut, successFlagForward);
    
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
      return xVal;
      //	compCase = 2;
    } 
    
    Double diffRem = diff.NormL2();
    
    /*
     * Solve nonlinear problem with help of Newton-based scheme
     * 
     * residual/target function: 
     *    F(xVal) = xVal + mu_inv( hystVal(xVal) - yVal )
     * 
     * Jacobian of residual:
     *    Jac(xVal)_ji = delta_ji + mu_inv[j][:]*dhystVal/dxVal_i
     * 
     * Newton:
     *    F(xNew) = F(xVal) + Jac(xVal) deltaX + O(deltaX^2)
     * 
     *  > solve:
     *      Jac(xVal) deltaX = -F(xVal)
     *  > perform linesearch to find s in [0,1]
     *  > update:
     *      xNew = xVal + s deltaX
     * 
     * >> pure Newton not working as Jacobian cannot be computed exactly
     * >> approximated Jacobian might cause convergence issues
     * 
     * Jacobian-Free-Newton-Krylov:
     *  > as Newton, but update constructed from Krylov subspace
     *  
     * -- solve step --   
     *  > take initial guess deltaX0
     *  > compute inital r0 = -F(xVal) - Jac(xVal) deltaX0
     *  > cnt = 1
     *  > while()
     *    > deltaX = deltaX0 + sum_i=0^cnt-1 beta_i [Jac(xVal)]^i r0
     *    > deltaX converged > stop
     *    > deltaX not converged > continue
     * ----------------
     *  > perform linesearch to find s in [0,1]
     *  > update:
     *      xNew = xVal + s deltaX
     * 
     * >> big advantage: Jac is not required directly but only Jac*someVector
     * >> alternative:
     *    move linesearch to computation of deltaX to determine beta_i
     * 
     * Levenberg-Marquardt:
     *  > try to compensate issues with Jacobian by regularization
     *  > no linesearch, instead optimize alpha
     * 
     *  > while()
     *    > solve:
     *       (Jac(xVal)^T Jac(xVal) + alpha^2 I) deltaX = -Jac(xVal)^T F(xVal)
     *    > check deltaX by trust region method
     *    > deltaX ok > stop
     *    > deltaX not ok > adept alpha
     *  > update:
     *      xNew = xVal + deltaX
     * 
     * >> works in many cases but still has issues due to approximation of Jac
     * >> adaption of alpha might require various iterations
     */
    
    /*
     * obtain start value
     */
    Vector<Double> xStart = obtainStartVector(prevXval, prevYval, yVal, diffRem, diffSat, stayBelowSat);
    
    traceMsg << "--D-- Inversion: Use Levenberg Marquart" << std::endl;
    
    // tolerance wrt y > 1e-10 or 1e-12 seems good > takes 2-3 its
    // only problem: y-x-loops look ugly as x can be quite off!
    //Double tolError = 1e-11;  
    // tolerance for reevalution
    
    UInt itCnt = 0;
    Double alpha = INV_alphaLSStart_;
    Double alphaMin = INV_alphaLSMin_;//1.0/256.0;
    Double alphaMax = INV_alphaLSMax_;//8192;//512.0;
    
    Vector<Double> hystVal = computeValue_vec(xStart, operatorIndex, overwriteMemory, debugOut, successFlagForward);
    xVal = xStart;
    
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
    //    std::cout << "stayBelowSat: " << stayBelowSat << std::endl;
    Vector<Double> xVal_krylov = Vector<Double>(dim_);
    xVal_krylov = xVal;
    while(true){ 
      itCnt++;
      totalNumberOfLMIterations++;
      //      std::cout << "itCnt: " << itCnt << std::endl;
      if(debug){
        traceMsg << "--OUTER ITERATION-- (" << itCnt << ")" << std::endl;
        traceMsg << "Current solution: " << xVal.ToString() << std::endl;
      }
      // do not override hyst memory here
      //          ClipDirection(xVal);
      hystVal = computeValue_vec(xVal, operatorIndex, overwriteMemory, debugOut, successFlagForward);
      
      res = computeResidual(xVal,yVal,hystVal,mu_inv);
      jac = computeJacobian(xVal,hystVal, mu_inv, 
              operatorIndex, sign, implementation, overwriteMemory, stayBelowSat);
      
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
              operatorIndex, overwriteMemory);
      
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
            successFlag = -1;
            
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
      
      discardUpdate = computeUpdate_LM_full(xVal, yVal, res, xUpdate, jac, jacT, mu, mu_inv, 
              operatorIndex, overwriteMemory, alpha, alphaMin, alphaMax, 
              numberOfIterations,xStart,factorToSat,stayBelowSat,sol);
      
      Vector<Double> xUpdate_Krylov = computeUpdate_Krylov(xVal, yVal, hystVal, mu_inv, operatorIndex);
      Vector<Double> xUpdate_Newton = computeUpdate_Newton(xVal, yVal, hystVal, mu_inv, operatorIndex);
      //      std::cout << "-- xUpdate - LM: " << xUpdate.ToString() << std::endl;
      //      std::cout << "-- xUpdate - Krylov: " << xUpdate_Krylov.ToString() << std::endl;
      //      std::cout << "-- xUpdate - Newton: " << xUpdate_Newton.ToString() << std::endl;
      xVal_krylov.Add(xUpdate_Krylov);
      
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
      
      LOG_DBG(vecpreisachInversion) << "Computed update: " << xUpdate.ToString();
      if(!discardUpdate){
        xVal = xVal+xUpdate;
      } else {
        LOG_DBG(vecpreisachInversion) << "Discard update";
        //std::cout << "Discard update; reset to best solution so far" << std::endl;
        alphaMin = alphaMin/2.0;
        alphaMax = alphaMax*2.0;
      }
      sign = sign*(-1.0);
    }
    //    std::cout << "xVal - LM: " << xVal.ToString() << std::endl;
    //    std::cout << "xVal - Krylov: " << xVal_krylov.ToString() << std::endl;
    
    
    //      std::cout << "xVal.NormL2(): " << xVal.NormL2() << std::endl;
    //      std::cout << "XSaturated_: " << XSaturated_ << std::endl;
    //      
    avgAlphaStatistics /= totalNumberOfLMIterations;
    
    if((xVal.NormL2() > XSaturated_)&&(stayBelowSat==1)){
      EXCEPTION("LM lead xVal into saturation although input is below > must not be the case!");
    }
    
    if((xVal.NormL2() < XSaturated_)&&(stayBelowSat==-1)){
      EXCEPTION("LM lead xVal below saturation although input is above> must not be the case!");
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
      checkInversionOutput(xVal, yVal, mu, INV_resTolB_, resYNorm, operatorIndex, overwriteMemory,true);
      LOG_TRACE(vecpreisachInversion) << traceMsg.str();
    }
    
    LOG_TRACE(vecpreisachInversion) << " --------- END IVERSION --------- ";
    
    return xVal;
  } 
  
	// returns true, if yTarget can be reached by mu*xComputed + hyst(xComputed)
	bool Hysteresis::checkInversionOutput(Vector<Double>& xComputed, Vector<Double>& yTarget, 
          Matrix<Double>& mu, Double tol, Double& resYNorm, Integer operatorIdx, bool overwriteMemory, bool output){
		
		Vector<Double> yCheck = Vector<Double>(dim_);
    bool debugOut = false;
    int successCode = 0;
    Vector<Double> hCheck = computeValue_vec(xComputed, operatorIdx, overwriteMemory, debugOut, successCode);
    
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
