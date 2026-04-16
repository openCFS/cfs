#include <iostream>
#include <fstream>
#include <iomanip>

#include "Hysteresis.hh"

#include "Utils/ToolsFull.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

namespace CoupledField
{ 
  DEFINE_LOG(hysteresis_inversion, "hysteresis_inversion")
  DEFINE_LOG(hysteresis_linesearch, "hysteresis_linesearch")

  /*
   * Base routines - Constructors, Destructor
   */
  Hysteresis::Hysteresis(Integer numElem)
  {
    Hysteresis(numElem, 1.0,1.0,1.0,0.0,0.0,0.0, false,true);
  }

  Hysteresis::Hysteresis(Integer numElem, Double XSaturated, Double YSaturated, Double hystSaturated,
          Double anhystA, Double anhystB, Double anhystC, Double anhystD, bool anhystOnly){

    EXCEPTION("Hysteresis - old type constructor - Should not be used anymore");
    numElements_ = numElem;
    anhyst_A_ = anhystA;
    anhyst_B_ = anhystB;
    anhyst_C_ = anhystC;
    anhyst_D_ = anhystD;
    anhystOnly_ = anhystOnly;

    XSaturated_ = XSaturated;
    PSaturated_ = YSaturated;
    hystSaturated_ = hystSaturated;
    dim_ = 0;
    invParamsSet_ = false;
  }

  Hysteresis::Hysteresis(Integer numElem, ParameterPreisachOperators operatorParams,
          ParameterPreisachWeights weightParams)
  {
    numElements_ = numElem;
    anhyst_A_ = weightParams.anhysteretic_a_;
    anhyst_B_ = weightParams.anhysteretic_b_;
    anhyst_C_ = weightParams.anhysteretic_c_;
    anhyst_D_ = weightParams.anhysteretic_d_;
    anhystOnly_ = weightParams.anhystOnly_;

    XSaturated_ = operatorParams.inputSat_;
    XForAlignment_ = operatorParams.inputForAlignment_;
//    std::cout << "3. and 5. - HysteresisConstructor - New parameter XForAlignment_ = " << XForAlignment_ << std::endl;
    PSaturated_ = operatorParams.outputSat_;
    hystSaturated_ = operatorParams.preisachSat_;
    checkInversionResult_ = operatorParams.checkInversionResult_;
    printWarnings_ = operatorParams.printWarnings_;
    dim_ = 0;
    invParamsSet_ = false;
  }

  Hysteresis::~Hysteresis()
  {
  }

	Double Hysteresis::bisectForAnhyst(Double Ytarget,
          Double Xdown, Double Xup, Double Poffset, Double eps_mu, Double tol, Vector<Double> dir, UInt idx, int& successFlag){

//    std::cout << "Biscet before normalization" << std::endl;
//      std::cout << " Ytarget: " << Ytarget << std::endl;
//      std::cout << " Xdown: " << Xdown << std::endl;
//      std::cout << " Xup: " << Xup << std::endl;
//      std::cout << " Poffset: " << Poffset << std::endl;
//      std::cout << " XSaturated_: " << XSaturated_ << std::endl;
//      std::cout << " PSaturated_: " << PSaturated_ << std::endl;
//      std::cout << " hystSaturated_: " << hystSaturated_ << std::endl;

		return XSaturated_*bisectForAnhyst_normalized(Ytarget/PSaturated_,
            Xdown/XSaturated_, Xup/XSaturated_, Poffset/PSaturated_, XSaturated_*eps_mu/PSaturated_, tol/XSaturated_, dir, idx, successFlag);

	}

	Double Hysteresis::bisectForAnhyst_normalized(Double Ytarget_normalized,
          Double Xdown_normalized, Double Xup_normalized, Double Poffset_normalized,
          Double eps_mu_normalized, Double tol, Vector<Double> dir, UInt idx, int& successFlag){

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
    /*
     * 8.6.2020
     *  Success flag added
     *    N = success after itartion N
     *    0 = success due to close enough initial guess
     *   -1 = fail
     */

		Double Xout_normalized, Xmid_normalized;
		Double resUp, resDown, resMid;

    bool evalRequired = false;
    Vector<Double> Pout = Vector<Double>(dim_);
    Vector<Double> Pin = Vector<Double>(dim_);
    if(dir.NormL2() != 0){
      evalRequired = true;
    }
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
      std::stringstream errormsg;
      UInt precisionDigits = 12;
      errormsg << "Bisection for anhystPart -> Starting values not appropriate! Both residuals have same sign!" << std::endl;
      errormsg << " XSaturated_: " << XSaturated_ << std::endl;
      errormsg << " PSaturated_: " << PSaturated_ << std::endl;
      errormsg << " eps_mu_normalized: " << eps_mu_normalized << std::endl;
      errormsg << " Ytarget_normalized: " << Ytarget_normalized << std::endl;
      errormsg << std::setprecision(precisionDigits) << " Xdown_normalized: " << Xdown_normalized << std::endl;
      errormsg << std::setprecision(precisionDigits) << " Xup_normalized: " << Xup_normalized << std::endl;
      errormsg << std::setprecision(precisionDigits) << " Xdown_normalized-Xup_normalized: " << Xdown_normalized-Xup_normalized << std::endl;
      errormsg << " Poffset_normalized: " << Poffset_normalized << std::endl;
      errormsg << " resUp: " << resUp << std::endl;
      errormsg << " resDown: " << resDown << std::endl;
      errormsg << " evalRequired? " << evalRequired << std::endl;
      EXCEPTION(errormsg.str());
    }
    if(abs(resUp) < tol){
      Xout_normalized = Xup_normalized;
      successFlag = 0;
    } else if(abs(resDown) < tol){
      Xout_normalized = Xdown_normalized;
      successFlag = 0;
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
          successFlag = i+1;
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
      successFlag = -1;
      //      LOG_DBG(scalpreisachInversion) << "bisection failed; reamining residual: " << resMid;
    }
    return Xout_normalized;

	}

  bool Hysteresis::checkConvergence(Vector<Double>& jacTres, Double& errorNorm, Double tol){
    // According to Dahmen&Reusken - Numerik partieller DFG
    // the residual is a non-sufficient condition
    // instead we have to check the norm of jacT*res
    errorNorm = jacTres.NormL2();

    if(errorNorm <= tol){
      return true;
    } else {
      return false;
    }

  }

  Integer Hysteresis::checkIncrementTrustRegion(Vector<Double>& x_new,
          Vector<Double>& res_cur, Vector<Double>& res_new,
          Vector<Double>& jac_dx, Double& alpha, Double alphaStepUp, Double alphaStepDown, int stayBelowSat){
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
    beta0 = INV_params_.trustLow;
    betaLow = INV_params_.trustMid;
    betaUp = INV_params_.trustHigh;

    //    beta0 = 0.2;
    //    betaLow = 0.4;
    //    betaUp = 0.8;

    Double factorUp = alphaStepUp;
    Double factorDown = alphaStepDown;

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

      LOG_TRACE(hysteresis_inversion) << "Nominator: ||F(x)||^2 - ||F(x + increment)||^2 = " << nominator;
      LOG_TRACE(hysteresis_inversion) << "Denominator: ||F(x||^2 - ||F(x) + F'(x)*increment||^2 = " << denominator;

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
    
    if (YdiffToRemancence < INV_params_.tolB){
      // we are quite close to remanence, but have not yet reached it
      // start from 0
      startXVector.Init();
      //  startAtZero=true;
      //	compCase = 3;
      //    } else if (YdiffToSaturation < 100*INV_params_.tolB){
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
	  /*
	  * DEBUGGING 17.6.2020: the first case cannot occur as startXVector is a fresh vector at this point
	  * it is up to the numerical precision if we end up here!
	  *
	  * startXVector = previousXVector; should appear BEFORE the if
	  * inside if we should not check for != 0 but for < tolerance!
	  *
	  */
      // 19.6.2020 move line out of if clause
      // > works apparently better than starting from zero-vector!
      startXVector = previousXVector;
      if(startXVector.NormL2() != 0){
//        std::cout << "startXVector.NormL2() "<<startXVector.NormL2()<<" != 0" << std::endl;
//        startXVector = previousXVector;

        if(previousYVector.NormL2() != 0){
          factor = currentYVector.NormL2()/previousYVector.NormL2();
        } else {
          if(currentYVector.NormL2() > previousYVector.NormL2()){
            factor = 1.25;
          } else {
			  // this point can only be reached if previousYVector.NormL2() AND currentYVector.NormL2() = 0!
            factor = 1.0/1.25;
          }
        }
        startXVector.ScalarMult(factor);
      } else {
//        std::cout << "startXVector.NormL2() "<<startXVector.NormL2()<<" == 0" << std::endl;
		  // currently we always end up here if the constructor Vector<Double>(size) initializes the vector really with 0! 
        //            traceMsg << "X is zero; scaling will not help; try a scaled version of yVal instead" << std::endl;
        startXVector = currentYVector;
		// this can be problematic as currentYVector contains the anhysteretic part whereas PSaturated does not!
		// shouldn't we scale by XSaturated_/YSaturated_ instead?
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

//    std::cout << "startXVector: " << startXVector.ToString(6) << std::endl;
    //    std::cout << "currentYVector.NormL2(): " << currentYVector.NormL2() << std::endl;
    //    std::cout << "dir: " << dir.ToString() << std::endl;
    //    std::cout << "XSaturated_: " << XSaturated_ << std::endl;
    if((startXVector.NormL2() >= XSaturated_)&&(stayBelowSat==1)){
      //        traceMsg << "Reset xVal as its value is above Xsaturation but yVal is not" << std::endl;
      startXVector = dir;
      startXVector.ScalarMult(1.0/1.5*XSaturated_);
    }

    if((startXVector.NormL2() < XSaturated_)&&(stayBelowSat==-1)){
      //        traceMsg << "Reset xVal as its value is below Xsaturation but yVal is above" << std::endl;
      startXVector = dir;
      startXVector.ScalarMult(1.5*XSaturated_);
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
          Vector<Double>& hyst_x, Matrix<Double> mu_inv, Integer operatorIdx, Double scalingForJacDiagonal){

    /*
     * Jac*dX = -res
     * > solve directly
     */
    Matrix<Double> jac = computeJacobian(x, hyst_x,
            mu_inv, operatorIdx, 1.0, false, 0, scalingForJacDiagonal);
//    computeJacobian(Vector<Double>& xVal, Vector<Double>& hystVal,
//          Matrix<Double> mu_inv, Integer operatorIdx, Double sign,
//					bool overwriteMemory, int stayBelowSat, Double scalingForJacDiagonal) 
//    std::cout << "scalingForJacDiagonal: " << scalingForJacDiagonal << std::endl;
//    std::cout << "Computed Jacobian for Newton = " << jac.ToString() << std::endl;
    
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

  bool Hysteresis::CompareJacobianApproximations(Vector<Double>& xVal, Vector<Double>& y,
          Matrix<Double> mu_inv, Integer operatorIdx){

    int successFlag = 0;
    Vector<Double> hystVal = computeValue_vec(xVal, operatorIdx, false, false, successFlag);

    Matrix<Double> FDJac_Kelly_v1 = ApproximateFDJacobian(xVal, y, hystVal, mu_inv, operatorIdx, 1);
    Matrix<Double> FDJac_Kelly_v2 = ApproximateFDJacobian(xVal, y, hystVal, mu_inv, operatorIdx, 2);
    Matrix<Double> FDJac_Old = computeJacobian(xVal, hystVal, mu_inv, operatorIdx, 1.0, false, 0, 1.0);

    std::cout << "FDJac_Kelly_v1 = " << FDJac_Kelly_v1.ToString() << std::endl;
    std::cout << "FDJac_Kelly_v2 = " << FDJac_Kelly_v2.ToString() << std::endl;
    std::cout << "FDJac_Old = " << FDJac_Old.ToString() << std::endl;


    Vector<Double> v1 = Vector<Double>(dim_);
    v1[0] = 1.0;
    v1[1] = -934;

    Vector<Double> JacV1_Kelly = ApproximateDirectionalDerivative(xVal, v1, y, hystVal, mu_inv, operatorIdx);
    Vector<Double> JacV1_Old = computeJacobianTimesVector(xVal, v1, y, hystVal, mu_inv, operatorIdx);

    std::cout << "JacV1_Kelly = " << JacV1_Kelly.ToString() << std::endl;
    std::cout << "JacV1_Old = " << JacV1_Old.ToString() << std::endl;

    Vector<Double> v2 = Vector<Double>(dim_);
    v2[0] = 1.0;
    v2[1] = -0.934;

    Vector<Double> JacV2_Kelly = ApproximateDirectionalDerivative(xVal, v2, y, hystVal, mu_inv, operatorIdx);
    Vector<Double> JacV2_Old = computeJacobianTimesVector(xVal, v2, y, hystVal, mu_inv, operatorIdx);

    std::cout << "JacV2_Kelly = " << JacV2_Kelly.ToString() << std::endl;
    std::cout << "JacV2_Old = " << JacV2_Old.ToString() << std::endl;

    return true;
  }

  Matrix<Double> Hysteresis::ApproximateFDJacobian(Vector<Double>& x, Vector<Double>& y, Vector<Double>& hyst_x,
          Matrix<Double> mu_inv, Integer operatorIdx, int version){
    /*
     * Approximate Jacobian of residual at point x using Finite Differences
     *
     * Source: "Iterative Methods for Linear and Nonlinear Equations" - Kelly 1995
     * p.80
     */
    Double h = 1e-7; // approx. sqrt(double precision)
    UInt vecSize = x.GetSize();

    Matrix<Double> jac = Matrix<Double>(vecSize,vecSize);
    jac.Init();

    Vector<Double> resAtX = computeResidual(x,y,hyst_x,mu_inv);

    Vector<Double> xShifted = Vector<Double>(vecSize);
    Vector<Double> hystAtXshifted = Vector<Double>(vecSize);
    Vector<Double> resAtXshifted = Vector<Double>(vecSize);

    Double scaling = h;
    if(x.NormL2() != 0){
      scaling *= x.NormL2();
    }

    bool overwriteMem = false;
    bool debugOut = false;
    int successFlag = 0;

    if(version == 1){
      // follow approach of Kelly directly and use residual for computation of Jacobian
      for(UInt col = 0; col < vecSize; col++){
        xShifted.Init();
        xShifted.Add(1.0,x);
        xShifted[col] += scaling;

        hystAtXshifted = computeValue_vec(xShifted, operatorIdx, overwriteMem, debugOut, successFlag);
        resAtXshifted = computeResidual(xShifted,y,hystAtXshifted,mu_inv);

        resAtXshifted.Add(-1.0,resAtX);
        resAtXshifted.ScalarDiv(scaling);

        for(UInt row = 0; row < vecSize; row++){
          jac[row][col] = resAtXshifted[row];
        }
      }
    } else if(version == 2){
      // use old approach and compute only the Jacobian of the hyst operator; then add identity
      // note: res = xVal + hystVal(xVal)/mu - yVal/mu
      //        > dres/dx_j = delta_ij + mu_inv*dhystVal(xVal)/dx_j
      Vector<Double> column = Vector<Double>(vecSize);
      for(UInt col = 0; col < vecSize; col++){
        xShifted.Init();
        xShifted.Add(1.0,x);
        xShifted[col] += scaling;

        hystAtXshifted = computeValue_vec(xShifted, operatorIdx, overwriteMem, debugOut, successFlag);

        hystAtXshifted.Add(-1.0,hyst_x);
        hystAtXshifted.ScalarDiv(scaling);
        mu_inv.Mult(hystAtXshifted,column);

        for(UInt row = 0; row < vecSize; row++){
          jac[row][col] = column[row];
        }
        jac[col][col] += 1.0;
      }
    }
    return jac;
  }


  Vector<Double> Hysteresis::ApproximateDirectionalDerivative(Vector<Double>& x, Vector<Double>& v,
          Vector<Double>& y, Vector<Double>& hyst_x, Matrix<Double> mu_inv, Integer operatorIdx){

    /*
     * Approximate directional derivative at point x into direction v
     * <=> approximate product of Jacobian(x) * v
     *
     * Note: Jacobian is meant with respect to the residual
     *
     * Source: "Iterative Methods for Linear and Nonlinear Equations" - Kelly 1995
     * p.81 (5.15)
     */
    Double h = 1e-7; // approx. sqrt(double precision)
    UInt vecSize = v.GetSize();

    assert(vecSize == x.GetSize());

    Vector<Double> directionalDerivative = Vector<Double>(vecSize);
    directionalDerivative.Init();

    if(v.NormL2() == 0){
      return directionalDerivative;
    }

    Vector<Double> resAtX = computeResidual(x,y,hyst_x,mu_inv);

    Vector<Double> xShifted = Vector<Double>(vecSize);
    Vector<Double> hystAtXshifted = Vector<Double>(vecSize);
    Vector<Double> resAtXshifted = Vector<Double>(vecSize);
    Double scaling = v.NormL2()/h;

    // xShifted = h*v/|v|
    xShifted.Add(h/v.NormL2(),v);

    if(x.NormL2() != 0){
      // xShifted = x + h*|x|*v/|v|
      xShifted.ScalarMult(x.NormL2());
      xShifted.Add(1.0,x);
      // scaling = |v|/(h*|x|)
      scaling /= x.NormL2();
    }

    bool overwriteMem = false;
    bool debugOut = false;
    int successFlag = 0;
    hystAtXshifted = computeValue_vec(xShifted, operatorIdx, overwriteMem, debugOut, successFlag);
    resAtXshifted = computeResidual(xShifted,y,hystAtXshifted,mu_inv);

    directionalDerivative.Add(1.0,resAtXshifted,-1.0,resAtX);
    directionalDerivative.ScalarMult(scaling);

    return directionalDerivative;
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
              mu_inv, operatorIdx, 1.0, 0, false, 1.0);

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
          Matrix<Double> mu_inv, Integer operatorIdx, Double sign,
					bool overwriteMemory, int stayBelowSat, Double scalingForJacDiagonal) {

    LOG_DBG(hysteresis_inversion) << " --------- computeJacobian --------- ";
    //    LOG_DBG(vecpreisach) << "VecPreisach::computeJacobian";
    if((xVal.NormL2() > XSaturated_)&&(stayBelowSat==1)){
      EXCEPTION("xVal.NormL2() > XSaturated_");
    }
    if((xVal.NormL2() < XSaturated_)&&(stayBelowSat==-1)){
      EXCEPTION("xVal.NormL2() < XSaturated_");
    }

    Double deltaX = 0.0;
		Double scal = INV_params_.jacRes;
    Double deltaXmin = scal; // XSaturated_; 
    //  19.06.2020 - Ongoing issue: Inversion via Newton fails bad at some points as Jacobian does not fit at all
    // Question:
    // > is the used implementation, especially the resolution, suitable at all?
    // in coef function hyst (where Jacobian of material relation is computed, i.e., of B = muH + Jm)
    // the stepping length is set to 
    //     steppingDistance = steppingSign*std::max(scaling,scaling*E_H.NormL2());
    // where scaling = 1e-7 (=default value of INV_params_.jacRes!);
    // here, we set 
    //        if( xVal[i] < 0 ){
    //          deltaX = sign*std::min( scal*xVal[i], -deltaXmin );
    //        } else {
    //          deltaX = sign*std::max( scal*xVal[i], deltaXmin );
    //        }
    // i.e. scal*xVal[i] is comparable to scaling*E_H.NormL2()
    // but deltaXmin = scal*XSaturated_ is way larger than scal!
    // > idea:
    //  set deltaX = sign(std::max( scal*xVal.NormL2(), scal)) 
    //  and check if this works better!
    // > note: error was found elsewhere! mayergoyz model had issues that had nothing to do with 
    //  this routine; nevertheless, set deltaXmin to INV_params_.jacRes;

    //		Double deltaXmin = 1e-8*XSaturated_;
    //    Double scal = 1e-8;//1e-5;
    //    bool overwriteMemory = false;
    Vector<Double> xShifted;
    Vector<Double> hystShifted;
    Vector<Double> deltaHyst;
    Matrix<Double> jac = Matrix<Double>(dim_,dim_);
    jac.Init();

    if(INV_params_.jacImplementation == 2){
      /*
       * Full Jacobian using forward/backward differences
       * BUT different approach on scaling factor
       * > taken from "Jacobian-free Newton-Krylov methods: a survey ..." - Knoll
       * > J_ij = ( F_i(u + eps_j*e_j) - F_i(u) )/eps_j
       * > eps_j = b*u_j + b with b approx sqrt(machine precision)
       */
      Double b = std::sqrt(INV_params_.jacRes);
      Double eps;
      for(UInt j = 0; j < dim_; j++){
        eps = sign*(b*abs(xVal[j]) + b);

        xShifted = xVal;

        xShifted[j] += eps;

        int successFlag = 0;
        bool debugOut = false;
        hystShifted = computeValue_vec(xShifted, operatorIdx, overwriteMemory, debugOut, successFlag);
        /*
         * Compute Jacobian for residual wrt x
         *
         * jac_ij = + delta_ij + mu_inv[i][:]*dhystVal/dxVal_j
         */
        deltaHyst = hystShifted;
        deltaHyst.Add(-1.0,hystVal);

        Vector<Double> curCol = Vector<Double>(dim_);
        mu_inv.Mult(deltaHyst,curCol);
        curCol.ScalarDiv(eps);

        // should be 1 (due to eps*1/eps, but apparently a higher value helps sometimes
        // idea: try to solve with 1, if this does not work increase step by step
        // > some sort of regularization
        // > DEPRECATED! Resulting H will not fit to input B!
        jac[j][j] = scalingForJacDiagonal;
        for(UInt i = 0; i < dim_; i++){
          jac[i][j] += curCol[i];
        }
      }

      // for testing of jacobianTimesVector, we can compute the Jacobian and multiply it by the vector
      // and compare; use std forward/backward differences for comparison
      // note: INV_params_.jacImplementation = -1 if jacobianTimesVector should be used
    } else if(INV_params_.jacImplementation <= 0){
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
        LOG_DBG(hysteresis_inversion) << " xVal " << xVal.ToString();
        LOG_DBG(hysteresis_inversion) << " xShifted " << xShifted.ToString();
        LOG_TRACE(hysteresis_inversion) << " dXvec " << dXvec.ToString();
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

        LOG_DBG(hysteresis_inversion) << " hystVal " << hystVal.ToString();
        LOG_DBG(hysteresis_inversion) << " hystShifted " << hystShifted.ToString();
        LOG_TRACE(hysteresis_inversion) << " deltaHyst " << deltaHyst.ToString();

        Vector<Double> curCol = Vector<Double>(dim_);
        mu_inv.Mult(deltaHyst,curCol);
        curCol.ScalarDiv(deltaX);

        jac[i][i] = 1.0;
        for(UInt j = 0; j < dim_; j++){
          jac[j][i] += curCol[j];
        }
      }      
    } else if(INV_params_.jacImplementation == 1){
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

  Vector<Double> Hysteresis::computeUpdate_LM_withFixedAlpha(Vector<Double> jacTres_neg,
          Matrix<Double>& jacTjac, Vector<Double>& res, Double alphaFix){

    /*
     * Compute update dx by solving
     *
     *  [Jac^T(x_k)*Jac(x_k) + alpha_k*Identity] * dx = -Jac^T(x_k)*res(x_k)
     *
     * with alpha_k = alpha*||res(x_k)||
     *
     */
    Matrix<Double> matToInvert = Matrix<Double>(dim_,dim_);
    Matrix<Double> matInverted = Matrix<Double>(dim_,dim_);
    Vector<Double> dx = Vector<Double>(dim_);

    Double alpha = alphaFix*res.NormL2();
    Double minDeterminant = 1e-12;
    Double detMatToInvert;

    matToInvert = jacTjac;
    /*
     * Add regularization
     */
    for(UInt i = 0; i < dim_; i++){
      matToInvert[i][i] += alpha;
    }

    matToInvert.Determinant(detMatToInvert);

    if(abs(detMatToInvert)<minDeterminant){
      EXCEPTION("LM with fixed alpha: Alpha*||res|| not sufficient to make matrix invertible!");
    }

    matToInvert.Invert(matInverted);
    matInverted.Mult(jacTres_neg,dx);

    assert(!dx.ContainsNaN() && !dx.ContainsInf());

    return dx;
  }

  Vector<Double> Hysteresis::computeUpdate_LM(Vector<Double> jacTres_neg,
          Matrix<Double>& jacTjac, Double& alpha, Double& alphaAcc, Double& alphaMinReg, Double alphaMaxReg){

    /*
     * Compute update s_k by solving
     *
     *  [Jac^T(x_k)*Jac(x_k) + alpha_k*Identity] * s_k = -Jac^T(x_k)*res(x_k)
     *
     * with alpha_k in [alphaMinReg, alphaMaxReg]
     * alpha_k is iteratively increased until [Jac^T(x_k)*Jac(x_k) + alpha_k*Identity]
     * becomes invertible
     *
     */
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
        if( (cnt > 200)||(alpha > alphaMaxReg) ){
          EXCEPTION("LM: Cannot find alpha, such that jacTjac becomes invertible!");
        }
      }
    }

    alphaMinReg = alpha;
    matToInvert.Invert(matInverted);
    matInverted.Mult(jacTres_neg,dx);

    assert(!dx.ContainsNaN() && !dx.ContainsInf());

    return dx;
  }

  bool Hysteresis::computeUpdateNewton_Linesearch(Vector<Double>& xStart, Vector<Double>& xCurrent, Vector<Double>& xUpdate,
          Vector<Double>& hystCurrent, Vector<Double>& resCurrent, Vector<Double>& yTarget,
          Matrix<Double>& mu_inv, Matrix<Double>& jacCurrent, Vector<Double>& jacTresCurrent,
          int operatorIdx, int stayBelowSat,
          Double& alpha, Double& alphaMin, Double& alphaMax,
          bool stopLineSearchAtLocalMin, Double scalingForJacDiagonal, UInt& numberOfIterations){

    /*
     * Initial checks
     */
    if((xCurrent.NormL2() > XSaturated_)&&(stayBelowSat==1)){
      EXCEPTION("xInput to Linesearch already above saturation > must not be the case!");
    }
    if((xCurrent.NormL2() < XSaturated_)&&(stayBelowSat==-1)){
      EXCEPTION("xInput to Linesearch already below saturation > must not be the case!");
    }

    /*
     * Compute update using trial and error linesearch
     * > use alpha that minimizes jacT*res
     */
    UInt numAlphas = INV_params_.maxNumLSIts;
    Double optAlpha = 1.0;
    Double deltaAlpha = (alphaMax-alphaMin)/(numAlphas-1);
    Double curAlpha;

    bool allAlphasOutOfBound = true;
    int successFlag = 0;
    bool debugOut = false;
    bool overwriteMemory = false;

    Double minErrorNorm = 1e18;
    Vector<Double> xUpdateStart;

    if(INV_params_.inversionMethod == 1){
      xUpdateStart = computeUpdate_Newton(xCurrent, yTarget, hystCurrent, mu_inv, operatorIdx, scalingForJacDiagonal);
//      std::cout << "xCurrent = " << xCurrent.ToString() << std::endl;
//      std::cout << "xUpdateStart = " << xUpdateStart.ToString() << std::endl;
    } else if(INV_params_.inversionMethod == 2){
      xUpdateStart = computeUpdate_Krylov(xCurrent, yTarget, hystCurrent, mu_inv, operatorIdx);
    } else {
      EXCEPTION("UpdateImplementation must be 1 or 2 for std. linesearch");
    }
    Vector<Double> hystNew;
    Vector<Double> resNew;
    Vector<Double> xNew = Vector<Double>(dim_);
    Vector<Double> jacTresNew = Vector<Double>(dim_);
    numberOfIterations = 0;
    for(UInt k = 0; k < numAlphas; k++){
      // for statistics; if stopLineSearchAtLocalMin = false, numberOfIterations will always reach numAlphas
      // > this is quite inefficient! a suitable stopping criterion like the Armijo-rule would help here
      // or the quick and dirty check for locally optimal stepping (stopLineSearchAtLocalMin = true)
      numberOfIterations++;
      curAlpha = alphaMax - deltaAlpha*k;

      xNew = xCurrent;
      xNew.Add(curAlpha,xUpdateStart);

      /*
       * Check if update fits bounds
       */
      if((xNew.NormL2() < XSaturated_)&&(stayBelowSat==-1)){
//        std::cout << "Continue as xNew < XSaturated but stayBelowSat==-1!" << std::endl;
        continue;
      }
      if((xNew.NormL2() > XSaturated_)&&(stayBelowSat==1)){
//        std::cout << "Continue as xNew > XSaturated but stayBelowSat==+1!" << std::endl;
        continue;
      }
      allAlphasOutOfBound = false;

      hystNew = computeValue_vec(xNew, operatorIdx, overwriteMemory, debugOut, successFlag);
      resNew = computeResidual(xNew,yTarget,hystNew,mu_inv);

      if(INV_params_.jacImplementation == -1){
        /*
         * compute jac_dx without setting up jacobian; needs current x and current hyst values NOT the new ones!
         */
        jacTresNew = computeJacobianTimesVector(xCurrent, resNew, yTarget, hystCurrent, mu_inv, operatorIdx);
      } else {
        jacCurrent.Mult(resNew,jacTresNew);
      }

      if(jacTresNew.NormL2() < minErrorNorm){
        minErrorNorm = jacTresNew.NormL2();
        optAlpha = curAlpha;
      } else {
        /*
         * if we assume that residual behaves somehow quadratically
         * an increase in error means, that optimal alpha has already been passed
         * > break
         */
        if((stopLineSearchAtLocalMin)&&(jacTresNew.NormL2() > 1.01*minErrorNorm)){
          /*
           * much faster but does not work perfectly all the time
           */
          break;
        }
      }
    }
    
    // return for statistics
    // the last found alpha=alphaCur may not be the best one; take optAlpha
    alpha = optAlpha;
    
    xUpdate.Init();
    if(allAlphasOutOfBound){
      // update failed > discard
      return true;
    } else {
      xUpdate.Add(optAlpha,xUpdateStart);
      return false;
    }

  }

  Vector<Double> Hysteresis::projectToSolutionSpace(Vector<Double> x, int stayBelowSat){
    // stayBelowSat:
    //  1  => X = [0, XSaturated]*arbitrary direction
    // -1  => X = [XSaturated,infty]*arbitrary direction
    //  0  => X = R^2 / R^3
    Vector<Double> xProjected = Vector<Double>(dim_);
    xProjected = x;

    Double factor = 0.99999999995;
    if(stayBelowSat == 1){
      // cut length of vector down
      if(x.NormL2() >= XSaturated_){
        xProjected.ScalarMult(factor*XSaturated_/x.NormL2());
      }
    } else if(stayBelowSat == -1){
      if(x.NormL2() <= XSaturated_){
        xProjected.ScalarMult(XSaturated_/x.NormL2()/factor);
      }
    }
    return xProjected;
  }

  bool Hysteresis::computeUpdateLM_Projected(Vector<Double>& xk, Vector<Double>& resk, Vector<Double>& y,
    Matrix<Double>& mu_inv, Matrix<Double>& jac, Vector<Double>& jacTres,
    Vector<Double>& dx, UInt operatorIdx, int stayBelowSat){

    /*
     * Based on
     * "Levenberg-Marquardt methods for constrained nonlinear equations with strong local convergence properties"
     * -Kanzow,Yamashita,Fukushima
     * Algorithm 3.12 + modification as described in section 4 (use of Armijo-type linesearch)
     */
    // todo: make the following parameter accesible from mat file
    // mu, rho > 0
    Double mu = INV_params_.projLM_mu;
    Double rho = INV_params_.projLM_rho;

    // beta,sigma,gamma,tau,c in (0,1)
    Double beta = INV_params_.projLM_beta;
    Double sigma = INV_params_.projLM_sigma;
    Double gamma = INV_params_.projLM_gamma;
    Double tau = INV_params_.projLM_tau;
    Double c = INV_params_.projLM_c;

    // p > 1
    Double p = INV_params_.projLM_p;

    // parameter common to other linesearch methods
    Double alphaStart = INV_params_.alphaLSMax; // > 0
    Double alphaMin = INV_params_.alphaLSMin;
    UInt maxkk = INV_params_.maxNumLSIts;
    UInt maxjj = INV_params_.maxNumLSIts;

    /*
     * Note: no loop here! Actual outer loop can be found in computeInput_vec
     * > remark 17.12.: shouldn't this mean, that the outer loop should be allowed more iterations?
     *
     * this function evaluates steps S.1, S.2, S.3 and S.4 for a given k
     */
     /*
      * S.1 > check residual
      *
      * ||F(x_k)|| == 0 ?
      *
      * F(x_k) = current residual
      */
    if(resk.NormL2() == 0){
      dx.Init();
      // this function returns true if update shall be discarded
      // and false if it is accepted
      return false;
    }

    /*
     * S.2 > Compute dx by solving regularized system
     *
     * (H_k^T H_k) + mu_k I) dx = -H_k^T F(x_k)
     *
     * H_k = current Jacobian
     * F(x_k) = current residual
     * mu_k = mu*||F(x_k)||
     */
    Matrix<Double> jacTjac = Matrix<Double>(dim_,dim_);
    Vector<Double> jacTres_neg = Vector<Double>(dim_);

    jac.MultT(jac,jacTjac);
    jacTres_neg.Init();
    jacTres_neg.Add(-1.0,jacTres);

    dx = computeUpdate_LM_withFixedAlpha(jacTres_neg, jacTjac, resk, mu);

    /*
     * S.3 > check if update is ok
     *
     * ||F( P(x_k + dx) )|| <= gamma ||F(x_k)||
     *
     * P(x) = projection of x to allowd solution space X
     */
    /*
     * substep 1: compute projection
     */
    Vector<Double> xNew = Vector<Double>(dim_);
    xNew.Add(1.0,xk,1.0,dx);

    Vector<Double> xProjected = projectToSolutionSpace(xNew, stayBelowSat);

    /*
     * substep 2: compute residual of projected value
     */
    bool overwriteMemory = false;
    bool debugOut = false;
    int successFlag = 0;
    Vector<Double> hystNew = computeValue_vec(xProjected, operatorIdx, overwriteMemory, debugOut, successFlag);
    Vector<Double> resNew = computeResidual(xProjected,y,hystNew,mu_inv);
    /*
     * substep 3: ||F( P(x_k + dx) )|| <= gamma ||F(x_k)|| ?
     */
    if(resNew.NormL2() <= gamma*resk.NormL2()){
//      std::cout << "---------" << std::endl;
//      std::cout << "LM step: " << std::endl;
//      std::cout << "resNew.NormL2(): " << resNew.NormL2() << std::endl;
//      std::cout << "gamma*resk.NormL2(): " << gamma*resk.NormL2() << std::endl;
      // update appropriate > leave function with discardUpdate = false
      // but: as we projected xNew (and therewith our opdate), we have to retrieve it from
      // dx = xProjected-xk
      dx.Init();
      dx.Add(1.0,xProjected,-1.0,xk);

      // this function returns true if update shall be discarded
      // and false if it is accepted
      return false;
    }

    /*
     * Additional step as described in section 4
     *
     * > check if
     *        s_k = P(x_k + dx) - x_k
     *   is a valied descent direction, i.e.
     *        grad f(x_k)^T s_k <= -rho||s_k||^p
     *
     *   with f(x_k) = ||F(x_k)||^2
     *    > grad f(x_k) = 2*H_k^T F(x_k)
     *    H_k = current Jacobian
     *    F(x_k) = current residual
     *
     *   if true:
     *      > perform Armijo-type linesearch
     *   if false:
     *      > perform S.4
     */
    /*
     * substep 1: compute s_k
     */
    Vector<Double> sk = Vector<Double>(dim_);
    sk.Add(1.0,xProjected,-1.0,xk);
    /*
     * substep 2: compute grad f(x_k)^T s_k
     *
     * Note: grad f(x_k) = 2*H_k^T F(x_k)
     *                   = 2*Jac^T * res = -2*jacTres_neg
     */
    Double gradfTsk;
    jacTres_neg.Inner(sk,gradfTsk);
    gradfTsk *= -2.0;

    /*
     * substep 3: grad f(x_k)^T s_k <= -rho||s_k||^p ?
     */
    Double rhs = -rho*std::pow(sk.NormL2(),p);
    if(gradfTsk <= rhs){
//      std::cout << "---------" << std::endl;
//      std::cout << "LS step: " << std::endl;
//      std::cout << "sk: " << sk.ToString() << std::endl;
//      std::cout << "sk.NormL2(): " << sk.NormL2() << std::endl;
      /*
       * Perform Armijo-type linesearch
       * > see e.g. https://en.wikipedia.org/wiki/Backtracking_line_search
       *
       * optimize for f(x) = ||F(x)||^2
       * > p^T grad f(x) = gradfTsk from previous computation
       */
      Double criterion;
      Double t = -c*gradfTsk;
      Double alphakk = alphaStart;
      Vector<Double> hystkk;
      Vector<Double> reskk;
      Vector<Double> xkk = Vector<Double>(dim_);
      Vector<Double> xProjected_kk;

      bool successArmijo = false;
      for(UInt kk = 0; kk < maxkk; kk++){
        if(alphakk < alphaMin){
//          std::cout << "alpha < alphaMin" << std::endl;
          break;
        }
        xkk.Init();
        xkk.Add(1.0,xk,alphakk,sk);

        xProjected_kk = projectToSolutionSpace(xkk, stayBelowSat);

        hystkk = computeValue_vec(xProjected_kk, operatorIdx, overwriteMemory, debugOut, successFlag);
        reskk = computeResidual(xProjected_kk,y,hystkk,mu_inv);

        criterion = alphakk*t;
        if((resk.NormL2_squared() - reskk.NormL2_squared()) >= criterion ){
          successArmijo = true;
          break;
        }
        alphakk *= tau;
      }
      if(successArmijo){
//        std::cout << "resk.NormL2_squared(): " << resk.NormL2_squared() << std::endl;
//        std::cout << "reskk.NormL2_squared(): " << reskk.NormL2_squared() << std::endl;
//        std::cout << "alphakk: " << alphakk << std::endl;
//        std::cout << "alphakk*t: " << alphakk*t << std::endl;
//
        // update appropriate > leave function with discardUpdate = false
        // but: as we projected xNew (and therewith our opdate), we have to retrieve it from
        // dx = xProjected-xk
        dx.Init();
        dx.Add(1.0,xProjected_kk,-1.0,xk);
        // this function returns true if update shall be discarded
        // and false if it is accepted
//        std::cout << "dx: " << dx.ToString() << std::endl;
        return false;
      } else {
//        std::cout << "resk.NormL2_squared(): " << resk.NormL2_squared() << std::endl;
//        std::cout << "reskk.NormL2_squared(): " << reskk.NormL2_squared() << std::endl;
//        std::cout << "alphakk: " << alphakk << std::endl;
//        std::cout << "alphakk*t: " << alphakk*t << std::endl;
//
        dx.Init();
        dx.Add(1.0,xProjected_kk,-1.0,xk);
        // this function returns true if update shall be discarded
        // and false if it is accepted
//        std::cout << "dx: " << dx.ToString() << std::endl;

//        WARN("Armijo linsearch was not successful!");
      }

    } else {
      /*
       * Proceed with step S.4 from algorithm 3.12
       * > in fact very similar to Armijo-type linsearch above
       * > Difference:
       *    Armijo: f(x + alpha sk) <= f(x) + alpha*c*gradf(x)*sk
       *          with sk = P(x + dx) - x
       *    > f(x + alpha (P(x + dx) - x)) <= f(x) + alpha*c*gradf(x)^T*(P(x + dx) - x)
       *
       *    Here: f(x_t) <= f(x) + sigma*gradf(x)^T (x_t - x)
       *          with x_t = P(x - t*gradf(x))
       *    > f( P(x - t*gradf(x)) ) <= f(x) + sigma*gradf(x)^T(P(x - t*gradf(x))-x)
       */
      Double criterion;
      Double lStart = 0;
      Double l = lStart;
      Double t;
      Vector<Double> gradf = Vector<Double>(dim_);
      Vector<Double> xjj = Vector<Double>(dim_);
      Vector<Double> xProjected_jj = Vector<Double>(dim_);
      Vector<Double> hystjj;
      Vector<Double> resjj;
      Vector<Double> dx_jj = Vector<Double>(dim_);
      gradf.Init();
      gradf.Add(-1.0,jacTres_neg);

      bool successProjectedIt = false;
      for(UInt jj = 0; jj < maxjj; jj++){

        t = std::pow(beta,l);

        xjj.Init();
        xjj.Add(1.0,xk,-t,gradf);

        // xProjected_jj = x_t from above description
        xProjected_jj = projectToSolutionSpace(xjj, stayBelowSat);

        hystjj = computeValue_vec(xProjected_jj, operatorIdx, overwriteMemory, debugOut, successFlag);
        resjj = computeResidual(xProjected_jj,y,hystjj,mu_inv);

        dx_jj.Init();
        dx_jj.Add(1.0,xProjected_jj,-1.0,xk);

        gradf.Inner(dx_jj,criterion);
        criterion*=sigma;

        if((resk.NormL2_squared() - resjj.NormL2_squared()) >= criterion ){
          successProjectedIt = true;
          break;
        }
        l++;
      }
      if(successProjectedIt){
        // update appropriate > leave function with discardUpdate = false
        // but: as we projected xNew (and therewith our opdate), we have to retrieve it from
        // dx = xProjected-xk
        dx.Init();
        dx.Add(1.0,xProjected_jj,-1.0,xk);
        // this function returns true if update shall be discarded
        // and false if it is accepted
        return false;
      } else {
        dx.Init();
        dx.Add(1.0,xProjected_jj,-1.0,xk);
//        EXCEPTION("Projected linsearch was not successful!");
      }

    }
    // if we end up here (something went probably wrong)
    return true;
  }

  bool Hysteresis::computeUpdateLM_TrustRegion(Vector<Double>& xStart, Vector<Double>& xCurrent, Vector<Double>& xUpdate,
          Vector<Double>& hystCurrent, Vector<Double>& resCurrent, Vector<Double>& yTarget,
          Matrix<Double>& mu_inv, Matrix<Double>& jacCurrent, Vector<Double>& jacTresCurrent,
          int operatorIdx, int stayBelowSat, Double& alpha, Double alphaMin, Double alphaMax, UInt& numberOfIterations){

    /*
     * Initial checks
     */
    if((xCurrent.NormL2() > XSaturated_)&&(stayBelowSat==1)){
      EXCEPTION("xInput to Linesearch already above saturation > must not be the case!");
    }
    if((xCurrent.NormL2() < XSaturated_)&&(stayBelowSat==-1)){
      EXCEPTION("xInput to Linesearch already below saturation > must not be the case!");
    }

    UInt maxIter = INV_params_.maxNumRegIts;
    UInt itCnt = 0;
    Double alphaMinLocal = alphaMin;

    /*
     * new: success is now an integer
     * <0 no success
     * >0 success but alpha was too large
     * =0 success alpha good
     */
    Integer success = -1;
    bool discard = false;
    Vector<Double> xNew = Vector<Double>(dim_);
    xNew.Init();
    xNew.Add(xCurrent);

    /*
     * updateImplementation:
     * 0 = LevenbergMarquardt > only one which is working proper at the moment
     * 1 = Newton, direct solve
     * 2 = Newton, Krylov space solution
     */
    Matrix<Double> jacTjac = Matrix<Double>(dim_,dim_);
    Vector<Double> jacTres_neg = Vector<Double>(dim_);

    Double alphaStepUp, alphaStepDown;

    /*
     * alphaMin and Max are taken from input
     */
    alpha = alphaMinLocal;
    jacCurrent.MultT(jacCurrent,jacTjac);
    jacTres_neg.Init();
    jacTres_neg.Add(-1.0,jacTresCurrent);
    alphaStepUp = 2.0;
    alphaStepDown = 0.5;

    /*
     * Accumulated alpha
     */
    Double alphaAcc = alpha*alpha;

    /*
     * Newton and Newton-Krylov compute update only once, then use
     * alpha for stepping
     */
    int successFlag = 0;
    bool debugOut = false;
    bool overwriteMemory = false;

    Vector<Double> xUpdateStart = Vector<Double>(dim_);
    Vector<Double> hystNew = Vector<Double>(dim_);
    Vector<Double> resNew = Vector<Double>(dim_);
    Vector<Double> jac_dx = Vector<Double>(dim_);

    while(true){
      itCnt++;
			// for statistics
			numberOfIterations=itCnt;

      /*
       * LM computes new update in every iteration as it depends on alpha
       */
      xUpdate = computeUpdate_LM(jacTres_neg,
              jacTjac, alpha, alphaAcc, alphaMinLocal, 1e16);

      xNew.Init();
      xNew.Add(1.0,xCurrent,1.0,xUpdate);

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
      hystNew = computeValue_vec(xNew, operatorIdx, overwriteMemory, debugOut, successFlag);
      resNew = computeResidual(xNew,yTarget,hystNew,mu_inv);

      Vector<Double> jac_dx = Vector<Double>(dim_);

      jacCurrent.Mult(xUpdate,jac_dx);

      success = checkIncrementTrustRegion(xNew, resCurrent, resNew, jac_dx, alpha, alphaStepUp, alphaStepDown, stayBelowSat);

      //      success = checkIncrement(xNew, xUpdate, res, resNew, jac, alpha,0);
      LOG_DBG(hysteresis_linesearch) << "Check trust region - return code: " << success;

      if(success >= 0){
        LOG_DBG(hysteresis_linesearch) << "Update ok according to trust region check; now check if resulting vector may be out of bounds";
        bool wasCut = false;
        if((xNew.NormL2() > XSaturated_)&&(stayBelowSat==1)){
          // too large update
          // > set xNew back to saturation
          LOG_DBG(hysteresis_linesearch) << "2: xNew.NormL2() > Saturation but stayBelowSat = +1 > cut down xNew";

          // restrict; take amplitude of Y into account
          Double factor = 0.99999/PSaturated_*yTarget.NormL2();
          //0.99999
          xNew.ScalarMult(factor*XSaturated_/xNew.NormL2());
          wasCut = true;
        }

        if((xNew.NormL2() < XSaturated_)&&(stayBelowSat==-1)){
          // too large update
          // > set xNew back to saturation
          LOG_DBG(hysteresis_linesearch) << "3: xNew.NormL2() < Saturation but stayBelowSat = -1 (i.e. stayAboveSat = true) > cut down xNew";

          Double factor = 1.0/0.99999/PSaturated_*yTarget.NormL2();
          //1.0/0.99999
          xNew.ScalarMult(factor*XSaturated_/xNew.NormL2());
          wasCut = true;
        }

        if(wasCut){
          // retrieve updated
          xUpdate = xNew;
          xUpdate -= xCurrent;
          LOG_DBG(hysteresis_linesearch) << "Cut down xUpdate = " << xUpdate.ToString();
          // check again
          hystNew = computeValue_vec(xNew, operatorIdx, overwriteMemory, debugOut, successFlag);
          resNew = computeResidual(xCurrent,yTarget,hystNew,mu_inv);
          //          resNew = computeResidual(xNew,yVal,hystNew,mu,mu_inv,wrtX,relative);
          jacCurrent.Mult(xUpdate,jac_dx);

          success = checkIncrementTrustRegion(xNew, resCurrent, resNew, jac_dx, alpha, alphaStepUp, alphaStepDown, stayBelowSat);

          //          success = checkIncrement(xNew, xUpdate, res, resNew, jac, alpha,stayBelowSat);
          LOG_DBG(hysteresis_linesearch) << "xNew (after cut): " << xNew.ToString();
          LOG_DBG(hysteresis_linesearch) << "hystNew (after cut): " << hystNew.ToString();
          LOG_DBG(hysteresis_linesearch) << "New res vector (after cut): " << resNew.ToString();
          LOG_DBG(hysteresis_linesearch) << "Check trust region - return code (after cut): " << success;
          if(success >= 0){
            LOG_DBG(hysteresis_linesearch) << "Update still ok > take it";
            break;
          } else {
            LOG_DBG(hysteresis_linesearch) << "Update no longer ok > go to next iteration except alpha out of limits";
            // any further
            // works better with alphaAcc > alphaMax instead of alpha > alphaMax
            Double alphaCrit;
            alphaCrit = alphaAcc;

            if(alphaCrit > alphaMax){
              LOG_TRACE(hysteresis_linesearch) << "Alpha max reached; take update > stop";
              alphaCrit = alphaMax;
              // LM: large alpha > small update; adding it might not hurt to much > discard = false
              discard = false;
              alphaAcc = alphaCrit;

              break;
            }

            if(alpha < alphaMinLocal){
              LOG_TRACE(hysteresis_linesearch) << "Alpha min reached; discard update > stop";
              alpha = alphaMinLocal;
              // LM: small alpha > large update; might be risky > discard
              discard = true;

              break;
            }
            if(itCnt >= maxIter){
              LOG_TRACE(hysteresis_linesearch) << "Linesearch was not successful; Discard update.";
              discard = true;
              break;
            }
          }
        } else {
          LOG_DBG(hysteresis_linesearch) << "Update not out of bounds > take it";
          break;
        }
      } else {
        // works better with alphaAcc > alphaMax instead of alpha > alphaMax
        Double alphaCrit;
        alphaCrit = alphaAcc;

        if(alphaCrit > alphaMax){
          LOG_TRACE(hysteresis_linesearch) << "Alpha max reached; take update > stop";
          alphaCrit = alphaMax;
          // LM: large alpha > small update; BUT: here it might still be too much > discard, too
          discard = true;

          alphaAcc = alphaCrit;

          break;
        }

        if(alpha < alphaMinLocal){
          LOG_TRACE(hysteresis_linesearch) << "Alpha min reached; discard update > stop";
          alpha = alphaMinLocal;
          // LM: small alpha > large update; might be risky > discard
          discard = true;
          break;
        }

        if(itCnt >= maxIter){
          //          std::cout << "Max number of iterations used" << std::endl;
          LOG_TRACE(hysteresis_linesearch) << "Linesearch was not successful; Discard update.";
          discard = true;
          break;
        }
      }
    }

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
    LOG_TRACE(hysteresis_inversion) << " --------- START IVERSION --------- ";
    LOG_TRACE(hysteresis_inversion) << " yVal = " << yVal.ToString();
    bool debug = !true;

    if(invParamsSet_ == false){
      EXCEPTION("Hysteresis.cc: Parameter for inversion have not been set yet!")
    }
    
    /*
     * updateImplementation
     *
     * 0: Levenberg Marquardt
     * 1: Standard Newton
     * 2: Jacobian-Free-Newton-Krylov
     * 3: projected LM
     */
    UInt updateImplementation = INV_params_.inversionMethod;
    bool stopLineSearchAtLocalMin = INV_params_.stopLineSearchAtLocalMin;
    bool tolHRel = INV_params_.tolH_useAsRelativeNorm;
    bool tolBRel = INV_params_.tolB_useAsRelativeNorm;

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
      traceMsg << "hystSat: " << hystSaturated_ << std::endl;
      traceMsg << "xSat: " << XSaturated_ << std::endl;
      traceMsg << "INV_params_.maxNumIts: " << INV_params_.maxNumIts << std::endl;
      traceMsg << "INV_params_.tolH (relative? = "<<tolHRel<< "): " << INV_params_.tolH << std::endl;
      if(tolHRel){
        traceMsg << "FIXED reference value for relative tolerance => norm of xVal_old = "<< prevXval.NormL2() << std::endl;
      }
      traceMsg << "INV_params_.tolB (relative? = "<<tolBRel<< "): " << INV_params_.tolB << std::endl;
      if(tolBRel){
        traceMsg << "FIXED reference value for relative tolerance => norm of yVal = "<< yVal.NormL2() << std::endl;
      }
    }

    // for statistics
    int successFlagForward = 0;
    successFlag = -1;
    // -1 = fail
    //  0 = reuse value
    //  1 = anhyst only
    //  2 = bisection
    //  3-6 only for vector jacobianImplementation using Levenberg Marquardt, Newton or fixpoint iteration
    //  3 = remanence
    //  4 = passed dut to error tolerance
    //  5 = passed due to tolerance wrt x
    //  6 = passed due to tolerance wrt y
    //  7 = passed fixpoint due to tolerance wrt to x
    //  8 = passed fixpoint due to tolerance wrt to y

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

    Double tolB = INV_params_.tolB;
    Double referenceForRelativeError = yVal.NormL2();
    bool checkRelativeError = tolBRel;
    if(checkRelativeError && (referenceForRelativeError != 0)){
      tolB /= referenceForRelativeError;
    }
    Double tolH = INV_params_.tolH;
    referenceForRelativeError = prevXval.NormL2();
    checkRelativeError = tolHRel;
    if(checkRelativeError && (referenceForRelativeError != 0)){
      tolH /= referenceForRelativeError;
    }

    if(diff.NormL2() < tolB){
      traceMsg << "--A-- Inversion: Reuse old value" << std::endl;
      xVal = prevXval;
      LOG_TRACE(hysteresis_inversion) << "Reused value xVal: " << xVal.ToString();
      successFlag = 0;

      return xVal;
    }

    /*
     * Invert eps/mu for later usage
     *
     * Remark regarding mu:
     *  According to B = mu_0*H + J we should always use mu_0 as mu.
     *  However, the more dominant J_sat/H_sat becomes compared to mu_0
     *  the worse the inversion will get.
     *  Tests showed for example, that J_sat = 1.07, H_sat = 110e3 does
     *  have troubles with mu_0 but works very well with 100*mu_0. If
     *  J_sat = 0.107 and H_sat = 110e3, mu_0 works fine.
     * > how to fix this?
     *
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
      xVal.Init();
      if(yNorm == 0){
        // anhyst part has no remanence
        successFlag = 1;
        return xVal;
      }
      //
//      Double tol = INV_params_.tolB; //1e-12;
      
      Double Xup, Xdown, Poffset, xScal;
      Xup = yNorm/eps_mu;
      Xdown = 0;
      Poffset = 0;
      
      int successFlagBisect = -1;
      xScal = bisectForAnhyst(yNorm, Xdown, Xup, Poffset, eps_mu, tolB,successFlagBisect);
      if(successFlagBisect >= 0){
        traceMsg << "-- bisection (anhyst only) successful after " << successFlagBisect << " iterations " << std::endl;
        successFlag = 1;
      } else {
        traceMsg << "-- bisection (anhyst only) exceeded maximal number of iterations! no success!" << std::endl;
        successFlag = -1;
      }
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
     * NOTE 6.6.2020:
     *  The above is not true unfortunately! Even above saturation, Mayergoyz model (and also Sutor revised model
     *  with rotResistance (= kappa) < 1) we do not have perfect alignment (this by the way leads to non-vanishing
     *  rotational losses); Thus, the pure linesearch case would not be valid!
     * TODO: check if this is handeled correctly by flag fieldsAlignedAboveSat
     *  
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
    // note 6.6.2020: revised sutor model with rotres (=kappa) < 1 is not necessarily aligned!
    //
    // note 22.05.2018: in case of sutor model, PSaturated might be unreachable if rotRes < 1 (for revised model)
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
    
    // make use of new parameter XForAlignment_
    // see coeffunction hyst ReadAndSetParamsForHystOperator for explanation
    Vector<Double> alignmentInput = Vector<Double>(dim_);
    Vector<Double> hystValAtAlignment = Vector<Double>(dim_);
    Double anhystPartAtAlginment = 0.0;
    alignmentInput.Init();
    if(XForAlignment_ < std::numeric_limits<Double>::max()){
      alignmentInput.Add(XForAlignment_,yDir);
      hystValAtAlignment = computeValue_vec(alignmentInput, operatorIndex, false, false, successCodeTmp);
      anhystPartAtAlginment = PSaturated_*evalAnhystPart_normalized(XForAlignment_/XSaturated_);
    }
 
    Double anhystPartPosSat = PSaturated_*evalAnhystPart_normalized(1.0);
//         std::cerr << "yNorm: " << yNorm << std::endl;
//          std::cerr << "hystValAtXSat: " << hystValAtXSat.ToString(12) << std::endl;
//          std::cerr << "XSaturated_: " << XSaturated_ << std::endl;
//          std::cerr << "PSaturated_: " << PSaturated_ << std::endl;
//          std::cerr << "eps_mu: " << eps_mu << std::endl;
    //
    // check saturation in direction of yIn might solve system
    //Double diffSat = yNorm - (PSaturated_ + eps_mu*XSaturated_ + anhystPartPosSat);
    Double diffSat = yNorm - (hystValAtXSat.NormL2() + eps_mu*XSaturated_);
    // note: both diffSat and diffSatAlt should be the same if fields are perfectly aligned, but this
    // is never really the case due to roundings; diffSat is more precise; diffSatAlt is tested during
    // bisection as it does not require a reevaluation of the hysteresis operator
    // usually the bisection will work quite well but not if solution is very close to saturation;
    // (found out via test of hyst operator)
    // hotfix solution: allow for a slightly larger tolerance for diffSat; 1e-8 seems to work fine
    // alternatively, one could try to set Poffset further below not to hystSaturated_ but to 
    // hystValAtXSat.NormL2()-anhystPartPosSat > this works!
    //Double diffSatAlt = yNorm - (hystSaturated_ + eps_mu*XSaturated_ + anhystPartPosSat);

    //std::cerr << "diffSat = " << diffSat << std::endl;
    LOG_TRACE(hysteresis_inversion) << "yNorm: " << yNorm;
    if(yNorm == 0){
      stayBelowSat = 1;
    } else {
      if(abs(diffSat) < tolB){
      //if(abs(diffSat) < 1e-8){//INV_params_.tolB){
        traceMsg << "--B Special-- Inversion: (Nearly) exact Saturation found" << std::endl;
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
          if(yNorm >= (hystValAtAlignment.NormL2() + eps_mu*XForAlignment_)){
//            std::cout << "Hysteresis.cc: yNorm >= (hystValAtAlignment.NormL2() + eps_mu*XForAlignment_)" << std::endl;
//            std::cout << "Use bisection!" << std::endl;
            useBisection = true;
          }
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
        traceMsg << "--I 3-- Inversion: Unclear if input leads to X above or below Saturation" << std::endl;
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
        std::stringstream convergenceStream;
        convergenceStream << "Use bisection" << std::endl;
        if(anhystPartPosSat == 0){
          convergenceStream << "--B1-- Inversion: Anhysteretic part zero > solve by simple division" << std::endl;
          //xScal = (yNorm - PSaturated_)/eps_mu;
          // > has to be changed to value at which fields align
          // xScal = (yNorm - hystValAtXSat.NormL2())/eps_mu;
          // will not change behavior of classical model and clipped revised model as here
          // the input for alignmet = Xsaturated_
          // for revised model with kappa_rev aka. rotress < 1 however we previously never came
          // to this point; the new treatment allows for it but requires higher field
          // (for revised XForAlignment_ = Xsaturated_/kappa_rev if kappa_rev < 1)
          xScal = (yNorm - hystValAtAlignment.NormL2())/eps_mu;
          xVal.Init();
          xVal.Add(xScal,yDir);
          successFlag = 2;
        } else {
          convergenceStream << "--B2-- Inversion: Anhysteretic non-zero > solve by bisection" << std::endl;
//          Double tol = tolB; //1e-12;
          Double Xup, Xdown, Poffset;
          //Xup = (yNorm - PSaturated_)/eps_mu;
          //Poffset = PSaturated_;
          Xup = (yNorm - hystSaturated_)/eps_mu;
          //Xup = (yNorm - hystValAtXSat.NormL2())/eps_mu;
          Xdown = XForAlignment_;
          //Xdown = XSaturated_;
          //Poffset = hystSaturated_;//hystValAtXSat.NormL2();
          // subtract anhyst part to get a better approximation 
          // of the Poffset than just hystSaturated_
          Poffset = hystValAtAlignment.NormL2()-anhystPartAtAlginment;
          //Poffset = hystValAtXSat.NormL2()-anhystPartPosSat; 
          
//          std::cerr << "hystValAtXSat.NormL2(): " << hystValAtXSat.NormL2() << std::endl;
//          std::cerr << "hystSaturated_: " << hystSaturated_ << std::endl;
//          std::cerr << "hystSaturated_-hystValAtXSat.NormL2(): " << hystSaturated_-hystValAtXSat.NormL2() << std::endl;
          int successFlagBisect = -1;
          xScal = bisectForAnhyst(yNorm, Xdown, Xup, Poffset, eps_mu, tolB,successFlagBisect);
          if(successFlagBisect >= 0){
            convergenceStream << "-- bisection successful after " << successFlagBisect << " iterations " << std::endl;
            successFlag = 2;
          } else {
            convergenceStream << "-- bisection exceeded maximal number of iterations! no success!" << std::endl;
            successFlag = -1;
          }
          xVal.Init();
          xVal.Add(xScal,yDir);
          successFlag = 2;
        }
        
        if(xVal.NormL2() > 1e20){
          std::cout << convergenceStream.str() << std::endl;
          if(xVal.NormL2() > 1e20){
            std::cout << "Norm of xVal" << xVal.NormL2() << std::endl;
            EXCEPTION("Bisection went terribly wrong.");
          }
        }

        return xVal;
      }
    } // yNorm != 0

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

    LOG_DBG(hysteresis_inversion) << "Check difference between yVal and hystVal(0); remanence?";
    LOG_DBG(hysteresis_inversion) << "Diff Vector: " << diff.ToString();
    LOG_DBG(hysteresis_inversion) << "Norm of diff: " << diff.NormL2();

    if(diff.NormL2() < tolB){
      traceMsg << "--C-- Inversion: Remanence detected" << std::endl;
      xVal = xTMP;
      LOG_DBG(hysteresis_inversion) << "Set xVal to 0: " << xVal.ToString();
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

     /* Feb 2019:
      * Fixpoint x-Based
      *  yVal = mu*xVal + hystVal(xVal) = mu_FP*xVal + hystVal_FP(xVal)
      *  > xVal = mu_inv(yVal - hystVal(xVal))
      *
      * FP-iteration
      *   xVal^k+1 = 1/mu_FP(yVal - hystVal_FP(xVal^k))
      *            = 1/mu_FP(yVal - hystVal(xVal^k) - (mu - mu_FP)*xVal^k)
      *            = xVal^k + 1/mu_FP(yVal - hystVal(xVal^k) - mu xVal^k)
      *            = xVal^k + deltaX
      *
      *   mu_FP = 0.5*(min_x dy/dx + max_x dy/dx)
      *
      */
    bool warnAtNonConvergence_ = INV_params_.printWarnings;
    if(INV_params_.inversionMethod == 5){
      // TEST inversion via GLOBAL FP approach
      Double slope1,slope2,muFP;
      Matrix<Double> diagForm = Matrix<Double>(dim_,dim_);
      Vector<Double> deltaY = Vector<Double>(dim_);
      UInt maxIter = INV_params_.maxNumIts; // must be large here!

      // we need mu_FP, i.e. dB/dH
      // from tracing hyst operator we got dP/dH so we have to add mu
      slope1 = INV_params_.minSlopeHystOperator;
      slope2 = INV_params_.maxSlopeHystOperator;

      mu.GetDiagInMatrix(diagForm);
      slope1 += diagForm.GetMin();
      slope2 += diagForm.GetMax();

      // should also work for H-version with globalFPFactor_C = 1 but apparently it does not > maybe slope estimate is wrong?!
      muFP = INV_params_.safetyFactor_C*0.5*(slope1 + slope2);
//      std::cout << "slope1: " << slope1 << std::endl;
//      std::cout << "slope2: " << slope2 << std::endl;
      Vector<Double> hystVal = Vector<Double>(dim_);
      std::stringstream convergenceStream;
      convergenceStream << "Start FP inversion" << std::endl;
      for(UInt i = 0; i < maxIter; i++){
        totalNumberOfLMIterations++; // for statistics
        hystVal = computeValue_vec(xVal, operatorIndex, overwriteMemory, debugOut, successFlagForward);
        mu.Mult(xVal,deltaY);
        deltaY.Add(1.0,hystVal);
        deltaY.Add(-1.0,yVal);
        convergenceStream << "Iteration/maxNumIterations = " << i << "/" << maxIter << std::endl;
        convergenceStream << "deltaY: " << deltaY.ToString(TS_PLAIN) << std::endl;
        convergenceStream << "deltaY NormL2: " << deltaY.NormL2() << std::endl;
        // test 27.10.2020 > should we add this before returing xVal or afterwards
        xVal.Add(-1.0/muFP,deltaY);
        convergenceStream << "xVal: " << xVal.ToString(TS_PLAIN,", ", 8) << std::endl;
        convergenceStream << "hystVal: " << hystVal.ToString(TS_PLAIN,", ", 8) << std::endl;
        convergenceStream << "yVal: " << yVal.ToString(TS_PLAIN,", ", 8) << std::endl;
        if(deltaY.NormL2()/abs(mu.GetMax()) < tolH){
          successFlag = 7;
          return xVal;
        }
      }

      if(deltaY.NormL2() < tolB){
        // Failback with tolB criterion
        successFlag = 8;
        return xVal;
      } 
      else {
        if(xVal.NormL2() > 1e20){
          std::cout << "FP inversion did not succeed! here is the log:" << std::endl;
          std::cout << convergenceStream.str() << std::endl;
          std::cout << "slope1: " << slope1 << std::endl;
          std::cout << "slope2: " << slope2 << std::endl;
          std::cout << "muFP: " << muFP << std::endl;
          std::cout << "Norm of xVal" << xVal.NormL2() << std::endl;
          EXCEPTION("FP-Inversion went terribly wrong.");
        }
      }

      if(warnAtNonConvergence_){
        std::stringstream warnmsg;
        warnmsg << "Fixpoint (H-version) was not successful on element level" << std::endl;
        warnmsg << "Final H: " << xVal << std::endl;
        warnmsg << "Final deltaB.NormL2()/abs(mu.GetMax()): " << deltaY.NormL2()/abs(mu.GetMax()) << std::endl;
        warnmsg << "Final deltaB.NormL2(): " << deltaY.NormL2() << std::endl;
        warnmsg << "Try increasing the convergenceFactor to a value larger 1.0 (currently "<<INV_params_.safetyFactor_C << ")";
        warnmsg << "and/or increase the number of iterations (currently " << INV_params_.maxNumIts << ")" << std::endl;
//        WARN(warnmsg.str());
        std::cout << warnmsg.str() << std::endl;
      }
      successFlag = -1;
      return xVal;
    }

    /*
     * obtain start value
     */
    Vector<Double> xStart = obtainStartVector(prevXval, prevYval, yVal, diffRem, diffSat, stayBelowSat);
    std::string inversionName = "";
    if(INV_params_.inversionMethod == 0){
      traceMsg << "--D-- Inversion: Use Levenberg Marquardt" << std::endl;
      inversionName = "LevenbergMarquardt";
//      std::cout << "--D-- Inversion: Use Levenberg Marquardt" << std::endl;
    }
    if(INV_params_.inversionMethod == 1){
      traceMsg << "--D-- Inversion: Use Newton" << std::endl;
      inversionName = "Newton";
//      std::cout << "--D-- Inversion: Use Newton" << std::endl;
    }
    if(INV_params_.inversionMethod == 2){
      traceMsg << "--D-- Inversion: Use JacobiFreeNewtonKrylov" << std::endl;
      inversionName = "JacobiFreeNewtonKrylov";
//      std::cout << "--D-- Inversion: Use JacobiFreeNewtonKrylov" << std::endl;
    }
    if(INV_params_.inversionMethod == 3){
      traceMsg << "--D-- Inversion: Use ProjectedLMWithLinesearch" << std::endl;
      inversionName = "ProjectedLMWithLinesearch";
//      std::cout << "--D-- Inversion: Use JacobiFreeNewtonKrylov" << std::endl;
    }

    // tolerance wrt y > 1e-10 or 1e-12 seems good > takes 2-3 its
    // only problem: y-x-loops look ugly as x can be quite off!
    //Double tolError = 1e-11;
    // tolerance for reevalution

    UInt itCnt = 0;
    Double alpha;
    Double alphaMin;
    Double alphaMax;

    if(INV_params_.inversionMethod == 0){
      alpha = INV_params_.alphaRegStart;
      alphaMin = INV_params_.alphaRegMin;
      alphaMax = INV_params_.alphaRegMax;
    } else {
      alpha = 1.0;
      alphaMin = INV_params_.alphaLSMin;
      alphaMax = INV_params_.alphaLSMax;
    }

    //    Vector<Double> hystVal = computeValue_vec(xStart, operatorIndex, overwriteMemory, debugOut, successFlagForward);
    xVal = xStart;

    Double sign = 1.0;
    bool successError = false;
    bool successX = false;
    bool successY = false;
    bool discardUpdate = false;

    Vector<Double> hystVal = Vector<Double>(dim_);
    Vector<Double> xUpdate = Vector<Double>(dim_);
    Vector<Double> res = Vector<Double>(dim_);
    Matrix<Double> jac = Matrix<Double>(dim_,dim_);
    Matrix<Double> jacT = Matrix<Double>(dim_,dim_);
    Vector<Double> jacTres = Vector<Double>(dim_);

    Double errorNorm;
    Double errorNormResX;
    Double errorNormResY;

    //    INV_params_.maxNumIts *= 10;

    Vector<Double> bestSol = Vector<Double>(dim_);
    Double bestErrorNorm = 1e16;

    Vector<Double> bestSolRes = Vector<Double>(dim_);
    Double bestErrorNormRes = 1e16;

    /*
     * Start outer iteration of inversion
     */
    //    std::cout << "stayBelowSat: " << stayBelowSat << std::endl;
    UInt discardCnt = 0;

//    xVal.Init();
    Vector<Double> startSol = Vector<Double>(dim_);
    startSol = xVal;
    Double scalingForJacDiagonal = 1.0;

//    std::cout << "---------------------" << std::endl;
//    std::cout << "Start outer iteration" << std::endl;
//    std::cout << "---------------------" << std::endl;

    UInt maxIts;
    if(updateImplementation == 3){
      // in the case of projected LM linsearch, there are no inner iterations during linsearch
      // > we just check at one position and if this does not work, we get back to main iteration
      // > here we should allow for more outer iterations
      maxIts = INV_params_.maxNumLSIts*INV_params_.maxNumIts;
    } else {
      maxIts = INV_params_.maxNumIts;
    }

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

      if(INV_params_.jacImplementation != -1){
        /*
         * Compute Jacobian and its transpose
         */
        jac = computeJacobian(xVal,hystVal, mu_inv,
                operatorIndex, sign, overwriteMemory, stayBelowSat,scalingForJacDiagonal);
        jac.Transpose(jacT);
        /*
         * Compute actual error criterion (Jacobian*residual) from Jacobain
         */
        jacT.Mult(res,jacTres);
      } else {
        /*
         * Compute jacTres without setting up the actual Jacobian
         */
        jacTres = computeJacobianTimesVector(xVal, res,
                yVal, hystVal, mu_inv, operatorIndex);
      }

      /*
       * Determine three error criteria:
       * 1. jacT*residual < tolerance
       *    > hardest to satisfy
       * 2. residual wrt x, i.e. x + muInv*(hyst - y) < tolerance
       *    > should be main criterion
       * 3. residual wrt y, i.e. y - hyst - mu*x < tolerance
       *    > some sort of failback criterion
       */
      successError = checkConvergence(jacTres,errorNorm,tolH);
      errorNormResX = res.NormL2();
      successX = (errorNormResX <= tolH);
      Vector<Double> yObtained = Vector<Double>(dim_);
      successY = checkInversionOutput(xVal, yVal, yObtained, mu, tolB, errorNormResY,operatorIndex, overwriteMemory);

      if(errorNormResX < bestErrorNormRes){
        bestErrorNormRes = errorNormResX;
        bestSolRes = xVal;
      }

      if(errorNorm < bestErrorNorm){
        bestErrorNorm = errorNorm;
        bestSol = xVal;
      }

      if(successError){
        // main criterion > would be best if this one could be satisfied
        traceMsg << "Success! Error estimate |jacT*ResX| = " << errorNorm << " < " << tolH << std::endl;
        successFlag = 4;
        LOG_TRACE(hysteresis_inversion) << "Success! Error estimate |jacT*ResX| = " << errorNorm << " < " << tolH << std::endl;
        break;
      } else if(successX){
        // failback; still top if this works
        traceMsg << "Success! Residual norm wrt X = |ResX| = " << errorNormResX << " < " << tolH << std::endl;
        successFlag = 5;
        LOG_TRACE(hysteresis_inversion) << "Success! Residual norm wrt X = |ResX| = " << errorNormResX << " < " << tolH << std::endl;
        break;
      } else {
        scalingForJacDiagonal = 1.0;
        if( itCnt > 0.3*INV_params_.maxNumIts){
          // no solution has been found yet; try larger scalingForJacDiagonal
          scalingForJacDiagonal = 1.5;
        }
        if( itCnt > 0.45*INV_params_.maxNumIts){
          // no solution has been found yet; try larger scalingForJacDiagonal
          scalingForJacDiagonal = 2.0;
        }
        if( itCnt > 0.6*INV_params_.maxNumIts){
          // no solution has been found yet; try larger scalingForJacDiagonal
          scalingForJacDiagonal = 4.0;
        }
        if( itCnt > 0.75*INV_params_.maxNumIts){
          // no solution has been found yet; try larger scalingForJacDiagonal
          scalingForJacDiagonal = 8.0;
        }
        if( itCnt > 0.9*INV_params_.maxNumIts){
          // no solution has been found yet; try larger scalingForJacDiagonal
          scalingForJacDiagonal = 16.0;
        }

//        std::cout << "INV_params_.maxNumIts: " << INV_params_.maxNumIts << std::endl;
//        std::cout << "Current scaling: " << scalingForJacDiagonal << std::endl;
        if( (totalNumberOfLMIterations%10 == 0)||(itCnt >= maxIts) ){
          if(operatorIndex == 0){
              LOG_DBG(hysteresis_inversion) << "Check failback after " << itCnt << " of " << maxIts;
              LOG_DBG(hysteresis_inversion) << "Remaining norm wrt x: " << errorNormResX;
              LOG_DBG(hysteresis_inversion) << "Remaining norm wrt y: " << errorNormResY;
          }

          if(successY){
            // failback; might still have large res error in x buz might be the best we can find
            // check only each 10th iteration
            LOG_TRACE(hysteresis_inversion) << "Success! Residual norm wrt Y = |ResY| = " << errorNormResY << " < " << tolB << std::endl;
            traceMsg << "Success! Residual norm wrt Y = |ResY| = " << errorNormResY << " < " << tolB << std::endl;

//            int code = 0;
//            Vector<Double> hSol = computeValue_vec(sol, operatorIndex, false, false, code);
//            Vector<Double> ySol = Vector<Double>(dim_);
//            mu.Mult(sol,ySol);
//
//            std::cout << "xSol = " << sol.ToString() << std::endl;
//            std::cout << "h(xSol) = " << hSol.ToString() << std::endl;
//            std::cout << "mu*xSol = " << ySol.ToString() << std::endl;
//            ySol.Add(hSol);
//            std::cout << "mu*xSol + h(xSol) = " << ySol.ToString() << std::endl;
//
//            Vector<Double> hCalc = computeValue_vec(xVal, operatorIndex, false, false, code);
//            Vector<Double> yCalc = Vector<Double>(dim_);
//            mu.Mult(xVal,yCalc);
//
//            std::cout << "xCalc = " << xVal.ToString() << std::endl;
//            std::cout << "h(xCalc) = " << hCalc.ToString() << std::endl;
//            std::cout << "mu*xCalc = " << yCalc.ToString() << std::endl;
//            yCalc.Add(hCalc);
//            std::cout << "mu*xCalc + h(xCalc) = " << yCalc.ToString() << std::endl;
//
            if(operatorIndex == 0){
              LOG_DBG(hysteresis_inversion) << "Failback after " << itCnt << " of " << maxIts << " satisfied";
              LOG_DBG(hysteresis_inversion) << "Input:  " << yVal.ToString();
              LOG_DBG(hysteresis_inversion) << "Obtained output:  " << yObtained.ToString();
              LOG_DBG(hysteresis_inversion) << "Remaining norm wrt x: " << errorNormResX;
              LOG_DBG(hysteresis_inversion) << "Remaining norm wrt y: " << errorNormResY;

            }
            successFlag = 6;
            break;
          } else if (itCnt >= maxIts) {

            LOG_TRACE(hysteresis_inversion) << "NO Success! Remaining residual norm wrt Y = |ResY| = " << errorNormResY << " < " << tolB << std::endl;
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
//            std::cout << "start Sol: " << startSol.ToString() << std::endl;
//            std::cout << "last found Sol: " << xVal.ToString() << std::endl;
//            std::cout << "best found Sol wrt error Norm: " << bestSol.ToString() << std::endl;
//            std::cout << "best found Sol wrt residual Norm: " << bestSolRes.ToString() << std::endl;
//            std::cout << "actual Sol: " << sol.ToString() << std::endl;
//            std::cout << "|jacT*ResX|: " << errorNorm << std::endl;
//            std::cout << "errorNormResX: " << errorNormResX << std::endl;
//            std::cout << "errorNormResY: " << errorNormResY << std::endl;

            if(warnAtNonConvergence_){
              std::stringstream warnmsg;
              warnmsg << inversionName << " was not successful on element level" << std::endl;
              warnmsg << "Final H: " << xVal << std::endl;
              warnmsg << "Final deltaB.NormL2()/abs(mu.GetMax()): " << errorNormResX << std::endl;
              warnmsg << "Final deltaB.NormL2(): " << errorNormResY << std::endl;
              std::cout << warnmsg.str() << std::endl;
//              WARN(warnmsg.str());
            }
//            if(operatorIndex == 0){
//              std::cout << "++FAIL++ Inversion of operator " << operatorIndex << " failed after " << maxIts << " iterations" << std::endl;
//              LOG_DBG(hysteresis_inversion) << "++FAIL++ Inversion of operator " << operatorIndex << " failed after " << maxIts << " iterations";
//              Vector<Double> tmp = Vector<Double>(dim_);
//              tmp = yVal;
//              tmp.Add(-1.0,hystVal);
//              Vector<Double> tmp2 = Vector<Double>(dim_);
//              mu.Mult(xVal,tmp2);
//              tmp.Add(-1.0,tmp2);
//              std::cout << "YIn - mu*XFound - P(XFound) = " << std::endl;
//              std::cout << yVal[0] << " - (" << mu[0][0] << " " << mu[0][1] << ")*" << xVal[0] << " - " << hystVal[0] << " = " << tmp[0] << std::endl;
//              std::cout << yVal[1] << "   (" << mu[1][0] << " " << mu[1][1] << ") " << xVal[1] << "   " << hystVal[1] << "   " << tmp[1] << std::endl;
//              std::cout << "res wrt y = " << tmp.NormL2() << std::endl;
//
//              std::cout << "For comparison: |YIn| = " << yVal.NormL2() << " / |hystValAtXSat| = " << hystValAtXSat.NormL2() << " / eps_mu*XSaturated_ = " << eps_mu*XSaturated_ << std::endl;
//              if(stayBelowSat != 0){
//                std::cout << "Solution was forced to be ";
//                if(stayBelowSat == -1){
//                  std::cout << " above saturation" << std::endl;
//                } else if(stayBelowSat == 1){
//                  std::cout << " below saturation" << std::endl;
//                }
//              } else {
//                std::cout << "Solution was not forced to be above or below saturation" << std::endl;
//              }
//
              LOG_DBG(hysteresis_inversion) << "YIn - mu*XFound - P(XFound) = ";
              LOG_DBG(hysteresis_inversion) << yVal[0] << " - (" << mu[0][0] << " " << mu[0][1] << ")*" << xVal[0] << " - " << hystVal[0] << " = " << tmp[0];
              LOG_DBG(hysteresis_inversion) << yVal[1] << "   (" << mu[1][0] << " " << mu[1][1] << ") " << xVal[1] << "   " << hystVal[1] << "   " << tmp[1];
              LOG_DBG(hysteresis_inversion) << "res wrt y = " << tmp.NormL2();

              LOG_DBG2(hysteresis_inversion) << "For comparison: |YIn| = " << yVal.NormL2() << " / |hystValAtXSat| = " << hystValAtXSat.NormL2() << " / eps_mu*XSaturated_ = " << eps_mu*XSaturated_;
              if(stayBelowSat != 0){
                LOG_DBG2(hysteresis_inversion) << "Solution was forced to be ";
                if(stayBelowSat == -1){
                  LOG_DBG2(hysteresis_inversion) << " above saturation";
                } else if(stayBelowSat == 1){
                  LOG_DBG2(hysteresis_inversion) << " below saturation";
                }
              } else {
                LOG_DBG2(hysteresis_inversion) << "Solution was not forced to be above or below saturation";
              }

//              if( INV_params_.jacImplementation == 2 ){
//                Vector<Double> tmp = Vector<Double>(dim_);
//                tmp = yVal;
//                tmp.Add(-1.0,hystVal);
//                Vector<Double> tmp2 = Vector<Double>(dim_);
//                mu.Mult(xVal,tmp2);
//                tmp.Add(-scalingForJacDiagonal,tmp2);
//                std::cout << "YIn - lastscaling*mu*XFound - P(XFound) = " << std::endl;
//                std::cout << yVal[0] << " - "<< scalingForJacDiagonal << "\t*(" << mu[0][0] << " " << mu[0][1] << ")*" << xVal[0] << " - " << hystVal[0] << " = " << tmp[0] << std::endl;
//                std::cout << yVal[1] << "     \t (" << mu[1][0] << " " << mu[1][1] << ") " << xVal[1] << "   " << hystVal[1] << "   " << tmp[1] << std::endl;
//                std::cout << "res wrt y = " << tmp.NormL2() << std::endl;
//              }
//            }

            break;
          }
        }
      }

      LOG_DBG(hysteresis_inversion) << "Solution not appropriate yet > Perform linesearch";
      UInt numberOfIterations = 0;

      /*
       * Compute new update
       */
      //      discardUpdate = computeUpdate(xVal, yVal, res, xUpdate, jac, jacT, mu, mu_inv,
      //              operatorIndex, overwriteMemory, alpha, alphaMin, alphaMax,
      //              numberOfIterations,xStart,stayBelowSat,sol);
      //
      if(updateImplementation == 0){
        discardUpdate = computeUpdateLM_TrustRegion(xStart, xVal, xUpdate, hystVal, res, yVal, mu_inv,
                jac, jacTres, operatorIndex, stayBelowSat,
                alpha, alphaMin, alphaMax, numberOfIterations);
      } else if(updateImplementation == 3){
        discardUpdate = computeUpdateLM_Projected(xVal, res, yVal, mu_inv, jac, jacTres, xUpdate, operatorIndex, stayBelowSat);
      } else {
        discardUpdate = computeUpdateNewton_Linesearch(xStart, xVal, xUpdate, hystVal, res, yVal, mu_inv,
                jac, jacTres, operatorIndex, stayBelowSat,
                alpha, alphaMin, alphaMax, stopLineSearchAtLocalMin, scalingForJacDiagonal, numberOfIterations);
      }

      LOG_DBG(hysteresis_inversion) << "Computed update: " << xUpdate.ToString();
      if(!discardUpdate){
        xVal = xVal+xUpdate;
      } else {
        discardCnt++;
        LOG_DBG(hysteresis_inversion) << "Discard update";
        bool useFailback = false;
        if(useFailback){
          //std::cout << "Discard update; reset to best solution so far" << std::endl;
          //        std::cout << "Update not usable: use failback inversion method = Newton" << std::endl;
          updateImplementation = 1;
          stopLineSearchAtLocalMin = false;
          alpha = 1;
          alphaMin = 0.00025;
          alphaMax = 1;
          if(discardCnt!=1){
            //          std::cout << "Still no success: change alpha min and alpha max" << std::endl;
            alphaMin = alphaMin/2.0;
            alphaMax = alphaMax*2.0;
          }
        } else {
          alphaMin = alphaMin/2.0;
          alphaMax = alphaMax*2.0;
        }
      }
      sign = sign*(-1.0);

      /*
       * Collect some statistics
       */
      if(alpha < minAlphaStatistics){
        minAlphaStatistics = alpha;
      }
      if(alpha > maxAlphaStatistics){
        maxAlphaStatistics = alpha;
      }
      avgAlphaStatistics += alpha;

      totalNumberOfLinesearchIterations += numberOfIterations;

      if(numberOfIterations > maximalNumberOfLinesearchIterations){
        maximalNumberOfLinesearchIterations = numberOfIterations;
      }
    }

    avgAlphaStatistics /= totalNumberOfLMIterations;

    /*
     * Final checks
     */
    if((xVal.NormL2() > XSaturated_)&&(stayBelowSat==1)){
      EXCEPTION("LM lead xVal into saturation although input is below > must not be the case!");
    }

    if((xVal.NormL2() < XSaturated_)&&(stayBelowSat==-1)){
      EXCEPTION("LM lead xVal below saturation although input is above> must not be the case!");
    }

    bool abortOnFail = false;
    if(checkInversionResult_){
      Double resYNorm = 0.0;
      Vector<Double> yObtained = Vector<Double>(dim_);
      bool finalCheckPassed = checkInversionOutput(xVal, yVal, yObtained, mu, tolB, resYNorm, operatorIndex, false,true);
      if(!finalCheckPassed){
        std::stringstream exceptionmsg;
        exceptionmsg << inversionName << "-Inversion failed final test for requested input: " << yVal.ToString() <<  std::endl;
        exceptionmsg << "retrievedOutput: " << yObtained.ToString() << " -  remaining residual norm wrt y: " << resYNorm;
        if(abortOnFail){
          EXCEPTION(exceptionmsg.str());
        } else {
          //WARN(exceptionmsg.str());
          std::cout << exceptionmsg.str() << std::endl;
        }
      }
    }

    LOG_TRACE(hysteresis_inversion) << " --------- END IVERSION --------- ";
//    std::cout << "Finished local inversion after " << itCnt << " terations" << std::endl;
    return xVal;
  }

	// returns true, if yTarget can be reached by mu*xComputed + hyst(xComputed)
	bool Hysteresis::checkInversionOutput(Vector<Double>& xComputed, Vector<Double>& yTarget, Vector<Double>& yObtained,
          Matrix<Double>& mu, Double tol, Double& resYNorm, Integer operatorIdx, bool overwriteMemory, bool output){

		Vector<Double> yCheck = Vector<Double>(dim_);
    bool debugOut = false;
    int successCode = 0;
    Vector<Double> hCheck = computeValue_vec(xComputed, operatorIdx, overwriteMemory, debugOut, successCode);

		mu.Mult(xComputed,yCheck);
		yCheck.Add(hCheck);
    yObtained = yCheck;

		Vector<Double> diff;
		diff = yCheck;
		diff -= yTarget;
    resYNorm = diff.NormL2();
    if(output){
      LOG_TRACE(hysteresis_inversion) << "CheckInversionOutput - norm of residual wrt y-vector: " << diff.NormL2();
    }
		if(diff.NormL2() < tol){
      if(output){
        LOG_TRACE(hysteresis_inversion) << "Solution OK; |diff| = " << diff.NormL2() << " < " << tol;
      }
			return true;
		} else {
      if(output){
        LOG_TRACE(hysteresis_inversion) << "Solution NOT OK; |diff| = " << diff.NormL2() << " > " << tol;
        LOG_TRACE(hysteresis_inversion) << "yRetrieved = " << yCheck.ToString();
        LOG_TRACE(hysteresis_inversion) << "yIn = " << yTarget.ToString();
      }
			return false;
		}
	}

}
