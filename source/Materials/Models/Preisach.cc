#include "Preisach.hh"

#include <iostream>
#include <fstream>
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Utils/Timer.hh"

namespace CoupledField
{ 
  class Preisach;
  
  DEFINE_LOG(scalpreisach, "scalpreisach")
  DEFINE_LOG(scalpreisachInversion, "scalpreisachInversion")
  DEFINE_LOG(scalpreisachVecExtension, "scalpreisachVecExtension")

  Preisach::Preisach(Integer numElem, ParameterPreisachOperators operatorParams,
          ParameterPreisachWeights weightParams, bool isVirgin, bool ignoreAnhystPart)
  : Hysteresis(numElem,operatorParams,weightParams){

//  Preisach::Preisach(Integer numElem, Double xSat, Double ySat,
//          Matrix<Double>& preisachWeight, bool isVirgin, Double anhyst_A, Double anhyst_B, Double anhyst_C,bool anhystOnly)
//  : Hysteresis(numElem,xSat,ySat,anhyst_A,anhyst_B,anhyst_C,anhystOnly)
//  {
    LOG_TRACE(scalpreisach) << "ScalarPreisach: Using Everett function";
    //tol_ = 1e-5;
    tol_ = 1e-15;

    // set in base class
//    if (xSat > 0 ) {
//      XSaturated_  = xSat;
//    }
//    else {
//      XSaturated_  = 1.0;
//    }
//
//    PSaturated_  = ySat;
    
    fixDirection_ = operatorParams.fixDirection_;
    dim_ =  fixDirection_.GetSize();

    isVirgin_    = isVirgin;
    
    preisachWeights_ = weightParams.weightTensor_;

    preisachSum_.Resize(numElem);
    preisachSum_.Init();
    StringLength_.Resize(numElem);

    // needed in combination with Mayergoyz model
    // > only pure hysteretic part wanted in that case
    if(ignoreAnhystPart){
      anhyst_A_ = 0;
      anhyst_B_ = 0;
      anhyst_C_ = 0;
    }
    
    strings_     = new Vector<Double>[numElem];
    helpStrings_ = new Vector<Double>[numElem];
    minmaxtype_ = new Vector<Integer>[numElem];
    evaluatedEverettPixel_ = new Vector<Double>[numElem];
    maxStringLength_ = 200;

    for (Integer el=0; el<numElem; el++) {
      strings_[el].Resize(maxStringLength_);
			minmaxtype_[el].Resize(maxStringLength_);
			evaluatedEverettPixel_[el].Resize(maxStringLength_);
			
      helpStrings_[el].Resize(maxStringLength_+1);
      
      StringLength_[el] = 1;
      for ( UInt i=0; i<maxStringLength_; i++) {
        strings_[el][i] = 0.0;
				minmaxtype_[el][i] = 0; // neither min nor max
				evaluatedEverettPixel_[el][i] = 0.0;
      }
    }
    
    // for inversion
    previousXval_ = Vector<Double>(numElem);
    previousXval_.Init();
    
    previousPval_ = Vector<Double>(numElem);
    previousPval_.Init();
  }
  
  Preisach::~Preisach()
  {
    delete [] strings_;
    delete [] helpStrings_;
		delete [] evaluatedEverettPixel_;
		delete [] minmaxtype_;
  }
  
  Double Preisach::getValue( Integer idx ) 
  {
    return ( preisachSum_[idx]*PSaturated_ );
  }
  

  Double Preisach::bisect(Double dY,Double xMin,Double xMax, Double xFixed, Double eps_mu, Double tol){
		LOG_DBG(scalpreisachInversion) << "perform bisection for dY: " << dY;
//    std::cout << "perform bisection" << std::endl;
//    std::cout << "dY: " << dY << std::endl;
//    std::cout << "xMin: " << xMin << std::endl;
//    std::cout << "xMax: " << xMax << std::endl;
//    std::cout << "xFixed: " << xFixed << std::endl;
    
    /*
     *  Try to find xOpt, such that
     * 
     *  abs(dY - dP - dX_to_dY*dX) = min
     *  with dP = everettPixel(xFixed,xOpt)
     *  and  dX = xOpt-xFixed
     * 
     * if xFixed < -1, find xOpt, such that
     * 
     *  abs(dY) = abs(everettPixel(-xOpt,xOpt)) = min
     * 
     * > use bisection to find xOpt between xMin and xMax
     */ 
    UInt maxIter = 50;
    Double xMiddle;
    Double diffMin, diffMiddle;
    Double everettMin, everettMiddle;
    Double dX_to_dY = eps_mu*XSaturated_/PSaturated_;

    // new: as preisach operator alone leads to hystSaturated instead of PSaturated
    //        but dY is normalized to PSaturated, we have to add the scaling dHyst_to_dY to all
    //        everett related parts (as those return +/-1 which corresponds to +/-hystSaturated
    Double dHyst_to_dY = hystSaturated_/PSaturated_;
    Double dX = 0.0;
    Double dAnhyst = 0.0;
    
    if(xFixed >= -1){
      LOG_DBG(scalpreisachInversion) << "perform bisection with fixed edge";
      //std::cout << "bisection with fixed edge" << std::endl;
      for(UInt i = 0; i<maxIter; i++){
        xMiddle = 0.5*(xMax+xMin);
        dX = xMiddle - xFixed;
        dAnhyst = evalAnhystPart_normalized(xMiddle) - evalAnhystPart_normalized(xFixed);
        // important: here we need a factor of 2.0 as we have an update
        // to the previous state (i.e. we have to revert from +1 to -1 or
        // vice versa)
        everettMiddle = 2.0*everettPixel(xFixed,xMiddle);
        // dY = dP + dX_to_dY * dX
        diffMiddle = dY - dHyst_to_dY*everettMiddle - dX*dX_to_dY - dAnhyst;
        
        if(abs(diffMiddle) < tol){
					LOG_DBG(scalpreisachInversion) << "Remaining diff (normalized): " << abs(diffMiddle);
          LOG_DBG(scalpreisachInversion) << "Remaining diff: " << PSaturated_*abs(diffMiddle);
          return xMiddle;
        }
        
        dX = xMin - xFixed;
        dAnhyst = evalAnhystPart_normalized(xMin) - evalAnhystPart_normalized(xFixed);

        // syntax for everettPixel: 1. parameter: older/previous entry, 2. parameter new/next entry
        everettMin = 2.0*everettPixel(xFixed,xMin);

        diffMin = dY - dHyst_to_dY*everettMin - dX*dX_to_dY - dAnhyst;
        
//				std::cout << "xMiddle: " << xMiddle << std::endl;
//        std::cout << "everettMiddle: " << everettMiddle << std::endl;
//        std::cout << "everettMin: " << everettMin << std::endl;
//        std::cout << "diffMiddle: " << diffMiddle << std::endl;
//        std::cout << "diffMin: " << diffMin << std::endl;
        if( ((diffMiddle < 0)&&(diffMin < 0)) || ((diffMiddle > 0)&&(diffMin > 0)) ){
          // search between middle and max
          xMin = xMiddle;
        } else {
          // search between min and middle
          xMax = xMiddle;
        }
      }
      LOG_DBG(scalpreisachInversion) << "Maxnumber of iteations reached ( case: xFixed >= -1 )";
			LOG_DBG(scalpreisachInversion) << "Remaining diff (normalized): " << abs(diffMiddle);
      LOG_DBG(scalpreisachInversion) << "Remaining diff: " << PSaturated_*abs(diffMiddle);
			return xMiddle;
    } else {
      LOG_DBG(scalpreisachInversion) << "perform bisection without fixed edge";
      // for the open edge case, we have to be careful as the everettPixel
      // inside an empty Preisach plane will always be positive for x in 0,1;
      // to get a negative dY, we have to search in -1,0 instead
      if(dY < 0){
        xMax = -xMax;
        xMin = -xMin;
      }
      for(UInt i = 0; i<maxIter; i++){
        //std::cout << "bisection with open edge" << std::endl;
        xMiddle = 0.5*(xMax+xMin);
        dX = xMiddle; // here we have dX towards 0
        dAnhyst = evalAnhystPart_normalized(xMiddle);
        // here we do not need a factor of 2.0 as we do not have to
        // revert a previous state
        everettMiddle = everettPixel(-xMiddle,xMiddle);

        diffMiddle = dY - dHyst_to_dY*everettMiddle - dX*dX_to_dY - dAnhyst;
        
        if(abs(diffMiddle) < tol){
          LOG_DBG(scalpreisachInversion) << "Remaining diff (normalized): " << abs(diffMiddle);
          LOG_DBG(scalpreisachInversion) << "Remaining diff: " << PSaturated_*abs(diffMiddle);
          return xMiddle;
        }
        
        // note that everettPixel always returns a positive value (the area)
        // and a sign that can be used to find out if this area is to be subtracted
        // or to be added
        dX = xMin; // here we have dX towards 0
        dAnhyst = evalAnhystPart_normalized(xMin);
        everettMin = everettPixel(-xMin,xMin);
        //everettMax = everettPixel(-xMax,xMax,sign);
				
        diffMin = dY - dHyst_to_dY*everettMin - dX*dX_to_dY - dAnhyst;
        //diffMin = dYabs - abs(everettMin + dX*dX_to_dY);
        //diffMax = dYabs - everettMax;
				
//        std::cout << "xMiddle: " << xMiddle << std::endl;
//        std::cout << "everettMiddle: " << everettMiddle << std::endl;
//        std::cout << "everettMin: " << everettMin << std::endl;
//        std::cout << "diffMiddle: " << diffMiddle << std::endl;
//        std::cout << "diffMin: " << diffMin << std::endl;
        if( ((diffMiddle < 0)&&(diffMin < 0)) || ((diffMiddle > 0)&&(diffMin > 0)) ){
          // search between middle and max
          xMin = xMiddle;
        } else {
          // search between min and middle
          xMax = xMiddle;
        }
      }
      LOG_DBG(scalpreisachInversion) << "Maxnumber of iteations reached ( case: xFixed < -1 )";
		  LOG_DBG(scalpreisachInversion) << "Remaining diff (normalized): " << abs(diffMiddle);
      LOG_DBG(scalpreisachInversion) << "Remaining diff: " << PSaturated_*abs(diffMiddle);
			return xMiddle;
    }
  }

  Double Preisach::computeInputAndUpdate(Double Yin, Double eps_mu, Integer idx, bool overwrite, int& successFlag){
    /*
     * NEW implementation; this shall replace class SimplePreisachInv that
     * apparently is not working
     * 
     * 
     * Source: "Generalisiertes Preisach-Modell für die Simulation uund Kompensation
     *         piezokeramischer Aktoren" - Dissertation, Felix Wolf, p. 127ff
     */
    
		LOG_TRACE(scalpreisachInversion) << "Compute inverse of Preisach operator for Yin = " << Yin;
    LOG_TRACE(scalpreisachInversion) << "Index = " << idx;
		LOG_DBG(scalpreisachInversion) << "Compute inverse of Preisach operator for Yin_normalized = " << Yin/PSaturated_;
    /*
     * 0. Check if Input drives system into saturation
     *    > make sure to not only compare with YSaturated but also 
     *      include the contribution of XSaturated!
     *    > if eps_mu*XSaturated_ is not considered
     *      we could return XSaturated as Yin > YSaturated but
     *      when recomputing Yin via YSaturated + eps_mu*XSaturated
     *      we can get a significant difference to the actual Yin
     */
    Double Xout;
		Double x1,x2; // search interval for later bisection
		x1 = 0.0;
		x2 = 0.0;
    Double tol = 1e-8/XSaturated_;
    // dX_to_dY is the factor which is needed to transfer a normalized dX to
    // a normalized dY
    // dY_norm = dY/YSaturated = eps_mu*dX/YSaturated = eps_mu*dX*XSaturated/XSaturated/YSaturated
    //         = eps_mu*XSaturated/YSaturated * dX_norm
    Double dX_to_dY = eps_mu*XSaturated_/PSaturated_;
    // new: as preisach operator alone leads to hystSaturated instead of PSaturated
    //        but dY is normalized to PSaturated, we have to add the scaling dHyst_to_dY to all
    //        everett related parts (as those return +/-1 which corresponds to +/-hystSaturated
    Double dHyst_to_dY = hystSaturated_/PSaturated_;
    
    UInt invcase = 0;
		UInt subcase = 0;
    //std::cout << "Inversion" << std::endl;
    //std::cout << "Yin: " << Yin << std::endl;
    //std::cout << "YSaturated: " << PSaturated_ << std::endl;
    //std::cout << "eps_mu: " << eps_mu << std::endl;
    // anhysteretic parameter are meant to be combined with saturated X > mult by 1.0 here
    Double anhystPartPosSat = PSaturated_*evalAnhystPart_normalized(1.0);
    Double anhystPartNegSat = PSaturated_*evalAnhystPart_normalized(-1.0);

    if(anhystOnly_ == true){
      // actual hyst operator is 0; we just solve for the anhyst part via bisection
      Double Xup, Xdown, Poffset;
      Xup = Yin/eps_mu;
      Xdown = 0.0;
      Poffset = 0.0;
      Xout = bisectForAnhyst(Yin, Xdown, Xup, Poffset, eps_mu, tol);
      //Xout = bisectForSaturation(Yin, eps_mu, tol, false);
      invcase = 9;
      successFlag = 1;
    } else {
      // we just invert the actual hysteresis operator, i.e we have to subtract all anhyst and reversible parts
      //  before starting the inversion
      if(Yin >= (hystSaturated_ + eps_mu*XSaturated_ + anhystPartPosSat) ){
        LOG_DBG(scalpreisachInversion) << "Pos saturation";
        //std::cout << "pos saturation" << std::endl;
        /*
         * System is actually in positive saturation
         *
         * yIn = Ysaturated_ + eps_mu*xOut + PSaturated_*(anhyst_A_*std::atan(anhyst_B_*xOut/Xsaturated_) + anhyst_C_*xOut/Xsaturated_)
         *     >= PSaturated_ + eps_mu*XSaturated_ + PSaturated_*(anhyst_A_*std::atan(anhyst_B_) + anhyst_C_);
         *
         * > inversion without anhysteretic part
         * xOut = (yIn - Ysaturated_)/eps_mu
         *
         * > with anhysteretic part, inversion more complicated!
         */
        if(anhystPartPosSat == 0){
          Xout = (Yin - hystSaturated_)/eps_mu;
          invcase = 1;
        } else {
          Double Xup, Xdown, Poffset;
          Xup = (Yin - hystSaturated_)/eps_mu;
          Xdown = XSaturated_;
          Poffset = hystSaturated_;
          Xout = bisectForAnhyst(Yin, Xdown, Xup, Poffset, eps_mu, tol);
          //Xout = bisectForSaturation(Yin, eps_mu, tol, false);
          invcase = 11;
        }
        successFlag = 2;
      } else if (Yin <= (-hystSaturated_ - eps_mu*XSaturated_ + anhystPartNegSat) ){
        LOG_DBG(scalpreisachInversion) << "Neg saturation";
        //std::cout << "neg saturation" << std::endl;
        /*
         * System is actually in negative saturation
         *
         * yIn = -Ysaturated_ + eps_mu*xOut
         *
         * xOut = (yIn + Ysaturated_)/eps_mu
         */
        if(anhystPartPosSat == 0){
          Xout = (Yin + hystSaturated_)/eps_mu;
          invcase = 2;
        } else {
          Double Xup, Xdown, Poffset;
          Xup = (Yin + hystSaturated_)/eps_mu;
          Xdown = -XSaturated_;
          Poffset = -hystSaturated_;
          Xout = bisectForAnhyst(Yin, Xdown, Xup, Poffset, eps_mu, tol);
          //Xout = bisectForSaturation(Yin, eps_mu, tol, true);
          invcase = 21;
        }
        successFlag = 2;
      } else {
        invcase = 3;
        /*
         * 1. Compute difference between requested yVal and previously computed value
         */
        // factor for relating a change in x to a change in y for
        // the non-hysteretic part
        //  > dY = dP + dX_to_dY*dX
        // for dY, dP and dX all normalized

        // TODO:; should be normalized AFTER subtracting x
        // do not cap to +/-1 here as Y is not the polarization (that is actually capped!)
        // EDIT: 26.3.2018
        //  we have to cap, but not to +/-1 (that would be for polarization) but
        //  to +/- (1 + dX_to_dY*Xsaturated)
        // > reason: in the following we compute dY as the difference between oldthat can
        //           be represented by the Preisach plane is
        //           2*PSaturated_ + 2*eps_mu*XSaturated_
        //           (considering the reversible part);
        //           If we start in saturation, however, we have a previous |Y| > YSaturated.
        //           This leads to a dY that can be larger than 2*PSaturated_ + 2*eps_mu*XSaturated_.
        //           Even if this is not the case, the later finesearch will try to compensate this
        //           dY by adapting the Preisach plane. This works quite well actually, but is wrong
        //           apparently as the finesearch computes dY from YSaturated+eps_mu*XSaturated and not from the actual
        //           Yold.
        //            > restrict Yold to (YSaturated+eps_mu*XSaturated)/YSaturated
        //

        // OLD without capping > works well unless we come from a saturated state
  //      Double Y_normalized = Yin/PSaturated_;
  //      Double P_old_normalized = previousPval_[idx];
  //      Double Y_old_normalized = P_old_normalized+dX_to_dY*previousXval_[idx];
  //      LOG_DBG(scalpreisach) << "yOld: " << Y_old_normalized*PSaturated_;
  //      Double dY = Y_normalized - Y_old_normalized;

        // NEW with capping
        // NEW 25.4.2018 with anhyst part
        Double Y_normalized = Yin/PSaturated_;
        Double P_old_normalized = previousPval_[idx];
        // IMPORTANT NOTE: previousXval_[idx] is normalized but NOT clipped > exactly what we need!
        // NOTE: previousPval_ contains already the anhysteretic part and in total is normalized to PSaturated_!
        Double Y_old_normalized = P_old_normalized + dX_to_dY*previousXval_[idx];// + evalAnhystPart_normalized(previousXval_[idx]);
        LOG_DBG(scalpreisachInversion) << "Y_old: " << Y_old_normalized*PSaturated_;
        LOG_DBG(scalpreisachInversion) << "Y_old_normalized: " << Y_old_normalized;
        if(Y_old_normalized > (dHyst_to_dY + dX_to_dY + evalAnhystPart_normalized(1.0))){
          Y_old_normalized = dHyst_to_dY + dX_to_dY + evalAnhystPart_normalized(1.0);
        } else if (Y_old_normalized < -(dHyst_to_dY + dX_to_dY)+evalAnhystPart_normalized(-1.0)){
          Y_old_normalized = -dHyst_to_dY - dX_to_dY + evalAnhystPart_normalized(-1.0);
        }
        LOG_DBG(scalpreisachInversion) << "Y_old_normalized (clipped): " << Y_old_normalized;

        Double dY = Y_normalized - Y_old_normalized;


  //     dY = dY - dX_to_dY*previousXval_[idx];
        //std::cout << "Y_normalized - P_old_normalized - dX_to_dY*previousXval_[idx]: " << dY << std::endl;
        LOG_TRACE(scalpreisachInversion) << "Starting value for dY: " << dY*PSaturated_;
        LOG_DBG(scalpreisachInversion) << "Starting value for dY (normalized): " << dY;
        Integer minmaxcur = 0.0;
        if(dY > 0){
          // larger value than before > input has to be larger than
          // previous one (this is only true if all Preisach weights >= 0!)
          // > retrieved input will be a maximum
          minmaxcur = 1.0;
        } else if(dY < 0){
          minmaxcur = -1.0;
        }
        //std::cout << "minmaxcur: " << minmaxcur << std::endl;
        /*
         * 2. Check if difference to previous value is relevant or not
         */
        if(abs(dY) < 1e-16){
          // difference is negligible
          LOG_TRACE(scalpreisachInversion) << "take previous xvalue" << std::endl;
          invcase = 31;
          successFlag = 0;
          // attention: previousXval_ is normalized by XSaturated_
          Xout = previousXval_[idx];
        } else{
          invcase = 32;
          Vector<Double> &stringEl = strings_[idx];
          Vector<int> &minmaxEl = minmaxtype_[idx];
          Vector<Double> &everettEl = evaluatedEverettPixel_[idx];

          x1 = -1.0; // lower bound for max, left bound for min
          x2 = 1.0; // upper bound for max, right bound for min

          Double xfix;
          UInt actLength = StringLength_[idx];

          /*
           * first step: list has to be cut down such that it ends with
           *						 the opposite minmaxtype (as it is easier to add
           *						 a max after a min than adding a max after a max)
           */
          if( (actLength >= 1)&&(minmaxcur == minmaxEl[actLength-1]) ){
            LOG_DBG(scalpreisachInversion) << "remove last entry of list: " << stringEl[actLength-1];
            //std::cout << "remove last entry of list: " << stringEl[actLength-1] << std::endl;
            // delete last entry (not actually here but just virtually)
            // the value of the corresponding everett pixel now increases
            // the difference dY of course
            // ATTENTION: if eps_mu*XSaturated_ comes into the range of 1e-3 or larger
            //            we might find an everett pixel that fits dY but the result
            //            will still be wrong as dY actually computes to be
            //              Y = everettPixel(X) + eps_mu*X
            // Conclusion: we have to subtract the influence of eps_mu*X from the
            //             searched area update
            Double dP = everettEl[actLength-1];
            //std::cout << "dP (=everett[actLength-1]): " << dP << std::endl;

            dY += dHyst_to_dY*dP;
            //std::cout << "dY+dP: " << dY << std::endl;

            // we search space for dY not dP, so we have to add the contribution of
            // the reversible part
            Double dX = stringEl[actLength-1];
            Double dAnhyst = evalAnhystPart_normalized(stringEl[actLength-1]);

            LOG_DBG(scalpreisachInversion) << "Invert Preisach - Add anhystPart " << evalAnhystPart_normalized(stringEl[actLength-1]);
            LOG_DBG(scalpreisachInversion) << "for X/norm = " << (stringEl[actLength-1]);

            if(actLength >= 2){
              LOG_DBG(scalpreisachInversion) << "remove a second entry: " << stringEl[actLength-2];
              dX -= stringEl[actLength-2];
              dAnhyst -= evalAnhystPart_normalized(stringEl[actLength-2]);
            }
            //std::cout << "dX = " << dX << std::endl;
            dY += dX*dX_to_dY + dAnhyst;
            //std::cout << "dY+dP+dX_to_dY*dX: " << dY << std::endl;

            if(minmaxcur > 0){
              // max> max was removed; lower bound can be increased to removed entry
              x1 = stringEl[actLength-1];
            } else {
              // min> min was removed; right bound can be decreased to removed entry
              x2 = stringEl[actLength-1];
            }
            LOG_DBG(scalpreisachInversion) << "Updated dY: " << dY;
            actLength--;
          }

          /*
           * second step: coarse search > find everettPixel that is large enought
           *							to hold dY
           */
          while(true){
            if(actLength <= 0){
              //std::cout << "list empty" << std::endl;
              // list empty

              // in bisection, search for x such that triangle
              // everett(-x,x) = dY
              // to indicate that we want to search for a triangle with
              // alpha and beta to be equal and variable, set xfix (the third parameter
              // to the bisect function to -2.0
              xfix = -2.0;
              // searching range between 0 and maximal value
              x1 = 0;
              x2 = 1;
              subcase = 0;
              break; // go to fine search
            } else {
              // here we need to check if our dY would lead to a new list entry
              // after the current one or if we have to step further backwards to find
              // a fitting place
              // Idea:
              //  everettEl[actLength-1] corresponds to the increment of the polarization, i.e. dP_old
              //  if abs(dP_new) < abs(dP_old), the new change in polarization (dP_new) could be
              //  achieved by inserting a new extremum to the end of the list
              //  if abs(dP_new) > abs(dP_old), the new input would wipe out older entries as the space is
              //  not large enougth; in this case, we go back by 2 positions (to get the same extremum again)
              //  and check again
              // Problem:
              //  aboves procedure only works for updates dP; however, we have dY = dP + eps*dX given
              //  if eps is negligble (i.e. for electrostatics with eps = eps0) this normally is no issue
              //  in case of magnetics and large values of dX, the contribution eps*dX might be significant, though
              //  > we would have to compute dP_new first but this would require the knowledge of dX_new (which of course is unknonw)
              //  > instead we compute dY_old = dP_old * eps*dX_old and compare dY_new with dY_old
              //    by aboves steps; at this point we can no longer compare the areas directly (i.e. abs(dP_new), abs(dP_old))
              //    this leads to problems like the following:
              //      material in saturation, dP (from 0 to saturation) = 1, eps*dX = 1 > dY_old = 2
              //      material shall now have a negative y-value, e.g. -0.1 > dY_new = -0.1-2 = -2.1
              //      > the value -0.1 is larger than neg. saturation (here -1) and thus should lead to a
              //        a corresponding x_new that is a minimum after the previous maximum of +1 (=xSat)
              //      > if we check for dY_old and dY_new, we see however, that dY_new does not fit into the space
              //        of dY_old (which obviously is wrong)
              //    > solution idea: take the first everett entry twice in terms of space computation
              //              this seems legit as the full preisach plane can swtich from +1 to -1 and thus
              //              has an effetive space of 2
              //

              // the availableSpace consists of the space for dP, i.e. the size of the
              // everett pixel
              Double availSpace = dHyst_to_dY*everettEl[actLength-1];

              Double dX = stringEl[actLength-1];
              Double dAnhyst = evalAnhystPart_normalized(stringEl[actLength-1]);
  //            std::cout << "stringEl[actLength-1]: "<< stringEl[actLength-1] << std::endl;
  //            std::cout << "actLength: "<< actLength << std::endl;
              if(actLength >= 2){
                //std::cout << "stringEl[actLength-2]: "<< stringEl[actLength-2] << std::endl;
                dX -= stringEl[actLength-2];
                dAnhyst -= evalAnhystPart_normalized(stringEl[actLength-2]);
                subcase = 2;
              } else {
                // last entry; everettEl counts twice
                availSpace += dHyst_to_dY*everettEl[actLength-1];
                // we also have to take the difference towards the other edge of the triangle, so add last entry another time
                // dX += stringEl[actLength-1];
                // > instead: add dX*dX_to_dY to availspace instead of increasing dX!
                // dX is simply stringEL here
                availSpace += dX_to_dY*dX + dAnhyst;
                subcase = 1;
                //
              }
              availSpace += dX_to_dY*dX + dAnhyst;

  //            std::cout << "everettEl[actLength-1]: " << everettEl[actLength-1] << std::endl;
  //            std::cout << "dX_to_dY*dX: " << dX_to_dY*dX << std::endl;
  //            std::cout << "abs(everettEl[actLength-1] + dX_to_dY*dX): " << abs(everettEl[actLength-1] + dX_to_dY*dX) << std::endl;
  //            std::cout << "abs(dY): " << abs(dY) << std::endl;
              // if( abs(everettEl[actLength-1] + dX_to_dY*dX) >= abs(dY) ){
              if( abs(availSpace) >= abs(dY) ){
                // std::cout << "enough space" << std::endl;
                // > enough space

                // fix edge of the everett triangle
                xfix = stringEl[actLength-1];
                if(actLength >= 2){
                  // reduce searching range (if possible due to previous entries)
                  if(minmaxcur > 0){
                    // entry at position actLength-2 is a maximum again and
                    // marks the upper bound for the everett triangle
                    x2 = stringEl[actLength-2];
                  } else {
                    // entry at position actLength-2 is a minimum again and
                    // marks the left bound for the everett triangle
                    x1 = stringEl[actLength-2];
                  }
                }

                break; // go to fine search
              } else {
                // std::cout << "not enough space" << std::endl;

                // not enough space; decrease list by 2 (it has to end with the
                // same minmaxtype again!) and check again
                // the removed entries have to be added to dY (as the difference
                // to the last computed value increases by the removed entries)
                // and adapt x1,x2 if possible

                // everettEl = dP
                // dY = dP + dX_to_dY*dX
                // we actually try to find a fitting dP here, so we must subtract
                // the influence of dX
                dY += dHyst_to_dY*everettEl[actLength-1]+dX_to_dY*dX + dAnhyst;
                actLength--;
                if(actLength >= 1){
                  dX = stringEl[actLength-1];
                  dAnhyst = evalAnhystPart_normalized(stringEl[actLength-1]);
                  if(actLength >= 2){
                    dX -= stringEl[actLength-2];
                    dAnhyst -= evalAnhystPart_normalized(stringEl[actLength-2]);
                  }
                  dY += dHyst_to_dY*everettEl[actLength-1] + dX_to_dY*dX + dAnhyst;
                  if(minmaxcur > 0){
                    // max> max was removed; lower bound can be increased to removed entry
                    x1 = stringEl[actLength-1];
                  } else {
                    // min> min was removed; right bound can be decreased to removed entry
                    x2 = stringEl[actLength-1];
                  }
                  LOG_DBG(scalpreisachInversion) << "Updated dY: " << dY;
                  actLength--;
                }
              }
            }
          } // while loop

          /*
           * thrid step: fine search via bisection (xOut in range -1,1)
           */
          Xout = bisect(dY,x1,x2,xfix,eps_mu,tol);
        } // reuse old value

        Xout *= XSaturated_;
      } // pos/neg/no saturation
      // rescale to -xSat to +xSat
      //LOG_TRACE(scalpreisachInversion) << "Found Xout (normalized): " << Xout;
      LOG_TRACE(scalpreisachInversion) << "Found Xout: " << Xout;

    } // anhyst only

		/*
     * final step: if overwrite == true > compute forward step to
     *						 update the list (this was not done yet!)
     * 
     * > has to be done for ALL cases
     */	
    bool debug = !true;
    successFlag = 4;
    if(overwrite || debug){
      LOG_DBG(scalpreisachInversion) << "overwrite? " << overwrite;
      int successFlagForward = 0;
      Double yRetrieved = computeValueAndUpdate( Xout, idx, overwrite, successFlagForward );
      //std::cout << "Xout: " << Xout << std::endl;
      //std::cout << "pRetrieved: " << yRetrieved << std::endl;
			LOG_DBG(scalpreisachInversion) << "Found Xout: " << Xout;
      yRetrieved+=Xout*eps_mu;
      LOG_TRACE(scalpreisachInversion) << "yRequested: " << Yin;
      LOG_TRACE(scalpreisachInversion) << "Found Yout: " << yRetrieved;
      LOG_TRACE(scalpreisachInversion) << "InversionCase: " << invcase;
      LOG_TRACE(scalpreisachInversion) << "SubCase: " << subcase;
      if(abs(yRetrieved-Yin) > tol){
        successFlag = -1;
      }
//        LOG_TRACE(scalpreisachInversion) << "Difference: " << yRetrieved-Yin;
//
//        Double P_old_normalized = previousPval_[idx];
//        Double Y_old_normalized = P_old_normalized+dX_to_dY*previousXval_[idx];
//
//        LOG_TRACE(scalpreisachInversion) << "yOld: " << Y_old_normalized*PSaturated_;
//        LOG_TRACE(scalpreisachInversion) << "yRequested-yOld: " << Yin-Y_old_normalized*PSaturated_;
//        LOG_TRACE(scalpreisachInversion) << "yRequested-yOld (normalized): " << Yin-Y_old_normalized;
//        LOG_TRACE(scalpreisachInversion) << "x1, x2, xOut: " << x1 << ", " << x2 << ", " << Xout;
//      }
    } else if(checkInversionResult_){
      int successFlagForward = 0;
      bool abortOnFail = false;
      Double yRetrieved = computeValueAndUpdate( Xout, idx, false, successFlagForward );
      yRetrieved+=Xout*eps_mu;
      if(abs(yRetrieved-Yin) > tol){
        std::stringstream exceptionmsg;
        exceptionmsg << "Everett-based-Inversion failed final test; remaining residual norm wrt y: " << abs(yRetrieved-Yin);
        if(abortOnFail){
          EXCEPTION(exceptionmsg.str());
        } else {
          //WARN(exceptionmsg.str());
          std::cout << exceptionmsg.str() << std::endl;
        }
      }

    } else {
      //store value for the case that we want to reuse it later
      // only needed if computeValueAndUpdate is not exectued with overwrite = true
      //NOTE: we must not execute this line before computeValueAndUpdate as this
      // function would return without overwriting the hyst memory
      //previousXval_[idx] = Xout/XSaturated_;
    }

		return Xout;
	}
  
  //  
  //  Double Preisach::computeValue(Double& Xin, Integer idx, bool overwrite) 
  //  {
  //
  //    /*
  //     * What is this function used for?
  //     * It does not update the list of minima and maxima but only evaluates
  //     * the everett function using the current list checking three cases:
  //     * a) new input wipes out complete list -> only one Everett pixel needed
  //     * b) new input replaces the currently last entry of list
  //     * c) new input attached to end of list
  //     *
  //     * What about cases in between?
  //     * No memory set.
  //     */
  //
  //    Vector<Double> &stringEl = strings_[idx];
  //    UInt& actLength = StringLength_[idx];
  //
  //    //normalize input
  //    Double newX, Yval;
  //    newX = normalizeInput(Xin);
  //
  //    Yval = 0.0;
  //    if ( abs(abs(newX)+eps_) > abs(stringEl[0]) || actLength == 0 ) {
  //      Yval =  everettPixel( -newX, newX );
  //    }
  //    else {
  //      Yval =  everettPixel(-stringEl[0],stringEl[0]);
  //      if ( abs(abs(newX)+eps_) > abs(stringEl[actLength-1]) ) {
  //        for ( UInt i=0; i<actLength-2; i++ ) {
  //          Yval +=  2.0*everettPixel(stringEl[i],stringEl[i+1]);
  //        }
  //        Yval +=  2.0*everettPixel(stringEl[actLength-2], newX);
  //      }
  //      else {
  //        for ( UInt i=0; i<actLength-1; i++ ) {
  //          Yval +=  2.0*everettPixel(stringEl[i],stringEl[i+1]);
  //        }
  //        Yval +=  2.0*everettPixel(stringEl[actLength-1], newX);
  //      }
  //    }
  //
  //    return ( Yval*PSaturated_ );
  //  }
  
  Vector<Double> Preisach::computeValue_vec(Vector<Double>& xVal, Integer idxElem, bool overwrite,
      bool debugOut, int& successFlag, bool skipAnhystPart) {
    
//    std::cout << "Scal Preisach: computeValue_vec; input: " << xVal.ToString() << std::endl;
    // apply the same steps as in CoefFunctionHyst
    Double scalInput;
    Double scalOutput;
    fixDirection_.Inner(xVal,scalInput);

    scalOutput = this->computeValueAndUpdate(scalInput,idxElem, overwrite,successFlag);
    Vector<Double> outputOfHystOperator = Vector<Double>(fixDirection_.GetSize());
    outputOfHystOperator.Init();
    outputOfHystOperator.Add(scalOutput,fixDirection_);

    return outputOfHystOperator;
  }

  Double Preisach::computeValueAndUpdate( Double Xin, Integer idx,
          bool overwrite, int& successFlag )
  {
    //do the deletion
    Double newY = 0.0;
    if(anhystOnly_ == false){
      // update minmaxlist will add anhyst part
      newY += updateMinMaxList(Xin, idx, overwrite, successFlag);
    } else  {
      Double X_norm_unclipped = Xin / XSaturated_;
      LOG_DBG(scalpreisachInversion) << "Eval Preisach - Add anhystPart " << evalAnhystPart_normalized(X_norm_unclipped);
      LOG_DBG(scalpreisachInversion) << "for X/norm = " << (X_norm_unclipped);
      newY += evalAnhystPart_normalized(X_norm_unclipped);
      successFlag = 1; // anhyst Only
    }
    return ( newY*PSaturated_ );
  }

  Double Preisach::computeValueAndUpdateMeasure( Double Xin, Integer idx,
          bool overwrite, int& successFlag, Double& time )
  {
    Timer* timer = new Timer();
    Double startTime = timer->GetCPUTime();
    timer->Start();

    Double newY = computeValueAndUpdate( Xin, idx,overwrite, successFlag );

    timer->Stop();
    Double endTime = timer->GetCPUTime();
    time = endTime-startTime;
    
    return newY;
  }
  
  

  Double Preisach::updateMinMaxList(Double Xin, Integer idx, 
          bool overwrite, int& successFlag )
  {
    LOG_TRACE(scalpreisachInversion) << "UpdateMinMaxList - Input: " << Xin;
    //std::cout << "UpdateMinMaxList - Input: " << Xin << std::endl;
    Double newY;
    
    //normalize input
    Double newX = normalizeAndClipInput(Xin);
    
    Vector<Double> &stringEl     = strings_[idx];
    Vector<Double> &helpStringEl = helpStrings_[idx];
		
    UInt& actLength = StringLength_[idx];
    UInt stringLength = actLength;
		
		// determine type of current input
		// only relevant if overwrite is true

    // NOTE: previousXval_[idx] = Xin/XSaturated_;
    // > previousXval_ is unclipped
    // > compare with clipped version of previousXval (better for saturation
    //    as > sat and >> sat will reuse max value
//    Double oldX = normalizeAndClipInput(previousXval_[idx]*XSaturated_);
//		Double diff = newX - oldX;


    Double diff = Xin-previousXval_[idx]*XSaturated_;
		int minmaxcur = 0;
		if(diff > 0){
			// new input is larger than last one > leads to a maximum
			minmaxcur = 1;
		} else if(diff < 0){
			// new input is smaller than last one > leads to a minimum
			minmaxcur = -1;
		} else if(diff == 0){
      // reuse old value but set previousXval anew
      // reason: above we compare the clipped values, Xin might be different from prevX
      LOG_TRACE(scalpreisachInversion) << "Reuse: " << preisachSum_[idx];
      successFlag = 0;
      if(overwrite){
        previousXval_[idx] = Xin/XSaturated_;
      }
      // note: preisachSum contains hyst part + anhyst part and is normalized to PSaturated_
      //  if anhystpart counts for PSaturated, preisachSum will be in range [-1,1] if input is below saturation
      //  if anhystpart does not count for PSaturated, preisachSum may be larger than 1 at saturation
			return preisachSum_[idx];
		}
    
    //    std::cout << "Element: " << idx << std::endl;
    //
    //    std::cout << "Input Xin: " << Xin << std::endl;
    
    //    std::cout << "Print out entries of min/max list before update" << std::endl;
    //    for ( UInt i=0; i<stringLength; i++ ) {
    //      std::cout << "index " << i << ": " << stringEl[i] << std::endl;
    //    }
    //    std::cout << "#############" << std::endl;
    
    //    std::cout << "Print out entries of helper min/max list before update" << std::endl;
    //    for ( UInt i=0; i<stringLength; i++ ) {
    //      std::cout << "index " << i << ": " << helpStringEl[i] << std::endl;
    //    }
    //    std::cout << "#############" << std::endl;
    
    if ( abs(newX) > abs(abs(stringEl[0]) - tol_) || stringLength == 0 ) {
      // reset
      successFlag = 4;
      stringLength = 1;
      if ( overwrite ){
				stringEl[0]  = newX;
        // at the beginning we check if input became larger or smaller and
        // according to that, we set minmaxcur to -1 (for min) or +1 (for max)
        // this type of extremum is correct if a new entry is added to the list
        // but might be wrong if it wipes out the list (e.g. if we stay in neg
        // saturation)
        // > reset minmaxcur here
        if(newX < 0){
          minmaxcur = -1;
        } else {
          minmaxcur = 1;
        }
				minmaxtype_[idx][0] = minmaxcur;
			} else {
        helpStringEl[0] = newX;
			}
    }
    else { 
      // check difference to last entry of list
      if ( abs( newX - stringEl[stringLength-1] ) > tol_ ) {
        // list modified
        successFlag = 5;
        helpStringEl[0] = -stringEl[0];
        for ( UInt i=1; i<=stringLength; i++ ) {
          helpStringEl[i] = stringEl[i-1];
        }
        
        UInt k = 0;
        
        Double a = helpStringEl[stringLength-1];
        Double b = helpStringEl[stringLength];
        
        //      std::cout << "a=" << a << "  b=" << b << std::endl;
        
        while ( ( k<stringLength-1) && 
                ( ( newX<=std::min(a,b) ) || ( newX>=std::max(a,b) ) ) ) {
          k = k + 1;
          a = helpStringEl[stringLength-k-1];
          b = helpStringEl[stringLength-k];
        }
        //        std::cout << "Old string length: " << stringLength << std::endl;
        stringLength = stringLength - k + 1;
        //        std::cout << "New string length: " << stringLength << std::endl;
        
        if (overwrite ) {
          //check, if capacity of string-arrays is too less
          if ( stringLength > maxStringLength_ ) {
            //resize the string-arrays
            maxStringLength_ += (UInt) round( (Double)maxStringLength_ / 2.0 );
            stringEl.Resize(maxStringLength_);
            
            //store the resulting strings
            for ( UInt i=0; i<stringLength-1; i++ ) {
              stringEl[i] = helpStringEl[i+1];
            }
            
            // resize help-String-array 
            helpStringEl.Resize(maxStringLength_+1);
          }
          else {
            //store the resulting strings
            for ( UInt i=0; i<stringLength-1; i++ ) {
              stringEl[i] = helpStringEl[i+1];
            }
          }
        }
        
        if ( overwrite ) {
          //store the new input
          stringEl[stringLength-1] = newX;
					// we only have to store the minmax type of the newly added entry
					// the others are already set (or will be overwritten)
					minmaxtype_[idx][stringLength-1] = minmaxcur;
        } else {
          
          //          std::cout << "Print out entries of helper min/max list before resorting" << std::endl;
          //
          //           for ( UInt i=0; i<stringLength+1; i++ ) {
          //             std::cout << "index " << i << ": " << helpStringEl[i] << std::endl;
          //           }
          //           std::cout << "#############" << std::endl;
          
          //correct storage in helpString
          for ( UInt i=0; i<stringLength-1; i++ ) {
            helpStringEl[i] = helpStringEl[i+1];
          }
          
          //          std::cout << "Print out entries of helper min/max list after resorting" << std::endl;
          //
          //           for ( UInt i=0; i<stringLength+1; i++ ) {
          //             std::cout << "index " << i << ": " << helpStringEl[i] << std::endl;
          //           }
          //           std::cout << "#############" << std::endl;
          
          helpStringEl[stringLength-1] = newX; 
        }
      }
      else {
        // append
        successFlag = 6;
        /*
         * Note: We check if a new input is large/small enough to modify the
         * actual storage stringEl;
         * if this is the case, we update
         * either stringEl (if overwrite = true) or helpStringEl (if overwrite = false)
         * and everything works just fine
         *
         * if this is not the case, we do neither change stringEl nor helpStringEl
         * if overwrite = true, this is no problem, as we evaluate stringEl which is
         * up to date (otherwise an update would have been needed)
         * if overwrite = false, we evaluate helpStringEl. This can be a serious problem
         * if helpStringEl is not stringEL
         * -> if no update of list was performed, copy stringEl (which is up to date)
         * to helpStringEl
         */
        //  std::cout << "No update needed; copy current stringEl to helpStringEl" << std::endl;
        for ( UInt i=0; i<=stringLength-1; i++ ) {
          helpStringEl[i] = stringEl[i];
        }
      }
    }
    
    if ( overwrite ) {
      actLength = stringLength;
      
      //      std::cout << "Overwrite = true" << std::endl;
      //      std::cout << "Print out entries of min/max list" << std::endl;
      //      for ( UInt i=0; i<stringLength; i++ ) {
      //        std::cout << "index " << i << ": " << stringEl[i] << std::endl;
      //      }
      //      std::cout << "#############" << std::endl;
      
      //      std::cout << "Updated non-temporary min/max list" << std::endl;
      //      std::cout << "Print out entries of min/max list" << std::endl;
      //      for ( UInt i=0; i<actLength; i++ ) {
      //        std::cout << "index " << i << ": " << stringEl[i] << std::endl;
      //      }
      //      std::cout << "#############" << std::endl;
      
      //compute preisach-sum
      Double pixelToAdd = everettPixel(-stringEl[0],stringEl[0]);
			evaluatedEverettPixel_[idx][0] = pixelToAdd;
      preisachSum_[idx] = pixelToAdd;
      
      for ( UInt i=0; i<actLength-1; i++ ) {
        pixelToAdd = 2.0*everettPixel(stringEl[i],stringEl[i+1]);
				evaluatedEverettPixel_[idx][i+1] = pixelToAdd;
        preisachSum_[idx] += pixelToAdd;   
      }

			
      // add optional anhysteretic part
      // > see e.g. "A preisach-based hysteresis model for magnetic and ferroelectric hysteresis" - A. Sutor 2010
      Double X_norm_unclipped = Xin / XSaturated_;
      LOG_DBG(scalpreisachInversion) << "Eval Preisach - Add anhystPart " << evalAnhystPart_normalized(X_norm_unclipped);
      LOG_DBG(scalpreisachInversion) << "for X/norm = " << (X_norm_unclipped);

      // 1. scale Preisach sum (till now only hysteretic part) by hystSaturated_
      // 2. add PSaturated_ * anhystPart
      // 3. normalize sum to PSaturated_

//      std::cout << "preisachSum_[idx]: " << preisachSum_[idx] << std::endl;

      preisachSum_[idx] *= hystSaturated_;
//      std::cout << "preisachSum_[idx]*= hystSaturated_: " << preisachSum_[idx] << std::endl;
      preisachSum_[idx] += PSaturated_*evalAnhystPart_normalized(X_norm_unclipped);
//      std::cout << "preisachSum_[idx]+= PSaturated_*evalAnhystPart_normalized(X_norm_unclipped): " << preisachSum_[idx] << std::endl;

      preisachSum_[idx] /= PSaturated_;
//      std::cout << "preisachSum_[idx]/= PSaturated_: " << preisachSum_[idx] << std::endl;

      //newY += anhyst_A_*std::atan(anhyst_B_*X_norm_unclipped) + anhyst_C_*X_norm_unclipped;
      newY = preisachSum_[idx];

			// store values for next evaluation
			// ONLY in case of overwrite
      // store the normalized values but Xin has to be unclipped!
      // (newX will be clipped to 1 if X > XSaturated)
			previousXval_[idx] = Xin/XSaturated_;
			previousPval_[idx] = newY;
			
    }
    else {
      /*
       * shouldn't we start at idx=1 here as helpStringEl[0] = -stringEl[0], helpStringEl[1] = stringEl[0]
       * by this we start with everettPixel(+stringEl[0],-stringEl[0])
       * then add 2*everettPixel(-stringEl[0],+stringEl[0])
       * ???
       * -> No, see line 181 -> elements get shifted to the right by 1
       *
       */
      
      //       std::cout << "Overwrite = false" << std::endl;
      //       std::cout << "Print out entries of helper min/max list" << std::endl;
      //
      //       for ( UInt i=0; i<stringLength; i++ ) {
      //         std::cout << "index " << i << ": " << helpStringEl[i] << std::endl;
      //       }
      //       std::cout << "#############" << std::endl;
      Double pixelToAdd = everettPixel(-helpStringEl[0], helpStringEl[0]);
      newY = pixelToAdd;
      for ( UInt i=0; i<stringLength-1; i++ ) {
				pixelToAdd = 2.0*everettPixel(helpStringEl[i],helpStringEl[i+1]);
        newY +=  pixelToAdd;
      }

      Double X_norm_unclipped = Xin / XSaturated_;
      LOG_DBG(scalpreisachInversion) << "Eval Preisach - Add anhystPart " << evalAnhystPart_normalized(X_norm_unclipped);
      LOG_DBG(scalpreisachInversion) << "for X/norm = " << (X_norm_unclipped);

      // 1. scale Preisach sum (till now only hysteretic part) by hystSaturated_
      // 2. add PSaturated_ * anhystPart
      // 3. normalize sum to PSaturated_
      newY *= hystSaturated_;
      newY += PSaturated_*evalAnhystPart_normalized(X_norm_unclipped);
      newY /= PSaturated_;

//      newY += evalAnhystPart_normalized(X_norm_unclipped);
      //newY += anhyst_A_*std::atan(anhyst_B_*X_norm_unclipped) + anhyst_C_*X_norm_unclipped;
    }
    LOG_TRACE(scalpreisachInversion) << "Computed new value: " << newY;
    //std::cout << "UpdateMinMaxList - Output: " << newY << std::endl;

//    if(overwrite){
//      successFlag = 2;
//    } else {
//      successFlag = 3;
//    }

    return newY;
  }
  
  
  Double Preisach::everettPixel(Double val1, Double val2)
  {
    
    Double X1 = std::max(val1,val2);
    Double X2 = std::min(val1,val2);
    
    UInt M = preisachWeights_.GetNumRows();
    Double delta = 2.0 / ( (Double) M );
    
    //  std::cout << "delta: " << delta << std::endl;
    //compute index for X1 (alpha)
    Integer idx1 = -1;
    Double alpha = -1.0;
    while ( alpha <= X1 ) {
      idx1++;
      alpha += delta;
    }
    
    if (alpha > 1.0) {
      alpha = 1.0;
      idx1 = M-1;
    }
    
    //compute index for X2 (beta)
    Integer idx2 = -1;
    Double  beta = -1.0;
    while ( beta <= X2 ) {
      idx2++;
      beta += delta;
    }
    beta -= delta;
    
    // X1 >= X2
    // -> mit alpha_0 = beta_0 = -1 und gleichem delta -> idx1 >= idx2
    
    //    std::cout << "idx1: " << idx1 << std::endl;
    //    std::cout << "idx2: " << idx2 << std::endl;
    
    Double area = 0.0;
    if ( idx1 >= 0 ) {
      UInt start = std::min(idx1,idx2);
      UInt stop  = std::max(idx1,idx2);
      
      // why do we iterate over the whole square and not over the triangle?
      for ( UInt i=start; i<=stop; i++ ) {
        for ( UInt j=start; j<=stop; j++ ) {
          area += preisachWeights_[i][j];
        }
      }
      
      area *= 0.5*delta*delta;
      
      //reduce the computed area
      Double diffX1 = alpha - X1;
      Double diffX2 = X2    - beta;
      Double minusArea;
      
      /*
       * The idea behind the minusArea is to handle input variation which are smaller than delta
       * Assume for example a min/max list with differences between mins and maxs < delta
       * In that case, the same preisach element would flip from +1 to -1 and the single steps would cancel
       * out in the sum. By reducing the weighting areas accordingly, we can also treat input lists of the described form
       * as now the summed up terms are weighted with different areas!
       * Commit out the section below and you will see, that it does not work anymore!
       */
      
      //check, if we are already on the diagonal!!
      if ( idx1 == idx2 ) {
        minusArea = (   diffX2 * (delta - 0.5*diffX2)
                + diffX1 * (delta - 0.5*diffX1)
                - diffX1*diffX2
                ) * preisachWeights_[idx1][idx2];
      }
      else {
        minusArea = ( (diffX1+diffX2 )*delta - diffX1*diffX2 )
                * preisachWeights_[idx1][idx2];
        
        Integer idx = idx1-1;
        while ( idx > idx2 ) {
          minusArea += diffX2 * delta * preisachWeights_[idx][idx2];
          idx--;
          //	  std::cout << "minusArea2=" << minusArea << std::endl;
        }
        minusArea += ( delta*diffX2 - 0.5*diffX2*diffX2 )
                * preisachWeights_[idx][idx2];
        
        idx = idx2 + 1;
        while ( idx < idx1 ) {
          minusArea += diffX1 * delta * preisachWeights_[idx1][idx];
          idx++;
        }
        minusArea += ( delta*diffX1 - 0.5*diffX1*diffX1 )
                * preisachWeights_[idx1][idx];
      }
      
      area -= minusArea;
    }
    
    //sgn-function
    if ( val2 < val1 ) {
      area *= -1.0;
    }
    
    return area;
  }
  
  // new name that better describes function
  // Xin gets not only normalized by dividing it with Xsaturated
  // but also gets clipped to +/-1 if Xsaturated is exceeded
  Double Preisach::normalizeAndClipInput(Double Xin)
  {
    
    Double Xout;
    
    if ( Xin > XSaturated_ ) {
      //saturation achieved!!
      Xout = 1.0;
    }
    else if ( Xin < -XSaturated_ ) {
      //saturation achieved!!
      Xout = -1.0;
    }
    else {
      //normalize input
      Xout = Xin / XSaturated_;
    }
    
    return Xout;
  }
  
}
