#include "Preisach.hh"

#include <iostream>
#include <fstream>



namespace CoupledField
{ 
  
  Preisach::Preisach(Integer numElem, Double xSat, Double ySat, 
          Matrix<Double>& preisachWeight, bool isVirgin) 
  : Hysteresis(numElem)
  {
    
    std::cout << "ScalarPreisach: Using Everett function" << std::endl;
    
    //tol_ = 1e-5;
    tol_ = 1e-15;
    
    if (xSat > 0 ) {
      XSaturated_  = xSat;
    }
    else {
      XSaturated_  = 1.0;
    }
    
    YSaturated_  = ySat;
    isVirgin_    = isVirgin;
    
    preisachWeights_ = preisachWeight;
    preisachSum_.Resize(numElem);
    preisachSum_.Init(0);
    StringLength_.Resize(numElem);
    
    strings_     = new Vector<Double>[numElem];
    helpStrings_ = new Vector<Double>[numElem];
    minmaxtype_     = new Vector<Integer>[numElem];
    evaluatedEverettPixel_ = new Vector<Double>[numElem];
    maxStringLength_ = 100;
    
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
    return ( preisachSum_[idx]*YSaturated_ );
  }
  
  Double Preisach::bisect(Double dY,Double xMin,Double xMax, Double xFixed, Double eps_mu, Double tol){
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
    Double dX_to_dY = eps_mu*XSaturated_/YSaturated_;
    Double dX = 0.0;
    
    if(xFixed >= -1){
      //std::cout << "bisection with fixed edge" << std::endl;
      for(UInt i = 0; i<maxIter; i++){
        xMiddle = 0.5*(xMax+xMin);
        dX = xMiddle - xFixed;
        // important: here we need a factor of 2.0 as we have an update
        // to the previous state (i.e. we have to revert from +1 to -1 or
        // vice versa)
        everettMiddle = 2.0*everettPixel(xFixed,xMiddle);
        // dY = dP + dX_to_dY * dX
        diffMiddle = dY - everettMiddle - dX*dX_to_dY;
        
        if(abs(diffMiddle) < tol){
          return xMiddle;
        }
        
        dX = xMin - xFixed;
        // syntax for everettPixel: 1. parameter: older/previous entry, 2. parameter new/next entry
        everettMin = 2.0*everettPixel(xFixed,xMin);
         
        diffMin = dY - everettMin - dX*dX_to_dY;
        
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
			return xMiddle;
    } else {
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
        // here we do not need a factor of 2.0 as we do not have to
        // revert a previous state
        everettMiddle = everettPixel(-xMiddle,xMiddle);

        diffMiddle = dY - everettMiddle - dX*dX_to_dY;
        
        if(abs(diffMiddle) < tol){
          return xMiddle;
        }
        
        // note that everettPixel always returns a positive value (the area)
        // and a sign that can be used to find out if this area is to be subtracted
        // or to be added
        dX = xMin; // here we have dX towards 0
        everettMin = everettPixel(-xMin,xMin);
        //everettMax = everettPixel(-xMax,xMax,sign);
				
        diffMin = dY - everettMin - dX*dX_to_dY;
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
			return xMiddle;
    }
  }
  
  Double Preisach::computeInputAndUpdate(Double Yin, Double eps_mu, Integer idx, bool overwrite){
    /*
     * NEW implementation; this shall replace class SimplePreisachInv that
     * apparently is not working
     * 
     * 
     * Source: "Generalisiertes Preisach-Modell für die Simulation uund Kompensation
     *         piezokeramischer Aktoren" - Dissertation, Felix Wolf, p. 127ff
     */
    
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
    Double tol = 1e-8;
    Double dX_to_dY = eps_mu*XSaturated_/YSaturated_;
    
    //std::cout << "Inversion" << std::endl;
    //std::cout << "Yin: " << Yin << std::endl;
    //std::cout << "YSaturated: " << YSaturated_ << std::endl;
    //std::cout << "eps_mu: " << eps_mu << std::endl;
    if(Yin >= (YSaturated_ + eps_mu*XSaturated_) ){
      //std::cout << "pos saturation" << std::endl;
      /*
       * System is actually in positive saturation
       * 
       * yIn = Ysaturated_ + eps_mu*xOut
       *     >= YSaturated_ + eps_mu*XSaturated_
       * 
       * xOut = (yIn - Ysaturated_)/eps_mu
       */
      Xout = (Yin - YSaturated_)/eps_mu;
    } else if (Yin <= (-YSaturated_ - eps_mu*XSaturated_) ){
      //std::cout << "neg saturation" << std::endl;
      /*
       * System is actually in negative saturation
       * 
       * yIn = -Ysaturated_ + eps_mu*xOut
       * 
       * xOut = (yIn + Ysaturated_)/eps_mu
       */
      Xout = (Yin + YSaturated_)/eps_mu;
    } else {
      /*
       * 1. Compute difference between requested yVal and previously computed value
       */
      // factor for relating a change in x to a change in y for
      // the non-hysteretic part
      //  > dY = dP + dX_to_dY*dX
      // for dY, dP and dX all normalized
      
      // TODO:; should be normalized AFTER subtracting x
      // do not cap to +/-1 here as Y is not the polarization (that is actually capped!)
      Double Y_normalized = Yin/YSaturated_;
      Double P_old_normalized = previousPval_[idx];
//      std::cout << "Y_normalized: " << Y_normalized << std::endl;
//      std::cout << "P_old_normalized: " << P_old_normalized << std::endl;
//      
      Double dY = Y_normalized - P_old_normalized;
//      std::cout << "Y_normalized - P_old_normalized: " << dY << std::endl;
//      std::cout << "previousXval_[idx]: " << previousXval_[idx] << std::endl;
//      std::cout << "dX_to_dY*previousXval_[idx]): " << dX_to_dY*previousXval_[idx] << std::endl;
//      
      dY = dY - dX_to_dY*previousXval_[idx];
      //std::cout << "Y_normalized - P_old_normalized - dX_to_dY*previousXval_[idx]: " << dY << std::endl;
      
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
      if(abs(dY) < tol){
        // difference is negligible
        //std::cout << "take previous xvalue" << std::endl;
        Xout = previousXval_[idx];
      } else{
        Vector<Double> &stringEl = strings_[idx];
				Vector<int> &minmaxEl = minmaxtype_[idx];
				Vector<Double> &everettEl = evaluatedEverettPixel_[idx];
        
				Double x1,x2; // search interval for later bisection
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
          
          dY += dP;
          //std::cout << "dY+dP: " << dY << std::endl;
          
          // we search space for dY not dP, so we have to add the contribution of
          // the reversible part
          Double dX = stringEl[actLength-1];
          if(actLength >= 2){
            dX -= stringEl[actLength-2];
          }
          //std::cout << "dX = " << dX << std::endl;
          dY += dX*dX_to_dY;
          //std::cout << "dY+dP+dX_to_dY*dX: " << dY << std::endl;
          
					if(minmaxcur > 0){
						// max> max was removed; lower bound can be increased to removed entry
						x1 = stringEl[actLength-1];
					} else {
						// min> min was removed; right bound can be decreased to removed entry
						x2 = stringEl[actLength-1];
					}
					
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
            Double availSpace = everettEl[actLength-1];
            
            Double dX = stringEl[actLength-1];
//            std::cout << "stringEl[actLength-1]: "<< stringEl[actLength-1] << std::endl;
//            std::cout << "actLength: "<< actLength << std::endl;
            if(actLength >= 2){
              //std::cout << "stringEl[actLength-2]: "<< stringEl[actLength-2] << std::endl;
              dX -= stringEl[actLength-2];
            } else {
              // last entry; everettEl counts twice
              availSpace += everettEl[actLength-1];
              // we also have to take the difference towards the other edge of the triangle, so add last entry another time
              dX += stringEl[actLength-1];
            }
            availSpace += dX_to_dY*dX;
            
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
							dY += everettEl[actLength-1]+dX_to_dY*dX;
							actLength--;
							if(actLength >= 1){
                dX = stringEl[actLength-1];
                if(actLength >= 2){
                  dX -= stringEl[actLength-2];
                }
								dY += everettEl[actLength-1]+dX_to_dY*dX;
								if(minmaxcur > 0){
									// max> max was removed; lower bound can be increased to removed entry
									x1 = stringEl[actLength-1];
								} else {
									// min> min was removed; right bound can be decreased to removed entry
									x2 = stringEl[actLength-1];
								}
								
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
      
      // rescale to -xSat to +xSat
      Xout *= XSaturated_;
		} // pos/neg/no saturation
		
		/*
     * final step: if overwrite == true > compute forward step to
     *						 update the list (this was not done yet!)
     * 
     * > has to be done for ALL cases
     */	
		//if(overwrite){
      computeValueAndUpdate( Xout, idx, overwrite );
      Double yRetrieved = computeValueAndUpdate( Xout, idx, overwrite );
      //std::cout << "Xout: " << Xout << std::endl;
      //std::cout << "pRetrieved: " << yRetrieved << std::endl;
      yRetrieved+=Xout*eps_mu;
      //std::cout << "yRetrieved: " << yRetrieved << std::endl;
      //std::cout << "yRequested: " << Yin << std::endl;
		//}
		
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
  //    return ( Yval*YSaturated_ );
  //  }
  
  Double Preisach::computeValueAndUpdate( Double Xin, Integer idx,
          bool overwrite )  
  {
    
    //do the deletion
    Double newY = updateMinMaxList(Xin, idx, overwrite);
    
    return ( newY*YSaturated_ );
  }
  
  
  Double Preisach::updateMinMaxList(Double Xin, Integer idx, 
          bool overwrite )
  {
    
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
		Double diff = newX - previousXval_[idx];
		int minmaxcur = 0;
		if(diff > 0){
			// new input is larger than last one > leads to a maximum
			minmaxcur = 1;
		} else if(diff < 0){
			// new input is smaller than last one > leads to a minimum
			minmaxcur = -1;
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
      if ( abs( newX - stringEl[stringLength-1] ) > tol_ ) {
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
        }
        else {
          
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
        /*
         * Note: We check if a new input is large/small enough to modify the
         * actual storage stringEl;
         * if this is the case, we update
         * either stringEl (if update = true) or helpStringEl (if update = false)
         * and everything works just fine
         *
         * if this is not the case, we do neither change stringEl nor helpStringEl
         * if update = true, this is no problem, as we evaluate stringEl which is
         * up to date (otherwise an update would have been needed)
         * if update = false, we evaluate helpStringEl. This can be a serious problem
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
    }
    
    //std::cout << "UpdateMinMaxList - Output: " << newY << std::endl;
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
