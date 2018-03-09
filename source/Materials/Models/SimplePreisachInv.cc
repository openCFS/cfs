#include <iostream>
#include <fstream>

#include "SimplePreisachInv.hh"

namespace CoupledField
{ 
  
  SimplePreisachInv::SimplePreisachInv(Integer numElem, 
          Double xSat, 
          Double ySat, 
          Matrix<Double>& preisachWeight, 
          bool isVirgin) 
  : Hysteresis(numElem) {
    
    eps_ = 1e-5;
    nu0_ = 7.9577e+05;
    
    if (xSat > 0 ) {
      Xsaturated_  = xSat;
    }
    else {
      Xsaturated_  = 1.0;
    }
    
    Ysaturated_ = ySat;
    isVirgin_   = isVirgin;
    
    //    std::cout << "xSat: " << Xsaturated_ << "  Ysat=" << YSaturated_ << std::endl;
    preisachWeights_ = preisachWeight;
    
    previousYval_.Resize(numElem);
    previousYval_.Init();
    
    actXval_.Resize(numElem);
    actXval_.Init();
    
    preisachSum_.Resize(numElem);
    preisachSum_.Init(0);
    
    StringLength_.Resize(numElem);
    
    strings_     = new Vector<Double>[numElem];
    helpStrings_ = new Vector<Double>[numElem];
    
    maxStringLength_ = 100;
    
    for (Integer el=0; el<numElem; el++) {
      strings_[el].Resize(maxStringLength_);
      helpStrings_[el].Resize(maxStringLength_+1);
      
      StringLength_[el] = 1;
      for ( UInt i=0; i<maxStringLength_; i++) {
        strings_[el][i] = 0.0;
      }
    }
    
    //compute grid points
    M_ = preisachWeight.GetNumRows(); 
    epts_.Resize(M_+1);
    Double dx = 2.0 / (Double)M_; 
    for (UInt i=0; i<=M_; i++ ) {
      epts_[i] = -1.0 + (Double)i * dx;
    }
    std::cout << "epts: \n" << epts_.ToString() << std::endl << std::endl;
    
    //    computePreisachWeights();
  }
  
  SimplePreisachInv::~SimplePreisachInv()
  {
    delete [] strings_;
    delete [] helpStrings_;
  }
    
  Double SimplePreisachInv::computeValue(Double& Yin, 
          Integer idxEl, 
          bool overwrite) 
  {
    /*
     * mnierla:
     * Quelle > Dissertation Felix Wolf?
     * > wohl eher nein, da zu unterschiedlich
     * 
     */
    
    Vector<Double> &stringEl = strings_[idxEl];
    UInt& actLength = StringLength_[idxEl];
    
    //normalize input
    Double newY, newX, dY, yStart, xStart;
    UInt j;
    
    //     if ( Yin > Ysaturated_ ) {
    //       newX = Xsaturated_ + nu0_ * (Yin - Ysaturated_);
    //       newX /= Xsaturated_;
    //     }
    //     else {
    
    std::cout << "Yin: " << Yin << std::endl;
    
    newY = normalizeYval(Yin);
    /*
     * previousYval_ wird nie gesetzt
     */
    dY   = newY - previousYval_[idxEl];
    
    std::cout << "newY: " << newY << std::endl;
    std::cout << "dY = " << dY << std::endl;
    std::cout << "actLength: " << actLength << std::endl;
    if ( (actLength <= 1)  && (dY > 1e-10) ) {
      std::cout << "Branche 1: ComputeInv: Yin: " << newY << std::endl;
      std::cout << "Branche 1: ComputeInv: Yprev: " << previousYval_[idxEl] << std::endl;
      j = 0;
      
      /*
       * only in Branche 1 we use everettPixel without the factor of 2
       * why?
       */
//      everettPixel( -1,+1) = +1
//      everettPixel( 1,-1) = -1  
//      evalEverett = everettPixel*ySat        

      // yStart = everettPixel( -1,+1) = ySat
      // yStart = everettPixel( 1,-1) = -ySat
      std::cout << "test everettPixel  " << std::endl;
      std::cout << "everettPixel( -1,+1) " << everettPixel( -1,+1) << std::endl;
      std::cout << "everettPixel( 1,-1) " << everettPixel( 1,-1) << std::endl;
      std::cout << "ySat: " << Ysaturated_ << std::endl;
      yStart = everettPixel( -epts_[j],epts_[j]);
      
      // yStart = -ySat, newY >= -ySat
      // > loop will be entered at least once
      // full Everett triangles will be evaluated, i.e.
      // the presiach plane must be in virgin state
      while ( yStart < (newY*0.999) ) {
        j += 1;
        yStart = everettPixel( -epts_[j],epts_[j]);
      }
      
      j -= 1;
      yStart = everettPixel( -epts_[j],epts_[j]);
      xStart = epts_[j];
      
      // linear interpolation between next larger everettPixel (j+1) and next smaller one (j) 
      newX = xStart + (2.0/(Double)M_) * ( newY - yStart ) / ( everettPixel( -epts_[j+1],epts_[j+1]) - yStart );
      
      if ( overwrite ) {
        actLength   = 1;
        stringEl[0] = newX;
      }
    }
    else {
      if ( actLength <= 1 ) {
        std::cout << "Branche 2,1: ComputeInv: Yin: " << newY << std::endl;
        
        Double yAdd, p1, p2;
        j = M_-1;
        yAdd =  2.0*everettPixel( stringEl[0], epts_[j] ); 
        while ( yAdd > dY ) {
          j -= 1;
          yAdd = 2.0*everettPixel( stringEl[0], epts_[j] ); 
        }
        xStart = epts_[j];
        p1 = 2.0*everettPixel( stringEl[0], epts_[j] );
        p2 = 2.0*everettPixel( stringEl[0], epts_[j+1] );
        if ( p2 > 0 ) {
          Double delta =  2.0*everettPixel( epts_[j+1], epts_[j] );
          newX = xStart + (2.0/(Double)M_) * (p1 - dY ) / delta;
        }
        else {
          newX = xStart + (2.0/(Double)M_) * (p1 - dY ) / (p1 -p2);
        }
      }
      else {
        Double yDiff, yDiffOld, p1, p2;
        j = M_-1;
        std::cout << "j " << j << std::endl;
        yDiffOld = 2.0*everettPixel( stringEl[0], stringEl[1] );
        yDiff = 2.0*everettPixel( stringEl[0], epts_[j] ) - yDiffOld; 
        std::cout << "yDiffOld: "<<yDiffOld<<std::endl;
        std::cout << "yDiff: "<<yDiff<<std::endl;
        std::cout << "stringEl[0]: " << stringEl[0] << std::endl;
        while ( yDiff > dY ) {
          j -= 1;
          std::cout << "j " << j << std::endl;
          std::cout << "epts_[j] " << epts_[j] << std::endl;
          yDiff = 2.0*everettPixel( stringEl[0], epts_[j] ) - yDiffOld; 
          std::cout << "yDiff: "<<yDiff<<std::endl;
        }
        xStart = epts_[j];
        p1 = 2.0*everettPixel( stringEl[0], epts_[j] );
        p2 = 2.0*everettPixel( stringEl[0], epts_[j+1] );
        if ( p2 > 0 ) {
          Double delta =  2.0*everettPixel( epts_[j+1], epts_[j] );
          newX = xStart + (2.0/(Double)M_) * (p1 - yDiffOld - dY ) / delta;
        }
        else {
          newX = xStart + (2.0/(Double)M_) * (p1 - yDiffOld - dY ) / (p1 -p2);
        }
      }
      if ( overwrite ) {
        /*
         * maximal 2 entries?
         */
        actLength   = 2;
        stringEl[1] = newX;
      }
    }
    
    std::cout << "ComputeInv: Xout: " << newX << std::endl;
    
    actXval_[idxEl] = newX;
    
    return ( newX*Xsaturated_ );
  }
  
  
  Double SimplePreisachInv::EvalEverett(Double xVal1, Double xVal2, Integer idx)
  {
    
    //normalize input
    Double X1 = normalizeXval( xVal1 );
    Double X2 = normalizeXval( xVal2 );
    
    
    //    std::cout << "StrLength: " <<  StringLenght_[idx] << std::endl;
    Double pPixel;
    if (  StringLength_[idx] < 2 ) {
      //    Double pPixel =  everettPixel(X1, X2);
      pPixel =  everett(-X2, X2);
    }
    else {
      pPixel =  everett(X1, X2);
    }
    
    std::cout << "X1=" << X1 << " X2=" << X2 << " pPixel=" << pPixel << std::endl;
    pPixel *= Ysaturated_;
    
    return pPixel;
  }
  
  
  Double SimplePreisachInv::everett(Double X1, Double X2)
  {
    
    //     Double newY;
    //     Double diffX = X2 - X1;
    
    //     if ( diffX > 0) {
    //       newY = 0.25*diffX*diffX;
    //     }
    //     else {
    //       newY = -0.25*diffX*diffX;
    //     }
    //return newY;
    
    Double pPixel =  everettPixel(X1, X2);
    return pPixel;
  }
  
  
  Double SimplePreisachInv::everettPixel(Double val1, Double val2)
  {
    
    Double X1 = std::max(val1,val2);
    Double X2 = std::min(val1,val2);
    
    UInt M = preisachWeights_.GetNumRows();
    Double delta = 2.0 / ( (Double) M );
    //     std::cout << "delta=" << delta << std::endl;
    //     std::cout << "e1=" << val1 << "  e2=" << val2 << std::endl;
    //     std::cout << "X1=" << X1 << "  X2=" << X2 << std::endl;
    
    // we compute the upper left corner within the preisach-domain
    // alpha > X1; beta < X2;
    
    //compute index for X1 (alpha)
    Integer idx1 = -1;
    Double alpha = -1.0;
    while ( alpha <= X1 ) {
      idx1++;
      alpha += delta;
      //std::cout << "idx1=" << idx1 << "  alpha=" << alpha << std::endl;
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
    
    std::cout << "idx1=" << idx1 << "  idx2=" << idx2 << std::endl;
    std::cout << "alpha=" << alpha << "  beta=" << beta << std::endl;
    
    Double area = 0.0;
    if ( idx1 >= 0 ) {
      UInt start = std::min(idx1,idx2);
      UInt stop  = std::max(idx1,idx2);
      for ( UInt i=start; i<=stop; i++ ) {
        for ( UInt j=start; j<=stop; j++ ) {
          area += preisachWeights_[i][j];
        }
      }
      
      area *= 0.5*delta*delta;
      
      //      std::cout << "areaToMuch=" << area << std::endl;
      
      //reduce the computed area
      Double diffX1 = alpha - X1;
      Double diffX2 = X2    - beta;
      //      std::cout << " diffX1=" << diffX1 << "  diffX2=" << diffX2 << std::endl;
      
      Double minusArea;
      
      //check, if we are already on the diagonal!!
      if ( idx1 == idx2 ) {
        minusArea = (   diffX2 * (delta - 0.5*diffX2) 
                + diffX1 * (delta - 0.5*diffX1)
                - diffX1*diffX2 
                ) * preisachWeights_[idx1][idx2];
        //	std::cout << "minusArea1end=" << minusArea << std::endl;
      }
      else {
        minusArea = ( (diffX1+diffX2 )*delta - diffX1*diffX2 ) 
                * preisachWeights_[idx1][idx2];
        
        //	std::cout << "minusArea1=" << minusArea << std::endl;
        Integer idx = idx1-1;
        while ( idx > idx2 ) {
          minusArea += diffX2 * delta * preisachWeights_[idx][idx2];
          idx--;
          //	  std::cout << "minusArea2=" << minusArea << std::endl;
        }
        minusArea += ( delta*diffX2 - 0.5*diffX2*diffX2 )
                * preisachWeights_[idx][idx2]; 
        
        // 	std::cout << "minusArea3=" <<  ( delta*diffX2 - 0.5*diffX2*diffX2 )
        // 	  * preisachWeights_[idx][idx2] << std::endl;
        
        idx = idx2 + 1;
        while ( idx < idx1 ) {
          minusArea += diffX1 * delta * preisachWeights_[idx1][idx];
          idx++;
        }
        minusArea += ( delta*diffX1 - 0.5*diffX1*diffX1 )  
                * preisachWeights_[idx1][idx];
      }
      
      
      std::cout << "Area: " << area << std::endl;
      std::cout << "minusArea: " << minusArea << std::endl;
      area -= minusArea; 
    }
    
    //sgn-function
    if ( val2 < val1 ) {
      //      std::cout << "val2 < val1" << std::endl;
      area *= -1.0; 
    }
    
    return area;
  }
  
  
  void SimplePreisachInv::computePreisachWeights()
  {
    
    UInt dim = 6;
    preisachWeights_.Resize(dim,dim);
    for ( UInt i=0; i<dim; i++) {
      for ( UInt j=0; j<dim; j++) {
        preisachWeights_[i][j] = 0.5;
      }
    }
  }
  
  Double SimplePreisachInv::normalizeXval(Double Xin)
  {
    
    Double Xout;
    
    if ( Xin > Xsaturated_ ) {
      //saturation achieved!!
      Xout = 1.0;
    }
    else if ( Xin < -Xsaturated_ ) {
      //saturation achieved!!
      Xout = -1.0;
    }
    else {
      //normalize input
      Xout = Xin / Xsaturated_;
    }
    
    return Xout;
  }
  
  Double SimplePreisachInv::normalizeYval(Double& Yin)
  {
    
    Double Yout;
    
    if ( Yin > Ysaturated_ ) {
      //saturation achieved!!
      Yout = 1.0;
      Yin = Ysaturated_;
    }
    else if ( Yin < -Ysaturated_ ) {
      //saturation achieved!!
      Yout = -1.0;
      Yin = -Ysaturated_;
    }
    else {
      //normalize input
      Yout = Yin / Ysaturated_;
    }
    
    return Yout;
  }
  
}
