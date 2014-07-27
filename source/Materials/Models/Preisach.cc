#include "Preisach.hh"

#include <iostream>
#include <fstream>



namespace CoupledField
{ 

  Preisach::Preisach(Integer numElem, Double xSat, Double ySat, 
		       Matrix<Double>& preisachWeight, bool isVirgin) 
    : Hysteresis(numElem)
  {

    eps_ = 1e-5;

    if (xSat > 0 ) {
      Xsaturated_  = xSat;
    }
    else {
      Xsaturated_  = 1.0;
    }

    YSaturated_  = ySat;
    isVirgin_    = isVirgin;

    preisachWeights_ = preisachWeight;

    preisachSum_.Resize(numElem);
    preisachSum_.Init(0);

    StringLenght_.Resize(numElem);

    strings_     = new Vector<Double>[numElem];
    helpStrings_ = new Vector<Double>[numElem];

    maxStringLength_ = 100;

    for (Integer el=0; el<numElem; el++) {
      strings_[el].Resize(maxStringLength_);
      helpStrings_[el].Resize(maxStringLength_+1);

      StringLenght_[el] = 1;
      for ( UInt i=0; i<maxStringLength_; i++) {
	strings_[el][i] = 0.0;
      }
    }

  }

  Preisach::~Preisach()
  {
    delete [] strings_;
    delete [] helpStrings_;

  }

  Double Preisach::getValue( Integer idx ) 
  {
    return ( preisachSum_[idx]*YSaturated_ );
  }

  Double Preisach::computeValue(Double& Xin, Integer idx, bool overwrite) 
  {

    Vector<Double> &stringEl     = strings_[idx];
    UInt& actLength = StringLenght_[idx];

    //normalize input
    Double newX, Yval;
    newX = normalizeInput(Xin);

    Yval = 0.0;
    if ( abs(abs(newX)+eps_) > abs(stringEl[0]) || actLength == 0 ) {
      Yval =  everettPixel( -newX, newX );
    }
    else {
      Yval =  everettPixel(-stringEl[0],stringEl[0]);
      if ( abs(abs(newX)+eps_) > abs(stringEl[actLength-1]) ) {
        for ( UInt i=0; i<actLength-2; i++ ) {
          Yval +=  2.0*everettPixel(stringEl[i],stringEl[i+1]);
        }
        Yval +=  2.0*everettPixel(stringEl[actLength-2], newX);
      }
      else {
        for ( UInt i=0; i<actLength-1; i++ ) {
          Yval +=  2.0*everettPixel(stringEl[i],stringEl[i+1]);
        }
        Yval +=  2.0*everettPixel(stringEl[actLength-1], newX);
      }
    }

    return ( Yval*YSaturated_ );
  }


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

    Double newY;

    //normalize input
    Double newX = normalizeInput(Xin);

    Vector<Double> &stringEl     = strings_[idx];
    Vector<Double> &helpStringEl = helpStrings_[idx];

    UInt& actLength = StringLenght_[idx];
    UInt stringLength = actLength;

    if ( abs(newX) > abs(abs(stringEl[0]) - eps_) || stringLength == 0 ) {
      stringLength = 1;
      if ( overwrite ) 
        stringEl[0]  = newX;
      else
        helpStringEl[0] = newX;
    }
    else { 
      if ( abs( newX - stringEl[stringLength-1] ) > eps_ ) {
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
        
        stringLength = stringLength - k + 1;
        
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
        }
        else {
          //correct storage in helpString
          for ( UInt i=0; i<stringLength-1; i++ ) {
              helpStringEl[i] = helpStringEl[i+1];
          }
          helpStringEl[stringLength-1] = newX; 
        }
      }
    }
    
    if ( overwrite ) {
      actLength = stringLength;

      //compute preisach-sum
      preisachSum_[idx] =  everettPixel(-stringEl[0],stringEl[0]);
      for ( UInt i=0; i<actLength-1; i++ ) {
        preisachSum_[idx] +=  2.0*everettPixel(stringEl[i],stringEl[i+1]);
      }
      newY = preisachSum_[idx]; 
    }
    else {
      newY = everettPixel(-helpStringEl[0], helpStringEl[0]);
      for ( UInt i=0; i<stringLength-1; i++ ) {
        newY +=  2.0*everettPixel(helpStringEl[i],helpStringEl[i+1]);
      }
    }

    return newY;
  }


  Double Preisach::everettPixel(Double val1, Double val2)
  {

    Double X1 = std::max(val1,val2);
    Double X2 = std::min(val1,val2);

    UInt M = preisachWeights_.GetNumRows();
    Double delta = 2.0 / ( (Double) M );

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

      //reduce the computed area
      Double diffX1 = alpha - X1;
      Double diffX2 = X2    - beta;
      Double minusArea;
 
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


  Double Preisach::normalizeInput(Double Xin)
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

}
