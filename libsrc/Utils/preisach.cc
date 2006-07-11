#include <iostream>
#include <fstream>

#include "preisach.hh"

namespace CoupledField
{ 

  Preisach :: Preisach(Integer numElem, Double xSat, Double ySat, 
		       Matrix<Double>& preisachWeight, bool isVirgin) 
    : Hysteresis(numElem)
  {
    ENTER_FCN("Preisach::Preisach" );

    if (xSat > 0 ) {
      Xsaturated_  = xSat;
    }
    else {
      Xsaturated_  = 1.0;
    }

    YSaturated_  = ySat;
    //    YRemnant_    = xRem;
    isVirgin_    = isVirgin;

    preisachWeights_ = preisachWeight;

    lastVal_.Resize(numElem);
    preisachSum_.Resize(numElem);
    preisachSum_.Init(0);

    StringLenght_.resize(numElem);

    strings_     = new std::vector<Double>[numElem];
    helpStrings_ = new std::vector<Double>[numElem];

    maxStringLength_ = 4;

    for (Integer el=0; el<numElem; el++) {
      strings_[el].resize(maxStringLength_);
      helpStrings_[el].resize(maxStringLength_+1);

      StringLenght_[el] = 1;
      for ( UInt i=0; i<maxStringLength_; i++) {
	strings_[el][i] = 0.0;
      }
    }

    //allocate memory for previous results, needed for the
    //effective material parameter formulation
    Xprevious_.Resize(numElem);
    Yprevious_.Resize(numElem);
    Xprevious_.Init(0.0);
    Yprevious_.Init(0.0);

    computePreisachWeights();
  }

  Preisach :: ~Preisach()
  {
    delete [] strings_;
    delete [] helpStrings_;

  }

  Double Preisach :: computeValue(Double Xin, Integer idxElem) 
  {
    ENTER_FCN( "Preisach::computeValue" );

    Integer idx = idxElem - 1;
  
    //normalize input
    Double newX = normalizeInput(Xin);

    //do the deletion
    updateMinMaxList(Xin, idxElem);

    return ( preisachSum_[idx]*YSaturated_ );
  }


  void Preisach :: updateMinMaxList(Double Xin, Integer nrEl)
  {
    ENTER_FCN( "Preisach::updateMinMaxList" );

    Integer idx = nrEl - 1;

    //normalize input
    Double newX = normalizeInput(Xin);

    std::vector<Double> &stringEl     = strings_[idx];
    std::vector<Double> &helpStringEl = helpStrings_[idx];

    UInt& actLength = StringLenght_[idx];

    Double eps = 1.0e-10;
    if ( abs(newX) > abs(stringEl[0] - eps) || actLength == 0 ) {
      actLength   = 1;
      stringEl[0] = newX;
    }
    else {
      helpStringEl[0] = -stringEl[0];
      for ( UInt i=1; i<=actLength; i++ ) {
	helpStringEl[i] = stringEl[i-1];
      }

      UInt k = 0;

      Double a = helpStringEl[actLength-1];
      Double b = helpStringEl[actLength];

      //      std::cout << "a=" << a << "  b=" << b << std::endl;

      while ( ( k<actLength-1) && 
	      ( ( newX<=std::min(a,b) ) || ( newX>=std::max(a,b) ) ) ) {
	k = k + 1;
	a = helpStringEl[actLength-k-1];
	b = helpStringEl[actLength-k];
      }

      actLength = actLength - k + 1;

      std::cout << "actLength= " << actLength << std::endl;
      //check, if capacity of string-arrays is too less
      if ( actLength > maxStringLength_ ) {
	//resize the string-arrays
	maxStringLength_ += (UInt) round( (Double)maxStringLength_ / 2.0 );
	stringEl.resize(maxStringLength_);

	//store the resulting strings
	for ( UInt i=0; i<actLength-1; i++ ) {
	  stringEl[i] = helpStringEl[i+1];
	}

	// resize help-String-array 
	helpStringEl.resize(maxStringLength_+1);
      }
      else {
	//store the resulting strings
	for ( UInt i=0; i<actLength-1; i++ ) {
	  stringEl[i] = helpStringEl[i+1];
	}
      }
      //store the new input
      stringEl[actLength-1] = newX;
    }

    //compute preisach-sum
    preisachSum_[idx] =  everett(-stringEl[0],stringEl[0]);
    for ( UInt i=0; i<actLength-1; i++ ) {
      preisachSum_[idx] +=  2.0*everett(stringEl[i],stringEl[i+1]);
    }

  }



  Double Preisach :: everett(Double X1, Double X2)
  {
    ENTER_FCN( "Preisach:: everett" );

    Double newY;
    Double diffX = X2 - X1;

    if ( diffX > 0) {
      newY = 0.25*diffX*diffX;
    }
    else {
      newY = -0.25*diffX*diffX;
    }

    Double pPixel =  everettPixel(X1, X2);

    std::cout << "P=" << newY << "  Ppixel=" << pPixel << std::endl;

    //return newY;
    return pPixel;
  }


  Double Preisach :: everettPixel(Double val1, Double val2)
  {
    ENTER_FCN( "Preisach:: everettPixel" );

    Double X1 = std::max(val1,val2);
    Double X2 = std::min(val1,val2);

    UInt M = preisachWeights_.GetSizeRow();
    Double delta = 2.0 / ( (Double) M );
    std::cout << "delta=" << delta << std::endl;
    std::cout << "e1=" << val1 << "  e2=" << val2 << std::endl;
    std::cout << "X1=" << X1 << "  X2=" << X2 << std::endl;

    // we compute the upper left corner within the preisach-domain
    // alpha > X1; beta < X2;

    //compute index for X1 (alpha)
    Integer idx1 = -1;
    Double alpha = -1.0;
    while ( alpha <= X1 ) {
      idx1++;
      alpha += delta;
      std::cout << "idx1=" << idx1 << "  alpha=" << alpha << std::endl;
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

      std::cout << "areaToMuch=" << area << std::endl;

      //reduce the computed area
      Double diffX1 = alpha - X1;
      Double diffX2 = X2    - beta;
      std::cout << " diffX1=" << diffX1 << "  diffX2=" << diffX2 << std::endl;

      Double minusArea;
 
      //check, if we are already on the diagonal!!
      if ( idx1 == idx2 ) {
	minusArea = (   diffX2 * (delta - 0.5*diffX2) 
		      + diffX1 * (delta - 0.5*diffX1)
                      - diffX1*diffX2 
		     ) * preisachWeights_[idx1][idx2];
	std::cout << "minusArea1end=" << minusArea << std::endl;
      }
      else {
	minusArea = ( (diffX1+diffX2 )*delta - diffX1*diffX2 ) 
	  * preisachWeights_[idx1][idx2];

	std::cout << "minusArea1=" << minusArea << std::endl;
	UInt idx = idx1-1;
	while ( idx > idx2 ) {
	  minusArea += diffX2 * delta * preisachWeights_[idx][idx2];
	  idx--;
	  std::cout << "minusArea2=" << minusArea << std::endl;
	}
	minusArea += ( delta*diffX2 - 0.5*diffX2*diffX2 )
	  * preisachWeights_[idx][idx2]; 

	std::cout << "minusArea3=" <<  ( delta*diffX2 - 0.5*diffX2*diffX2 )
	  * preisachWeights_[idx][idx2] << std::endl;
	
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
      std::cout << "val2 < val1" << std::endl;
      area *= -1.0; 
    }

    return area;
  }


  void Preisach :: computePreisachWeights()
  {
    ENTER_FCN( "Preisach::computePreisachWeights" );

    UInt dim = 6;
    preisachWeights_.Resize(dim,dim);
    for ( UInt i=0; i<dim; i++) {
      for ( UInt j=0; j<dim; j++) {
	preisachWeights_[i][j] = 0.5;
      }
    }
  }

  Double Preisach :: normalizeInput(Double Xin)
  {
    ENTER_FCN( "Preisach::normalizeInput" );

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
