#include <iostream>
#include <fstream>

#include "bubbleStiffInt.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"

namespace CoupledField {



  BubbleStiffInt::BubbleStiffInt(const Double aDensity,
				 Double densityforbubble,
				 Double sonicVel, 
				 Double viscosity,
				 Double bubbleDensity,
				 Double frequency,
				 bool axi )
    : BaseForm(NULL), 
      density_(aDensity)

      
  {
    ENTER_FCN( "BubbleStiffInt::BubbleStiffInt" );
    name_ = "BubbleStiffInt";
    factor_ = 1.0;
    isaxi_ = axi;
    baseType_ = STIFFNESS;
    isSolDependent_ = true;

    densityforbubble_ = densityforbubble;
    sonicVel_ = sonicVel;
    viscosity_ = viscosity;
    bubbleDensity_ = bubbleDensity;
    frequency_ = frequency;

    // Set Expression for parser
    mParser_->SetExpr( mHandle_, "t" );
  }

 
  BubbleStiffInt::~BubbleStiffInt()
  {
    ENTER_FCN( "BubbleStiffInt::~BubbleStiffInt" );
  }

 

  void BubbleStiffInt::CalcElementMatrix( Matrix<Double>& elemMat,
					 EntityIterator& ent1, 
					 EntityIterator& ent2  )  {
    ENTER_FCN( "BubbleStiffInt::CalcElemMatrix" );
  
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );
    
    const UInt nrIntPts= ptelem->GetNumIntPoints();
    const UInt nrNodes = ptelem->GetNumNodes();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Double jacDet;

    Vector<Double> shapeFncAtIp;
    Matrix<Double> partElemMat;
    Vector<Double> CoordAtIP;

    // set matrix to desired size and set all elements to zero
    //    partElemMat.Resize(nrNodes);
    elemMat.Resize(nrNodes);
    elemMat.Init();
    
    // get value of radius und its derivative
    Double radius =  (*radius_)[indexMap_[ent1.GetElem()->elemNum]];
    Double radiusDeriv =  (*radiusDeriv_)[indexMap_[ent1.GetElem()->elemNum]];


    if (mParser_->Eval( mHandle_ ) <= 1.0 / frequency_ ){
      if (indexMap_[ent1.GetElem()->elemNum] == 0)
	std::cerr<< "Faktor 0  in bubbleStiffIt" <<std::endl;
      factor_ = 0.0;
    }
    else{
      if (indexMap_[ent1.GetElem()->elemNum] == 0)
	std::cerr <<"Faktor computed  in bubbleStiffIt" <<std::endl;
      factor_ = 4.0 * PI * bubbleDensity_ * radius * radius;
      factor_*= (1.0 + radiusDeriv /sonicVel_) ;  
      factor_ /= ((1.0- (radiusDeriv /sonicVel_)) * radius 
		  + 4.0 * viscosity_ / densityforbubble_ /sonicVel_);     
    }


    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {

      jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_);
        
      ptelem-> GetShFncAtIp(shapeFncAtIp, actIntPt);
        
      partElemMat.DyadicMult(shapeFncAtIp);
        

     

      if (isaxi_) {
        CoordAtIP = ptCoord_ * shapeFncAtIp;
        partElemMat *= 2 * PI * intWeights[actIntPt-1] * density_ * factor_* jacDet * CoordAtIP[0];
      }
      else 
        partElemMat *= intWeights[actIntPt-1] * density_ * factor_ * jacDet;
        
      elemMat += partElemMat;
    }
  
  }


  void  BubbleStiffInt::SetValues( StdVector<UInt> elemNumbers, Vector<Double> * radius,
				   Vector<Double> * radiusDeriv ) {
    ENTER_FCN( "BubbleStiffInt::SetValues" );
    
    elemNumbers_ = elemNumbers;
    radius_ = radius;
    std::cerr<<"SetValues bubbleStiff size radius"<< (*radius).GetSize()<<std::endl;
    radiusDeriv_ = radiusDeriv;


    // map index positions
    for( UInt i = 0; i < elemNumbers.GetSize(); i++ ) {
      indexMap_[elemNumbers[i]] = i;
    }

  }

  
} // end namespace CoupledField



