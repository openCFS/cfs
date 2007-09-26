// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include "bubbleDampInt.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"

namespace CoupledField {



  BubbleDampInt::BubbleDampInt(const Double aDensity,
			       Double densityforbubble,
			       Double sonicVel, 
			       Double viscosity,
			       Double bubbleDensity,
			       Double frequency,
			       bool axi )
    : BaseForm(NULL), 
      density_(aDensity)     
  {
    name_ = "BubbleDampInt";
    isaxi_ = axi;
    baseType_ = DAMPING;
    isSolDependent_ = true;

    densityforbubble_ = densityforbubble;
    sonicVel_ = sonicVel;
    viscosity_ = viscosity;
    bubbleDensity_ = bubbleDensity;
    frequency_ = frequency;

    // Set Expression for parser
    mParser_->SetExpr( mHandle_, "t" );

  }

 
  BubbleDampInt::~BubbleDampInt()
  {
  }

 

  void BubbleDampInt::CalcElementMatrix( Matrix<Double>& elemMat,
					 EntityIterator& ent1, 
					 EntityIterator& ent2  )  {
  

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );
    
    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt nrIntPts= ptelem->GetNumIntPoints();
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Double jacDet;

    Vector<Double> shapeFncAtIp;
    Matrix<Double> partElemMat;
    Vector<Double> CoordAtIP;

    // set matrix to desired size and set all elements to zero
    //    partElemMat.Resize(nrNodes);
    elemMat.Resize( numFncs );
    elemMat.Init();
    
    // get value of radius und its derivative
    Double radius =  (*radius_)[indexMap_[ent1.GetElem()->elemNum]];
    Double radiusDeriv =  (*radiusDeriv_)[indexMap_[ent1.GetElem()->elemNum]];
    Double factor;


    Double dampFact = 0.0;
    if (mParser_->Eval( mHandle_ ) < 10.0 / frequency_ ){
//       if (indexMap_[ent1.GetElem()->elemNum] == 0)
//        	std::cerr<<"Faktor 0  in bubbleDampIt" <<std::endl;
      dampFact = 4.0 * PI * bubbleDensity_ * radius * radius * radius / sonicVel_ /((1.0 - (radiusDeriv /sonicVel_)) * radius + 4.0 * viscosity_ / densityforbubble_ /sonicVel_); 
      //     std::cout<<dampFact<<std::endl;

//       if( ent1.GetElem()->elemNum == 294 ){
//      	std::cout<<mParser_->Eval( mHandle_ ) <<"   "<< dampFact<<std::endl; 
//       }

    }
    else{   
//       if (indexMap_[ent1.GetElem()->elemNum] == 0)
//        	std::cerr<<"Faktor computed  in bubbleDampIt" <<std::endl;
      factor= 4.0 * PI * bubbleDensity_ * radius * radius * radius / sonicVel_ /((1.0 - (radiusDeriv /sonicVel_)) * radius + 4.0 * viscosity_ / densityforbubble_ /sonicVel_); 
    }

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {

      jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, ent1.GetElem() );
        
      ptelem->GetShFncAtIp(shapeFncAtIp, actIntPt, ent1.GetElem() );
        
      partElemMat.DyadicMult(shapeFncAtIp);
     

      if (isaxi_) {
        CoordAtIP = ptCoord_ * shapeFncAtIp;
        partElemMat *= 2 * PI * intWeights[actIntPt-1] * density_ * factor* jacDet * CoordAtIP[0];
      }
      else 
        partElemMat *= intWeights[actIntPt-1] * density_ * factor * jacDet;
        
      elemMat += partElemMat;
    }
  
  }


  void  BubbleDampInt::SetValues( StdVector<UInt> elemNumbers, Vector<Double> * radius,
				   Vector<Double> * radiusDeriv ) {
    
    elemNumbers_ = elemNumbers;
    radius_ = radius;
    radiusDeriv_ = radiusDeriv;


    // map index positions
    for( UInt i = 0; i < elemNumbers.GetSize(); i++ ) {
      indexMap_[elemNumbers[i]] = i;
    }

  }

  
} // end namespace CoupledField



