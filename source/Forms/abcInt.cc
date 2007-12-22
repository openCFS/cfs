// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "abcInt.hh"


namespace CoupledField {
  
  
  AbsorbingBCsInt::AbsorbingBCsInt(bool isAxi) {

    name_ = "AbsorbingBCsInt";
    isaxi_ = isAxi;
    
  }

  AbsorbingBCsInt::~AbsorbingBCsInt() {
  }


  void AbsorbingBCsInt::SetFactor( const std::string& factor ) {
    mParser_->SetExpr( mHandle_, factor );
  }

  void AbsorbingBCsInt::CalcElementMatrix( Matrix<Double>& elemMat,
                                           EntityIterator& ent1, 
                                           EntityIterator& ent2 ) {

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );

    UInt j = 0;
    Double jacDet, factor, density, compressibility;
    
    Vector<Double> shapeFncAtIp;
    Matrix<Double> partElemMat;
    Double coordAtIp;

    ptelem->SetAnsatzFct( ansatzFct1_ );
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const UInt nrIntPts= ptelem->GetNumIntPoints();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    std::map<RegionIdType, BaseMaterial*> * acouMaterials;    

    // 1) determine correct material list and regionIds
    if ( firstPDEName_ == "acoustic" ) {
      acouMaterials = &firstMaterials_;
    } else {      
      acouMaterials = &secondMaterials_;
    }
    
    it = acouMaterials->find(actElem_->ptVolElem1->regionId);
    if ( it == acouMaterials->end() ) {
      it = acouMaterials->find(actElem_->ptVolElem2->regionId);
    } 

    it->second->GetScalar(density,DENSITY,REAL);
    it->second->GetScalar(compressibility,ACOU_BULK_MODULUS,REAL);
    factor = mParser_->Eval( mHandle_ );
    factor *= density / sqrt( compressibility / density );

    // 2) Calculate a normal mass matrix
    elemMat.Resize(numFncs);
    elemMat.Init();

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {

      jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, 
                                           ent1.GetElem() );
        
      ptelem->GetShFncAtIp(shapeFncAtIp, actIntPt, ent1.GetElem() );
        
      partElemMat.DyadicMult(shapeFncAtIp);
        
      if (isaxi_) {
        // For the entry 1/r things are more complicated
        coordAtIp = 0.0;
        for( j = 0; j < numFncs; j++ ) {
          coordAtIp += ptCoord_[0][j] * shapeFncAtIp[j];
        }
        partElemMat *= 2 * PI * intWeights[actIntPt-1] 
          * factor * jacDet * coordAtIp;
      }
      else 
        partElemMat *= intWeights[actIntPt-1] * factor * jacDet;
        
      elemMat += partElemMat;
    }
  }
  
  
  
  
  
  HeatFluxInt::HeatFluxInt( const std::string& factor, bool isaxi ) {
    
    name_ = "HeatFluxInt";
    isaxi_ = isaxi;
    mParser_->SetExpr( mHandle_, factor );
  }
      
  HeatFluxInt::~HeatFluxInt() {
  }
  
  void HeatFluxInt::CalcElementMatrix( Matrix<Double>& elemMat,
                                       EntityIterator& ent1, 
                                       EntityIterator& ent2 ) {

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );

    Double jacDet, factor;  
    Vector<Double> shapeFncAtIp;
    Matrix<Double> partElemMat;
    Vector<Double> coordAtIp;

    ptelem->SetAnsatzFct( ansatzFct1_ );
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const UInt nrIntPts= ptelem->GetNumIntPoints();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    
    factor = mParser_->Eval( mHandle_ );

    // Calculate a normal mass matrix
    elemMat.Resize(numFncs);
    elemMat.Init();

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {

      jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, 
                                           ent1.GetElem() );
        
      ptelem->GetShFncAtIp(shapeFncAtIp, actIntPt, ent1.GetElem() );
        
      partElemMat.DyadicMult(shapeFncAtIp);
        
      if (isaxi_) {
        coordAtIp = ptCoord_ * shapeFncAtIp;
        partElemMat *= 2 * PI * intWeights[actIntPt-1] * factor * jacDet * coordAtIp[0];
      }
      else {
        partElemMat *= intWeights[actIntPt-1] * factor * jacDet;
      }
      
      elemMat += partElemMat;
    }
  }
  

}
