#include "abcInt.hh"


namespace CoupledField {
  
  
  AbsorbingBCsInt::AbsorbingBCsInt(bool isAxi) {
    ENTER_FCN( "AbsorbingBCsInt::AbsorbingBCsInt" );

    name_ = "AbsorbingBCsInt";
    isaxi_ = isAxi;
    factor_ = 1.0;
    
  }

  AbsorbingBCsInt::~AbsorbingBCsInt() {
    ENTER_FCN( "AbsorbingBCsInt::~AbsorbingBCsInt" );
  }

  void AbsorbingBCsInt::CalcElementMatrix( Matrix<Double>& elemMat,
                                           EntityIterator& ent1, 
                                           EntityIterator& ent2 ) {
    ENTER_FCN( "AbsorbingBCsInt::CalcElementMatrix" );

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

    // determine correct material list and regionIds
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
    factor = factor_ * density / sqrt( compressibility / density );

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
  

}
