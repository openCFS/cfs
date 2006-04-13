#include "abcInt.hh"


namespace CoupledField {
  
  
  AbsorbingBCsInt::AbsorbingBCsInt(Boolean isAxi) {
    ENTER_FCN( "AbsorbingBCsInt::AbsorbingBCsInt" );

    isaxi_ = isAxi;
    factor_ = 1.0;
    
  }

  AbsorbingBCsInt::~AbsorbingBCsInt() {
    ENTER_FCN( "AbsorbingBCsInt::~AbsorbingBCsInt" );
  }

  void AbsorbingBCsInt::CalcElementMatrix(Matrix<Double>& ptCoord, 
                                          Matrix<Double> & elemMat) {
    ENTER_FCN( "AbsorbingBCsInt::CalcElementMatrix" );
    
    UInt j = 0;
    Double jacDet, factor, density, compressibility;
    
    Vector<Double> shapeFncAtIp;
    Matrix<Double> partElemMat;
    Double coordAtIp;
    
    const UInt nrIntPts= ptelem->GetNumIntPoints();
    const UInt nrNodes = ptelem->GetNumNodes();
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
    elemMat.Resize(nrNodes);

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {

      jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord);
        
      ptelem->GetShFncAtIp(shapeFncAtIp, actIntPt);
        
      partElemMat.DyadicMult(shapeFncAtIp);
        
      if (isaxi_) {
        // For the entry 1/r things are more complicated
        coordAtIp = 0.0;
        for( j = 0; j < nrNodes; j++ ) {
          coordAtIp += ptCoord[0][j] * shapeFncAtIp[j];
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
