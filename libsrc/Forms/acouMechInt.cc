#include "acouMechInt.hh"


namespace CoupledField {
  
  
  AcouMechInt::AcouMechInt( UInt dofsPerNode, Boolean isAxi) {
    ENTER_FCN( "AcouMechInt::AcouMechInt" );

    isaxi_ = isAxi;
    dofs_ = dofsPerNode;
    
  }

  AcouMechInt::~AcouMechInt() {
    ENTER_FCN( "AcouMechInt::~AcouMechInt" );
  }

  void AcouMechInt::CalcElementMatrix(Matrix<Double>& ptCoord, 
                                           Matrix<Double> & elemMat) {
    ENTER_FCN( "AcouMechInt::CalcElementMatrix" );
    
    UInt j = 0;
    Double jacDet, density;
    
    Vector<Double> shapeFncAtIp;
    Matrix<Double> partHelpMat, helpMat;
    Double coordAtIp;
    
    const UInt nrIntPts= ptelem->GetNumIntPoints();
    const UInt nrNodes = ptelem->GetNumNodes();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    
    Integer index = -1;
    
    const StdVector<RegionIdType> * acouRegionIds;
    const MaterialData * acouMaterials;
    
    // determine correct material list and regionIds
    if ( firstPDEName_ == "acoustic" ) {
      acouRegionIds = &firstRegionIds_; 
      acouMaterials = firstMaterials_;
    } else {
      acouRegionIds = &secondRegionIds_; 
      acouMaterials = secondMaterials_;
    }
    
    index = acouRegionIds->Find(actElem_->ptVolElem1->regionId);
    if ( index == -1 ) {
      index = acouRegionIds->Find(actElem_->ptVolElem2->regionId);
    } else {
      normal_ *= -1.0;
    }

    density = acouMaterials[index].GetDensity();
    if (formulation_ == ACOU_PRESSURE && firstPDEName_ == "acoustic" ) {
      //multiplicative factor in case of pressure formulation
      density *= -density;
    }
    else if (formulation_ == ACOU_PRESSURE && firstPDEName_ == "mechanic" ) {
      density = 1.0;
    }
    
    // 2) Calculate a normal mass matrix
    helpMat.Resize(nrNodes);
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {

      jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord);
        
      ptelem->GetShFncAtIp(shapeFncAtIp, actIntPt);
        
      partHelpMat.DyadicMult(shapeFncAtIp);
        
      if (isaxi_) {
        // For the entry 1/r things are more complicated
        coordAtIp = 0.0;
        for( j = 0; j < nrNodes; j++ ) {
          coordAtIp += ptCoord[0][j] * shapeFncAtIp[j];
        }
        partHelpMat *= 2 * PI * intWeights[actIntPt-1] 
          * density * jacDet * coordAtIp;
      }
      else 
        partHelpMat *= intWeights[actIntPt-1] * density * jacDet;
        
      helpMat += partHelpMat;
    }


    // 3) Create a multi-dof matrix, multiplied by normal vector
    elemMat.Resize( nrNodes*dofs_, nrNodes );

    //for pressure forumation
    if (formulation_ == ACOU_PRESSURE && firstPDEName_ == "acoustic" ) {
      for ( UInt iRow = 0; iRow < nrNodes; iRow++ ) {
	for ( UInt iCol = 0; iCol < nrNodes; iCol++ ) {
	  for ( UInt iDof = 0; iDof < dofs_; iDof++ ) {
	    elemMat[iRow][iCol*dofs_+iDof] = 
	      normal_[iDof] * helpMat[iRow][iCol];
	  }
	}
      }
    }
    else {
      for ( UInt iRow = 0; iRow < nrNodes; iRow++ ) {
	for ( UInt iCol = 0; iCol < nrNodes; iCol++ ) {
	  for ( UInt iDof = 0; iDof < dofs_; iDof++ ) {
	    elemMat[iRow*dofs_+iDof][iCol] = 
	      normal_[iDof] * helpMat[iRow][iCol];
	  }
	}
      }
    }
  }
  

}
