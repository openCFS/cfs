// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <stddef.h>
#include <map>
#include <ostream>
#include <string>
#include <utility>

#include "Domain/elem.hh"
#include "Domain/entityList.hh"
#include "Domain/surfElem.hh"
#include "Elements/basefe.hh"
#include "General/exception.hh"
#include "MatVec/exprt/xpr1.hh"
#include "MatVec/exprt/xpr2.hh"
#include "MatVec/matrix.hh"
#include "MatVec/vector.hh"
#include "Materials/baseMaterial.hh"
#include "acouMechInt.hh"


namespace CoupledField {
  
  
  AcouMechInt::AcouMechInt( UInt dofsPerNode, bool isAxi) {

    name_ = "AcouMechInt";
    isaxi_ = isAxi;
    dofs_ = dofsPerNode;
    formulation_ = NO_SOLUTION_TYPE; 
    
    // this bilinearform is never symmetric
    isSymmetric_ = false;
    
  }

  AcouMechInt::~AcouMechInt() {
  }
  
  void AcouMechInt::SetFormulation( SolutionType aformulation) {
    
    formulation_ = aformulation;
  }
    
  void AcouMechInt::CalcElementMatrix( Matrix<Double>& elemMat,
                                       EntityIterator& ent1, 
                                       EntityIterator& ent2 ) {
    
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );

    UInt j = 0;
    Double jacDet, density;
    
    Vector<Double> shapeFncAtIp;
    Matrix<Double> partHelpMat, helpMat;
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
    
    if ( actElem_->ptVolElem1 == NULL ) {
      EXCEPTION("Cannot apply mechanic-acoustic coupling on surface element "
          << actElem_->elemNum
          << ", because it has no adjacent volume element.\n"
          << "Go check your mesh file!");
    }
    it = acouMaterials->find(actElem_->ptVolElem1->regionId);
    if ( it == acouMaterials->end() ) {
      if ( actElem_->ptVolElem2 == NULL ) {
        EXCEPTION("Cannot apply mechanic-acoustic coupling on surface element "
            << actElem_->elemNum
            << ", because it has no adjacent volume element in an acoustic region.\n"
            << "Go check your mesh file!");
      }
      it = acouMaterials->find(actElem_->ptVolElem2->regionId);
      if ( it == acouMaterials->end()) {
        EXCEPTION("Acoustic parent region of surface element "
            << actElem_->elemNum << " could not be determined.");

      }
    } else {
      normal_ *= -1.0;
    }

    it->second->GetScalar(density,DENSITY,Global::REAL);
    if (formulation_ == ACOU_PRESSURE && firstPDEName_ == "acoustic" ) {
      //multiplicative factor in case of pressure formulation
      density *= -density;
    }
    else if (formulation_ == ACOU_PRESSURE && firstPDEName_ == "mechanic" ) {
      density = 1.0;
    }
    
    // 2) Calculate a normal mass matrix
    helpMat.Resize(numFncs);
    helpMat.Init();
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {

      jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_,
                                           ent1.GetElem() );
        
      ptelem->GetShFncAtIp(shapeFncAtIp, actIntPt, ent1.GetElem() );
        
      partHelpMat.DyadicMult(shapeFncAtIp);
        
      if (isaxi_) {
        // For the entry 1/r things are more complicated
        coordAtIp = 0.0;
        for( j = 0; j < numFncs; j++ ) {
          coordAtIp += ptCoord_[0][j] * shapeFncAtIp[j];
        }
        partHelpMat *= 2 * PI * intWeights[actIntPt-1] 
          * density * jacDet * coordAtIp;
      }
      else 
        partHelpMat *= intWeights[actIntPt-1] * density * jacDet;
        
      helpMat += partHelpMat;
    }


    // 3) Create a multi-dof matrix, multiplied by normal vector


    //for pressure forumation
    if (formulation_ == ACOU_PRESSURE && firstPDEName_ == "acoustic" ) {
      elemMat.Resize( numFncs, numFncs*dofs_ );
      for ( UInt iRow = 0; iRow < numFncs; iRow++ ) {
        for ( UInt iCol = 0; iCol < numFncs; iCol++ ) {
          for ( UInt iDof = 0; iDof < dofs_; iDof++ ) {
            elemMat[iRow][iCol*dofs_+iDof] = 
                normal_[iDof] * helpMat[iCol][iRow];
          }
        }
      }
    }
    else {
      elemMat.Resize( numFncs*dofs_, numFncs );
      for ( UInt iRow = 0; iRow < numFncs; iRow++ ) {
        for ( UInt iCol = 0; iCol < numFncs; iCol++ ) {
          for ( UInt iDof = 0; iDof < dofs_; iDof++ ) {
            elemMat[iRow*dofs_+iDof][iCol] = 
                normal_[iDof] * helpMat[iRow][iCol];
          }
        }
      }
    }
//     std::cerr<<"Density :"<<density<<std::endl;
//     std::cerr<<"Matrix\n"<<elemMat<<std::endl;
  }
  

}
