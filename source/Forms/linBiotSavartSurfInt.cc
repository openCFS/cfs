// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "linBiotSavartSurfInt.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "Utils/biotSavart.hh"
#include "Domain/domain.hh"

DECLARE_LOG(forms)

namespace CoupledField {

  // =================================================================
  // LinBiotSavartSurfInt 
  // =================================================================


  LinBiotSavartSurfInt::LinBiotSavartSurfInt( shared_ptr<BiotSavart> bisa, 
                                              bool isAxi )
    : LinearSurfForm() {

    name_ = "LinBiotSavartSurfInt";
    isaxi_ = isAxi;

    biotSavart_ = bisa;
    mParser_->SetExpr( mHandle_, "1.0" );      
  }

  LinBiotSavartSurfInt::~LinBiotSavartSurfInt() {
  }

  void LinBiotSavartSurfInt::SetFactor( const std::string& value ) {
    mParser_->SetExpr( mHandle_, value ); 
  }
  void LinBiotSavartSurfInt::CalcElemVector( Vector<Double> & elemVec,
                                             EntityIterator& ent ) {
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );
    ptelem->SetAnsatzFct( ansatzFct1_ );

    const UInt nrIntPts = ptelem->GetNumIntPoints();
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const Vector<Double> * intPoints = ptelem->GetIntPoints();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    
    // Calculate element vector  
    elemVec.Resize(numFncs);
    elemVec.Init(0.0);

    Double value = 0.0;    
    Double hNormal = 0.0;
    const Double mu0=1.25664e-6;

    Vector<Double> hField,  shapeFnc, CoordAtIP, normalOutwards;
    
    /* Determine correct normal:
     * By definition, the normal on the boundary points out of the
     * simulation domain. 
     * To guarantee this, we look for the correct volume neighbour of the 
     * surface element, i.e. the one, which belongs to the magnetic domain.
     * The nomalSign_ of the surface element is defined to point into the
     * first volume element neighbour.
     */
    
    // first check, if the current element has a volume element at all
    if ( actElem_->ptVolElem1 == NULL ) {
      EXCEPTION("Cannot apply absorbing BC on surface element "
          << actElem_->elemNum
          << ", because it has no adjacent volume element.\n"
          << "Go check your mesh file!");
    }

    std::map<RegionIdType, BaseMaterial*>::iterator it;
    it = materials_.find(actElem_->ptVolElem1->regionId);
    if ( it == materials_.end() ) {
      if ( actElem_->ptVolElem2 == NULL ) {
        EXCEPTION("Cannot apply absorbing BC on surface element "
            << actElem_->elemNum
            << ", because it has no adjacent volume element in an acoustic region.\n"
            << "Go check your mesh file!");
      }
      it = materials_.find(actElem_->ptVolElem2->regionId);
      if ( it == materials_.end()) {
        EXCEPTION("Acoustic parent region of surface element "
            << actElem_->elemNum << " could not be determined.");
      } else {
        normalOutwards = normal_ * -1.0;  
      }
    } else {
      normalOutwards = normal_;
    }

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {
      ptelem->GetShFncAtIp(shapeFnc,actIntPt,ent.GetElem() );
      value = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, ent.GetElem()) 
              * intWeights[actIntPt-1] * mParser_->Eval(mHandle_);
      
      // calculate Biot-Savart field
      ptelem->Local2GlobalCoord( CoordAtIP, intPoints[actIntPt-1],ptCoord_, ent.GetElem() );
      biotSavart_->CalcFieldAtPoint( hField, CoordAtIP );
      
      hNormal =  normalOutwards * hField;
      elemVec += shapeFnc * ( hNormal * mu0 * value );
      
    } // loop over integration points

  }
  
} // end of namespace
