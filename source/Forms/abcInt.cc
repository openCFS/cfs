// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "abcInt.hh"


namespace CoupledField {


  AbsorbingBCsInt::AbsorbingBCsInt( bool optimalImpedance,
                                    const std::string &imp_magn,
                                    const std::string &imp_phase,
                                    bool isAxi )
  {

    name_ = "AbsorbingBCsInt";
    optImp_ = optimalImpedance;
    impedMagn_ = imp_magn;
    impedPhase_ = imp_phase;
    isaxi_ = isAxi;
    
    // determine if impedance is always real or must be considered complex
    if ( impedPhase_.size() == 0 )
    {
      isComplex_ = false;
    }
    else
    {
      mParser_->SetExpr( mHandle_, impedPhase_ );
      if ( mParser_->IsExprConstant( mHandle_ ) )
      {
        isComplex_ = ( mParser_->Eval(mHandle_) != 0.0 );
      }
      else
      {
        isComplex_ = true;
      }
    }

  }

  AbsorbingBCsInt::~AbsorbingBCsInt() {
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
    const UInt nrIntPts = ptelem->GetNumIntPoints();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();

    std::map<RegionIdType, BaseMaterial*>::iterator it;
    std::map<RegionIdType, BaseMaterial*> * acouMaterials;

    // 1) determine correct material list and regionIds
    if ( firstPDEName_ == "acoustic" ) {
      acouMaterials = &firstMaterials_;
    } else {
      acouMaterials = &secondMaterials_;
    }

    if ( actElem_->ptVolElem1 == NULL ) {
      EXCEPTION("Cannot apply absorbing BC on surface element "
          << actElem_->elemNum
          << ", because it has no adjacent volume element.\n"
          << "Go check your mesh file!");
    }
    it = acouMaterials->find(actElem_->ptVolElem1->regionId);
    if ( it == acouMaterials->end() ) {
      if ( actElem_->ptVolElem2 == NULL ) {
        EXCEPTION("Cannot apply absorbing BC on surface element "
            << actElem_->elemNum
            << ", because it has no adjacent volume element in an acoustic region.\n"
            << "Go check your mesh file!");
      }
      it = acouMaterials->find(actElem_->ptVolElem2->regionId);
      if ( it == acouMaterials->end()) {
        EXCEPTION("Acoustic parent region of surface element "
            << actElem_->elemNum << " could not be determined.");
      }
    }

    it->second->GetScalar( density, DENSITY, Global::REAL );
    if ( ! impedMagn_.size() == 0 )
    {
      mParser_->SetExpr( mHandle_, impedMagn_ );
      factor = mParser_->Eval( mHandle_ );
    }
    else
    {
      factor = 1.0;
    }
    if ( optImp_ ) // compute optimal impedance rho/c ?
    {
      it->second->GetScalar(compressibility, ACOU_BULK_MODULUS, Global::REAL);
      factor *= density / sqrt( compressibility / density );
    }
    else
    {
      factor = density * density / factor;
    }

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

  void AbsorbingBCsInt::CalcElementMatrix( Matrix<Complex>& elemMat,
                                           EntityIterator& ent1,
                                           EntityIterator& ent2 ) {
    if ( optImp_ )
    {
      EXCEPTION("An optimally absorbing impedance can never be complex!");
    }
    
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );

    UInt j = 0;
    Double jacDet, density, magnitude, phase, realPart, imagPart;

    Vector<Double> shapeFncAtIp;
    Matrix<Double> partElemMat;
    Double coordAtIp;

    ptelem->SetAnsatzFct( ansatzFct1_ );
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const UInt nrIntPts = ptelem->GetNumIntPoints();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();

    std::map<RegionIdType, BaseMaterial*>::iterator it;
    std::map<RegionIdType, BaseMaterial*> * acouMaterials;

    // 1) determine correct material list and regionIds
    if ( firstPDEName_ == "acoustic" ) {
      acouMaterials = &firstMaterials_;
    } else {
      acouMaterials = &secondMaterials_;
    }

    if ( actElem_->ptVolElem1 == NULL ) {
      EXCEPTION("Cannot apply absorbing BC on surface element "
          << actElem_->elemNum
          << ", because it has no adjacent volume element.\n"
          << "Go check your mesh file!");
    }
    it = acouMaterials->find(actElem_->ptVolElem1->regionId);
    if ( it == acouMaterials->end() ) {
      if ( actElem_->ptVolElem2 == NULL ) {
        EXCEPTION("Cannot apply absorbing BC on surface element "
            << actElem_->elemNum
            << ", because it has no adjacent volume element in an acoustic region.\n"
            << "Go check your mesh file!");
      }
      it = acouMaterials->find(actElem_->ptVolElem2->regionId);
      if ( it == acouMaterials->end()) {
        EXCEPTION("Acoustic parent region of surface element "
            << actElem_->elemNum << " could not be determined.");
      }
    }

    it->second->GetScalar( density, DENSITY, Global::REAL );
    
    mParser_->SetExpr( mHandle_, impedMagn_ );
    magnitude = mParser_->Eval( mHandle_ );
    mParser_->SetExpr( mHandle_, impedPhase_ );
    phase = mParser_->Eval( mHandle_ );
    
    realPart = magnitude * cos( phase / 180 * PI );
    imagPart = magnitude * sin( phase / 180 * PI );
    
    Complex impedance( realPart, imagPart );
    Complex factor = density * density / impedance;

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
        partElemMat *= 2 * PI * intWeights[actIntPt-1] * jacDet * coordAtIp;
      }
      else
        partElemMat *= intWeights[actIntPt-1] * jacDet;

      elemMat += partElemMat * factor;
    }
  }


//=============================== HeatFluxInt ================================
  

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
