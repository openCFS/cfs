// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "ansatzFct.hh"
#include "General/environment.hh"

namespace CoupledField {


  AnsatzFct::AnsatzFct() {
    ENTER_FCN( "AnsatzFct::AnsatzFct" );
    
    isIsotropic_ = true;

  }

  AnsatzFct::~AnsatzFct() {
    ENTER_FCN( "AnsatzFct::~AnsatzFct" );
  }

  ConstFct::ConstFct() {
    ENTER_FCN( "ConstFct::ConstFct" );
    type_ = CONST;
  }

  LagrangeFct::LagrangeFct() {
    ENTER_FCN( "LagrangeFct::LagrangeFct" );
    type_ = LAGRANGE;
  }

  LegendreFct::LegendreFct() {
    ENTER_FCN( "LegendreFct::LegendreFct" );
    type_ = LEGENDRE;
    isoOrder_ = 0;
    maxOrder_ = 0;
    isIsotropic_ = true;
    subSpace_ = PRODUCT;
  }
  
  void LegendreFct::SetIsoOrder( UInt order ) {
    ENTER_FCN( "LegendreFct::SetIsoOrder" );

    isoOrder_ = order;
    isIsotropic_ = true;
    maxOrder_ = order;
  }
  
  
  UInt LegendreFct::GetIsoOrder() const {
    ENTER_FCN( " LegendreFct::GetIsoOrder" );
    if( !isIsotropic_) {
      Error( "Approximation is anisotropic!", __FILE__, __LINE__ );
    }
      
    return isoOrder_;
  }
  
  
  void LegendreFct::SetAnisoOrder( Matrix<UInt>& order ) {
    ENTER_FCN( "LegendreFct::SetAnsioOrder" );

    anOrder_ = order;
    isIsotropic_ = false;
    maxOrder_ = 0;

    // determine maximum w.r.t. to local coordinate directions
    maxOrderLocDir_.Resize( anOrder_.GetSizeRow() );
    maxOrderDof_.Resize( anOrder_.GetSizeCol() );
    maxOrderLocDir_.Init( 0 );
    maxOrderDof_.Init( 0 );
    for( UInt iLoc = 0; iLoc < anOrder_.GetSizeRow(); iLoc++ ) {
      for( UInt iDof = 0; iDof < anOrder_.GetSizeCol(); iDof++ ) {
        if( anOrder_[iLoc][iDof] > maxOrderLocDir_[iLoc] ) {
          maxOrderLocDir_[iLoc] = anOrder_[iLoc][iDof];
        }
        if( anOrder_[iLoc][iDof] > maxOrderDof_[iDof] ) {
          maxOrderDof_[iDof] = anOrder_[iLoc][iDof];
        }
        if( anOrder_[iLoc][iDof] > maxOrder_ ) {
          maxOrder_ = anOrder_[iLoc][iDof];
        }
      }
    }
    
    std::cerr << "order = \n" << order << std::endl;
    std::cerr << "maxOrderDof: " << maxOrderDof_.Serialize() << std::endl;
    std::cerr << "maxOrderLocDir: " << maxOrderLocDir_.Serialize() << std::endl;
  }
  
  
  const Matrix<UInt>& LegendreFct::GetAnisoOrder() const {
    ENTER_FCN( "LegendreFct::GetAnisoOrder" );
    if( isIsotropic_) {
      Error( "Approximation is isotropic!", __FILE__, __LINE__ );
    }
    return anOrder_;
  }
  
  
  UInt LegendreFct::GetMaxOrderLocDir( UInt locDir ) const {
    ENTER_FCN( "LegendreFct::GetMaxOrderLocDir" );
    if( isIsotropic_ ) {
      return isoOrder_;
    } else {
      return maxOrderLocDir_[locDir];
    }
  }
  
  
  UInt LegendreFct::GetMaxOrderDof( UInt dof ) const {
    ENTER_FCN( "LegendreFct::GetMaxOrderDof" );
    if( isIsotropic_ ) {
      return isoOrder_;
    } else {
      return maxOrderDof_[dof];
    }
  }

  UInt LegendreFct::GetMaxOrder( ) const {
    return maxOrder_;
  }
  
}
