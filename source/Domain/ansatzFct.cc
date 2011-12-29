// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "General/exception.hh"
#include "ansatzFct.hh"

namespace CoupledField {


  AnsatzFct::AnsatzFct() {

    isIsotropic_ = true;
    isDiscontinuous_ = false;

  }

  AnsatzFct::~AnsatzFct() {
  }

  bool AnsatzFct::IsDiscontinuous() const { 
    return isDiscontinuous_; 
  }

  void AnsatzFct::SetDiscontinuity(bool discontinuous){
    isDiscontinuous_ = discontinuous;
  }

  ConstFct::ConstFct() {
    type_ = CONST;
  }

  LagrangeFct::LagrangeFct() {
    type_ = LAGRANGE;
  }

//======================================
//Spectral Ansatz function
//======================================
// Define supporting points for the polinomial
    Double SpectralFct::l1[][1]  = {
      {-1},
      {1}
    };

    Double SpectralFct::l2[][1] = {
      {-1.000000000000000e+00},
      { 0.000000000000000e+00},
      { 1.000000000000000e+00}
    };

    Double SpectralFct::l3[][1] = {
      {-1.000000000000000e+00},
      {-4.472135954999579e-01},
      { 4.472135954999579e-01},
      { 1.000000000000000e+00}
    };

    Double SpectralFct::l4[][1] = {
      {-1.000000000000000e+00},
      {-6.546536707079771e-01},
      { 0.000000000000000e+00},
      { 6.546536707079771e-01},
      { 1.000000000000000e+00}
    };

    Double SpectralFct::l5[][1] = {
      {-1.000000000000000e+00},
      {-7.650553239294647e-01},
      {-2.852315164806451e-01},
      { 2.852315164806451e-01},
      { 7.650553239294647e-01},
      { 1.000000000000000e+00}
    };

    Double SpectralFct::l6[][1] = {
      {-1.000000000000000e+00},
      {-8.302238962785670e-01},
      {-4.688487934707142e-01},
      { 0.0000000000000000000},
      { 4.688487934707142e-01},
      { 8.302238962785670e-01},
      { 1.000000000000000e+00}
    };
    Double SpectralFct::l7[][1] = {
      {-1.000000000000000e+00},
      {-8.717401485096066e-01},
      {-5.917001814331423e-01},
      {-2.092992179024789e-01},
      { 2.092992179024789e-01},
      { 5.917001814331423e-01},
      { 8.717401485096066e-01},
      { 1.000000000000000e+00}
    };

  SpectralFct::SpectralFct() {
    type_ = SPECTRAL;
    Order_ = 1;
    //unlike for hirarchical functions...
    isIsotropic_ = false;
  }

  void SpectralFct::SetOrder( UInt order ){
    Order_ = order;
    supportingPoints_.Resize(Order_+1);
    switch(Order_){
      case 1:
        for ( UInt i = 0;i<=Order_ ;i++ )
          supportingPoints_[i] = l1[i][0];
        break;
      case 2:
        for ( UInt i = 0;i<=Order_ ;i++ )
          supportingPoints_[i] = l2[i][0];
        break;
      case 3:
        for ( UInt i = 0;i<=Order_ ;i++ )
          supportingPoints_[i] = l3[i][0];
        break;
      case 4:
         for ( UInt i = 0;i<=Order_ ;i++ )
           supportingPoints_[i] = l4[i][0];
        break;
      case 5:
        for ( UInt i = 0;i<=Order_ ;i++ )
          supportingPoints_[i] = l5[i][0];
        break;
      case 6:
        for ( UInt i = 0;i<=Order_ ;i++ )
          supportingPoints_[i] = l6[i][0];
        break;
      case 7:
        for ( UInt i = 0;i<=Order_ ;i++ )
          supportingPoints_[i] = l7[i][0];
        break;
      default:
        EXCEPTION( "Supplied Order of approximation not supported" );
        break;
    }


  }

  void SpectralFct::EvaluatePolynomial( Vector<Double> & shape, Double coord )
  {
    shape.Resize(Order_+1);
    shape.Init();
    //get iutegration Pointes
    // From Zienkiewicz, The Finite Element Method. Vol 1, page 122.
    // corner nodes
     for( UInt i=0; i<Order_+1; i++)
     {
      shape[i] = 1.0;
      for( UInt p = 0;p< Order_ + 1; p++)
      {
        if(p==i)
          continue;
        else
          shape[i] *= (coord - supportingPoints_[p]) / (supportingPoints_[i] - supportingPoints_[p]);
      }
    }
  }

  void SpectralFct::EvaluateDerivPolynomial( Vector<Double> & deriv, Double coord )
  {
    deriv.Resize(Order_+1);
    deriv.Init();
    for ( UInt i = 0; i<=Order_  ; i++)
    {
      Double sum = 1.0;
      for ( UInt k=0;k<=Order_ ; k++ )
      {
        if ( k != i )
        {
          sum=1.0;
          for ( UInt l=0;l<=Order_ ; l++)
            if ( (l!=i) && (l != k))
              sum *= (coord - supportingPoints_[l]);
          deriv[i] += sum;
        }
      }
    }
    for( UInt i=0; i< Order_+1; i++)
    {
      for( UInt p = 0;p< Order_+1; p++)
      {
        if(p==i)
          continue;
        else
          deriv[i] *= 1.0 / (supportingPoints_[i] - supportingPoints_[p]);
      }
    }
  }

  UInt SpectralFct::GetOrder( ) const{
    return Order_;
  }

  NedelecFct::NedelecFct() {
    type_ = NEDELEC;
  }

  LegendreFct::LegendreFct() {
    type_ = LEGENDRE;
    isoOrder_ = 0;
    maxOrder_ = 0;
    isIsotropic_ = true;
    subSpace_ = PRODUCT;
  }

  void LegendreFct::SetIsoOrder( UInt order ) {

    isoOrder_ = order;
    isIsotropic_ = true;
    maxOrder_ = order;
  }


  UInt LegendreFct::GetIsoOrder() const {
    if( !isIsotropic_) {
      EXCEPTION( "Approximation is anisotropic!" );
    }

    return isoOrder_;
  }


  void LegendreFct::SetAnisoOrder( Matrix<UInt>& order ) {

    anOrder_ = order;
    isIsotropic_ = false;
    maxOrder_ = 0;

    // determine maximum w.r.t. to local coordinate directions
    maxOrderLocDir_.Resize( anOrder_.GetNumRows() );
    maxOrderDof_.Resize( anOrder_.GetNumCols() );
    maxOrderLocDir_.Init( 0 );
    maxOrderDof_.Init( 0 );
    for( UInt iLoc = 0; iLoc < anOrder_.GetNumRows(); iLoc++ ) {
      for( UInt iDof = 0; iDof < anOrder_.GetNumCols(); iDof++ ) {
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

//    std::cerr << "order = \n" << order << std::endl;
//    std::cerr << "maxOrderDof: " << maxOrderDof_.Serialize() << std::endl;
//    std::cerr << "maxOrderLocDir: " << maxOrderLocDir_.Serialize() << std::endl;
  }


  const Matrix<UInt>& LegendreFct::GetAnisoOrder() const {
    if( isIsotropic_) {
      EXCEPTION( "Approximation is isotropic!" );
    }
    return anOrder_;
  }


  UInt LegendreFct::GetMaxOrderLocDir( UInt locDir ) const {
    if( isIsotropic_ ) {
      return isoOrder_;
    } else {
      return maxOrderLocDir_[locDir];
    }
  }


  UInt LegendreFct::GetMaxOrderDof( UInt dof ) const {
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
