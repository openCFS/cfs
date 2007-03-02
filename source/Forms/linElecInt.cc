#include <iostream>
#include <fstream>

#include "General/environment.hh"
#include "linElecInt.hh"

namespace CoupledField {


  // ============
  //   calcBMat
  // ============
  void linElecInt::calcBMat( Matrix<Double> &bMat, UInt ip,
			     Matrix<Double> &ptCoord ) {

    ENTER_FCN( "linElecInt::calcBMat" );

    // Obtain info on number of elements' funtions
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    
    // Set correct size of matrix B and initialise with zeros
    bMat.Resize( dim_, numFncs );
    bMat.Init();

    // Get derivatives of local shape functions with respect to global
    // coords (format: nrNodes x spaceDim)
    Matrix<Double> xiDx;
    ptelem->SetAnsatzFct( ansatzFct1_ );
    ptelem->GetGlobDerivShFncAtIp( xiDx, ip, ptCoord, it1_.GetElem() );

    if ( subTensorType_ == FULL ) {
      // The matrix bMat can be seen as a 3 x numFncs block-vector.
      // The k-th entry of this block vector corresponds to the matrix
      // B of the BDB product evaluated at the k-th node of the finite
      // element. 
      for( UInt actNode = 0; actNode < numFncs; actNode++ ) {
	bMat[0][actNode] = xiDx[actNode][0];
	bMat[1][actNode] = xiDx[actNode][1];
	bMat[2][actNode] = xiDx[actNode][2];
      }
    }
    else  {
      // The matrix bMat can be seen as a 1 x numFncs block-vector.
      // The k-th entry of this block vector corresponds to the matrix
      // B of the ADB product evaluated at the k-th node of the finite
      // element. We assume that the first coordinate equals y and the
      // second z.
      for( UInt actNode = 0; actNode < numFncs; actNode++ ) {
        bMat[0][actNode] = xiDx[actNode][0];
        bMat[1][actNode] = xiDx[actNode][1];
      }
    }

    
  }


  // ============
  //   calcDMat
  // ============
  void linElecInt::calcDMat( Matrix<Double> &dMat ) {

    ENTER_FCN( "linElecInt::calcDMat" );
    ptMaterial->GetTensor(dMat,ELEC_PERMITTIVITY,matDataType_,subTensorType_);
    dMat *= mParser_->Eval( mHandle_ );
  }


  void linElecInt::SetFactor( const std::string& factor ) {
    ENTER_FCN( "linElecInt::SetFactor" );
    
    mParser_->SetExpr( mHandle_, factor );
  }
  

  // ================
  //   Constructors
  // ================


  linElecInt::linElecInt(BaseMaterial* matData, SubTensorType type,
                         bool coordUpdate ) 
    : BDBInt(matData, type, coordUpdate ) {

    ENTER_FCN( "linElecInt::linElecInt" );

    name_ = "linElecInt";
    if ( type == FULL ) {
      dim_ = 3;
    }
    else {
      dim_ = 2;
    }

    if ( type == AXI ) {
      isaxi_     = true;
    }
  }


}
