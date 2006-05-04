#include <iostream>
#include <fstream>

#include "General/environment.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "linElecInt.hh"

namespace CoupledField {


  // ============
  //   calcBMat
  // ============
  void linElecInt::calcBMat( Matrix<Double> &bMat, UInt ip,
			     Matrix<Double> &ptCoord ) {

    ENTER_FCN( "linElecInt::calcAMat" );

    // Obtain info on number of element's nodes
    const UInt numNodes = ptelem->GetNumNodes();

    // Set correct size of matrix B and initialise with zeros
    bMat.Resize( dim_, numNodes );

    // Get derivatives of local shape functions with respect to global
    // coords (format: nrNodes x spaceDim)
    Matrix<Double> xiDx;
    ptelem->GetGlobDerivShFncAtIp( xiDx, ip, ptCoord );

    if ( subTensorType_ == FULL ) {
      // The matrix bMat can be seen as a 3 x numNodes block-vector.
      // The k-th entry of this block vector corresponds to the matrix
      // B of the BDB product evaluated at the k-th node of the finite
      // element. 
      for( UInt actNode = 0; actNode < numNodes; actNode++ ) {
	bMat[0][actNode] = xiDx[actNode][0];
	bMat[1][actNode] = xiDx[actNode][1];
	bMat[2][actNode] = xiDx[actNode][2];
      }
    }
    else  {
      // The matrix bMat can be seen as a 1 x numNodes block-vector.
      // The k-th entry of this block vector corresponds to the matrix
      // B of the ADB product evaluated at the k-th node of the finite
      // element. We assume that the first coordinate equals y and the
      // second z.
      for( UInt actNode = 0; actNode < numNodes; actNode++ ) {
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
    dMat *= factor_;
  }



  // ================
  //   Constructors
  // ================


  linElecInt::linElecInt(SubTensorType type) {
    ENTER_FCN( "linElecInt::linElecInt" );

    subTensorType_ = type;
    factor_  = 1.0;

    if ( type == FULL ) {
      dim_ = 3;
    }
    else {
      dim_ = 2;
    }

    if ( type == AXI ) {
      isaxi_     = TRUE;
    }
  }

  linElecInt::linElecInt(BaseMaterial* matData, SubTensorType type) 
             : BDBInt(matData) {

    ENTER_FCN( "linElecInt::linElecInt" );

    subTensorType_ = type;
    factor_ = 1.0;

    if ( type == FULL ) {
      dim_ = 3;
    }
    else {
      dim_ = 2;
    }

    if ( type == AXI ) {
      isaxi_     = TRUE;
    }
  }


}
