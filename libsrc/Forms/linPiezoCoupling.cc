#include <iostream>
#include <fstream>

#include "General/environment.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "linPiezoCoupling.hh"

namespace CoupledField {


  // ============
  //   calcAMat
  // ============
  void linPiezoCoupling::calcAMat( Matrix<Double> &aMat, UInt ip,
				   const Matrix<Double> &ptCoord ) {

    ENTER_FCN( "linPiezoCoupling::calcAMat" );

    // Obtain info on problem sizes
    const UInt numNodes = ptelem->GetNumNodes();
    //    const UInt nDofMech = 3;
    //    const UInt nRowsD   = 6;

    // Set correct size of matrix A and initialise with zeros
    aMat.Resize( numNodes * numDofsA_, matDimRow_ );

    // Get derivatives of local shape functions with respect to global coords
    // (format: nrNodes x spaceDim)
    Matrix<Double> xiDx;
    ptelem->GetGlobDerivShFncAtIp( xiDx, ip, ptCoord );

    // The matrix aMat can be seen as a numNodes x 1 block-vector.
    // The k-th entry of this block vector corresponds to the matrix
    // A of the ADB product evaluated at the k-th node of the finite
    // element. We assemble aMat in a top-down fashion.
    UInt actInd = 0;
    UInt actNode;

    if ( subTensorType_ == FULL ) {
      for( actNode = 0; actNode < numNodes; actNode++ ) {

	// 1st row of sub-matrix A(actNode)
	aMat[actInd][0] = xiDx[actNode][0];
	aMat[actInd][4] = xiDx[actNode][2];
	aMat[actInd][5] = xiDx[actNode][1];
	actInd++;
	
	// 2nd row of sub-matrix A(actNode)
	aMat[actInd][1] = xiDx[actNode][1];
	aMat[actInd][3] = xiDx[actNode][2];
	aMat[actInd][5] = xiDx[actNode][0];
	actInd++;
	
	// 3rd row of sub-matrix A(actNode)
	aMat[actInd][2] = xiDx[actNode][2];
	aMat[actInd][3] = xiDx[actNode][1];
	aMat[actInd][4] = xiDx[actNode][0];
	actInd++;
      }
    }
    else if ( subTensorType_ == AXI ) {
      // we assume that the first coordinate equals r and 
      // the second z.
      
      UInt j;
      Double coordAtIp;
      Vector<Double> ShpFncAtIp;
      ptelem->GetShFncAtIp( ShpFncAtIp, ip );
      
      for( actNode = 0; actNode < numNodes; actNode++ ) {
	
	// 1st row of sub-matrix A(actNode)
	aMat[actInd][0] = xiDx[actNode][0];    // dN/dr
	aMat[actInd][2] = xiDx[actNode][1];    // dN/dz
	
	// For the entry 1/r things are more complicated
	coordAtIp = 0.0;
	for( j = 0; j < numNodes; j++ ) {
	  coordAtIp += ptCoord[0][j] * ShpFncAtIp[j];
	}
	aMat[actInd][3] = ShpFncAtIp[actNode] / coordAtIp;    // 1/r
	actInd++;
	
	// 2nd row of sub-matrix A(actNode)
	aMat[actInd][1] = xiDx[actNode][1];   // dN/dz
	aMat[actInd][2] = xiDx[actNode][0];   // dN/dr
	
	actInd++;
	
      }
    }

    else if ( subTensorType_ == PLANE_STRAIN ) {
      // we assume that the first coordinate equals y and the second z.
      
      for( actNode = 0; actNode < numNodes; actNode++ ) {
	
	// 1st row of sub-matrix A(actNode)
	aMat[actInd][0] = xiDx[actNode][0];   // dN/dy
	aMat[actInd][2] = xiDx[actNode][1];   // dN/dz
	actInd++;
	
	// 2nd row of sub-matrix A(actNode)
	aMat[actInd][1] = xiDx[actNode][1];    // dN/dz
	aMat[actInd][2] = xiDx[actNode][0];    // dN/dy
	actInd++;
      }
    }
    
  }


  // ============
  //   calcBMat
  // ============
  void linPiezoCoupling::calcBMat( Matrix<Double> &bMat, UInt ip,
				   const Matrix<Double> &ptCoord ) {

    ENTER_FCN( "linPiezoCoupling::calcBMat" );

    // Obtain info on number of element's nodes
    const UInt numNodes = ptelem->GetNumNodes();

    // Set correct size of matrix B and initialise with zeros
    bMat.Resize( numDofsA_, numNodes );

    // Get derivatives of local shape functions with respect to global coords
    // (format: nrNodes x spaceDim)
    Matrix<Double> xiDx;
    ptelem->GetGlobDerivShFncAtIp( xiDx, ip, ptCoord );


    // The matrix bMat can be seen as a 1 x numNodes block-vector.
    // The k-th entry of this block vector corresponds to the matrix
    // B of the ADB product evaluated at the k-th node of the finite
    // element. We simply must transpose xiDx.
    if ( subTensorType_ == FULL ) {
      for( UInt actNode = 0; actNode < numNodes; actNode++ ) {
	bMat[0][actNode] = xiDx[actNode][0];
	bMat[1][actNode] = xiDx[actNode][1];
	bMat[2][actNode] = xiDx[actNode][2];
      }
    }
    else if ( subTensorType_ == AXI ) {
      // element. We assume that the first coordinate equals r and the
      // second z.
      for( UInt actNode = 0; actNode < numNodes; actNode++ ) {
	bMat[0][actNode] = xiDx[actNode][0];
	bMat[1][actNode] = xiDx[actNode][1];
      }
    }
   else if ( subTensorType_ == PLANE_STRAIN ) {
     // We assume that the first coordinate equals y and the
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
  void linPiezoCoupling::calcDMat( Matrix<Double> &dMat ) {

    ENTER_FCN( "linPiezoCoupling::calcDMat" );
    Matrix<Double> matMatrix;
    ptMaterial->GetTensor(matMatrix,PIEZO_TENSOR,matDataType_,subTensorType_);
    matMatrix.Transpose(dMat);
  }


  // ============
  //   Constructor
  // ============
  linPiezoCoupling::linPiezoCoupling(SubTensorType type) {
    ENTER_FCN( "linPiezoCoupling::linPiezoCoupling" );

    subTensorType_ = type;

    if ( type == FULL ) {
      numDofsA_  = 3;
      numDofsB_  = 1;
      matDimRow_ = 6;
      matDimCol_ = 3;
    }
    else if ( type == PLANE_STRAIN || type == PLANE_STRESS ) {
      numDofsA_  = 2;
      numDofsB_  = 1;
      matDimRow_ = 3;
      matDimCol_ = 2;
    }
    else if ( type == AXI ) {
      numDofsA_  = 2;
      numDofsB_  = 1;
      matDimRow_ = 4;
      matDimCol_ = 2;
      isaxi_     = TRUE;
    }
  }


}
