#include <iostream>
#include <fstream>

#include "General/environment.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "linPiezoCoupling3D.hh"

namespace CoupledField {


  // ============
  //   calcAMat
  // ============
  void linPiezoCoupling3D::calcAMat( Matrix<Double> &aMat, UInt ip,
                                     const Matrix<Double> &ptCoord ) {

    ENTER_FCN( "linPiezoCoupling3D::calcAMat" );

    // Obtain info on problem sizes
    const UInt numNodes = ptelem->GetNumNodes();
    const UInt nDofMech = 3;
    const UInt nRowsD   = 6;

    // Set correct size of matrix A and initialise with zeros
    aMat.Resize( numNodes * nDofMech, nRowsD );

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


  // ============
  //   calcBMat
  // ============
  void linPiezoCoupling3D::calcBMat( Matrix<Double> &bMat, UInt ip,
                                     const Matrix<Double> &ptCoord ) {

    ENTER_FCN( "linPiezoCoupling3D::calcBMat" );

    // Obtain info on number of element's nodes
    const UInt numNodes = ptelem->GetNumNodes();

    // Set correct size of matrix B and initialise with zeros
    bMat.Resize( 3, numNodes );

    // Get derivatives of local shape functions with respect to global coords
    // (format: nrNodes x spaceDim)
    Matrix<Double> xiDx;
    ptelem->GetGlobDerivShFncAtIp( xiDx, ip, ptCoord );

    // The matrix bMat can be seen as a 1 x numNodes block-vector.
    // The k-th entry of this block vector corresponds to the matrix
    // B of the ADB product evaluated at the k-th node of the finite
    // element. We simply must transpose xiDx.
    for( UInt actNode = 0; actNode < numNodes; actNode++ ) {
      bMat[0][actNode] = xiDx[actNode][0];
      bMat[1][actNode] = xiDx[actNode][1];
      bMat[2][actNode] = xiDx[actNode][2];
    }
  }


  // ============
  //   calcDMat
  // ============
  void linPiezoCoupling3D::calcDMat( Matrix<Double> &dMat ) {

    ENTER_FCN( "linPiezoCoupling3D::calcDMat" );

    dMat.Resize( 6, 3 );
    Matrix<Double> *matMatrix = ptMaterial->GetMatrix();

    // Extract e^T part from piezo-material matrix
    for( UInt i = 0; i < 6; i++ ) {
      for ( UInt j = 0; j < 3; j++ ) {
        dMat[i][j] = (*matMatrix)[i][j+6];
      }
    }

  }

}
