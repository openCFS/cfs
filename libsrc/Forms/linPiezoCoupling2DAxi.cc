#include <iostream>
#include <fstream>

#include "General/environment.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "linPiezoCoupling2DAxi.hh"

namespace CoupledField {


  // ============
  //   calcAMat
  // ============
  void linPiezoCoupling2DAxi::calcAMat( Matrix<Double> &aMat, Integer ip,
                                     const Matrix<Double> &ptCoord ) {

    ENTER_FCN( "linPiezoCoupling2DAxi::calcAMat" );

    // Obtain info on problem sizes
    const Integer numNodes = ptelem->GetNumNodes();
    const Integer nDofMech = 2;
    const Integer nRowsD   = 4;

    // Set correct size of matrix A and initialise with zeros
    aMat.Resize( numNodes * nDofMech, nRowsD );

    // Get local shape functions and their derivatives with respect to global
    // coords (format for the latter: nrNodes x spaceDim)
    Matrix<Double> xiDx;
    Vector<Double> ShpFncAtIp;
    ptelem->GetGlobDerivShFncAtIp( xiDx, ip, ptCoord );
    ptelem->GetShFncAtIp( ShpFncAtIp, ip );

    // The matrix aMat can be seen as a numNodes x 1 block-vector.
    // The k-th entry of this block vector corresponds to the matrix
    // A of the ADB product evaluated at the k-th node of the finite
    // element. We assemble aMat in a top-down fashion and assume that
    // the first coordinate equals r and the second z.
    Integer actInd = 0;
    Integer j;
    Integer actNode;
    Double coordAtIp;

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


  // ============
  //   calcBMat
  // ============
  void linPiezoCoupling2DAxi::calcBMat( Matrix<Double> &bMat, Integer ip,
                                     const Matrix<Double> &ptCoord ) {

    ENTER_FCN( "linPiezoCoupling2DAxi::calcBMat" );

    // Obtain info on number of element's nodes
    const Integer numNodes = ptelem->GetNumNodes();

    // Set correct size of matrix B and initialise with zeros
    bMat.Resize( 2, numNodes );

    // Get derivatives of local shape functions with respect to global coords
    // (format: nrNodes x spaceDim)
    Matrix<Double> xiDx;
    ptelem->GetGlobDerivShFncAtIp( xiDx, ip, ptCoord );

    // The matrix bMat can be seen as a 1 x numNodes block-vector.
    // The k-th entry of this block vector corresponds to the matrix
    // B of the ADB product evaluated at the k-th node of the finite
    // element. We assume that the first coordinate equals r and the
    // second z.
    for( Integer actNode = 0; actNode < numNodes; actNode++ ) {
      bMat[0][actNode] = xiDx[actNode][0];
      bMat[1][actNode] = xiDx[actNode][1];
    }
  }


  // ============
  //   calcDMat
  // ============
  void linPiezoCoupling2DAxi::calcDMat( Matrix<Double> &dMat ) {

    ENTER_FCN( "linPiezoCoupling2DAxi::calcDMat" );

    dMat.Resize( 4, 2 );
    Matrix<Double> *matMatrix = ptMaterial->GetMatrix();

    dMat[0][0] = (*matMatrix)[1][7];
    dMat[1][0] = (*matMatrix)[2][7];
    dMat[2][0] = (*matMatrix)[3][7];
    dMat[3][0] = (*matMatrix)[0][7];

    dMat[0][1] = (*matMatrix)[1][8];
    dMat[1][1] = (*matMatrix)[2][8];
    dMat[2][1] = (*matMatrix)[3][8];
    dMat[3][1] = (*matMatrix)[0][8];
  }

}
