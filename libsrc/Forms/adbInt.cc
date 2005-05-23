#include <iostream>
#include <fstream>

#include "adbInt.hh"


namespace CoupledField {


  // =====================
  //   CalcElementMatrix
  // =====================
  void ADBInt::CalcElementMatrix( Matrix<Double> &ptCoord,
                                  Matrix<Double> &elemMat ) {

    ENTER_FCN( "ADBInt::CalcElementMatrix" );

    const Integer nrIntPts = ptelem->GetNumIntPoints(); 
    const Integer nrNodes  = ptelem->GetNumNodes();   
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    double jacDet;

    Matrix<Double> aMat;
    Matrix<Double> bMat;
    Matrix<Double> dMat;
    Matrix<Double> dbMat;
    Double aux;

    elemMat.Resize( nrNodes * getNumDofsA(), nrNodes * getNumDofsB() );


    // **************************************************
    //  Material matrix independent of integration point
    // **************************************************

    // Setup material matrix once and for all
    calcDMat( dMat );

    // Loop over all integration points
    for ( Integer actIntPt = 1; actIntPt <= nrIntPts; actIntPt++ ) {

      // Setup the A matrix for current integration point
      calcAMat( aMat, actIntPt, ptCoord );

      // Setup the B matrix for current integration point
      calcBMat( bMat, actIntPt, ptCoord );

      // Compute Jacobian for integration point
      jacDet = ptelem->CalcJacobianDetAtIp( actIntPt, ptCoord );

      // Perform a safety check
      if ( jacDet < 0.0 ) {
        (*error) << "ADBInt::CalcElementMatrix: Encountered "
                 << "negative Jacobian determinant!";
        Error( __FILE__, __LINE__ );
      }

      // Special things must be done in the axi-symmetric case
      // We need to additionally scale with 2 pi r.
      //
      // NOTE: We assume here that computation is in they-z plane
      //       with z being the axis of symmetry and that y is
      //       represented by the first co-ordinate, thus
      //       2 pi r = "2 pi x"
      if ( isaxi_ ) {
        Vector<Double> ShpFncAtIp;
        ptelem->GetShFncAtIp( ShpFncAtIp, actIntPt );
        Double aux = 0.0;
        
        for ( Integer i = 0; i < nrNodes; i++ ) {
          aux += ptCoord[0][i] * ShpFncAtIp[i];
        }
        
        jacDet *= 2.0 * PI * aux;
      }

      // Compute the matrix product D * B and store as intermediate matrix
      dbMat.Resize( dMat.GetSizeRow(), bMat.GetSizeCol() );
      dMat.Mult( bMat, dbMat );

      // We now compute A * D * B and scale it by the determinant
      // of the Jacobian and the weight of the current integration
      // point. The result is added to the element matrix
      for ( Integer i = 0; i < aMat.GetSizeRow(); i++ ) {
        for ( Integer j = 0; j < dbMat.GetSizeCol(); j++ ) {

          // Compute entry (i,j) of A * D * B
          aux = 0.0;
          for ( Integer k = 0; k < aMat.GetSizeCol(); k++ ) {
            aux += aMat[i][k] * dbMat[k][j];
          }

          // Scale result and add to corresponding entry
          // of element matrix
          elemMat[i][j] += aux * jacDet * intWeights[actIntPt-1];
        }
      }
    }
  }

}
