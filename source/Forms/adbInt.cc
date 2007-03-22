// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include "adbInt.hh"


namespace CoupledField {


  // =====================
  //   CalcElementMatrix
  // =====================
  void ADBInt::CalcElementMatrix( Matrix<Double>& elemMat,
                                  EntityIterator& ent1, 
                                  EntityIterator& ent2 ) {
    ENTER_FCN( "ADBInt::CalcElementMatrix" );

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );

    ptelem->SetAnsatzFct( ansatzFct1_ );
    UInt numFncs1 = ptelem->GetNumFncs( ansatzFct1_ );
    UInt numFncs2 = ptelem->GetNumFncs( ansatzFct2_ );
    
    const UInt nrIntPts = ptelem->GetNumIntPoints(); 
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    double jacDet;

    Matrix<Double> aMat;
    Matrix<Double> bMat;
    Matrix<Double> dMat;
    Matrix<Double> dbMat;
    Double aux;

    elemMat.Resize( numFncs1 * getNumDofsA(), numFncs2 * getNumDofsB() );
    elemMat.Init();


    // **************************************************
    //  Material matrix independent of integration point
    // **************************************************

    // Setup material matrix once and for all
    calcDMat( dMat );

    // Loop over all integration points
    for ( UInt actIntPt = 1; actIntPt <= nrIntPts; actIntPt++ ) {
      
      // Check if material has to be rotated for each integration point
      if( ptMaterial->GetCoordSys() != NULL ) {
        // Get global coordinates
        Vector<Double> * intPoints = ptelem->GetIntPoints();
        Vector<Double> globIntPoint;
        
        ptelem->Local2GlobalCoord(globIntPoint, intPoints[actIntPt-1], 
                                  ptCoord_, ent1.GetElem() );
        ptMaterial->RotateTensorByPointCoord( globIntPoint, getDMaterialType() );
        calcDMat( dMat );
      }
      
      //std::cerr << "*** Calculating A ****\n";
      // Setup the A matrix for current integration point
      calcAMat( aMat, actIntPt, ptCoord_ );

      //std::cerr << "*** Calculating B ****\n";
      // Setup the B matrix for current integration point
      calcBMat( bMat, actIntPt, ptCoord_ );

      // Compute Jacobian for integration point
      jacDet = ptelem->CalcJacobianDetAtIp( actIntPt, ptCoord_, 
                                            ent1.GetElem() );

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
        ptelem->GetShFncAtIp( ShpFncAtIp, actIntPt, ent1.GetElem() );
        Double aux = 0.0;
        
        for ( UInt i = 0; i < numFncs1; i++ ) {
          aux += ptCoord_[0][i] * ShpFncAtIp[i];
        }
        
        jacDet *= 2.0 * PI * aux;
      }

      // Compute the matrix product D * B and store as intermediate matrix
      dbMat.Resize( dMat.GetSizeRow(), bMat.GetSizeCol() );
      dMat.Mult( bMat, dbMat );

      // We now compute A * D * B and scale it by the determinant
      // of the Jacobian and the weight of the current integration
      // point. The result is added to the element matrix
      for ( UInt i = 0; i < aMat.GetSizeRow(); i++ ) {
        for ( UInt j = 0; j < dbMat.GetSizeCol(); j++ ) {

          // Compute entry (i,j) of A * D * B
          aux = 0.0;
          for ( UInt k = 0; k < aMat.GetSizeCol(); k++ ) {
            aux += aMat[i][k] * dbMat[k][j];
          }

          // Scale result and add to corresponding entry
          // of element matrix
          elemMat[i][j] += aux * jacDet * intWeights[actIntPt-1];
        }
      }
    }


    //std::cerr << "elemMat = \n" << elemMat << std::endl;
  }

}
