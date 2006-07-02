#include <iostream>
#include <fstream>

#include "bdbInt.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"

namespace CoupledField {

  // New version seems to be buggy
#define BDB_NEW_VERSION


#ifndef BDB_NEW_VERSION

  // *********************
  //   CalcElementMatrix
  // *********************


  void BDBInt::CalcElementMatrix( Matrix<Double> &ptCoord,
                                  Matrix<Double> &elemMat ) {

    ENTER_FCN( "BDBInt::CalcElementMatrix" );

    const UInt nrIntPts = ptelem->GetNumIntPoints(); 
    const UInt nrNodes  = ptelem->GetNumNodes();   
    const UInt nrDofs   = getNrDofs();  
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    double jacDet;

    Matrix<Double> bMat; 
    Matrix<Double> dMat; 
    Matrix<Double> dB; 
    Matrix<Double> bTrans; 
    Matrix<Double> partElemMat;
  
    elemMat.Resize(nrNodes * nrDofs);
    elemMat.Init();
    dB.Resize( getDimD(), nrNodes * nrDofs );

    // If the material parameters are constant within the element
    // we can compute the D matrix once and for all
    if ( updateDMatInEveryIP_ == false ) {

      calcDMat( dMat );
    }

    // Loop over all integration points
    for ( UInt actIntPt = 1; actIntPt <= nrIntPts; actIntPt++ ) {

      // Check if D matrix must be re-determined for
      // the current integration point
      if ( updateDMatInEveryIP_ == true )
        calcDMat( dMat, actIntPt, ptCoord );
      }

      // Setup the B matrix for this integration point
      calcBMat( bMat, actIntPt, ptCoord );

      dMat.Mult( bMat, dB );
      // dB = dMat * bMat;
      bMat.Transpose(bTrans);
      partElemMat = bTrans * dB;

      jacDet = ptelem->CalcJacobianDetAtIp( actIntPt, ptCoord );

      // Perform a safety check
      if ( jacDet < 0.0 ) {
        (*error) << "BDBInt::CalcElementMatrix: Encountered "
                 << "negative Jacobian determinant!";
        Error( __FILE__, __LINE__ );
      }

      // Special things must be done in the axi-symmetric case
      if ( isaxi_ ) {
        Vector<Double> ShpFncAtIp;
        Vector<Double> CoordAtIP;
        ptelem->GetShFncAtIp( ShpFncAtIp, actIntPt );

        CoordAtIP = ptCoord * ShpFncAtIp;
        jacDet *= 2 * PI * CoordAtIP[0];
      }

      elemMat += partElemMat * jacDet * intWeights[actIntPt-1];
    }

  }

#else
  void BDBInt::CalcElementMatrix( Matrix<Double>& elemMat,
                                  EntityIterator& ent1, 
                                  EntityIterator& ent2 ) {

    ENTER_FCN( "BDBInt::CalcElementMatrix" );

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );

    const UInt nrNodes  = ptelem->GetNumNodes();   
    const UInt nrDofs   = getNrDofs();  

    double jacDet;

    Matrix<Double> bMat; 
    Matrix<Double> dMat; 
    Matrix<Double> dbMat; 
    Double aux, fac, *ptr1, *ptr2;

    elemMat.Resize( nrNodes * nrDofs);
    elemMat.Init();

    dbMat.Resize( getDimD(), nrNodes * nrDofs);

    //if softening, get maximal/minimal edge lenght
    if ( softeningPart_ == "bendingBK1" ) {
      ptelem->GetMaxMinEdgeLength(ptCoord_,maxEdgeLength_,minEdgeLength_);
    }

    if ( softeningPart_ == "shearBK1" ||  softeningPart_ == "shearSRI" ) {
      //do reduced order of integration
      ptelem->SetReducedIntegration();
    }

    //get integration points
    const UInt nrIntPts = ptelem->GetNumIntPoints(); 
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  

    // **************************************************
    //  Material matrix independent of integration point
    // **************************************************
    if ( updateDMatInEveryIP_ == false ) {
      

      // // Check if material has to be rotated
      if( ptMaterial->GetCoordSys() == NULL ) {
        calcDMat( dMat );
      }
      // Loop over all integration points
      for ( UInt actIntPt = 1; actIntPt <= nrIntPts; actIntPt++ ) {

        // Check if material has to be rotated for each integration point
        if( ptMaterial->GetCoordSys() != NULL ) {
          // Get global coordinates
          Vector<Double> * intPoints = ptelem->GetIntPoints();
          Vector<Double> globIntPoint;
          
          ptelem->Local2GlobalCoord(globIntPoint, intPoints[actIntPt-1], ptCoord_);
          ptMaterial->RotateTensorByPointCoord( globIntPoint,MECH_STIFFNESS_TENSOR );
          calcDMat( dMat );
        }
        

        // Setup the B matrix for current integration point
        calcBMat( bMat, actIntPt, ptCoord_ );

        // Compute Jacobian for integration point
        jacDet = ptelem->CalcJacobianDetAtIp( actIntPt, ptCoord_ );

        // Perform a safety check
        if ( jacDet < 0.0 ) {
          (*error) << "BDBInt::CalcElementMatrix: Encountered "
                   << "negative Jacobian determinant!";
          Error( __FILE__, __LINE__ );
        }

        // Special things must be done in the axi-symmetric case
        if ( isaxi_ ) {
          Vector<Double> ShpFncAtIp;
          Vector<Double> CoordAtIP;
          ptelem->GetShFncAtIp( ShpFncAtIp, actIntPt );

          CoordAtIP = ptCoord_ * ShpFncAtIp;
          jacDet *= 2 * PI * CoordAtIP[0];
        }

        // Compute the matrix product D * B and store as intermediate matrix
        dMat.Mult( bMat, dbMat );

        // We now compute B^T * D * B and scale it by the determinant
        // of the Jacobian and the weight of the current integration
        // point. The result is added to the element matrix.
        fac = jacDet * intWeights[actIntPt-1];
        for ( UInt k = 0; k < bMat.GetSizeRow(); k++ ) {
          ptr1 =  bMat[k];
          ptr2 = dbMat[k];
          for ( UInt i = 0; i < bMat.GetSizeCol(); i++ ) {
            aux = fac * ptr1[i];
            for ( UInt j = 0; j < dbMat.GetSizeCol(); j++ ) {
              elemMat[i][j] += aux * ptr2[j];
            }
          }
        }
      }
    }


    // **********************************************
    //  Material matrix depends on integration point
    // **********************************************
    else {

      // Loop over all integration points
      for ( UInt actIntPt = 1; actIntPt <= nrIntPts; actIntPt++ ) {

        // Setup material matrix for current integration point
        calcDMat( dMat, actIntPt, ptCoord_ );

        // Setup the B matrix for current integration point
        calcBMat( bMat, actIntPt, ptCoord_ );

        // Compute Jacobian for integration point
        jacDet = ptelem->CalcJacobianDetAtIp( actIntPt, ptCoord_ );

        // Perform a safety check
        if ( jacDet < 0.0 ) {
          (*error) << "BDBInt::CalcElementMatrix: Encountered "
                   << "negative Jacobian determinant!";
          Error( __FILE__, __LINE__ );
        }

        // Special things must be done in the axi-symmetric case
        if ( isaxi_ ) {
          Vector<Double> ShpFncAtIp;
          Vector<Double> CoordAtIP;
          ptelem->GetShFncAtIp( ShpFncAtIp, actIntPt );

          CoordAtIP = ptCoord_ * ShpFncAtIp;
          jacDet *= 2 * PI * CoordAtIP[0];
        }

        // Compute the matrix product D * B and store as intermediate matrix
        dbMat.Resize( dMat.GetSizeRow(), bMat.GetSizeCol() );
        dMat.Mult( bMat, dbMat );

        // We now compute B^T * D * B and scale it by the determinant
        // of the Jacobian and the weight of the current integration
        // point. The result is added to the element matrix.
        fac = jacDet * intWeights[actIntPt-1];
        for ( UInt k = 0; k < bMat.GetSizeRow(); k++ ) {
          ptr1 =  bMat[k];
          ptr2 = dbMat[k];
          for ( UInt i = 0; i < bMat.GetSizeCol(); i++ ) {
            aux = fac * ptr1[i];
            for ( UInt j = 0; j < dbMat.GetSizeCol(); j++ ) {
              elemMat[i][j] += aux * ptr2[j];
            }
          }
        }
      }
    }

    if ( softeningPart_ == "shearSRI" || softeningPart_ == "shearBK1" ) {
      //set back to standard integration
      ptelem->SetStandardIntegration();
    }

  }

#endif



  // ****************************
  //   CalcComplexElementMatrix
  // ****************************
  void BDBInt:: CalcComplexElementMatrix( Matrix<Complex> & elemMat,
                                          EntityIterator& ent1, 
                                          EntityIterator& ent2,
                                          Double & beta, Double & omega) {
    
    ENTER_FCN( "BDBInt::CalcComplexElementMatrix" );
    
    // Get pointer to reference element and its coordinates
    ExtractElemInfo( ent1 );
    
    const UInt nrIntPts = ptelem->GetNumIntPoints(); 
    const UInt nrNodes  = ptelem->GetNumNodes();   
    const UInt nrDofs   = getNrDofs();  
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    double jacDet;
    
    Matrix<Double> bMat; 
    Matrix<Complex> dMat; 
    Matrix<Complex> dB; 
    Matrix<Double> bTrans; 
    Matrix<Complex> partElemMat;
    
    elemMat.Resize(nrNodes * nrDofs,true);
    elemMat.Init(); 
    
    calcDMaterialMatWithComplexDamping( dMat, beta, omega );
    
    for ( UInt actIntPt = 1; actIntPt <= nrIntPts; actIntPt++ ) {
      
      if (updateDMatInEveryIP_) {
        calcDMaterialMatWithComplexDamping(dMat,beta,omega);
      }
      
      calcBMat(bMat, actIntPt, ptCoord_);
      
      //   hardcoded dB = dMat * bMat;
      dB.Resize(dMat.GetSizeRow(), bMat.GetSizeCol());
      Complex a;
      
      for ( UInt i = 0; i < dMat.GetSizeRow(); i++ ) {
        for ( UInt j = 0; j < bMat.GetSizeCol(); j++ ) {       
          a = dMat[i][0] * bMat[0][j];
          for ( UInt k = 1; k < bMat.GetSizeRow(); k++ ) {
            a += dMat[i][k] * bMat[k][j];
          }
          dB[i][j] = a;
        }
      }

      bMat.Transpose(bTrans);

      // hardcoded: partElemMat = bTrans * dB;
      partElemMat.Resize(bTrans.GetSizeRow(), dB.GetSizeCol());
      for ( UInt i = 0; i < bTrans.GetSizeRow(); i++ ) {
        for ( UInt j = 0; j < dB.GetSizeCol(); j++ ) {       
          a = bTrans[i][0] *dB[0][j];
          for ( UInt k = 1; k < dB.GetSizeRow(); k++ ) {
            a += bTrans[i][k] * dB[k][j];
          }
          partElemMat[i][j] = a;
        }
      }

      jacDet = ptelem->CalcJacobianDetAtIp(actIntPt,ptCoord_);

      // Perform a safety check
      if ( jacDet < 0.0 ) {
        (*error) << "BDBInt::CalcElementMatrix: Encountered "
                 << "negative Jacobian determinant!";
        Error( __FILE__, __LINE__ );
      }

      if (isaxi_) {
        Vector<Double> ShpFncAtIp;
        Vector<Double> CoordAtIP;
        ptelem->GetShFncAtIp(ShpFncAtIp,actIntPt);

        CoordAtIP = ptCoord_ * ShpFncAtIp;
        jacDet *= 2 * PI * CoordAtIP[0];
      }

      for ( UInt i = 0; i < elemMat.GetSizeRow(); i++ ) {
        for ( UInt j = 0; j < elemMat.GetSizeCol(); j++ ) {
          elemMat[i][j] += partElemMat[i][j] * jacDet *
            intWeights[actIntPt-1] ;
        }
      }
    }
  }


  void BDBInt::GetDMat(Matrix<Double> &dMat) {
    ENTER_FCN( "BDBInt::GetDMat" );
    calcDMat(dMat);
  }

  void BDBInt::GetBMat(Matrix<Double> &bMat, Matrix<Double> & ptCoord_) {
    ENTER_FCN( "BDBInt::GetBMat" );
    const Integer nrIntPts = ptelem->GetNumIntPoints(); 
   
    for (Integer actIntPt=1; actIntPt<=nrIntPts; actIntPt++) {
      calcBMat(bMat, actIntPt, ptCoord_);     
    }
  }



  // ***************
  //   Constructor
  // ***************
  BDBInt::BDBInt( BaseMaterial* matData, SubTensorType type,
                  bool geoUpdate ) 
    : BaseForm(matData,type, geoUpdate ), updateDMatInEveryIP_(0) {
    ENTER_FCN( "BDBInt::BDBInt" );
    name_ = "BDBInt";
    baseType_ = STIFFNESS;
  }

  
  // **************
  //   Destructor
  // **************
  BDBInt::~BDBInt() {
    ENTER_FCN( "BDBInt::~BDBInt" );
  }


}
