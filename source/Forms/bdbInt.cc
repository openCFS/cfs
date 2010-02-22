// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include "bdbInt.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "Optimization/DesignElement.hh"
#include "Elements/integrationScheme.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"
#include "Elements/H1Elems.hh"
#include "Elements/fespace.hh"

namespace CoupledField {

  void BDBInt::CalcElementMatrix( Matrix<Double>& elemMat,
                                  EntityIterator& ent1,
                                  EntityIterator& ent2,
                                  DesignElement::Type direction ) {
    
    // Extract physical element
    const Elem* ptElem = ent1.GetElem();

    Matrix<Double> dMat, bMat, dbMat;
    Double fac = 0.0;

    // Obtain FE element from feSpace
    BaseFE* ptFe = ptFeSpace1_->GetFe( ent1 ); 
    UInt nrFncs = ptFe->GetNumFncs();

    // Get shape map from grid
    shared_ptr<ElemShapeMap> esm = domain->GetGrid()->GetElemShapeMap( ptElem );

    // Get integration points (shortcut: from basefe instead of 
    // IntegrationScheme class)
    StdVector<LocPoint> intPoints;
    StdVector<Double> weights;
    intScheme_->GetIntPoints( ptElem->type, intPoints, weights );
    
    const UInt nrDofs   = getNrDofs();
    elemMat.Resize( nrFncs * nrDofs);
    elemMat.Init();
    
    // Loop over all integration points
    LocPointMapped lp;
    for( UInt i = 0; i < intPoints.GetSize(); i++  ) {

      // Calculate for each integration point the LocPointMapped
      // ... implement me
      lp.Set( intPoints[i], esm );

      // Call the CalcBMat()-method
      calcBMat( bMat, lp, ptFe);
      
      // Calculate D-Mat
      calcDMat(dMat, ent1.GetElem());

      fac = lp.jacDet * weights[i];
      // Compute the matrix product D * B and store as intermediate matrix
      // resize dbMat to handle SurfaceNortmalInt
      dbMat.Resize(bMat.GetNumRows(), bMat.GetNumCols());
      dbMat = (dMat * bMat) * fac;
      elemMat += Transpose(bMat) * dbMat;
    }

//    // Extract pointer to reference element and get coordinates
//    ExtractElemInfo( ent1 );
//
//
//    // First of all, set ansatz function to element
//    ptelem->SetAnsatzFct( ansatzFct1_ );
//
//    UInt nrFncs = ptelem->GetNumFncs( ansatzFct1_ );
//    const UInt nrDofs   = getNrDofs();
//
//    double jacDet;
//
//    Matrix<Double> bMat;
//    Matrix<Double> dMat;
//    Matrix<Double> dbMat;
//    Double aux, fac, *ptr1, *ptr2;
//
//    elemMat.Resize( nrFncs * nrDofs);
//    elemMat.Init();
//
//    dbMat.Resize( getDimD(), nrFncs * nrDofs);
//
//    //if softening, get maximal/minimal edge lenght
//    if ( softeningPart_ == "bendingBK1" ) {
//      ptelem->GetMaxMinEdgeLength(ptCoord_,maxEdgeLength_,minEdgeLength_);
//    }
//
//    if ( softeningPart_ == "shearBK1" ||  softeningPart_ == "shearSRI" ) {
//      //do reduced order of integration
//      ptelem->SetReducedIntegration();
//    }
//
//    //get integration points
//
//    const UInt nrIntPts = ptelem->GetNumIntPoints();
//    const Vector<Double> & intWeights = ptelem->GetIntWeights();
//
//    // **************************************************
//    //  Material matrix independent of integration point
//    // **************************************************
//    if ( updateDMatInEveryIP_ == false ) {
//
//
//      // // Check if material has to be rotated
//      if( ptMaterial->GetCoordSys() == NULL ) {
//        if(direction == DesignElement::NO_DERIVATIVE){
//          calcDMat(dMat, ent1.GetElem());
//        }else{
//          calcDMat(dMat, ent1.GetElem(), direction);
//        }
//      }
//
//      // Loop over all integration points
//      for ( UInt actIntPt = 1; actIntPt <= nrIntPts; actIntPt++ ) {
//
//        // Check if material has to be rotated for each integration point
//        if( ptMaterial->GetCoordSys() != NULL ) {
//          // Get global coordinates
//          Vector<Double> * intPoints = ptelem->GetIntPoints();
//          Vector<Double> globIntPoint;
//
//          ptelem->Local2GlobalCoord(globIntPoint, intPoints[actIntPt-1],
//                                    ptCoord_, ent1.GetElem() );
//          ptMaterial->RotateTensorByPointCoord( globIntPoint,getDMaterialType() );
//          if(direction == DesignElement::NO_DERIVATIVE){
//            calcDMat(dMat, ent1.GetElem());
//          }else{
//            calcDMat(dMat, ent1.GetElem(), direction);
//          }
//        }
//
//
//        // Setup the B matrix for current integration point
//        calcBMat( bMat, actIntPt, ptCoord_ );
//
//        // Compute Jacobian for integration point
//        jacDet = ptelem->CalcJacobianDetAtIp( actIntPt, ptCoord_, ent1.GetElem() );
//
//        // Perform a safety check
//        if ( jacDet < 0.0 ) {
//          EXCEPTION( "BDBInt::CalcElementMatrix: Encountered "
//                   << "negative Jacobian determinant!" );
//        }
//
//        // Special things must be done in the axi-symmetric case
//        if ( isaxi_ ) {
//          Vector<Double> ShpFncAtIp;
//          Vector<Double> CoordAtIP;
//          ptelem->GetShFncAtIp( ShpFncAtIp, actIntPt, ent1.GetElem() );
//
//          CoordAtIP = ptCoord_ * ShpFncAtIp;
//          jacDet *= 2 * PI * CoordAtIP[0];
//        }
//
//        // Compute the matrix product D * B and store as intermediate matrix
//        // resize dbMat to handle SurfaceNortmalInt
//        dbMat.Resize(bMat.GetNumRows(), bMat.GetNumCols());
//        dMat.Mult( bMat, dbMat );
//
//        // We now compute B^T * D * B and scale it by the determinant
//        // of the Jacobian and the weight of the current integration
//        // point. The result is added to the element matrix.
//        fac = jacDet * intWeights[actIntPt-1];
//        for ( UInt k = 0; k < bMat.GetNumRows(); k++ ) {
//          ptr1 =  bMat[k];
//          ptr2 = dbMat[k];
//          for ( UInt i = 0; i < bMat.GetNumCols(); i++ ) {
//            aux = fac * ptr1[i];
//            for ( UInt j = 0; j < dbMat.GetNumCols(); j++ ) {
//              elemMat[i][j] += aux * ptr2[j];
//            }
//          }
//        }
//      }
//
//    }
//
//    // **********************************************
//    //  Material matrix depends on integration point
//    // **********************************************
//    else {
//
//      // Loop over all integration points
//      for ( UInt actIntPt = 1; actIntPt <= nrIntPts; actIntPt++ ) {
//
//        // Setup material matrix for current integration point
//        calcDMat(dMat, actIntPt, ptCoord_);
//
//        // Setup the B matrix for current integration point
//        calcBMat( bMat, actIntPt, ptCoord_ );
//
//        // Compute Jacobian for integration point
//        jacDet = ptelem->CalcJacobianDetAtIp( actIntPt, ptCoord_, ent1.GetElem() );
//
//        // Perform a safety check
//        if ( jacDet < 0.0 ) {
//          EXCEPTION( "BDBInt::CalcElementMatrix: Encountered "
//                   << "negative Jacobian determinant!" );
//        }
//
//        // Special things must be done in the axi-symmetric case
//        if ( isaxi_ ) {
//          Vector<Double> ShpFncAtIp;
//          Vector<Double> CoordAtIP;
//          ptelem->GetShFncAtIp( ShpFncAtIp, actIntPt, ent1.GetElem() );
//
//          CoordAtIP = ptCoord_ * ShpFncAtIp;
//          jacDet *= 2 * PI * CoordAtIP[0];
//        }
//
//        // Compute the matrix product D * B and store as intermediate matrix
//        dbMat.Resize( dMat.GetNumRows(), bMat.GetNumCols() );
//        dMat.Mult( bMat, dbMat );
//
//        // We now compute B^T * D * B and scale it by the determinant
//        // of the Jacobian and the weight of the current integration
//        // point. The result is added to the element matrix.
//        fac = jacDet * intWeights[actIntPt-1];
//        for ( UInt k = 0; k < bMat.GetNumRows(); k++ ) {
//          ptr1 =  bMat[k];
//          ptr2 = dbMat[k];
//          for ( UInt i = 0; i < bMat.GetNumCols(); i++ ) {
//            aux = fac * ptr1[i];
//            for ( UInt j = 0; j < dbMat.GetNumCols(); j++ ) {
//              elemMat[i][j] += aux * ptr2[j];
//            }
//          }
//        }
//      }
//    }
//
//    if ( softeningPart_ == "shearSRI" || softeningPart_ == "shearBK1" ) {
//      //set back to standard integration
//      ptelem->SetStandardIntegration();
//    }
  }


  

  // ****************************
  //   CalcComplexElementMatrix
  // ****************************
  void BDBInt::CalcComplexElementMatrix( Matrix<Complex> & elemMat,
                                          EntityIterator& ent1,
                                          EntityIterator& ent2,
                                          Double & beta, Double & omega) {

      EXCEPTION("Implement me");
//
//    // Get pointer to reference element and its coordinates
//    ExtractElemInfo( ent1 );
//
//    const UInt nrIntPts = ptelem->GetNumIntPoints();
//    const UInt nrNodes  = ptelem->GetNumNodes();
//    const UInt nrDofs   = getNrDofs();
//    const Vector<Double> & intWeights = ptelem->GetIntWeights();
//    double jacDet;
//
//    Matrix<Double> bMat;
//    Matrix<Complex> dMat;
//    Matrix<Complex> dB;
//    Matrix<Double> bTrans;
//    Matrix<Complex> partElemMat;
//
//    elemMat.Resize(nrNodes * nrDofs,true);
//    elemMat.Init();
//
//    calcDMaterialMatWithComplexDamping( dMat, beta, omega );
//
//    for ( UInt actIntPt = 1; actIntPt <= nrIntPts; actIntPt++ ) {
//
//      if (updateDMatInEveryIP_) {
//        calcDMaterialMatWithComplexDamping(dMat,beta,omega);
//      }
//
//      calcBMat(bMat, actIntPt, ptCoord_);
//
//      //   hardcoded dB = dMat * bMat;
//      dB.Resize(dMat.GetNumRows(), bMat.GetNumCols());
//      Complex a;
//
//      for ( UInt i = 0; i < dMat.GetNumRows(); i++ ) {
//        for ( UInt j = 0; j < bMat.GetNumCols(); j++ ) {
//          a = dMat[i][0] * bMat[0][j];
//          for ( UInt k = 1; k < bMat.GetNumRows(); k++ ) {
//            a += dMat[i][k] * bMat[k][j];
//          }
//          dB[i][j] = a;
//        }
//      }
//
//      bMat.Transpose(bTrans);
//
//      // hardcoded: partElemMat = bTrans * dB;
//      partElemMat.Resize(bTrans.GetNumRows(), dB.GetNumCols());
//      for ( UInt i = 0; i < bTrans.GetNumRows(); i++ ) {
//        for ( UInt j = 0; j < dB.GetNumCols(); j++ ) {
//          a = bTrans[i][0] *dB[0][j];
//          for ( UInt k = 1; k < dB.GetNumRows(); k++ ) {
//            a += bTrans[i][k] * dB[k][j];
//          }
//          partElemMat[i][j] = a;
//        }
//      }
//
//      jacDet = ptelem->CalcJacobianDetAtIp(actIntPt,ptCoord_, ent1.GetElem());
//
//      // Perform a safety check
//      if ( jacDet < 0.0 ) {
//        EXCEPTION( "BDBInt::CalcElementMatrix: Encountered "
//                 << "negative Jacobian determinant!" );
//      }
//
//      if (isaxi_) {
//        Vector<Double> ShpFncAtIp;
//        Vector<Double> CoordAtIP;
//        ptelem->GetShFncAtIp(ShpFncAtIp, actIntPt, ent1.GetElem());
//
//        CoordAtIP = ptCoord_ * ShpFncAtIp;
//        jacDet *= 2 * PI * CoordAtIP[0];
//      }
//
//      for ( UInt i = 0; i < elemMat.GetNumRows(); i++ ) {
//        for ( UInt j = 0; j < elemMat.GetNumCols(); j++ ) {
//          elemMat[i][j] += partElemMat[i][j] * jacDet *
//            intWeights[actIntPt-1] ;
//        }
//      }
//    }
  }


  void BDBInt::calcBMat(EntityIterator it, Matrix<Double>& bMat ) {

      EXCEPTION("Implement me");
//    // get midpoint
//    ExtractElemInfo( it );
//    Vector<Double> midPoint;
//    it.GetElem()->ptElem->GetCoordMidPoint(midPoint);
//
//    // Set integration point to midpont
//    SetIntPoint( midPoint);
//    calcBMat( bMat, 1, ptCoord_ );
//    UnsetIntPoint();
  }

  void BDBInt::calcDBMat(EntityIterator it, Matrix<Double>& dbMat ) {
    EXCEPTION( "Implement me");
//    // get midpoint
//    ExtractElemInfo( it );
//    Vector<Double> midPoint;
//    it.GetElem()->ptElem->GetCoordMidPoint(midPoint);
//
//    // Set integration point to midpont
//    SetIntPoint( midPoint);
//    Matrix<Double> temp1 , temp2;
//    calcBMat( temp1, 1, ptCoord_ );
//    calcDMat( temp2, 1, ptCoord_ );
//    dbMat = temp1*temp2;
//    UnsetIntPoint();
  }

  // ***************
  //   Constructor
  // ***************
  BDBInt::BDBInt( BaseMaterial* matData, SubTensorType type,
                  bool geoUpdate )
    : BaseForm(matData,type, geoUpdate ), updateDMatInEveryIP_(0) {
    name_ = "BDBInt";
    baseType_ = STIFFNESS;
  }


  // **************
  //   Destructor
  // **************
  BDBInt::~BDBInt() {
  }


}
