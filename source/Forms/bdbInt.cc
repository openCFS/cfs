// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <stddef.h>
#include <complex>
#include <fstream>
#include <string>

#include "Domain/elem.hh"
#include "Domain/entityList.hh"
#include "Elements/basefe.hh"
#include "Forms/baseForm.hh"
#include "MatVec/exprt/xpr2.hh"
#include "MatVec/matrix.hh"
#include "MatVec/vector.hh"
#include "Materials/baseMaterial.hh"
#include "Optimization/Design/DesignElement.hh"
#include "bdbInt.hh"


namespace CoupledField {

  void BDBInt::CalcElementMatrix( Matrix<Double>& elemMat,
                                  EntityIterator& ent1,
                                  EntityIterator& ent2,
                                  DesignElement::Type direction ) {


    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );


    // First of all, set ansatz function to element
    ptelem->SetAnsatzFct( ansatzFct1_ );

    UInt nrFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const UInt nrDofs   = getNrDofs();

    double jacDet;

    Matrix<Double> bMat;
    Matrix<Double> dMat;
    Matrix<Double> dbMat;
    Double aux, fac, *ptr1, *ptr2, *ptr3;

    UInt bRows(getDimD());
    UInt bCols(nrFncs * nrDofs);

    elemMat.Resize(bCols);
    elemMat.Init();

    dbMat.Resize( bRows, bCols);

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
    const Vector<Double> * intPoints = ptelem->GetIntPoints();

    // **************************************************
    //  Material matrix independent of integration point
    // **************************************************
    if ( updateDMatInEveryIP_ == false ) {


      // // Check if material has to be rotated
      if( ptMaterial->GetCoordSys() == NULL )
        calcDMat(dMat, ent1.GetElem(), direction);
      /* shall be obsolete, may be deleted - Fabian
      {
        if(direction == DesignElement::NO_DERIVATIVE){
          calcDMat(dMat, ent1.GetElem());
        }else{
          calcDMat(dMat, ent1.GetElem(), direction);
        }
      }*/

      Vector<Double> globIntPoint;

      // Loop over all integration points
      for ( UInt actIntPt = 1; actIntPt <= nrIntPts; actIntPt++ ) {

        // Check if material has to be rotated for each integration point
        if( ptMaterial->GetCoordSys() != NULL ) {
          // Get global coordinates
          Vector<Double> * intPoints = ptelem->GetIntPoints();

          ptelem->Local2GlobalCoord(globIntPoint, intPoints[actIntPt-1],
                                    ptCoord_, ent1.GetElem() );
          ptMaterial->RotateTensorByPointCoord( globIntPoint,getDMaterialType() );
          if(direction == DesignElement::NO_DERIVATIVE){
            calcDMat(dMat, ent1.GetElem());
          }else{
            calcDMat(dMat, ent1.GetElem(), direction);
          }
        }


        // Setup the B matrix for current integration point
        CalcBMat(bMat, actIntPt, ptCoord_ );

        // Compute Jacobian for integration point
        jacDet = ptelem->CalcJacobianDetAtIp( actIntPt, ptCoord_, ent1.GetElem() );

        // Perform a safety check
        if ( jacDet < 0.0 ) {
          EXCEPTION( "BDBInt::CalcElementMatrix: Encountered "
                   << "negative Jacobian determinant!" );
        }

        // Special things must be done in the axi-symmetric case
        if ( isaxi_ ) {
          Vector<Double> CoordAtIP;
          ptelem->Local2GlobalCoord( CoordAtIP, intPoints[actIntPt-1],
                                   ptCoord_, ent1.GetElem() );
          //ptelem->GetShFncAtIp( ShpFncAtIp, actIntPt, ent1.GetElem() );
          jacDet *= 2 * PI * CoordAtIP[0];
        }

        // Compute the matrix product D * B and store as intermediate matrix
        dMat.Mult( bMat, dbMat );

        // We now compute B^T * D * B and scale it by the determinant
        // of the Jacobian and the weight of the current integration
        // point. The result is added to the element matrix.
        fac = jacDet * intWeights[actIntPt-1];
        for ( UInt k = 0; k < bRows; ++k ) {
          ptr1 =  bMat[k];
          ptr2 = dbMat[k];
          for ( UInt i = 0; i < bCols; ++i ) {
            ptr3 = elemMat[i];
            aux = fac * ptr1[i];
            for ( UInt j = 0; j < bCols; ++j ) {
              ptr3[j] += aux * ptr2[j];
            }
          }
        }
      } // end of loop over integration points

    }

    // **********************************************
    //  Material matrix depends on integration point
    // **********************************************
    else {

      // Loop over all integration points
      for ( UInt actIntPt = 1; actIntPt <= nrIntPts; ++actIntPt ) {

        // Setup material matrix for current integration point
        calcDMat(dMat, actIntPt, ptCoord_);

        // Setup the B matrix for current integration point
        CalcBMat( bMat, actIntPt, ptCoord_ );

        // Compute Jacobian for integration point
        jacDet = ptelem->CalcJacobianDetAtIp( actIntPt, ptCoord_, ent1.GetElem() );

        // Perform a safety check
        if ( jacDet < 0.0 ) {
          EXCEPTION( "BDBInt::CalcElementMatrix: Encountered "
                   << "negative Jacobian determinant!" );
        }

        // Special things must be done in the axi-symmetric case
        if ( isaxi_ ) {
          Vector<Double> CoordAtIP;
          ptelem->Local2GlobalCoord( CoordAtIP, intPoints[actIntPt-1],
                                     ptCoord_, ent1.GetElem() );
          //ptelem->GetShFncAtIp( ShpFncAtIp, actIntPt, ent1.GetElem() );
          jacDet *= 2 * PI * CoordAtIP[0];
        }
        
        // Compute the matrix product D * B and store as intermediate matrix
        dbMat.Resize( dMat.GetNumRows(), bMat.GetNumCols() );
        dMat.Mult( bMat, dbMat );

        // We now compute B^T * D * B and scale it by the determinant
        // of the Jacobian and the weight of the current integration
        // point. The result is added to the element matrix.
        fac = jacDet * intWeights[actIntPt-1];
        bRows = bMat.GetNumRows();
        bCols = bMat.GetNumCols();
        for ( UInt k = 0; k < bRows; ++k ) {
          ptr1 =  bMat[k];
          ptr2 = dbMat[k];
          for ( UInt i = 0; i < bCols; ++i ) {
            ptr3 = elemMat[i];
            aux = fac * ptr1[i];
            for ( UInt j = 0; j < bCols; ++j ) {
              ptr3[j] += aux * ptr2[j];
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



  // ****************************
  //   CalcComplexElementMatrix
  // ****************************
  void BDBInt::CalcComplexElementMatrix( Matrix<Complex> & elemMat,
                                          EntityIterator& ent1,
                                          EntityIterator& ent2,
                                          Double & beta, Double & omega) {


    // Get pointer to reference element and its coordinates
    ExtractElemInfo( ent1 );

    const UInt nrIntPts = ptelem->GetNumIntPoints();
    const UInt nrNodes  = ptelem->GetNumNodes();
    const UInt nrDofs   = getNrDofs();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();
    const Vector<Double> * intPoints = ptelem->GetIntPoints();
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

      CalcBMat(bMat, actIntPt, ptCoord_);

      //   hardcoded dB = dMat * bMat;
      dB.Resize(dMat.GetNumRows(), bMat.GetNumCols());
      Complex a;

			{
      const unsigned int brows(bMat.GetNumRows());
      const unsigned int drows(dMat.GetNumRows());
      const unsigned int bcols(bMat.GetNumCols());
      for ( UInt i = 0; i < drows; i++ ) {
        for ( UInt j = 0; j < bcols; j++ ) {
          a = dMat[i][0] * bMat[0][j];
          for ( UInt k = 1; k < brows; k++ ) {
            a += dMat[i][k] * bMat[k][j];
          }
          dB[i][j] = a;
        }
      }
			}

      bMat.Transpose(bTrans);

      // hardcoded: partElemMat = bTrans * dB;
      partElemMat.Resize(bTrans.GetNumRows(), dB.GetNumCols());
			{
      const unsigned int bTrows(bTrans.GetNumRows());
      const unsigned int dBcols(dMat.GetNumCols());
      const unsigned int dBrows(dB.GetNumRows());
      for ( UInt i = 0; i < bTrows; i++ ) {
        for ( UInt j = 0; j < dBcols; j++ ) {
          a = bTrans[i][0] *dB[0][j];
          for ( UInt k = 1; k < dBrows; k++ ) {
            a += bTrans[i][k] * dB[k][j];
          }
          partElemMat[i][j] = a;
        }
      }
			}

      jacDet = ptelem->CalcJacobianDetAtIp(actIntPt,ptCoord_, ent1.GetElem());

      // Perform a safety check
      if ( jacDet < 0.0 ) {
        EXCEPTION( "BDBInt::CalcElementMatrix: Encountered "
                 << "negative Jacobian determinant!" );
      }

      if ( isaxi_ ) {
        Vector<Double> CoordAtIP;
        ptelem->Local2GlobalCoord( CoordAtIP, intPoints[actIntPt-1],
                                   ptCoord_, ent1.GetElem() );
        //ptelem->GetShFncAtIp( ShpFncAtIp, actIntPt, ent1.GetElem() );
        jacDet *= 2 * PI * CoordAtIP[0];
      }
      
      const unsigned int ecols(elemMat.GetNumCols());
      const unsigned int erows(elemMat.GetNumRows());
      for ( UInt i = 0; i < erows; i++ ) {
        for ( UInt j = 0; j < ecols; j++ ) {
          elemMat[i][j] += partElemMat[i][j] * jacDet *
            intWeights[actIntPt-1] ;
        }
      }
    }
  }





  void BDBInt::calcBMat(EntityIterator it, Matrix<Double>& bMat ) {

    // get midpoint
    ExtractElemInfo( it );
    Vector<Double> midPoint;
    it.GetElem()->ptElem->GetCoordMidPoint(midPoint);

    // Set integration point to midpont
    SetIntPoint( midPoint);
    CalcBMat( bMat, 1, ptCoord_);

    //virtual void calcBMat(Matrix<Double> &bMat, UInt ip, const Matrix<Double> &ptCoord) {
    UnsetIntPoint();
  }

  void BDBInt::calcDBMat(EntityIterator it, Matrix<Double>& dbMat ) {

    // get midpoint
    ExtractElemInfo( it );
    Vector<Double> midPoint;
    it.GetElem()->ptElem->GetCoordMidPoint(midPoint);

    // Set integration point to midpont
    SetIntPoint( midPoint);
    Matrix<Double> temp1 , temp2;
    CalcBMat( temp1, 1, ptCoord_ );
    calcDMat( temp2, 1, ptCoord_ );
    dbMat = temp1*temp2;
    UnsetIntPoint();
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
