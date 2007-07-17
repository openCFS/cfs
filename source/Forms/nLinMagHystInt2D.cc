// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include "nLinMagHystInt2D.hh"
#include "Utils/nodestoresol.hh"
#include "Forms/curlCurlNodeInt.hh"

namespace CoupledField
{

  void nLinMagHystInt2D::CalcElementMatrix( Matrix<Double>& elemMat,
                                          EntityIterator& ent1, 
                                          EntityIterator& ent2 ) {

    ENTER_FCN( "nLinMagHystInt2D::CalcElementMatrix" );

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );


    // First of all, set ansatz function to element
    ptelem->SetAnsatzFct( ansatzFct1_ );

    const UInt nrFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const UInt spaceDim = ptelem->GetDim();  
    const UInt elNr = ent1.GetElem()->elemNum;

    Double jacDet;
    Matrix<Double> bMat; 
    Matrix<Double> dMat; 
    Matrix<Double> dbMat; 
    Double aux, fac, *ptr1, *ptr2;

    elemMat.Resize( nrFncs);
    elemMat.Init();

    dMat.Resize( spaceDim, spaceDim );
    dbMat.Resize( spaceDim, nrFncs );

    //get integration points   
    const UInt nrIntPts = ptelem->GetNumIntPoints(); 
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  

    // get solution of element
    Matrix<Double> temp;
    sol_->GetElemSolutionAsMatrix( temp, ent1 );
    temp.ConvertToVec_AppendRows( magPot_ );

    // Loop over all integration points
    for ( UInt actIntPt = 1; actIntPt <= nrIntPts; actIntPt++ ) {
    
      // Setup the B matrix for current integration point
      // and compute the B-field
      calcMyBMat( bMat, actIntPt, ptCoord_ );

      //compute material tensor as a function of the B-field      
      calcDMat( dMat, elNr );

      // Compute Jacobian for integration point
      jacDet = ptelem->CalcJacobianDetAtIp( actIntPt, ptCoord_, ent1.GetElem() );
      
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
        ptelem->GetShFncAtIp( ShpFncAtIp, actIntPt, ent1.GetElem() );
        
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
    std::cout << "elenMat:\n" << elemMat << std::endl;
  }
  
  
  // returns B - matrix for BDB
  void nLinMagHystInt2D::calcMyBMat( Matrix<Double> &bMat, UInt ip,
                                 Matrix<Double> &ptCoord ) {

    ENTER_FCN( "nLinMagHystInt2D::calcBMat" );

    const UInt numFncs  = ptelem->GetNumFncs( ansatzFct1_ );
    const UInt spaceDim = ptelem->GetDim();  
    const UInt nrDofs   = getNrDofs();  

    UInt actNode;
    
    
    bMat.Resize( spaceDim, numFncs * nrDofs);
    bMat.Init();

    // local shape functions derived after global coords
    // (format: numFncs x spaceDim)
    Matrix<Double> xiDx;
    ptelem->SetAnsatzFct( ansatzFct1_ );
      
    if (isSetIntPoint_) 
      ptelem->GetGlobDerivShFnc(xiDx, intPoint_, ptCoord, it1_.GetElem() );
    else
      ptelem->GetGlobDerivShFncAtIp(xiDx, ip, ptCoord, it1_.GetElem() );

    switch(spaceDim)
      {
      case 2:
        for (actNode = 0; actNode < numFncs; actNode++) {
          for ( UInt dim=0; dim < spaceDim; dim++ ) {
            bMat[dim][actNode] = xiDx[actNode][0];
            bMat[dim][actNode] = xiDx[actNode][1];
          }
        }

        if (isaxi_) {
          Error("nLinMagHystInt::calcBMat for axisymmetric case not implemented",
                __FILE__, __LINE__ );
        }
        break;


      case 3:
        Error("nLinMagHystInt::calcBMat for 3D case not implemented",
              __FILE__, __LINE__ );
        break;
      }

    // here we also compute B, which we will need for the hysteresis
    Bfield_.Resize(spaceDim);
    Bfield_.Init();
    for( UInt i=0; i<spaceDim; i++ )
      for( UInt j=0; j<numFncs; j++ )
        Bfield_[i] += xiDx[j][i] * magPot_[j];

    // Account for curl 
    Double temp = Bfield_[0];
    if ( isaxi_ ) {
      Bfield_[0] = -Bfield_[1];
      Bfield_[1] = temp;
    } else {
      Bfield_[0] = Bfield_[1];
      Bfield_[1] = -temp;
    }
  }


  void nLinMagHystInt2D::calcDMat(Matrix<Double> & dMat, UInt elNr)
  {
    ENTER_FCN( "nLinMagHystInt::calcDMat" );

    Vector<Double> scalarVals(2);

    ptMaterial->ComputeScalarDiffValues( elNr, Bfield_, scalarVals );
    std::cout << "scalarVals:\n" << scalarVals << std::endl;

    dMat[0][0] = 1 + scalarVals[0];
    dMat[1][1] = 1 + scalarVals[0];
    dMat[0][1] = scalarVals[1];
    dMat[1][0] = scalarVals[1];

    dMat *=  reluctivity0_;
    std::cout << "dMat:\n" << dMat << std::endl;
  }

  // returns B for postprocessing
  void nLinMagHystInt2D::calcBMat( Matrix<Double> &bMat, UInt ip,
                                     Matrix<Double> &ptCoord ) {
    ENTER_FCN( "nLinMagHystInt::calcBMat" );

    ptBMat_->ExtractElemInfo( it1_ );
    ptBMat_->SetIntPoint( intPoint_ );
    ptBMat_->calcBMat( bMat, ip, ptCoord );
    ptBMat_->UnsetIntPoint(); 

  }


  // =========================================================================
  // =============== standard con- and destructors (just for tracing) ========
  // =========================================================================

  // NOTE: We hardcode the orientation for 2D simulations to use the yz plane!

  nLinMagHystInt2D::nLinMagHystInt2D( BaseMaterial* matData, bool axi, bool coordUpdate) :
    BDBInt(matData) {
    ENTER_FCN( "nLinMagHystInt2D::nLinMagHystInt2D" );

    name_        = "nLinMagHystInt2D";
    isaxi_       = axi;
    coordUpdate_ = coordUpdate;

    isSolDependent_ = true;

    ptMaterial->GetScalar( startmatVal_, MAG_RELUCTIVITY,REAL);

    reluctivity0_ = 7.9577E+05;

    ptBMat_ = new  CurlCurlNode2DInt( matData, axi, coordUpdate);
  
  }


  nLinMagHystInt2D::~nLinMagHystInt2D() {
    ENTER_FCN( "nLinMagHystInt2D::~nLinMagHystInt2D" );
  }

} // end namespace CoupledField
