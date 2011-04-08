// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include "curlCurlNodeInt.hh"

namespace CoupledField
{

  CurlCurlNode2DInt::CurlCurlNode2DInt( BaseMaterial* matData, bool axi,
                                        bool coordUpdate )
    : BaseForm( matData, FULL, coordUpdate )
  {

    name_ = "CurlCurlNode2DInt";
    isaxi_ = axi;
    if ( matData != NULL )
      ptMaterial->GetScalar( matVal_, MAG_RELUCTIVITY, Global::REAL);
  }


 
  CurlCurlNode2DInt::~CurlCurlNode2DInt()
  {
  }



  void CurlCurlNode2DInt::CalcElementMatrix( Matrix<Double>& elemMat,
                                             EntityIterator& ent1, 
                                             EntityIterator& ent2  ) {
  
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );

    ptelem->SetAnsatzFct( ansatzFct1_ );
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    

    const UInt nrIntPts= ptelem->GetNumIntPoints();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Double jacDet;  


    // derivation of shape functions after global coordinates 
    Matrix<Double> xiDx;
    Matrix<Double> xiDxTransp;
    Matrix<Double> partElemMat;
    Matrix<Double> partElemMatAxi;
    Vector<Double> ShpFncAtIp;
    Vector<Double> CoordAtIP;
    Vector<Double> drAtIp;


    // set matrix to desired size and set all elements to zero
    elemMat.Resize(numFncs); 
    elemMat.Init();
    
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {
      jacDet = 0;
      
      ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord_, 
                                    jacDet, ent1.GetElem() );
      
      if (isaxi_) {
        ptelem->GetShFncAtIp(ShpFncAtIp,actIntPt, ent1.GetElem() );
        CoordAtIP = ptCoord_ * ShpFncAtIp;
        for (UInt i=0; i<numFncs; i++)
          xiDx[i][0] += ShpFncAtIp[i] / CoordAtIP[0];
        
        jacDet *= 2 * PI * CoordAtIP[0];
      }
      
      xiDx.Transpose(xiDxTransp);
      partElemMat = xiDx * xiDxTransp;
      partElemMat *= intWeights[actIntPt-1] * jacDet * matVal_;
      elemMat += partElemMat;
    }
  
  }

  
  void CurlCurlNode2DInt::CalcBMat( Matrix<Double> &bMat, UInt ip, const Matrix<Double> &ptCoord )
  {
    Matrix<Double> xiDx;
    if (isSetIntPoint_) {
      ptelem->GetGlobDerivShFnc(xiDx, intPoint_, ptCoord, 
                                it1_.GetElem() );
    } else { 
      ptelem->GetGlobDerivShFncAtIp(xiDx, ip, ptCoord, 
                                    it1_.GetElem() );
    }
    if (isaxi_) {
      Vector<Double> ShpFncAtIp;
      Vector<Double> CoordAtIP;
      if( isSetIntPoint_ ) {
        ptelem->GetShFnc(ShpFncAtIp,intPoint_,it1_.GetElem());
      } else {
        ptelem->GetShFncAtIp(ShpFncAtIp,ip,it1_.GetElem());
      }
      CoordAtIP = ptCoord * ShpFncAtIp;
      for (UInt i=0; i<ShpFncAtIp.GetSize(); i++)
        xiDx[i][0] += ShpFncAtIp[i] / CoordAtIP[0];
    }
    xiDx.Transpose(bMat);

  }
  
  //============================Curl-Curl-3D ====================================
           
  CurlCurlNode3DInt::CurlCurlNode3DInt( BaseMaterial* matData, bool coordUpdate )
    : BaseForm( matData, FULL, coordUpdate )

  {

    name_   = "CurlCurlNode3DInt";
    isaxi_  = false;
    nrDofs_ = 3;

    isOrthotropic_ = false;
    if ( matData != NULL ) {
      if ( matData->GetSymmetryType() == BaseMaterial::ORTHOTROPIC ) {
        isOrthotropic_ = true;
        reluctivityVec_.Resize(3);
        ptMaterial->GetScalar( matVal_, MAG_PERMEABILITY_1, Global::REAL);
        reluctivityVec_[0] = 1.0 / matVal_;
        ptMaterial->GetScalar( matVal_, MAG_PERMEABILITY_2, Global::REAL);
        reluctivityVec_[1] = 1.0 / matVal_;
        ptMaterial->GetScalar( matVal_, MAG_PERMEABILITY_3, Global::REAL);
        reluctivityVec_[2] = 1.0 / matVal_;
        //        std::cout << "Orthotropic: \n" << reluctivityVec_ << std::endl;
      }
      else {
        ptMaterial->GetScalar( matVal_, MAG_RELUCTIVITY, Global::REAL);
        //std::cout << "Isotropic: mu=" << matVal_ << std::endl;
      }
    }
  }


 
  CurlCurlNode3DInt::~CurlCurlNode3DInt()
  {
  }



  void CurlCurlNode3DInt::CalcElementMatrix( Matrix<Double>& elemMat,
                                             EntityIterator& ent1, 
                                             EntityIterator& ent2  ) {
  
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );

    ptelem->SetAnsatzFct( ansatzFct1_ );
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    
    Double jacDet;  

    Matrix<Double> bMatCurl, bMatDiv;
    Double aux1, aux2, fac, *ptr1, *ptr2, *ptr3, *ptr4;

    elemMat.Resize( numFncs * nrDofs_ );
    elemMat.Init();

    //get integration points
    const UInt nrIntPts = ptelem->GetNumIntPoints(); 
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  

    // Loop over all integration points
    for ( UInt actIntPt = 1; actIntPt <= nrIntPts; actIntPt++ ) {

      // Setup the B matrix for current integration point
      calcBMat( bMatCurl, bMatDiv, actIntPt, ptCoord_ );

      // Compute Jacobian for integration point
      jacDet = ptelem->CalcJacobianDetAtIp( actIntPt, ptCoord_, ent1.GetElem() );

      // Perform a safety check
      if ( jacDet < 0.0 ) {
        EXCEPTION( "CurlCurlNode3DInt::CalcElementMatrix: Encountered "
            << "negative Jacobian determinant!");
      }

      // We now compute B^T * D * B and scale it by the determinant
      // of the Jacobian and the weight of the current integration
      // point. The result is added to the element matrix.
      fac = jacDet * intWeights[actIntPt-1] * matVal_;

//       bMatCurl.Transpose(bMatCurlT);
//       partElemMat1 = bMatCurlT * bMatCurl;

//       bMatDiv.Transpose(bMatDivT);
//       partElemMat2 = bMatDivT * bMatDiv;

//       partMat += partElemMat1;
//       //      partMat += partElemMat2;
//       partMat *= fac;
 
//       elemMat += partMat;

      for ( UInt k = 0; k < bMatCurl.GetNumRows(); k++ ) {
        if ( isOrthotropic_ ) 
          fac =  jacDet * intWeights[actIntPt-1] * reluctivityVec_[k];

        ptr1 = bMatCurl[k];
        ptr2 = bMatCurl[k];
        ptr3 = bMatDiv[k];
        ptr4 = bMatDiv[k];
        for ( UInt i = 0; i < bMatCurl.GetNumCols(); i++ ) {
          aux1 = ptr1[i];
          aux2 = ptr3[i];
          for ( UInt j = 0; j < bMatCurl.GetNumCols(); j++ ) {
            elemMat[i][j] += ( aux1 * ptr2[j] + aux2 * ptr4[j] ) * fac;
          }
        }
      }
    }
    //std::cout << "StiffMat:\n" << elemMat << std::endl;
  }


  // returns curl and div - matrix
  void CurlCurlNode3DInt::calcBMat( Matrix<Double> &bMatCurl, Matrix<Double> &bMatDiv,
				    UInt ip, Matrix<Double> &ptCoord ) {


    const UInt numFncs  = ptelem->GetNumFncs( ansatzFct1_ );

    UInt actNode;
        
    bMatCurl.Resize( nrDofs_, numFncs * nrDofs_);
    bMatCurl.Init();

    bMatDiv.Resize( nrDofs_, numFncs * nrDofs_);
    bMatDiv.Init();

    // local shape functions derived after global coords
    // (format: numFncs x spaceDim)
    Matrix<Double> xiDx;
    ptelem->SetAnsatzFct( ansatzFct1_ );
      
    if (isSetIntPoint_) 
      ptelem->GetGlobDerivShFnc(xiDx, intPoint_, ptCoord, it1_.GetElem() );
    else {
      //std::cout << "ip : " << ip << std::endl;
      ptelem->GetGlobDerivShFncAtIp(xiDx, ip, ptCoord, it1_.GetElem() );
    }

    for(actNode=0; actNode < numFncs; actNode++) {
      //see M. Kaltenbacher 1.st edition, pp. 92
      //first row
      bMatCurl[0][actNode * nrDofs_ + 1] = -xiDx[actNode][2];
      bMatCurl[0][actNode * nrDofs_ + 2] =  xiDx[actNode][1];

      //second row
      bMatCurl[1][actNode * nrDofs_]     =  xiDx[actNode][2];
      bMatCurl[1][actNode * nrDofs_ + 2] = -xiDx[actNode][0];

      //third row
      bMatCurl[2][actNode * nrDofs_]     = -xiDx[actNode][1];
      bMatCurl[2][actNode * nrDofs_ + 1] =  xiDx[actNode][0];

//       bMatCurl[0][actNode * nrDofs_]       +=  xiDx[actNode][0];
//       bMatCurl[1][actNode * nrDofs_ + 1]   +=  xiDx[actNode][1];
//       bMatCurl[2][actNode * nrDofs_ + 2]   +=  xiDx[actNode][2];

      //first row
      bMatDiv[0][actNode * nrDofs_]      = xiDx[actNode][0];

      //second row
      bMatDiv[1][actNode * nrDofs_ + 1]  = xiDx[actNode][1];

      //third row
      bMatDiv[2][actNode * nrDofs_ + 2]  = xiDx[actNode][2];

    }
    //bMatDiv.Init();
//     std::cout << "bCurl:\n" << bMatCurl << std::endl;
//     std::cout << "bDiv:\n" << bMatDiv << std::endl;


    isSetIntPoint_ = false;
  }


  //============================coupling of vector and scalar potential ====================================

  MagCoupVectorScalarPotentialInt::MagCoupVectorScalarPotentialInt(Double aVal, bool coordUpdate )
    : BaseForm(NULL,FULL,coordUpdate ),matVal_ (aVal)
  {

    name_   = "MagCoupVectorScalarPotentialInt";
    isaxi_  = false;
    nrDofsVec_ = 3;

    // this bilinearform is never symmetric
    isSymmetric_ = false;
  }


 
  MagCoupVectorScalarPotentialInt::~MagCoupVectorScalarPotentialInt()
  {
  }



  void MagCoupVectorScalarPotentialInt::CalcElementMatrix( Matrix<Double>& elemMat,
                                             EntityIterator& ent1, 
                                             EntityIterator& ent2  ) {
  
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );

    ptelem->SetAnsatzFct( ansatzFct1_ );
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    elemMat.Resize( numFncs * nrDofsVec_, numFncs );
    elemMat.Init();

    Matrix<Double> xiDx;
    Vector<Double> ShpFncAtIp;
    Double jacDet, factor, shpFnc;  
    
    //get integration points
    const UInt nrIntPts = ptelem->GetNumIntPoints(); 
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  

    // Loop over all integration points
    for ( UInt actIntPt = 1; actIntPt <= nrIntPts; actIntPt++ ) {

      jacDet = 0;
      ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord_, 
                                    jacDet, ent1.GetElem());
      ptelem->GetShFncAtIp(ShpFncAtIp,actIntPt,ent1.GetElem());

      // Perform a safety check
      if ( jacDet < 0.0 ) {
        EXCEPTION("MagCoupVectorScalarPotentialInt::CalcElementMatrix: Encountered "
            << "negative Jacobian determinant!");
      }

      factor = jacDet * intWeights[actIntPt-1] * matVal_;
      for ( UInt k=0; k<numFncs; k++ ) {
        shpFnc = ShpFncAtIp[k] * factor;
        for ( UInt  i= 0; i < nrDofsVec_; i++ ) {
          for ( UInt j = 0; j < numFncs; j++ ) {
            elemMat[i+k*nrDofsVec_][j] += shpFnc * xiDx[j][i];
          }
        }
      }
    }

  }


} // end namespace CoupledField
