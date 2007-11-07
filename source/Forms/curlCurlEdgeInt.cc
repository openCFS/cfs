// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include "curlCurlEdgeInt.hh"

namespace CoupledField
{

  CurlCurlEdgeInt::CurlCurlEdgeInt( BaseMaterial* matData, bool coordUpdate )
    : BaseForm( matData, FULL, coordUpdate )
  {
    name_ = "CurlCurlEdgeInt";

    isaxi_  = false;
    nrDofs_ = 1;

    isOrthotropic_ = false;
    if ( matData != NULL ) {
      if ( matData->GetSymmetryType() == BaseMaterial::ORTHOTROPIC ) {
        isOrthotropic_ = true;
        reluctivityVec_.Resize(3);
        ptMaterial->GetScalar( matVal_, MAG_PERMEABILITY_1, REAL);
        reluctivityVec_[0] = 1.0 / matVal_;
        ptMaterial->GetScalar( matVal_, MAG_PERMEABILITY_2, REAL);
        reluctivityVec_[1] = 1.0 / matVal_;
        ptMaterial->GetScalar( matVal_, MAG_PERMEABILITY_3, REAL);
        reluctivityVec_[2] = 1.0 / matVal_;
        //        std::cout << "Orthotropic: \n" << reluctivityVec_ << std::endl;
      }
      else {
        ptMaterial->GetScalar( matVal_, MAG_RELUCTIVITY, REAL);
        //std::cout << "Isotropic: mu=" << matVal_ << std::endl;
      }
    }
  }


 
  CurlCurlEdgeInt::~CurlCurlEdgeInt()
  {
  }



  void CurlCurlEdgeInt::CalcElementMatrix( Matrix<Double>& elemMat,
                                           EntityIterator& ent1, 
                                           EntityIterator& ent2 ) {

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );
  
    const UInt nrIntPts= ptelem->GetNumIntPoints();
    const UInt nrEdges = ptelem->GetNumEdges();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Double jacDet;  
  
  
    // derivation of shape functions after global coordinates 
    StdVector< Matrix<Double>* > xiDx;
    xiDx.Resize(nrEdges);
    for (UInt i=0; i<nrEdges; i++)
      xiDx[i] = new Matrix<Double>;
  
    Matrix<Double> curl;
    Matrix<Double> curlTransp;
    Matrix<Double> partElemMat;
  
    // set matrix to desired size and set all elements to zero
    elemMat.Resize(nrEdges,true); 
    elemMat.Init();
  
  
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {
        // calc glob derivs of shape functions and jacobian determinante
        ptelem->GetEdgeGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord_);
      
        CalcEdgeCurl(curl, xiDx);
      
        curl.Transpose(curlTransp);
      
        partElemMat = curlTransp * curl;
      
        jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, ent1.GetElem());
      
        partElemMat *= intWeights[actIntPt-1] * jacDet * matVal_;
      
        elemMat += partElemMat;
      }
  
  }


  // calculates the curl, if the global derivates are already given in shapeDeriv
  void CurlCurlEdgeInt::CalcEdgeCurl(Matrix<Double>& curl, 
                                     const StdVector<Matrix<Double>*>& shapeDeriv)
  {
    UInt nrEdges = shapeDeriv.GetSize();
    UInt dim = shapeDeriv[0]->GetSizeRow();
    
    curl.Resize(dim, nrEdges);
    
    for (UInt actEdge=0; actEdge < nrEdges; actEdge++)
      for (UInt actDim=0; actDim < dim; actDim++)
        curl[actDim][actEdge] = 
          (*shapeDeriv[actEdge])[(actDim+2)%dim][(actDim+1)%dim] -
          (*shapeDeriv[actEdge])[(actDim+1)%dim][(actDim+2)%dim];
    
  }
  

} // end namespace CoupledField
