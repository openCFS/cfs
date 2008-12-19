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
     Double jacDet, aux1, fac, *ptr1, *ptr2;
   
   
     // derivation of shape functions after global coordinates 
     StdVector< Matrix<Double> > xiDx;
     xiDx.Resize(nrEdges);
   
     Matrix<Double> curl;
   
     // set matrix to desired size and set all elements to zero
     elemMat.Resize(nrEdges); 
     elemMat.Init();
   
     const Elem* geoElem = ent1.GetElem();

     for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
       {
         // calc glob derivs of shape functions and jacobian determinante
         ptelem->GetEdgeGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord_,
                                           geoElem);
       
         CalcEdgeCurl(curl, xiDx);

         jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, 
                                              ent1.GetElem());
         fac = jacDet * intWeights[actIntPt-1] * matVal_;
         // We now compute B^T * D * B and scale it by the determinant
         // of the Jacobian and the weight of the current integration
         // point. The result is added to the element matrix.

         for ( UInt k = 0; k < curl.GetSizeRow(); k++ ) {
           if ( isOrthotropic_ ) 
             fac =  jacDet * intWeights[actIntPt-1] * reluctivityVec_[k];

           ptr1 = curl[k];
           ptr2 = curl[k];
           for ( UInt i = 0; i < curl.GetSizeCol(); i++ ) {
             aux1 = fac * ptr1[i];
             for ( UInt j = 0; j < curl.GetSizeCol(); j++ ) {
               elemMat[i][j] += aux1 * ptr2[j];
             }
           }
         }
       }
   }


   // calculates the curl; the global derivates are already given in shapeDeriv
   void CurlCurlEdgeInt::CalcEdgeCurl(Matrix<Double>& curl, 
                                      const StdVector<Matrix<Double> >& shapeDeriv)
   {
     UInt nrEdges = shapeDeriv.GetSize();
     UInt dim = shapeDeriv[0].GetSizeRow();
     
     curl.Resize(dim, nrEdges);

     for (UInt actEdge=0; actEdge < nrEdges; actEdge++) {
       for (UInt actDim=0; actDim < dim; actDim++) {
         curl[actDim][actEdge] = 
           (shapeDeriv[actEdge])[(actDim+2)%dim][(actDim+1)%dim] -
           (shapeDeriv[actEdge])[(actDim+1)%dim][(actDim+2)%dim];
       }
     }
     
   }
} // end namespace CoupledField
