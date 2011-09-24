// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include "curlCurlEdgeInt.hh"
#include "Elements/HCurlElems.hh"
#include "Elements/fespace.hh"

namespace CoupledField
{

  CurlCurlEdgeInt::CurlCurlEdgeInt( BaseMaterial* matData, bool coordUpdate )
    : BaseForm( matData, FULL, coordUpdate )
  {
    name_ = "CurlCurlEdgeInt";

    isaxi_  = false;
    nrDofs_ = 1;
    logging_ = true;

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


 
  CurlCurlEdgeInt::~CurlCurlEdgeInt()
  {
  }
  
  void CurlCurlEdgeInt::CalcElementMatrix( Matrix<Double>& elemMat,
                                           EntityIterator& ent1,
                                           EntityIterator& ent2 ) {
    // Extract physical element
    const Elem* ptElem = ent1.GetElem();
    // Obtain FE element from feSpace
    shared_ptr<IntScheme> intScheme;
    FeHCurl *ptFe = 
        dynamic_cast<FeHCurl*>(ptFeSpace1_->GetFe( ent1, intScheme ));
    UInt nrFncs = ptFe->GetNumFncs();

    // Get shape map from grid
    shared_ptr<ElemShapeMap> esm = domain->GetGrid()->GetElemShapeMap( ptElem );

    // Get integration points (shortcut: from basefe instead of 
    // IntegrationScheme class)
    StdVector<LocPoint> intPoints;
    StdVector<Double> weights;
    intScheme->GetIntPoints( Elem::GetShapeType(ptElem->type), intPoints, weights );
    
    elemMat.Resize( nrFncs );
    elemMat.Init();

    
#define USE_BLAS_VERSION      
#ifdef USE_BLAS_VERSION
//    Matrix<Double> bdbMat;
//    bdbMat.Resize(nrFncs);
#endif
    
    // Loop over all integration points
    LocPointMapped lp;
    Matrix<Double> bMat;
    Double fac = 0.0;
    for( UInt i = 0; i < intPoints.GetSize(); i++  ) {

      // Calculate for each integration point the LocPointMapped
      lp.Set( intPoints[i], esm );

      // Call the CalcBMat()-method
      CalcBMat( bMat, lp, ptFe);

      fac = lp.jacDet * weights[i] * matVal_;
      
#ifdef USE_BLAS_VERSION
      bMat.Mult_Blas(bMat,elemMat,true,false,fac,1);
#else 
      // Compute the matrix product D * B and store as intermediate matrix
      elemMat += Transpose(bMat) * bMat * fac;
#endif
    }
  }
  
  void CurlCurlEdgeInt::CalcBMat( Matrix<Double>& bMat, 
                                  LocPointMapped& lp, FeHCurl* ptFe ) {
   
    
    ptFe->GetCurlShFnc(bMat, lp, lp.shapeMap->GetElem(), 1);
  }
  
  void CurlCurlEdgeInt::ApplyBMat( Vector<Double>& retVec,  
                                   LocPointMapped& lpm, BaseFE* ptFe, 
                                   const Vector<Double>& solVec ) {
    FeHCurl *fe = dynamic_cast<FeHCurl*>(ptFe);
    Matrix<Double> bMat;
    CalcBMat( bMat, lpm, fe);
    retVec = bMat * solVec;
  }
  void CurlCurlEdgeInt::ApplyBMat( Vector<Complex>& retVec,  
                                     LocPointMapped& lpm, BaseFE* ptFe, 
                                     const Vector<Complex>& solVec ) {
      
    FeHCurl *fe = dynamic_cast<FeHCurl*>(ptFe);
      Matrix<Double> bMat;
      CalcBMat( bMat, lpm, fe);
      retVec = bMat * solVec;
    }


//  void CurlCurlEdgeInt::CalcElementMatrix( Matrix<Double>& elemMat,
//                                            EntityIterator& ent1, 
//                                            EntityIterator& ent2 ) {
//
//     // Extract pointer to reference element and get coordinates
//     ExtractElemInfo( ent1 );
//   
//     const UInt nrIntPts= ptelem->GetNumIntPoints();
//     const UInt nrEdges = ptelem->GetNumEdges();
//     const Vector<Double> & intWeights = ptelem->GetIntWeights();  
//     Double jacDet, aux1, fac, *ptr1, *ptr2;
//   
//   
//     // derivation of shape functions after global coordinates 
//     StdVector< Matrix<Double> > xiDx;
//     xiDx.Resize(nrEdges);
//   
//     Matrix<Double> curl;
//   
//     // set matrix to desired size and set all elements to zero
//     elemMat.Resize(nrEdges); 
//     elemMat.Init();
//   
//     const Elem* geoElem = ent1.GetElem();
//
//     for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
//       {
//         // calc glob derivs of shape functions and jacobian determinante
//         ptelem->GetEdgeGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord_,
//                                           geoElem);
//       
//         CalcEdgeCurl(curl, xiDx);
//
//         jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, 
//                                              ent1.GetElem());
//         fac = jacDet * intWeights[actIntPt-1] * matVal_;
//         // We now compute B^T * D * B and scale it by the determinant
//         // of the Jacobian and the weight of the current integration
//         // point. The result is added to the element matrix.
//
//         for ( UInt k = 0; k < curl.GetNumRows(); k++ ) {
//           if ( isOrthotropic_ ) 
//             fac =  jacDet * intWeights[actIntPt-1] * reluctivityVec_[k];
//
//           ptr1 = curl[k];
//           ptr2 = curl[k];
//           for ( UInt i = 0; i < curl.GetNumCols(); i++ ) {
//             aux1 = fac * ptr1[i];
//             for ( UInt j = 0; j < curl.GetNumCols(); j++ ) {
//               elemMat[i][j] += aux1 * ptr2[j];
//             }
//           }
//         }
//       }
//   }


//   // calculates the curl; the global derivates are already given in shapeDeriv
//   void CurlCurlEdgeInt::CalcEdgeCurl(Matrix<Double>& curl, 
//                                      const StdVector<Matrix<Double> >& shapeDeriv)
//   {
//     UInt nrEdges = shapeDeriv.GetSize();
//     UInt dim = shapeDeriv[0].GetNumRows();
//     
//     curl.Resize(dim, nrEdges);
//
//     for (UInt actEdge=0; actEdge < nrEdges; actEdge++) {
//       for (UInt actDim=0; actDim < dim; actDim++) {
//         curl[actDim][actEdge] = 
//           (shapeDeriv[actEdge])[(actDim+2)%dim][(actDim+1)%dim] -
//           (shapeDeriv[actEdge])[(actDim+1)%dim][(actDim+2)%dim];
//       }
//     }
//     
//   }
} // end namespace CoupledField
