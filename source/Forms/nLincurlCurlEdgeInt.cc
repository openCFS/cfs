// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include "Utils/nodestoresol.hh"
#include "nLincurlCurlEdgeInt.hh"

namespace CoupledField
{

  nLinCurlCurlEdgeInt::nLinCurlCurlEdgeInt( BaseMaterial* matData, bool coordUpdate )
    : CurlCurlEdgeInt( matData, coordUpdate )
  {
    name_ = "nLinCurlCurlEdgeInt";

    isSolDependent_ = true;
    nonLinType_ = NEWTON;

    // get pointer to nonlinear BH curve approximation
    ptMaterial->NeedApproxMatCurve( MAG_PERMEABILITY );
    
    // fetch real start values for BH-curve
    if( isOrthotropic_ ) {
      EXCEPTION("Not implemented")
//      reluctivityVec_.Resize(3);
//      ptMaterial->GetScalar( matVal_, MAG_PERMEABILITY_START_1, REAL);
//      reluctivityVec_[0] = 1.0 / matVal_;
//      ptMaterial->GetScalar( matVal_, MAG_PERMEABILITY_START_2, REAL);
//      reluctivityVec_[1] = 1.0 / matVal_;
//      ptMaterial->GetScalar( matVal_, MAG_PERMEABILITY_START_3, REAL);
//      reluctivityVec_[2] = 1.0 / matVal_;
//
//      // resize for some additional arrays
//      currReluctivityVec_.Resize(3);
//      currDerivReluctivityVec_.Resize(3);
    } else {
      //ptMaterial->GetScalar( matVal_, MAG_PERMEABILITY_START, Global::REAL);
      //matVal_ = 1.0 / matVal_;
      ptMaterial->GetScalar( matVal_, MAG_RELUCTIVITY,Global::REAL);
    }
  }


 
  nLinCurlCurlEdgeInt::~nLinCurlCurlEdgeInt()
  {
  }



  void nLinCurlCurlEdgeInt::CalcElementMatrix( Matrix<Double>& elemMat,
                                           EntityIterator& ent1, 
                                           EntityIterator& ent2 ) {

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );
    const UInt nrIntPts= ptelem->GetNumIntPoints();
    const UInt nrEdges = ptelem->GetNumEdges();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  

    Double jacDet = 0.0;
    Double reluctivity = 0.0;
    Double derivReluctivity = 0.0;  
  
    // get pointer to nonlinear BH curve approximation
    ApproxData* nlinFnc = ptMaterial->GetNonlinFnc(MAG_PERMEABILITY);

 
    // get solution of current element
    sol_->GetElemSolution( magPot_, ent1 );

    // derivation of shape functions after global coordinates 
    StdVector< Matrix<Double> > xiDx;
    xiDx.Resize(nrEdges);
  
    Matrix<Double> curl;
    Vector<Double> elemFlux(3);
    Vector<Double> help( nrEdges );
    Vector<Double> helpAI( nrEdges );
    Double aux1, fac, *ptr1, *ptr2;

    // set matrix to desired size and set all elements to zero
    elemMat.Resize(nrEdges); 
    elemMat.Init();
  
    const Elem* geoElem = ent1.GetElem();

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {
      // calc glob derivs of shape functions and jacobian determinant
      ptelem->GetEdgeGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord_,
                                        geoElem);
      
      CalcEdgeCurl(curl, xiDx);
         
      //compute magnetic flux density
      elemFlux = curl * magPot_;
      Double Babs = elemFlux.NormL2();
      
      if ( Babs == 0 ) {
        if ( isOrthotropic_ ) 
          currReluctivityVec_ = reluctivityVec_;
        else 
          reluctivity = matVal_;
      }
      else {
        if ( isOrthotropic_ ) {
//          for ( UInt i=0; i<3; i++ ) {
//            currReluctivityVec_[i] = 
//              nlinFnc_[i]->EvaluateFuncNu( abs(elemFlux[i]) );
//          }
        }
        else {
          reluctivity = nlinFnc->EvaluateFuncNu(Babs);
        }
      }
        
      jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, 
                                           ent1.GetElem());
      
      // We now compute B^T * D * B and scale it by the determinant
      // of the Jacobian and the weight of the current integration
      // point. The result is added to the element matrix.
      fac = jacDet * intWeights[actIntPt-1] * reluctivity;
      for ( UInt k = 0; k < curl.GetNumRows(); k++ ) {
        if ( isOrthotropic_ ) 
          fac =  jacDet * intWeights[actIntPt-1] * currReluctivityVec_[k];
        
        ptr1 = curl[k];
        ptr2 = curl[k];
        for ( UInt i = 0; i < curl.GetNumCols(); i++ ) {
          aux1 = fac * ptr1[i];
          for ( UInt j = 0; j < curl.GetNumCols(); j++ ) {
            elemMat[i][j] += aux1 * ptr2[j];
          }
        }
      }
//      std::cerr << "\n\n-------------------------------\n";
//      std::cerr << "matrix before:\n" << elemMat << std::endl; 
      if ( nonLinType_ == NEWTON ) {
        if ( Babs == 0.0 ) 
          derivReluctivity = 0;
        else {          
          //Newton method
          if ( isOrthotropic_ ) {
            //fac = jacDet * intWeights[actIntPt-1] * Babs;
            //          if ( isOrthotropic_ ) {
            //            for ( UInt i=0; i<3; i++ ) {
            //              currDerivReluctivityVec_[i] = 
            //                nlinFnc_[i]->EvaluatePrimeNu( abs(elemFlux[i]) );
            //           }
          } else {
            derivReluctivity =  nlinFnc->EvaluatePrimeNu(Babs);
            fac = jacDet * intWeights[actIntPt-1] * derivReluctivity * Babs;
          }

          Vector<Double> eB(3); eB = elemFlux /Babs;
          for ( UInt k = 0; k < curl.GetNumCols(); k++ ) 
            for ( UInt i = 0; i < curl.GetNumRows(); i++ ) 
              help[k] =  curl[i][k] * eB[i];
          
          if ( isOrthotropic_ ) {
            for ( UInt k = 0; k < curl.GetNumCols(); k++ ) 
              for ( UInt i = 0; i < curl.GetNumRows(); i++ ) 
                helpAI[k] =  curl[i][k] * eB[i] 
                  * currDerivReluctivityVec_[i];
          }
          
          if ( isOrthotropic_ ) {
            for ( UInt i = 0; i< curl.GetNumCols(); i++ ) 
              for ( UInt j = 0; j< curl.GetNumCols(); j++ ) 
                elemMat[i][j] += fac * helpAI[i] * help[j];
          }
          else {
            for ( UInt i = 0; i< curl.GetNumCols(); i++ ) 
              for ( UInt j = 0; j< curl.GetNumCols(); j++ ) 
                elemMat[i][j] += fac * help[i] * help[j];
          }
        }
      //std::cerr << "matrix after:\n" << elemMat << "\n\n";
      }
    }

#ifndef NDEBUG 
    (*OutInfo::debug) << "nLinCurCurlEdge Matrix:  "  << std::endl
             << elemMat << std::endl << std::endl;
#endif
  
  }

  void  nLinCurlCurlEdgeInt::SetNonLinMethod(NonLinMethodType atype) {
    nonLinType_ = atype;
  }

} // end namespace CoupledField
