// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <string>

#include "Domain/entityList.hh"
#include "Elements/basefe.hh"
#include "Forms/curlCurlEdgeInt.hh"
#include "General/exception.hh"
#include "MatVec/exprt/xpr1.hh"
#include "MatVec/exprt/xpr2.hh"
#include "MatVec/matrix.hh"
#include "Materials/baseMaterial.hh"
#include "Utils/ApproxData.hh"
#include "Utils/StdVector.hh"
#include "Utils/nodestoresol.hh"
#include "nLincurlCurlEdgeInt.hh"

namespace CoupledField {
struct Elem;
}  // namespace CoupledField

namespace CoupledField
{

  nLinCurlCurlEdgeInt::nLinCurlCurlEdgeInt( BaseMaterial* matData, bool coordUpdate )
    : CurlCurlEdgeInt( matData, coordUpdate )
  {
    name_ = "nLinCurlCurlEdgeInt";

    isSolDependent_ = true;
    nonLinType_ = NEWTON;

    //check for anisotropic / isotropic material
    isAnisotropic_ = false;    
    if ( ptMaterial->IsNonLinMagBHcurves() )
      isAnisotropic_ = true;

    // get pointer to nonlinear BH curve approximation
    ptMaterial->NeedApproxMatCurve( MAG_PERMEABILITY );
    
    // fetch real start values for BH-curve
    ptMaterial->GetScalar( matVal_, MAG_RELUCTIVITY,Global::REAL);

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
          reluctivity = matVal_;
      }
      else if ( isAnisotropic_ ) {
        ComputeAnisotropicReluctivityOrDeriv( reluctivity, 
                                              elemFlux, true);
      }
      else {
        reluctivity = nlinFnc->EvaluateFuncNu(Babs);
      }
              
      jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, 
                                           ent1.GetElem());
      
      // We now compute B^T * D * B and scale it by the determinant
      // of the Jacobian and the weight of the current integration
      // point. The result is added to the element matrix.
      fac = jacDet * intWeights[actIntPt-1] * reluctivity;
      for ( UInt k = 0; k < curl.GetNumRows(); k++ ) {
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
//      std::cerr << "matrix before:\n" << elemMat << reluctivity = std::endl; 
      if ( nonLinType_ == NEWTON ) {
        if ( Babs == 0.0 ) {
          derivReluctivity = 0;
        }
        else {
          if ( isAnisotropic_ ) {
            ComputeAnisotropicReluctivityOrDeriv( derivReluctivity,
                                                  elemFlux, false);
          }
          else {
            //Newton method
            derivReluctivity =  nlinFnc->EvaluatePrimeNu(Babs);
          }

          fac = jacDet * intWeights[actIntPt-1] * derivReluctivity * Babs;
          
          Vector<Double> eB(3); eB = elemFlux * (1/Babs);

          Matrix<Double> dMat;
          dMat.DyadicMult( eB, eB );
          dMat *= fac;

          Matrix<Double> dbMat;
          dbMat = dMat * curl;

          for ( UInt k = 0; k < curl.GetNumRows(); k++ ) {
            ptr1 = curl[k];
            ptr2 = dbMat[k];
            for ( UInt i = 0; i < curl.GetNumCols(); i++ ) {
              aux1 = ptr1[i];
              for ( UInt j = 0; j < dbMat.GetNumCols(); j++ ) {
                elemMat[i][j] += aux1 * ptr2[j];
              }
            }
          }

// ==============  old implementation: wrong!! =================
//           for ( UInt k = 0; k < curl.GetNumCols(); k++ ) 
//             for ( UInt i = 0; i < curl.GetNumRows(); i++ ) 
//               help[k] =  curl[i][k] * eB[i];
          
//           for ( UInt i = 0; i< curl.GetNumCols(); i++ ) 
//             for ( UInt j = 0; j< curl.GetNumCols(); j++ ) 
//               elemMat[i][j] += fac * help[i] * help[j];

        }
      }

    }

    //#ifndef NDEBUG 
    //    (*OutInfo::debug) << "nLinCurCurlEdge Matrix:  "  << std::endl
    //             << elemMat << std::endl << std::endl;
    //#endif
  
  }

  void  nLinCurlCurlEdgeInt::SetNonLinMethod(NonLinMethodType atype) {
    nonLinType_ = atype;
  }


  void nLinCurlCurlEdgeInt::ComputeAnisotropicReluctivityOrDeriv( Double& val,
                                                                  Vector<double>& vecB,
                                                                  bool isReluctivity ) {

    //get pointers to the nlFunctions
    StdVector<ApproxData*> nlinFncs = 
     ptMaterial->GetNonlinFncs(  MAG_PERMEABILITYCURVES );

    //compute angle 
    Double angleB;
    if ( abs(vecB[0]) > 1e-5 ) {
      angleB = abs( std::atan( vecB[1] / vecB[0] ) );
      angleB *= 180.0/3.1416;
    }
    else {
      angleB = 90.0;
    }

    StdVector<Double>& anglesCurve = ptMaterial->GetAnisotropicAngles();
    UInt pos = 0;
    Double dist, minDist;
    minDist = abs( anglesCurve[0] - angleB );
    for (UInt i=1; i<anglesCurve.GetSize(); i++ ) {
      dist = abs( anglesCurve[i] - angleB );
      if ( dist < minDist ) {
        pos = i;
        minDist = dist;
      }
    }

    Double absB = vecB.NormL2();
    if ( isReluctivity ) 
      val = nlinFncs[pos]->EvaluateFuncNu( absB );
    else
      val = nlinFncs[pos]->EvaluatePrimeNu( absB );

  }

} // end namespace CoupledField
