// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include "Utils/nodestoresol.hh"
#include "nLincurlCurlEdgeInt.hh"
#include "Elements/HCurlElems.hh"

namespace CoupledField
{

  nLinCurlCurlEdgeInt::nLinCurlCurlEdgeInt( BaseMaterial* matData, bool coordUpdate )
    : CurlCurlEdgeInt( matData, coordUpdate )
  {
    name_ = "nLinCurlCurlEdgeInt";
    isSolDependent_ = true;
    nonLinType_ = NEWTON;

    // get pointer to nonlinear BH curve approximation
    ptMaterial->NeedApproxMatCurve( magBH );
    
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
    
    // Initialize file streams with elements to observe
    fluxElems_ =  1043, 3152, 3656, 3396;

  }


 
  nLinCurlCurlEdgeInt::~nLinCurlCurlEdgeInt()
  {
  }



  void nLinCurlCurlEdgeInt::CalcElementMatrix( Matrix<Double>& elemMat,
                                           EntityIterator& ent1, 
                                           EntityIterator& ent2 ) {
    // Extract physical element
    const Elem* ptElem = ent1.GetElem();
    FeHCurl* ptFe = dynamic_cast<FeHCurl*>(ptFeSpace1_->GetFe( ent1 )); 
    UInt nrFncs = ptFe->GetNumFncs();
    shared_ptr<ElemShapeMap> esm = domain->GetGrid()->GetElemShapeMap( ptElem );

    // Get integration points (shortcut: from basefe instead of 
    // IntegrationScheme class)
    StdVector<LocPoint> intPoints;
    StdVector<Double> weights;
    intScheme_->GetIntPoints( Elem::GetShapeType(ptElem->type), intPoints, weights );

    elemMat.Resize( nrFncs );
    elemMat.Init();
    
    // get pointer to nonlinear BH curve approximation
    ApproxData* nlinFnc = ptMaterial->GetNonlinFncBH(MAG_PERMEABILITY);
    // get solution of current element
    sol_->GetElemSolution( magPot_, ent1 );
        
    // Loop over all integration points
    Double reluctivity = 0.0;
    Double derivReluctivity = 0.0;
    Vector<Double> elemFlux(3);
    Vector<Double> help( nrFncs );
    Double aux1, fac, *ptr1, *ptr2;
    Matrix<Double> curl;
    LocPointMapped lp;
    
    Double elemFluxAvg = 0.0,  elemFluxMin = 1e50, elemFluxMax = 0.0;
    Double elemPotAvg = 0.0, elemPotMin = 1e150, elemPotMax = 0.0;
    Double elemNuAvg = 0.0, elemNuMin = 1e150, elemNuMax = 0.0;
    for( UInt i = 0; i < intPoints.GetSize(); i++  ) {
      
      lp.Set( intPoints[i], esm );
      ptFe->GetCurlShFnc(curl, lp, lp.shapeMap->GetElem(), 1);
    
      //compute magnetic flux density
      elemFlux = curl * magPot_;
      Double Babs = elemFlux.NormL2();
      
      elemFluxAvg += Babs / (Double) intPoints.GetSize();
      if (Babs < elemFluxMin) elemFluxMin = Babs;
      if (Babs > elemFluxMax) elemFluxMax = Babs;
      
      Double aAbs = magPot_.NormL2(); 
      elemPotAvg += aAbs / (Double) intPoints.GetSize();
      if (aAbs < elemPotMin) elemPotMin = aAbs;
      if (aAbs > elemPotMax) elemPotMax = aAbs;

      if ( Babs == 0 ) {
          reluctivity = matVal_;
      }
      else {
        reluctivity = nlinFnc->EvaluateFuncNu(Babs);
      }
              
      // We now compute B^T * D * B and scale it by the determinant
      // of the Jacobian and the weight of the current integration
      // point. The result is added to the element matrix.
      fac = lp.jacDet * weights[i] * reluctivity;
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
      elemNuAvg += reluctivity / (Double) intPoints.GetSize();
      if (reluctivity < elemNuMin) elemNuMin = reluctivity;
      if (reluctivity > elemNuMax) elemNuMax = reluctivity;
      
      if ( nonLinType_ == NEWTON ) {
        if ( Babs == 0.0 ) 
          derivReluctivity = 0;
        else {          
          //Newton method
          derivReluctivity =  nlinFnc->EvaluatePrimeNu(Babs);
          fac = lp.jacDet * weights[i] * derivReluctivity * Babs;

          Vector<Double> eB(3); 
          eB = elemFlux /Babs;
          for ( UInt k = 0; k < curl.GetNumCols(); k++ ) 
            for ( UInt i = 0; i < curl.GetNumRows(); i++ ) 
              help[k] =  curl[i][k] * eB[i];

          for ( UInt i = 0; i< curl.GetNumCols(); i++ ) 
            for ( UInt j = 0; j< curl.GetNumCols(); j++ ) 
              elemMat[i][j] += fac * help[i] * help[j];
        }
      }
    }

    // write flux to file, if required
    if( fluxElems_.Find( ptElem->elemNum) != -1 && 
        logging_ == true ) {
      std::string fileName = "elemFlux-" + lexical_cast<std::string>(ptElem->elemNum);
      std::cerr << "fileName is " << fileName << std::endl;
      std::ofstream out(fileName.c_str(), std::ios::out | std::ios::app );
      out << elemFluxAvg << "\t"
          << elemFluxMin << "\t"
          << elemFluxMax << "\t"
          << elemPotAvg << "\t"
          << elemPotMin << "\t"
          << elemPotMax << "\t"
          << elemNuAvg << "\t"
          << elemNuMin << "\t"
          << elemNuMax << "\n";
      out.close();
    }
  }

  void  nLinCurlCurlEdgeInt::SetNonLinMethod(NonLinMethodType atype) {
    nonLinType_ = atype;
  }

} // end namespace CoupledField
