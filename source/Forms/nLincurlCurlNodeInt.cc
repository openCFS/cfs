// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <string>

#include "Domain/entityList.hh"
#include "Elements/basefe.hh"
#include "Forms/curlCurlNodeInt.hh"
#include "General/exception.hh"
#include "MatVec/exprt/xpr1.hh"
#include "MatVec/exprt/xpr2.hh"
#include "MatVec/matrix.hh"
#include "Materials/baseMaterial.hh"
#include "Utils/ApproxData.hh"
#include "Utils/nodestoresol.hh"
#include "nLincurlCurlNodeInt.hh"

namespace CoupledField
{


  nLinCurlCurlNode2DInt::
  nLinCurlCurlNode2DInt( BaseMaterial* matData, 
                         bool axi, bool coordUpdate )
    : CurlCurlNode2DInt( matData, axi, coordUpdate )
  {

    name_ = "nLinCurlCurlNode2DInt";
    isSolDependent_ = true;
    nonLinType_ = NEWTON;

    ptMaterial->GetScalar( startmatVal_, MAG_RELUCTIVITY,Global::REAL);

    isHysteresis_ = false;
    if ( ptMaterial->IsSetHysteresis() ) {
      isHysteresis_ = true;
    }
  }

 
  nLinCurlCurlNode2DInt::~nLinCurlCurlNode2DInt()
  {
  }


  void nLinCurlCurlNode2DInt::CalcElementMatrix( Matrix<Double>& elemMat,
                                                  EntityIterator& ent1, 
                                                  EntityIterator& ent2 )
  {

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );  
    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt nrIntPts= ptelem->GetNumIntPoints();
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Double jacDet;  


    // get solution of current element
    Matrix<Double> temp;
    sol_->GetElemSolutionAsMatrix( temp, ent1 );
    temp.ConvertToVec_AppendRows( magPot_ );


    // derivation of shape functions after global coordinates 
    Matrix<Double> xiDx;
    Matrix<Double> xiDxTransp;
    Matrix<Double> partElemMat;
    Matrix<Double> partElemMatAxi;
    Vector<Double> ShpFncAtIp;
    Vector<Double> CoordAtIP;
    Vector<Double> drAtIp;

    Double reluctivity, derivReluctivity = 0.0;

    // set matrix to desired size and set all elements to zero
    elemMat.Resize(numFncs); elemMat.Init();

    // get pointer to nonlinear BH curve approximation
    nlinFnc_ = ptMaterial->GetNonlinFnc(MAG_PERMEABILITY);

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {
        jacDet = 0;
        
        ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord_, 
                                      jacDet, ent1.GetElem() );

        if (isaxi_)
          {
            ptelem->GetShFncAtIp(ShpFncAtIp,actIntPt, ent1.GetElem() );
            CoordAtIP = ptCoord_ * ShpFncAtIp;
            for (UInt i=0; i<numFncs; i++)
              xiDx[i][0] += ShpFncAtIp[i] / CoordAtIP[0];
            
            jacDet *= 2 * PI * CoordAtIP[0];
          }
  
        xiDx.Transpose(xiDxTransp);
        partElemMat = xiDx * xiDxTransp;

        //compute value for nonlinear reluctivity
        Vector<Double> B(2);
        UInt dim = 2;
        for( UInt i=0; i<dim; i++ )
          for( UInt j=0; j<numFncs; j++ )
            B[i] += xiDx[j][i] * magPot_[j];


        Double Babs = B.NormL2();

        if ( isHysteresis_ ) {
          EXCEPTION("Magnetics with hysteresis currently not supported");
//           //hysteresis modeling

//           // Account for curl 
//           Double temp = B[0];
//           if ( isaxi_ ) {
//             B[0] = -B[1];
//             B[1] = temp;
//           } else {
//             B[0] = B[1];
//             B[1] = -temp;
//           }
//           UInt nrEl = ent1.GetElem()->elemNum;
//           reluctivity = ComputeDiffReluctivity( nrEl, B );
//           //std::cout << "Bfield:\n" << B << "\n nu=" << reluctivity << std::endl;
        }
        else {
          //nonlinear BH curve
          if (Babs ==0) 
            reluctivity = startmatVal_;
          else {
            reluctivity = nlinFnc_->EvaluateFuncNu(Babs);
          }
        }

        //        std::cout << "b=" << Babs << "  nu=" << reluctivity << std::endl;
        partElemMat *= reluctivity;
        
        if ( !isHysteresis_ && nonLinType_ == NEWTON) {
          if (Babs ==0) 
            derivReluctivity = 0;
          else {          
            //Newton method
            Vector<Double> eB(2); eB = B * (1/Babs);
            derivReluctivity =  nlinFnc_->EvaluatePrimeNu(Babs);

            //            std::cout << "Newton, B=" << Babs << "   nuprime=" << derivReluctivity << std::endl;
            for (UInt p=0;  p<numFncs; p++)
              for (UInt q=0; q<numFncs; q++) {               
                partElemMat[p][q] +=  derivReluctivity * 
                  (eB[0]*eB[0]*xiDx[p][1]*xiDx[q][1] +
                   eB[1]*eB[1]*xiDx[p][0]*xiDx[q][0] -
                   eB[0]*eB[1]*xiDx[p][1]*xiDx[q][0] -
                   eB[1]*eB[0]*xiDx[p][0]*xiDx[q][1] ) * Babs;
              }
          }
        }

        partElemMat *= intWeights[actIntPt-1] * jacDet;
        elemMat += partElemMat;
      }
  }


  void nLinCurlCurlNode2DInt::SetNonLinMethod(NonLinMethodType atype) {
      nonLinType_ = atype;
    }


  Double nLinCurlCurlNode2DInt::ComputeDiffReluctivity( UInt nrEl, Vector<Double>& Bvec )
  {

    Double diffRelucVal = 1.0;

    //    diffRelucVal = ptMaterial->ComputeScalarDiffVal( nrEl, Bvec );

    if (  diffRelucVal <= 0.0 ) 
      EXCEPTION("Negative effective permeability");


    return diffRelucVal;

  }


  //============================Curl-Curl-3D ====================================

  nLinCurlCurlNode3DInt::
  nLinCurlCurlNode3DInt( BaseMaterial* matData, bool coordUpdate )
    : CurlCurlNode3DInt( matData, coordUpdate )
  {

    name_ = "nLinCurlCurlNode3DInt";
    isaxi_      = false;
    coordUpdate_ = coordUpdate;
    isSolDependent_ = true;
    nonLinType_ = NEWTON;

    ptMaterial->GetScalar( startmatVal_, MAG_RELUCTIVITY,Global::REAL);

  }

 
  nLinCurlCurlNode3DInt::~nLinCurlCurlNode3DInt()
  {
  }


  void nLinCurlCurlNode3DInt::CalcElementMatrix( Matrix<Double>& elemMat,
                                                  EntityIterator& ent1, 
                                                  EntityIterator& ent2 )
  {

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );

    ptelem->SetAnsatzFct( ansatzFct1_ );
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    
    Double jacDet, reluctivity, derivReluctivity;  

    Matrix<Double> bMatCurl, bMatDiv; 
    Double aux1, aux2, fac, *ptr1, *ptr2, *ptr3, *ptr4;

    elemMat.Resize( numFncs * nrDofs_ );
    elemMat.Init();

    Vector<Double> elemFlux(3);
    Vector<Double> help( numFncs * nrDofs_ );

//     Matrix<Double> temp;
//     sol_->GetElemSolutionAsMatrix( temp, ent1 );
//     temp.ConvertToVec_AppendRows( magPot_ );

    // get solution of current element
    sol_->GetElemSolution( magPot_, ent1 );

    //get integration points
    const UInt nrIntPts = ptelem->GetNumIntPoints(); 
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  

    // get pointer to nonlinear BH curve approximation
    nlinFnc_ = ptMaterial->GetNonlinFnc(MAG_PERMEABILITY);

    // Loop over all integration points
    for ( UInt actIntPt = 1; actIntPt <= nrIntPts; actIntPt++ ) {

      // Setup the B matrix for current integration point
      calcBMat( bMatCurl, bMatDiv, actIntPt, ptCoord_ );

      // Compute Jacobian for integration point
      jacDet = ptelem->CalcJacobianDetAtIp( actIntPt, ptCoord_, ent1.GetElem() );

      // Perform a safety check
      if ( jacDet < 0.0 ) {
        EXCEPTION("CurlCurlNode3DInt::CalcElementMatrix: Encountered "
            << "negative Jacobian determinant!");
      }

      //compute magnetic flux density
      elemFlux = bMatCurl * magPot_;
      Double Babs = elemFlux.NormL2();

      if ( Babs == 0 ) 
        reluctivity = startmatVal_;
      else {
        reluctivity = nlinFnc_->EvaluateFuncNu(Babs);
      }
           
      // We now compute B^T * D * B and scale it by the determinant
      // of the Jacobian and the weight of the current integration
      // point. The result is added to the element matrix.
      fac = jacDet * intWeights[actIntPt-1] * reluctivity;
      for ( UInt k = 0; k < bMatCurl.GetNumRows(); k++ ) {
	ptr1 = bMatCurl[k];
	ptr2 = bMatCurl[k];
	ptr3 = bMatDiv[k];
	ptr4 = bMatDiv[k];
	for ( UInt i = 0; i < bMatCurl.GetNumCols(); i++ ) {
	  aux1 = fac * ptr1[i];
	  aux2 = fac * ptr3[i];
	  for ( UInt j = 0; j < bMatCurl.GetNumCols(); j++ ) {
	    elemMat[i][j] += aux1 * ptr2[j] + aux2 * ptr4[j];
	  }
	}
      }

      if ( nonLinType_ == NEWTON ) {
        if ( Babs == 0.0 ) 
          derivReluctivity = 0;
        else {          
          //Newton method
          Vector<Double> eB(3); eB = elemFlux * (1/Babs);
          derivReluctivity =  nlinFnc_->EvaluatePrimeNu(Babs);
          fac = jacDet * intWeights[actIntPt-1] * derivReluctivity * Babs;

          for ( UInt k = 0; k < bMatCurl.GetNumCols(); k++ ) 
            for ( UInt i = 0; i < bMatCurl.GetNumRows(); i++ ) 
              help[k] =  bMatCurl[i][k] * eB[i];

          for ( UInt i = 0; i< bMatCurl.GetNumCols(); i++ ) 
            for ( UInt j = 0; j< bMatCurl.GetNumCols(); j++ ) 
              elemMat[i][j] += fac * help[i] * help[j];
        }
      }
    }
    
  }


  void nLinCurlCurlNode3DInt::SetNonLinMethod(NonLinMethodType atype) {
    nonLinType_ = atype;
  }



} // end namespace CoupledField
