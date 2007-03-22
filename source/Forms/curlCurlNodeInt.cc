// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include "curlCurlNodeInt.hh"

namespace CoupledField
{

  CurlCurlNode2DInt::CurlCurlNode2DInt(Double aVal, bool axi,
                                       bool coordUpdate )
    : BaseForm(NULL,FULL,coordUpdate ),matVal_ (aVal)
  {
    ENTER_FCN( "CurlCurlNode2DInt::CurlCurlNode2DInt" );

    name_ = "CurlCurlNode2DInt";
    isaxi_ = axi;
  }


 
  CurlCurlNode2DInt::~CurlCurlNode2DInt()
  {
    ENTER_FCN( "CurlCurlNode2DInt::~CurlCurlNode2DInt" );
  }



  void CurlCurlNode2DInt::CalcElementMatrix( Matrix<Double>& elemMat,
                                             EntityIterator& ent1, 
                                             EntityIterator& ent2  ) {
    ENTER_FCN( "CurlCurlNode2DInt::CalcElementMatrix" );
  
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
        partElemMat *= intWeights[actIntPt-1] * jacDet * matVal_;
        elemMat += partElemMat;
      }
  
  }

} // end namespace CoupledField
