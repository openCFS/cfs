// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include "laplaceInt.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"

namespace CoupledField
{


  LaplaceInt::LaplaceInt(Double aVal, bool axi, bool coordUpdate )
    : BaseForm(NULL),laplVal_ (aVal)
  {
    name_ = "LaplaceInt";
    isaxi_ = axi;
    coordUpdate_ = coordUpdate;
    if ( coordUpdate )
      isSolDependent_ = true;
  }

  LaplaceInt::LaplaceInt(BaseMaterial* mat, const MaterialDescriptor& md, bool axi, bool coordUpdate )
    : BaseForm(mat)
  {
    name_ = "LaplaceInt";
    isaxi_ = axi;
    coordUpdate_ = coordUpdate;
    if ( coordUpdate ) 
      isSolDependent_ = true;

    md_ = md;
    laplVal_ = md_.GetScalar(mat); // we won't use it but do it always again and again :)
  }

 
  LaplaceInt::~LaplaceInt()
  {
  }



  void LaplaceInt::CalcElementMatrix( Matrix<Double>& elemMat,
                                      EntityIterator& ent1, 
                                      EntityIterator& ent2 ) {
  
   
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );

    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt nrFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const UInt nrIntPts= ptelem->GetNumIntPoints();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Double jacDet;  


    // derivation of shape functions after global coordinates 
    Matrix<Double> xiDx;
    Matrix<Double> xiDxTransp;
    Matrix<Double> partElemMat;
    Vector<Double> ShpFncAtIp;
    Vector<Double> CoordAtIP;

    // set matrix to desired size and set all elements to zero
    elemMat.Resize(nrFncs); 
    elemMat.Init();
    
    // we call the laplace value density and it might contain an topology optimization
    // ersatz material factor
    Double density = md_.GetErsatzMaterial(this, ent1.GetElem(), laplVal_);

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
    {
        jacDet = 0;
        
        ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord_, 
                                      jacDet, ent1.GetElem());

        xiDx.Transpose(xiDxTransp);

        partElemMat = xiDx * xiDxTransp;
        
        if (isaxi_)
          {
            ptelem->GetShFncAtIp(ShpFncAtIp,actIntPt,ent1.GetElem());
            CoordAtIP = ptCoord_ * ShpFncAtIp;
            partElemMat *= 2 * PI * intWeights[actIntPt-1] 
              * jacDet * density * CoordAtIP[0];

          }
        else 
          partElemMat *= intWeights[actIntPt-1] * jacDet * density;

        elemMat += partElemMat;
    }
  }



} // end namespace CoupledField
