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

    

    //     //check for material value
    //     if (materialArray_ != NULL) {
    //       laplVal_ = (*materialArray_)[actSD_][actElemNr_];
    //     }


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
              * jacDet * laplVal_ * CoordAtIP[0];

          }
        else 
          partElemMat *= intWeights[actIntPt-1] * jacDet * laplVal_;

        elemMat += partElemMat;
      }

    //    std::cout << "ElemMatLaplace:\n" << elemMat << std::endl;
  }



} // end namespace CoupledField
