// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include "nlinLaplaceInt.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"
#include "Utils/nodestoresol.hh"

namespace CoupledField
{


  nLinLaplaceInt::nLinLaplaceInt(Double aVal, bool axi,  bool isSpeedVariable, 
                                 bool coordUpdate )
    : BaseForm(NULL),laplVal_ (aVal)
  {

    name_ = "nLinLaplaceInt";

    isSolDependent_ = false;
    if ( isSpeedVariable ) 
      isSolDependent_ = true;

    isaxi_ = axi;
    coordUpdate_ = coordUpdate;
  }


 
  nLinLaplaceInt::~nLinLaplaceInt()
  {
  }



  void nLinLaplaceInt::CalcElementMatrix( Matrix<Double>& elemMat,
                                      EntityIterator& ent1, 
                                      EntityIterator& ent2 ) {
  
   
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );
   
    const UInt nrIntPts= ptelem->GetNumIntPoints();
    const UInt nrNodes = ptelem->GetNumNodes();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Double jacDet;  


    // derivation of shape functions after global coordinates 
    Matrix<Double> xiDx;
    Matrix<Double> xiDxTransp;
    Matrix<Double> partElemMat;
    Vector<Double> ShpFncAtIp;
    Vector<Double> CoordAtIP;

    // set matrix to desired size and set all elements to zero
    elemMat.Resize(nrNodes); 
    elemMat.Init();

    //     //check for material value
    //     if (materialArray_ != NULL) {
    //       laplVal_ = (*materialArray_)[actSD_][actElemNr_];
    //     }

    // get the variable speed of sound, stored in the solution object
    Matrix<Double> temp;
    Vector<Double> sos; // speed of sound
    Double factor = 0.0;

    sol_->GetElemSolutionAsMatrix( temp, ent1 );
    temp.ConvertToVec_AppendRows( sos );
    for ( UInt i=0; i<sos.GetSize(); i++ ) {
       factor += sos[i];
    }
    factor /= (Double)sos.GetSize(); 
//    laplVal_ *= factor * factor;
    laplVal_ = factor*factor;

    //Info->PrintF( "", "SoS = %e\n", laplVal_);

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {
        jacDet = 0;
        
        ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord_, 
                                      jacDet, ent1.GetElem() );

        xiDx.Transpose(xiDxTransp);

        partElemMat = xiDx * xiDxTransp;
        
        if (isaxi_) {
	  ptelem->GetShFncAtIp( ShpFncAtIp, actIntPt, ent1.GetElem() );
	  CoordAtIP = ptCoord_ * ShpFncAtIp;
	  partElemMat *= 2 * PI * intWeights[actIntPt-1] 
	    * jacDet * laplVal_ * CoordAtIP[0];
	}
        else 
          partElemMat *= intWeights[actIntPt-1] * jacDet * laplVal_;

        elemMat += partElemMat;
      }

    //    std::cout << "nlinLaplace" << std::endl;
    //    std::cout << "ElemMatnLinnLinLaplace:\n" << elemMat << std::endl;
  }



} // end namespace CoupledField
