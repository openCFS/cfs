// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include "laplaceXYZInt.hh"

namespace CoupledField
{

  LaplaceXYZInt::LaplaceXYZInt(Double aVal, const UInt nrDofsPerNode, bool axi)
    : BaseForm( NULL ),nrDofsPerNode_(nrDofsPerNode),laplVal_ (aVal)
  {
    
    name_ = "LaplaceXYZInt";
    
    isaxi_ = axi;
  }


 
  LaplaceXYZInt::~LaplaceXYZInt()
  {
  }



  void LaplaceXYZInt::CalcElementMatrix( Matrix<Double>& elemMat,
                                         EntityIterator& ent1, 
                                         EntityIterator& ent2 )
  {

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );
    
    ptelem->SetAnsatzFct( ansatzFct1_ );
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const UInt nrIntPts= ptelem->GetNumIntPoints();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Double jacDet, factor;  


    // derivation of shape functions after global coordinates 
    Matrix<Double> xiDx;
    Matrix<Double> elemMatXYZ;
    Vector<Double> ShpFncAtIp;
    Vector<Double> CoordAtIP;


    // set matrix to desired size and set all elements to zero
    elemMat.Resize(numFncs*nrDofsPerNode_,true); 
    elemMat.Init();

    //(Kx, Ky, Kz)^T
    elemMatXYZ.Resize(numFncs*nrDofsPerNode_,numFncs); 
    elemMatXYZ.Init();

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {
        jacDet = 0;
        
        ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord_, 
                                      jacDet, ent1.GetElem());
      
        if (isaxi_) {
	  ptelem->GetShFncAtIp(ShpFncAtIp,actIntPt,ent1.GetElem());
	  CoordAtIP = ptCoord_ * ShpFncAtIp;
	  factor = 2 * PI * intWeights[actIntPt-1] * jacDet * laplVal_ * CoordAtIP[0];
	}
        else {
          factor = intWeights[actIntPt-1] * jacDet * laplVal_;
	}

	for (UInt i=0; i<numFncs; i++) {
	  for (UInt j=0; j<numFncs; j++) {
	    for (UInt k=0; k<nrDofsPerNode_; k++) {
	      elemMatXYZ[i+k*numFncs][j] += xiDx[i][k]*xiDx[j][k]*factor;
	    }
	  }
	}

      }

    for (UInt i=0; i<numFncs; i++) {
      for (UInt j=0; j<numFncs; j++) {
	for (UInt k=0; k<nrDofsPerNode_; k++) {
	  for (UInt l=0; l<nrDofsPerNode_; l++) {
	    elemMat[i*nrDofsPerNode_+k][j*nrDofsPerNode_+l]   = elemMatXYZ[i+k*numFncs][j];
	  }
	}
      }
    }

    //    std::cout << "ElemMatLaplace:\n" << elemMatXYZ << std::endl;
    //    std::cout << "ElemMatLaplaceXYZ:\n" << elemMat << std::endl;
  }

} // end namespace CoupledField
