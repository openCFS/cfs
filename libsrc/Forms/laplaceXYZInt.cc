#include <iostream>
#include <fstream>

#include "laplaceXYZInt.hh"

namespace CoupledField
{

  LaplaceXYZInt::LaplaceXYZInt(Double aVal, const UInt nrDofsPerNode, bool axi)
    : BaseForm( NULL ),nrDofsPerNode_(nrDofsPerNode),laplVal_ (aVal)
  {
    ENTER_FCN( "LaplaceXYZInt::LaplaceXYZInt" );
    
    name_ = "LaplaceXYZInt";
    
    isaxi_ = axi;
  }


 
  LaplaceXYZInt::~LaplaceXYZInt()
  {
    ENTER_FCN( "LaplaceXYZInt::~LaplaceXYZInt" );
  }



  void LaplaceXYZInt::CalcElementMatrix( Matrix<Double>& elemMat,
                                         EntityIterator& ent1, 
                                         EntityIterator& ent2 )
  {
    ENTER_FCN( "LaplaceXYZInt::CalcElementMatrix" );

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );
    
    const UInt nrIntPts= ptelem->GetNumIntPoints();
    const UInt nrNodes = ptelem->GetNumNodes();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Double jacDet, factor;  


    // derivation of shape functions after global coordinates 
    Matrix<Double> xiDx;
    Matrix<Double> elemMatXYZ;
    Vector<Double> ShpFncAtIp;
    Vector<Double> CoordAtIP;


    // set matrix to desired size and set all elements to zero
    elemMat.Resize(nrNodes*nrDofsPerNode_,true); elemMat.Init();

    //(Kx, Ky, Kz)^T
    elemMatXYZ.Resize(nrNodes*nrDofsPerNode_,nrNodes); elemMatXYZ.Init();

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {
        jacDet = 0;
        
        ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord_, jacDet);
      
        if (isaxi_) {
	  ptelem->GetShFncAtIp(ShpFncAtIp,actIntPt);
	  CoordAtIP = ptCoord_ * ShpFncAtIp;
	  factor = 2 * PI * intWeights[actIntPt-1] * jacDet * laplVal_ * CoordAtIP[0];
	}
        else {
          factor = intWeights[actIntPt-1] * jacDet * laplVal_;
	}

	for (UInt i=0; i<nrNodes; i++) {
	  for (UInt j=0; j<nrNodes; j++) {
	    for (UInt k=0; k<nrDofsPerNode_; k++) {
	      elemMatXYZ[i+k*nrNodes][j] += xiDx[i][k]*xiDx[j][k]*factor;
	    }
	  }
	}

      }

    for (UInt i=0; i<nrNodes; i++) {
      for (UInt j=0; j<nrNodes; j++) {
	for (UInt k=0; k<nrDofsPerNode_; k++) {
	  for (UInt l=0; l<nrDofsPerNode_; l++) {
	    elemMat[i*nrDofsPerNode_+k][j*nrDofsPerNode_+l]   = elemMatXYZ[i+k*nrNodes][j];
	  }
	}
      }
    }

    //    std::cout << "ElemMatLaplace:\n" << elemMatXYZ << std::endl;
    //    std::cout << "ElemMatLaplaceXYZ:\n" << elemMat << std::endl;
  }

} // end namespace CoupledField
