#include <iostream>
#include <fstream>

#include "laplaceXYZInt.hh"

namespace CoupledField
{

  LaplaceXYZInt::LaplaceXYZInt(BaseFE * aptelem, Double aVal, Boolean axi)
    : BaseForm(aptelem),laplVal_ (aVal)
  {
    ENTER_FCN( "LaplaceXYZInt::LaplaceXYZInt");

    Error("This LaplaceXYZInt-Constructor is no longer allowed",__FILE__,__LINE__);
  }


  LaplaceXYZInt::LaplaceXYZInt(Double aVal, const UInt nrDofsPerNode, Boolean axi)
    : BaseForm(),laplVal_ (aVal),nrDofsPerNode_(nrDofsPerNode)
  {
    ENTER_FCN( "LaplaceXYZInt::LaplaceXYZInt" );

    
    isaxi_ = axi;
  }


 
  LaplaceXYZInt::~LaplaceXYZInt()
  {
    ENTER_FCN( "LaplaceXYZInt::~LaplaceXYZInt" );
  }



  void LaplaceXYZInt::CalcElementMatrix(Matrix<Double> & ptCoord, Matrix<Double> & elemMat)
  {
    ENTER_FCN( "LaplaceXYZInt::CalcElementMatrix" );
  
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
    elemMat.Resize(nrNodes*nrDofsPerNode_); elemMat.Init();

    //(Kx, Ky, Kz)^T
    elemMatXYZ.Resize(nrNodes*nrDofsPerNode_,nrNodes); elemMatXYZ.Init();

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {
        jacDet = 0;
        
        ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord, jacDet);
      
        if (isaxi_) {
	  ptelem->GetShFncAtIp(ShpFncAtIp,actIntPt);
	  CoordAtIP = ptCoord * ShpFncAtIp;
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



  void LaplaceXYZInt::Print(std::ostream * out, const Matrix<Double> Result) const
  {
    ENTER_FCN( "LaplaceXYZInt::Print"); 
    (*out)<< "LaplaceXYZ stiffness matrix:" << std::endl << Result;
  }

} // end namespace CoupledField
