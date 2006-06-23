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

    const UInt nrIntPts= ptelem->GetNumIntPoints();
    const UInt nrNodes = ptelem->GetNumNodes();
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
    elemMat.Resize(nrNodes); 
    elemMat.Init();
    
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {
        jacDet = 0;
        
        ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord_, jacDet);

        if (isaxi_)
          {
            ptelem->GetShFncAtIp(ShpFncAtIp,actIntPt);
            CoordAtIP = ptCoord_ * ShpFncAtIp;
            for (UInt i=0; i<nrNodes; i++)
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
