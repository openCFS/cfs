#include <iostream>
#include <fstream>

#include "curlCurlNodeInt.hh"

namespace CoupledField
{

  CurlCurlNode2DInt::CurlCurlNode2DInt(BaseFE * aptelem, Double aVal, Boolean axi)
    : BaseForm(aptelem),matVal_ (aVal), isaxi_(axi)
  {
#ifdef TRACE
    (*trace) << "entering CurlCurlNode2DInt::CurlCurlNode2DInt" << std::endl;
#endif
  }


  CurlCurlNode2DInt::CurlCurlNode2DInt(Double aVal, Boolean axi)
    : BaseForm(),matVal_ (aVal), isaxi_(axi)
  {
#ifdef TRACE
    (*trace) << "entering CurlCurlNode2DInt::CurlCurlNode2DInt" << std::endl;
#endif
  }


 
  CurlCurlNode2DInt::~CurlCurlNode2DInt()
  {
#ifdef TRACE
    (*trace) << "entering CurlCurlNode2DInt::~CurlCurlNode2DInt" << std::endl;
#endif
  }



  void CurlCurlNode2DInt::CalcElementMatrix(Matrix<Double> & ptCoord, Matrix<Double> & elemMat)
  {
#ifdef TRACE
    (*trace) << "entering CurlCurlNode2DInt::CalcElementMatrix" << std::endl;
#endif
  
    const Integer nrIntPts= ptelem->GetNumIntPoints();
    const Integer nrNodes = ptelem->GetNumNodes();
    const std::vector<Double> & intWeights = ptelem->GetIntWeights();  
    Double jacDet;  


    // derivation of shape functions after global coordinates 
    Matrix<Double> xiDx;
    Matrix<Double> xiDxTransp;
    Matrix<Double> partElemMat;
    Matrix<Double> partElemMatAxi;
    std::vector<Double> ShpFncAtIp;
    std::vector<Double> CoordAtIP;
    std::vector<Double> drAtIp;


    // set matrix to desired size and set all elements to zero
    elemMat.Resize(nrNodes); elemMat.Init();
    
    for (Integer actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {
	jacDet = 0;
	
	ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord, jacDet);

	if (isaxi_)
	  {
	    ptelem->GetShFncAtIp(ShpFncAtIp,actIntPt);
	    CoordAtIP = ptCoord * ShpFncAtIp;
            for (Integer i=0; i<nrNodes; i++)
		xiDx[i][0] += ShpFncAtIp[i] / CoordAtIP[0];
            
            jacDet *= 2 * PI * CoordAtIP[0];
	}
  
	xiDx.Transpose(xiDxTransp);
	partElemMat = xiDx * xiDxTransp;
	partElemMat *= intWeights[actIntPt-1] * jacDet * matVal_;
	elemMat += partElemMat;
      }
  

#ifdef TRACE
    (*trace) << "leaving CurlCurlNode2DInt::CalcElemMatrix" << std::endl;
#endif
  }



  void CurlCurlNode2DInt::Print(std::ostream * out, const Matrix<Double> Result) const
  {
#ifdef TRACE
    (*trace) << "entering CurlCurlNode2DInt::Print" << std::endl;
#endif
    (*out)<< "Laplace stiffness matrix:" << std::endl << Result;
  }

} // end namespace CoupledField
