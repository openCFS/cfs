#include <iostream>
#include <fstream>

#include "laplaceInt.hh"

namespace CoupledField
{

  LaplaceInt::LaplaceInt(BaseFE * aptelem, Double aVal, Boolean axi)
    : BaseForm(aptelem),laplVal_ (aVal), isaxi_(axi)
  {
#ifdef TRACE
    (*trace) << "entering LaplaceInt::LaplaceInt" << std::endl;
#endif
  }


  LaplaceInt::LaplaceInt(Double aVal, Boolean axi)
    : BaseForm(),laplVal_ (aVal), isaxi_(axi)
  {
#ifdef TRACE
    (*trace) << "entering LaplaceInt::LaplaceInt" << std::endl;
#endif
  }


 
  LaplaceInt::~LaplaceInt()
  {
#ifdef TRACE
    (*trace) << "entering LaplaceInt::~LaplaceInt" << std::endl;
#endif
  }



  void LaplaceInt::CalcElementMatrix(Matrix<Double> & ptCoord, Matrix<Double> & elemMat)
  {
#ifdef TRACE
    (*trace) << "entering LaplaceInt::CalcElementMatrix" << std::endl;
#endif
  
    const Integer nrIntPts= ptelem->GetNumIntPoints();
    const Integer nrNodes = ptelem->GetNumNodes();
    const std::vector<Double> & intWeights = ptelem->GetIntWeights();  
    Double jacDet;  


    // derivation of shape functions after global coordinates 
    Matrix<Double> xiDx;
    Matrix<Double> xiDxTransp;
    Matrix<Double> partElemMat;
    std::vector<Double> ShpFncAtIp;
    std::vector<Double> CoordAtIP;

    // set matrix to desired size and set all elements to zero
    elemMat.Resize(nrNodes); elemMat.Init();
    
    for (Integer actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {
	jacDet = 0;
	
	ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord, jacDet);

	xiDx.Transpose(xiDxTransp);

	partElemMat = xiDx * xiDxTransp;
	
	if (isaxi_)
	  {
	    ptelem->GetShFncAtIp(ShpFncAtIp,actIntPt);
	    CoordAtIP = ptCoord * ShpFncAtIp;
	    partElemMat *= 2 * PI * intWeights[actIntPt-1] * jacDet * laplVal_ * CoordAtIP[0];

	  }
	else 
	  partElemMat *= intWeights[actIntPt-1] * jacDet * laplVal_;

	elemMat += partElemMat;
      }
  

#ifdef TRACE
    (*trace) << "leaving LaplaceInt::CalcElemMatrix" << std::endl;
#endif
  }



  void LaplaceInt::Print(std::ostream * out, const Matrix<Double> Result) const
  {
#ifdef TRACE
    (*trace) << "entering LaplaceInt::Print" << std::endl;
#endif
    (*out)<< "Laplace stiffness matrix:" << std::endl << Result;
  }

} // end namespace CoupledField
