#include <iostream>
#include <fstream>

#include "laplaceInt.hh"

namespace CoupledField
{

  LaplaceInt::LaplaceInt(BaseFE * aptelem, Double aVal)
    : BaseForm(aptelem),laplVal_ (aVal)
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
    double jacDet;  


    // derivation of shape functions after global coordinates 
    Matrix<Double> xiDx;
    Matrix<Double> xiDxTransp;
    Matrix<Double> partElemMat;
  


    // set matrix to desired size and set all elements to zero
    //    partElemMat.Resize(nrNodes);
    elemMat.Resize(nrNodes); elemMat.Init();
    

    for (Integer actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {
	jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord);

	ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord);

	xiDx.Transpose(xiDxTransp);

	partElemMat = xiDx * xiDxTransp;
	
	partElemMat *= intWeights[actIntPt-1] * jacDet * laplVal_;

#ifdef DEBUG 
	(*debug) << "Partelemmat on intPt " << actIntPt << std::endl
		 << partElemMat << std::endl
		 << "xiDx \n" << xiDx
		 << "\n xiDxTransp \n " << xiDxTransp 
		 << "\n intWeights " << intWeights[actIntPt-1] << std::endl;
	
#endif

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
