#include <iostream>
#include <fstream>

#include "laplaceInt.hh"

namespace CoupledField
{

  LaplaceInt::LaplaceInt(BaseFE * aptelem, Double aVal, Boolean axi)
    : BaseForm(aptelem),laplVal_ (aVal)
  {
    ENTER_FCN( "LaplaceInt::LaplaceInt");

    isaxi_ = axi;
  }


  LaplaceInt::LaplaceInt(Double aVal, Boolean axi)
    : BaseForm(),laplVal_ (aVal)
  {
    ENTER_FCN( "LaplaceInt::LaplaceInt" );

    isaxi_ = axi;
  }


 
  LaplaceInt::~LaplaceInt()
  {
    ENTER_FCN( "LaplaceInt::~LaplaceInt" );
  }



  void LaplaceInt::CalcElementMatrix(Matrix<Double> & ptCoord, Matrix<Double> & elemMat)
  {
    ENTER_FCN( "LaplaceInt::CalcElementMatrix" );
  
    const Integer nrIntPts= ptelem->GetNumIntPoints();
    const Integer nrNodes = ptelem->GetNumNodes();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Double jacDet;  


    // derivation of shape functions after global coordinates 
    Matrix<Double> xiDx;
    Matrix<Double> xiDxTransp;
    Matrix<Double> partElemMat;
    Vector<Double> ShpFncAtIp;
    Vector<Double> CoordAtIP;

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
  
  }



  void LaplaceInt::Print(std::ostream * out, const Matrix<Double> Result) const
  {
    ENTER_FCN( "LaplaceInt::Print"); 
    (*out)<< "Laplace stiffness matrix:" << std::endl << Result;
  }

} // end namespace CoupledField
