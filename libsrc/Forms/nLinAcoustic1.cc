#include <iostream>
#include <fstream>

#include "Utils/SmoothSpline.hh"
#include "nLinAcoustic1.hh"

namespace CoupledField {

nLinAcoustic1::nLinAcoustic1(Double factor, Boolean axi)
  : BaseForm()
{
  ENTER_FCN( "nLinAcoustic1::nLinAcoustic1");

  factor_     = factor;
  isaxi_      = axi;
  nonLinType_ = FIXEDPOINT;
}


nLinAcoustic1::~nLinAcoustic1()
{
  ENTER_FCN( "nLinAcoustic1::~nLinAcoustic1");
}


void nLinAcoustic1::CalcElementMatrix(Matrix<Double> & ptCoord, Matrix<Double> & elemMat)
{
  ENTER_FCN( "nLinAcoustic1::CalcElementMatrix");
  
  const Integer nrIntPts= ptelem->GetNumIntPoints();
  const Integer nrNodes = ptelem->GetNumNodes();
  const Vector<Double> & intWeights = ptelem->GetIntWeights();  
  Double jacDet;  

  Matrix<Double> partElemMat;
  Vector<Double> ShpFncAtIp;
  Vector<Double> CoordAtIP;

  Double solderiv1AtIP;

  // set matrix to desired size and set all elements to zero
  elemMat.Resize(nrNodes); elemMat.Init(0);

  for (Integer actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {
	
	jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord);
	ptelem->GetShFncAtIp(ShpFncAtIp,actIntPt);

	//get 1st derivartive of solution at integration point
	solderiv1AtIP = solderiv1_*ShpFncAtIp;
	
	partElemMat.DyadicMult(ShpFncAtIp);

	if (isaxi_) {
	  CoordAtIP = ptCoord * ShpFncAtIp;
	  partElemMat *= 2 * PI * intWeights[actIntPt-1] * factor_
                     * jacDet * CoordAtIP[0] * solderiv1AtIP;
	}
	else 
	  partElemMat *= intWeights[actIntPt-1] * factor_ 
                     * jacDet * solderiv1AtIP;
	
	elemMat += partElemMat;
  }
}


void nLinAcoustic1::SetNonLinMethod(std::string atype)
{
  ENTER_FCN( "nLinAcoustic1::SetNonLinMethod");
    
  if (atype == "fixPoint")
	nonLinType_ = FIXEDPOINT;
    
}


void nLinAcoustic1::Print(std::ostream * out, const Matrix<Double> Result) const
{
  ENTER_FCN( "nLinAcoustic1::Print");

  (*out)<< "nLinAcoustic1 matrix:" << std::endl << Result;
}

} // end namespace CoupledField
