#include <iostream>
#include <fstream>

#include "Utils/SmoothSpline.hh"
#include "nLinAcoustic2.hh"

namespace CoupledField {

nLinAcoustic2::nLinAcoustic2(Double factor, Boolean axi)
  : BaseForm()
{
  ENTER_FCN( "nLinAcoustic2::nLinAcoustic2");

  factor_     = factor;
  isaxi_      = axi;
  nonLinType_ = FIXEDPOINT;
}

 
nLinAcoustic2::~nLinAcoustic2()
{
  ENTER_FCN( "nLinAcoustic2::~nLinAcoustic2");
}


void nLinAcoustic2::CalcElementMatrix(Matrix<Double> & ptCoord, Matrix<Double> & elemMat)
{
  ENTER_FCN( "nLinAcoustic2::CalcElementMatrix");
  
  const Integer nrIntPts= ptelem->GetNumIntPoints();
  const Integer nrNodes = ptelem->GetNumNodes();
  const Vector<Double> & intWeights = ptelem->GetIntWeights();  
  Double jacDet;  


  // derivation of shape functions after global coordinates 
  Matrix<Double> xiDx;
  Matrix<Double> xiDxTransp;
  //  Matrix<Double> partElemMat(nrNodes,nrNodes);
  //  partElemMat.Init(0);
  Vector<Double> ShpFncAtIp;
  Vector<Double> CoordAtIP;
  Vector<Double> solGradAtIP;
  Vector<Double> vec1(nrNodes);

  // set matrix to desired size and set all elements to zero
  elemMat.Resize(nrNodes); elemMat.Init(0);

  for (Integer actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {
	 
	jacDet = 0;
	ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord, jacDet);
	ptelem->GetShFncAtIp(ShpFncAtIp,actIntPt);

	if (isaxi_) {

	  CoordAtIP = ptCoord * ShpFncAtIp;
	  for (Integer i=0; i<nrNodes; i++)
		xiDx[i][0] += ShpFncAtIp[i] / CoordAtIP[0];
            
	  jacDet *= 2 * PI * CoordAtIP[0];
	}

	xiDx.Transpose(xiDxTransp);

	//compute gradient of solution at integration point
	solGradAtIP = xiDxTransp * sol_;

	for (Integer i=0; i< nrNodes; i++)
	  for (Integer j=0; j<xiDx.GetSizeCol(); j++)
		vec1[i] += xiDx[i][j]*solGradAtIP[j];

	Matrix<Double> mat2;
	mat2.DyadicMult(ShpFncAtIp,vec1);
	
	elemMat += mat2 * jacDet;

	//	std::cout << "partElemMat:" << std::endl << partElemMat << std::endl;

  }
}


void nLinAcoustic2::SetNonLinMethod(std::string atype)
{
  ENTER_FCN( "nLinAcoustic2::SetNonLinMethod");
    
  if (atype == "fixPoint")
	nonLinType_ = FIXEDPOINT;
    
}


void nLinAcoustic2::Print(std::ostream * out, const Matrix<Double> Result) const
{
  ENTER_FCN( "nLinAcoustic2::Print");

  (*out)<< "nLinAcoustic2 matrix:" << std::endl << Result;
}

} // end namespace CoupledField
