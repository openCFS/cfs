#include <iostream>
#include <fstream>

#include "Utils/SmoothSpline.hh"
#include "nLinAcoustic1.hh"

namespace CoupledField {

nLinAcoustic1::nLinAcoustic1(Boolean axi)
   : BaseForm()
{
   ENTER_FCN( "nLinAcoustic1::nLinAcoustic1");

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


   // derivation of shape functions after global coordinates 
   Matrix<Double> xiDx;
   Matrix<Double> xiDxTransp;
   Matrix<Double> partElemMat;
   Matrix<Double> partElemMatAxi;
   Vector<Double> ShpFncAtIp;
   Vector<Double> CoordAtIP;
   Vector<Double> drAtIp;

   // set matrix to desired size and set all elements to zero
   elemMat.Resize(nrNodes); elemMat.Init(0);

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
