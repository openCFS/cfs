#include <iostream>
#include <fstream>

#include "massInt.hh"

namespace CoupledField
{


  void MassInt::CalcElementMatrix(Matrix<Double> & ptCoord, Matrix<Double> & elemMat)
  {
#ifdef TRACE
    (*trace) << "entering MassInt::CalcElemMatrix" << std::endl;
#endif
  
    const Integer nrIntPts= ptelem->GetNumIntPoints();
    const Integer nrNodes = ptelem->GetNumNodes();
    const std::vector<Double> & intWeights = ptelem->GetIntWeights();  
    Double jacDet;
    const Double density = ptMaterial->GetDensity();
    

    // derivation of shape functions after global coordinates 

    std::vector<Double> shapeFnc;
    Matrix<Double> partElemMat;


    // set matrix to desired size and set all elements to zero
    //    partElemMat.Resize(nrNodes);
    elemMat.Resize(nrNodes);
    

    for (Integer actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {
	jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord);

	ptelem-> GetShFncAtIp(shapeFnc, actIntPt);

	partElemMat.DyadicMult(shapeFnc);

	partElemMat *= intWeights[actIntPt-1] * density * jacDet;
      
	elemMat += partElemMat;
      }
  

#ifdef TRACE
    (*trace) << "leaving MassInt::CalcElemMatrix" << std::endl;
#endif
  }






  void MassInt::Print(std::ostream * out, const Matrix<Double> Result) const
  {
#ifdef TRACE
    (*trace) << "entering MassInt::Print" << std::endl;
#endif
    (*out)<< "Mass matrix:" << std::endl << Result;
  }


  MassInt::MassInt(BaseFE * aptelem, MaterialData & aMat): BaseForm(aptelem, aMat)
  {
#ifdef TRACE
    (*trace) << "entering MassInt::MassInt" << std::endl;
#endif
  }


 
  MassInt::~MassInt()
  {
#ifdef TRACE
    (*trace) << "entering MassInt::~MassInt" << std::endl;
#endif
  }





} // end namespace CoupledField



