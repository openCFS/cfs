#include <iostream>
#include <fstream>

#include "massInt.hh"

namespace CoupledField
{


  void MassInt::CalcElementMatrix(Matrix<Double> & ptCoord,
				  Matrix<Double> & elemMat)
  {
#ifdef TRACE
    (*trace) << "entering MassInt::CalcElemMatrix" << std::endl;
#endif
  
    const Integer nrIntPts= ptelem->GetNumIntPoints();
    const Integer nrNodes = ptelem->GetNumNodes();
    const std::vector<Double> & intWeights = ptelem->GetIntWeights();  
    Double jacDet;

    // derivation of shape functions after global coordinates 

    std::vector<Double> shapeFnc;
    Matrix<Double> partElemMat;
    std::vector<Double> ShpFncAtIp;
    std::vector<Double> CoordAtIP;

    // set matrix to desired size and set all elements to zero
    //    partElemMat.Resize(nrNodes);
    elemMat.Resize(nrNodes);
    

    for (Integer actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {
	jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord);

	ptelem-> GetShFncAtIp(shapeFnc, actIntPt);

	partElemMat.DyadicMult(shapeFnc);

	if (isaxi_)
	  {
	    ptelem->GetShFncAtIp(ShpFncAtIp,actIntPt);
	    CoordAtIP = ptCoord * ShpFncAtIp;
	    partElemMat *= 2 * PI * intWeights[actIntPt-1] * density_ * jacDet * CoordAtIP[0];
	  }
	else 
	  partElemMat *= intWeights[actIntPt-1] * density_ * jacDet;
      
	elemMat += partElemMat;
      }

    if (nrDofsPerNode_ > 1)
      {
	Matrix <Double> multDofMass;
	MassMultiDof(multDofMass, elemMat, nrDofsPerNode_);
	elemMat = multDofMass;
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


  MassInt::MassInt(BaseFE * aptelem, Double aDensity, Boolean axi)
    : BaseForm(aptelem), 
      density_(aDensity), 
      isaxi_(axi),
      nrDofsPerNode_(1)
  {
#ifdef TRACE
    (*trace) << "entering MassInt::MassInt" << std::endl;
#endif
  }


  MassInt::MassInt(const Double aDensity,  const Integer nrDofsPerNode, Boolean axi)
    : BaseForm(), 
      density_(aDensity), 
      isaxi_(axi), 
      nrDofsPerNode_(1)
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



  void MassInt::MassMultiDof(Matrix<Double>& massMultDof, const Matrix<Double>& massMatSingleDof,  const Integer nrDofs)
  {
    
#ifdef TRACE
    (*trace) << "entering MassInt::MassMultiDof" << std::endl;
#endif
    
    const Integer singleDofSize = massMatSingleDof.getSize();
    const Integer multDofSize   = singleDofSize * nrDofs;
    
    Integer i, j, actDof;
    
    massMultDof.Resize(multDofSize);
    massMultDof.Init();
    
    for (i=0; i < singleDofSize; i++)
      for (j=0; j < singleDofSize; j++)
	for (actDof=0; actDof < nrDofs; actDof++)
	  massMultDof[i*nrDofs + actDof][j*nrDofs + actDof] = massMatSingleDof[i][j]; 
}


} // end namespace CoupledField



