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

    std::vector<Double> shapeFncAtIp;
    Matrix<Double> partElemMat;
    std::vector<Double> CoordAtIP;

    // set matrix to desired size and set all elements to zero
    //    partElemMat.Resize(nrNodes);
    elemMat.Resize(nrNodes);
    

    for (Integer actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {
	jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord);

	ptelem-> GetShFncAtIp(shapeFncAtIp, actIntPt);

	partElemMat.DyadicMult(shapeFncAtIp);

	if (isaxi_)
	  {
	    CoordAtIP = ptCoord * shapeFncAtIp;
	    partElemMat *= 2 * PI * intWeights[actIntPt-1] * density_ * jacDet * CoordAtIP[0];
	  }
	else 
	  partElemMat *= intWeights[actIntPt-1] * density_ * jacDet;
      
	elemMat += partElemMat;
      }

    if (nrDofsPerNode_ > 1 && dofzero_ > 0 )
      {
	Matrix <Double> multDofMassZero;
	MassMultiDofZero(multDofMassZero, elemMat);
	elemMat = multDofMassZero;
      }
    
    else if (nrDofsPerNode_ > 1 && dofzero_ == 0 )
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
  }


  MassInt::MassInt(BaseFE * aptelem, Double aDensity, Boolean axi)
    : BaseForm(aptelem), 
      density_(aDensity), 
      nrDofsPerNode_(1)
  {
#ifdef TRACE
    (*trace) << "entering MassInt::MassInt" << std::endl;
#endif
    isaxi_ = axi;
    dofzero_ = 0;
  }



  MassInt::MassInt(const Double aDensity,  const Integer nrDofsPerNode, Boolean axi)
    : BaseForm(), 
      density_(aDensity), 
      nrDofsPerNode_(nrDofsPerNode)
  {
#ifdef TRACE
    (*trace) << "entering MassInt::MassInt" << std::endl;
#endif
    isaxi_ = axi;
    dofzero_ = 0;
    
  }


 
  MassInt::~MassInt()
  {
#ifdef TRACE
    (*trace) << "entering MassInt::~MassInt" << std::endl;
#endif
  }



  void MassInt::MassMultiDof(Matrix<Double>& massMultDof, 
			     const Matrix<Double>& massMatSingleDof,  const Integer nrDofs)
  {
    
#ifdef TRACE
    (*trace) << "entering MassInt::MassMultiDof" << std::endl;
#endif
    
    const Integer singleDofSize = massMatSingleDof.GetSizeRow();

    Integer i, j, actDof;
    
    massMultDof.Resize(nrDofsPerNode_*singleDofSize);
    massMultDof.Init();
    
    for (i=0; i < singleDofSize; i++)
      for (j=0; j < singleDofSize; j++)
	for (actDof=0; actDof < nrDofs; actDof++)
	  massMultDof[i*nrDofs + actDof][j*nrDofs + actDof] = massMatSingleDof[i][j]; 
}


  MassInt::MassInt(const Double aDensity,  const Integer nrDofsPerNode, Integer dofzero,
                   Boolean axi)
    : BaseForm(), 
      density_(aDensity), 
      nrDofsPerNode_(nrDofsPerNode)
  {
#ifdef TRACE
    (*trace) << "entering MassInt::MassInt" << std::endl;
#endif
    isaxi_ = axi;
    dofzero_ = dofzero;
  }

  void MassInt::MassMultiDofZero(Matrix<Double>& massMultDofZero, const Matrix<Double>& massMatSingleDof)
  {
#ifdef TRACE
    (*trace) << "entering MassInt::MassMultiDofZero" << std::endl;
#endif
    
    Integer nrDofs = nrDofsPerNode_ -1;
    const Integer singleDofSize = massMatSingleDof.GetSizeRow();
    
    Integer i, j, actDof;
    
    massMultDofZero.Resize(nrDofsPerNode_*singleDofSize);
    massMultDofZero.Init();
    
    for (i=0; i < singleDofSize; i++)
      for (j=0; j < singleDofSize; j++)
	for (actDof=0; actDof < nrDofsPerNode_ ; actDof++)
	  {
	    if (actDof+1 != dofzero_)
	      massMultDofZero[i*nrDofsPerNode_ + actDof][j*nrDofsPerNode_ + actDof] = 
		massMatSingleDof[i][j]; 
	  }
}

} // end namespace CoupledField



