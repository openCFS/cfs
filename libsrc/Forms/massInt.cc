#include <iostream>
#include <fstream>

#include "massInt.hh"

namespace CoupledField {

  MassInt::MassInt(BaseFE * aptelem, Double aDensity, Boolean axi)
    : BaseForm(aptelem), 
      density_(aDensity), 
      nrDofsPerNode_(1)
  {
    ENTER_FCN( "MassInt::MassInt" );

    factor_ = 1.0;
    isaxi_ = axi;
    dofzero_ = 0;
    baseType_ = MASS;
  }

  MassInt::MassInt(const Double aDensity,  const UInt nrDofsPerNode, Boolean axi)
    : BaseForm(), 
      density_(aDensity), 
      nrDofsPerNode_(nrDofsPerNode)
  {
    ENTER_FCN( "MassInt::MassInt" );

    factor_ = 1.0;
    isaxi_ = axi;
    dofzero_ = 0;
    baseType_ = MASS;
  }

  MassInt::MassInt(const Double aDensity,  const UInt nrDofsPerNode, UInt dofzero,
                   Boolean axi)
    : BaseForm(), 
      density_(aDensity), 
      nrDofsPerNode_(nrDofsPerNode)
  {
    ENTER_FCN( "MassInt::MassInt" );

    factor_ = 1.0;
    isaxi_ = axi;
    dofzero_ = dofzero;
    baseType_ = MASS;
  }
 
  MassInt::~MassInt()
  {
    ENTER_FCN( "MassInt::~MassInt" );
  }

  void MassInt::CalcElementMatrix(Matrix<Double> & ptCoord,
                                  Matrix<Double> & elemMat)
  {
    ENTER_FCN( "MassInt::CalcElemMatrix" );
  
    const UInt nrIntPts= ptelem->GetNumIntPoints();
    const UInt nrNodes = ptelem->GetNumNodes();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Double jacDet;

    Vector<Double> shapeFncAtIp;
    Matrix<Double> partElemMat;
    Vector<Double> CoordAtIP;

    // set matrix to desired size and set all elements to zero
    //    partElemMat.Resize(nrNodes);
    elemMat.Resize(nrNodes);
    

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {

      jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord);
        
      ptelem-> GetShFncAtIp(shapeFncAtIp, actIntPt);
        
      partElemMat.DyadicMult(shapeFncAtIp);
        
      if (isaxi_) {
        CoordAtIP = ptCoord * shapeFncAtIp;
        partElemMat *= 2 * PI * intWeights[actIntPt-1] * density_ * factor_* jacDet * CoordAtIP[0];
      }
      else 
        partElemMat *= intWeights[actIntPt-1] * density_ * factor_ * jacDet;
        
      elemMat += partElemMat;
    }
  
    if (nrDofsPerNode_ > 1 && dofzero_ > 0 ) {
      Matrix <Double> multDofMassZero;
      MassMultiDofZero(multDofMassZero, elemMat);
      elemMat = multDofMassZero;
    }
    
    else if (nrDofsPerNode_ > 1 && dofzero_ == 0 ) {
      Matrix <Double> multDofMass;
      MassMultiDof(multDofMass, elemMat, nrDofsPerNode_);       
      elemMat = multDofMass;
    }

  }

  void MassInt::Print(std::ostream * out, const Matrix<Double> Result) const
  {
    ENTER_FCN( "MassInt::Print" );
  }


  void MassInt::MassMultiDof(Matrix<Double>& massMultDof, 
                             const Matrix<Double>& massMatSingleDof,  const UInt nrDofs)
  {
    ENTER_FCN( "MassInt::MassMultiDof" );
    
    const UInt singleDofSize = massMatSingleDof.GetSizeRow();

    UInt i, j, actDof;
    
    massMultDof.Resize(nrDofsPerNode_*singleDofSize);
    massMultDof.Init();
    
    for (i=0; i < singleDofSize; i++)
      for (j=0; j < singleDofSize; j++)
        for (actDof=0; actDof < nrDofs; actDof++)
          massMultDof[i*nrDofs + actDof][j*nrDofs + actDof] = massMatSingleDof[i][j]; 
  }


  void MassInt::MassMultiDofZero(Matrix<Double>& massMultDofZero, const Matrix<Double>& massMatSingleDof)
  {
    ENTER_FCN( "MassInt::MassMultiDofZero" );
    
    const UInt singleDofSize = massMatSingleDof.GetSizeRow();
    
    UInt i, j, actDof;
    
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



