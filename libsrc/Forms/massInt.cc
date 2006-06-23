#include <iostream>
#include <fstream>

#include "massInt.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"

namespace CoupledField {



  MassInt::MassInt(const Double aDensity,  const UInt nrDofsPerNode, 
                   bool axi, bool coordUpdate )
    : BaseForm(NULL), 
      density_(aDensity), 
      nrDofsPerNode_(nrDofsPerNode)
      
  {
    ENTER_FCN( "MassInt::MassInt" );
    name_ = "MassInt";
    factor_ = 1.0;
    isaxi_ = axi;
    coordUpdate_ = coordUpdate;
    baseType_ = MASS;
  }

 
  MassInt::~MassInt()
  {
    ENTER_FCN( "MassInt::~MassInt" );
  }

 

  void MassInt::CalcElementMatrix( Matrix<Double>& elemMat,
                                   EntityIterator& ent1, 
                                   EntityIterator& ent2  )  {
    ENTER_FCN( "MassInt::CalcElemMatrix" );
  
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );
    
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
    elemMat.Init();
    

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {

      jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_);
        
      ptelem-> GetShFncAtIp(shapeFncAtIp, actIntPt);
        
      partElemMat.DyadicMult(shapeFncAtIp);
        
      if (isaxi_) {
        CoordAtIP = ptCoord_ * shapeFncAtIp;
        partElemMat *= 2 * PI * intWeights[actIntPt-1] * density_ * factor_* jacDet * CoordAtIP[0];
      }
      else 
        partElemMat *= intWeights[actIntPt-1] * density_ * factor_ * jacDet;
        
      elemMat += partElemMat;
    }
  
    if (nrDofsPerNode_ > 1 ) {
      Matrix <Double> multDofMass;
      MassMultiDof(multDofMass, elemMat, nrDofsPerNode_);       
      elemMat = multDofMass;
    }

    //    std::cout << "ElemMatMass:\n" << elemMat << std::endl;
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
            massMultDofZero[i*nrDofsPerNode_ + actDof][j*nrDofsPerNode_ + actDof] = 
              massMatSingleDof[i][j]; 
          }
  }
  
} // end namespace CoupledField



