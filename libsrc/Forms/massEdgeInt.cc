#include <iostream>
#include <fstream>

#include "massEdgeInt.hh"

namespace CoupledField
{


  void MassEdgeInt::CalcElementMatrix(Matrix<Double> & ptCoord, Matrix<Double> & elemMat)
  {
    ENTER_FCN( "MassEdgeInt::CalcElemMatrix" );
    const UInt nrIntPts = ptelem->GetNumIntPoints();
    const UInt nrEdges  = ptelem->GetNumEdges();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Double jacDet;
  
    // derivation of shape functions after global coordinates 
  
    Matrix<Double> shapeEdge;
    Matrix<Double> shapeEdgeTransp;
    Matrix<Double> partElemMat;
    //    partElemMat.Resize(nrNodes);
  
    // set matrix to desired size and set all elements to zero
    elemMat.Resize(nrEdges);
  
  
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {
        jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord);
      
        ptelem->CalcEdgeShapeFncAtIp(shapeEdge, actIntPt, ptCoord);
      
        shapeEdge.Transpose(shapeEdgeTransp);
      
        partElemMat = shapeEdge * shapeEdgeTransp;
      
        partElemMat *= intWeights[actIntPt-1] * conductivity_ * jacDet;
      
        elemMat += partElemMat;
      }
  
  }


  void MassEdgeInt::Print(std::ostream * out, const Matrix<Double> Result) const
  {
    ENTER_FCN( "MassEdgeInt::Print" );
    (*out)<< "Mass matrix:" << std::endl << Result;
  }


  MassEdgeInt::MassEdgeInt(BaseFE * aptelem, Double acond)
    : BaseForm(aptelem), conductivity_(acond)
  {
    ENTER_FCN( "MassEdgeInt::MassEdgeInt" );
  }


 
  MassEdgeInt::~MassEdgeInt()
  {
    ENTER_FCN( "MassEdgeInt::~MassEdgeInt" );
  }

} // end namespace CoupledField



