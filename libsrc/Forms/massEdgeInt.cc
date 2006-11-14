#include <iostream>
#include <fstream>

#include "massEdgeInt.hh"

namespace CoupledField
{


  void MassEdgeInt::CalcElementMatrix( Matrix<Double>& elemMat,
                                       EntityIterator& ent1, 
                                       EntityIterator& ent2 ) {
    ENTER_FCN( "MassEdgeInt::CalcElemMatrix" );

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );

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
    elemMat.Init();
  
  
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {
        jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, 
                                             ent1.GetElem() );
      
        ptelem->CalcEdgeShapeFncAtIp(shapeEdge, actIntPt, ptCoord_);
      
        shapeEdge.Transpose(shapeEdgeTransp);
      
        partElemMat = shapeEdge * shapeEdgeTransp;
      
        partElemMat *= intWeights[actIntPt-1] * conductivity_ * jacDet;
      
        elemMat += partElemMat;
      }
  
  }



  MassEdgeInt::MassEdgeInt( Double acond )
    : BaseForm( NULL ), conductivity_(acond)
  {
    ENTER_FCN( "MassEdgeInt::MassEdgeInt" );

    name_ = "MassEdgeInt";
  }


 
  MassEdgeInt::~MassEdgeInt()
  {
    ENTER_FCN( "MassEdgeInt::~MassEdgeInt" );
  }

} // end namespace CoupledField



