#include <iostream>
#include <fstream>

#include "massEdgeInt.hh"

namespace CoupledField
{


  void MassEdgeInt::CalcElementMatrix(Matrix<Double> & ptCoord, Matrix<Double> & elemMat)
  {
#ifdef TRACE
    (*trace) << "entering MassEdgeInt::CalcElemMatrix" << std::endl;
#endif
  
    const Integer nrIntPts = ptelem->GetNumIntPoints();
    const Integer nrNodes  = ptelem->GetNumNodes();
    const Integer nrEdges  = ptelem->GetNumEdges();
    const std::vector<Double> & intWeights = ptelem->GetIntWeights();  
    Double jacDet;

    // derivation of shape functions after global coordinates 

    Matrix<Double> shapeEdge;
    Matrix<Double> shapeEdgeTransp;
    Matrix<Double> partElemMat;
    //    partElemMat.Resize(nrNodes);

    // set matrix to desired size and set all elements to zero
     elemMat.Resize(nrEdges);
    

    for (Integer actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {
	jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord);

	ptelem->CalcEdgeShapeFncAtIp(shapeEdge, actIntPt, ptCoord);

	shapeEdge.Transpose(shapeEdgeTransp);

	partElemMat = shapeEdge * shapeEdgeTransp;

	partElemMat *= intWeights[actIntPt-1] * conductivity_ * jacDet;
      
	elemMat += partElemMat;
      }
  

#ifdef TRACE
    (*trace) << "leaving MassEdgeInt::CalcElemMatrix" << std::endl;
#endif
  }


  void MassEdgeInt::Print(std::ostream * out, const Matrix<Double> Result) const
  {
#ifdef TRACE
    (*trace) << "entering MassEdgeInt::Print" << std::endl;
#endif
    (*out)<< "Mass matrix:" << std::endl << Result;
  }


  MassEdgeInt::MassEdgeInt(BaseFE * aptelem, Double acond)
    : BaseForm(aptelem), conductivity_(acond)
  {
#ifdef TRACE
    (*trace) << "entering MassEdgeInt::MassEdgeInt" << std::endl;
#endif
  }


 
  MassEdgeInt::~MassEdgeInt()
  {
#ifdef TRACE
    (*trace) << "entering MassEdgeInt::~MassEdgeInt" << std::endl;
#endif
  }





} // end namespace CoupledField



