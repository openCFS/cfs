#include <iostream>
#include <fstream>
#include <math.h>

#include "linearForm.hh"


namespace CoupledField
{


  LinearForm::LinearForm(BaseFE * aptelem) : BaseForm(aptelem)
  {
#ifdef TRACE
    (*trace) << "entering LinearForm::LinearForm" << std::endl;
#endif
  }


  
  LinearForm ::~LinearForm()
  {
#ifdef TRACE
    (*trace) << "entering LinearForm::~LinearForm" << std::endl;
#endif
  }



  void LinearForm::CalcElemVector(Matrix<Double>& ptCoord, std::vector<Double> & Result)
  {
#ifdef TRACE
    (*trace) << "entering LinearForm::CalcElemVector" << std::endl;
#endif
  }
  


  LinearEdgeInt::LinearEdgeInt(BaseFE * aptelem, Double aVal, Integer aDirection) 
    : LinearForm(aptelem), val_(aVal), direction_(aDirection)
  {
#ifdef TRACE
    (*trace) << "entering LinearEdgeInt::LinearEdgeInt" << std::endl;
#endif
  }


  
  LinearEdgeInt ::~LinearEdgeInt()
  {
#ifdef TRACE
    (*trace) << "entering LinearEdgeInt::~LinearEdgeInt" << std::endl;
#endif
  }



  void LinearEdgeInt::CalcElemVector(Matrix<Double>& ptCoord, std::vector<Double> & elemVec)
  {
#ifdef TRACE
    (*trace) << "entering LinearEdgeInt::CalcElemVector" << std::endl;
#endif

    const Integer nrIntPts = ptelem->GetNumIntPoints();
    const Integer nrEdges  = ptelem->GetNumEdges();
    const Integer dim      = ptelem->GetDim();
    const std::vector<Double> & intWeights = ptelem->GetIntWeights();  
    Double jacDet;  


    // derivation of shape functions after global coordinates 
    Matrix<Double> shapeEdge;
    std::vector<Double> partElemVec;
    partElemVec.resize(nrEdges);


    // set vector to desired size and set all elements to zero    
    elemVec.resize(nrEdges); 
    for (Integer i=0; i<nrEdges; i++)
      elemVec[i] = 0;
    
    std::vector<Double> currentVec(dim, 0);
    
    // define direction of current
    currentVec[direction_-1] = val_;	

    for (Integer actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {
	ptelem->CalcEdgeShapeFncAtIp(shapeEdge, actIntPt, ptCoord);

	partElemVec = shapeEdge * currentVec;
	
	jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord);

	for(Integer i=0; i<partElemVec.size(); i++)
	  elemVec[i] += partElemVec[i] * intWeights[actIntPt-1] * jacDet;
      }
  

#ifdef DEBUG 
	(*debug) << "CalcElemVector:  "  << std::endl
		 << partElemVec << std::endl
		 << "\n jacDet " << jacDet << std::endl;
#endif


#ifdef TRACE
    (*trace) << "leaving linearEdgeInt::CalcElementVector" << std::endl;
#endif

  }
  



  } // end of namespace
