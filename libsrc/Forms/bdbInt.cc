#include <iostream>
#include <fstream>

#include "bdbInt.hh"

namespace CoupledField
{


  void BDBInt::CalcElementMatrix(Matrix<Double>& ptCoord, Matrix<Double> & elemMat)
  {
#ifdef TRACE
    (*trace) << "entering BDBInt::CalcElementMatrix" << std::endl;
#endif

    const Integer nrIntPts = ptelem->GetNumIntPoints(); 
    const Integer nrNodes  = ptelem->GetNumNodes();   
    const Integer nrDofs   = getNrDofs();  
    const std::vector<Double> & intWeights = ptelem->GetIntWeights();  
    double jacDet;  


    Matrix<Double> bMat; 
    Matrix<Double> dMat; 
    Matrix<Double> dB; 
    Matrix<Double> bTrans; 
    Matrix<Double> partElemMat;
  
    elemMat.Resize(nrNodes * nrDofs);
 
    calcDMat(dMat);

#ifdef DEBUG
    *debug << std::endl << "d-Matrix of BDB integrator: " << std::endl 
	   << dMat << std::endl;
#endif

    for (Integer actIntPt=1; actIntPt<=nrIntPts; actIntPt++)
      {
	calcBMat(bMat, actIntPt, ptCoord);  

#ifdef DEBUG
	*debug << std::endl << "B-Matrix of BDB integrator: at intPt " << actIntPt << std::endl 
	       << "coords: " << ptCoord << std::endl
	       << bMat << std::endl;
#endif

	dB = dMat * bMat;

	bMat.Transpose(bTrans);
	
	partElemMat = bTrans * dB;
	
	jacDet = ptelem->CalcJacobianDetAtIp(actIntPt,ptCoord);

#ifdef DEBUG
	*debug << std::endl << "BDBInt::CalcElementMatrix: partelemmat at intPt " << actIntPt << std::endl 
	       << partElemMat << std::endl
	       << "jacDet " << jacDet << std::endl
	       << "intWeights[actIntPt]: " << intWeights[actIntPt-1] << std::endl;
#endif
	
	elemMat += partElemMat * jacDet * intWeights[actIntPt-1] ;
      }
  }




  BDBInt::BDBInt(BaseFE * aptelem, MaterialData & matData) : BaseForm(aptelem, matData)
  {
#ifdef TRACE
    (*trace) << "entering BDBInt::BDBInt" << std::endl;
#endif
  }

  
  BDBInt::~BDBInt()
  {
#ifdef TRACE
    (*trace) << "entering BDBInt::~BDBInt" << std::endl;
#endif

    ;
  }



} // namespace CoupledField
