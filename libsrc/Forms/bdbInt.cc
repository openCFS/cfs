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

    for (Integer actIntPt=1; actIntPt<=nrIntPts; actIntPt++)
      {
	calcBMat(bMat, actIntPt, ptCoord);  

	dB = dMat * bMat;

	bMat.Transpose(bTrans);
	
	partElemMat = bTrans * dB;
	
	jacDet = ptelem->CalcJacobianDetAtIp(actIntPt,ptCoord);
	
	elemMat += partElemMat * jacDet * intWeights[actIntPt] ;
      }
  }




  BDBInt::BDBInt(BaseFE * aptelem) : BaseForm(aptelem)
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
