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
    elemMat.Init();
 
    calcDMat(dMat);

#ifdef DEBUG
    *debug << std::endl << "d-Matrix of BDB integrator: " << std::endl 
	   << dMat << std::endl;
#endif

    for (Integer actIntPt=1; actIntPt<=nrIntPts; actIntPt++)
      {
	calcBMat(bMat, actIntPt, ptCoord);  

	*debug << std::endl << "b-Matrix of BDB integrator: " << std::endl 
	   << bMat << std::endl;


	dB = dMat * bMat;

	*debug << "b * d =  " << std::endl 
	   << dB << std::endl;

	bMat.Transpose(bTrans);

	partElemMat = bTrans * dB;

	*debug << "bTranspMalDB = [ " << std::endl 
	   << partElemMat << std::endl;
	
	
	jacDet = ptelem->CalcJacobianDetAtIp(actIntPt,ptCoord);

	elemMat += partElemMat * jacDet * intWeights[actIntPt-1] ;

	*debug << "jacDet =  " << jacDet << std::endl
	       << "intWeights " << intWeights[actIntPt-1] << std::endl
	       << "elemMat " << std::endl << elemMat << std::endl;
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
