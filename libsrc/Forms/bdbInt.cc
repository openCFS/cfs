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
 
    if (!updateDMatInEveryIP_)
      calcDMat(dMat);
    

    for (Integer actIntPt=1; actIntPt<=nrIntPts; actIntPt++)
      {
	if (updateDMatInEveryIP_)
	  calcDMat(dMat, actIntPt, ptCoord);
    
	calcBMat(bMat, actIntPt, ptCoord);

	dB = dMat * bMat;

	bMat.Transpose(bTrans);

	partElemMat = bTrans * dB;

	jacDet = ptelem->CalcJacobianDetAtIp(actIntPt,ptCoord);

	if (jacDet < 0)
	  Error("Negative Jacobian determinant!", __FILE__, __LINE__);

	if (isaxi_)
	  {
	    std::vector<Double> ShpFncAtIp;
	    std::vector<Double> CoordAtIP;
	    ptelem->GetShFncAtIp(ShpFncAtIp,actIntPt);
	    CoordAtIP = ptCoord * ShpFncAtIp;
            jacDet *= 2 * PI * CoordAtIP[0];
	  }

	elemMat += partElemMat * jacDet * intWeights[actIntPt-1] ;
      }
  }




  BDBInt::BDBInt(BaseFE * aptelem, MaterialData & matData) 
    : BaseForm(aptelem, matData), updateDMatInEveryIP_(0)
  {
#ifdef TRACE
    (*trace) << "entering BDBInt::BDBInt" << std::endl;
#endif
  }


  BDBInt::BDBInt(MaterialData & matData) 
    : BaseForm(matData), updateDMatInEveryIP_(0)
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
