#include <iostream>
#include <fstream>

#include "bdbInt.hh"

namespace CoupledField
{

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

  void BDBInt::CalcElementMatrix(Matrix<Double>& ptCoord, Matrix<Double> & elemMat)
  {

    Integer nrIntPts = ptelem->GetNumIntPoints(); // l - number of integration points
    Integer nrNodes  = ptelem->GetNumNodes();   // n - number of nodes

    Matrix<Double> bMat; 
    Matrix<Double> dMat; 
    Matrix<Double> dB; 
    Matrix<Double> bTrans; 
    Matrix<Double> partElemMat; 

    Integer const spaceDim = 2;  // has to be read with an apropriate method on the future !!!!!

    std::vector<Double> intWeights=ptelem->GetIntWeights();

  

  
    elemMat.Resize(nrNodes*spaceDim);
 
    getBMat(bMat);  
    getDMat(dMat);

    for (int actIntPt=1; actIntPt<=nrIntPts; actIntPt++)
      {
	dB = dMat * bMat;
	bMat.Transpose(bTrans);
	
	partElemMat = bTrans * dB;
	
	elemMat += partElemMat * ptelem->CalcJacobianDetAtIp(actIntPt,ptCoord) * intWeights[actIntPt] ;
      }
  }
}
