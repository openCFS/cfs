#include <iostream>
#include <fstream>

#include "bdbInt.hh"

namespace CoupledField
{


  void BDBInt::CalcElementMatrix(Matrix<Double>& ptCoord, Matrix<Double> & elemMat)
  {
    ENTER_FCN( "BDBInt::CalcElementMatrix" );

    const Integer nrIntPts = ptelem->GetNumIntPoints(); 
    const Integer nrNodes  = ptelem->GetNumNodes();   
    const Integer nrDofs   = getNrDofs();  
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
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
	    Vector<Double> ShpFncAtIp;
	    Vector<Double> CoordAtIP;
	    ptelem->GetShFncAtIp(ShpFncAtIp,actIntPt);

	    CoordAtIP = ptCoord * ShpFncAtIp;
            jacDet *= 2 * PI * CoordAtIP[0];
	  }

	elemMat += partElemMat * jacDet * intWeights[actIntPt-1] ;
      }
  }


  void BDBInt::CalcComplexElementMatrix(Matrix<Double> & ptCoord, Matrix<Complex> & elemMat,Double & beta, Double & omega)
  {
    ENTER_FCN( "BDBInt::CalcComplexElementMatrix" );

    const Integer nrIntPts = ptelem->GetNumIntPoints(); 
     const Integer nrNodes  = ptelem->GetNumNodes();   
     const Integer nrDofs   = getNrDofs();  
     const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    double jacDet;
 
    Matrix<Double> bMat; 
    Matrix<Complex> dMat; 
    Matrix<Complex> dB; 
    Matrix<Double> bTrans; 
    Matrix<Complex> partElemMat;
  
    elemMat.Resize(nrNodes * nrDofs);
    elemMat.Init();
    
    calcDMaterialMatWithComplexDamping(dMat,beta,omega);
	//      calcDMat(dMat);
    
    for (Integer actIntPt=1; actIntPt<=nrIntPts; actIntPt++)
      {
	if (updateDMatInEveryIP_)
	  calcDMaterialMatWithComplexDamping(dMat,beta,omega);
	  //	  calcDMat(dMat, actIntPt, ptCoord);
    
	calcBMat(bMat, actIntPt, ptCoord);
      
	//	std::cout<<"BMat= "<< bMat.GetSizeRow()<< " x " << bMat.GetSizeCol()<< " DMat = "<< dMat.GetSizeRow()<< " x " << dMat.GetSizeCol() <<  std::endl;

	//   hardcoded dB = dMat * bMat;
	dB.Resize(dMat.GetSizeRow(), bMat.GetSizeCol());
	Complex a;

	for (int i = 0; i < dMat.GetSizeRow(); i++)
	  for (int j = 0; j < bMat.GetSizeCol(); j++)
	    {       
	      a = dMat[i][0] * bMat[0][j];
	      for (Integer k = 1; k <bMat.GetSizeRow(); k++)
		a += dMat[i][k] * bMat[k][j];
	      dB[i][j] = a;
	    }  

	bMat.Transpose(bTrans);

	//	hardcoded: partElemMat = bTrans * dB;
	partElemMat.Resize(bTrans.GetSizeRow(), dB.GetSizeCol());
	for (int i = 0; i < bTrans.GetSizeRow(); i++)
	  for (int j = 0; j < dB.GetSizeCol(); j++)
	    {       
	      a = bTrans[i][0] *dB[0][j];
	      for (Integer k = 1; k <dB.GetSizeRow(); k++)
		a += bTrans[i][k] * dB[k][j];
	     partElemMat[i][j] = a;
	    }
        
	jacDet = ptelem->CalcJacobianDetAtIp(actIntPt,ptCoord);

	if (jacDet < 0)
	  Error("Negative Jacobian determinant!", __FILE__, __LINE__);

	if (isaxi_)
	  {
	    Vector<Double> ShpFncAtIp;
	    Vector<Double> CoordAtIP;
	    ptelem->GetShFncAtIp(ShpFncAtIp,actIntPt);

	    CoordAtIP = ptCoord * ShpFncAtIp;
            jacDet *= 2 * PI * CoordAtIP[0];
	  }

	for (int i = 0; i < elemMat.GetSizeRow(); i++)
	  for (int j = 0; j < elemMat.GetSizeCol(); j++)
       	elemMat[i][j] = elemMat[i][j]+ partElemMat[i][j] * jacDet * intWeights[actIntPt-1] ;
      }
  }




  BDBInt::BDBInt(BaseFE * aptelem, MaterialData & matData) 
    : BaseForm(aptelem, matData), updateDMatInEveryIP_(0)
  {
    ENTER_FCN( "BDBInt::BDBInt" );
  }


  BDBInt::BDBInt(MaterialData & matData) 
    : BaseForm(matData), updateDMatInEveryIP_(0)
  {
    ENTER_FCN( "BDBInt::BDBInt" );
  }

  
  BDBInt::~BDBInt()
  {
    ENTER_FCN( "BDBInt::~BDBInt" );
    ;
  }


} // namespace CoupledField
