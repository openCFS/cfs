// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include "bdInt.hh"

namespace CoupledField
{


  void BDInt::calcElementVector(Vector<Double> & resultStressVector,
                                EntityIterator& ent,
                                Vector<Double> & fracDerivStress)  {
    
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );

    const Integer nrIntPts = ptelem->GetNumIntPoints();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();
    double jacDet;

    Matrix<Double> bMat;
    Matrix<Double> bTrans;
    Matrix<Double> aMat;    
    Matrix<Double> alphaMat;
    Vector<Double> temp;

    temp.Resize(getDim());

    calcAMat(aMat);
    calcAlphaMat(alphaMat);

    fracDerivStress = alphaMat * fracDerivStress;
    fracDerivStress = aMat * fracDerivStress;


    for (Integer actIntPt=1; actIntPt<=nrIntPts; actIntPt++)
      {
	calcBMat(bMat, actIntPt, ptCoord_);	
	bMat.Transpose(bTrans);
	
	temp = bTrans * fracDerivStress;

	jacDet = ptelem->CalcJacobianDetAtIp(actIntPt,ptCoord_,
                                             ent.GetElem());

	if (jacDet < 0)
	  Error("Negative Jacobian determinant!", __FILE__, __LINE__);

	if (isaxi_)
	  {
	    Vector<Double> ShpFncAtIp;
	    Vector<Double> CoordAtIP;
	    ptelem->GetShFncAtIp(ShpFncAtIp,actIntPt, it1_.GetElem() );

	    CoordAtIP = ptCoord_ * ShpFncAtIp;
            jacDet *= 2 * PI * CoordAtIP[0];
	  }	

	resultStressVector +=  temp  * jacDet * intWeights[actIntPt-1];
      }
  }

  void BDInt::calcBMat(Matrix<Double> & bMat, Integer ip, Matrix<Double> & ptCoord_)
  {

    ptelem->SetAnsatzFct( ansatzFct1_ );
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const UInt spaceDim = ptelem->GetDim();  
    const UInt nrDofs   = getNrDofs();  
   
    
    UInt actDim, actNode, j, k;
    
    
    bMat.Resize(getDim(), numFncs * nrDofs);
    bMat.Init();
    
    // local shape functions derived after global coords (format: numFncs x spaceDim)
    Matrix<Double> xiDx;

    if (isSetIntPoint_) 
      ptelem->GetGlobDerivShFnc(xiDx, intPoint_, ptCoord_, it1_.GetElem() );
    else
      ptelem->GetGlobDerivShFncAtIp(xiDx, ip, ptCoord_, it1_.GetElem() );

    for(actDim=0; actDim < spaceDim; actDim++)
      for(actNode=0; actNode < numFncs; actNode++)
	bMat[actDim][actNode * spaceDim + actDim] = xiDx[actNode][actDim];

    switch(spaceDim)
      {
      case 2:
	j = 1;
	k = 0;
	
	for (actNode = 0; actNode < numFncs; actNode++)
	  {
	    bMat[spaceDim][actNode * spaceDim + 1] = xiDx[actNode][0];
	    bMat[spaceDim][actNode * spaceDim]     = xiDx[actNode][1];
	  }

	if (isaxi_)
	  {
	    Integer idxtheta = getDimD();
	    Vector<Double> ShpFncAtIp;
	    Vector<Double> CoordAtIP;

	    if (isSetIntPoint_) 
	      ptelem->GetShFnc(ShpFncAtIp,intPoint_, it1_.GetElem() );
	    else
	      ptelem->GetShFncAtIp(ShpFncAtIp, ip, it1_.GetElem() );

	    CoordAtIP = ptCoord_ * ShpFncAtIp;

	    for (actNode = 0; actNode < numFncs; actNode++)	     
	      bMat[idxtheta-1][actNode * spaceDim] = ShpFncAtIp[actNode] / CoordAtIP[0];
	  }

	break;

      case 3:
	Integer actDim=spaceDim;
	for (actNode = 0; actNode < numFncs; actNode++)
	  {
	    bMat[actDim][actNode * spaceDim + 1] = xiDx[actNode][2];
	    bMat[actDim][actNode * spaceDim + 2] = xiDx[actNode][1];
	  }

	actDim++;
	for (actNode = 0; actNode < numFncs; actNode++)
	  {
	    bMat[actDim][actNode * spaceDim]     = xiDx[actNode][2];
	    bMat[actDim][actNode * spaceDim + 2] = xiDx[actNode][0];
	  }

	actDim++;
	for (actNode = 0; actNode < numFncs; actNode++)
	  {
	    bMat[actDim][actNode * spaceDim]     = xiDx[actNode][1];
	    bMat[actDim][actNode * spaceDim + 1] = xiDx[actNode][0];
	  }
	break;
      }

    isSetIntPoint_ = false;
  }

  void BDInt::calcAlphaMat(Matrix<Double> & aMat)
  {    
    double val = 0.0;

    // compute matrix A,  same entries, 
    aMat.Resize(getDim());
    aMat.Init();
    
    // set entries on the diagonal
    val = timeStepPowerFracDeriv_ * dampAlpha_ ;
	
    // set the emtries on the diagonal matrix, to get the inverse, 
    // the value on the diagonal are 1/val
    for(int i=0;i<getDim();i++)
      {
	aMat.SetEntry(i,i,val);
      }  
  }


  void BDInt::calcAMat(Matrix<Double> & aMat)
  {    
    
    double val = 0.0;
    // compute matrix A,  same entries, 
    aMat.Resize(getDim());
    aMat.Init();
    
    // set entries on the diagonal
    val = timeStepPowerFracDeriv_ * dampAlpha_ ;
	
    // set the emtries on the diagonal matrix, to get the inverse, 
    // the value on the diagonal are 1/val
    for(int i=0;i<getDim();i++)
      {
	aMat.SetEntry(i,i,1/(val+1));
      }  
  }

 Integer BDInt::getDim()
  {
    
    if(geomType_ == "axi")
      {
    	return 4;
      }
    else if(geomType_ == "planeStrain")
      {
    	return 3;   	
      }
    else if(geomType_ == "3d")
      {
    	return 6;   	
      }
    else
      {
	Error("wrongh geomType(axi,planeStrein,3d) specified", __FILE__, __LINE__);  
      }
    return 0;
  }

  /// returns nr. of degrees of freedom
  Integer BDInt::getNrDofs()
  {
   
    if(geomType_ == "axi" || geomType_ == "planeStrain")
      {
    	return 2;
      }
    else if(geomType_ == "3d")
      {
    	return 3;   	
      }
    else 
      {
	Error("wrongh geomType(axi,planeStrein,3d) specified", __FILE__, __LINE__);
	return -1;
      }    
  }


  BDInt::BDInt(BaseMaterial* matData,std::string geomType, Double timeStep)
    : BaseForm(matData), updateDMatInEveryIP_(0) {

    *warning << "This class is a candidate for major reimplementation and therefore "
             << "not working -> Contact Gerhard and Andreas";
    Warning( __FILE__, __LINE__ );
  //   name_ = "BDInt";
//     geomType_ = geomType;
//     ptMaterial->GetScalar(dampAlpha_,RAYLEIGH_ALPHA,REAL);

//    StdVector<Double> fracDerivList_;
//    params->GetList( "fracDeriv", fracDerivList_, "mechanic", "damping" );
//    Double fracDeriv_ = fracDerivList_[0];
//    timeStepPowerFracDeriv_ = std::pow(timeStep,-fracDeriv_);
  }


  BDInt::~BDInt()
  {
  }


} // namespace CoupledField
