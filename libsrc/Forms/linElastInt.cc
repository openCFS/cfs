#include <iostream>
#include <fstream>

#include "linElastInt.hh"

namespace CoupledField
{

  // returns B - matrix for BDB
  void linElastInt::calcBMat(Matrix<Double> & bMat, Integer ip, Matrix<Double> & ptCoord)
  {
#ifdef TRACE
  (*trace) << "entering linElastInt::calcBMat " << std::endl;
#endif
    const Integer nrNodes  = ptelem->GetNumNodes();
    const Integer spaceDim = ptelem->GetDim();  
    const Integer nrDofs   = getNrDofs();  

    Integer actDim, actNode, help, j, k;
    
    
    bMat.Resize(getDimD(), nrNodes * nrDofs);
    
    // local shape functions derived after global coords (format: nrNodes x spaceDim)
    Matrix<Double> xiDx;

    ptelem->GetGlobDerivShFncAtIp(xiDx, ip, ptCoord);


    for(actDim=0; actDim < spaceDim; actDim++)
      for(actNode=0; actNode < nrNodes; actNode++)
	bMat[actDim][actNode * spaceDim + actDim] = xiDx[actNode][actDim];

    switch(spaceDim)
      {
      case 2:
	j = 1;
	k = 0;
	
	for (actNode = 0; actNode < nrNodes; actNode++)
	  {
	    bMat[spaceDim][actNode * spaceDim + 1] = xiDx[actNode][0];
	    bMat[spaceDim][actNode * spaceDim]     = xiDx[actNode][1];
	  }
	break;


      case 3:
	Integer actDim=spaceDim;
	for (actNode = 0; actNode < nrNodes; actNode++)
	  {
	    bMat[actDim][actNode * spaceDim + 1] = xiDx[actNode][2];
	    bMat[actDim][actNode * spaceDim + 2] = xiDx[actNode][1];
	  }

	actDim++;
	for (actNode = 0; actNode < nrNodes; actNode++)
	  {
	    bMat[actDim][actNode * spaceDim]     = xiDx[actNode][2];
	    bMat[actDim][actNode * spaceDim + 2] = xiDx[actNode][0];
	  }

	actDim++;
	for (actNode = 0; actNode < nrNodes; actNode++)
	  {
	    bMat[actDim][actNode * spaceDim]     = xiDx[actNode][1];
	    bMat[actDim][actNode * spaceDim + 1] = xiDx[actNode][0];
	  }
	break;
      }
  }
  
  
  void mechPlainStrainInt::calcDMat(Matrix<Double> & dMat)
  {
#ifdef TRACE
  (*trace) << "entering linElastInt::calcDMat " << std::endl;
#endif
    CalcPlainStrainMat(dMat);
  }
  

  // a 2d-problem is calculated in the x-y-plane
  void mechPlainStrainInt::CalcPlainStrainMat(Matrix<Double> & matDat)
  {
    const Integer nrElems2d = getDimD();
    
    Integer rowPtrXY[] = {1,2,6,7,8};  // indices of rows and lines for xy-plane
    Integer rowPtrYZ[] = {2,3,4,8,9};  // indices of rows and lines for yz-plane
    Integer rowPtrXZ[] = {1,3,5,7,9};  // indices of rows and lines for xz-plane
    Integer * rowPtr;
    Integer i,j;

    switch(actOrientation)
      {	
      case xy: 
	{
	  rowPtr = rowPtrXY;
	  break;
	}
      case xz: 
	{
	  rowPtr = rowPtrXZ;
	  break;
	}

      case yz: 
	{
	  rowPtr = rowPtrYZ;    
	  break;
	}
      }    
	
    Matrix<Double> * matMatrix =  ptMaterial->GetMatrix();
    
    matDat.Resize(nrElems2d);

    for (i=0; i<nrElems2d; i++)
      for (j=0; j<nrElems2d; j++)
	matDat[i][j] = (*matMatrix)[rowPtr[i]-1][rowPtr[j]-1];	
}




  // ===================================================================================
  // =================== standard con- and destructors (just for tracing) ==============
  // ===================================================================================


  linElastInt::linElastInt(BaseFE * aptelem, MaterialData & matData) 
    : BDBInt(aptelem, matData)
  {
#ifdef TRACE
    (*trace) << "entering linElastInt::linElastInt" << std::endl;
#endif

    ptelem=aptelem;
  }
 

  linElastInt::~linElastInt()
  {
#ifdef TRACE
    (*trace) << "entering linElastInt::~linElastInt" << std::endl;
#endif
  }



  mechPlainStrainInt::mechPlainStrainInt(BaseFE * aptelem, MaterialData & matData) 
    : linElastInt(aptelem, matData), actOrientation(xy)
  {
#ifdef TRACE
    (*trace) << "entering mechPlainStrainInt::mechPlainStrainInt" << std::endl;
#endif

    ptelem=aptelem;
  }
 

  mechPlainStrainInt::~mechPlainStrainInt()
  {
#ifdef TRACE
    (*trace) << "entering mechPlainStrainInt::~mechPlainStrainInt" << std::endl;
#endif
  }

}
