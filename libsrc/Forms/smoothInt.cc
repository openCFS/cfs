#include <iostream>
#include <fstream>

#include "smoothInt.hh"

namespace CoupledField
{



  // returns B - matrix for BDB
  void SmoothInt::calcBMat(Matrix<Double> & bMat, Integer ip, Matrix<Double> & ptCoord)
  {
    ENTER_FCN( "SmoothInt::calcBMat" );
    
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
  

  
  // calculated the D-matrix for the plain strain state
  void smoothPlainStrainInt::calcDMat(Matrix<Double> & dMat, 
				      Integer ip, 
				      Matrix<Double> & ptCoord)
  {
    ENTER_FCN( "smoothPlainStrainInt::calcDMat" );

    const Integer nrElems2d = getDimD();
    
    Integer rowPtrXY[] = {1,2,6,7,8};  // indices of rows and lines for xy-plane
    Integer rowPtrYZ[] = {2,3,4,8,9};  // indices of rows and lines for yz-plane
    Integer rowPtrXZ[] = {1,3,5,7,9};  // indices of rows and lines for xz-plane
    Integer * rowPtr;
    Integer i,j;

    switch(actOrientation_)
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
    dMat.Resize(nrElems2d);

    Double jacDetInv;
    jacDetInv = 1.0/ptelem->CalcJacobianDetAtIp(ip,ptCoord);

    for (i=0; i<nrElems2d; i++)
      for (j=0; j<nrElems2d; j++)
	dMat[i][j] = (*matMatrix)[rowPtr[i]-1][rowPtr[j]-1]*jacDetInv;	

}



  // calculated the D-matrix for the axisymmetric state
  void SmoothAxiInt::calcDMat(Matrix<Double> & dMat, Integer ip, Matrix<Double> & ptCoord)
  {
    ENTER_FCN( "SmoothAxiInt::calcDMat" );
    
    const Integer nrElemsAxi = 4;
    
    Integer rowPtrXY[] = {1,2,6,3};  // indices of rows and lines for xy-plane
    Integer rowPtrYZ[] = {2,3,4,1};  // indices of rows and lines for yz-plane
    Integer rowPtrXZ[] = {1,3,5,2};  // indices of rows and lines for xz-plane
    Integer * rowPtr;

    switch(actOrientation_)
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
    
    dMat.Resize(nrElemsAxi);

    Double jacDetInv;
    jacDetInv = 1.0/ptelem->CalcJacobianDetAtIp(ip,ptCoord);

    for (Integer i=0; i<nrElemsAxi; i++)
      for (Integer j=0; j<nrElemsAxi; j++)
	dMat[i][j] = (*matMatrix)[rowPtr[i]-1][rowPtr[j]-1] * jacDetInv;
  }



  // calculates the D-matrix of a 3d-problem 
  void smooth3DInt::calcDMat(Matrix<Double> & dMat, Integer ip, Matrix<Double> & ptCoord)
  {
    ENTER_FCN( "smooth3DInt::calcDMat" );

    const Integer nrElems3d = getDimD();
    
    Matrix<Double> * matMatrix =  ptMaterial->GetMatrix();
    dMat.Resize(nrElems3d);

    Double jacDetInv;
    jacDetInv = 1.0/ptelem->CalcJacobianDetAtIp(ip,ptCoord);

    for (Integer i=0; i<nrElems3d; i++)
      for (Integer j=0; j<nrElems3d; j++)
	dMat[i][j] = (*matMatrix)[i][j] * jacDetInv;	
  }




  // ===================================================================================
  // =================== standard con- and destructors (just for tracing) ==============
  // ===================================================================================

  // calculate (for 2D problems) by default in the xy-plane
  SmoothInt::SmoothInt(BaseFE * aptelem, MaterialData & matData) 
    : BDBInt(aptelem, matData), actOrientation_(xy)
  {
    ENTER_FCN( "SmoothInt::SmoothInt" );
    updateDMatInEveryIP_ = 1;
  }


  // calculate (for 2D problems) by default in the xy-plane
  SmoothInt::SmoothInt(MaterialData & matData) 
    : BDBInt(matData), actOrientation_(xy)
  {
    ENTER_FCN( "SmoothInt::SmoothInt" );
    updateDMatInEveryIP_ = 1;
  }
 

  SmoothInt::~SmoothInt()
  {
    ENTER_FCN( "SmoothInt::~SmoothInt" );
  }



  smoothPlainStrainInt::smoothPlainStrainInt(BaseFE * aptelem, MaterialData & matData) 
    : SmoothInt(aptelem, matData)
  {
    ENTER_FCN( "smoothPlainStrainInt::smoothPlainStrainInt" );

    ptelem=aptelem;
  }


  smoothPlainStrainInt::smoothPlainStrainInt(MaterialData & matData) 
    : SmoothInt(matData)
  {
    ENTER_FCN( "smoothPlainStrainInt::smoothPlainStrainInt" );
  }
 

  smoothPlainStrainInt::~smoothPlainStrainInt()
  {
    ENTER_FCN( "smoothPlainStrainInt::~smoothPlainStrainInt" );
  }


  SmoothAxiInt::SmoothAxiInt(BaseFE * aptelem, MaterialData & matData) 
    : SmoothInt(aptelem, matData)
  {
    ENTER_FCN( "SmoothAxiInt::SmoothAxiInt" );
    isaxi_ = TRUE;
    ptelem=aptelem;
  }


  SmoothAxiInt::SmoothAxiInt(MaterialData & matData) 
    : SmoothInt(matData)
  {
    ENTER_FCN( "SmoothAxiInt::SmoothAxiInt" );
    isaxi_ = TRUE;
  }
 

  SmoothAxiInt::~SmoothAxiInt()
  {
    ENTER_FCN( "SmoothAxiInt::~SmoothAxiInt" );
  }





  smooth3DInt::smooth3DInt(BaseFE * aptelem, MaterialData & matData) 
    : SmoothInt(aptelem, matData)
  {
    ENTER_FCN( "smooth3DInt::smooth3DInt" );
  }


  smooth3DInt::smooth3DInt(MaterialData & matData) 
    : SmoothInt(matData)
  {
    ENTER_FCN( "smooth3DInt::smooth3DInt" );
  }
 

  smooth3DInt::~smooth3DInt()
  {
    ENTER_FCN( "smooth3DInt::~smooth3DInt" );
  }

} // end namespace CoupledField
