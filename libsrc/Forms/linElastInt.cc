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

	if (isaxi_)
	  {
	    Integer idxtheta = getDimD();
	    std::vector<Double> ShpFncAtIp;
	    std::vector<Double> CoordAtIP;

	    ptelem->GetShFncAtIp(ShpFncAtIp,ip);
	    CoordAtIP = ptCoord * ShpFncAtIp;

	    for (actNode = 0; actNode < nrNodes; actNode++)	     
	      bMat[idxtheta-1][actNode * spaceDim] = ShpFncAtIp[actNode] / CoordAtIP[0];
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
  


  void linElastInt::
  CalcAxiMaterialMat(Matrix<Double> & dMat, enum orientation2D actOrientation)
  {
#ifdef TRACE
    (*trace) << "entering linElastInt::CalcAxiMaterialMat " << std::endl;
#endif
    
    const Integer nrElemsAxi = 4;
    
    Integer rowPtrXY[] = {1,2,6,3};  // indices of rows and lines for xy-plane
    Integer rowPtrYZ[] = {2,3,4,1};  // indices of rows and lines for yz-plane
    Integer rowPtrXZ[] = {1,3,5,2};  // indices of rows and lines for xz-plane
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
    
    dMat.Resize(nrElemsAxi);

    for (i=0; i<nrElemsAxi; i++)
      for (j=0; j<nrElemsAxi; j++)
	dMat[i][j] = (*matMatrix)[rowPtr[i]-1][rowPtr[j]-1];

}


  
  // calculated the D-matrix for the plain strain state
  void mechPlainStrainInt::calcDMat(Matrix<Double> & dMat)
  {
#ifdef TRACE
  (*trace) << "entering mechPlainStrainInt::calcDMat " << std::endl;
#endif
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
    
    dMat.Resize(nrElems2d);

    for (i=0; i<nrElems2d; i++)
      for (j=0; j<nrElems2d; j++)
	dMat[i][j] = (*matMatrix)[rowPtr[i]-1][rowPtr[j]-1];	
}




  // calculated the D-matrix for the axisymmetric state
  void mechAxiInt::calcDMat(Matrix<Double> & dMat)
  {
#ifdef TRACE
  (*trace) << "entering mechAxiInt::calcDMat " << std::endl;
#endif
  
  CalcAxiMaterialMat(dMat, actOrientation);
  
  }
  



  // calculates the D-matrix of a 3d-problem 
  void mech3DInt::calcDMat(Matrix<Double> & dMat)
  {
#ifdef TRACE
  (*trace) << "entering mech3DInt::calcDMat " << std::endl;
#endif
    const Integer nrElems3d = getDimD();
    
    Matrix<Double> * matMatrix =  ptMaterial->GetMatrix();
    
    dMat.Resize(nrElems3d);

    for (Integer i=0; i<nrElems3d; i++)
      for (Integer j=0; j<nrElems3d; j++)
	dMat[i][j] = (*matMatrix)[i][j];	
  }




  // ===================================================================================
  // =================== standard con- and destructors (just for tracing) ==============
  // ===================================================================================

  // calculate (for 2D problems) by default in the xy-plane
  linElastInt::linElastInt(BaseFE * aptelem, MaterialData & matData) 
    : BDBInt(aptelem, matData), actOrientation(xy)
  {
#ifdef TRACE
    (*trace) << "entering linElastInt::linElastInt" << std::endl;
#endif
  }


  // calculate (for 2D problems) by default in the xy-plane
  linElastInt::linElastInt(MaterialData & matData) 
    : BDBInt(matData), actOrientation(xy)
  {
#ifdef TRACE
    (*trace) << "entering linElastInt::linElastInt" << std::endl;
#endif
  }
 

  linElastInt::~linElastInt()
  {
#ifdef TRACE
    (*trace) << "entering linElastInt::~linElastInt" << std::endl;
#endif
  }



  mechPlainStrainInt::mechPlainStrainInt(BaseFE * aptelem, MaterialData & matData) 
    : linElastInt(aptelem, matData)
  {
#ifdef TRACE
    (*trace) << "entering mechPlainStrainInt::mechPlainStrainInt" << std::endl;
#endif

    ptelem=aptelem;
  }


  mechPlainStrainInt::mechPlainStrainInt(MaterialData & matData) 
    : linElastInt(matData)
  {
#ifdef TRACE
    (*trace) << "entering mechPlainStrainInt::mechPlainStrainInt" << std::endl;
#endif
  }
 

  mechPlainStrainInt::~mechPlainStrainInt()
  {
#ifdef TRACE
    (*trace) << "entering mechPlainStrainInt::~mechPlainStrainInt" << std::endl;
#endif
  }

  mechAxiInt::mechAxiInt(BaseFE * aptelem, MaterialData & matData) 
    : linElastInt(aptelem, matData)
  {
#ifdef TRACE
    (*trace) << "entering mechAxiInt::mechAxiInt" << std::endl;
#endif
    isaxi_ = TRUE;
    ptelem=aptelem;
  }


  mechAxiInt::mechAxiInt(MaterialData & matData) 
    : linElastInt(matData)
  {
#ifdef TRACE
    (*trace) << "entering mechAxiInt::mechAxiInt" << std::endl;
#endif
    isaxi_ = TRUE;
  }
 

  mechAxiInt::~mechAxiInt()
  {
#ifdef TRACE
    (*trace) << "entering mechAxiInt::~mechAxiInt" << std::endl;
#endif
  }


  mech3DInt::mech3DInt(BaseFE * aptelem, MaterialData & matData) 
    : linElastInt(aptelem, matData)
  {
#ifdef TRACE
    (*trace) << "entering mech3DInt::mech3DInt" << std::endl;
#endif
  }


  mech3DInt::mech3DInt(MaterialData & matData) 
    : linElastInt(matData)
  {
#ifdef TRACE
    (*trace) << "entering mech3DInt::mech3DInt" << std::endl;
#endif
  }
 

  mech3DInt::~mech3DInt()
  {
#ifdef TRACE
    (*trace) << "entering mech3DInt::~mech3DInt" << std::endl;
#endif
  }

} // end namespace CoupledField
