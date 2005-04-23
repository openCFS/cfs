#include <iostream>
#include <fstream>

#include "linElastInt.hh"

namespace CoupledField
{



  // returns B - matrix for BDB
  void linElastInt::calcBMat( Matrix<Double> &bMat, Integer ip,
                              Matrix<Double> &ptCoord ) {

    ENTER_FCN( "linElastInt::calcBMat" );

    const Integer nrNodes  = ptelem->GetNumNodes();
    const Integer spaceDim = ptelem->GetDim();  
    const Integer nrDofs   = getNrDofs();  

    Integer actDim, actNode, help, j, k;
    
    
    bMat.Resize(getDimD(), nrNodes * nrDofs);
    
    // local shape functions derived after global coords
    // (format: nrNodes x spaceDim)
    Matrix<Double> xiDx;

    if (isSetIntPoint_) 
      ptelem->GetGlobDerivShFnc(xiDx, intPoint_, ptCoord);
    else
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
	    Vector<Double> ShpFncAtIp;
	    Vector<Double> CoordAtIP;

	    if (isSetIntPoint_) 
	      ptelem->GetShFnc(ShpFncAtIp,intPoint_);
	    else
	      ptelem->GetShFncAtIp(ShpFncAtIp,ip);

	    CoordAtIP = ptCoord * ShpFncAtIp;

	    for (actNode = 0; actNode < nrNodes; actNode++)	     
	      bMat[idxtheta-1][actNode * spaceDim] =
                ShpFncAtIp[actNode] / CoordAtIP[0];
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

    isSetIntPoint_ = FALSE;
  }

  void linElastInt::
  CalcAxiMaterialMat(Matrix<Double> & dMat, enum orientation2D actOrientation)
  {
    ENTER_FCN( "linElastInt::CalcAxiMaterialMat" );
    
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
	
    Matrix<Double> const & matMatrix =  *(ptMaterial->GetMatrix());
    dMat.Resize(nrElemsAxi);

    for (i=0; i<nrElemsAxi; i++)
      for (j=0; j<nrElemsAxi; j++)
	dMat[i][j] = matMatrix[rowPtr[i]-1][rowPtr[j]-1];

  }

  // calculated the D-matrix for the plain strain state
  void linElastInt::
  CalcPlaneStrainMaterialMat( Matrix<Double> &dMat,
                              enum orientation2D actOrientation ){

    ENTER_FCN( "linElastInt::CalcPlaneStrainMaterialMat" );

    const Integer nrElems2d = getDimD();
    
    Integer rowPtrXY[] = {1,2,6};  // indices of rows and lines for xy-plane
    Integer rowPtrYZ[] = {2,3,4};  // indices of rows and lines for yz-plane
    Integer rowPtrXZ[] = {1,3,5};  // indices of rows and lines for xz-plane
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
	
    Matrix<Double> const  & matMatrix = *( ptMaterial->GetMatrix());
    
    dMat.Resize(nrElems2d);

    for (i=0; i<nrElems2d; i++)
      for (j=0; j<nrElems2d; j++)
	dMat[i][j] = matMatrix[rowPtr[i]-1][rowPtr[j]-1];	
  }


  // calculates the D-matrix of a 3d-problem 
  void linElastInt::Calc3DMaterialMat(Matrix<Double> & dMat)
  {
    ENTER_FCN( "linElastInt::Calc3DMaterialMat" );

    const Integer nrElems3d = getDimD();
    
    Matrix<Double> const & matMatrix =  *(ptMaterial->GetMatrix());
    
    dMat.Resize(nrElems3d);

    for (Integer i=0; i<nrElems3d; i++)
      for (Integer j=0; j<nrElems3d; j++)
	dMat[i][j] = matMatrix[i][j];	
  }

  
  // calculated the D-matrix for the plain strain state
  void mechPlainStrainInt::calcDMat(Matrix<Double> & dMat)
  {
    ENTER_FCN( "mechPlainStrainInt::calcDMat" );
    CalcPlaneStrainMaterialMat(dMat,actOrientation);
  }


  // calculated the D-matrix for the axisymmetric state
  void mechAxiInt::calcDMat(Matrix<Double> & dMat)
  {
    ENTER_FCN( "mechAxiInt::calcDMat" );
  
    CalcAxiMaterialMat(dMat, actOrientation);
  
  }
  
  // calculates the D-matrix of a 3d-problem 
  void mech3DInt::calcDMat(Matrix<Double> & dMat)
  {
    ENTER_FCN( "mech3DInt::calcDMat" );
    Calc3DMaterialMat(dMat);
  }




  // =========================================================================
  // =============== standard con- and destructors (just for tracing) ========
  // =========================================================================

  // NOTE: We hardcode the orientation for 2D simulations to use the yz plane!

  linElastInt::linElastInt( BaseFE *aptelem, MaterialData &matData ) :
    BDBInt(aptelem, matData), actOrientation(yz) {
    ENTER_FCN( "linElastInt::linElastInt" );
  }

  linElastInt::linElastInt( MaterialData & matData) :
    BDBInt(matData), actOrientation(yz) {
    ENTER_FCN( "linElastInt::linElastInt" );
  }

  linElastInt::~linElastInt() {
    ENTER_FCN( "linElastInt::~linElastInt" );
  }

  mechPlainStrainInt::mechPlainStrainInt(BaseFE * aptelem,
                                         MaterialData & matData) :
    linElastInt(aptelem, matData) {
    ENTER_FCN( "mechPlainStrainInt::mechPlainStrainInt" );

    ptelem=aptelem;
  }


  mechPlainStrainInt::mechPlainStrainInt(MaterialData & matData) 
    : linElastInt(matData)
  {
    ENTER_FCN(" mechPlainStrainInt::mechPlainStrainInt" );

  }
 

  mechPlainStrainInt::~mechPlainStrainInt()
  {
    ENTER_FCN( "mechPlainStrainInt::~mechPlainStrainInt" );
  }

  mechAxiInt::mechAxiInt(BaseFE * aptelem, MaterialData & matData) 
    : linElastInt(aptelem, matData)
  {
    ENTER_FCN( "mechAxiInt::mechAxiInt" );

    isaxi_ = TRUE;
    ptelem=aptelem;
  }


  mechAxiInt::mechAxiInt(MaterialData & matData) 
    : linElastInt(matData)
  {
    ENTER_FCN( "mechAxiInt::mechAxiInt" );

    isaxi_ = TRUE;
  }
 

  mechAxiInt::~mechAxiInt()
  {
    ENTER_FCN( "mechAxiInt::~mechAxiInt" );
  }


  mech3DInt::mech3DInt(BaseFE * aptelem, MaterialData & matData) 
    : linElastInt(aptelem, matData)
  {
    ENTER_FCN( "mech3DInt::mech3DInt" );

  }


  mech3DInt::mech3DInt(MaterialData & matData) 
    : linElastInt(matData)
  {
    ENTER_FCN( "mech3DInt::mech3DInt" );

  }
 

  mech3DInt::~mech3DInt()
  {
    ENTER_FCN( "mech3DInt::~mech3DInt" );

  }

} // end namespace CoupledField
