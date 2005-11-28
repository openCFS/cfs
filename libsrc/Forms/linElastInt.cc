#include <iostream>
#include <fstream>

#include "linElastInt.hh"

namespace CoupledField
{

  // =====================
  //   CalcElementMatrix
  // =====================
  void linElastInt::CalcElementMatrix( Matrix<Double> &ptCoord,
                                  Matrix<Double> &elemMat ) {

    ENTER_FCN( "linElastInt::CalcElementMatrix" );

    Boolean Softening;

    if (softeningModel_=="BK1") {
      softeningPart_ = "bendingBK1";
      BDBInt::CalcElementMatrix(ptCoord, elemMat);
      Matrix<Double> helpMat = elemMat;
      // std::cout << "Bending Mat:\n" << helpMat << std::endl;

      softeningPart_ = "shearBK1";
      BDBInt::CalcElementMatrix(ptCoord, elemMat);
      elemMat += helpMat;
      //     std::cout << "Total Mat:\n" << elemMat << std::endl;
    }
    else if (softeningModel_=="SRI") {
      softeningPart_ = "bendingSRI";
      BDBInt::CalcElementMatrix(ptCoord, elemMat);
      Matrix<Double> helpMat = elemMat;
      // std::cout << "Bending Mat:\n" << helpMat << std::endl;

      softeningPart_ = "shearSRI";
      BDBInt::CalcElementMatrix(ptCoord, elemMat);
      elemMat += helpMat;
      //     std::cout << "Total Mat:\n" << elemMat << std::endl;
    }
    else {
      BDBInt::CalcElementMatrix(ptCoord, elemMat);
    }
  }

  // returns B - matrix for BDB
  void linElastInt::calcBMat( Matrix<Double> &bMat, UInt ip,
                              Matrix<Double> &ptCoord ) {

    ENTER_FCN( "linElastInt::calcBMat" );

    const UInt nrNodes  = ptelem->GetNumNodes();
    const UInt spaceDim = ptelem->GetDim();  
    const UInt nrDofs   = getNrDofs();  

    UInt actDim, actNode, j, k;
    
    
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
            UInt idxtheta = getDimD();
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
        UInt actDim=spaceDim;
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
    
    const UInt nrElemsAxi = 4;
    
    UInt rowPtrXY[] = {1,2,6,3};  // indices of rows and lines for xy-plane
    UInt rowPtrYZ[] = {2,3,4,1};  // indices of rows and lines for yz-plane
    UInt rowPtrXZ[] = {1,3,5,2};  // indices of rows and lines for xz-plane
    UInt * rowPtr;
    UInt i,j;

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

    //check for softening model
    if (softeningPart_ == "bendingBK1" ) {
      Error("BK1-Softening-Model for 2D not implemented",__FILE__,__LINE__);
    }

    else if (softeningPart_ == "bendingSRI" ) {
      UInt idx = nrElemsAxi-2;
      dMat[idx][idx] = 0.0;
    }

    else if (softeningPart_ == "shearBK1" ) {
      Error("BK1-Softening-Model for 2D not implemented",__FILE__,__LINE__);
    }

    else if (softeningPart_ == "shearSRI" ) {
      for (UInt i=0; i<nrElemsAxi-2; i++) {
	for (UInt j=0; j<nrElemsAxi-2; j++) {
	  dMat[i][j] = 0.0;
	}
      }

      UInt rowidx=nrElemsAxi-1;
      for (UInt i=0; i<nrElemsAxi; i++) {
	dMat[rowidx][i] = 0.0;
      }

      UInt colidx=nrElemsAxi-1;
      for (UInt i=0; i<nrElemsAxi; i++) {
	dMat[i][colidx] = 0.0;
      }

    }

  }

  // calculated the D-matrix for the plain strain state
  void linElastInt::
  CalcPlaneStrainMaterialMat( Matrix<Double> &dMat,
                              enum orientation2D actOrientation ){

    ENTER_FCN( "linElastInt::CalcPlaneStrainMaterialMat" );

    const UInt nrElems2d = getDimD();
    
    UInt rowPtrXY[] = {1,2,6};  // indices of rows and lines for xy-plane
    UInt rowPtrYZ[] = {2,3,4};  // indices of rows and lines for yz-plane
    UInt rowPtrXZ[] = {1,3,5};  // indices of rows and lines for xz-plane
    UInt * rowPtr;
    UInt i,j;

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

    //check for softening model
    if (softeningPart_ == "bendingBK1" ) {
      Error("BK1-Softening-Model for 2D not implemented",__FILE__,__LINE__);
    }

    else if (softeningPart_ == "bendingSRI" ) {
      UInt idx = nrElems2d-1;
      dMat[idx][idx] = 0.0;
    }

    else if (softeningPart_ == "shearBK1" ) {
      Error("BK1-Softening-Model for 2D not implemented",__FILE__,__LINE__);
    }

    else if (softeningPart_ == "shearSRI" ) {
      for (UInt i=0; i<nrElems2d; i++) {
	for (UInt j=0; j<nrElems2d; j++) {
	  dMat[i][j] = 0.0;
	}
      }
    }
   
  }


  // calculates the D-matrix of a 3d-problem 
  void linElastInt::Calc3DMaterialMat(Matrix<Double> & dMat)
  {
    ENTER_FCN( "linElastInt::Calc3DMaterialMat" );

    const UInt nrElems3d = getDimD();
    
    Matrix<Double> const & matMatrix =  *(ptMaterial->GetMatrix());
    
    dMat.Resize(nrElems3d);

    for (UInt i=0; i<nrElems3d; i++) {
      for (UInt j=0; j<nrElems3d; j++) {
        dMat[i][j] = matMatrix[i][j];
      }
    }

    if (softeningPart_ == "bendingBK1" ) {
      
      if (  maxEdgeLength_ < 1e-15 ) {
	Error("linElastInt::Calc3DMaterialMat: maxEdgeLength_=0",__FILE__,__LINE__);
      }

      Double f1, factor;
      f1 = (minEdgeLength_* minEdgeLength_);
      factor =  2*f1 / ( 2*f1 + maxEdgeLength_ * maxEdgeLength_);

      dMat[3][3] *= factor;   
      dMat[4][4] *= factor;   
      dMat[5][5] *= factor;   
    }

    else if (softeningPart_ == "bendingSRI" ) {
      for (UInt i=3; i<nrElems3d; i++) {
	for (UInt j=3; j<nrElems3d; j++) {
	  dMat[i][j] = 0.0;
	}
      }
    }

    else if (softeningPart_ == "shearBK1" ) {
      for (UInt i=0; i<3; i++) {
	for (UInt j=0; j<3; j++) {
	  dMat[i][j] = 0.0;
	}
      }

      Double f1, factor;
      f1 = maxEdgeLength_ * maxEdgeLength_;
      factor = f1 / ( f1 + 2* minEdgeLength_* minEdgeLength_ ); 
      dMat[3][3] *= factor;   
      dMat[4][4] *= factor;   
      dMat[5][5] *= factor; 
    }

    else if (softeningPart_ == "shearSRI" ) {
      for (UInt i=0; i<3; i++) {
	for (UInt j=0; j<3; j++) {
	  dMat[i][j] = 0.0;
	}
      }
    }

    //   std::cout << "Dmat:\n" << dMat << std::endl;
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
