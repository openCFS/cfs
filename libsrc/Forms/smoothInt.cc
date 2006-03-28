#include <iostream>
#include <fstream>

#include "smoothInt.hh"

namespace CoupledField
{



  // returns B - matrix for BDB
  void SmoothInt::calcBMat(Matrix<Double> & bMat, UInt ip, Matrix<Double> & ptCoord)
  {
    ENTER_FCN( "SmoothInt::calcBMat" );
    
    const UInt nrNodes  = ptelem->GetNumNodes();
    const UInt spaceDim = ptelem->GetDim();  
    const UInt nrDofs   = getNrDofs();  

    UInt actDim, actNode, j, k;
    
    
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
  }
  

  // calculates the D-matrix 
  void SmoothInt::calcDMat(Matrix<Double> & dMat, UInt ip, Matrix<Double> & ptCoord)
  {
    ENTER_FCN( "smooth3DInt::calcDMat" );

    const UInt nrElems3d = getDimD();
    
    Matrix<Double> matMatrix;
    ptMaterial->GetTensor(matMatrix,MECH_STIFFNESS_TENSOR,REAL,subTensorType_);
    UInt dimRow = matMatrix.GetSizeRow();
    UInt dimCol = matMatrix.GetSizeCol();

    dMat.Resize(dimRow, dimCol);

    Double jacDetInv;
    jacDetInv = 1.0/ptelem->CalcJacobianDetAtIp(ip,ptCoord);

    for (UInt i=0; i<dimRow; i++)
      for (UInt j=0; j<dimCol; j++)
        dMat[i][j] = matMatrix[i][j] * jacDetInv;    
  }


  // ===================================================================================
  // =================== standard con- and destructors (just for tracing) ==============
  // ===================================================================================

  SmoothInt::SmoothInt(BaseMaterial* matData, SubTensorType type) 
    : BDBInt(matData, type) 
  {
    ENTER_FCN( "SmoothInt::SmoothInt" );
    updateDMatInEveryIP_ = 1;

    if ( type == FULL ) {
      dimD_   = 6;
      nrDofs_ = 3;
    }
    else if ( type == PLANE_STRAIN ) {
      dimD_   = 3;
      nrDofs_ = 2;
    }
    else if ( type == AXI ) {
      dimD_   = 4;
      nrDofs_ = 2;
      isaxi_  = TRUE;
    }
  }
 

  SmoothInt::~SmoothInt()
  {
    ENTER_FCN( "SmoothInt::~SmoothInt" );
  }


} // end namespace CoupledField
