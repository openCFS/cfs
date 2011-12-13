// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <string>

#include "Domain/entityList.hh"
#include "Elements/basefe.hh"
#include "Forms/bdbInt.hh"
#include "MatVec/matrix.hh"
#include "Materials/baseMaterial.hh"
#include "smoothInt.hh"

namespace CoupledField
{



  // returns B - matrix for BDB
  void SmoothInt::CalcBMat(Matrix<Double> & bMat, UInt ip, const Matrix<Double> & ptCoord)
  {
    
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt spaceDim = ptelem->GetDim();  
    const UInt nrDofs   = getNrDofs();  

    UInt actDim, actNode;
    // UInt j, k; // TODO: Unused variables j,k
    
    
    bMat.Resize(getDimD(), numFncs * nrDofs);
    bMat.Init();
    
    // local shape functions derived after global coords (format: numFncs x spaceDim)
    Matrix<Double> xiDx;

    ptelem->GetGlobDerivShFncAtIp(xiDx, ip, ptCoord, it1_.GetElem());


    for(actDim=0; actDim < spaceDim; actDim++)
      for(actNode=0; actNode < numFncs; actNode++)
        bMat[actDim][actNode * spaceDim + actDim] = xiDx[actNode][actDim];

    switch(spaceDim)
      {
      case 2:
        // j = 1;
        // k = 0;
        
        for (actNode = 0; actNode < numFncs; actNode++)
          {
            bMat[spaceDim][actNode * spaceDim + 1] = xiDx[actNode][0];
            bMat[spaceDim][actNode * spaceDim]     = xiDx[actNode][1];
          }
        break;


      case 3:
        UInt actDim=spaceDim;
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
  }
  

  // calculates the D-matrix 
  void SmoothInt::calcDMat(Matrix<Double> & dMat, UInt ip, Matrix<Double> & ptCoord)
  {

    Matrix<Double> matMatrix;
    ptMaterial->GetTensor(matMatrix,MECH_STIFFNESS_TENSOR,Global::REAL,subTensorType_);
    UInt dimRow = matMatrix.GetNumRows();
    UInt dimCol = matMatrix.GetNumCols();

    dMat.Resize(dimRow, dimCol);

    Double jacDetInv;
    jacDetInv = 1.0/ptelem->CalcJacobianDetAtIp(ip,ptCoord, it1_.GetElem() );

    for (UInt i=0; i<dimRow; i++)
      for (UInt j=0; j<dimCol; j++)
        dMat[i][j] = matMatrix[i][j] * jacDetInv;    
  }


  // ===================================================================================
  // =================== standard con- and destructors (just for tracing) ==============
  // ===================================================================================

  SmoothInt::SmoothInt(BaseMaterial* matData, SubTensorType type,
                       bool coordUpdate ) 
    : BDBInt(matData, type) 
  {

    name_ = "SmoothInt";
    updateDMatInEveryIP_ = 1;
    coordUpdate_ = coordUpdate;

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
      isaxi_  = true;
    }
  }
 

  SmoothInt::~SmoothInt()
  {
  }


} // end namespace CoupledField
