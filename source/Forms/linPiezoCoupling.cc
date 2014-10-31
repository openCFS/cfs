// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <stddef.h>
#include <fstream>
#include <string>

#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"
#include "Domain/domain.hh"
#include "Domain/elem.hh"
#include "Domain/entityList.hh"
#include "Elements/basefe.hh"
#include "General/environment.hh"
#include "MatVec/matrix.hh"
#include "MatVec/vector.hh"
#include "Materials/baseMaterial.hh"
#include "linPiezoCoupling.hh"
#include "Optimization/Design/DesignSpace.hh"

DECLARE_LOG(forms)

namespace CoupledField {


  // ============
  //   calcAMat
  // ============
  void linPiezoCoupling::calcAMat( Matrix<Double> &aMat, UInt ip,
				   const Matrix<Double> &ptCoord ) {


    // Obtain info on problem sizes
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );

    // Set type of ansatz function , but do not recalculate
    // integration points
    ptelem->SetAnsatzFct( ansatzFct1_, false );
    //    const UInt nDofMech = 3;
    //    const UInt nRowsD   = 6;

    // Set correct size of matrix A and initialise with zeros
    aMat.Resize( numFncs * numDofsA_, matDimRow_ );
    aMat.Init();

    // Get derivatives of local shape functions with respect to global coords
    // (format: numFncs x spaceDim)
    Matrix<Double> xiDx;
    ptelem->GetGlobDerivShFncAtIp( xiDx, ip, ptCoord, it1_.GetElem() );

    //std::cerr << "xiDx = \n" << xiDx << std::endl;

    // The matrix aMat can be seen as a numFncs x 1 block-vector.
    // The k-th entry of this block vector corresponds to the matrix
    // A of the ADB product evaluated at the k-th node of the finite
    // element. We assemble aMat in a top-down fashion.
    UInt actInd = 0;
    UInt actNode;

    if ( subTensorType_ == FULL ) {
      for( actNode = 0; actNode < numFncs; actNode++ ) {

        // 1st row of sub-matrix A(actNode)
        aMat[actInd][0] = xiDx[actNode][0];
        aMat[actInd][4] = xiDx[actNode][2];
        aMat[actInd][5] = xiDx[actNode][1];
        actInd++;

        // 2nd row of sub-matrix A(actNode)
        aMat[actInd][1] = xiDx[actNode][1];
        aMat[actInd][3] = xiDx[actNode][2];
        aMat[actInd][5] = xiDx[actNode][0];
        actInd++;

        // 3rd row of sub-matrix A(actNode)
        aMat[actInd][2] = xiDx[actNode][2];
        aMat[actInd][3] = xiDx[actNode][1];
        aMat[actInd][4] = xiDx[actNode][0];
        actInd++;
      }
    }
    else if ( subTensorType_ == AXI ) {
      // we assume that the first coordinate equals r and 
      // the second z.

      UInt j;
      Double coordAtIp;
      Vector<Double> ShpFncAtIp;
      ptelem->GetShFncAtIp( ShpFncAtIp, ip, it1_.GetElem() );

      for( actNode = 0; actNode < numFncs; actNode++ ) {

        // 1st row of sub-matrix A(actNode)
        aMat[actInd][0] = xiDx[actNode][0];    // dN/dr
        aMat[actInd][2] = xiDx[actNode][1];    // dN/dz

        // For the entry 1/r things are more complicated
        coordAtIp = 0.0;
        for( j = 0; j < numFncs; j++ ) {
          coordAtIp += ptCoord[0][j] * ShpFncAtIp[j];
        }
        aMat[actInd][3] = ShpFncAtIp[actNode] / coordAtIp;    // 1/r
        actInd++;

        // 2nd row of sub-matrix A(actNode)
        aMat[actInd][1] = xiDx[actNode][1];   // dN/dz
        aMat[actInd][2] = xiDx[actNode][0];   // dN/dr

        actInd++;

      }
    }

    else if ( subTensorType_ == PLANE_STRAIN || subTensorType_ == PLANE_STRESS) {
      // we assume that the first coordinate equals y and the second z.

      for( actNode = 0; actNode < numFncs; actNode++ ) {

        // 1st row of sub-matrix A(actNode)
        aMat[actInd][0] = xiDx[actNode][0];   // dN/dy
        aMat[actInd][2] = xiDx[actNode][1];   // dN/dz
        actInd++;

        // 2nd row of sub-matrix A(actNode)
        aMat[actInd][1] = xiDx[actNode][1];    // dN/dz
        aMat[actInd][2] = xiDx[actNode][0];    // dN/dy
        actInd++;
      }
    }
   // std::cerr << "linPiezoCoupling::aMat = \n" << aMat << std::endl;
  }


  // ============
  //   calcBMat
  // ============
  void linPiezoCoupling::CalcBMat( Matrix<Double> &bMat, UInt ip,
				   const Matrix<Double> &ptCoord ) {


    // Set type of ansatz function , but do not recalculate
    // integration points
    ptelem->SetAnsatzFct( ansatzFct2_, false );

    // Obtain info on number of element's nodes
    UInt numFncs = ptelem->GetNumFncs( ansatzFct2_ );

    // Set correct size of matrix B and initialise with zeros
    bMat.Resize( numDofsA_, numFncs );
    bMat.Init();

    // Get derivatives of local shape functions with respect to global coords
    // (format: numFncs x spaceDim)
    Matrix<Double> xiDx;
    ptelem->GetGlobDerivShFncAtIp( xiDx, ip, ptCoord, it1_.GetElem() );


    // The matrix bMat can be seen as a 1 x numFncs block-vector.
    // The k-th entry of this block vector corresponds to the matrix
    // B of the ADB product evaluated at the k-th node of the finite
    // element. We simply must transpose xiDx.
    if ( subTensorType_ == FULL ) {
      for( UInt actNode = 0; actNode < numFncs; actNode++ ) {
        bMat[0][actNode] = xiDx[actNode][0];
        bMat[1][actNode] = xiDx[actNode][1];
        bMat[2][actNode] = xiDx[actNode][2];
      }
    }
    else if ( subTensorType_ == AXI ) {
      // element. We assume that the first coordinate equals r and the
      // second z.
      for( UInt actNode = 0; actNode < numFncs; actNode++ ) {
        bMat[0][actNode] = xiDx[actNode][0];
        bMat[1][actNode] = xiDx[actNode][1];
      }
    }
   else if ( subTensorType_ == PLANE_STRAIN || subTensorType_ == PLANE_STRESS) {
     // We assume that the first coordinate equals y and the
     // second z.
     for( UInt actNode = 0; actNode < numFncs; actNode++ ) {
       bMat[0][actNode] = xiDx[actNode][0];
       bMat[1][actNode] = xiDx[actNode][1];
     }
   }

    
    //std::cerr << "linPiezoCoupling::bMat = \n" << bMat << std::endl;
  }


  // ============
  //   calcDMat
  // ============
  void linPiezoCoupling::calcDMat(Matrix<Double> &dMat, const Elem* elem, DesignElement::Type direction)
  {
    Matrix<Double> matMatrix;

    // consider standard FEM, SIMP and FMO
    if(elem != NULL &&
       domain->HasNonDensityDesignMaterial() &&
       domain->GetErsatzMaterial()->GetPiezoCouplingTensor(matMatrix, elem, direction))
    {
      matMatrix.Transpose(dMat); // this is the FMO tensor
      LOG_DBG3(forms) << GetName() << "::calcDMat("  << (elem != NULL ? Integer(elem->elemNum) : -1)  << ") -> FMO=" << dMat.ToString();
    }
    else
    {
      ptMaterial->GetTensor(matMatrix,PIEZO_TENSOR,matDataType_,subTensorType_);
      matMatrix.Transpose(dMat);

      // possibly SIMP?
      Double density = elem != NULL ? GetErsatzMaterialFactor(elem) : 1.0;
      if(density != 1.0) dMat *= density;
      LOG_DBG3(forms) << GetName() << "::calcDMat("  << (elem != NULL ? Integer(elem->elemNum) : -1)  << ") density=" << density << " -> " << dMat.ToString();
    }
  }


  // ============
  //   Constructor
  // ============
  linPiezoCoupling::linPiezoCoupling(BaseMaterial* matData,
                                     SubTensorType type) 
    : ADBInt( matData, type ) {

    name_ = "linPiezoCoupling";

    if ( type == FULL ) {
      numDofsA_  = 3;
      numDofsB_  = 1;
      matDimRow_ = 6;
      matDimCol_ = 3;
    }
    else if ( type == PLANE_STRAIN || type == PLANE_STRESS ) {
      numDofsA_  = 2;
      numDofsB_  = 1;
      matDimRow_ = 3;
      matDimCol_ = 2;
    }
    else if ( type == AXI ) {
      numDofsA_  = 2;
      numDofsB_  = 1;
      matDimRow_ = 4;
      matDimCol_ = 2;
      isaxi_     = true;
    }
  }
}
