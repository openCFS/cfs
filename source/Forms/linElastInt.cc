// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <fstream>
#include <string>

#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"
#include "Domain/domain.hh"
#include "Domain/elem.hh"
#include "Domain/entityList.hh"
#include "Elements/basefe.hh"
#include "Forms/bdbInt.hh"
#include "General/environment.hh"
#include "General/exception.hh"
#include "MatVec/exprt/xpr2.hh"
#include "MatVec/matrix.hh"
#include "MatVec/vector.hh"
#include "Materials/baseMaterial.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "linElastInt.hh"

DECLARE_LOG(forms)

namespace CoupledField {

// =====================
//   CalcElementMatrix
// =====================
void linElastInt::CalcElementMatrix( Matrix<Double>& elemMat,
    EntityIterator& ent1,
    EntityIterator& ent2,
    const DesignElement::Type direction) {

  //  softeningModel_ = "ICM_Taylor";

  if (softeningModel_=="BK1") {
    softeningPart_ = "bendingBK1";
    //std::cout << "Bending:" <<  std::endl;
    BDBInt::CalcElementMatrix( elemMat, ent1, ent2 );
    Matrix<Double> helpMat = elemMat;
    // std::cout << "Bending Mat:\n" << helpMat << std::endl;

    softeningPart_ = "shearBK1";
    //std::cout << "Shear:" <<  std::endl;
    CalcElementMatrixShearBK1( elemMat, ent1, ent2 );
    //std::cout << "Shear Mat:\n" << elemMat << std::endl;
    elemMat += helpMat;
    //     std::cout << "Total Mat:\n" << elemMat << std::endl;
  }
  else if (softeningModel_=="SRI") {
    softeningPart_ = "bendingSRI";
    BDBInt::CalcElementMatrix( elemMat, ent1, ent2 );
    Matrix<Double> helpMat = elemMat;
    // std::cout << "Bending Mat:\n" << helpMat << std::endl;

    softeningPart_ = "shearSRI";
    BDBInt::CalcElementMatrix( elemMat, ent1, ent2);
    elemMat += helpMat;
    //     std::cout << "Total Mat:\n" << elemMat << std::endl;
  }
  else if (softeningModel_=="icModesTW") {
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );
    if ( ( ptelem->GetNumFncs( ansatzFct1_) == 4 &&
        ptelem->feType() == Elem::QUAD4 ) ||
        ( ptelem->GetNumFncs( ansatzFct1_) == 8 &&
            ptelem->feType() == Elem::HEXA8 )  ) {
      //just for quadrilaterals and hexahedrals we can do the softening with
      //incompatible modes
      CalcElementMatrixICM( elemMat, ent1, ent2);
    }
    else {
      BDBInt::CalcElementMatrix( elemMat, ent1, ent2, direction);
    }
  }
  else {
    // No softening
    BDBInt::CalcElementMatrix( elemMat, ent1, ent2, direction);
  }
}

void linElastInt::ReorderBLikeMatrix(Matrix<Double>& in, Matrix<Double>& out, UInt ip, BaseFE* elem, const Matrix<Double>& ptCoord){

  const UInt numFncs  = elem->GetNumFncs( ansatzFct1_ );
  const UInt spaceDim = elem->GetDim();  
  const UInt nrDofs   = getNrDofs();  
  const Vector<Double> * intPoints = elem->GetIntPoints();

  UInt actDim, actNode;
  // UInt j, k; // TODO: Unused variables j,k

  out.Resize(getDimD(), numFncs * nrDofs);
  out.Init();

  for(actDim=0; actDim < spaceDim; actDim++)
    for(actNode=0; actNode < numFncs; actNode++)
      out[actDim][actNode * spaceDim + actDim] = in[actNode][actDim];

  switch(spaceDim)
  {
  case 2:
    // j = 1;
    // k = 0;

    for (actNode = 0; actNode < numFncs; actNode++)
    {
      out[spaceDim][actNode * spaceDim + 1] = in[actNode][0];
      out[spaceDim][actNode * spaceDim]     = in[actNode][1];
    }

    if (isaxi_)
    {
      UInt idxtheta = getDimD();
      Vector<Double> ShpFncAtIp;
      Vector<Double> CoordAtIP;

      if (isSetIntPoint_) {
        elem->GetShFnc(ShpFncAtIp,intPoint_, it1_.GetElem());
        elem->Local2GlobalCoord( CoordAtIP, intPoint_,
                                 ptCoord, it1_.GetElem() );
      } else {
        elem->GetShFncAtIp(ShpFncAtIp,ip, it1_.GetElem() );
        elem->Local2GlobalCoord( CoordAtIP, intPoints[ip-1],
                                 ptCoord, it1_.GetElem() );
      }


      for (actNode = 0; actNode < numFncs; actNode++)          
        out[idxtheta-1][actNode * spaceDim] =
          ShpFncAtIp[actNode] / CoordAtIP[0];
    }

    break;


  case 3:
    UInt actDim=spaceDim;
    for (actNode = 0; actNode < numFncs; actNode++)
    {
      out[actDim][actNode * spaceDim + 1] = in[actNode][2];
      out[actDim][actNode * spaceDim + 2] = in[actNode][1];
    }

    actDim++;
    for (actNode = 0; actNode < numFncs; actNode++)
    {
      out[actDim][actNode * spaceDim]     = in[actNode][2];
      out[actDim][actNode * spaceDim + 2] = in[actNode][0];
    }

    actDim++;
    for (actNode = 0; actNode < numFncs; actNode++)
    {
      out[actDim][actNode * spaceDim]     = in[actNode][1];
      out[actDim][actNode * spaceDim + 1] = in[actNode][0];
    }
    break;
  }
}

// returns B - matrix for BDB
void linElastInt::CalcBMat( Matrix<Double> &bMat, UInt ip,
    const Matrix<Double> &ptCoord ) {

  // local shape functions derived after global coords
  // (format: numFncs x spaceDim)
  Matrix<Double> xiDx;
  ptelem->SetAnsatzFct( ansatzFct1_ );

  if (isSetIntPoint_) 
    ptelem->GetGlobDerivShFnc(xiDx, intPoint_, ptCoord, it1_.GetElem() );
  else
    ptelem->GetGlobDerivShFncAtIp(xiDx, ip, ptCoord, it1_.GetElem() );
  
  // LOG_DBG3(forms) << "calcBMat: xiDx: " << xiDx.ToString() << std::endl;
  
  ReorderBLikeMatrix(xiDx, bMat, ip, ptelem, ptCoord);

  // LOG_DBG2(forms) << "calcBMat: bMat: " << bMat.ToString() << std::endl;

  isSetIntPoint_ = false;
}


void linElastInt::calcDMat(Matrix<Double> & dMat, const Elem* elem, const DesignElement::Type direction, double force_factor)
{
  // Bastian's stuff. If not applicable we might consider doing SIMP
  if(!GetErsatzMaterialTensor(dMat, elem, direction))
  {
    // check if we do SIMP optimization
    Double pseudo_density = elem != NULL ? GetErsatzMaterialFactor(elem) : 1.0;
    // pseudo_density can be overwritten by force_factor
    assert(!(elem != NULL && force_factor != 0.0));
    if(force_factor != 0.0) pseudo_density = force_factor;

    // do SIMP?
    if(pseudo_density != 1.0)
    {
      LOG_DBG3(forms) << GetName() << "::calcDMat(Matrix<Double>, " << (elem != NULL ? Integer(elem->elemNum) : -1)  << ") -> density=" << pseudo_density;

      // is there any chance we do bimaterial stuff at all?
      BaseMaterial* bm = elem != NULL && pseudo_density != 1.0 ? domain->GetErsatzBiMaterial(elem,  MECHANIC) : NULL;
      // Get the material tensor, in the standard case simply the tensor
      GetScaledMaterial(pseudo_density, false, bm, dMat);
    }
    else
      ptMaterial->GetTensor(dMat, MECH_STIFFNESS_TENSOR, matDataType_, subTensorType_);
  }

  //check for softening model
  if ( subTensorType_ == AXI ) {
    if (softeningPart_ == "bendingBK1" ) {
      EXCEPTION("BK1-Softening-Model for 2D not implemented");
    }

    else if (softeningPart_ == "bendingSRI" ) {
      UInt idx = dimD_-2;
      dMat[idx][idx] = 0.0;
    }

    else if (softeningPart_ == "shearBK1" ) {
      EXCEPTION("BK1-Softening-Model for 2D not implemented");
    }

    else if (softeningPart_ == "shearSRI" ) {
      for (UInt i=0; i<dimD_-2; i++) {
        for (UInt j=0; j<dimD_-2; j++) {
          dMat[i][j] = 0.0;
        }
      }

      UInt rowidx=dimD_-1;
      for (UInt i=0; i<dimD_; i++) {
        dMat[rowidx][i] = 0.0;
      }

      UInt colidx=dimD_-1;
      for (UInt i=0; i<dimD_; i++) {
        dMat[i][colidx] = 0.0;
      }

    }
  }

  else if ( subTensorType_ ==  PLANE_STRAIN || subTensorType_ ==  PLANE_STRESS ) {
    //check for softening model
    if (softeningPart_ == "bendingBK1" ) {
      EXCEPTION("BK1-Softening-Model for 2D not implemented");
    }

    else if (softeningPart_ == "bendingSRI" ) {
      UInt idx = dimD_-1;
      dMat[idx][idx] = 0.0;
    }

    else if (softeningPart_ == "shearBK1" ) {
      EXCEPTION("BK1-Softening-Model for 2D not implemented");
    }

    else if (softeningPart_ == "shearSRI" ) {
      for (UInt i=0; i<dimD_; i++) {
        for (UInt j=0; j<dimD_; j++) {
          dMat[i][j] = 0.0;
        }
      }
    }
  }

  else if (subTensorType_ == FULL ) {
    if (softeningPart_ == "bendingBK1" ) {

      if (  maxEdgeLength_ < 1e-15 ) {
        EXCEPTION("linElastInt::Calc3DMaterialMat: maxEdgeLength_=0");
      }

      Double f1, factor;
      f1 = (minEdgeLength_* minEdgeLength_);
      factor =  2.0*f1 / ( 2.0*f1 + maxEdgeLength_ * maxEdgeLength_);
      //	factor = 1.0;

      dMat[3][3] *= factor;
      dMat[4][4] *= factor;
      //	dMat[5][5] *= factor;
    }

    else if (softeningPart_ == "bendingSRI" ) {
      for (UInt i=3; i<dimD_; i++) {
        for (UInt j=3; j<dimD_; j++) {
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
      factor = f1 / ( f1 + 2.0*minEdgeLength_* minEdgeLength_ );
      //factor = 1.0;
      dMat[3][3] *= factor;
      dMat[4][4] *= factor;
      dMat[5][5]  = 0.0;

      //	dMat[5][5] *= factor;
    }

    else if (softeningPart_ == "shearSRI" ) {
      for (UInt i=0; i<3; i++) {
        for (UInt j=0; j<3; j++) {
          dMat[i][j] = 0.0;
        }
      }
    }
  }
}


void linElastInt::calcDMat(Matrix<Complex> & dMat, const Elem* elem)
{
  //check for softening model
  if (softeningModel_ != "no" )
    EXCEPTION("Softening Model just supported for real valed material data");

  ptMaterial->GetTensor(dMat,MECH_STIFFNESS_TENSOR,Global::COMPLEX,subTensorType_);
      
  Double density = elem != NULL ? GetErsatzMaterialFactor(elem) : 1.0;
  if(density != 1.0) dMat *= density;
  LOG_DBG3(forms) << GetName() << "->linElastInt::calcDMat(<Complex>, "
                  << (elem != NULL ? Integer(elem->elemNum) : -1)
                  << ") -> density=" << density;
}


// returns G - matrix for GDG (incompatible modes)
void linElastInt::calcGMat( Matrix<Double> &gMat, UInt ip,
                            Matrix<Double> &ptCoord ) {


  const UInt spaceDim = ptelem->GetDim();
  const UInt nrDofs   = getNrDofs();

  UInt actDim, actNode;
  // UInt j, k;  // TODO: Unused variables j,k

  gMat.Resize(getDimD(), nrICModes_ * nrDofs);
  gMat.Init();

  // local shape functions derived after global coords
  // (format: nrICModes x spaceDim)
  Matrix<Double> xiDx;
  ptelem->SetAnsatzFct( ansatzFct1_ );

  ptelem->GetGlobDerivShFncICModesAtIp(xiDx, ip, ptCoord, it1_.GetElem());

  for(actDim=0; actDim < spaceDim; actDim++)
    for(actNode=0; actNode < nrICModes_; actNode++)
      gMat[actDim][actNode * spaceDim + actDim] = xiDx[actNode][actDim];

  switch(spaceDim)
  {
  case 2:
    // j = 1;
    // k = 0;

    for (actNode = 0; actNode < nrICModes_; actNode++)
    {
      gMat[spaceDim][actNode * spaceDim + 1] = xiDx[actNode][0];
      gMat[spaceDim][actNode * spaceDim]     = xiDx[actNode][1];
    }

    if (isaxi_)
    {
      UInt idxtheta = getDimD();
      Vector<Double> ShpFncAtIp;
      Vector<Double> CoordAtIP;

      //compute coordinates at integration point with standad
      //shape functions
      ptelem->GetShFncAtIp(ShpFncAtIp,ip, it1_.GetElem() );
      CoordAtIP = ptCoord * ShpFncAtIp;

      //get incompatible modes hape function
      ptelem->GetShFncICModesAtIp(ShpFncAtIp,ip, it1_.GetElem() );

      // std::cout << "N=" <<   ShpFncAtIp[actNode] / CoordAtIP[0] << std::endl;
      for (actNode = 0; actNode < nrICModes_; actNode++)
        gMat[idxtheta-1][actNode * spaceDim] =
          ShpFncAtIp[actNode] / CoordAtIP[0];
    }

    break;


  case 3:
    UInt actDim=spaceDim;
    for (actNode = 0; actNode < nrICModes_; actNode++)
    {
      gMat[actDim][actNode * spaceDim + 1] = xiDx[actNode][2];
      gMat[actDim][actNode * spaceDim + 2] = xiDx[actNode][1];
    }

    actDim++;
    for (actNode = 0; actNode < nrICModes_; actNode++)
    {
      gMat[actDim][actNode * spaceDim]     = xiDx[actNode][2];
      gMat[actDim][actNode * spaceDim + 2] = xiDx[actNode][0];
    }

    actDim++;
    for (actNode = 0; actNode < nrICModes_; actNode++)
    {
      gMat[actDim][actNode * spaceDim]     = xiDx[actNode][1];
      gMat[actDim][actNode * spaceDim + 1] = xiDx[actNode][0];
    }
    break;

  }
}

void linElastInt::CalcElementMatrixICM( Matrix<Double>& elemMat,
    EntityIterator& ent1,
    EntityIterator& ent2 ) {

  // Extract pointer to reference element and get coordinates
  ExtractElemInfo( ent1 );


  // First of all, set ansatz function to element
  ptelem->SetAnsatzFct( ansatzFct1_ );

  UInt nrFncs = ptelem->GetNumFncs( ansatzFct1_ );
  const UInt nrDofs   = getNrDofs();

  double jacDet;

  Matrix<Double> bMat;
  Matrix<Double> dMat;
  Matrix<Double> dbMat;
  Matrix<Double> gMat;
  Matrix<Double> dgMat;

  Double aux, fac, *ptr1, *ptr2;

  Matrix<Double> elemMat12;  //mixed element stiffness matrix
  Matrix<Double> elemMat22;  //element stiffness matrix due to incompatible modes

  elemMat.Resize(nrFncs * nrDofs);
  elemMat.Init();
  //    elemMat11.Resize(nrFncs * nrDofs);
  //     elemMat11.Init();
  elemMat12.Resize(nrICModes_ * nrDofs, nrFncs * nrDofs);
  elemMat12.Init();
  elemMat22.Resize(nrICModes_ * nrDofs);
  elemMat22.Init();

  dbMat.Resize( getDimD(), nrFncs * nrDofs);
  dgMat.Resize( getDimD(), nrICModes_ * nrDofs);

  //get integration points
  const UInt nrIntPts = ptelem->GetNumIntPoints();
  const Vector<Double> & intWeights = ptelem->GetIntWeights();

  // **************************************************
  //  Material matrix independent of integration point
  // **************************************************
  if ( updateDMatInEveryIP_ == true ) {
    EXCEPTION("CalcElementMatrixICM does not support material dependent integration points");
  }

  // Check if material has to be rotated
  if( ptMaterial->GetCoordSys() == NULL ) {
    calcDMat(dMat, ent1.GetElem());
  }

  // Loop over all integration points
  for ( UInt actIntPt = 1; actIntPt <= nrIntPts; actIntPt++ ) {

    // Check if material has to be rotated for each integration point
    if( ptMaterial->GetCoordSys() != NULL ) {
      // Get global coordinates
      Vector<Double> * intPoints = ptelem->GetIntPoints();
      Vector<Double> globIntPoint;

      ptelem->Local2GlobalCoord(globIntPoint, intPoints[actIntPt-1],
          ptCoord_, ent1.GetElem() );
      ptMaterial->RotateTensorByPointCoord( globIntPoint,MECH_STIFFNESS_TENSOR );
      calcDMat(dMat, ent1.GetElem());
    }


    // Setup the B matrix for current integration point
    CalcBMat( bMat, actIntPt, ptCoord_ );
    //      std::cout << "bMat:\n" << bMat << std::endl;

    //incompatible modes
    calcGMat( gMat, actIntPt, ptCoord_ );

    // Compute Jacobian for integration point
    jacDet = ptelem->CalcJacobianDetAtIp( actIntPt, ptCoord_, ent1.GetElem() );

    // Perform a safety check
    if ( jacDet < 0.0 )
      EXCEPTION("BDBInt::CalcElementMatrixICM: Encountered "
                << "negative Jacobian determinant!");

    Double jacDetInv = 1.0 / jacDet;
    if ( isaxi_ ) {
      for ( UInt i=0; i<gMat.GetNumCols(); i++)
        for ( UInt j=0; j<gMat.GetNumRows()-1; j++ )
          gMat[j][i] *=  jacDetInv;
    }
    else
      gMat *=   jacDetInv;

    //      std::cout << "gMat:\n" << gMat << std::endl;

    // Special things must be done in the axi-symmetric case
    if ( isaxi_ ) {
      Vector<Double> ShpFncAtIp;
      Vector<Double> CoordAtIP;
      ptelem->GetShFncAtIp( ShpFncAtIp, actIntPt, ent1.GetElem() );

      CoordAtIP = ptCoord_ * ShpFncAtIp;
      jacDet *= 2 * PI * CoordAtIP[0];
    }

    // now we devide by the actual determinant, so that by the later multiplication
    // the determinant cancles out; with this trick, the element passes the patch test
    // see paper Taylor
    // gMat /=  sqrt(jacDet);

    // Compute the matrix product D * B and store as intermediate matrix
    dMat.Mult( bMat, dbMat );
    //    std::cout << "dMat:\n" << dMat << std::endl;

    // Compute the matrix product D * G and store as intermediate matrix
    dMat.Mult( gMat, dgMat );
    //std::cout << "dgMat:\n" << dgMat << std::endl;

    // We now compute B^T * D * B and scale it by the determinant
    // of the Jacobian and the weight of the current integration
    // point. The result is added to the element matrix.
    fac = jacDet * intWeights[actIntPt-1];
    for ( UInt k = 0; k < bMat.GetNumRows(); k++ ) {
      ptr1 =  bMat[k];
      ptr2 = dbMat[k];
      for ( UInt i = 0; i < bMat.GetNumCols(); i++ ) {
        aux = fac * ptr1[i];
        for ( UInt j = 0; j < dbMat.GetNumCols(); j++ ) {
          elemMat[i][j] += aux * ptr2[j];
        }
      }
    }

    // We now compute G^T * D * G and scale it by the determinant
    // of the Jacobian and the weight of the current integration
    // point. The result is added to the element matrix.
    fac = jacDet * intWeights[actIntPt-1];
    for ( UInt k = 0; k < gMat.GetNumRows(); k++ ) {
      ptr1 =  gMat[k];
      ptr2 = dgMat[k];
      for ( UInt i = 0; i < gMat.GetNumCols(); i++ ) {
        aux = fac * ptr1[i];
        for ( UInt j = 0; j < dgMat.GetNumCols(); j++ ) {
          elemMat22[i][j] += aux * ptr2[j];
        }
      }
    }

    // We now compute G^T * D * B and scale it by the determinant
    // of the Jacobian and the weight of the current integration
    // point. The result is added to the element matrix.
    for ( UInt k = 0; k < gMat.GetNumRows(); k++ ) {
      ptr1 =  gMat[k];
      ptr2 = dbMat[k];
      for ( UInt i = 0; i < gMat.GetNumCols(); i++ ) {
        aux = fac * ptr1[i];
        for ( UInt j = 0; j < dbMat.GetNumCols(); j++ ) {
          elemMat12[i][j] += aux * ptr2[j];
        }
      }
    }
  }

  //    std::cout << "elemMat11:\n" << elemMat11 << std::endl;
  //    std::cout << "elemMat22:\n" << elemMat22 << std::endl;
  //    std::cout << "elemMat12:\n" << elemMat12 << std::endl;

  //invert k22
  Matrix<Double> invEelemMat22;
  elemMat22.Invert(invEelemMat22);

  Matrix<Double> part1(nrICModes_ * nrDofs, nrFncs * nrDofs);
  invEelemMat22.Mult( elemMat12, part1 );

  // We now compute elemMat - k12^T * k22^(-1) * k12
  for ( UInt k = 0; k < elemMat12.GetNumRows(); k++ ) {
    ptr1 = elemMat12[k];
    ptr2 = part1[k];
    for ( UInt i = 0; i < elemMat12.GetNumCols(); i++ ) {
      aux = ptr1[i];
      for ( UInt j = 0; j < part1.GetNumCols(); j++ ) {
        elemMat[i][j] -= aux * ptr2[j];
      }
    }
  }

  //    std::cout << "softedMat:\n" << elemMat << std::endl;

}

void linElastInt::CalcElementMatrixShearBK1( Matrix<Double>& elemMat,
    EntityIterator& ent1,
    EntityIterator& ent2 )
{
  // Extract pointer to reference element and get coordinates
  ExtractElemInfo( ent1 );

  // First of all, set ansatz function to element
  ptelem->SetAnsatzFct( ansatzFct1_ );

  UInt nrFncs = ptelem->GetNumFncs( ansatzFct1_ );
  const UInt nrDofs   = getNrDofs();

  double jacDet, jacDetP, jacDetN;

  Matrix<Double> bMatP, bMatN, partElemMat;
  Matrix<Double> dMat;
  Matrix<Double> dbMat;
  Double aux, fac, *ptr1, *ptr2;

  elemMat.Resize( nrFncs * nrDofs);
  elemMat.Init();
  partElemMat.Resize( nrFncs * nrDofs);
  partElemMat.Init();

  dbMat.Resize( getDimD(), nrFncs * nrDofs);

  ptelem->GetMaxMinEdgeLength(ptCoord_,maxEdgeLength_,minEdgeLength_);

  // first compute the odd part: standard integration

  UInt nrIntPts = ptelem->GetNumIntPoints();
  const Vector<Double> & intWeights = ptelem->GetIntWeights();

  // // Check if material has to be rotated
  if( ptMaterial->GetCoordSys() == NULL )
    calcDMat(dMat, ent1.GetElem());

  if ( nrIntPts != 14 )
    EXCEPTION("For BK1 formulation we need ECONOMICAL, 4, 14 int. points");

  UInt intPtStart  = 5;
  UInt intPtEnd    = 9;
  UInt intPtOffset = 5;

  // Loop over just the intergration points located at z not equal zero:
  for ( UInt actIntPt = intPtStart; actIntPt <= intPtEnd; actIntPt++ ) {
    // Check if material has to be rotated for each integration point
    if( ptMaterial->GetCoordSys() != NULL ) {
      // Get global coordinates
      Vector<Double> * intPoints = ptelem->GetIntPoints();
      Vector<Double> globIntPoint;

      ptelem->Local2GlobalCoord(globIntPoint, intPoints[actIntPt-1],
          ptCoord_, ent1.GetElem() );
      ptMaterial->RotateTensorByPointCoord( globIntPoint,getDMaterialType() );
      calcDMat(dMat, ent1.GetElem());
    }


    // Setup the B matrix for negative z position
    CalcBMat( bMatN, actIntPt, ptCoord_ );
    CalcBMat( bMatP, actIntPt+intPtOffset, ptCoord_ );


    // Compute Jacobian for integration point
    jacDetN = ptelem->CalcJacobianDetAtIp( actIntPt, ptCoord_, ent1.GetElem() );
    jacDetP = ptelem->CalcJacobianDetAtIp( actIntPt+intPtOffset, ptCoord_,
        ent1.GetElem() );
    //     jacDet = 0.5 * (jacDetP + jacDetN);

    bMatP *=  sqrt(jacDetP * intWeights[actIntPt+intPtOffset-1]);
    bMatN *=  sqrt(jacDetN * intWeights[actIntPt-1]);

    //take the difference!!
    bMatP -= bMatN;

    // Perform a safety check
    if ( jacDetP < 0.0 )
      EXCEPTION("BDBInt::CalcElementMatrix: Encountered negative Jacobian determinant!");

    // Special things must be done in the axi-symmetric case
    if ( isaxi_ ) {
      Vector<Double> ShpFncAtIp;
      Vector<Double> CoordAtIP;
      ptelem->GetShFncAtIp( ShpFncAtIp, actIntPt, ent1.GetElem() );

      CoordAtIP = ptCoord_ * ShpFncAtIp;
      jacDet *= 2 * PI * CoordAtIP[0];
    }

    // Compute the matrix product D * B and store as intermediate matrix
    dMat.Mult( bMatP, dbMat );

    // We now compute B^T * D * B and scale it by the determinant
    // of the Jacobian and the weight of the current integration
    // point. The result is added to the element matrix.
    fac = 1.0; //jacDet * intWeights[actIntPt-1];
    for ( UInt k = 0; k < bMatP.GetNumRows(); k++ ) {
      ptr1 =  bMatP[k];
      ptr2 = dbMat[k];
      for ( UInt i = 0; i < bMatP.GetNumCols(); i++ ) {
        aux = fac * ptr1[i];
        for ( UInt j = 0; j < dbMat.GetNumCols(); j++ ) {
          partElemMat[i][j] += aux * ptr2[j];
        }
      }
    }
  }

  elemMat += partElemMat * 0.5;

  // do the reduced integration with the even part
  ptelem->SetReducedIntegration();

  nrIntPts = ptelem->GetNumIntPoints();
  if ( nrIntPts != 3 )
    EXCEPTION("For BK1 formulation we need for reduced integration SPECIAL,1");

  partElemMat.Init();

  //reduced integration: just 3 points
  intPtEnd    = 2;
  intPtOffset = 1;

  // Loop over the intergration points located at z equal zero and not equal zero:
  for ( UInt actIntPt = 1; actIntPt <= intPtEnd; actIntPt++ ) {

    // Check if material has to be rotated for each integration point
    if( ptMaterial->GetCoordSys() != NULL ) {
      // Get global coordinates
      Vector<Double> * intPoints = ptelem->GetIntPoints();
      Vector<Double> globIntPoint;

      ptelem->Local2GlobalCoord(globIntPoint, intPoints[actIntPt-1],
          ptCoord_, ent1.GetElem() );
      ptMaterial->RotateTensorByPointCoord( globIntPoint,getDMaterialType() );
      calcDMat(dMat, ent1.GetElem());
    }

    if (  actIntPt < 2 ) {
      //      if (  actIntPt < 5 ) {
      //	std::cout << "Do mult 4, z=0" << std::endl;
      // Setup the B matrix for z=0 position
      CalcBMat( bMatP, actIntPt, ptCoord_ );
      jacDet = 2.0*ptelem->CalcJacobianDetAtIp( actIntPt, ptCoord_, ent1.GetElem() );
    }
    else {
      // Setup the B matrix for negative z position
      CalcBMat( bMatN, actIntPt, ptCoord_ );
      // Setup the B matrix for positive z position
      CalcBMat( bMatP, actIntPt+intPtOffset, ptCoord_ );

      //take the sum!!
      bMatP += bMatN;

      // Compute Jacobian for integration point
      jacDetN = ptelem->CalcJacobianDetAtIp( actIntPt, ptCoord_, ent1.GetElem() );
      jacDetP = ptelem->CalcJacobianDetAtIp( actIntPt+intPtOffset, ptCoord_,
          ent1.GetElem() );
      jacDet = 0.5 * (jacDetP + jacDetN);
    }

    // Perform a safety check
    if ( jacDet < 0.0 )
      EXCEPTION("BDBInt::CalcElementMatrix: Encountered negative Jacobian determinant!");

    // Special things must be done in the axi-symmetric case
    if ( isaxi_ ) {
      Vector<Double> ShpFncAtIp;
      Vector<Double> CoordAtIP;
      ptelem->GetShFncAtIp( ShpFncAtIp, actIntPt, ent1.GetElem() );

      CoordAtIP = ptCoord_ * ShpFncAtIp;
      jacDet *= 2 * PI * CoordAtIP[0];
    }

    // Compute the matrix product D * B and store as intermediate matrix
    dMat.Mult( bMatP, dbMat );

    // We now compute B^T * D * B and scale it by the determinant
    // of the Jacobian and the weight of the current integration
    // point. The result is added to the element matrix.
    fac = jacDet * intWeights[actIntPt-1];
    for ( UInt k = 0; k < bMatP.GetNumRows(); k++ ) {
      ptr1 =  bMatP[k];
      ptr2 = dbMat[k];
      for ( UInt i = 0; i < bMatP.GetNumCols(); i++ ) {
        aux = fac * ptr1[i];
        for ( UInt j = 0; j < dbMat.GetNumCols(); j++ ) {
          partElemMat[i][j] += aux * ptr2[j];
        }
      }
    }
  }

  elemMat += partElemMat * 0.5;

  //set back to standard integration
  ptelem->SetStandardIntegration();
}


void linElastInt::SetDimensions(SubTensorType type)
{
  switch(type)
  {
  case FULL:
    dimD_      = 6;
    nrDofs_    = 3;
    nrICModes_ = 3;
    break;

  case PLANE_STRAIN:
  case PLANE_STRESS:
    dimD_      = 3;
    nrDofs_    = 2;
    nrICModes_ = 2;
    break;

  case AXI:
    dimD_      = 4;
    nrDofs_    = 2;
    nrICModes_ = 2;
    isaxi_     = true;
    break;

  case PLANE:
    EXCEPTION("SubTensorType PLANE not handled by linElastInt::SetDimensions");
    break;
  default:
    assert(false);
  }
}

linElastInt::linElastInt( BaseMaterial* matData, SubTensorType type) :
  BDBInt(matData,type)
{
  name_ = "linElastInt";
  subTensorType_ = type;
  SetDimensions(type);
}


linElastInt::~linElastInt()
{
}

bool linElastInt::GetErsatzMaterialTensor(Matrix<double>& t, const Elem* elem, DesignElement::Type direction){
  // it is possible that a region is not subject to optimization. such a region is only identified in GetErsatzMaterialTensor
  return (elem != NULL
      && domain->HasErsatzMaterialTensor()
      && domain->GetErsatzMaterial()->GetErsatzMaterialTensor(t, subTensorType_, elem, direction));
}

} // end of namespace
