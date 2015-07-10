// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <stddef.h>
#include <cmath>
#include <iostream>
#include <string>

#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"
#include "Domain/domain.hh"
#include "Domain/elem.hh"
#include "Domain/entityList.hh"
#include "Domain/grid.hh"
#include "Domain/ncElem.hh"
#include "Domain/surfElem.hh"
#include "Elements/basefe.hh"
#include "General/environment.hh"
#include "General/exception.hh"
#include "MatVec/exprt/xpr1.hh"
#include "MatVec/matrix.hh"
#include "MatVec/vector.hh"
#include "Utils/StdVector.hh"
#include "nonConformingInt.hh"


namespace CoupledField {

  DECLARE_LOG(ncint)
  DEFINE_LOG(ncint, "forms.nonconformingint")

  NonConformingInt::NonConformingInt( UInt dofsPerNode, bool isAxi,
                                      bool coordUpdate ) {

    name_ = "NonConformingInt";
    isaxi_ = isAxi;
    dofs_ = dofsPerNode;
    coordUpdate_ = coordUpdate;

  }

  NonConformingInt::~NonConformingInt() {
  }

  void NonConformingInt::CalcElementMatrix(  Matrix<Double>& stiffMat,
                                             EntityIterator& ent1,
                                             EntityIterator& ent2 ) {

    SurfElem *surfElem[2];
    Matrix<Double> coordMat;
    Matrix<Double> coordMatLG, coordMatSurf;
    Vector< Double >* intPoints;
    Matrix<Double> globalIntPoints, localIntPointsLG, localIntPointsSurf;
    Vector< Double > point;
    Vector< Double > shapeFuncsLG, shapeFuncsSurf;
    UInt dim;
    Double jacDet;
    Vector<Double> intWeights;
    UInt shapeFuncsSurfSize, shapeFuncsLGSize, numIntPoints;

    Grid* ptGrid = domain->GetGrid();

    const NCElem *elem = dynamic_cast< const NCElem* >(ent1.GetElem());
    if(elem == NULL)
    {
      EXCEPTION("NonConformingInt can only handle NCElems!");
    }

    ptGrid->GetElemNodesCoord( ptCoord_,
                               elem->connect,
                               coordUpdate_ );

    LOG_DBG3(ncint) << "On NC-elem #" << elem->elemNum
                    << " with connect = " << elem->connect.Serialize()
                    << " with coordinates \n"
                    << ptCoord_;

    surfElem[0] = elem->ptSurfParent;
    surfElem[1] = elem->ptLagrangeParent;

    dim = ptGrid->GetDim();

    // calculate surface normal
    ptGrid->CalcSurfNormal(normal_, *surfElem[1]);
    normal_ *= (Double) surfElem[1]->normalSign;

    surfElem[0]->ptElem->SetAnsatzFct(ansatzFct1_);
    shapeFuncsSurfSize = surfElem[0]->ptElem->GetNumFncs(ansatzFct1_);
    surfElem[1]->ptElem->SetAnsatzFct(ansatzFct2_);
    shapeFuncsLGSize = surfElem[1]->ptElem->GetNumFncs(ansatzFct2_);

    // obtain number of integration points for the NCElem
    numIntPoints = elem->ptElem->GetNumIntPoints();
    intPoints = elem->ptElem->GetIntPoints();
    intWeights = elem->ptElem->GetIntWeights();

    // resize element matrix to size
    // num_degrees(SurfParent) x num_degrees(LagrangeParent)
    stiffMat.Resize( shapeFuncsSurfSize * dofs_, shapeFuncsLGSize * dofs_ );
    stiffMat.Init();

    // local matrix:
    //
    //           slave side dofs -->
    //
    //   master  * * * * * ...
    //    side   * * * * * ...
    //    dofs   * * * * * ...
    //     |     * * * * * ...
    //     |     . . . . . ...
    //     v     . . . . . ...
    //

    // get corner coordinates of NCElem and resize
    // the intPoints matrices to dim x numIntPoints
    globalIntPoints.Resize(dim, numIntPoints);
    localIntPointsSurf.Resize(dim-1, numIntPoints);
    localIntPointsLG.Resize(dim-1, numIntPoints);

    point.Resize(dim);

    // Project local int points from ncelem to surfelem
    // and lagrange elem
    //
    // x---------------o------o----x
    //                 |      |
    //            x----o------o----x
    //                 |      |
    //            x----o------o---------------x

    for(UInt i=0; i < numIntPoints; i++)
    {
      elem->ptElem->Local2GlobalCoord(point,
                                      intPoints[i],
                                      ptCoord_,
                                      NULL);
      for(UInt j=0; j < dim; j++)
      {
        globalIntPoints[j][i] = point[j];
        LOG_DBG3(ncint) << "Global ip "<< i << ": " << point[j];
      }
    }

    ptGrid->GetElemNodesCoord( coordMat, surfElem[1]->connect, coordUpdate_);
    surfElem[1]->ptElem->Global2LocalCoords(localIntPointsLG,
                                            globalIntPoints,
                                            coordMat);

    LOG_DBG3(ncint) << "Local ip Lagrange: " << surfElem[1]->elemNum
                    << ": " << std::endl << localIntPointsLG;

    ptGrid->GetElemNodesCoord( coordMat, surfElem[0]->connect, coordUpdate_);
    if (!ptGrid->IsNcInterfaceCoplanar(elem->regionId)) { // projection is necessary
      Double scale, dist, sign;
      Vector<Double> normal, p0, line;

      // calculate surface normal of master element
      ptGrid->CalcSurfNormal(normal, *surfElem[0], coordUpdate_);

      // get first node of master element
      p0.Resize(dim);
      for (UInt i = 0; i < dim; ++i)
        p0[i] = coordMat[i][0];

      // project integration points onto master element
      for (UInt i = 0; i < numIntPoints; ++i) {
        // get global coordinates of i-th integration point
        for (UInt j = 0; j < dim; ++j)
          point[j] = globalIntPoints[j][i];

        line = point - p0;

        // calculate distance of int. point to master element
        normal.Inner(line, dist);

        // determine orientation of normal on slave element
        normal_.Inner(line, sign);
        sign /= fabs(sign);

        if(std::isnan(sign))
          continue;
        
        // scale the distance for the normal projection
        normal_.Inner(normal, scale);

        // do the projection
        for (UInt j = 0; j < dim; ++j) {
          globalIntPoints[j][i] -= sign * normal_[j] * fabs(dist) * fabs(scale);
        }
      }
    }
    surfElem[0]->ptElem->Global2LocalCoords(localIntPointsSurf,
                                            globalIntPoints,
                                            coordMat);

    LOG_DBG3(ncint) << "Local ip Surf: " << surfElem[0]->elemNum
                    << ": " << std::endl << localIntPointsSurf;

    point.Resize(dim-1);

    for(UInt i=0; i < numIntPoints; i++)
    {
      for(UInt j=0; j < dim-1; j++)
        point[j] = localIntPointsSurf[j][i];

      try
      {
          surfElem[0]->ptElem->GetShFnc(shapeFuncsSurf, point, surfElem[0]);
      }
      catch (Exception& ex)
      {
          std::cerr << "In ncElem " << elem->elemNum
                    << ", slaveElem " << surfElem[1]->elemNum
                    << ", masterElem " << surfElem[0]->elemNum << ": "
                    << std::endl;
          std::cerr << ex.what();
          return;
      }

      for(UInt j=0; j < dim-1; j++)
        point[j] = localIntPointsLG[j][i];

      try
      {
          surfElem[1]->ptElem->GetShFnc(shapeFuncsLG, point, surfElem[1]);
      }
      catch (Exception& ex)
      {
          std::cerr << "In ncElem " << elem->elemNum
                    << ", slaveElem " << surfElem[1]->elemNum
                    << ", masterElem " << surfElem[0]->elemNum << ": "
                    << std::endl;
          std::cerr << ex.what();
          return;
      }

      jacDet = std::fabs(elem->ptElem->CalcJacobianDet(intPoints[i],
                                                       ptCoord_, elem));

      jacDet = -jacDet;

      if (isaxi_) {
        for(UInt s=0; s < shapeFuncsSurfSize; s++) {
          for(UInt l=0; l < shapeFuncsLGSize; l++) {
            for(UInt dof=0; dof < dofs_; dof++) {
            // TODO: strieben - Take a closer look!
            Double radius = globalIntPoints[0][i];
              stiffMat[s*dofs_+dof][l*dofs_+dof] += 
                2.0 * PI * radius
                * jacDet * intWeights[i]
                * shapeFuncsSurf[s] * shapeFuncsLG[l];
            }
          }
        }
      }
      else {
        for(UInt s=0; s < shapeFuncsSurfSize; s++) {
          for(UInt l=0; l < shapeFuncsLGSize; l++) {
            for(UInt dof=0; dof < dofs_; dof++) {
              // stiffMat[s][l] += jacDet * intWeights[i] *
              //                   shapeFuncsSurf[shapeFuncsSurfSize - s - 1] *
              //                   shapeFuncsLG[shapeFuncsLGSize - l - 1 ];

              // This destroyed the circular pattern for the l2d example.
              // TODO: strieben - Take a closer look!
              stiffMat[s*dofs_+dof][l*dofs_+dof] += jacDet * intWeights[i] *
                shapeFuncsSurf[s] *
                shapeFuncsLG[l];
            }
          }
        }
      }
    }

    LOG_DBG2(ncint) << "element matrix of NCElem #" << elem->elemNum
                    << ":" << std::endl << stiffMat << std::endl;
  }


}
