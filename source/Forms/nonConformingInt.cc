// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "DataInOut/Logging/cfslog.hh"
#include "Domain/grid.hh"
#include "Domain/ncElem.hh"
#include "Domain/domain.hh"
#include "Elements/elements_header.hh"
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
      Error("NonConformingInt can only handle NCElems!",__FILE__,__LINE__);
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
    stiffMat.Resize( shapeFuncsSurfSize, shapeFuncsLGSize );
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

    ptGrid->GetElemNodesCoord( coordMat, surfElem[0]->connect);
    surfElem[0]->ptElem->Global2LocalCoords(localIntPointsSurf,
                                            globalIntPoints,
                                            coordMat);
    
    LOG_DBG3(ncint) << "Local ip Surf: " << surfElem[0]->elemNum
                    << ": " << std::endl << localIntPointsSurf;

    ptGrid->GetElemNodesCoord( coordMat, surfElem[1]->connect);
    surfElem[1]->ptElem->Global2LocalCoords(localIntPointsLG,
                                            globalIntPoints,
                                            coordMat);

    LOG_DBG3(ncint) << "Local ip Lagrange: " << surfElem[1]->elemNum
                    << ": " << std::endl << localIntPointsLG;

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
        
      for(UInt s=0; s < shapeFuncsSurfSize; s++)
      {
        for(UInt l=0; l < shapeFuncsLGSize; l++)
        {
          // stiffMat[s][l] += jacDet * intWeights[i] * 
          //                   shapeFuncsSurf[shapeFuncsSurfSize - s - 1] *
          //                   shapeFuncsLG[shapeFuncsLGSize - l - 1 ];

          // This destroyed the circular pattern for the l2d example.
          // TODO: strieben - Take a closer look!
          stiffMat[s][l] += jacDet * intWeights[i] * 
                            shapeFuncsSurf[s] *
                            shapeFuncsLG[l];
        }
      }
    }
  }
      

}
