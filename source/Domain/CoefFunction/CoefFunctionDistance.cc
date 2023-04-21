/*
 * CoefFunctionDistance.cc
 *
 *  Created on: 30 Nov 2022
 *      Author: Dominik Mayrhofer
 */

#include "CoefFunctionDistance.hh"

#include "DataInOut/Logging/LogConfigurator.hh"

DEFINE_LOG(coeffctdistance, "coeffctdistance")

namespace CoupledField {
  CoefFunctionDistance::CoefFunctionDistance(shared_ptr<BaseFeFunction> feFct,
                                                        std::string surf1,
                                                        std::string surf2,
                                                        std::string vol,
                                                        bool useSurfaceMidpoints) : CoefFunction() {
    dimType_ = VECTOR; // should be determined from input
    dependType_ = SOLUTION; // should be determined from input
    isAnalytic_ = false;
    isComplex_  = false;

    // Use NN search on elements to determine the distance between a source surface and a target reion (nodes or surface midpoints can be used for te evaluation)
    feFct_ = feFct;
    surf1_ = surf1;
    surf2_ = surf2;
    vol_ = vol;
    useSurfaceMidpoints_ = useSurfaceMidpoints;

    volRegion_ = feFct_->GetGrid()->GetRegionId(vol);
    surfRegion1_ = feFct_->GetGrid()->GetRegionId(surf1);
    surfRegion2_ = feFct_->GetGrid()->GetRegionId(surf2);

  }

  CoefFunctionDistance::~CoefFunctionDistance() {
    ;
  }

  void CoefFunctionDistance::GetScalar(Double& scal, const LocPointMapped& surfLpm ) {
    LOG_DBG(coeffctdistance) << "+++++ CoefFunctionDistance::GetScalar ++++++";
    
    // Check that we are actually dealing with a surface lpm
    assert(surfLpm.isSurface);

    // get the updated query point
    Vector<Double> queryPoint;
    //surfLpm.GetShapeMap()->GetGlobMidPoint(queryPoint);
    const SurfElem* srcElem;
    srcElem = surfLpm.GetShapeMap()->GetSurfElem();

    StdVector<UInt> sourceNodes;
    StdVector<UInt> sourceNodesCorr;
    StdVector<Vector<Double>> sourceNodeCoords;
    Vector<Double> sourceVecSum;
    // get the global node numbers of the current element
    feFct_->GetGrid()->GetElemNodes(sourceNodes,srcElem->elemNum);
    // get the nodal coordinates belonging to the node numbers
    feFct_->GetGrid()->GetNodeCoordinates(sourceNodeCoords, sourceNodes, true); // always use the updated geometry!
    //Vector<Double> curMidPoint = targetElems[i]->GetShape().midPointCoord;
    
    // calculate the barycenter
    sourceVecSum.Resize(sourceNodeCoords.GetSize());
    sourceVecSum.Init(0);
    queryPoint.Resize(sourceNodeCoords.GetSize());
    queryPoint.Init(0);
    for (UInt u = 0; u < sourceNodeCoords.GetSize(); ++u) {
      sourceVecSum+=sourceNodeCoords[u];
    }
    for (UInt u = 0; u < sourceVecSum.GetSize(); ++u) {
      queryPoint[u] = sourceVecSum[u]/(sourceNodeCoords.GetSize());
    }


    // get the correct target region
    RegionIdType surfaceRegion = surfLpm.ptEl->regionId;

    RegionIdType targetRegion = surfaceRegion==surfRegion1_  ? surfRegion2_ : surfRegion1_; // target surface for the NN search

    StdVector<SurfElem*> targetElems;
    feFct_->GetGrid()->GetSurfElems(targetElems, targetRegion);
    
    // Delete all points before setting them
    ResetPoints();

    // loop over the target and add points
    // switch if we want to use just the node list or calculate the barycenters and use them
    if( useSurfaceMidpoints_ ) {
      // use surface midpoints (calculate barycenters and use them for the NN search)
      StdVector<UInt> targetNodes;
      StdVector<UInt> targetNodesCorr;
      StdVector<Vector<Double>> nodeCoords;
      Vector<Double> vecSum;
      Vector<Double> curMidPoint;

      // loop over all elements
      LOG_DBG(coeffctdistance) << "Looping over all elements using surface midpoints";
      for (UInt i = 0; i < targetElems.GetSize(); ++i) {
        // get the global node numbers of the current element
        feFct_->GetGrid()->GetElemNodes(targetNodes,targetElems[i]->elemNum);
        // get the nodal coordinates belonging to the node numbers
        feFct_->GetGrid()->GetNodeCoordinates(nodeCoords, targetNodes, true); // always use the updated geometry!
        
        // calculate the barycenter
        vecSum.Resize(nodeCoords.GetSize());
        vecSum.Init(0);
        curMidPoint.Resize(nodeCoords.GetSize());
        curMidPoint.Init(0);
        for (UInt u = 0; u < nodeCoords.GetSize(); ++u) {
          vecSum+=nodeCoords[u];
        }
        for (UInt u = 0; u < vecSum.GetSize(); ++u) {
          curMidPoint[u] = vecSum[u]/(nodeCoords.GetSize());
        }
        AddPoint(curMidPoint);
        // For debugging purposes
        LOG_DBG2(coeffctdistance) << "curMidPoint[0] " << i << ": " << curMidPoint[0]<< std::endl;
        LOG_DBG2(coeffctdistance) << "curMidPoint[1]: " << curMidPoint[1]<< std::endl;
      }
    } else {
      // use the nodes of the surface list

      // get the nodes based on the region itself
      StdVector<UInt> targetNodes;
      feFct_->GetGrid()->GetNodesByRegion(targetNodes, targetRegion);

      // get the nodal coordinates belonging to the node numbers
      StdVector<Vector<Double>> nodeCoords;
      feFct_->GetGrid()->GetNodeCoordinates(nodeCoords, targetNodes, true); // always use the updated geometry!

      // loop over all nodes and add the points
      for (UInt u = 0; u < targetNodes.GetSize(); ++u) {
        AddPoint(nodeCoords[u]);
      }
    }
    // For debugging purposes
    LOG_DBG(coeffctdistance) << "queryPoint[0]: " << queryPoint[0]<< std::endl;

    // use the query point to evaluate the distance vector
    QueryPoint(scal, queryPoint);

    // For debugging purposes
    LOG_DBG(coeffctdistance) << "Calculated distance " << scal << std::endl;
  }

  void CoefFunctionDistance::AddPoint(Vector<Double> targetPoint) {
    targetPoints_.Push_back(targetPoint);
  }

  void CoefFunctionDistance::QueryPoint(Double& dist, Vector<Double> queryPoint) {
    Vector<Double> diffVec;
    Double curDiff;
    Double minDiff = std::numeric_limits<double>::infinity(); // init to infinity so that we do not run into problems with large distances
    Vector<Double> distVec;
    for (UInt i = 0; i < targetPoints_.GetSize(); ++i) {
      diffVec = targetPoints_[i]-queryPoint;
      curDiff = diffVec.NormL2();
      if (curDiff<=minDiff) {
        minDiff = curDiff;
      }
      //std::cout << "Diff vec " << i << ": " << diffVec << std::endl;
    }
    dist = minDiff;
  }

  void CoefFunctionDistance::ResetPoints() {
    targetPoints_.Clear();
  }

}
