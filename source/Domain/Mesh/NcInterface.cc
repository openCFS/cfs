/*
 * NcInterface.cc
 *
 *  Created on: 20.02.2013
 *      Author: jens
 */

#include "NcInterface.hh"
#include "Grid.hh"

namespace CoupledField {

  void NcInterface::NcInterface(Grid* grid, PtrParamNode nciNode) {
    UInt numIfaces;
    StdVector<SurfElem*> ifaceElems;
    StdVector<SurfElem*> masterElems;
    StdVector<SurfElem*> slaveElems;
    //StdVector<SurfElem*> surfElems;
    std::string ncRegionName;
    //RegionIdType slaveId;
    //RegionIdType masterId;
    PtrParamNode pNode;
    ParamNodeList ncIfaceNodes;

    ptGrid_ = grid;
    
    name_ = nciNode->Get("name")->As<std::string>();

    masterSurfRegion_ = ptGrid_->GetRegion().Parse(nciNode->Get("masterSide")
        ->As<std::string>());
    slaveSurfRegion_ = ptGrid_->GetRegion().Parse(nciNode->Get("slaveSide")
        ->As<std::string>());
    
    if ( masterSurfRegion_ == -1 || slaveSurfRegion_ == -1 ) {
      EXCEPTION("Cannot find master/slave regions of ncInterface '"
          << name_ << "'.");
    }
    
    ptGrid_->GetSurfElems(masterElems, masterSurfRegion_);
    ptGrid_->GetSurfElems(slaveElems, slaveSurfRegion_);

    if ( masterElems.GetSize() == 0 || slaveElems.GetSize() == 0 ) {
      EXCEPTION("Cannot find surface elements in master/slave regions of "
          << "ncInterface '" << name_ << "'.");
    }

    masterVolRegion_ = masterElems[0]->ptVolElems[0]->regionId;
    slaveVolRegion_ = slaveElems[0]->ptVolElems[0]->regionId;
    
    if (ptGrid_->GetDim() == 2) {
      intersectAlgo_ = NCI_INTERSECT_LINE;
    }
    else {
      std::string isecCalc;
      nciNode->GetValue("isecCalculation", isecCalc, ParamNode::PASS);
      if (isecCalc == "coaxi") {
        intersectAlgo_ = NCI_INTERSECT_RECT;
      }
      else {
        intersectAlgo_ = NCI_INTERSECT_POLYGON;
      }
    }

    nciNode->GetValue("tolAbs", tolAbs_, ParamNode::PASS);
    nciNode->GetValue("tolRel", tolRel_, ParamNode::PASS);

    PtrParamNode rotNode = nciNode->Get("rotation", ParamNode::PASS);
    if (rotNode) {
      rotNode->GetValue("coordSysId", coordSysId_, ParamNode::PASS);
      rotNode->GetValue("rpm", rpm_, ParamNode::PASS);
      motionType_ = NCI_ROTATION;
    }
    else {
      rpm_ = 0.0;
      motionType_ = NCI_NOT_MOVING;
    }

    if (rpm_ != 0.0) {
      coplanar_ = false;
    }
    else {
      UInt numMasterElems = masterElems.GetSize(),
           numSlaveElems = slaveElems.GetSize();
      
      ifaceElems.Reserve(numMasterElems + numSlaveElems);
      
      for (UInt i = 0; i < numMasterElems; ++i) {
        ifaceElems.Push_back(masterElems[i]);
      }

      for (UInt i = 0; i < numSlaveElems; ++i) {
        ifaceElems.Push_back(slaveElems[i]);
      }

      coplanar_ = IsSurfacePlanar(ifaceElems);
    }

    // create a surface region for the intersection elements
    // TODO: Do we really need this?
    region_ = ptGrid_->AddSurfaceRegion(name_);

    // initialize nodeOffsets if there is a moving region
    if (motionType_ != NCI_NOT_MOVING) {
      StdVector<UInt> nodeNums;
      Vector<Double> nullOffsets;
      ptGrid_->GetNodesByRegion(nodeNums, slaveVolRegion_);
      nullOffsets.Resize(nodeNums.GetSize()*ptGrid_->GetDim(), 0.0);
      ptGrid_->SetNodeOffset(nodeNums, nullOffsets);
    }

    // Calculate the intersection
    UpdateIntersection();

  }

} // namespace CoupledField
