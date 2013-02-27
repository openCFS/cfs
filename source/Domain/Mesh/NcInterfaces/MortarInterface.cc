/*
 * MortarInterface.cc
 *
 *  Created on: 23.02.2013
 *      Author: jens
 */

#include "MortarInterface.hh"
#include "Domain/ElemMapping/SurfElem.hh"

namespace CoupledField {

MortarInterface::MortarInterface(Grid* grid, PtrParamNode nciNode) :
  BaseNcInterface(grid)
{
  StdVector<SurfElem*> masterElems;
  StdVector<SurfElem*> slaveElems;

  SetName( nciNode->Get("name")->As<std::string>() );

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

  std::string isecCalc;
  nciNode->GetValue("isecCalculation", isecCalc, ParamNode::PASS);
  if (ptGrid_->GetDim() == 2) {
    intersectAlgo_ = NCI_INTERSECT_LINE;
    if ( isecCalc != "" ) {
      WARN("ncInterface '" << name_ << "': intersection algorithm '"
           << isecCalc << "' is not valid in a 2D simulation.");
    }
  }
  else {
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
    SetRotation(
        rotNode->Get("coordSysId", ParamNode::PASS)->As<std::string>(),
        rotNode->Get("rpm", ParamNode::PASS)->As<Double>());
  }

  // initialize nodeOffsets if there is a moving region
  if ( isMoving_ ) {
    StdVector<UInt> nodeNums;
    Vector<Double> nullOffsets;
    ptGrid_->GetNodesByRegion(nodeNums, slaveVolRegion_);
    nullOffsets.Resize(nodeNums.GetSize()*ptGrid_->GetDim(), 0.0);
    ptGrid_->SetNodeOffset(nodeNums, nullOffsets);
  } else {
    // Calculate the intersection, if interface is stationary.
    // If there is motion, UpdateInterface will be called by MoveInterface.
    UpdateInterface();
  }
}

MortarInterface::~MortarInterface() : ~BaseNcInterface() {}

void MortarInterface::UpdateInterface() {
  StdVector<SurfElem*> masterElems;
  StdVector<SurfElem*> slaveElems;
  StdVector<SurfElem*> ifaceElems;
  StdVector<NcSurfElem*> ncElems;
  StdVector<SurfElem*> ncElemsHelper;
  StdVector<UInt> masterNodes;
  StdVector<UInt> ncElemIds;
  StdVector<UInt> newNodes;
  StdVector<std::string> listNodeNames;
  std::string newNodesName = name_ + "_nodes";

  ptGrid_->ClearRegion(region_);
  ptGrid_->GetListNodeNames(listNodeNames);
  if ( listNodeNames.Find(newNodesName) != -1 ) {
    ptGrid_->DeleteNamedNodes(newNodesName);
  }

  ptGrid_->GetSurfElems(masterElems, masterSurfRegion_);
  ptGrid_->GetSurfElems(slaveElems, slaveSurfRegion_);

  // check if interface is coplanar
  if ( isMoving_ && coordSys_ && coordSys_->GetDofName(2) == "phi" ) {
    // rotational motion cannot create coplanar interfaces
    coplanar_ = false;
  } else {
    UInt numMasterElems = masterElems.GetSize(),
        numSlaveElems = slaveElems.GetSize();

    ifaceElems.Reserve(numMasterElems + numSlaveElems);

    for (UInt i = 0; i < numMasterElems; ++i) {
      ifaceElems.Push_back(masterElems[i]);
    }

    for (UInt i = 0; i < numSlaveElems; ++i) {
      ifaceElems.Push_back(slaveElems[i]);
    }

    coplanar_ = ptGrid_->IsSurfacePlanar(ifaceElems);
  }

  switch (intersectAlgo_) {
    case NCI_INTERSECT_LINE:
      for (UInt i = 0; i < masterElems.GetSize(); ++i) {
        for (UInt j = 0; j < slaveElems.GetSize(); ++j) {
          SideOnSide(masterElems[i], slaveElems[j], coplanar_,
              needsUpdate, ncElems, newNodes );
        }
      }
      break;
    case NCI_INTERSECT_RECT:
      if (!coplanar_) {
        EXCEPTION("Only coplanar interfaces are supported with coaxial "
            << "rectangle algorithm.");
      }
      for (UInt i = 0; i < masterElems.GetSize(); ++i) {
        for(UInt j = 0; j < slaveElems.GetSize(); ++j) {
          SurfElem* m_el = masterElems[i];
          SurfElem* s_el = slaveElems[j];
          if( (m_el->type != Elem::ET_QUAD4 )
              || s_el->type != Elem::ET_QUAD4 )
          {
            EXCEPTION("Only quadrilaterals can be intersected with coaxial "
                << "rectangle algorithm.");
          }

          if(RectangleOnRectangle(masterElems[i], slaveElems[j],
              ncIf.needsUpdate, ncElems))
          {
            /*LOG_DBG3(grid) << "Intersection between "
                << masterElems[i]->elemNum << " and "
                << slaveElems[j]->elemNum << std::endl;*/
          }
        }
      }
      break;
    case NCI_INTERSECT_POLYGON:
      for (UInt i = 0; i < masterElems.GetSize(); ++i) {
        for (UInt j = 0; j < slaveElems.GetSize(); ++j) {
          SurfElem* m_el = masterElems[i];
          SurfElem* s_el = slaveElems[j];
          if( (m_el->type != Elem::ET_QUAD4 )
              || s_el->type != Elem::ET_QUAD4 )
          {
            EXCEPTION("Only triangles and quadrilaterals can be intersected"
                << " with polygon algorithm.");
          }

          if (PolygonOnPolygon(masterElems[i], slaveElems[j], coplanar_,
              ncIf.needsUpdate, ncElems,
              tolAbs_, tolRel_))
          {
            /*LOG_DBG3(grid) << "Intersection between "
                << masterElems[i]->elemNum << " and "
                << slaveElems[j]->elemNum << std::endl;*/
          }
        }
      }
      break;
    default:
      EXCEPTION("Intersection algorithm is not implemented");
      break;
  }

  if(ncElems.GetSize() > 0)
  {
    UInt numElems = ncElems.GetSize();

    ncElemsHelper.Resize(numElems);

    for(UInt i=0; i<numElems; ++i)
    {
      ncElemsHelper[i] = ncElems[i];
    }

    ptGrid_->AddSurfaceElems(region_, ncElemsHelper, ncElemIds);

    if ( newNodes.GetSize() > 0 )
    {
      ptGrid_->AddNamedNodes(newNodesName, newNodes);
    }
  }
  else {
    EXCEPTION("No intersection elements were computed for non-conforming"
        << " interface '" << name_
        << "'. Please check your mesh file.");
  }
}

} /* namespace CoupledField */
