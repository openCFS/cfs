/*
 * MortarInterface.cc
 *
 *  Created on: 23.02.2013
 *      Author: jens
 */

#include "MortarInterface.hh"
#include "PolygonIterators.hh"
#include "Domain/Domain.hh"
#include "Domain/ElemMapping/SurfElem.hh"
#include "Domain/ElemMapping/EntityLists.hh"
#include "Domain/CoordinateSystems/CoordSystem.hh"
#include "Driver/TransientDriver.hh"
#include "Utils/StdVector.hh"
#include "MatVec/Vector.hh"
#include "Utils/mathParser/mathParser.hh"

#include <sstream>
#include <boost/shared_ptr.hpp>

namespace CoupledField {

MortarInterface::MortarInterface(Grid* grid, PtrParamNode nciNode) :
  BaseNcInterface(grid),
  isCoplanar_(false),
  isMoving_(false),
  moveMaster_(false),
  exportToGrid_(true),
  geoWarn_(true),
  coordSysId_(""),
  coordSys_(NULL),
  mParser_(NULL),
  tolAbs_(1e-12),
  tolRel_(1e-4),
  region_(NO_REGION_ID),
  isReset_(false)
{
  name_ = nciNode->Get("name")->As<std::string>();
  elemList_->SetName(name_);
  
  StdVector<SurfElem*> masterElems;
  StdVector<SurfElem*> slaveElems;

  masterSurfRegion_ = ptGrid_->GetRegion().Parse(nciNode->Get("masterSide")
      ->As<std::string>(), NO_REGION_ID);
  slaveSurfRegion_ = ptGrid_->GetRegion().Parse(nciNode->Get("slaveSide")
      ->As<std::string>(), NO_REGION_ID);

  if (masterSurfRegion_ == NO_REGION_ID || slaveSurfRegion_ == NO_REGION_ID) {
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
  nciNode->GetValue("intersectionMethod", isecCalc, ParamNode::PASS);
  if (ptGrid_->GetDim() == 2) {
    intersectAlgo_ = NCI_INTERSECT_LINE;
    if ( !isecCalc.empty() ) {
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
  
  nciNode->GetValue("storeIntegrationGrid", exportToGrid_, ParamNode::PASS);

  nciNode->GetValue("geometryWarnings", geoWarn_, ParamNode::PASS);
  
  PtrParamNode motionNode = nciNode->Get("rotation", ParamNode::PASS);
  if (motionNode) {
    SetRotation( motionNode->Get("coordSysId")->As<std::string>(),
                 motionNode->Get("rpm")->As<Double>());
    moveMaster_ = (motionNode->Get("movingSide", ParamNode::INSERT)
                    ->As<std::string>() == "master");
  }
  
  motionNode = nciNode->Get("generalMotion", ParamNode::PASS);
  if (motionNode) {
    std::string coordSysId = "default";
    StdVector<std::string> displaceExpr;
    
    if (motionNode->Has("displace3")) {
      displaceExpr.Resize(3);
      displaceExpr[2] = motionNode->Get("displace3")->As<std::string>();
    } else {
      displaceExpr.Resize(2);
    }
    displaceExpr[0] = motionNode->Get("displace1")->As<std::string>();
    displaceExpr[1] = motionNode->Get("displace2")->As<std::string>();
    
    motionNode->GetValue("coordSysId", coordSysId, ParamNode::INSERT);
    
    SetMotion(displaceExpr, coordSysId);

    moveMaster_ = (motionNode->Get("movingSide", ParamNode::INSERT)
                    ->As<std::string>() == "master");
  }

  //if ( !isMoving_ ) {
    // Calculate the intersection, if interface is stationary.
    // If there is motion, UpdateInterface will be called by TransientDriver.
    UpdateInterface();
  //}
}

MortarInterface::~MortarInterface() {
  ptGrid_ = NULL;
  coordSys_ = NULL;
  
  if ( mParser_ ) {
    for ( UInt i=0; i<3; ++i ) {
      if (mphOffset_[i] != MathParser::GLOB_HANDLER)
        mParser_->ReleaseHandle(mphOffset_[i]);
    }
    mParser_ = NULL;
  }
}


void MortarInterface::SetRotation(const std::string &coordSysId, Double rpm) {
  if ( rpm == 0.0 ) return;
  
  coordSys_ = domain->GetCoordSystem(coordSysId);
  if ( coordSys_->GetDofName(2) != "phi" ) {
    EXCEPTION("For a rotating ncInterface the coordinate system must be "
        << "either polar (in 2D) of cylindrical (3D).");
  }
  
  std::ostringstream sstr("");
  sstr << std::setprecision(std::numeric_limits<Double>::digits10)
      // 1 rpm = 360° / 60 s = 6 °/s
       << std::scientific << (6.0*rpm) << "*t";
  offsetExpr_.Resize(coordSys_->GetDim());
  offsetExpr_[1] = sstr.str();

  SetMotion(offsetExpr_, coordSysId);
}

void MortarInterface::SetMotion(const StdVector<std::string> &offsetExpr,
                                const std::string &coordSysId)
{
  offsetExpr_ = offsetExpr;
  coordSys_ = domain->GetCoordSystem(coordSysId);
  coordSysId_ = coordSysId;
  
  UInt dim = coordSys_->GetDim(); 
  if ( dim != offsetExpr.GetSize() ) {
    EXCEPTION("You must provide exactly as many math parser expressions as "
              << "there are dimensions in the coordinate system.")
  }

  mParser_ = domain->GetMathParser();
  for ( UInt i=0; i<3; ++i ) {
    if ( i < dim && !offsetExpr_[i].empty() ) {
      mphOffset_[i] = mParser_->GetNewHandle(true);
      mParser_->SetExpr(mphOffset_[i], offsetExpr_[i]);
      if ( mParser_->IsExprVariable( mphOffset_[i], "t") ) 
        isMoving_ = true;
    }
    else {
      mphOffset_[i] = MathParser::GLOB_HANDLER;
    }
  }
  
  if ( isMoving_ ) {
    // initialize nodeOffsets if there is a moving region
    StdVector<UInt> nodeNums;
    Vector<Double> nullOffsets;
    ptGrid_->GetNodesByRegion(nodeNums, slaveVolRegion_);
    nullOffsets.Resize(nodeNums.GetSize()*dim, 0.0);
    ptGrid_->SetNodeOffset(nodeNums, nullOffsets);
  } else {
    WARN("You supplied constant expressions as time-dependent "
         << "displacements for moving ncInterface '" << name_
         << "'. Interface is assumed stationary.");
  }
}

void MortarInterface::MoveInterface() {
  
  if ( !isMoving_ ) return;
  
  const bool useGlobalCoords = (coordSysId_ == "default");
  UInt dim = coordSys_->GetDim(), numNodes;
  StdVector<UInt> nodeNums;
  Vector<Double> coordOrig, coordTmp, coordNew, nodeOffsets;

  if (moveMaster_) {
    ptGrid_->GetNodesByRegion(nodeNums, masterVolRegion_);
  }
  else {
    ptGrid_->GetNodesByRegion(nodeNums, slaveVolRegion_);
  }

  numNodes = nodeNums.GetSize();
  nodeOffsets.Resize(numNodes * dim);

  if ( useGlobalCoords ) {
    for (UInt i = 0; i < numNodes; ++i) {
      for (UInt j = 0; j < dim; ++j) {
        nodeOffsets[i*dim+j] = mParser_->Eval(mphOffset_[j]);
      }
    }
  }
  else {
    for (UInt i = 0; i < numNodes; ++i) {
      ptGrid_->GetNodeCoordinate(coordOrig, nodeNums[i], false);

      coordSys_->Global2LocalCoord(coordTmp, coordOrig);
      for (UInt j = 0; j < dim; ++j) {
        coordTmp[j] += mParser_->Eval(mphOffset_[j]);
      }
      coordSys_->Local2GlobalCoord(coordNew, coordTmp);

      nodeOffsets[i * dim    ] = coordNew[0] - coordOrig[0];
      nodeOffsets[i * dim + 1] = coordNew[1] - coordOrig[1];
      if (dim == 3)
        nodeOffsets[i*dim + 2] = coordNew[2] - coordOrig[2];
    }
  }

  ptGrid_->SetNodeOffset(nodeNums, nodeOffsets);
}

void MortarInterface::ResetInterface(){

  if ( !isMoving_ && elemList_->GetSize() > 0 ) return;
  StdVector<std::string> listNodeNames;
  std::string newNodesName = name_ + "_nodes";

  elemList_->Clear();
  if ( region_ != NO_REGION_ID ) {
    ptGrid_->ClearRegion(region_);
  }
  ptGrid_->GetListNodeNames(listNodeNames);
  if ( listNodeNames.Find(newNodesName) != -1 ) {
    ptGrid_->DeleteNamedNodes(newNodesName);
  }
  isReset_ = true;
}

void MortarInterface::UpdateInterface() {
  
  if ( !isMoving_ && (elemList_->GetSize() > 0 || isReset_) ) return;
  
  isReset_ = false;
  StdVector<SurfElem*> masterElems;
  StdVector<SurfElem*> slaveElems;
  StdVector<SurfElem*> ifaceElems;
  StdVector<SurfElem*> ncElemsHelper;
  StdVector<UInt> masterNodes;
  StdVector<UInt> ncElemIds;
  StdVector<UInt> newNodes;
  StdVector<std::string> listNodeNames;
  std::string newNodesName = name_ + "_nodes";

  MoveInterface();

  ptGrid_->GetSurfElems(masterElems, masterSurfRegion_);
  ptGrid_->GetSurfElems(slaveElems, slaveSurfRegion_);

  // check if interface is coplanar
  if ( isMoving_ && coordSys_ && coordSys_->GetDofName(2) == "phi" ) {
    // rotational motion cannot create coplanar interfaces
    isCoplanar_ = false;
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

    isCoplanar_ = ptGrid_->IsSurfacePlanar(ifaceElems);
  }

  switch (intersectAlgo_) {
    case NCI_INTERSECT_LINE:
#pragma omp parallel for
      for (UInt i = 0; i < masterElems.GetSize(); ++i) {
        for (UInt j = 0; j < slaveElems.GetSize(); ++j) {
          IntersectLines(masterElems[i], slaveElems[j], newNodes );
        }
      }
      break;
    case NCI_INTERSECT_RECT:
      if (!isCoplanar_) {
        EXCEPTION("Only coplanar interfaces are supported with coaxial "
            << "rectangle algorithm.");
      }
#pragma omp parallel for
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

          if(IntersectRects( masterElems[i], slaveElems[j], newNodes))
          {
            /*LOG_DBG3(grid) << "Intersection between "
                << masterElems[i]->elemNum << " and "
                << slaveElems[j]->elemNum << std::endl;*/
          }
        }
      }
      break;
    case NCI_INTERSECT_POLYGON:
#pragma omp parallel for
      for (UInt i = 0; i < masterElems.GetSize(); ++i) {
        for (UInt j = 0; j < slaveElems.GetSize(); ++j) {
          SurfElem* m_el = masterElems[i];
          SurfElem* s_el = slaveElems[j];
          if ( (m_el->type != Elem::ET_QUAD4 && m_el->type != Elem::ET_QUAD8
                && m_el->type != Elem::ET_QUAD9 && m_el->type != Elem::ET_TRIA3
                && m_el->type != Elem::ET_TRIA6)
              || (s_el->type != Elem::ET_QUAD4 && s_el->type != Elem::ET_QUAD8
                  && s_el->type != Elem::ET_QUAD9 && s_el->type != Elem::ET_TRIA3
                  && s_el->type != Elem::ET_TRIA6) )
          {
            EXCEPTION("Only triangles and quadrilaterals can be intersected"
                << " with polygon algorithm.");
          }

          if (IntersectPolygons( masterElems[i], slaveElems[j], newNodes ))
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

  if ( newNodes.GetSize() > 0 ) {
    ptGrid_->AddNamedNodes(newNodesName, newNodes);
  }

  UInt numElems = elemList_->GetSize();
  
  if( numElems > 0 ) {
    UpdateIntegrators();
    
    if ( exportToGrid_ ) {
      if ( region_ == NO_REGION_ID ) {
        region_ = ptGrid_->AddSurfaceRegion(name_);
      }

      ncElemsHelper.Resize(numElems);

      for ( UInt i=0; i<numElems; ++i ) {
        // We need to make explicit copies of the NcSurfElems, because the
        // Grid deletes all its elements when it gets destroyed.
        ncElemsHelper[i] = new SurfElem(*(elemList_->GetSurfElem(i)));
      }

      ptGrid_->AddSurfaceElems(region_, ncElemsHelper, ncElemIds);

      for ( UInt i=0; i<numElems; ++i ) {
        std::map<std::string, UInt>::iterator it = ptGrid_->entityDim_.find(name_);

        if( it != ptGrid_->entityDim_.end() ) {
          if( it->second != Elem::shapes[elemList_->GetSurfElem(i)->type].dim ) {
            EXCEPTION( "Region '" << name_
                << "' contains elements of different dimensions!");
          }
        } else {
          ptGrid_->entityDim_[name_] =
              Elem::shapes[elemList_->GetSurfElem(i)->type].dim;
        }
      }
    }
  }
  else {
    EXCEPTION("No intersection elements were computed for non-conforming"
        << " interface '" << name_
        << "'. Please check your mesh file.");
  }
}

/****************************************************************************
 **
 ** IntersectLines
 **
 **   computes the local coordinates of the overlap of the master and slave
 **   element with respect to the master side in the order of the
 **   orientation of the slave side element. It pushes back the intersection
 **   element to elemList.
 **
 ** Input Parameters:
 **   ifaceElem1:  Master Side
 **   ifaceElem2:  Slave Side
 **
 ** Output Parameters:
 **   elemList: the found intersection NCElems will be pushed
 **                     back to this vector
 **
 */

bool MortarInterface::IntersectLines( SurfElem *ifaceElem1,
                                      SurfElem *ifaceElem2,
                                      StdVector<UInt> &newNodes )
{
  // c0, c1, d0 and d1 are the endpoints of the two line elements
  //
  //           d0 x-----------+--------------x d1
  // c0 x---------+-----------x c1

  Vector<Double> c0, c1, d0, d1, tmp; // endpoint-coordinates of the two lines
  Vector<Double> diff0, diff1;
  Vector<Double> s, t;
  Vector<Double> normal;
  StdVector<UInt> connect2;
  Double dist, fac;
  UInt nodenum_c0, nodenum_c1, nodenum_d0, nodenum_d1;
  Double relativeElemVol;

  s.Resize(2);
  t.Resize(2);
  connect2.Resize(2);


  // Get coordinates of the endpoints
  nodenum_c0 = ifaceElem1->connect[0];
  nodenum_c1 = ifaceElem1->connect[1];
  nodenum_d0 = ifaceElem2->connect[0];
  nodenum_d1 = ifaceElem2->connect[1];
  ptGrid_->GetNodeCoordinate(c0, nodenum_c0, isMoving_);
  ptGrid_->GetNodeCoordinate(c1, nodenum_c1, isMoving_);
  ptGrid_->GetNodeCoordinate(d0, nodenum_d0, isMoving_);
  ptGrid_->GetNodeCoordinate(d1, nodenum_d1, isMoving_);
  
  // Project master nodes onto slave element, if interface is not coplanar
  if ( !isCoplanar_ ) {
    shared_ptr<ElemShapeMap> sm = ptGrid_->GetElemShapeMap(ifaceElem2, isMoving_);
    LocPoint lp = Elem::shapes[ifaceElem2->type].midPointCoord;

    // compute maximal allowed distance as sum of lengths of both lines
    tmp = c1 - c0;
    Double maxDist = tmp.NormL2();
    tmp = d1 - d0;
    maxDist += tmp.NormL2();
    // compute normal vector of slave element
    sm->CalcNormal(normal, lp);
    // compute distance of c0 to plane of slave element
    tmp = c0 - d0;
    fac = normal.Inner(tmp);
    // make sure that distance does not exceed maximum distance
    if (fabs(fac) > maxDist)
      return false;
    // do the projection if necessary
    if ( fabs(fac) > 1e-12 )
    {
      c0 -= normal * fac;
      // add new node later, if necessary
      nodenum_c0 = 0;
    }

    // do the same for c1
    tmp = c1 - d0;
    fac = normal.Inner(tmp);
    if ( fabs(fac) > 1e-12 )
    {
    c1 -= normal * fac;
      nodenum_c1 = 0;
    }
  }

  // Compute and normalize vector from c0 to c1.
  // This becomes the new x-unit vector.
  diff1 = c1 - c0;
  dist = diff1.NormL2();
  // check if both elements were perpendicular
  if ( dist < 1e-12 )
    return false;
  fac = 1.0 / dist;
  diff1 *= fac;

  // Compute x1 coordinate of line2 in respect to line1.
  diff0 = d0 - c0;
  diff0.Inner(diff1, s[0]);
  s[0] *= fac;

  // Compute x2 coordinate of line2 in respect to line1.
  diff0 = d1 - c0;
  diff0.Inner(diff1, s[1]);
  s[1] *= fac;

  // Bring line2's endpoints into ascending order.
  if(s[1] < s[0])
  {
    t[0] = s[1];
    t[1] = s[0];
    connect2[0] = nodenum_d1;
    connect2[1] = nodenum_d0;
  }
  else
  {
    t[0] = s[0];
    t[1] = s[1];
    connect2[0] = nodenum_d0;
    connect2[1] = nodenum_d1;
  }

  // Check if an intersection between line1 and line2 exists.
  if(t[0] >= 1.0)
    return false;

  if(t[1] <= 0.0)
    return false;

  shared_ptr<MortarNcSurfElem> ncElem(new MortarNcSurfElem());
  ncElem->connect.Resize(2);

  relativeElemVol = t[1] - t[0];

  // If an intersection exists, we must distinguish 4 different cases.
  if(t[0] <= 0)
  {
    if (nodenum_c0 == 0) {
      // create a new node for the projection of c0 onto d
      Vector<Double> new_node;
      new_node.Resize(3);
      new_node[0] = c0[0];
      new_node[1] = c0[1];
      if (c0.GetSize() == 2)
        new_node[2] = 0.0;
      else
        new_node[2] = c0[2];
      ptGrid_->AddNode(new_node, nodenum_c0);
      newNodes.Push_back(nodenum_c0);
    }
    ncElem->connect[0] = nodenum_c0;

    if(t[1] >= 1)
    {
      // connect2[0] x--------|--------------------|-----x connect2[1]
      //                   c0 x--------------------x c1

      if (nodenum_c1 == 0) {
        // create a new node for the projection of c1 onto d
        Vector<Double> new_node;
        new_node.Resize(3);
        new_node[0] = c1[0];
        new_node[1] = c1[1];
        if (c1.GetSize() == 2)
          new_node[2] = 0.0;
        else
          new_node[2] = c1[2];
        ptGrid_->AddNode(new_node, nodenum_c1);
        newNodes.Push_back(nodenum_c1);
      }
      ncElem->connect[1] = nodenum_c1;
      relativeElemVol = 1;
    }
    else
    {
      // connect2[0] x--------|---------x connect2[1]
      //                   c0 x---------|-----------x c1

      relativeElemVol = t[1];
      ncElem->connect[1] = connect2[1];
    }

  }
  else
  {
    ncElem->connect[0] = connect2[0];

    if(t[1] >= 1)
    {
      // connect2[0] x----------------|------x connect2[1]
      //      c0 x---|----------------x c1

      if (nodenum_c1 == 0) {
        // create a new node for the projection of c1 onto d
        Vector<Double> new_node;
        new_node.Resize(3);
        new_node[0] = c1[0];
        new_node[1] = c1[1];
        if (c1.GetSize() == 2)
          new_node[2] = 0.0;
        else
          new_node[2] = c1[2];
        ptGrid_->AddNode(new_node, nodenum_c1);
        newNodes.Push_back(nodenum_c1);
      }
      ncElem->connect[1] = nodenum_c1;
      relativeElemVol = 1-t[0];
    }
    else
    {
      // connect2[0] x------------x connect2[1]
      //      c0 x---|------------|---x c1

      ncElem->connect[1] = connect2[1];
    }
  }

  if (relativeElemVol < 1e-6) {
    if (geoWarn_) {
      WARN("Rejecting ncElem due to a relative volume of " << relativeElemVol
          << std::endl
          << "  for intersection of elements " << ifaceElem1->elemNum
          // << " (" << region_.ToString(ifaceElem1->regionId) << ")"
          << " and " << ifaceElem2->elemNum);
      // << " (" << this->region_.ToString(ifaceElem2->regionId) << ") ");
    }
    return false;
  }

  // In case of a curved interface store the projected master element.
  // This is needed for coordinate transform of integration points.
  shared_ptr<SurfElem> projMaster( new SurfElem() );
  if ( !isCoplanar_ ) {
    projMaster->type = Elem::ET_LINE2;
    projMaster->connect.Resize(2);
    
    if (nodenum_c0 == 0) {
      // create a new node for the projection of c0 onto d
      Vector<Double> new_node;
      new_node.Resize(3);
      new_node[0] = c0[0];
      new_node[1] = c0[1];
      if (c0.GetSize() == 2)
        new_node[2] = 0.0;
      else
        new_node[2] = c0[2];
      ptGrid_->AddNode(new_node, nodenum_c0);
      newNodes.Push_back(nodenum_c0);
    }
    projMaster->connect[0] = nodenum_c0;
    
    if (nodenum_c1 == 0) {
      // create a new node for the projection of c0 onto d
      Vector<Double> new_node;
      new_node.Resize(3);
      new_node[0] = c1[0];
      new_node[1] = c1[1];
      if (c1.GetSize() == 2)
        new_node[2] = 0.0;
      else
        new_node[2] = c1[2];
      ptGrid_->AddNode(new_node, nodenum_c1);
      newNodes.Push_back(nodenum_c1);
    }
    projMaster->connect[1] = nodenum_c1;
    
    // projection might be identical to ncElem itself
    if ( projMaster->connect != ncElem->connect ) {
      ncElem->projectedMaster = projMaster;
    } // else case: NULL pointer signals projectedMaster==ncElem
  }
  
  ncElem->type = Elem::ET_LINE2;
  ncElem->ptMaster = ifaceElem1;
  ncElem->ptSlave = ifaceElem2;

  elemList_->AddElement(ncElem);

  return true;
}

bool MortarInterface::IntersectRects( SurfElem *ifaceElem1,
                                      SurfElem *ifaceElem2,
                                      StdVector<UInt> &newNodes )
{
  Vector<Double> c0, c1, c2, d0, d1, d2;
  Vector<Double> diffS, diffX, diffY, diffX2;
  Vector<Double> s, t;
  StdVector<UInt> connect2;
  Double distX, distY, distX2, facX, facY, r;
  UInt nodeNr;
  // Introduce a tolerance to account for roundoff errors during the calculation of
  // normed new x base vector. 
  Double tol_r;

  s.Resize(4);
  t.Resize(4);
  connect2.Resize(4);

  // The meaning of the points c0, c1, c2, d0 and d2
  // is as follows:
  //                x------------------x d2
  //                |                  |
  //                |                  |
  // c2 x-----------+----------x       |
  //    |           |          |       |
  //    |           |          |       |
  //    |        d0 x----------+-------x d1
  //    |                      |
  //    |                      |
  // c0 x----------------------x c1


  // Get coordinates of the endpoints
  ptGrid_->GetNodeCoordinate(c0, ifaceElem1->connect[0], isMoving_);
  ptGrid_->GetNodeCoordinate(c1, ifaceElem1->connect[1], isMoving_);
  ptGrid_->GetNodeCoordinate(c2, ifaceElem1->connect[3], isMoving_);
  ptGrid_->GetNodeCoordinate(d0, ifaceElem2->connect[0], isMoving_);
  ptGrid_->GetNodeCoordinate(d1, ifaceElem2->connect[1], isMoving_);
  ptGrid_->GetNodeCoordinate(d2, ifaceElem2->connect[2], isMoving_);

  // Compute and normalize vector from c0 to c1.
  // This becomes the new x-unit vector.
  diffX = c1 - c0;
  distX = diffX.NormL2();
  facX = 1.0 / distX;
  diffX *= facX;

  // Compute and normalize vector from c0 to c2
  // This becomes the new y-unit vector.
  diffY = c2 - c0;
  distY = diffY.NormL2();
  facY = 1.0 / distY;
  diffY *= facY;

  // Now compute vector from c0 to d0 and project
  // the result onto the new x- and y-axis.
  diffS = d0 - c0;
  diffS.Inner(diffX, s[0]);
  diffS.Inner(diffY, s[1]);
  s[0] *= facX;
  s[1] *= facY;

  // Now compute vector from c0 to d2 and project
  // the result onto the new x- and y-axis.
  diffS = d2 - c0;
  diffS.Inner(diffX, s[2]);
  diffS.Inner(diffY, s[3]);
  s[2] *= facX;
  s[3] *= facY;

  // Determine the orientation of the second rectangle
  // to make sure that the edges which connect c0 and c1
  // are parallel to the edges which connect d0 and d1.
  diffX2 = d1 - d0;
  distX2 = diffX2.NormL2();
  diffX.Inner(diffX2, r);

  // Set the tolerance for determining if the edges 
  // mentioned in the last comment are parallel.
  tol_r = distX2 < distX ? distX2 / 10 : distX / 10;

  // Bring the x- and y-coordinates of the intersection
  // into an order, where the smaller coordinates come
  // first.
  if(s[2] < s[0])
  {
    t[0] = s[2];
    t[2] = s[0];

    if(s[3] < s[1])
    {
      t[1] = s[3];
      t[3] = s[1];
      connect2[0] = ifaceElem2->connect[2];
      connect2[2] = ifaceElem2->connect[0];
      if (fabs(r) < tol_r) {
        connect2[1] = ifaceElem2->connect[1];
        connect2[3] = ifaceElem2->connect[3];
      }
      else {
        connect2[1] = ifaceElem2->connect[3];
        connect2[3] = ifaceElem2->connect[1];
      }
    }
    else
    {
      t[1] = s[1];
      t[3] = s[3];
      connect2[1] = ifaceElem2->connect[0];
      connect2[3] = ifaceElem2->connect[2];
      if (fabs(r) < tol_r) {
        connect2[0] = ifaceElem2->connect[3];
        connect2[2] = ifaceElem2->connect[1];
      }
      else {
        connect2[0] = ifaceElem2->connect[1];
        connect2[2] = ifaceElem2->connect[3];
      }
    }

  }
  else
  {
    t[0] = s[0];
    t[2] = s[2];

    if(s[3] < s[1])
    {
      t[1] = s[3];
      t[3] = s[1];
      connect2[1] = ifaceElem2->connect[2];
      connect2[3] = ifaceElem2->connect[0];
      if (fabs(r) < tol_r) {
        connect2[0] = ifaceElem2->connect[1];
        connect2[2] = ifaceElem2->connect[3];
      }
      else {
        connect2[0] = ifaceElem2->connect[3];
        connect2[2] = ifaceElem2->connect[1];
      }
    }
    else
    {
      t[1] = s[1];
      t[3] = s[3];
      connect2[0] = ifaceElem2->connect[0];
      connect2[2] = ifaceElem2->connect[2];
      if (fabs(r) < tol_r) {
        connect2[1] = ifaceElem2->connect[3];
        connect2[3] = ifaceElem2->connect[1];
      }
      else {
        connect2[1] = ifaceElem2->connect[1];
        connect2[3] = ifaceElem2->connect[3];
      }
    }
  }

  // Check if an intersection between rectangle1
  // and rectangle2 exists.
  if(t[0] >= 1.0)
    return false;

  if(t[2] <= 0.0)
    return false;

  if(t[1] >= 1.0)
    return false;

  if(t[3] <= 0.0)
    return false;

  shared_ptr<MortarNcSurfElem> ncElem(new MortarNcSurfElem());
  ncElem->connect.Resize(4);

  diffX *= distX;
  diffY *= distY;

  // If an intersection actually exist, we eventually
  // have to compute the intersection points.
  // There exist 16 different cases how two axiparallel
  // rectangles can intersect each other.

  Vector<Double> tmp;

  if(t[0] <= 0)
  {
    if(t[2] >= 1)
    {
      if(t[1] <= 0)
      {
        ncElem->connect[0] = ifaceElem1->connect[0];
        ncElem->connect[1] = ifaceElem1->connect[1];

        if(t[3] >= 1)
        {
          ncElem->connect[2] = ifaceElem1->connect[2];
          ncElem->connect[3] = ifaceElem1->connect[3];
        }
        else
        {
          tmp = c0 + diffX     + diffY*t[3];
          ptGrid_->AddNode(tmp, nodeNr);
          newNodes.Push_back(nodeNr);
          ncElem->connect[2] = nodeNr;
          tmp = c0 + diffX*0.0 + diffY*t[3];
          ptGrid_->AddNode(tmp, nodeNr);
          newNodes.Push_back(nodeNr);
          ncElem->connect[3] = nodeNr;
        }

      }
      else
      {
        if(t[3] >= 1)
        {
          tmp = c0 + diffX*0.0 + diffY*t[1];
          ptGrid_->AddNode(tmp, nodeNr);
          newNodes.Push_back(nodeNr);
          ncElem->connect[0] = nodeNr;
          tmp = c0 + diffX     + diffY*t[1];
          ptGrid_->AddNode(tmp, nodeNr);
          newNodes.Push_back(nodeNr);
          ncElem->connect[1] = nodeNr;
          ncElem->connect[2] = ifaceElem1->connect[2];
          ncElem->connect[3] = ifaceElem1->connect[3];
        }
        else
        {
          tmp = c0 + diffX*0.0 + diffY*t[1];
          ptGrid_->AddNode(tmp, nodeNr);
          newNodes.Push_back(nodeNr);
          ncElem->connect[0] = nodeNr;
          tmp = c0 + diffX     + diffY*t[1];
          ptGrid_->AddNode(tmp, nodeNr);
          newNodes.Push_back(nodeNr);
          ncElem->connect[1] = nodeNr;
          tmp = c0 + diffX     + diffY*t[3];
          ptGrid_->AddNode(tmp, nodeNr);
          newNodes.Push_back(nodeNr);
          ncElem->connect[2] = nodeNr;
          tmp = c0 + diffX*0.0 + diffY*t[3];
          ptGrid_->AddNode(tmp, nodeNr);
          newNodes.Push_back(nodeNr);
          ncElem->connect[3] = nodeNr;
        }
      }
    }
    else
    {
      if(t[1] <= 0)
      {
        if(t[3] >= 1)
        {
          ncElem->connect[0] = ifaceElem1->connect[0];
          tmp = c0 + diffX*t[2] + diffY*0.0;
          ptGrid_->AddNode(tmp, nodeNr);
          newNodes.Push_back(nodeNr);
          ncElem->connect[1] = nodeNr;
          tmp = c0 + diffX*t[2] + diffY;
          ptGrid_->AddNode(tmp, nodeNr);
          newNodes.Push_back(nodeNr);
          ncElem->connect[2] = nodeNr;
          ncElem->connect[3] = ifaceElem1->connect[3];
        }
        else
        {
          ncElem->connect[0] = ifaceElem1->connect[0];
          tmp = c0 + diffX*t[2] + diffY*0.0;
          ptGrid_->AddNode(tmp, nodeNr);
          newNodes.Push_back(nodeNr);
          ncElem->connect[1] = nodeNr;
          ncElem->connect[2] = connect2[2];
          tmp = c0 + diffX*0.0  + diffY*t[3];
          ptGrid_->AddNode(tmp, nodeNr);
          newNodes.Push_back(nodeNr);
          ncElem->connect[3] = nodeNr;
        }

      }
      else
      {
        if(t[3] >= 1)
        {
          tmp = c0 + diffX*0.0  + diffY*t[1];
          ptGrid_->AddNode(tmp, nodeNr);
          newNodes.Push_back(nodeNr);
          ncElem->connect[0] = nodeNr;
          ncElem->connect[1] = connect2[1];
          tmp = c0 + diffX*t[2] + diffY;
          ptGrid_->AddNode(tmp, nodeNr);
          newNodes.Push_back(nodeNr);
          ncElem->connect[2] = nodeNr;
          ncElem->connect[3] = ifaceElem1->connect[3];
        }
        else
        {
          tmp = c0 + diffX*0.0  + diffY*t[1];
          ptGrid_->AddNode(tmp, nodeNr);
          newNodes.Push_back(nodeNr);
          ncElem->connect[0] = nodeNr;
          ncElem->connect[1] = connect2[1];
          ncElem->connect[2] = connect2[2];
          tmp = c0 + diffX*0.0  + diffY*t[3];
          ptGrid_->AddNode(tmp, nodeNr);
          newNodes.Push_back(nodeNr);
          ncElem->connect[3] = nodeNr;
        }
      }
    }
  }
  else
  {
    if(t[2] >= 1)
    {
      if(t[1] <= 0)
      {
        if(t[3] >= 1)
        {
          tmp = c0 + diffX*t[0] + diffY*0.0;
          ptGrid_->AddNode(tmp, nodeNr);
          newNodes.Push_back(nodeNr);
          ncElem->connect[0] = nodeNr;
          ncElem->connect[1] = ifaceElem1->connect[1];
          ncElem->connect[2] = ifaceElem1->connect[2];
          tmp = c0 + diffX*t[0] + diffY;
          ptGrid_->AddNode(tmp, nodeNr);
          newNodes.Push_back(nodeNr);
          ncElem->connect[3] = nodeNr;
        }
        else
        {
          tmp = c0 + diffX*t[0] + diffY*0.0;
          ptGrid_->AddNode(tmp, nodeNr);
          newNodes.Push_back(nodeNr);
          ncElem->connect[0] = nodeNr;
          ncElem->connect[1] = ifaceElem1->connect[1];
          tmp = c0 + diffX      + diffY*t[3];
          ptGrid_->AddNode(tmp, nodeNr);
          newNodes.Push_back(nodeNr);
          ncElem->connect[2] = nodeNr;
          ncElem->connect[3] = connect2[3];
        }

      }
      else
      {
        if(t[3] >= 1)
        {
          ncElem->connect[0] = connect2[0];
          tmp = c0 + diffX      + diffY*t[1];
          ptGrid_->AddNode(tmp, nodeNr);
          newNodes.Push_back(nodeNr);
          ncElem->connect[1] = nodeNr;
          ncElem->connect[2] = ifaceElem1->connect[2];
          tmp = c0 + diffX*t[0] + diffY;
          ptGrid_->AddNode(tmp, nodeNr);
          newNodes.Push_back(nodeNr);
          ncElem->connect[3] = nodeNr;
        }
        else
        {
          ncElem->connect[0] = connect2[0];
          tmp = c0 + diffX      + diffY*t[1];
          ptGrid_->AddNode(tmp, nodeNr);
          newNodes.Push_back(nodeNr);
          ncElem->connect[1] = nodeNr;
          tmp = c0 + diffX      + diffY*t[3];
          ptGrid_->AddNode(tmp, nodeNr);
          newNodes.Push_back(nodeNr);
          ncElem->connect[2] = nodeNr;
          ncElem->connect[3] = connect2[3];
        }
      }
    }
    else
    {
      if(t[1] <= 0)
      {
        if(t[3] >= 1)
        {
          tmp = c0 + diffX*t[0] + diffY*0.0;
          ptGrid_->AddNode(tmp, nodeNr);
          newNodes.Push_back(nodeNr);
          ncElem->connect[0] = nodeNr;
          tmp = c0 + diffX*t[2] + diffY*0.0;
          ptGrid_->AddNode(tmp, nodeNr);
          newNodes.Push_back(nodeNr);
          ncElem->connect[1] = nodeNr;
          tmp = c0 + diffX*t[2] + diffY;
          ptGrid_->AddNode(tmp, nodeNr);
          newNodes.Push_back(nodeNr);
          ncElem->connect[2] = nodeNr;
          tmp = c0 + diffX*t[0] + diffY;
          ptGrid_->AddNode(tmp, nodeNr);
          newNodes.Push_back(nodeNr);
          ncElem->connect[3] = nodeNr;
        }
        else
        {
          tmp = c0 + diffX*t[0] + diffY*0.0;
          ptGrid_->AddNode(tmp, nodeNr);
          newNodes.Push_back(nodeNr);
          ncElem->connect[0] = nodeNr;
          tmp = c0 + diffX*t[2] + diffY*0.0;
          ptGrid_->AddNode(tmp, nodeNr);
          newNodes.Push_back(nodeNr);
          ncElem->connect[1] = nodeNr;
          ncElem->connect[2] = connect2[2];
          ncElem->connect[3] = connect2[3];
        }

      }
      else
      {
        if(t[3] >= 1)
        {
          ncElem->connect[0] = connect2[0];
          ncElem->connect[1] = connect2[1];
          tmp = c0 + diffX*t[2] + diffY;
          ptGrid_->AddNode(tmp, nodeNr);
          newNodes.Push_back(nodeNr);
          ncElem->connect[2] = nodeNr;
          tmp = c0 + diffX*t[0] + diffY;
          ptGrid_->AddNode(tmp, nodeNr);
          newNodes.Push_back(nodeNr);
          ncElem->connect[3] = nodeNr;
        }
        else
        {
          ncElem->connect[0] = connect2[0];
          ncElem->connect[1] = connect2[1];
          ncElem->connect[2] = connect2[2];
          ncElem->connect[3] = connect2[3];
        }
      }
    }
  }
  
  ncElem->type = Elem::ET_QUAD4;
  ncElem->ptMaster = ifaceElem1;
  ncElem->ptSlave = ifaceElem2;

  elemList_->AddElement(ncElem);
  
  return true;
}

bool MortarInterface::IntersectPolygons( SurfElem *ifElem1, SurfElem *ifElem2,
                                        StdVector<UInt> &newNodes )
{
  UInt i, j, n, p1Size, p2Size;
  StdVector< Vector<Double> > p1, p2, r;

  switch(ifElem1->type) {
    case Elem::ET_TRIA3:
    case Elem::ET_TRIA6:
      p1Size = 3;
      break;

    case Elem::ET_QUAD4:
    case Elem::ET_QUAD8:
    case Elem::ET_QUAD9:
      p1Size = 4;
      break;

    default:
      EXCEPTION("First argument to PolygonOnPolygon may not be of type '"
                << Elem::feType.ToString(ifElem1->type) << "!");
      break;
  }
  
  p1.Resize(p1Size);
  for (i = 0; i < p1Size; ++i)
    ptGrid_->GetNodeCoordinate(p1[i], ifElem1->connect[i], isMoving_);

  switch(ifElem2->type) {
    case Elem::ET_TRIA3:
    case Elem::ET_TRIA6:
      p2Size = 3;
      break;

    case Elem::ET_QUAD4:
    case Elem::ET_QUAD8:
    case Elem::ET_QUAD9:
      p2Size = 4;
      break;

    default:
      EXCEPTION("Second argument to PolygonOnPolygon may not be of type '"
          << Elem::feType.ToString(ifElem2->type) << "'!");
      break;
  }

  p2.Resize(p2Size);
  for (i = 0; i < p2Size; ++i)
    ptGrid_->GetNodeCoordinate(p2[i], ifElem2->connect[i], isMoving_);

  if (CutPolys(p1, p2, isCoplanar_, r))
  {
    StdVector<MortarNcSurfElem*> newElems;
    
    UInt nodeStart = newNodes.GetSize();
    UInt elemStart = TriangulatePoly(r, newElems, newNodes);

    // store projected master element
    shared_ptr<SurfElem> projMaster;
    
    if ( !isCoplanar_ ) {
      Vector<Double> nodeCoord;
      projMaster.reset( new SurfElem() );
      projMaster->type = ifElem1->type;
      projMaster->connect.Resize(p1Size);
      
      n = newNodes.GetSize();
      for ( i = 0; i < p1Size; ++i ) {
        Vector<Double> & p1Coord = p1[i];
        for ( j = nodeStart; j < n; ++j ) {
          ptGrid_->GetNodeCoordinate( nodeCoord, newNodes[j], isMoving_ );
          if ( nodeCoord == p1Coord ) {
            projMaster->connect[i] = newNodes[j];
            break;
          }
        }
        if ( j == n ) {
          ptGrid_->AddNode(p1[i], projMaster->connect[i]);
          newNodes.Push_back(projMaster->connect[i]);
        }
      }
      
      if ( newElems.GetSize() - elemStart == 1 ) { // there is only 1 element
        if ( projMaster->connect == newElems[elemStart]->connect ) {
          // projected master is identical to Mortar element
          projMaster.reset(); // Reset to NULL pointer
        }
      }
    }

    n = newElems.GetSize();
    for ( i = elemStart; i < n; ++i)
    {
      shared_ptr<MortarNcSurfElem> newElem(newElems[i]);
      newElem->ptMaster = ifElem1;
      newElem->ptSlave = ifElem2;
      newElem->projectedMaster = projMaster;
      elemList_->AddElement(newElem);
    }

    return true;
  }
  return false;
}

MortarInterface::LineIntersectType MortarInterface::CutLines(const Vector<Double> &a,
                                   const Vector<Double> &b, const Vector<Double> &c,
                                   const Vector<Double> &d, Vector<Double> &e) const
{
  Double l1, l2;
  Vector<Double> v1, v2, temp;

#ifdef CHECK_INDEX
  if ((a.GetSize() != 3) || (b.GetSize() != 3) || (c.GetSize() != 3) ||
      (d.GetSize() != 3)) {
    EXCEPTION("Points must be given as 3D coordinates");
    return INTERSECT_NONE;
  }
#endif

  v1.Resize(3);
  v2.Resize(3);
  e.Resize(3);

  // calculate vectors of both lines
  v1 = b - a;
  v2 = d - c;
  // calculate lengths of both lines
  l1 = v1.NormL2();
  l2 = v2.NormL2();

  if ( v1.Collinear(v2) ) { // lines are parallel
    // if line from a to d is also parallel then lines may intersect
    e = d - a;
    if ( v1.Collinear(e) ) {
      Double l_ac, l_ad, l_bc, l_bd;
      // calculate distances between points
      temp = (c - a);
      l_ac = temp.NormL2();

      temp = (d - a);
      l_ad = temp.NormL2();

      temp = (c - b);
      l_bc = temp.NormL2();

      temp = (d - b);
      l_bd = temp.NormL2();

      // does a lie on [c,d]?
      if (fabs(l_ac + l_ad - l2) < tolAbs_) {
        e = a;
        if (l_ac < tolAbs_) // is a=c?
          return INTERSECT_A_EQ_C;
        if (l_ad < tolAbs_) { // is a=d?
          // usually a cut at d wins over a cut at a,
          // but c might be an endpoint here, too
          if (fabs(l_ac + l_bc - l1) < tolAbs_) {
            e = c;
            return INTERSECT_IN_C;
          }
          return INTERSECT_IN_D;
        }
        // does c lie on [a,b]?
        if (fabs(l_ac + l_bc - l1) < tolAbs_)
          return INTERSECT_A_AND_C; // intersection is [a,c]
        
        return INTERSECT_IN_A;
      }
      if (fabs(l_bc + l_bd - l2) < tolAbs_) {
        e = b;
        return INTERSECT_IN_B;
      }
      if (fabs(l_ac + l_bc - l1) < tolAbs_) {
        e = c;
        return INTERSECT_IN_C;
      }
      if (fabs(l_ad + l_bd - l1) < tolAbs_) {
        e = d;
        return INTERSECT_IN_D;
      }

      // both lines lie on one infinite virtual line
      if (l_bc < l_ac) { // we only consider [a,inf) for type 0
        // return point closest to a
        if (l_ad < l_ac)
          e = d;
        else
          e = c;
        return INTERSECT_OUTSIDE;
      }
      return INTERSECT_NONE; // line a->b does not point to line [c,d]
    }
    return INTERSECT_NONE; // lines are parallel but do not intersect
  }

  /* At this point we know that the lines are not parallel,
   * so compute intersection.
   *
   * a + h * v1 = c + k * v2
   *
   * This is a system with 2 unknowns (h,k) and 3 equations. Compute k1
   * from equations 1 and 2, k2 from equations 1 and 3, and k3 from
   * equations 2 and 3. Depending on the orientation of the lines in 3D
   * space, up to two values out of (k1,k2,k3) may be undefined, because
   * the denominator is zero. Therefore we need to select the right k.
   */

  Double h, k, k1 = 0.0, k2 = 0.0, k3 = 0.0, denom1, denom2, denom3;

  denom1 = v1[1] * v2[0] - v1[0] * v2[1];
  if (fabs(denom1) > tolAbs_)
    k1 = (v1[0] * (c[1] - a[1]) + v1[1] * (a[0] - c[0])) / denom1;
  denom2 = v1[2] * v2[0] - v1[0] * v2[2];
  if (fabs(denom2) > tolAbs_)
    k2 = (v1[0] * (c[2] - a[2]) + v1[2] * (a[0] - c[0])) / denom2;
  denom3 = v1[2] * v2[1] - v1[1] * v2[2];
  if (fabs(denom3) > tolAbs_)
    k3 = (v1[1] * (c[2] - a[2]) + v1[2] * (a[1] - c[1])) /denom3;

  // If this system has no solution, lines do not intersect.
  if ((fabs(denom1) <= tolAbs_)
      && (fabs(denom2) <= tolAbs_)
      && (fabs(denom3) <= tolAbs_))
    return INTERSECT_NONE;

  /* TODO: jens
   * This check makes no sense for 3 k's. Maybe add a check based on the
   * standard deviation of k.
   */
  /*if ((fabs(denom1) > tolAbs_)
      && (fabs(denom2) > tolAbs_)) {
    if (fabs(k1 - k2) > tolRel_)
      return INTERSECT_NONE;
  }*/

  // select the right one out of (k1,k2,k3)
  if (fabs(denom1) > tolAbs_)
    k = k1;
  else if (fabs(denom2) > tolAbs_)
    k = k2;
  else
    k = k3;

  // compute second unknown
  if (fabs(v1[0]) > tolAbs_)
    h = (c[0] - a[0] + v2[0] * k) / v1[0];
  else if (fabs(v1[1]) > tolAbs_)
    h = (c[1] - a[1] + v2[1] * k) / v1[1];
  else
    h = (c[2] - a[2] + v2[2] * k) / v1[2];

  // compute point of intersection
  e = c + v2 * k; // do not use h, because it was computed from k

  if (h > -tolRel_) { // we consider only [a,inf)
    if ((k > -tolRel_) && (k < 1.0 + tolRel_)) { // intersection on [c,d]?
      if (h < 1.0 + tolRel_) { // intersection on [a,b]?
        // treat special cases
        if (fabs(h - 1.0) < tolRel_) // h=1 means intersection in b
          return INTERSECT_IN_B;
        if (fabs(k - 1.0) < tolRel_) // k=1 means intersection in d
          return INTERSECT_IN_D;
        if (fabs(k) < tolRel_) { // k=0 means intersection in c
          if (fabs(h) < tolRel_) // h=0 means intersection in a
            return INTERSECT_A_EQ_C;
          return INTERSECT_IN_C;
        }
        if (fabs(h) < tolRel_) // h=0 means intersection in a
          return INTERSECT_IN_A;
        return INTERSECT_CROSS; // X intersection
      }
      return INTERSECT_ON_LINE2; // [a,inf) with [c,d]
    }
    return INTERSECT_OUTSIDE; // intersection not on any line
  }

  return INTERSECT_NONE; // no intersection (with [a,inf))
}

bool MortarInterface::CutPolys(StdVector< Vector<Double> > &p1,
                    StdVector< Vector<Double> > &p2, const bool coplanar,
                    StdVector< Vector<Double> > &r) const
{
  Double r1, r2;
  UInt i, inside = 0, nCuts = 0, start_cur = p1.GetSize();
  Vector<Double> c1, c2, e;
  Vector<Double> temp1, temp2;
  struct Intersection {
    UInt index;
    UInt type;
    bool swap;
    Vector<Double> loc;
  } cuts[2];

#ifdef CHECK_INDEX
  // check that we have actually polygons
  if ((p1.GetSize() < 3) || (p2.GetSize() < 3))
  {
    EXCEPTION("A polygon must consist of 3 points at least");
    return false;
  }
#endif

  // compute surrounding circles of both polygons
  r1 = PolyCentroid(p1, c1);
  r2 = PolyCentroid(p2, c2);

  // quit, if surrounding circles do not intersect
  temp1 = (c1 - c2);
  if (r1 + r2 < temp1.NormL2())
    return false;

  // if interface is not coplanar then project p1 onto p2
  if (!coplanar) {
    Double scale;
    Vector<Double> n;

    // compute surface normal of p2
    temp1 = p2[1]- p2[0];
    temp2 = p2[2] - p2[0];
    temp1.CrossProduct(temp2, n);
    n.Normalize();

    // project each point of p1
    for (i = 0; i < p1.GetSize(); ++i) {
      temp1 = p1[i] - p2[0];
      scale = n.Inner( temp1 );
      p1[i] -= n * scale;
    }
  }

  // Count those points of p1 that are contained in p2. Choose a point
  // that lies outside of p2 as starting point.
  for (i = 0; i < p1.GetSize(); ++i) {
    if (PointInsidePoly(p1[i], p2, &c2))
      ++inside;
    else if ((inside == 0) || (start_cur == p1.GetSize()))
      start_cur = i;

  }
  // Is p1 contained completely in p2?
  if (inside == p1.GetSize()) {
    r = p1; // intersection is p1 itself (for convex polygons)
    return true;
  }
  if (inside == 0) { // no points of p1 inside p2?
    for (i = 0; i < p2.GetSize(); ++i) {
      if (PointInsidePoly(p2[i], p1, &c1))
        ++inside;
    }
    // Is p2 contained completely in p1?
    if (inside == p2.GetSize()) {
      r = p2; // intersection is p2 itself (for convex polygons)
      return true;
    }
  }
  // WARNING: One can not conclude that two polygons do not intersect from
  // the fact that no point lies inside the other polygon.

  // make sure that both polygons have the same orientation
  PolygonIterator pi1(p1, start_cur), pi2(p2);
  temp1 = p1[1] - p1[0];
  temp2 = p1[2] - p1[0];
  temp1.CrossProduct(temp2, c1);

  temp1 = p2[1] - p2[0];
  temp2 = p2[2] - p2[0];
  temp1.CrossProduct(temp2, c2);
  if (c1 * c2 < 0.0)
    pi2.Reverse();

  // find the first cut of two edges of the polygons
  do {
    do {
      LineIntersectType cuttype
        = CutLines(*pi1, pi1.Next(), *pi2, pi2.Next(), e);
      Intersection cut = {pi2.GetPos(), cuttype, false, e};
      // See what kind of cut we have found.
      // This section is different from the main loop, because we do
      // not know if we cut from outside into p2 or vice versa. This
      // can happen despite the starting point lying outside, because
      // an edge of p1 might cut p2 into halves.
      switch (cuttype) {
        case INTERSECT_CROSS: // lines cross each other
          break; // always store the cut
        case INTERSECT_IN_A:
          // see if [a,b] lies inside of p2 or not
          if ((CutLines(*pi1, pi1.Next(), *pi2, pi2.Next(2), e)
                >= INTERSECT_ON_LINE2) ||
              (CutLines(*pi1, pi1.Next(), pi2.Prev(),
                        pi2.Next(), e) >= INTERSECT_ON_LINE2))
            break;
          continue; // [a,b] lies outside of p2 => no cut
        case INTERSECT_IN_C:
        case INTERSECT_A_EQ_C:
          // does [a,b] cut into p2?
          if (CutLines(*pi1, pi1.Next(), pi2.Prev(), pi2.Next(),
                e) >= INTERSECT_ON_LINE2)
            break; // yes, store cut
          // does [c,d] lie inside of p1?
          if (CutLines(*pi2, pi2.Next(), pi1.Prev(), pi1.Next(),
                e) >= INTERSECT_ON_LINE2) {
            cut.swap = true; // continue with p2
            break; // add cut
          }
          if (cuttype != INTERSECT_A_EQ_C) { // not for a=c
            if (CutLines(*pi2, pi2.Next(), *pi1, pi1.Next(2), e)
                >= INTERSECT_ON_LINE2) {
              cut.swap = true; // continue with p2
              break; // add cut
            }
          }
          continue; // polygons touch in c only
        case INTERSECT_A_AND_C:
          nCuts = 2;
          cuts[0] = cut;
          cuts[0].type = INTERSECT_IN_A;
          cuts[1].index = cuts[0].index;
          cuts[1].type = INTERSECT_IN_C;
          cuts[1].swap = true;
          cuts[1].loc = *pi2;
          continue;
        default:
          // cases for cuts in b and d are not stored, because
          // they would give duplicate cuts (polygons are closed!)
          continue;
      }

      if (nCuts < 2)
        cuts[nCuts] = cut; // store the cut
      ++nCuts;

      // test next line of passive polygon
    } while ( ! (++pi2).AtBegin() );

    // exit loop if first active line with cut is found
    if (nCuts > 0)
      break;

    // next line of active polygon
  } while ( ! (++pi1).AtBegin() );

  // do not proceed if there is no cut to start with
  if (nCuts == 0)
    return false;
  // make sure there are not more cuts than possible
  if (nCuts > 2) {
    EXCEPTION("A line cannot cut more than two edges of a convex polygon");
  }

  // save the position of the first cut in the active polygon
  pi1.SetBegin(pi1.GetPos());

  if (nCuts == 2) { // two cuts
    // make sure we do not treat a duplicate cut
    temp1 = (cuts[1].loc - cuts[0].loc);
    if (temp1.NormL2() < tolAbs_) {
      nCuts = 1;
    } else {
      // Here we can assume that we have found two "real" cuts. In
      // this case [a,b] runs completely through p2.
      // => sort cuts by distance to a
      temp1 = (cuts[0].loc - *pi1);
      temp2 = (cuts[1].loc - *pi1);
      if (temp1.NormL2() < temp2.NormL2())
      {
        r.Push_back(cuts[0].loc);
        r.Push_back(cuts[1].loc);

        pi2.Seek(cuts[0].index);
        if ((cuts[0].type != INTERSECT_IN_C) &&
            (cuts[0].type != INTERSECT_A_EQ_C)) ++pi2;
        pi2.SetBegin();

        pi2.Seek(cuts[1].index);
      } else {
        r.Push_back(cuts[1].loc);
        r.Push_back(cuts[0].loc);

        pi2.Seek(cuts[1].index);
        if ((cuts[1].type != INTERSECT_IN_C) &&
            (cuts[1].type != INTERSECT_A_EQ_C)) ++pi2;
        pi2.SetBegin();

        pi2.Seek(cuts[0].index);
      }
      ++pi1; // avoid finding the same cut twice
      pi1.Swap(pi2); // continue with p2
    }
  } // NO else clause here in order to catch duplicate cuts
  if (nCuts == 1) { // one cut with line [a,b]
    // save the position of the first cut with the passive polygon

    pi2.Seek(cuts[0].index);
    if ((cuts[0].type != INTERSECT_IN_C) &&
        (cuts[0].type != INTERSECT_A_EQ_C)) ++pi2;
    pi2.SetBegin();

    // store first point of intersection polygon
    r.Push_back(cuts[0].loc);

    // avoid finding the same cut twice
    ++pi1;
    // continue with p2, if indicated
    if (cuts[0].swap) {
      pi1.Swap(pi2);
    } else {// [a,b] cuts into p2, so add b
      r.Push_back(*pi1);
    }
  }

  bool swap, swapped = false;
  UInt start_act = pi1.GetPos(), start_pas = pi2.GetPos();

  // main loop
  do {
    swap = false;
    if ( ! pi2.AtBegin() || !swapped ) {
      do {
        switch (CutLines(*pi1, pi1.Next(), *pi2, pi2.Next(), e)) {
          case INTERSECT_CROSS:
          case INTERSECT_IN_C:
            r.Push_back(e);
            swap = true;
            break;
          case INTERSECT_IN_A:
          case INTERSECT_A_AND_C:
            if (CutLines(*pi1, pi1.Next(), *pi2, pi2.Next(2), e)
                >= INTERSECT_ON_LINE2)
              continue;
          case INTERSECT_A_EQ_C:
            if (CutLines(*pi1, pi1.Next(), pi2.Prev(), pi2.Next(),
                  e) >= INTERSECT_ON_LINE2)
              continue;
            swap = true;
            break;
          default:
            break;
        }
        if (swap) break;
      } while ( ! (++pi2).AtBegin() );
    }
    ++pi1;
    if (swap) {
      pi1.Swap(pi2);
      start_act = pi1.GetPos();
      pi1.Seek(start_act);
      start_pas = pi2.GetPos();
      swapped = true;
    } else {
      r.Push_back(*pi1);
      // Return to the point directly after the last cut (we can do
      // this due to the polygons being convex and having the same
      // orientation).
      pi2.Seek(start_pas);
    }
  } while ( !pi1.AtEnd() );

  r.Erase(r.GetSize() - 1);
  return (r.GetSize() > 2);
}

bool MortarInterface::PointInsidePoly(const Vector<Double> &p,
                           const StdVector< Vector<Double> > &poly,
                           const Vector<Double> *const c) const
{
  bool result = false;
  LineIntersectType s;
  Vector<Double> center, e, temp;
  ConstPolygonIterator pi(poly);

  // compute centroid of polygon, if not given
  if (c == NULL)
    PolyCentroid(poly, center);
  else
    center = *c;

  // Test if p is the centroid of the polygon (should always lie inside of a
  // convex polygon). In this case the algorithm below will not work.
  temp = (p - center);
  if ( temp.NormL2() < tolAbs_)
    return true;

  // try intersecting [c,p] with each edge of the polygon
  do {
    s = CutLines(center, p, *pi, pi.Next(), e);
    if (s <= INTERSECT_OUTSIDE)
      continue;
    if ((s == INTERSECT_ON_LINE2) || (s == INTERSECT_IN_B)) {
      result = true;
      break;
    }
    if ((s == INTERSECT_CROSS) || (s >= INTERSECT_IN_C)) {
      result = false;
      break;
    }
  } while ( ! (++pi).AtBegin() );

  return result;
}

Double MortarInterface::PolyCentroid(const StdVector< Vector<Double> > &p,
                          Vector<Double> &c) const
{
  UInt i, n = p.GetSize();
  Double r, r_max = 0.0;
  Vector<Double> temp;

  // set c to 0
  c.Resize(3);
  c.Init();

  // compute center of gravity
  for (i = 0; i < n; ++i)
    c += p[i];
  c /= (Double) n;

  // find point with maximum distance from centroid
  for (i = 0; i < n; ++i) {
    temp = (p[i] - c);
    r = temp.NormL2();
    if (r > r_max)
      r_max = r;
  }

  return r_max;
}

UInt MortarInterface::TriangulatePoly(const StdVector< Vector<Double> > &p,
                           StdVector<MortarNcSurfElem*> &tri,
                           StdVector<UInt> &newNodes )
{
  UInt nodeNo, firstNo = tri.GetSize();
  MortarNcSurfElem *ncElem;
  UInt numPoints = p.GetSize();
  Vector<Double> temp1, temp2;

  switch (numPoints) {
    case 3:
      ncElem = new MortarNcSurfElem;
      ncElem->type = Elem::ET_TRIA3;
      ncElem->connect.Resize(3);

      ptGrid_->AddNode(p[0], nodeNo);
      newNodes.Push_back(nodeNo);
      ncElem->connect[0] = nodeNo;
      ptGrid_->AddNode(p[1], nodeNo);
      newNodes.Push_back(nodeNo);
      ncElem->connect[1] = nodeNo;
      ptGrid_->AddNode(p[2], nodeNo);
      newNodes.Push_back(nodeNo);
      ncElem->connect[2] = nodeNo;

      tri.Push_back(ncElem);
      break;
    case 4:
      ncElem = new MortarNcSurfElem;
      ncElem->type = Elem::ET_QUAD4;
      ncElem->connect.Resize(4);

      ptGrid_->AddNode(p[0], nodeNo);
      newNodes.Push_back(nodeNo);
      ncElem->connect[0] = nodeNo;
      ptGrid_->AddNode(p[1], nodeNo);
      newNodes.Push_back(nodeNo);
      ncElem->connect[1] = nodeNo;
      ptGrid_->AddNode(p[2], nodeNo);
      newNodes.Push_back(nodeNo);
      ncElem->connect[2] = nodeNo;
      ptGrid_->AddNode(p[3], nodeNo);
      newNodes.Push_back(nodeNo);
      ncElem->connect[3] = nodeNo;

      tri.Push_back(ncElem);

      break;
    default:
      UInt i, centerNode, firstNode;
      Vector<Double> c;

      PolyCentroid(p, c);
      ptGrid_->AddNode(c, centerNode);
      newNodes.Push_back(centerNode);
      ptGrid_->AddNode(p[0], nodeNo);
      newNodes.Push_back(nodeNo);
      firstNode = nodeNo;

      for (i = 1; i < p.GetSize(); ++i)
      {
        ncElem = new MortarNcSurfElem;
        ncElem->type = Elem::ET_TRIA3;
        ncElem->connect.Resize(3);

        ncElem->connect[0] = nodeNo;
        ptGrid_->AddNode(p[i], nodeNo);
        newNodes.Push_back(nodeNo);
        ncElem->connect[1] = nodeNo;
        ncElem->connect[2] = centerNode;

        tri.Push_back(ncElem);
      }

      ncElem = new MortarNcSurfElem;
      ncElem->type = Elem::ET_TRIA3;
      ncElem->connect.Resize(3);

      ncElem->connect[0] = nodeNo;
      ncElem->connect[1] = firstNode;
      ncElem->connect[2] = centerNode;

      tri.Push_back(ncElem);
      break;
  }
  return firstNo;
}

} /* namespace CoupledField */
