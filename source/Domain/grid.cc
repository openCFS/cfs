// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <cmath>
#include <string>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/programOptions.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/WriteInfo.hh"
#include "Elements/elements_header.hh"
#include "elem.hh"
#include "domain.hh"
#include "General/exception.hh"
#include "Utils/mathfunctions.hh"
#include "Utils/coordSystem.hh"
#include "polygonIterator.hh"

#include "grid.hh"

#include <limits>

namespace CoupledField
{

  // declare class specific logging stream
  DECLARE_LOG(grid)
  DEFINE_LOG(grid, "grid")

  Grid::Grid()
  {
    isInitialized_ = false; // set by FinishInit()

    if (!ptQ1)     ptQ1     = new Quad1FE();
    if (!ptQ2)     ptQ2     = new Quad2FE();
    if (!ptQ9)     ptQ9     = new Quad9FE();
    if (!ptTet1)   ptTet1   = new Tetra1FE();
    if (!ptTet2)   ptTet2   = new Tetra2FE();
    if (!ptL1)     ptL1     = new Line1FE();
    if (!ptL2)     ptL2     = new Line2FE();
    if (!ptTr1)    ptTr1    = new Triangle1FE();
    if (!ptTr2)    ptTr2    = new Triangle2FE();
    if (!ptHexa1)  ptHexa1  = new Hexa1FE();
    if (!ptHexa2)  ptHexa2  = new Hexa2FE();
    if (!ptHexa27) ptHexa27 = new Hexa27FE();
    if (!ptPyra1)  ptPyra1  = new Pyra1FE();
    if (!ptPyra2)  ptPyra2  = new Pyra2FE();
    if (!ptWedge1) ptWedge1 = new Wedge1FE();
    if (!ptWedge2) ptWedge2 = new Wedge2FE();

 #ifdef USE_SCRIPTING
    // Register functions
    RegisterFunctions();
 #endif

    region_.SetName("Grid::region");
    region_.Add(ALL_REGIONS, "all");
  }



  Grid::~Grid()
  {
    delete ptQ1;      ptQ1 = NULL;
    delete ptQ2;      ptQ2 = NULL;
    delete ptQ9;      ptQ9 = NULL;
    delete ptTet1;    ptTet1 = NULL;
    delete ptTet2;    ptTet2 = NULL;
    delete ptL1;      ptL1 = NULL;
    delete ptL2;      ptL2 = NULL;
    delete ptTr1;     ptTr1 = NULL;
    delete ptTr2;     ptTr2 = NULL;
    delete ptHexa1;   ptHexa1 = NULL;
    delete ptHexa2;   ptHexa2 = NULL;
    delete ptHexa27;  ptHexa27 = NULL;
    delete ptPyra1;   ptPyra1 = NULL;
    delete ptPyra2;   ptPyra2 = NULL;
    delete ptWedge1;  ptWedge1 = NULL;
    delete ptWedge2;  ptWedge2 = NULL;
  }

  Grid::RegionData::RegionData()
  {
    name = "";
    id   = -1;
    type = NOT_SET;
    type_idx = -1;
    regular = false;
    homogeneous = true;
    barycenters = false;
  }

  RegionIdType Grid::AddRegion(const std::string& name)
  {
    RegionData rd;
    rd.name = name;
    rd.id   = regionData.GetSize();
    regionData.Push_back(rd);
    region_.Add(rd.id, rd.name);

    nameTypeMap_[name] = EntityList::REGION;

    return rd.id;
  }

  void Grid::AddRegions(const StdVector<std::string> & names, StdVector<RegionIdType> & ids)
  {
    ids.Resize(names.GetSize());

    for(unsigned int i = 0; i < names.GetSize(); i++)
      ids[i] = AddRegion(names[i]);
  }

  RegionIdType Grid::AddRegion(const std::string& name, RegionType type)
  {
    // Check if entities with given name exist already
    if( nameTypeMap_.find( name) != nameTypeMap_.end() )
      EXCEPTION("Entities with name " << name << " are already defined");

    if(!isInitialized_)
      EXCEPTION("Cannot add a region to an uninitialized grid!");

    RegionIdType id = AddRegion(name);
    regionData[id].type = type;

    StdVector<Elem*> dummy_elems;
    std::set<UInt> dummy_nodes;

    if(type == SURFACE_REGION)
    {
      regionData[id].type_idx = surfRegionIds_.GetSize();
      surfRegionIds_.Push_back(id);
      surfElems_.Push_back(dummy_elems);
      surfElemNodes_.Push_back(dummy_nodes);
    }
    else
    {
      regionData[id].type_idx = volRegionIds_.GetSize();
      volRegionIds_.Push_back(id);
      volElems_.Push_back(dummy_elems);
      volElemNodes_.Push_back(dummy_nodes);
    }

    return id;
  }


  UInt Grid::GetNumVolRegions()
  {
    return volRegionIds_.GetSize();
  }

  UInt Grid::GetNumSurfRegions()
  {
    return surfRegionIds_.GetSize();
  }

  bool Grid::IsRegionRegular(StdVector<RegionIdType>& regions) const
  {
    bool regular = true;

    for(unsigned int i = 0; i < regions.GetSize(); i++)
      if(!regionData[regions[i]].regular)
        regular = false;

    return regular;
  }

  void Grid::GetRegionNames( StdVector<std::string>& regionNames )
  {
    regionNames.Resize(regionData.GetSize());

    for(UInt i = 0; i < regionData.GetSize(); i++)
      regionNames[i] = regionData[i].name;
  }


  void Grid::GetVolRegionIds( StdVector<RegionIdType> & volRegions ) {
    volRegions = volRegionIds_;
  }

  void Grid::GetSurfRegionIds( StdVector<RegionIdType> & surfRegions ) {

    surfRegions = surfRegionIds_;
  }

  UInt Grid::SetElementBarycenters(RegionIdType reg, bool updated)
  {
    RegionData& rd = regionData[reg];

    if(rd.barycenters) return 0;

    // common for all elements
    Matrix<Double>  coords;

    // our operation target
    StdVector<Elem*>& elems = rd.type == VOLUME_REGION ? volElems_[rd.type_idx] : surfElems_[rd.type_idx];
    for(UInt i = 0;  i < elems.GetSize(); i++)
    {
      Elem* elem = elems[i];

      StdVector<UInt>& connect = elem->connect;

      GetElemNodesCoord(coords, connect, updated);

      // a barycenter is simply the average of all coordinates
      // TODO: handle axis symmetry!!
      BaseFE::CalcBarycenter(coords, elem->barycenter);
    }

    rd.barycenters = true; // don't do it again!

    return elems.GetSize();
  }

  shared_ptr<EntityList> Grid::GetEntityList( EntityList::ListType listType,
                                              const std::string& name,
                                              EntityList::DefineType defineType ) {

    // First check, if there any entites with this name at all
    if( nameTypeMap_.find( name) == nameTypeMap_.end() ) {
      EXCEPTION( "There are no entities with name '" << name
                 << "' in the mesh" ) ;
    }

    EntityList::DefineType entityType = nameTypeMap_[name];

    shared_ptr<EntityList> ret;

    if( listType == EntityList::ELEM_LIST ) {
      shared_ptr<ElemList> eList  = shared_ptr<ElemList>( new ElemList(this) );
      if( entityType == EntityList::REGION ) {
        RegionIdType regionId = GetRegion().Parse( name );
        eList->SetRegion( regionId);
      } else {
        eList->SetNamedElems( name );
      }
      ret = eList;

    } else if( listType == EntityList::SURF_ELEM_LIST ) {
      shared_ptr<SurfElemList> surfList  =
        shared_ptr<SurfElemList>( new SurfElemList(this) );
      if( entityType == EntityList::REGION ) {
        RegionIdType regionId = GetRegion().Parse( name );
        surfList->SetRegion( regionId);
      } else {
        surfList->SetNamedElems( name );
      }
      ret = surfList;

    } else if( listType == EntityList::NODE_LIST ) {
      shared_ptr<NodeList> nodeList = shared_ptr<NodeList>( new NodeList(this) );
      // Check if name describes a nodeList
      if( entityType == EntityList::NAMED_NODES ) {
        StdVector<std::string> nodeNames;
        GetListNodeNames( nodeNames );
        nodeList->SetNamedNodes( name );
      } else if( entityType == EntityList::NAMED_ELEMS ) {
          nodeList->SetNamedNodes( name );
      } else if( entityType == EntityList::REGION ) {
        RegionIdType regionId = GetRegion().Parse( name );
        nodeList->SetNodesOfRegion( regionId );
      } else {
        EXCEPTION("GetEntityList with NODE_LIST works only with regions"
                  << " and named nodes!" );
      }
      ret = nodeList;
    } else if( listType == EntityList::REGION_LIST ) {
      shared_ptr<RegionList> regionList =
        shared_ptr<RegionList>( new RegionList(this) );
      if( entityType == EntityList::REGION ) {
        RegionIdType regionId = GetRegion().Parse( name );
        regionList->SetRegionId( regionId );
      } else {
        EXCEPTION( "GetEntityList with REGION_LIST works only with regions!" );
      }
      ret = regionList;
    } else {
      EXCEPTION( "Type '" << listType << "' describes no EntityList which is created "
                 << "by the grid-class." );
    }

    return ret;

  }

   void Grid::Dump()
   {
     StdVector<Elem*>   elems;

     std::cout << "Grid: elements=" << GetNumElems() << " nodes=" << GetNumNodes() << std::endl;

     for(UInt i = 0; i < regionData.GetSize(); i++)
     {
       GetElems(elems, i);

       std::cout << "region: " << regionData[i].name << " id=" << i << " elements=" << elems.GetSize() <<  std::endl;
     }
   }


  bool Grid::InitNonmatchingInterfaces() {
    UInt numIfaces;
    StdVector<SurfElem*> ifaceElems;
    StdVector<SurfElem*> masterElems;
    StdVector<SurfElem*> slaveElems;
    std::string ncRegionName;
    RegionIdType slaveId;
    RegionIdType masterId;
    ParamNode* pNode;
    StdVector<ParamNode*> ncIfaceNodes;

    pNode = param->Get("domain");
    pNode = pNode->Get("ncInterfaceList",false);
    if(!pNode)
    {
      Info->PrintF("","\nNo non-conforming interfaces found!\n");
      LOG_DBG2(grid) << "NonMatching: No non-conforming interfaces found!";
      return true;
    }

    ncIfaceNodes = pNode->GetList("ncInterface");
    numIfaces = ncIfaceNodes.GetSize();

    ncInterfaces_.Resize(numIfaces);

    for (UInt i = 0; i < numIfaces; ++i) {
      ParamNode* ncIfaceNode = ncIfaceNodes[i];
      ParamNode* rotNode = ncIfaceNode->Get("rotation", false);

      ncInterfaces_[i].name = ncIfaceNode->Get("name")->AsString();

      masterId = GetRegion().Parse(ncIfaceNode->Get("masterSide")->AsString());
      slaveId = GetRegion().Parse(ncIfaceNode->Get("slaveSide")->AsString());
      if ((masterId == -1) || (slaveId == -1)) {
        EXCEPTION("Cannot find master/slave regions of ncInterface '"
            << ncInterfaces_[i].name << "'.");
      }
      ncInterfaces_[i].masterRegion = masterId;
      ncInterfaces_[i].slaveRegion = slaveId;

      if (GetDim() == 2) {
        ncInterfaces_[i].intersectAlgo = LINE_INTERSECT;
      }
      else {
        std::string isecCalc;
        ncIfaceNode->Get("isecCalculation", isecCalc, false);
        if (isecCalc == "coaxi")
          ncInterfaces_[i].intersectAlgo = RECT_INTERSECT;
        else
          ncInterfaces_[i].intersectAlgo = POLYGON_INTERSECT;
      }

      ncIfaceNode->Get("tolAbs", ncInterfaces_[i].tolAbs, false);
      ncIfaceNode->Get("tolRel", ncInterfaces_[i].tolRel, false);

      if (rotNode) {
        rotNode->Get("coordSysId", ncInterfaces_[i].coordSysId, false);
        rotNode->Get("angleStep", ncInterfaces_[i].rotationAngle, false);
      }
      else {
        ncInterfaces_[i].rotationAngle = 0.0;
      }

      if (ncInterfaces_[i].rotationAngle > 0.0) {
        ncInterfaces_[i].coplanar = false;
      }
      else {
        UInt numElems;
        GetSurfElems(masterElems, masterId);
        GetSurfElems(slaveElems, slaveId);

        if ((masterElems.GetSize() == 0) || (slaveElems.GetSize() == 0)) {
          EXCEPTION("Cannot find surface elements in master/slave regions of "
              << "ncInterface '" << ncInterfaces_[i].name << "'.");
        }

        numElems = masterElems.GetSize();
        ifaceElems.Reserve(numElems + slaveElems.GetSize());
        for (UInt j = 0; j < numElems; ++j) {
          ifaceElems.Push_back(masterElems[j]);
        }
        masterElems.Clear();

        numElems = slaveElems.GetSize();
        for (UInt j = 0; j < numElems; ++j) {
          ifaceElems.Push_back(slaveElems[j]);
        }
        slaveElems.Clear();

        ncInterfaces_[i].coplanar = IsSurfacePlanar(ifaceElems);
        ifaceElems.Clear();
      }

      ncInterfaces_[i].region = AddSurfaceRegion(ncInterfaces_[i].name);

      //if (ncInterfaces_[i].rotationAngle == 0.0)
      UpdateNcIntersection(ncInterfaces_[i]);
    }

    return true;
  }


  void Grid::UpdateNcIntersection(const ncInterface& ncIf) {
    StdVector<SurfElem*> masterElems;
    StdVector<SurfElem*> slaveElems;
    StdVector<NCElem*> ncElems;
    StdVector<SurfElem*> ncElemsHelper;
    StdVector<UInt> masterNodes;
    StdVector<UInt> ncElemIds;
    std::string isecCalc;

    GetSurfElems(masterElems, ncIf.masterRegion);
    GetSurfElems(slaveElems, ncIf.slaveRegion);

    switch (ncIf.intersectAlgo) {
      case LINE_INTERSECT:
        LOG_DBG3(grid) << "\ncomputing integration grid of ncInterfaces...";
        for (UInt i = 0; i < masterElems.GetSize(); ++i) {
          for (UInt j = 0; j < slaveElems.GetSize(); ++j) {
            SideOnSide(masterElems[i], slaveElems[j], ncIf.coplanar,
                       ncIf.rotationAngle != 0.0, ncElems);
          }
        }
        break;
      case RECT_INTERSECT:
        if (!ncIf.coplanar) {
          EXCEPTION("Only coplanar interfaces are supported with coaxial "
                    << "rectangle algorithm.");
        }
        for (UInt i = 0; i < masterElems.GetSize(); ++i) {
          for(UInt j = 0; j < slaveElems.GetSize(); ++j) {
            BaseFE* m_el = masterElems[i]->ptElem;
            BaseFE* s_el = slaveElems[j]->ptElem;
            if ( (m_el != ptQ1 && m_el != ptQ2 && m_el != ptQ9)
               ||(s_el != ptQ1 && s_el != ptQ2 && s_el != ptQ9))
            {
              EXCEPTION("Only quadrilaterals can be intersected with coaxial "
                        << "rectangle algorithm.");
            }

            if(RectangleOnRectangle(masterElems[i], slaveElems[j],
                ncIf.coplanar, ncElems))
            {
              LOG_DBG3(grid) << "Intersection between "
                             << masterElems[i]->elemNum << " and "
                             << slaveElems[j]->elemNum << std::endl;
            }
          }
        }
        break;
      case POLYGON_INTERSECT:
        for (UInt i = 0; i < masterElems.GetSize(); ++i) {
          for (UInt j = 0; j < slaveElems.GetSize(); ++j) {
            BaseFE* m_el = masterElems[i]->ptElem;
            BaseFE* s_el = slaveElems[j]->ptElem;
            if ( (m_el != ptTr1 && m_el != ptQ1 && m_el != ptQ2 && m_el != ptQ9 && m_el != ptTr2)
               ||(s_el != ptTr1 && s_el != ptQ1 && s_el != ptQ2 && s_el != ptQ9 && s_el != ptTr2))
            {
              EXCEPTION("Only triangles and quadrilaterals can be intersected"
                        << " with polygon algorithm.");
            }

            if (PolygonOnPolygon(masterElems[i], slaveElems[j], ncIf.coplanar,
                                 ncIf.rotationAngle != 0.0, ncElems,
                                 ncIf.tolAbs, ncIf.tolRel))
            {
              LOG_DBG3(grid) << "Intersection between "
                << masterElems[i]->elemNum << " and "
                << slaveElems[j]->elemNum << std::endl;
            }
          }
        }
    }

    if(ncElems.GetSize() > 0)
    {
      ClearRegion(ncIf.region);

      ncElemsHelper.Resize(ncElems.GetSize());

      for(UInt i=0; i<ncElems.GetSize(); i++)
      {
        ncElemsHelper[i] = ncElems[i];
      }

      AddSurfaceElems(ncIf.region, ncElemsHelper, ncElemIds);
    }
    else {
      EXCEPTION("No intersection elements were computed for non-conforming"
                << " interface '" << ncIf.name
                << "'. Please check your mesh file.");
    }
  }


  bool Grid::IsSurfacePlanar(const StdVector<SurfElem*>& ifaceElems)
  {
    std::set<Integer> ifaceNodes;
    std::set<Integer>::iterator it,end;
    Vector<Double> pv1, pv2, pv3;
    Vector<Double> v1, v2, normal, n;
    Double innerProd = 0.0, norm1 = 0.0, norm2 = 0.0;
    Double eps = 1e-15;
    UInt pnum=0; // number of point in ifaceNodes
    UInt dim = GetDim();

    // Determine set of points on interface.
    for(UInt i=0; i<ifaceElems.GetSize(); i++)
    {
      for(UInt j=0; j<ifaceElems[i]->connect.GetSize(); j++)
      {
        ifaceNodes.insert(ifaceElems[i]->connect[j]);
      }
    }

    normal.Resize(3);
    normal[0] = 0.0;
    normal[1] = 0.0;
    normal[2] = 0.0;

    // Loop through all interface points and determine if they are coplanar.
    for(it=ifaceNodes.begin(), end=ifaceNodes.end(); it!=end; it++)
    {
      if(pnum == 0)
        // Get coordinates of first point.
        GetNodeCoordinate(pv1, (*it));
      else if(pnum == 1)
      {
        // Get coordinates of second point, compute vector
        // from first to second point and normalize it.
        GetNodeCoordinate(pv2, (*it));
        v1 = pv2 - pv1;
        // Normalize v1.
        norm1 = Normalize(v1);
        if(norm1 < eps)
          norm1 = 0.0;
      }
      else
      {
        // Get coordinate of point pnum, compute vector
        // from first to point pnum and normalize it.
        GetNodeCoordinate(pv3, (*it));
        v2 = pv3 - pv1;

        // Normalize v2.
        norm2 = Normalize(v2);
        if(norm2 < eps)
          norm2 = 0.0;

        // If point pnum and the first point coincide
        // both points are on the interface -> continue.
        if(norm2 == 0.0)
          continue;

        // If the first and second point coincide we need
        // a new vector v1 to get the direction of the interface.
        if(norm1 == 0.0)
        {
          v1 = v2;
          norm1 = norm2;
          // We don't want to compare the vector to itself!
          continue;
        }

        switch(dim)
        {
        case 2:
          // Compute the inner product of v1 and v2. If the vectors
          // point to the same direction the result must be 1
          v1.Inner(v2, innerProd);

          if((1.0 - std::fabs(innerProd)) >= eps)
            return false;

          break;

        case 3:
          if(normal.NormL2() == 0)
          {
            CrossProd(v1, v2, normal);

            // We may not use a zero-normal vector to perform our
            // test for coplanarity.
            if(normal.NormL2() < eps)
            {
              normal[0] = 0.0;
              normal[1] = 0.0;
              normal[2] = 0.0;
              continue;
            }
            Normalize(normal);
            continue;
          }
          CrossProd(v1, v2, n);

          // If the norm of n is smaller than eps we have multiplied
          // linearly dependant vectors -> a point in the plane
          // has been found.
          if(n.NormL2() < eps)
          {
            continue;
          }

          Normalize(n);

          n.Inner(normal, innerProd);

          // At this place we should have either linearly dependant normal
          // and n vectors or they face to different directions.
          // The value of innerProd indicates what is the case.
          if((1.0 - std::fabs(innerProd)) >= eps)
            return false;

          break;
        }

      }
      pnum++;
    }

    return true;
  }

  bool Grid::IsNcInterfaceCoplanar(const std::string &ncIfaceName) {
    return IsNcInterfaceCoplanar(GetRegion().Parse(ncIfaceName));
  }

  bool Grid::IsNcInterfaceCoplanar(RegionIdType regionId) {
    UInt numIfaces = ncInterfaces_.GetSize();

    for (UInt i = 0; i < numIfaces; ++i) {
      if (ncInterfaces_[i].region == regionId)
        return ncInterfaces_[i].coplanar;
    }

    EXCEPTION("Non-matching grid interface does not exist: " << region_.ToString(regionId));
  }

  /****************************************************************************
   **
   ** SideOnSide
   **
   **   computes the local coordinates of the overlap of the master and slave
   **   element with respect to the master side in the order of the
   **   orientation of the slave side element. It pushes back the intersection
   **   element to elemList.
   **
   ** Input Parameters:
   **   ifaceElemM:  Master Side
   **   ifaceElemS:  Slave Side
   **   collinear:   indicates if the interface is coplanar or not
   **   coordUpdate: indicates if updated coordinates should be used
   **
   ** Output Parameters:
   **   elemList: the found intersection NCElems will be pushed
   **                     back to this vector
   **
   */


  bool Grid::SideOnSide (SurfElem *ifaceElem1, SurfElem *ifaceElem2,
                         bool collinear, bool coordUpdate,
                         StdVector<NCElem*>& elemList)
  {
    // c0, c1, d0 and d1 are the endpoints of the two line elements
    //
    //           d0 x--------------------------x d1
    // c0 x---------|-----------x c1

    Vector<Double> c0, c1, d0, d1, tmp; // endpoint-coordinates of the two lines
    Vector<Double> diff0, diff1;
    Vector<Double> s, t;
    Vector<Double> normal;
    StdVector<UInt> connect2;
    Double dist, dist1, dist2, fac;
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
    GetNodeCoordinate(c0, nodenum_c0);
    GetNodeCoordinate(c1, nodenum_c1);
    GetNodeCoordinate(d0, nodenum_d0);
    GetNodeCoordinate(d1, nodenum_d1);

    // Project master nodes onto slave element, if interface is not coplanar
    if (!collinear) {
      UInt node_num;
      Vector<Double> new_node;
      // compute maximal allowed distance as sum of lengths of both lines
      tmp = c1 - c0;
      Double maxDist = tmp.NormL2();
      tmp = d1 - d0;
      maxDist += tmp.NormL2();
      // compute normal vector of slave element
      CalcSurfNormal(normal, *ifaceElem2);
      // compute distance of c0 to plane of slave element
      tmp = c0 - d0;
      normal.Inner(tmp, fac);
      // make sure that distance does not exceed maximum distance
      if (fabs(fac) > maxDist)
        return false;
      // do the projection
      c0 -= normal * fac;
      // add new node
      new_node.Resize(3);
      new_node[0] = c0[0];
      new_node[1] = c0[1];
      if (c0.GetSize() == 2) {
        new_node[2] = 0.0;
      } else {
        new_node[2] = c0[2];
      }
      AddNode(new_node, node_num);
      nodenum_c0 = node_num;

      // do the same for c1
      tmp = c1 - d0;
      normal.Inner(tmp, fac);
      c1 -= normal * fac;
      new_node[0] = c1[0];
      new_node[1] = c1[1];
      if (c1.GetSize() == 2) {
        new_node[2] = 0.0;
      } else {
        new_node[2] = c1[2];
      }
      AddNode(new_node, node_num);
      nodenum_c1 = node_num;
    }

    // Compute and normalize vector from c0 to c1.
    // This becomes the new x-unit vector.
    diff1 = c1 - c0;
    dist = diff1.NormL2();
    fac = 1.0 / dist;
    diff1 *= fac;

    // Compute x1 coordinate of line2 in respect
    // to line1.
    diff0 = d0 - c0;
    dist1 = diff0.NormL2();
    diff0.Inner(diff1, s[0]);
    s[0] *= fac;

    // Compute x2 coordinate of line2 in respect
    // to line1.
    diff0 = d1 - c0;
    dist2 = diff0.NormL2();
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

    // Check if an intersection between line1 and line2
    // exists.
    if(t[0] >= 1.0)
      return false;

    if(t[1] <= 0.0)
      return false;

    NCElem* ncElem = new NCElem;
    ncElem->connect.Resize(2);

    relativeElemVol = t[1] - t[0];

    // If an intersection exists, we must distinguish
    // 4 different cases.
    if(t[0] <= 0)
    {
      ncElem->connect[0] = nodenum_c0;

      if(t[1] >= 1)
      {
        // connect2[0] x--------|--------------------|-----x connect2[1]
        //                   c0 x--------------------x c1

        relativeElemVol = 1;
        ncElem->connect[1] = nodenum_c1;
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

    if(relativeElemVol < 1e-3) {
      std::stringstream sstr;
      sstr << "Rejecting ncElem due to a relative volume of " << relativeElemVol;
      sstr << std::endl;
      sstr << "  for intersection of elements " << ifaceElem1->elemNum;
      sstr << " (" << region_.ToString(ifaceElem1->regionId) << ") ";
      sstr << "and " << ifaceElem2->elemNum;
      sstr << " (" << region_.ToString(ifaceElem2->regionId) << ") ";
      Warning(sstr.str().c_str(), __FILE__, __LINE__);
      delete ncElem;
      return false;
    }

    ncElem->ptElem = ptL1;
    ncElem->ptLagrangeParent = ifaceElem2;
    ncElem->ptSurfParent = ifaceElem1;

    elemList.Push_back(ncElem);

    return true;
  }

  bool Grid::RectangleOnRectangle(SurfElem *ifaceElem1, SurfElem *ifaceElem2,
                                  bool coplanar, StdVector<NCElem*>& elemList)
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
    GetNodeCoordinate(c0, ifaceElem1->connect[0]);
    GetNodeCoordinate(c1, ifaceElem1->connect[1]);
    GetNodeCoordinate(c2, ifaceElem1->connect[3]);
    GetNodeCoordinate(d0, ifaceElem2->connect[0]);
    GetNodeCoordinate(d1, ifaceElem2->connect[1]);
    GetNodeCoordinate(d2, ifaceElem2->connect[2]);

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

    NCElem* ncElem = new NCElem;
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
            AddNode(tmp, nodeNr);
            ncElem->connect[2] = nodeNr;
            tmp = c0 + diffX*0.0 + diffY*t[3];
            AddNode(tmp, nodeNr);
            ncElem->connect[3] = nodeNr;
          }

        }
        else
        {
          if(t[3] >= 1)
          {
            tmp = c0 + diffX*0.0 + diffY*t[1];
            AddNode(tmp, nodeNr);
            ncElem->connect[0] = nodeNr;
            tmp = c0 + diffX     + diffY*t[1];
            AddNode(tmp, nodeNr);
            ncElem->connect[1] = nodeNr;
            ncElem->connect[2] = ifaceElem1->connect[2];
            ncElem->connect[3] = ifaceElem1->connect[3];
          }
          else
          {
            tmp = c0 + diffX*0.0 + diffY*t[1];
            AddNode(tmp, nodeNr);
            ncElem->connect[0] = nodeNr;
            tmp = c0 + diffX     + diffY*t[1];
            AddNode(tmp, nodeNr);
            ncElem->connect[1] = nodeNr;
            tmp = c0 + diffX     + diffY*t[3];
            AddNode(tmp, nodeNr);
            ncElem->connect[2] = nodeNr;
            tmp = c0 + diffX*0.0 + diffY*t[3];
            AddNode(tmp, nodeNr);
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
            AddNode(tmp, nodeNr);
            ncElem->connect[1] = nodeNr;
            tmp = c0 + diffX*t[2] + diffY;
            AddNode(tmp, nodeNr);
            ncElem->connect[2] = nodeNr;
            ncElem->connect[3] = ifaceElem1->connect[3];
          }
          else
          {
            ncElem->connect[0] = ifaceElem1->connect[0];
            tmp = c0 + diffX*t[2] + diffY*0.0;
            AddNode(tmp, nodeNr);
            ncElem->connect[1] = nodeNr;
            ncElem->connect[2] = connect2[2];
            tmp = c0 + diffX*0.0  + diffY*t[3];
            AddNode(tmp, nodeNr);
            ncElem->connect[3] = nodeNr;
          }

        }
        else
        {
          if(t[3] >= 1)
          {
            tmp = c0 + diffX*0.0  + diffY*t[1];
            AddNode(tmp, nodeNr);
            ncElem->connect[0] = nodeNr;
            ncElem->connect[1] = connect2[1];
            tmp = c0 + diffX*t[2] + diffY;
            AddNode(tmp, nodeNr);
            ncElem->connect[2] = nodeNr;
            ncElem->connect[3] = ifaceElem1->connect[3];
          }
          else
          {
            tmp = c0 + diffX*0.0  + diffY*t[1];
            AddNode(tmp, nodeNr);
            ncElem->connect[0] = nodeNr;
            ncElem->connect[1] = connect2[1];
            ncElem->connect[2] = connect2[2];
            tmp = c0 + diffX*0.0  + diffY*t[3];
            AddNode(tmp, nodeNr);
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
            AddNode(tmp, nodeNr);
            ncElem->connect[0] = nodeNr;
            ncElem->connect[1] = ifaceElem1->connect[1];
            ncElem->connect[2] = ifaceElem1->connect[2];
            tmp = c0 + diffX*t[0] + diffY;
            AddNode(tmp, nodeNr);
            ncElem->connect[3] = nodeNr;
          }
          else
          {
            tmp = c0 + diffX*t[0] + diffY*0.0;
            AddNode(tmp, nodeNr);
            ncElem->connect[0] = nodeNr;
            ncElem->connect[1] = ifaceElem1->connect[1];
            tmp = c0 + diffX      + diffY*t[3];
            AddNode(tmp, nodeNr);
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
            AddNode(tmp, nodeNr);
            ncElem->connect[1] = nodeNr;
            ncElem->connect[2] = ifaceElem1->connect[2];
            tmp = c0 + diffX*t[0] + diffY;
            AddNode(tmp, nodeNr);
            ncElem->connect[3] = nodeNr;
          }
          else
          {
            ncElem->connect[0] = connect2[0];
            tmp = c0 + diffX      + diffY*t[1];
            AddNode(tmp, nodeNr);
            ncElem->connect[1] = nodeNr;
            tmp = c0 + diffX      + diffY*t[3];
            AddNode(tmp, nodeNr);
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
            AddNode(tmp, nodeNr);
            ncElem->connect[0] = nodeNr;
            tmp = c0 + diffX*t[2] + diffY*0.0;
            AddNode(tmp, nodeNr);
            ncElem->connect[1] = nodeNr;
            tmp = c0 + diffX*t[2] + diffY;
            AddNode(tmp, nodeNr);
            ncElem->connect[2] = nodeNr;
            tmp = c0 + diffX*t[0] + diffY;
            AddNode(tmp, nodeNr);
            ncElem->connect[3] = nodeNr;
          }
          else
          {
            tmp = c0 + diffX*t[0] + diffY*0.0;
            AddNode(tmp, nodeNr);
            ncElem->connect[0] = nodeNr;
            tmp = c0 + diffX*t[2] + diffY*0.0;
            AddNode(tmp, nodeNr);
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
            AddNode(tmp, nodeNr);
            ncElem->connect[2] = nodeNr;
            tmp = c0 + diffX*t[0] + diffY;
            AddNode(tmp, nodeNr);
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

    ncElem->ptElem = ptQ1;
    ncElem->ptLagrangeParent = ifaceElem2;
    ncElem->ptSurfParent = ifaceElem1;

    elemList.Push_back(ncElem);

    return true;
  }

  bool Grid::PolygonOnPolygon(SurfElem *ifElem1, SurfElem *ifElem2,
                              bool coplanar, bool coordUpdate,
                              StdVector<NCElem*> &elemList,
                              Double absTol, Double relTol)
  {
    UInt i;
    StdVector< Vector<Double> > p1, p2, r;
	UInt p1Size, p2Size;

    polysectAbsTol_ = absTol;
    polysectRelTol_ = relTol;

    switch(ifElem1->ptElem->feType()) {
      case Elem::TRIA3:
      case Elem::TRIA6:
        p1Size = 3;
        break;

      case Elem::QUAD4:
      case Elem::QUAD8:
      case Elem::QUAD9:
        p1Size = 4;
        break;

      default:
    	EXCEPTION("First argument to PolygonOnPolygon may not be of type '"
    			<< Elem::feType.ToString(ifElem1->ptElem->feType()) << "!");
    }

    p1.Resize(p1Size);
    for (i = 0; i < p1Size; ++i)
      GetNodeCoordinate(p1[i], ifElem1->connect[i], coordUpdate);

    switch(ifElem2->ptElem->feType()) {
      case Elem::TRIA3:
      case Elem::TRIA6:
        p2Size = 3;
        break;

      case Elem::QUAD4:
      case Elem::QUAD8:
      case Elem::QUAD9:
        p2Size = 4;
        break;

      default:
    	EXCEPTION("Second argument to PolygonOnPolygon may not be of type '"
    			<< Elem::feType.ToString(ifElem2->ptElem->feType()) << "'!");
    }

    p2.Resize(p2Size);
    for (i = 0; i < p2Size; ++i)
      GetNodeCoordinate(p2[i], ifElem2->connect[i], coordUpdate);

    if (CutPolys(p1, p2, coplanar, r))
    {
      UInt start = TriangulatePoly(r, elemList);

      for (UInt i = start; i < elemList.GetSize(); ++i)
      {
        elemList[i]->ptLagrangeParent = ifElem2;
        elemList[i]->ptSurfParent = ifElem1;
      }

      return true;
    }
    return false;
  }

  Grid::LineIntersectType Grid::CutLines(const Vector<Double> &a,
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

    if (CoLinear(v1, v2)) { // lines are parallel
      // if line from a to d is also parallel then lines may intersect
      e = d - a;
      if (CoLinear(v1, e)) {
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
        if (fabs(l_ac + l_ad - l2) < polysectAbsTol_) {
          e = a;
          if (l_ac < polysectAbsTol_) // is a=c?
            return INTERSECT_A_EQ_C;
          if (l_ad < polysectAbsTol_) { // is a=d?
            // usually a cut at d wins over a cut at a,
            // but c might be an endpoint here, too
            if (fabs(l_ac + l_bc - l1) < polysectAbsTol_) {
              e = c;
              return INTERSECT_IN_C;
            }
            return INTERSECT_IN_D;
          }
          // does c lie on [a,b]?
          if (fabs(l_ac + l_bc - l1) < polysectAbsTol_)
            return INTERSECT_A_AND_C; // intersection is [a,c]
          
          return INTERSECT_IN_A;
        }
        if (fabs(l_bc + l_bd - l2) < polysectAbsTol_) {
          e = b;
          return INTERSECT_IN_B;
        }
        if (fabs(l_ac + l_bc - l1) < polysectAbsTol_) {
          e = c;
          return INTERSECT_IN_C;
        }
        if (fabs(l_ad + l_bd - l1) < polysectAbsTol_) {
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
    if (fabs(denom1) > polysectAbsTol_)
      k1 = (v1[0] * (c[1] - a[1]) + v1[1] * (a[0] - c[0])) / denom1;
    denom2 = v1[2] * v2[0] - v1[0] * v2[2];
    if (fabs(denom2) > polysectAbsTol_)
      k2 = (v1[0] * (c[2] - a[2]) + v1[2] * (a[0] - c[0])) / denom2;
    denom3 = v1[2] * v2[1] - v1[1] * v2[2];
    if (fabs(denom3) > polysectAbsTol_)
      k3 = (v1[1] * (c[2] - a[2]) + v1[2] * (a[1] - c[1])) /denom3;

    // If this system has no solution, lines do not intersect.
    if ((fabs(denom1) <= polysectAbsTol_)
        && (fabs(denom2) <= polysectAbsTol_)
        && (fabs(denom3) <= polysectAbsTol_))
      return INTERSECT_NONE;

    /* TODO: jens
     * This check makes no sense for 3 k's. Maybe add a check based on the
     * standard deviation of k.
     */
    /*if ((fabs(denom1) > polysectAbsTol_)
        && (fabs(denom2) > polysectAbsTol_)) {
      if (fabs(k1 - k2) > polysectRelTol_)
        return INTERSECT_NONE;
    }*/

    // select the right one out of (k1,k2,k3)
    if (fabs(denom1) > polysectAbsTol_)
      k = k1;
    else if (fabs(denom2) > polysectAbsTol_)
      k = k2;
    else
      k = k3;

    // compute second unknown
    if (fabs(v1[0]) > polysectAbsTol_)
      h = (c[0] - a[0] + v2[0] * k) / v1[0];
    else if (fabs(v1[1]) > polysectAbsTol_)
      h = (c[1] - a[1] + v2[1] * k) / v1[1];
    else
      h = (c[2] - a[2] + v2[2] * k) / v1[2];

    // compute point of intersection
    e = c + v2 * k; // do not use h, because it was computed from k

    if (h > -polysectRelTol_) { // we consider only [a,inf)
      if ((k > -polysectRelTol_) && (k < 1.0 + polysectRelTol_)) { // intersection on [c,d]?
        if (h < 1.0 + polysectRelTol_) { // intersection on [a,b]?
          // treat special cases
          if (fabs(h - 1.0) < polysectRelTol_) // h=1 means intersection in b
            return INTERSECT_IN_B;
          if (fabs(k - 1.0) < polysectRelTol_) // k=1 means intersection in d
            return INTERSECT_IN_D;
          if (fabs(k) < polysectRelTol_) { // k=0 means intersection in c
            if (fabs(h) < polysectRelTol_) // h=0 means intersection in a
              return INTERSECT_A_EQ_C;
            return INTERSECT_IN_C;
          }
          if (fabs(h) < polysectRelTol_) // h=0 means intersection in a
            return INTERSECT_IN_A;
          return INTERSECT_CROSS; // X intersection
        }
        return INTERSECT_ON_LINE2; // [a,inf) with [c,d]
      }
      return INTERSECT_OUTSIDE; // intersection not on any line
    }

    return INTERSECT_NONE; // no intersection (with [a,inf))
  }

  bool Grid::CutPolys(StdVector< Vector<Double> > &p1,
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
      CrossProd(temp1,temp2, n);
      Normalize(n);

      // project each point of p1
      for (i = 0; i < p1.GetSize(); ++i) {
        temp1 = p1[i] - p2[0];
        n.Inner( temp1, scale);
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
    CrossProd(temp1, temp2, c1);

    temp1 = p2[1] - p2[0];
    temp2 = p2[2] - p2[0];
    CrossProd(temp1, temp2, c2);
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
      if (temp1.NormL2() < polysectAbsTol_) {
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

  bool Grid::PointInsidePoly(const Vector<Double> &p,
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
    if ( temp.NormL2() < polysectAbsTol_)
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

  Double Grid::PolyCentroid(const StdVector< Vector<Double> > &p,
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

  UInt Grid::TriangulatePoly(const StdVector< Vector<Double> > &p,
                             StdVector<NCElem*> &tri)
  {
    UInt nodeNo, firstNo = tri.GetSize();
    NCElem *ncElem/*, *ncElem2*/;
    UInt numPoints = p.GetSize();
    Vector<Double> temp1, temp2;

    switch (numPoints) {
      case 3:
        ncElem = new NCElem;
        ncElem->ptElem = ptTr1;
        ncElem->connect.Resize(3);

        AddNode(p[0], nodeNo);
        ncElem->connect[0] = nodeNo;
        AddNode(p[1], nodeNo);
        ncElem->connect[1] = nodeNo;
        AddNode(p[2], nodeNo);
        ncElem->connect[2] = nodeNo;

        tri.Push_back(ncElem);
        break;
      case 4:
        ncElem = new NCElem;
        ncElem->ptElem = ptQ1;
        ncElem->connect.Resize(4);

        AddNode(p[0], nodeNo);
        ncElem->connect[0] = nodeNo;
        AddNode(p[1], nodeNo);
        ncElem->connect[1] = nodeNo;
        AddNode(p[2], nodeNo);
        ncElem->connect[2] = nodeNo;
        AddNode(p[3], nodeNo);
        ncElem->connect[3] = nodeNo;

        tri.Push_back(ncElem);

        /*ncElem = new NCElem;
        ncElem2 = new NCElem;
        ncElem->connect.Resize(3);
        ncElem2->connect.Resize(3);

        if ((p[0] - p[2]).NormL2() < (p[1] - p[3]).NormL2())
        {
          AddNode(p[0], nodeNo);
          ncElem->connect[0] = nodeNo;
          ncElem2->connect[0] = nodeNo;
          AddNode(p[1], nodeNo);
          ncElem->connect[1] = nodeNo;
          AddNode(p[2], nodeNo);
          ncElem->connect[2] = nodeNo;
          ncElem2->connect[1] = nodeNo;
          AddNode(p[3], nodeNo);
          ncElem2->connect[2] = nodeNo;
        } else {
          AddNode(p[0], nodeNo);
          ncElem->connect[0] = nodeNo;
          AddNode(p[1], nodeNo);
          ncElem->connect[1] = nodeNo;
          ncElem2->connect[0] = nodeNo;
          AddNode(p[2], nodeNo);
          ncElem2->connect[1] = nodeNo;
          AddNode(p[3], nodeNo);
          ncElem->connect[2] = nodeNo;
          ncElem2->connect[2] = nodeNo;
        }

        tri.Push_back(ncElem);
        tri.Push_back(ncElem2);*/
        break;
      default:
        UInt i, centerNode, firstNode;
        Vector<Double> c;

        PolyCentroid(p, c);
        AddNode(c, centerNode);
        AddNode(p[0], nodeNo);
        firstNode = nodeNo;

        for (i = 1; i < p.GetSize(); ++i)
        {
          ncElem = new NCElem;
          ncElem->ptElem = ptTr1;
          ncElem->connect.Resize(3);

          ncElem->connect[0] = nodeNo;
          AddNode(p[i], nodeNo);
          ncElem->connect[1] = nodeNo;
          ncElem->connect[2] = centerNode;

          tri.Push_back(ncElem);
        }

        ncElem = new NCElem;
        ncElem->ptElem = ptTr1;
        ncElem->connect.Resize(3);

        ncElem->connect[0] = nodeNo;
        ncElem->connect[1] = firstNode;
        ncElem->connect[2] = centerNode;

        tri.Push_back(ncElem);
    }
    return firstNo;
  }

  void Grid::SurfRegionFromVolRegions(
    const std::string& surfRegionName,
    const std::string& region1,
    const std::string& region2) {

    StdVector<SurfElem*> newSurfaceElems;
    StdVector<Elem*> region1Elems;
    StdVector<UInt> region2Nodes;
    StdVector<UInt> surfElemIds;
    RegionIdType region1Id;
    RegionIdType region2Id;
    RegionIdType surfRegionId;
    Elem* el;
    SurfElem* surfEl;

    if(GetDim() != 2)
      EXCEPTION("SurfRegionFromVolRegions is only implemented for 2D!");

    if(region2 == "NO_REGION") {
      SurfRegionFromSingleVolRegion(surfRegionName, region1);
      return;
    }

    region1Id = GetRegion().Parse(region1);
    region2Id = GetRegion().Parse(region2);

    this->GetElems(region1Elems, region1Id);
    this->GetNodesByRegion(region2Nodes, region2Id);

    UInt nElemsRegion1 = region1Elems.GetSize();
    UInt numCorners;
    UInt lastCornerInRegion2;

    for(UInt i=0; i<nElemsRegion1; i++) {
      el = region1Elems[i];
      numCorners = el->ptElem->GetNumCorners();
      lastCornerInRegion2 = 0;
      for(UInt n=0; n<numCorners; ) {
        if(region2Nodes.Find(el->connect[n]) < 0) {
          n++;
          continue;
        }

        if(region2Nodes.Find(el->connect[(n+1) % numCorners]) < 0) {
          n+=2;
          continue;
        }

        surfEl = new SurfElem();
        surfEl->connect.Resize(3);
        surfEl->connect[0] = el->connect[n];
        surfEl->connect[1] = el->connect[(n+1) % numCorners];
        surfEl->ptElem = ptL1;

        switch(el->ptElem->feType()) {
        case Elem::TRIA6:
        case Elem::QUAD8:
        case Elem::QUAD9:
          surfEl->connect[2] = el->connect[n+numCorners];
          surfEl->ptElem = ptL2;
          break;
        default:
          break;
        }

        newSurfaceElems.Push_back(surfEl);

        n++;
      }
    }

    if(newSurfaceElems.GetSize() != 0)
    {
      surfRegionId = AddSurfaceRegion(surfRegionName);
      AddSurfaceElems( surfRegionId, newSurfaceElems, surfElemIds);
    }

  }

  void Grid::SurfRegionFromSingleVolRegion(
    const std::string& surfRegionName,
    const std::string& region)
  {
    StdVector<SurfElem*> newSurfaceElems;
    StdVector<Elem*> regionElems;
    StdVector<UInt> regionNodes;
    StdVector<UInt> surfElemIds;
    std::map<UInt, UInt> edgeCounts;
    RegionIdType regionId;
    RegionIdType surfRegionId;
    Elem* el;
    SurfElem* surfEl;

    regionId = GetRegion().Parse(region);

    this->GetElems(regionElems, regionId);
    this->GetNodesByRegion(regionNodes, regionId);

    UInt nElemsRegion = regionElems.GetSize();
    UInt numEdges;
    //    UInt lastCornerOnBnd;

    for(UInt i=0; i<nElemsRegion; i++) {
      el = regionElems[i];
      // get number of edges
      numEdges = el->ptElem->GetNumEdges();
      for(UInt n=0; n<numEdges; n++) {
        UInt edgeNum = el->edges[n] < 0 ? -el->edges[n] : el->edges[n];
        edgeCounts[edgeNum]++;
      }
    }

    for(UInt i=0; i<nElemsRegion; i++) {
      el = regionElems[i];
      // get number of edges
      numEdges = el->ptElem->GetNumEdges();
      for(UInt n=0; n<numEdges; n++) {
        UInt edgeNum = el->edges[n] < 0 ? -el->edges[n] : el->edges[n];
        if(edgeCounts[edgeNum] != 1)
        {
          n++;
          continue;
        }


        /*
                  UInt cornerCount1 = cornerCounts[el->connect[n]];
                  UInt cornerCount2 = cornerCounts[el->connect[(n+1) % numCorners]];

                  switch(cornerCount1) {
                  case 1:
                  case 2:
                  if(cornerCount2 > 3) {
                  n+=2;
                  continue;
                  }
                  break;
                  case 3:
                  if(cornerCount2 > 2) {
                  n+=2;
                  continue;
                  }
                  break;
                  default:
                  n++;
                  continue;
                  }
                  /*
                  if(cornerCount1 > 2) {
                  n++;
                  continue;
                  }

                  if(cornerCount1 > 2) {
                  n+=2;
                  continue;
                  }
        */

        surfEl = new SurfElem();
        surfEl->connect.Resize(3);
        surfEl->connect[0] = el->connect[n];
        surfEl->connect[1] = el->connect[(n+1) % numEdges];
        surfEl->ptElem = ptL1;

        switch(el->ptElem->feType()) {
        case Elem::TRIA6:
        case Elem::QUAD8:
        case Elem::QUAD9:
          surfEl->connect[2] = el->connect[n+numEdges];
          surfEl->ptElem = ptL2;
          break;
        default:
          break;
        }

        newSurfaceElems.Push_back(surfEl);

        n++;
      }
    }

    if(newSurfaceElems.GetSize() != 0)
    {
      surfRegionId = AddSurfaceRegion(surfRegionName);
      AddSurfaceElems( surfRegionId, newSurfaceElems, surfElemIds);
    }

  }

  // =======================================================================
  // Method wrappers for scripting
  // =======================================================================

  void Grid::Wrap_GetNodesByName() {
    SCRIPT_GET(std::string, name);
    StdVector<UInt> nodeNrs;
    GetNodesByName( nodeNrs, name);
    nodeNrs.ToString(SCRIPT_RETVAL);
  }

  void Grid::Wrap_GetNodeCoordinate() {
    SCRIPT_GET(UInt, nodeNr);
    Vector<Double> coord;
    GetNodeCoordinate(coord, nodeNr );
    SCRIPT_RETVAL = coord.ToString( ' ');;
  }

  void Grid::Wrap_GetNodesByRegion() {
    SCRIPT_GET(std::string, name);
    StdVector<UInt> nodeNrs;
    GetNodesByRegion(nodeNrs, GetRegion().Parse(name) );
    nodeNrs.ToString(SCRIPT_RETVAL);
  }

  void Grid::Wrap_GetListNodeNames() {
    GetListNodeNames( SCRIPT_RETVAL );
  }

  void Grid::Wrap_GetListElemNames() {
    GetListElemNames( SCRIPT_RETVAL );
  }

  void Grid::Wrap_GetRegionNames() {
    StdVector<RegionIdType> regionIds;
    for(UInt i = 0; i < regionData.GetSize(); i++)
      regionIds[i] = i;
    region_.ToString( regionIds, SCRIPT_RETVAL );

  }

  void Grid::Wrap_GetNumNodes() {
    SCRIPT_RETVAL.Push_back( lexical_cast<std::string>( GetNumNodes() ) );
  }

  void Grid::Wrap_GetNumElems() {
    SCRIPT_RETVAL.Push_back( lexical_cast<std::string>( GetNumElems() ) );

  }

  void Grid::Wrap_GetNumSurfElems() {
    SCRIPT_RETVAL.Push_back( lexical_cast<std::string>( GetNumSurfElems() ) );

  }

  void Grid::Wrap_GetNumVolElems() {
    SCRIPT_RETVAL.Push_back( lexical_cast<std::string>( GetNumVolElems() ) );
  }

  void Grid::Wrap_GetNumNodesOfRegion() {
    SCRIPT_GET(std::string, name);
    StdVector<RegionIdType> ids;
    ids.Push_back( GetRegion().Parse (name ) );
    SCRIPT_RETVAL.Push_back( lexical_cast<std::string>( GetNumNodes ( ids ) ) );
  }

  void Grid::Wrap_GetNumElemsOfRegion() {
    SCRIPT_GET(std::string, name);
    StdVector<RegionIdType> ids;
    ids.Push_back( GetRegion().Parse (name ) );
    SCRIPT_RETVAL.Push_back( lexical_cast<std::string>( GetNumElems (ids ) ) );
  }


  void Grid::RegisterFunctions() {

    typedef FctPointer<Grid> FCPT;
    StdVector<ArgList> a;
    StdVector<FCPT*> pt;
    StdVector<std::string> name;

    // --- GetNodesByName ---
    a.Push_back();
    a.Last().RegisterParam("name", ArgList::STRING);
    pt.Push_back( new FCPT(this,&Grid::Wrap_GetNodesByName) );
    name.Push_back( "getNodesByName" );

    // --- GetNodeCoordinate ---
    a.Push_back();
    a.Last().RegisterParam("nodeNr", ArgList::UINT);
    pt.Push_back(new FCPT(this, &Grid::Wrap_GetNodeCoordinate) );
    name.Push_back( "getNodeCoordinate" );

    // --- GetNodesByRegion ---
    a.Push_back();
    a.Last().RegisterParam("name", ArgList::STRING);
    pt.Push_back(new FCPT(this, &Grid::Wrap_GetNodesByRegion) );
    name.Push_back( "getNodesByRegion" );

    // --- GetListNodeNames ---
    a.Push_back();
    pt.Push_back(new FCPT(this, &Grid::Wrap_GetListNodeNames) );
    name.Push_back( "getNodeNames" );

    // --- GetListElemNames ---
    a.Push_back();
    pt.Push_back(new FCPT(this, &Grid::Wrap_GetListElemNames) );
    name.Push_back( "getElemNames" );

    // --- GetRegionNames---
    a.Push_back();
    pt.Push_back(new FCPT(this, &Grid::Wrap_GetRegionNames) );
    name.Push_back( "getRegionNames" );

    // --- GetNumNodes ---
    a.Push_back();
    pt.Push_back(new FCPT(this, &Grid::Wrap_GetNumNodes) );
    name.Push_back( "getNumNodes" );

    // --- GetNumElems ---
    a.Push_back();
    pt.Push_back(new FCPT(this, &Grid::Wrap_GetNumElems) );
    name.Push_back( "getNumElems" );

    // --- GetNumSurfElems ---
    a.Push_back();
    pt.Push_back(new FCPT(this, &Grid::Wrap_GetNumSurfElems) );
    name.Push_back( "getNumSurfElems" );

    // --- GetNumVolElems ---
    a.Push_back();
    pt.Push_back(new FCPT(this, &Grid::Wrap_GetNumVolElems) );
    name.Push_back( "getNumVolElems" );

    // --- GetNumNodesOfRegion ---
    a.Push_back();
    a.Last().RegisterParam("name", ArgList::STRING);
    pt.Push_back(new FCPT(this, &Grid::Wrap_GetNumNodesOfRegion) );
    name.Push_back( "getNumNodesOfRegion" );

    // --- GetNumElemsOfRegion ---
    a.Push_back();
    a.Last().RegisterParam("name", ArgList::STRING);
    pt.Push_back( new FCPT(this, &Grid::Wrap_GetNumElemsOfRegion) );
    name.Push_back( "getNumElemsOfRegion" );

    // Now register all functions with scripting
    for (UInt i = 0; i < pt.GetSize(); i++ ) {
      Script_RegisterFct(name[i], pt[i], a[i] );
    }

  }

#ifdef USE_INTERPOLATION

  /*
  // coordinates for 9 boxes of a grid
  int p[9*4]   = { 0,0,1,1,  1,0,2,1,  2,0,3,1, // lower
                   0,1,1,2,  1,1,2,2,  2,1,3,2, // middle
                   0,2,1,3,  1,2,2,3,  2,2,3,3};// upper
  // 9 boxes
  Grid::Box boxes[9] = { Box( p,    p+ 2),  Box( p+ 4, p+ 6),  Box( p+ 8, p+10),
                   Box( p+12, p+14),  Box( p+16, p+18),  Box( p+20, p+22),
                   Box( p+24, p+26),  Box( p+28, p+30),  Box( p+32, p+34)};
  // 2 selected boxes as query; center and upper right
  Grid::Box query[2] = { Box( p+16, p+18),  Box( p+32, p+34)};

  // callback function object writing results to an output iterator
  template <class OutputIterator>
  struct Report {
    OutputIterator it;
    Report( OutputIterator i) : it(i) {} // store iterator in object
    // We write the id-number of box a to the output iterator assuming
    // that box b (the query box) is not interesting in the result.
    void operator()( const Grid::Box& a, const Grid::Box&) { *it++ = a.id(); }
  };
  template <class Iter> // helper function to create the function object
  Report<Iter> report( Iter it) { return Report<Iter>(it); }




  // callback function that reports all truly intersecting triangles
  void report_inters( const HandleBox& a, const HandleBox& b) {
    std::cout << "Box " << *a.handle() << " and "
              << *b.handle() << " intersect";
    std::cout << '.' << std::endl;
  }

  // callback function that reports all truly intersecting triangles
  void report_inters2( const Box& a, const Box& b) {
    std::cout << "Box " << a.id() << " and "
              << b.id() << " intersect";
    std::cout << '.' << std::endl;
  }

  void Grid::intersection() {
    double xmin, ymin, xmax, ymax, zmin, zmax;
    StdVector<Elem*> elems;
    Point p;
    UInt dim = GetDim();
    std::vector<HandleBox> elemBoxes;
    std::vector<HandleBox> elemBoxes2;

    std::cout << std::endl;
    std::cout << std::endl;

    GetVolElems(elems, ALL_REGIONS);
    elemBoxes.resize(elems.GetSize());

    for(UInt i = 0, m=elems.GetSize(); i < m; i++)
    {
      GetNodeCoordinate(p, elems[i]->connect[0]);

      xmin = xmax = p[0];
      ymin = ymax = p[1];
      zmin = zmax = p[2];

      for(UInt j = 1, n=elems[i]->connect.GetSize(); j < n; j++)
      {
        GetNodeCoordinate(p, elems[i]->connect[j]);
        xmin = p[0] < xmin ? p[0] : xmin;
        xmax = p[0] > xmax ? p[0] : xmax;
        ymin = p[1] < ymin ? p[1] : ymin;
        ymax = p[1] > ymax ? p[1] : ymax;
        zmin = p[2] < zmin ? p[2] : zmin;
        zmax = p[2] > zmax ? p[2] : zmax;
      }

      elemBoxes[i] = HandleBox(BBox3D(xmin, ymin, zmin, xmax, ymax, zmax),
                               &elems[i]->elemNum);

      std::cout << "element " << elems[i]->elemNum << " BBox3D (" << xmin << ", " << ymin << ", " << zmin << ") (" << xmax <<  ", " << ymax << ", " << zmax << ")" << std::endl;

    }

    UInt test = 2333;
    elemBoxes2.push_back(HandleBox(BBox3D(0.001, 0.0, 0.001, 0.001, 0.0, 0.001),
                                   &test));

    // run the intersection algorithm and store results in a vector
    std::vector<std::size_t> result;
    //    CGAL::box_self_intersection_d( boxes, boxes+9, report_inters2 );
    //    CGAL::box_self_intersection_d( boxes, boxes+9, report( std::back_inserter( result)) );

    //    CGAL::box_intersection_d( boxes, boxes+9,
    //                              query, query+2,
    //                              report( std::back_inserter( result)));

    //    CGAL::box_intersection_d( elemBoxes.begin(), elemBoxes.end(),
    //                              elemBoxes2.begin(), elemBoxes2.end(),
    //                              report_inters);

    CGAL::box_intersection_d( elemBoxes.begin(), elemBoxes.end(),
                              elemBoxes2.begin(), elemBoxes2.end(),
                              report2( std::back_inserter( result), *this));

    // sort, check, and show result
    std::sort( result.begin(), result.end());
    //    std::size_t check1[13] = {0,1,2,3,4,4,5,5,6,7,7,8,8};
    //    assert(result.size() == 13 && std::equal(check1,check1+13,result.begin()));

    std::copy( result.begin(), result.end(),
               std::ostream_iterator<std::size_t>( std::cout, " "));
    std::cout << std::endl;
    std::cout << std::endl;
  }
  */

  const Elem* Grid::GetElemAtGlobalCoord(const Vector<double>& globCoord,
                                         Vector<double>& localCoords)
  {
    double xmin, ymin, xmax, ymax, zmin, zmax;
    //UInt dim = GetDim();
    std::vector<HandleBox> elemBoxes2;

    // If we haven't initialized the grid bounding boxes yet, do so now!
    if(elemBoxes_.empty())
    {
      StdVector<Elem*> elems;
      Point p;

      GetVolElems(elems, ALL_REGIONS);
      UInt size = elems.GetSize();
      
      elemBoxes_.reserve( size );

      for(UInt i = 0, m=elems.GetSize(); i < m; i++)
      {
        GetNodeCoordinate(p, elems[i]->connect[0]);

        xmin = xmax = p[0];
        ymin = ymax = p[1];
        zmin = zmax = p[2];

        for(UInt j = 1, n=elems[i]->connect.GetSize(); j < n; j++)
        {
          GetNodeCoordinate(p, elems[i]->connect[j]);
          xmin = p[0] < xmin ? p[0] : xmin;
          xmax = p[0] > xmax ? p[0] : xmax;
          ymin = p[1] < ymin ? p[1] : ymin;
          ymax = p[1] > ymax ? p[1] : ymax;
          zmin = p[2] < zmin ? p[2] : zmin;
          zmax = p[2] > zmax ? p[2] : zmax;
        }

        elemBoxes_.push_back( HandleBox(BBox3D(xmin, ymin, zmin, xmax, ymax, zmax),
                                  &elems[i]->elemNum) );

        //        std::cout << "element " << elems[i]->elemNum << " BBox3D (" << xmin << ", " << ymin << ", " << zmin << ") (" << xmax <<  ", " << ymax << ", " << zmax << ")" << std::endl;

      }
    }

    UInt test = 0xFFFFFFFF;
    elemBoxes2.push_back(HandleBox(BBox3D(globCoord[0], globCoord[1], globCoord[2],
                                          globCoord[0], globCoord[1], globCoord[2]),
                                   &test));

    // run the intersection algorithm and store results in a vector
    std::vector< std::pair<const Elem*, Vector<double> > > result;

    CGAL::box_intersection_d( elemBoxes_.begin(), elemBoxes_.end(),
                              elemBoxes2.begin(), elemBoxes2.end(),
                              report2( std::back_inserter( result ), *this));

    if(!result.empty())
    {
      localCoords = result.begin()->second;
      return result.begin()->first;
    }
    else
      return NULL;

  }


  void Grid::ComputeConservativeInterpolationWeights(const ElemList& destElemList,
          const NodeList& sourceNodeList,
          const std::string& coordSysId,
          ciTolerance& globalEpsilon,
          ciTolerance& localEpsilon,
          Double z,
          Double zEpsilon,
          std::vector< std::map<UInt, Double> >& consInterpWeights,
          StdVector<UInt> &unmapped_nodes)
  {
    Double xmin, ymin, xmax, ymax, zmin, zmax;
    Double globEps, locEps;
    UInt i;
    UInt dim = GetDim();
    UInt numSourceNodes = sourceNodeList.GetSize(); // number of nodes in source region
    UInt numActualSourceNodes = 0; // actual number of nodes to be interpolated
    CoordSystem* coordSys = domain->GetCoordSystem(coordSysId);
    Grid* source = sourceNodeList.GetGrid();
    UInt srcDim = source->GetDim();
    Point p;
    StdVector<UInt> sourceNodeNumbers, sourceNodeIndices;
    Vector<Double> point;
    Vector<Double> globPoint;
    std::vector< Vector<Double> > nodeCoords;
    std::vector<HandleBox> elemBoxes2;

    // initialize memory of interpolation weights, if necessary
    if (consInterpWeights.empty())
      consInterpWeights.resize(numSourceNodes);

    // initialize memory for coordinates of source nodes
    nodeCoords.resize(numSourceNodes);
    for (UInt i=0; i<numSourceNodes; ++i) {
      nodeCoords[i].Clear();
    }

    // If we haven't initialized the grid bounding boxes yet, do so now!
    if (elemBoxes_.empty()) {
      
      elemBoxes_.reserve( destElemList.GetSize() );
      const Elem* elem = NULL;

      for(UInt i = 0, m=destElemList.GetSize(); i < m; ++i)
      {
        elem = destElemList.GetElem(i);
        GetNodeCoordinate(p, elem->connect[0]);

        xmin = xmax = p[0];
        ymin = ymax = p[1];
        zmin = zmax = p[2];

        for(UInt j = 1, n=elem->connect.GetSize(); j < n; ++j)
        {
          GetNodeCoordinate(p, elem->connect[j]);
          xmin = p[0] < xmin ? p[0] : xmin;
          xmax = p[0]> xmax ? p[0] : xmax;
          ymin = p[1] < ymin ? p[1] : ymin;
          ymax = p[1]> ymax ? p[1] : ymax;
          zmin = p[2] < zmin ? p[2] : zmin;
          zmax = p[2]> zmax ? p[2] : zmax;
        }

        elemBoxes_.push_back( HandleBox(BBox3D(xmin, ymin, zmin, xmax, ymax, zmax),
                              &elem->elemNum) );

        //std::cout << "element " << elems[i]->elemNum << " BBox3D (" << xmin
        //          << ", " << ymin << ", " << zmin << ") (" << xmax <<  ", "
        //          << ymax << ", " << zmax << ")" << std::endl;

      }
    }

    // check that tolerances make sense
    if (globalEpsilon.end < globalEpsilon.start)
      globalEpsilon.end = globalEpsilon.start;
    if (globalEpsilon.inc == 0.0) {
      if (globalEpsilon.start == 0.0)
        globalEpsilon.inc = 1.0e-6;
      else
        globalEpsilon.inc = globalEpsilon.start;
    }
    if (localEpsilon.end < localEpsilon.start)
      localEpsilon.end = localEpsilon.start;
    if (localEpsilon.inc == 0.0) {
      if (localEpsilon.start == 0.0)
        localEpsilon.inc = 1.0e-3;
      else
        localEpsilon.inc = localEpsilon.start;
    }

    // loop over tolerance ranges
    for (locEps = localEpsilon.start;
         locEps <= localEpsilon.end;
         locEps += localEpsilon.inc) {

      // global tolerance is inner loop, because
      // it doesn't cause numerical errors
      for (globEps = globalEpsilon.start;
           globEps <= globalEpsilon.end;
           globEps += localEpsilon.inc) {

        EntityIterator it = sourceNodeList.GetIterator();
        sourceNodeNumbers.Clear();
        sourceNodeIndices.Clear();
        elemBoxes2.clear();
        i=0;

        // create a list of nodes that still need to be interpolated
        while (!it.IsEnd())
        {
          if (!consInterpWeights.empty())
          {
            if (consInterpWeights[i].empty())
            {
              // If the source grid is 3D and the destination grid is 2D
              // we have to map the global source coordinates into the
              // local source coordinate sys and only use those nodes
              // with given z.
              if ( srcDim == 3 && dim == 2)
              {
                source->GetNodeCoordinate(point, it.GetNode(), true);
                coordSys->Global2LocalCoord(globPoint, point);

                if ( std::fabs(globPoint[2] - z) < zEpsilon )
                {
                  sourceNodeNumbers.Push_back(it.GetNode());
                  sourceNodeIndices.Push_back(it.GetPos());
                }
              }
              else
              {
                sourceNodeNumbers.Push_back(it.GetNode());
                sourceNodeIndices.Push_back(it.GetPos());
              }
            }
          }
          else
          {
            if ( srcDim == 3 && dim == 2)
            {
              source->GetNodeCoordinate(point, it.GetNode(), true);
              coordSys->Global2LocalCoord(globPoint, point);

              if ( std::fabs(globPoint[2] - z) < zEpsilon )
              {
                sourceNodeNumbers.Push_back(it.GetNode());
                sourceNodeIndices.Push_back(it.GetPos());
              }
            }
            else
            {
              sourceNodeNumbers.Push_back(it.GetNode());
              sourceNodeIndices.Push_back(it.GetPos());
            }

          }

          it++;
          ++i;
        }

        numActualSourceNodes = sourceNodeNumbers.GetSize();

        if (numActualSourceNodes == 0)
          return;

        // add global tolerance to bounding boxes
        for(UInt n=0; n<numActualSourceNodes; ++n)
        {
          source->GetNodeCoordinate(point, sourceNodeNumbers[n], true);
          coordSys->Global2LocalCoord(nodeCoords[sourceNodeIndices[n]], point);

          // subtract origin here!
          if(dim == 3)
            elemBoxes2.push_back(HandleBox(
                BBox3D(nodeCoords[sourceNodeIndices[n]][0]-globEps,
                       nodeCoords[sourceNodeIndices[n]][1]-globEps,
                       nodeCoords[sourceNodeIndices[n]][2]-globEps,
                       nodeCoords[sourceNodeIndices[n]][0]+globEps,
                       nodeCoords[sourceNodeIndices[n]][1]+globEps,
                       nodeCoords[sourceNodeIndices[n]][2]+globEps),
                  &sourceNodeIndices[n]));
          else
            elemBoxes2.push_back(HandleBox(
                BBox3D(nodeCoords[sourceNodeIndices[n]][0]-globEps,
                       nodeCoords[sourceNodeIndices[n]][1]-globEps,
                       0.0,
                       nodeCoords[sourceNodeIndices[n]][0]+globEps,
                       nodeCoords[sourceNodeIndices[n]][1]+globEps,
                       0.0),
                  &sourceNodeIndices[n]));
        }

        // run the intersection algorithm and store results in a vector

        CGAL::box_intersection_d( elemBoxes_.begin(), elemBoxes_.end(),
                                  elemBoxes2.begin(), elemBoxes2.end(),
                                  GenConsInterpReportFunctor(destElemList,
                                      sourceNodeList,
                                      nodeCoords,
                                      locEps,
                                      consInterpWeights));
      }
    }

    for (UInt i=0; i<numActualSourceNodes; ++i) {
      if (consInterpWeights[sourceNodeIndices[i]].empty()) {
        // just add the number of the unmapped node.
        // DO NOT CLEAR the vector before this loop,
        // because this function is called several times.
        unmapped_nodes.Push_back(sourceNodeNumbers[i]);
      }
    }

  }

  Grid::ConsInterpReportFunctor::ConsInterpReportFunctor(const ElemList& destElemList,
                                                         const NodeList& sourceNodeList,
                                                         const std::vector< Vector<Double> >& nodeCoords,
                                                         Double localEpsilon,
                                                         std::vector< std::map<UInt, Double> >& consInterpWeights)
      : //destElemList_(destElemList),
      //sourceNodeList_(sourceNodeList),
        nodeCoords_(nodeCoords),
        localEpsilon_(localEpsilon),
        consInterpWeights_(consInterpWeights),
        nodeCounter_(0),
        percentage_(0),
        oldPercentage_(9)
    {
      numSourceNodes_ = consInterpWeights_.size();
      connect_.resize(64);

      sourceGrid_ = sourceNodeList.GetGrid();
      destGrid_ = destElemList.GetGrid();

      StdVector<UInt> destNodeNumbers;
      NodeList destNodeList(destGrid_);
      destNodeList.SetNodesOfRegion(destElemList.GetRegion());
      EntityIterator it = destNodeList.GetIterator();
      while(!it.IsEnd())
      {
        destNodeNumToPosMap_[it.GetNode()] = it.GetPos();
        it++;
      }

      destNodeNumbers = destNodeList.GetNodes();
    } // store iterator in object

  void Grid::ConsInterpReportFunctor::operator()( const HandleBox& a,
                                                  const HandleBox& b) {
      UInt destElemNum = *a.handle();
      UInt sourceNodeIndex = *b.handle();
      UInt dim = destGrid_->GetDim();
      UInt localDim;
      UInt numElemNodes;
      Elem::FEType type;
      RegionIdType region;
      Matrix<Double> coordMat;
      Matrix<Double> globCoordMat;
      Matrix<Double> localCoords;
      Vector<Double> point;
      StdVector<bool> coordsInside;
      Vector<Double> locCoords;
      const Elem* elem = NULL;

      //      std::cout << "Elem Number " << elemNum << " <- " << sourceNodeNum << std::endl;

      // If source to destination node map already contains entries
      // for the source node, refuse to do any further calculations
      // because weights already exist for this node.
      if(!consInterpWeights_[sourceNodeIndex].empty())
      {
        //        std::cout << "Rejecting mapping from source node index " << sourceNodeIndex << " to element " << destElemNum << std::endl;
        return;
      }
      

      // Fill coordinate matrix with node coords of destination element
      destGrid_->GetElemData(destElemNum, type, region, &connect_[0]);
      numElemNodes = Elem::GetNumElemNodes(type);
      coordMat.Resize(dim, numElemNodes);
      for(UInt i=0; i<numElemNodes; i++)
      {
        destGrid_->GetNodeCoordinate(point, connect_[i], true);
        for(UInt j=0; j<dim; j++) {
          coordMat[j][i] = point[j];
        }
      }

      // Fill global coordinate matrix for node of source grid
      globCoordMat.Resize(dim, 1);
      for(UInt j=0; j<dim; j++) {
        globCoordMat[j][0] = nodeCoords_[sourceNodeIndex][j];
      }

      // Get local coordinate of source node in respect to potential
      // destination element
      elem = destGrid_->GetElem(destElemNum);
      elem->ptElem->Global2LocalCoords(localCoords, globCoordMat, coordMat);
      elem->ptElem->CoordsInsideElem(localCoords, localEpsilon_, coordsInside);

      // If node is inside potential destination element, calculate
      // conservative interpolation weights.
      if(coordsInside[0])
      {
        localDim = localCoords.GetNumRows();
        locCoords.Resize(localDim);

        for(UInt j=0; j<localDim; j++) {
          locCoords[j] = localCoords[j][0];
        }

        // The vector S contains the values of the shape functions
        // at the local coordinate of the source node. These values
        // serve also as the interpolation weights.
        Vector<double> S;
        elem->ptElem->GetShFnc(S, locCoords, elem );

        //        std::cout << "Local Coord: " << locCoords << std::endl;
        //        std::cout << "Shape functions: " << S << std::endl;

        // Put weights into correct position of source -> destination
        // node map.
        for(UInt i=0; i<numElemNodes; i++)
        {
          UInt pos = destNodeNumToPosMap_[connect_[i]];

          if(consInterpWeights_[sourceNodeIndex].find(pos) ==
             consInterpWeights_[sourceNodeIndex].end())
          {
            if(S[i] != 0.0)
            {
              consInterpWeights_[sourceNodeIndex][pos] = S[i];
            }
            //            std::cout << "Node: " << connect_[i] << ": " << S[i] << std::endl;
          }
        }

        nodeCounter_++;
        percentage_ = (UInt)(100*(Double)nodeCounter_ / (Double)numSourceNodes_);
        if(((percentage_ % 10) == 0) && ((oldPercentage_ % 10) == 9))
        {
          //std::cout << percentage_ << "% done... " << std::endl;
          std::cout << "."; // use a short status display
        }
        oldPercentage_ = percentage_;
      }
#if 0
      else
      {
        std::cout << sourceNodeNum << ": Local Coord: " << localCoords[0][0] << " " << localCoords[1][0] << std::endl;
      }
#endif
    }


#endif // USE_INTERPOLATION

} // end of namespace
