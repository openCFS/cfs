// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <string>
#include <fstream>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "Elements/elements_header.hh"
#include "Domain/elem.hh"
#include "grid.hh"
#include "General/exception.hh"
#include "Utils/mathfunctions.hh"
#include "polygonIterator.hh"

#ifdef USE_INTERPOLATION
#include <CGAL/box_intersection_d.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2_algorithms.h>
#endif

namespace CoupledField
{

  // declare class specific logging stream
  DECLARE_LOG(grid)
  DEFINE_LOG(grid, "grid")

  Grid::Grid()
  {
    ENTER_FCN( "Grid::Grid" );

    ptQ1     = new Quad1FE();
    ptQ2     = new Quad2FE();
    ptTet1   = new Tetra1FE();
    ptTet2   = new Tetra2FE();
    ptL1     = new Line1FE();
    ptL2     = new Line2FE();
    ptTr1    = new Triangle1FE();
    ptTr2    = new Triangle2FE();
    ptHexa1  = new Hexa1FE();
    ptHexa2  = new Hexa2FE();
    ptPyra1  = new Pyra1FE();
    ptPyra2  = new Pyra2FE();
    ptWedge1 = new Wedge1FE();
    ptWedge2 = new Wedge2FE();

 #ifdef USE_SCRIPTING
    // Register functions
    RegisterFunctions();
 #endif

  }

  

  Grid::~Grid()
  {
    ENTER_FCN( "Grid::~Grid" );
    if (ptQ1)     delete ptQ1;
    if (ptQ2)     delete ptQ2;
    if (ptTet1)   delete ptTet1;
    if (ptTet2)   delete ptTet2;
    if (ptL1)     delete ptL1;
    if (ptL2)     delete ptL2;
    if (ptTr1)    delete ptTr1;
    if (ptTr2)    delete ptTr2;
    if (ptHexa1)  delete ptHexa1;
    if (ptHexa2)  delete ptHexa2;
    if (ptPyra1)  delete ptPyra1;
    if (ptPyra2)  delete ptPyra2;
    if (ptWedge1) delete ptWedge1;
    if (ptWedge2) delete ptWedge2;
  
  }


  void Grid::AddRegion(const std::string regionName, RegionIdType & rId)
  {
    rId = regionNames_.GetSize();
    regionNames_.Push_back(regionName);
  }
    
  void Grid::AddRegions(const StdVector<std::string> & regionNames,
                        StdVector<RegionIdType> & rIds)
  {
    UInt numRegions = regionNames.GetSize();
    UInt newId = regionNames_.GetSize();

    rIds.Resize(numRegions);
    for(UInt i=0; i<numRegions; i++, newId++)
    {
      rIds[i] = newId;
      regionNames_.Push_back(regionNames[i]);
    }
  }
    
  void Grid::AddRegions(const std::vector<std::string> & regionNames,
                        std::vector<RegionIdType> & rIds)
  {
    UInt numRegions = regionNames.size();
    UInt newId = regionNames_.GetSize();

    rIds.resize(numRegions);
    for(UInt i=0; i<numRegions; i++, newId++)
    {
      rIds[i] = newId;
      regionNames_.Push_back(regionNames[i]);
    }
  }

  UInt Grid::GetNumRegions()
  {
    return regionNames_.GetSize();
  }
    
  UInt Grid::GetNumVolRegions()
  {
    return volRegionIds_.GetSize();
  }
    
  UInt Grid::GetNumSurfRegions()
  {
    return surfRegionIds_.GetSize();
  }
    
  void Grid::GetRegionIds( StdVector<RegionIdType> & regions ) {
    ENTER_FCN( "Grid::GetRegionIds" );
    
    UInt numRegions = volRegionIds_.GetSize() + surfRegionIds_.GetSize();
    
    regions.Clear();

    for(UInt i=0; i<numRegions; i++)
      regions.Push_back(i);
  }
  
  void Grid::GetRegionIds( std::vector<RegionIdType> & regions ) {
    ENTER_FCN( "Grid::GetRegionIds" );
        
    UInt numRegions = volRegionIds_.GetSize() + surfRegionIds_.GetSize();
        
    regions.clear();
        
    for(UInt i=0; i<numRegions; i++)
      regions.push_back(i);
  }

  
  void Grid::GetVolRegionIds( StdVector<RegionIdType> & volRegions ) {
    ENTER_FCN( "Grid::GetVolRegionIds" );
    
    volRegions = volRegionIds_;
  }
  
  void Grid::GetVolRegionIds( std::vector<RegionIdType> & volRegions ) {
    ENTER_FCN( "Grid::GetVolRegionIds" );
        
    UInt numRegions = volRegionIds_.GetSize();
        
    volRegions.resize(numRegions);
        
    for(UInt i=0; i<numRegions; i++)
      volRegions[i] = volRegionIds_[i];
  }

  
  void Grid::GetSurfRegionIds( StdVector<RegionIdType> & surfRegions ) {
    ENTER_FCN( "Grid::GetSurfRegionIds" );
    
    surfRegions = surfRegionIds_;
  }

  void Grid::GetSurfRegionIds( std::vector<RegionIdType> & surfRegions ) {
    ENTER_FCN( "Grid::GetSurfRegionIds" );
        
    UInt numRegions = surfRegionIds_.GetSize();
        
    surfRegions.resize(numRegions);
        
    for(UInt i=0; i<numRegions; i++)
      surfRegions[i] = surfRegionIds_[i];
  }

  RegionIdType Grid::RegionNameToId( const std::string & regionName ) {
    ENTER_FCN( "Grid::RegionNameToId" );

    RegionIdType ret = NO_REGION_ID;

    if (regionName == "all" ) {
      ret= ALL_REGIONS;
    } else {
      ret =  regionNames_.Find(regionName);
      if (ret == -1 ) {
        EXCEPTION( "The region with name '" << regionName 
                   << "' is not contained in the grid!" );
      }
    }
    return ret;
  }

  void Grid::RegionIdToName( StdVector<std::string> & regionNames,
                             const StdVector<RegionIdType> & regionId ) {
    ENTER_FCN( "Grid::RegionIdToName" );
  
    regionNames.Resize( regionId.GetSize() );
    for (UInt i=0; i<regionId.GetSize(); i++ ) {
      if ( regionId[i] == ALL_REGIONS )
        regionNames[i] = "all";
      else
        regionNames[i] = regionNames_[regionId[i]];
    }
  }

  void Grid::RegionNameToId( StdVector<RegionIdType> & regionIds,
                             const StdVector<std::string> 
                             & regionNames ) {
    ENTER_FCN( "Grid::RegionNameToId" );

    RegionIdType ret = NO_REGION_ID;
    
    regionIds.Resize( regionNames.GetSize() );
    for (UInt i=0; i<regionNames.GetSize(); i++ ) {
    
      if (regionNames[i] == "all" ) {
        ret = ALL_REGIONS;
      } else {
        ret = regionNames_.Find(regionNames[i]);
        if ( ret == -1 ) {
          EXCEPTION("The region with name '" << regionNames[i]
                    << "' is not contained in the grid!" );
        }
      }
      regionIds[i] = ret;
    }
  }
  
  std::string Grid::RegionIdToName( const RegionIdType regionId ) {
    ENTER_FCN( "Grid::RegionIdToName" );
  
    if ( regionId == ALL_REGIONS )
      return "all";
    else
      return regionNames_[regionId];
  }
  

  shared_ptr<EntityList> Grid::GetEntityList( EntityList::ListType listType, 
                                              const std::string& name,
                                              EntityList::DefineType defineType ) {
    ENTER_FCN( "Grid::GetEntityList" );

    shared_ptr<EntityList> ret;

    if( listType == EntityList::ELEM_LIST ) {
      shared_ptr<ElemList> eList  = shared_ptr<ElemList>( new ElemList(this) );
      if( defineType == EntityList::REGION ) {
        RegionIdType regionId = RegionNameToId( name );
        eList->SetRegion( regionId);
      } else {
        eList->SetNamedElems( name );
      }
      ret = eList;
      
    } else if( listType == EntityList::SURF_ELEM_LIST ) {
      shared_ptr<SurfElemList> surfList  = 
        shared_ptr<SurfElemList>( new SurfElemList(this) );
      if( defineType == EntityList::REGION ) {
        RegionIdType regionId = RegionNameToId( name );
        surfList->SetRegion( regionId);
      } else {
        surfList->SetNamedElems( name );
      }
      ret = surfList;

    } else if( listType == EntityList::NODE_LIST ) {
      shared_ptr<NodeList> nodeList = shared_ptr<NodeList>( new NodeList(this) );
      // Check if name describes a nodeList
      if( defineType == EntityList::NAMED_NODES ) {
        StdVector<std::string> nodeNames;
        GetListNodeNames( nodeNames );
        nodeList->SetNamedNodes( name );
      } else if( defineType == EntityList::REGION ) {
        RegionIdType regionId = RegionNameToId( name );
        nodeList->SetNodesOfRegion( regionId );
      } else {
        EXCEPTION("GetEntityList with NODE_LIST works only with regions"
                  << " and named nodes!" );
      }
      ret = nodeList;
    } else if( listType == EntityList::REGION_LIST ) {
      shared_ptr<RegionList> regionList = 
        shared_ptr<RegionList>( new RegionList(this) );
      if( defineType == EntityList::REGION ) {
        RegionIdType regionId = RegionNameToId( name );
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
     StdVector<Elem*>   elems;;
    
     std::cout << "Grid: elements=" << GetNumElems() << " nodes=" << GetNumNodes() << std::endl;             

     for(UInt i = 0; i < regionNames_.GetSize(); i++) 
     {
       std::string  region_name = regionNames_[i];
       RegionIdType region_id   = RegionNameToId(region_name);
           
       GetElems(elems, region_id);
                   
       std::cout << "region: " << region_name << " id=" << region_id << " elements=" << elems.GetSize() <<  std::endl;            
     }
   }


  bool Grid::InitNonmatchingInterfaces() {
    StdVector<SurfElem*> slaveElems;
    StdVector<SurfElem*> masterElems;
    StdVector<std::string> ncRegionNames;
    StdVector<std::string> ncSlaveIfaceRegionNames;
    StdVector<std::string> ncMasterIfaceRegionNames;
    std::string slaveName;
    std::string masterName;
    std::string ncRegionName;
    RegionIdType slaveId;
    RegionIdType masterId;
    UInt n;
    bool coPlanarInterface;

    StdVector<SurfElem*> ifaceElems;

#ifdef USE_INTERPOLATION
    intersection();
#endif
    
    ParamNode* pNode;
    StdVector<ParamNode*> ncIfaceNodes;
    pNode = param->Get("domain");
    pNode = pNode->Get("ncInterfaceList",false);
    if(!pNode)
    {
      Info->PrintF("","\nNo non-conforming interfaces found!\n");
      LOG_DBG2(grid) << "NonMatching: No non-conforming interfaces found!"
                     << slaveName;
      return true;
    }

    ncIfaceNodes = pNode->GetList("ncInterface");

    for( UInt i = 0; i < ncIfaceNodes.GetSize(); i++ ) {
      ParamNode* ncIfaceNode = ncIfaceNodes[i];
        
      pNode = ncIfaceNode->Get("name");
      ncRegionNames.Push_back( pNode->AsString() );
      pNode = ncIfaceNode->Get("slaveSide");
      ncSlaveIfaceRegionNames.Push_back( pNode->AsString() );
      pNode = ncIfaceNode->Get("masterSide");
      ncMasterIfaceRegionNames.Push_back( pNode->AsString() );
    }

    n = ncRegionNames.GetSize();
    for(UInt i=0; i<n; i++)
    {
      slaveName = ncSlaveIfaceRegionNames[i];
      masterName = ncMasterIfaceRegionNames[i];
      ncRegionName = ncRegionNames[i];
      
      slaveId = RegionNameToId(slaveName);
      masterId = RegionNameToId(masterName);

      if((slaveId == -1) ||
         (masterId == -1))
      {
        EXCEPTION( "Could not find region id(s) for master/slave interface!");
      }
      
      LOG_DBG2(grid) << "Trying to get SLAVE Elems: " << slaveName << std::endl;
      GetSurfElems(slaveElems, slaveId);
      LOG_DBG2(grid) << "Trying to get MASTER Elems: " << masterName << std::endl;
      GetSurfElems(masterElems, masterId);


      if((slaveElems.GetSize() == 0) ||
         (masterElems.GetSize() == 0))
      {
        EXCEPTION( "Could not find elements on master/slave interface!");
      }

      
      LOG_DBG3(grid) << "NonMatching: Slave Element numbers:";
      for(UInt i=0; i<slaveElems.GetSize(); i++)
      {
        LOG_DBG3(grid) << "NonMatching: " << slaveElems[i]->elemNum;
        ifaceElems.Push_back(slaveElems[i]);
      }

      LOG_DBG3(grid) << "NonMatching: Master Element numbers:";
      for(UInt i=0; i<masterElems.GetSize(); i++)
      {
        LOG_DBG3(grid) << "NonMatching: " << masterElems[i]->elemNum;
        ifaceElems.Push_back(masterElems[i]);
      }
      
      coPlanarInterface = AreInterfacesCoplanar(ifaceElems);

      ifaceElems.Clear();
      slaveElems.Clear();
      masterElems.Clear();

      
      InitMasterInterface(ncRegionName, slaveName, masterName, coPlanarInterface);
    }
    
    return true;
  }


  bool Grid::InitMasterInterface(const std::string & ncRegionBaseName,
                                 const std::string & slaveInterfaceName,
                                 const std::string & masterInterfaceName,
                                 const bool coPlanarIface)
  {
    StdVector<SurfElem*> masterElems;
    StdVector<SurfElem*> slaveElems;
    StdVector<NCElem*> ncElems;
    StdVector<SurfElem*> ncElemsHelper;
    RegionIdType masterId = RegionNameToId(masterInterfaceName);
    RegionIdType slaveId = RegionNameToId(slaveInterfaceName);
    RegionIdType ncRegionId;
    StdVector<UInt> masterNodes;
    StdVector< UInt > ncElemIds;
    std::string isecCalc;
    bool lostDOFs = false;

    GetSurfElems(masterElems, masterId);
    GetSurfElems(slaveElems, slaveId);

    if(regionNames_.Find(ncRegionBaseName) != -1)
    {
      ncRegionId = RegionNameToId(ncRegionBaseName);
      ClearRegion(ncRegionId);
    }
    else
      AddSurfaceRegion(ncRegionBaseName, ncRegionId);

    ParamNode* ncIfaceNode;
    ncIfaceNode = param->Get("domain")->Get("ncInterfaceList");
    ncIfaceNode = ncIfaceNode->Get("ncInterface", "name", ncRegionBaseName);
    isecCalc = "polygon";
    ncIfaceNode->Get("isecCalculation", isecCalc, false);
    
    if(GetDim() == 2)
    {
      if (!coPlanarIface)
      {
        EXCEPTION( "Only collinear interfaces are supported at the moment!");
      }
      for(UInt i=0; i<masterElems.GetSize(); i++)
      {
        for(UInt j=0; j<slaveElems.GetSize(); j++)
        {
          SideOnSide(masterElems[i], slaveElems[j],
                     coPlanarIface, ncElems);
        }
      }
    }
    else //GetDim == 2
    {
      for(UInt i=0; i<masterElems.GetSize(); i++)
      {
        for(UInt j=0; j<slaveElems.GetSize(); j++)
        {
          if((isecCalc == "dirty") || (isecCalc == ""))
          {
            EXCEPTION("ISNMG: Dirty Intersection not implemented anymore!");
          }

          if((isecCalc == "coaxi"))
          {
	    if (!coPlanarIface) {
	      EXCEPTION("Only coplanar interfaces are supported with coaxi algorithm!");
	    }
            if((masterElems[i]->ptElem != ptQ1) ||
               (slaveElems[j]->ptElem != ptQ1))
            {
              EXCEPTION("Only linear Quads can be intersected with coaxi algorithm!");
            }
            
            if(RectangleOnRectangle(masterElems[i], slaveElems[j], coPlanarIface, ncElems))
              LOG_DBG3(grid) << "Intersection between " << masterElems[i]->elemNum << " and " << slaveElems[j]->elemNum << std::endl;
          }
	  else if (isecCalc == "polygon")
	  {
            if (((masterElems[i]->ptElem != ptTr1) &&
                 (masterElems[i]->ptElem != ptQ1)) ||
                ((slaveElems[j]->ptElem != ptTr1) &&
                 (slaveElems[j]->ptElem != ptQ1)))
            {
              EXCEPTION("Only linear triangles and quads can be intersected with polygon algorithm.");
            }

            if (PolygonOnPolygon(masterElems[i], slaveElems[j],
                                 coPlanarIface, ncElems))
              LOG_DBG3(grid) << "Intersection between " <<
                masterElems[i]->elemNum << " and " <<
                slaveElems[j]->elemNum << std::endl;
	  }
        }
      }
    }

    if(lostDOFs)
    {
      Warning( "Not all lagrange DOFs might be incorporated in system matrix!\n"
               "Try to set isecRefines to a higher value.",
               __FILE__ ,__LINE__);
    }

    if(ncElems.GetSize() == 0)
    {
      ClearRegion(ncRegionId);
    }
    else
    {
      ncElemsHelper.Resize(ncElems.GetSize());
      
      for(UInt i=0; i<ncElems.GetSize(); i++)
      {
        ncElemsHelper[i] = ncElems[i];
      }

      AddSurfaceElems( ncRegionId, ncElemsHelper, ncElemIds);
    }
    
    return true;
  }

  
  bool Grid::AreInterfacesCoplanar(const StdVector<SurfElem*>& ifaceElems)
  {
    std::set<Integer> ifaceNodes;
    std::set<Integer>::iterator it,end;
    Vector<Double> pv1, pv2, pv3;
    Vector<Double> v1, v2, normal, n;
    Double innerProd, norm1, norm2;
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

  /****************************************************************************
   ** 
   ** SideOnSide
   **
   **		computes the local coordinates of the overlap of the master
   **		and slave element with respect to the master side in the 
   **		order of the orientation of the slave side element
   **		without projection to any hyperplane. It pushes back the
   **           intersection element to elemList.
   **           This functions assumes a coplanar interface!
   **
   ** Input Parameters:   
   **		ifaceElemM:  Master Side
   **		ifaceElemS: Slave Side
   **		coPlanarIface: this parameter indicates wether interface is
   **                          coplanar
   **
   ** Output Parameters:
   **		elemList: the found intersection NCElems will be pushed
   **                     back to this vector 
   **		
   */


  bool Grid::SideOnSide (SurfElem *ifaceElem1,
                         SurfElem *ifaceElem2,
                         const bool coPlanarIface,
                         StdVector<NCElem*>& elemList)
  {
    // c0, c1, d0 and d1 are the endpoints of the two line elements
    //
    //           d0 x--------------------------x d1
    // c0 x---------|-----------x c1

    Vector<Double> c0, c1, d0, d1; // endpoint-coordinates of the two lines
    Vector<Double> diff0, diff1;
    Vector<Double> s, t;
    StdVector<UInt> connect2;
    Double dist, dist1, fac;
    
    s.Resize(2);
    t.Resize(2);
    connect2.Resize(2);
    
    // Get coordinates of the endpoints
    GetNodeCoordinate(c0, ifaceElem1->connect[0]);
    GetNodeCoordinate(c1, ifaceElem1->connect[1]);
    GetNodeCoordinate(d0, ifaceElem2->connect[0]);
    GetNodeCoordinate(d1, ifaceElem2->connect[1]);

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
    dist1 = diff0.NormL2();
    diff0.Inner(diff1, s[1]);
    s[1] *= fac;

    // Bring line2's endpoints into ascending order.
    if(s[1] < s[0])
    {
      t[0] = s[1];
      t[1] = s[0];
      connect2[0] = ifaceElem2->connect[1];
      connect2[1] = ifaceElem2->connect[0];      
    }
    else
    {
      t[0] = s[0];
      t[1] = s[1];
      connect2[0] = ifaceElem2->connect[0];
      connect2[1] = ifaceElem2->connect[1];
    }

    // Check if an intersection between line1 and line2
    // exists.
    if(t[0] >= 1.0)
      return false;

    if(t[1] <= 0.0)
      return false;


    NCElem* ncElem = new NCElem;
    ncElem->connect.Resize(2);

    // If an intersection exists, we must distinguish
    // 4 different cases.
    if(t[0] <= 0)
    {
      ncElem->connect[0] = ifaceElem1->connect[0];

      if(t[1] >= 1)
      {
        // connect2[0] x--------|--------------------|-----x connect2[1]
        //                   c0 x--------------------x c1

        ncElem->connect[1] = ifaceElem1->connect[1];
      }      
      else
      {
        // connect2[0] x--------|---------x connect2[1]
        //                   c0 x---------|-----------x c1

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
        
        ncElem->connect[1] = ifaceElem1->connect[1];
      }
      else
      {
        // connect2[0] x------------x connect2[1]
        //      c0 x---|------------|---x c1

        ncElem->connect[1] = connect2[1];              
      }
    }

    ncElem->ptElem = ifaceElem1->ptElem;
    ncElem->ptLagrangeParent = ifaceElem2;
    ncElem->ptSurfParent = ifaceElem1;

    elemList.Push_back(ncElem);

    return true;
  }
  
  bool Grid::RectangleOnRectangle(SurfElem *ifaceElem1,
                                  SurfElem *ifaceElem2,
                                  const bool coPlanarIface,
                                  StdVector<NCElem*>& elemList)
  {
    Vector<Double> c0, c1, c2, d0, d1;
    Vector<Double> diffS, diffX, diffY;
    Vector<Double> s, t;
    StdVector<UInt> connect2;
    Double distX, distY, facX, facY;
    UInt nodeNr;
    
    
    s.Resize(4);
    t.Resize(4);
    connect2.Resize(4);

    // The meaning of the points c0, c1, c2, d0 and d1
    // is as follows:
    //                x------------------x d1
    //                |                  |
    //                |                  |
    // c2 x-----------+----------x       |
    //    |           |          |       |
    //    |           |          |       |
    //    |        d0 x----------+-------x
    //    |                      |
    //    |                      |
    // c0 x----------------------x c1

    
    // Get coordinates of the endpoints
    GetNodeCoordinate(c0, ifaceElem1->connect[0]);
    GetNodeCoordinate(c1, ifaceElem1->connect[1]);
    GetNodeCoordinate(c2, ifaceElem1->connect[3]);
    GetNodeCoordinate(d0, ifaceElem2->connect[0]);
    GetNodeCoordinate(d1, ifaceElem2->connect[2]);

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
    
    // Now compute vector from c0 to d1 and project
    // the result onto the new x- and y-axis.
    diffS = d1 - c0;
    diffS.Inner(diffX, s[2]);
    diffS.Inner(diffY, s[3]);
    s[2] *= facX;
    s[3] *= facY;

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
        connect2[1] = ifaceElem2->connect[3];
        connect2[2] = ifaceElem2->connect[0];
        connect2[3] = ifaceElem2->connect[1];        
      }
      else
      {
        t[1] = s[1];
        t[3] = s[3];
        connect2[0] = ifaceElem2->connect[1];
        connect2[1] = ifaceElem2->connect[0];
        connect2[2] = ifaceElem2->connect[3];
        connect2[3] = ifaceElem2->connect[2];        
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
        connect2[0] = ifaceElem2->connect[3];
        connect2[1] = ifaceElem2->connect[2];
        connect2[2] = ifaceElem2->connect[1];
        connect2[3] = ifaceElem2->connect[0];
      }
      else
      {
        t[1] = s[1];
        t[3] = s[3];
        connect2[0] = ifaceElem2->connect[0];
        connect2[1] = ifaceElem2->connect[1];
        connect2[2] = ifaceElem2->connect[2];
        connect2[3] = ifaceElem2->connect[3];        
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

    ncElem->ptElem = ifaceElem1->ptElem;
    ncElem->ptLagrangeParent = ifaceElem2;
    ncElem->ptSurfParent = ifaceElem1;

    elemList.Push_back(ncElem);

    return true;
  }

  const double TOL = 1e-8;

  bool Grid::PolygonOnPolygon(SurfElem *ifElem1, SurfElem *ifElem2,
                              const bool coPlanarIface, StdVector<NCElem*> &elemList)
  {
    UInt i;
    StdVector< Vector<Double> > p1, p2, r;
      
    p1.Resize(ifElem1->connect.GetSize());
    for (i = 0; i < p1.GetSize(); ++i)
      GetNodeCoordinate(p1[i], ifElem1->connect[i]);
      
    p2.Resize(ifElem2->connect.GetSize());
    for (i = 0; i < p2.GetSize(); ++i)
      GetNodeCoordinate(p2[i], ifElem2->connect[i]);

    if (CutPolys(p1, p2, coPlanarIface, r))
    {
      UInt start = TriangulatePoly(r, elemList);

      for (UInt i = start; i < elemList.GetSize(); ++i)
      {
        elemList[i]->ptElem = ifElem1->ptElem;
        elemList[i]->ptLagrangeParent = ifElem2;
        elemList[i]->ptSurfParent = ifElem1;
        elemList[i]->ptElem = ptTr1;
      }
	  
      return true;
    }
    return false;
  }

  Grid::IntersectType Grid::CutLines(const Vector<Double> &a,
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

        if (fabs(l_ac + l_ad - l2) < TOL) {
          e = a;
          if (l_ac < TOL)
            return INTERSECT_A_EQ_C;
          if (l_ad < TOL)
            return INTERSECT_IN_D;
          return INTERSECT_IN_A;
        }
        if (fabs(l_bc + l_bd - l2) < TOL) {
          e = b;
          return INTERSECT_IN_B;
        }
        if (fabs(l_ac + l_bc - l1) < TOL) {
          e = c;
          return INTERSECT_IN_C;
        }
        if (fabs(l_ad + l_bd - l1) < TOL) {
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

    // lines are not parallel, so compute intersection
    Double h, k, k1, k2, denom1, denom2;

    /* a + h * v1 = c + k * v2
     * This is a system with 2 unknowns and 3 equations. Compute k1 from
     * equations 1 and 2, k2 from equations 1 and 3.
     */
    denom1 = v1[1] * v2[0] - v1[0] * v2[1];
    if (fabs(denom1) > TOL)
      k1 = (v1[0] * (c[1] - a[1]) + v1[1] * (a[0] - c[0])) / denom1;
    denom2 = v1[2] * v2[0] - v1[0] * v2[2];
    if (fabs(denom2) > TOL)
      k2 = (v1[0] * (c[2] - a[2]) + v1[2] * (a[0] - c[0])) / denom2;
    // If this system has no solution, lines do not intersect.
    if ((fabs(denom1) <= TOL) && (fabs(denom2) <= TOL))
      return INTERSECT_NONE;
    if ((fabs(denom1) > TOL) && (fabs(denom2) > TOL)) {
      if (fabs(k1 - k2) > TOL)
        return INTERSECT_NONE;
    }

    // compute second unknown
    if (fabs(denom1) > TOL) k = k1; else k = k2;
    if (fabs(v1[0]) > TOL)
      h = (c[0] - a[0] + v2[0] * k) / v1[0];
    else if (fabs(v1[1]) > TOL)
      h = (c[1] - a[1] + v2[1] * k) / v1[1];
    else
      h = (c[2] - a[2] + v2[2] * k) / v1[2];

    e = a + v1 * h; // compute point of intersection

    if (h > -TOL) { // we consider only [a,inf)
      if ((k > -TOL) && (k < 1.0 + TOL)) { // intersection on [c,d]?
        if (h < 1.0 + TOL) { // intersection on [a,b]?
          // treat special cases
          if (fabs(h - 1.0) < TOL) // h=1 means intersection in b
            return INTERSECT_IN_B;
          if (fabs(k - 1.0) < TOL) // k=1 means intersection in d
            return INTERSECT_IN_D;
          if (fabs(k) < TOL) { // k=0 means intersection in c
            if (fabs(h) < TOL) // h=0 means intersection in a
              return INTERSECT_A_EQ_C;
            return INTERSECT_IN_C;
          }
          if (fabs(h) < TOL) // h=0 means intersection in a
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
                      StdVector< Vector<Double> > &p2, const bool coPlanarIface,
                      StdVector< Vector<Double> > &r) const
  {
    Double r1, r2;
    UInt i, inside = 0, nCuts = 0, start_cur = p1.GetSize();
    Vector<Double> c1, c2, e/*, last_cut*/;
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
      
    //last_cut.Resize(3);

    // compute surrounding circles of both polygons
    r1 = PolyCentroid(p1, c1);
    r2 = PolyCentroid(p2, c2);

    // quit, if surrounding circles do not intersect
    temp1 = (c1 - c2);
    if (r1 + r2 < temp1.NormL2())
      return false;

    // if interface is not coplanar then project p1 onto p2
    if (!coPlanarIface) {
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
        IntersectType cuttype
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
      EXCEPTION("A line can not cut more than two edges of a convex polygon");
    }

    // save the position of the first cut in the active polygon
    pi1.SetBegin(pi1.GetPos());

    if (nCuts == 2) { // two cuts
      // make sure we do not treat a duplicate cut
      temp1 = (cuts[1].loc - cuts[0].loc);
      if (temp1.NormL2() < TOL) {
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
          //last_cut = cuts[1].loc;
          r.Push_back(cuts[1].loc);
		  
          pi2.Seek(cuts[0].index);
          if ((cuts[0].type != INTERSECT_IN_C) &&
              (cuts[0].type != INTERSECT_A_EQ_C)) ++pi2;
          pi2.SetBegin();

          pi2.Seek(cuts[1].index);
        } else {
          r.Push_back(cuts[1].loc);
          //last_cut = cuts[0].loc;
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
      //last_cut = cuts[0].loc;
      r.Push_back(cuts[0].loc);

      // avoid finding the same cut twice
      ++pi1;
      // continue with p2, if indicated
      if (cuts[0].swap) {
        pi1.Swap(pi2);
      } else {// [a,b] cuts into p2, so add b
        /*r.Push_back(last_cut);
          last_cut = *pi1;*/
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
            /*r.Push_back(last_cut);
              last_cut = e;*/
            r.Push_back(e);
            swap = true;
            break;
          case INTERSECT_IN_A:
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
        /*r.Push_back(last_cut);
          last_cut = *pi1;*/
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
    IntersectType s;
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
    if ( temp.NormL2() < TOL)
      return TRUE;

    // try intersecting [c,p] with each edge of the polygon
    do {
      s = CutLines(center, p, *pi, pi.Next(), e);
      if (s <= INTERSECT_OUTSIDE)
        continue;
      if ((s == INTERSECT_ON_LINE2) || (s == INTERSECT_IN_B)) {
        result = TRUE;
        break;
      }
      if ((s == INTERSECT_CROSS) || (s >= INTERSECT_IN_C)) {
        result = FALSE;
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
    c.Init(0.0);

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
    NCElem *ncElem, *ncElem2;
    UInt numPoints = p.GetSize();
    Vector<Double> temp1, temp2;

    /*      
            Vector<Double> point, midPoint;
            Double dist, maxDist = 0.0;

            midPoint.Resize(3);
            midPoint[0] = midPoint[1] = midPoint[2] = 0;

            for(UInt i=0; i<numPoints; i++)
            {
            midPoint += p[i];
            }

            midPoint *= 1.0/numPoints;

            for(UInt i=0; i<numPoints; i++)
            {
            dist = p[i][0] - midPoint[0] +
            p[i][1] - midPoint[1] +
            p[i][2] - midPoint[2];
            maxDist = maxDist > dist ? maxDist : dist;
            }

            if(maxDist < 0.001)
            {
            std::cerr << "maxDist = " << maxDist <<std::endl;
            return firstNo;
            }
    */
      
    switch (numPoints) {
    case 3:
      ncElem = new NCElem;
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
      ncElem2 = new NCElem;
      ncElem->connect.Resize(3);
      ncElem2->connect.Resize(3);
	      
      temp1 = (p[0] - p[2]);
      temp2 = (p[1] - p[3]);
      if (temp1.NormL2() < temp2.NormL2())
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
      tri.Push_back(ncElem2);
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
        ncElem->connect.Resize(3);
		  
        ncElem->connect[0] = nodeNo;
        AddNode(p[i], nodeNo);
        ncElem->connect[1] = nodeNo;
        ncElem->connect[2] = centerNode;
		  
        tri.Push_back(ncElem);
      }
	      
      ncElem = new NCElem;
      ncElem->connect.Resize(3);

      ncElem->connect[0] = nodeNo;
      ncElem->connect[1] = firstNode;
      ncElem->connect[2] = centerNode;

      tri.Push_back(ncElem);
    }
    return firstNo;
  }

  // =======================================================================
  // Method wrappers for scripting
  // =======================================================================
    
  void Grid::Wrap_GetNodesByName() {
    SCRIPT_GET(std::string, name);
    StdVector<UInt> nodeNrs;
    GetNodesByName( nodeNrs, name);
    UInt2String(SCRIPT_RETVAL, nodeNrs); 
  }
  
  void Grid::Wrap_GetNodeCoordinate() {
    SCRIPT_GET(UInt, nodeNr);
    Vector<Double> coord;
    GetNodeCoordinate(coord, nodeNr );
    Double2String(SCRIPT_RETVAL, coord);
  }

  void Grid::Wrap_GetNodesByRegion() {
    SCRIPT_GET(std::string, name);
    StdVector<UInt> nodeNrs;
    GetNodesByRegion(nodeNrs, RegionNameToId(name) );
    UInt2String(SCRIPT_RETVAL, nodeNrs);
  }
  void Grid::Wrap_GetListNodeNames() {
    GetListNodeNames( SCRIPT_RETVAL );
  }

  void Grid::Wrap_GetListElemNames() {
    GetListElemNames( SCRIPT_RETVAL );
  }

  void Grid::Wrap_GetRegionNames() {
    StdVector<RegionIdType> regionIds;
    GetRegionIds( regionIds );
    RegionIdToName( SCRIPT_RETVAL, regionIds );

  }

  void Grid::Wrap_GetNumNodes() {
    SCRIPT_RETVAL.Push_back( GenStr( GetNumNodes() ) );
  }

  void Grid::Wrap_GetNumElems() {
    SCRIPT_RETVAL.Push_back( GenStr( GetNumElems() ) );

  }

  void Grid::Wrap_GetNumSurfElems() {
    SCRIPT_RETVAL.Push_back( GenStr( GetNumSurfElems() ) );

  }

  void Grid::Wrap_GetNumVolElems() {
    SCRIPT_RETVAL.Push_back( GenStr( GetNumVolElems() ) );
  }

  void Grid::Wrap_GetNumNodesOfRegion() {
    SCRIPT_GET(std::string, name);
    StdVector<RegionIdType> ids;
    ids.Push_back( RegionNameToId (name ) );
    SCRIPT_RETVAL.Push_back( GenStr( GetNumNodes ( ids ) ) );
  }

  void Grid::Wrap_GetNumElemsOfRegion() {
    SCRIPT_GET(std::string, name);             
    StdVector<RegionIdType> ids;
    ids.Push_back( RegionNameToId (name ) );
    SCRIPT_RETVAL.Push_back( GenStr( GetNumElems (ids ) ) );
  }
  
  
  void Grid::RegisterFunctions() {
    ENTER_FCN( "Grid::RegisterFunctions" );

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
  typedef CGAL::Box_intersection_d::Box_d<int,2> Box;

  //  typedef CGAL::Box_intersection_d::Box_d<double,2> Box2D;
  //  typedef CGAL::Box_intersection_d::Box_d<double,3> Box3D;
  //  typedef CGAL::Bbox_2                              BBox2D;
  typedef CGAL::Bbox_3                              BBox3D;
  typedef CGAL::Box_intersection_d::Box_with_handle_d<double,3,UInt*> HandleBox;

  typedef CGAL::Cartesian<double> K;
  typedef K::Point_2 Point2D;
  
  // coordinates for 9 boxes of a grid
  int p[9*4]   = { 0,0,1,1,  1,0,2,1,  2,0,3,1, // lower
                   0,1,1,2,  1,1,2,2,  2,1,3,2, // middle
                   0,2,1,3,  1,2,2,3,  2,2,3,3};// upper
  // 9 boxes
  Box boxes[9] = { Box( p,    p+ 2),  Box( p+ 4, p+ 6),  Box( p+ 8, p+10),
                   Box( p+12, p+14),  Box( p+16, p+18),  Box( p+20, p+22),
                   Box( p+24, p+26),  Box( p+28, p+30),  Box( p+32, p+34)};
  // 2 selected boxes as query; center and upper right
  Box query[2] = { Box( p+16, p+18),  Box( p+32, p+34)};

  // callback function object writing results to an output iterator
  template <class OutputIterator>
  struct Report {
    OutputIterator it;
    Report( OutputIterator i) : it(i) {} // store iterator in object
    // We write the id-number of box a to the output iterator assuming
    // that box b (the query box) is not interesting in the result.
    void operator()( const Box& a, const Box&) { *it++ = a.id(); }
  };
  template <class Iter> // helper function to create the function object
  Report<Iter> report( Iter it) { return Report<Iter>(it); }


  // callback function object writing results to an output iterator
  template <class OutputIterator, class Grid>
  struct Report2 {
    OutputIterator it;
    Grid& grid;
    Report2(OutputIterator i, Grid& g) : it(i), grid(g) {} // store iterator in object
    // We write the id-number of box a to the output iterator assuming
    // that box b (the query box) is not interesting in the result.
    void operator()( const HandleBox& a, const HandleBox& b) {
      UInt elemNum = *a.handle();
      UInt dim = grid.GetDim();
      UInt numElemNodes;
      FEType type;
      RegionIdType region;
      std::vector<UInt> connect;
      Matrix<Double> coordMat;
      Matrix<Double> globCoordMat;
      Matrix<Double> localCoords;
      Vector<Double> point;

      *it++ = elemNum;
      connect.resize(64);
      
      std::cout << "Elem Number " << elemNum << std::endl;

      grid.GetElemData(elemNum, type, region, &connect[0]);
      numElemNodes = NUM_ELEM_NODES[type];
      coordMat.Resize(dim, numElemNodes);
      globCoordMat.Resize(dim, 1);

      for(UInt i=0; i<numElemNodes; i++)
      {
        grid.GetNodeCoordinate(point, connect[i], true);
        // coordMat auffüllen!
        for(UInt j=0; j<dim; j++)
        {
          coordMat[j][i] = point[j];

          //          std::cout << "Corner " << (i+1) << " Coord "
          //                    << (j+1) << " " << point[j] << std::endl;
        }
      }

      for(UInt j=0; j<dim; j++)
      {
        globCoordMat[j][0] = b.min_coord(j);

        //        std::cout << "Glob Coord " << (j+1)
        //                  << " Coord " << b.min_coord(j) << std::endl;

      }

      grid.GetElem(elemNum)->ptElem->Global2LocalCoords(localCoords, globCoordMat, coordMat);

    }
  };
  template <class Iter, class Grid> // helper function to create the function object
  Report2<Iter, Grid> report2(Iter it, Grid& g) { return Report2<Iter, Grid>(it, g); }

  
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
  
#endif // USE_INTERPOLATION

} // end of namespace
