#include <string>
#include <fstream>

#include "Elements/elements_header.hh"
#include "Domain/elem.hh"
#include "grid.hh"
#include "General/exception.hh"

namespace CoupledField
{

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
  
  
} // end of namespace
