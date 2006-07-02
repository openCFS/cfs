#include <string>
#include <fstream>

#include "Elements/elements_header.hh"
#include "DataInOut/filetype.hh"
#include "Domain/elem.hh"
#include "grid.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"

namespace CoupledField
{

  Grid::Grid(FileType * aptFileType)
  {
    ENTER_FCN( "Grid::Grid" );

    ptFileType = aptFileType;


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

#ifdef TCL_INTERFACE
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




  RegionIdType Grid::RegionNameToId( const std::string & regionName ) {
    ENTER_FCN( "Grid::RegionNameToId" );

    RegionIdType ret = NO_REGION_ID;

    if (regionName == "all" ) {
      ret= ALL_REGIONS;
    } else {
      ret =  regionNames_.Find(regionName);
      if (ret == -1 ) {
        (*error) << "The region with name '" << regionName 
                 << "' is not contained in the grid!";
        Error( __FILE__, __LINE__ );
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
          (*error) << "The region with name '" << regionNames[i]
                   << "' is not contained in the grid!";
          Error( __FILE__, __LINE__ );
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
  

  shared_ptr<EntityList> Grid::GetEntityList( EntityList::Type type, 
                                             const std::string& name ) {
    ENTER_FCN( "Grid::GetEntityList" );

    shared_ptr<EntityList> ret;

    if( type == EntityList::ELEM_LIST ) {
      shared_ptr<ElemList> eList  = shared_ptr<ElemList>( new ElemList(this) );
      RegionIdType regionId = RegionNameToId( name );
      eList->SetRegion( regionId);
      ret = eList;
      
    } else if( type == EntityList::SURF_ELEM_LIST ) {
      shared_ptr<SurfElemList> surfList  = 
        shared_ptr<SurfElemList>( new SurfElemList(this) );
      RegionIdType regionId = RegionNameToId( name );
      surfList->SetRegion( regionId);
      ret = surfList;

    } else if( type == EntityList::NODE_LIST ) {
      shared_ptr<NodeList> nodeList = shared_ptr<NodeList>( new NodeList(this) );

      // Check if name describes a nodeList
      StdVector<std::string> nodeNames;
      GetListNodeNames( nodeNames );
      if( nodeNames.Find(name) < 0 ) {
        *error << "'" << name << "' describes no named nodes!";
        Error( __FILE__, __LINE__ );
      }
      nodeList->SetNodes( name );
      ret = nodeList;
    } else {
      *error << "Type '" << type << "' describes no EntityList which is created "
             << "by the grid-class.";
      Error( __FILE__, __LINE__ );
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
