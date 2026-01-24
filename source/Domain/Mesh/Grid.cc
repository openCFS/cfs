#include "Grid.hh"
#include "NcInterfaces/BaseNcInterface.hh"
#include "NcInterfaces/MortarInterface.hh"


#include <cmath>
#include <string>
#include <limits>
#include <boost/scoped_array.hpp>

#ifdef USE_LIBFBI
#include <fbi/tuplegenerator.h> //TraitsGenerator
#include <fbi/fbi.h> //SetA::intersect
#include <fbi/tuple.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include "General/Exception.hh"
#include "Utils/mathfunctions.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ProgramOptions.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Domain/Domain.hh"

#include "Domain/CoordinateSystems/CoordSystem.hh"

namespace CoupledField
{

  // declare class specific logging stream
  DEFINE_LOG(grid, "grid")

  Grid::Grid(PtrParamNode param, PtrParamNode infoNode)
  {
    isInitialized_ = false; // set by FinishInit()
    isAxi_ = false;
    depth2dPlane_ = 1.0;
    param_ = param;
    info_ = infoNode;
    
    region_.SetName("Grid::region");
    region_.Add(ALL_REGIONS, "all");
    
    integScheme_.reset( new IntScheme() );

    // in addition, add always the NO_REGION to the enum
    region_.Add( NO_REGION_ID, "_NO_REGION_");
  }

  Grid::~Grid()
  {
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

  Matrix<double> Grid::CalcRegionsBoundingBox(const StdVector<RegionIdType>& regs, CoordSystem* sys) 
  {
    Matrix<double> minMax;
    CalcBoundingBoxOfRegion(regs[0], minMax, sys);

    LOG_DBG(grid) << "CRBB regs=" << regs.ToString() << " minMax=" << minMax;

    for(unsigned int i = 1; i < regs.GetSize(); i++)
    {
      Matrix<double> tmp;
      CalcBoundingBoxOfRegion(regs[i], tmp, sys);
      for(unsigned int d = 0; d < dim_; d++ ) 
      {
        minMax[d][0] = std::min(minMax[d][0],tmp[d][0]);
        minMax[d][1] = std::max(minMax[d][1],tmp[d][1]);
      }
    }
    return minMax;
  }

  Matrix<double>& Grid::CalcGridBoundingBox(CoordSystem* sys, bool force_3D)
  {
    #pragma omp critical
    {
      Matrix<double>& box = grid_bounding_box_;
      if(box.GetNumRows() == 0)
      {
        // set the box ignoring force_3D!
        if(sys == NULL)
          sys = domain->GetCoordSystem();

        StdVector<RegionIdType> regs;
        GetVolRegionIds(regs);

        Matrix<double> tmp;

        for(unsigned int r = 0; r < regs.GetSize(); r++)
        {
          CalcBoundingBoxOfRegion(regs[r], tmp, sys);

          LOG_DBG(grid) << "CGBB: tmp rows= " << tmp.GetNumRows() << " cols = " << tmp.GetNumCols();
          LOG_DBG(grid) << "CGBB: " << r << " regs[r]reg=" << regs[r] << " = " << region_.ToString(regs[r]) << " bb=" << tmp.ToString();
          if(r == 0) // the first region is the first guess
            box = tmp;
          else
          {
            for(unsigned int d = 0; d < tmp.GetNumRows(); d++)
            {
              box[d][0] = std::min(box[d][0], tmp[d][0]);
              box[d][1] = std::max(box[d][1], tmp[d][1]);
            }
          }
        }
      }

      // now the box is set but it might be that force_3D is ignored
      // this also works if box was created in a previous call but with another force_3D parameter
      if(GetDim() == 2 && ((!force_3D && box.GetNumRows() == 3) || (force_3D && box.GetNumRows() == 2)))
      {
        Matrix<double> tmp((force_3D ? 3 : 2), 2);
        tmp.Assign(box, 1.0, true); // size tolerant
        box = tmp;
      }
    } // end of critical guard

    return grid_bounding_box_;
  }

  shared_ptr<ElemShapeMap> Grid::GetElemShapeMap(const Elem *ptElem, bool isUpdated)
  {
    shared_ptr<ElemShapeMap> ret(new LagrangeElemShapeMap(this));
    ret->SetElem(ptElem, isUpdated);
    return ret;
  }

  RegionIdType Grid::AddRegion(const std::string& name, bool reg)
  {
    RegionData rd;
    rd.name = name;
    rd.regular = reg;
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
      numSurfElemNodes_.Push_back(0);
    }
    else
    {
      regionData[id].type_idx = volRegionIds_.GetSize();
      volRegionIds_.Push_back(id);
      volElems_.Push_back(dummy_elems);
      numVolElemNodes_.Push_back(0);
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

  UInt Grid::GetNumNodes(const StdVector<RegionIdType>& regions) const
  {
    UInt numNodes = 0;

    for (UInt i=0; i < regions.GetSize(); i++)
      numNodes += GetNumNodes(regions[i]);
    return numNodes;
  }


  bool Grid::IsRegionRegular(StdVector<RegionIdType>& regions) const
  {
    for(unsigned int i = 0; i < regions.GetSize(); i++)
      if(!regionData[regions[i]].regular)
        return false;

    return true;
  }

  bool Grid::IsGridRegular() const
  {
    for(unsigned int i = 0; i < regionData.GetSize(); i++)
      if(!regionData[i].regular)
        return false;

    return true;
  }

  std::string Grid::GetRegionName(RegionIdType id )
  {
    for(UInt i = 0; i < regionData.GetSize(); i++) {
      if(regionData[i].id == id) {
        return regionData[i].name;
      }
    }
    EXCEPTION( "No Region name found for id " << id ) ;
    return std::string("");
  }

  void Grid::GetRegionNames( StdVector<std::string>& regionNames )
  {
    regionNames.Resize(regionData.GetSize());

    for(UInt i = 0; i < regionData.GetSize(); i++)
      regionNames[i] = regionData[i].name;
  }

  RegionIdType Grid::GetRegionId(const std::string name ){
      for(UInt i = 0; i < regionData.GetSize(); i++) {
          if(regionData[i].name == name) {
              return regionData[i].id;
          }
      }
      EXCEPTION( "No Region found with name '" << name << "'" ) ;
      return -1;
  }

  void Grid::GetVolRegionIds( StdVector<RegionIdType> & volRegions ) {
    volRegions = volRegionIds_;
  }

  void Grid::GetSurfRegionIds( StdVector<RegionIdType> & surfRegions ) {

    surfRegions = surfRegionIds_;
  }
  

  const Elem* Grid::GetElemAtGlobalCoord(const Vector<double>& globCoord,
                                         LocPoint& locCoord,
                                         const StdVector<shared_ptr<EntityList> >& srcEntities,
                                         bool printWarnings,
                                         bool updatedGeo) {
    
    StdVector<Vector<Double> > globCoords(1);
    StdVector<LocPoint> lps;
    StdVector<const Elem*> elems;
    globCoords[0] = globCoord;
    globCoords[0].Resize(GetDim());
    GetElemsAtGlobalCoords( globCoords, lps, elems, srcEntities, 1e-3, 1e-2, true, updatedGeo);
    if( elems.GetSize() == 0 && printWarnings ) {
      WARN( "Could not find element at global position " << globCoord.ToString() );
    }
    locCoord = lps[0];
    return elems[0];
  }

  void Grid::GetElemsAtGlobalCoords( const StdVector<Vector<Double> >& globCoords,
                                     StdVector< LocPoint >& localCoords,
                                     StdVector< const Elem* > & elems,
                                     const StdVector<shared_ptr<EntityList> >& srcEntities,
                                     Double globalTol, Double localTol,
                                     bool printWarnings,
                                     bool updatedGeo) {

    // 1) first, determine element candidates for each point, determined by
    //    intersection of bounding-boxes. The algorithm used depends on the
    //    library used (CGAL, own one, libfbi)
    const UInt numPts = globCoords.GetSize();
    StdVector<PointElemMatch> matches( numPts );
    for( UInt i = 0; i < numPts; ++i ) {
     matches[i].globCoord = globCoords[i];
    }
    
    MapPointsToBoundingBoxes( matches, srcEntities, globalTol, updatedGeo );

    // Debug information about found macthes
//    std::cerr << "Found the following matches:\n";
//    for( UInt i = 0; i < numPts; ++i ) {
//      std::cerr << "coord: " <<matches[i].globCoord.ToString() << std::endl;
//      std::cerr << "matches # " << matches[i].matches.size() << ":\n\t";
//      std::set<const Elem*>::const_iterator it = matches[i].matches.begin();
//      for( ; it != matches[i].matches.end(); ++it ) {
//        std::cerr << (*it)->elemNum << ", ";
//      }
//      std::cerr<< "\n";
//    }// loop over all matches
    
    // 2) Afterwards loop over all candidates 
    MapGlobPointsToLoc( matches, elems, localCoords, localTol, printWarnings, updatedGeo );
  }
  
  const Elem* Grid::GetElemAtNode( UInt nodeNum,
                                   LocPoint& locCoord,
                                   const std::set<RegionIdType>& srcRegions ) {
    const Elem* ret = NULL;
    // check if nodes were already mapped
    if( !midNodeProjections_.size() ) {
      // Perform mapping of mid-side nodes
      MapMidSideNodes();
    }
    boost::unordered_map<UInt, NodeElemMatch>::const_iterator it;
    it = midNodeProjections_.find(nodeNum); 
    if( it != midNodeProjections_.end() ) {
      const NodeElemMatch& matches = it->second;
      const UInt size = matches.GetSize();
      for( UInt i = 0; i < size; ++i ) {
        if( srcRegions.find(matches[i].first->regionId) != srcRegions.end() ) {
          ret = (matches[i]).first;
          locCoord = (matches[i]).second;
          break;
        }
      }
    }
    return ret;
  }
  
  
  UInt Grid::SetElementBarycenters(RegionIdType reg, bool updated)
  {
    RegionData& rd = regionData[reg];

    if(rd.barycenters)
      return 0;

    // our operation target
    StdVector<Elem*>& elems = rd.type == VOLUME_REGION ? volElems_[rd.type_idx] : surfElems_[rd.type_idx];
    for(UInt i = 0;  i < elems.GetSize(); i++)
      GetElemShapeMap(elems[i], updated)->CalcBarycenter(elems[i]->extended->barycenter);

    rd.barycenters = true; // don't do it again!

    return elems.GetSize();
  }

  shared_ptr<EntityList> Grid::GetEntityList( EntityList::ListType listType,
                                              const std::string& name ) {

    // First check, if there any entites with this name at all
    if( nameTypeMap_.find( name) == nameTypeMap_.end() ) {
      EXCEPTION( "There are no entities with name '" << name
                 << "' in the mesh" ) ;
    }

    EntityList::DefineType entityType = nameTypeMap_[name];

    shared_ptr<EntityList> ret;

    // check, if name denotes a list of surface elements
    bool isSurface = false;
    if( GetEntityDim(name) == GetDim()-1 ) {
      isSurface = true;
    }
    
    if( listType == EntityList::ELEM_LIST ) {
      shared_ptr<ElemList> eList;
      if( isSurface ) {
        eList.reset(new SurfElemList(this) );
      } else {
        eList.reset( new ElemList(this) );
      }
      if( entityType == EntityList::REGION ) {
        RegionIdType regionId = GetRegion().Parse( name );
        eList->SetRegion( regionId);
      } else {
        eList->SetNamedElems( name );
      }
      ret = eList;

    } else if( listType == EntityList::SURF_ELEM_LIST ) {
      shared_ptr<SurfElemList> surfList ( new SurfElemList(this) );
      if( entityType == EntityList::REGION ) {
        RegionIdType regionId = GetRegion().Parse( name );
        surfList->SetRegion( regionId);
      } else {
        surfList->SetNamedElems( name );
      }
      ret = surfList;

    } else if( listType == EntityList::NODE_LIST ) {
      shared_ptr<NodeList> nodeList ( new NodeList(this) );
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
      shared_ptr<RegionList> regionList ( new RegionList(this) );
      if( entityType == EntityList::REGION ) {
        RegionIdType regionId = GetRegion().Parse( name );
        regionList->SetRegion( regionId );
      } else if (entityType){
        EXCEPTION( "GetEntityList with REGION_LIST works only with regions!" );
      }
      ret = regionList;
    } else if( listType == EntityList::NAME_LIST ) {
      shared_ptr<NameList> nameList =
          shared_ptr<NameList>( new NameList(this) );
      nameList->SetName( name );
      ret = nameList;
    } else {
      EXCEPTION( "Type '" << listType << "' describes no EntityList which is created "
                 << "by the grid-class." );
    }

    return ret;

  }
  
  EntityList::DefineType Grid::GetEntityType( const std::string& name ) const {
    EntityList::DefineType ret = EntityList::NO_TYPE;
    std::map<std::string, EntityList::DefineType>::const_iterator it;
    it = nameTypeMap_.find(name);
    if( it != nameTypeMap_.end() ){
      ret = it->second;
    }
    return ret;  
  }
  
  UInt Grid::GetEntityDim( const std::string& name ) const {
    UInt dim = 0;
    std::map<std::string, UInt>::const_iterator it = entityDim_.find(name);
    if( it != entityDim_.end() ){
      dim = it->second;
    } else {
      EXCEPTION( "No entities with name '" << name << "' are defined in the grid")
    }
    return dim;  
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

  // =======================================================================
  // Automatic Layer Generation and Geometry Computation
  // =======================================================================
  void Grid::GetGeometryOnRegionNodes(shared_ptr<NodeGeometry>& geometry, const PtrParamNode layerGenNode, const RegionIdType& surfRegionId) {
    // check if geometry already exists for given region
    if (geometryRegionMap_.count(surfRegionId) == 0) {
      // directly compute geometry of passed surface region
      this->ComputeGeometryOnSurfaceRegionNodes(layerGenNode);
    }
    // pass geometry
    geometry = geometryRegionMap_[surfRegionId];
  };

  void Grid::GetConnectedSurfaceRegions(StdVector<RegionIdType>& connecedSurfRegionIds, const RegionIdType& volumeRegionId) {
    // check if the connection has been set
    if (volumeSurfaceRegionMap_.count(volumeRegionId) == 0) {
      EXCEPTION("Connection map must be set manually before calling Grid::GetConnectedSurfaceRegions().");
    } else {
      connecedSurfRegionIds = volumeSurfaceRegionMap_[volumeRegionId];
    }
  };

  // =======================================================================
  // FINITE VOLUME REPRESENTATION SECTION
  // =======================================================================

  Grid::FiniteVolumeRepresentation::FiniteVolumeRepresentation() {
    isSet = false;
  }

  Grid::FiniteVolumeRepresentation& Grid::GetFiniteVolumeRepresentation() {
    return fvr_;    
  }

  // =========================================================================
  // NONCONFORMING INTERFACES SECTION
  // =========================================================================

  void Grid::InitNcInterfacesFromXML() {
    // if no param object is present, just leave
    if (!param_) return;

    // check if there is a ncInterfaceList, if not just leave
    PtrParamNode nciListNode = param_->Get("domain")->Get("ncInterfaceList", ParamNode::PASS);
    if (!nciListNode) return;

    ParamNodeList nciList = nciListNode->GetList("ncInterface");
    UInt numNCIs = nciList.GetSize();
    ncInterfaces_.Reserve(numNCIs);

    // loop twice to ensure that static interfaces are added before actively moving ones,
    // before passively moving ones. This is to assure that changing interfaces are deleted in correct order.
    // In case of connectedRegions the passively moving interfaces must already contain the node offset when they
    // are updated.
    bool passiveUpdate = false;
    for ( UInt i=0; i<numNCIs; ++i ) {
      passiveUpdate = (nciList[i]->Get("passiveGeomUpdate")->As<std::string>()=="yes");
      if(!nciList[i]->Has("rotation") && !nciList[i]->Has("generalMotion") && !passiveUpdate)
        AddNcInterface(shared_ptr<BaseNcInterface>(new MortarInterface(this, nciList[i])));
    }
    for ( UInt i=0; i<numNCIs; ++i ) {
      passiveUpdate = (nciList[i]->Get("passiveGeomUpdate")->As<std::string>()=="yes");
      if((nciList[i]->Has("rotation") || nciList[i]->Has("generalMotion")) && !passiveUpdate)
        AddNcInterface(shared_ptr<BaseNcInterface>(new MortarInterface(this, nciList[i])));
    }
    for ( UInt i=0; i<numNCIs; ++i ) {
      passiveUpdate = (nciList[i]->Get("passiveGeomUpdate")->As<std::string>()=="yes");
      if(passiveUpdate)
        AddNcInterface(shared_ptr<BaseNcInterface>(new MortarInterface(this, nciList[i])));
    }
    // sanity check
    assert(ncInterfaces_.GetSize() == numNCIs);
  }

  shared_ptr<BaseNcInterface> Grid::GetNcInterface(NcInterfaceId ncId) const {
    if ( ncId < ncInterfaces_.GetSize() ) {
      return ncInterfaces_[ncId];
    } else {
      EXCEPTION("NcInterface with ID " << ncId << " is unknown.");
    }
  }

  Grid::NcInterfaceId Grid::GetNcInterfaceId(const std::string &name) const {
    std::map< std::string, NcInterfaceId >::const_iterator ncId
        =nciNameMap_.find(name);
    if ( ncId != nciNameMap_.end() ) {
      return ncId->second;
    } else {
      EXCEPTION("NcInterface with name '" << name << " is unknown.");
    }
  }

  Grid::NcInterfaceId Grid::AddNcInterface(shared_ptr<BaseNcInterface> ncIf) {
    ncInterfaces_.Push_back(ncIf);
    if ( ncIf->GetName().length() > 0 ) {
      nciNameMap_[ncIf->GetName()] = ncInterfaces_.GetSize()-1;
    }
    return ncInterfaces_.GetSize()-1;
  }

  void Grid::MoveNcInterfaces() {
    // First, we reset the moving interfaces. MortarInterface::ResetInterface() checks if the interface needs to be reset first.
    // We reset the interfaces in reverse order, so we always only delete nodes and elements that were created last
    // the double loop in Grid::InitNcInterfacesFromXML() assures that moving interfaces are updated last.
    // also, we need to assure that actively moving interfaces are updated before passively moving ones.
    for (auto it = ncInterfaces_.End(); it-- != ncInterfaces_.Begin();) {
      (*it)->ResetInterface();
    }
    for (auto it = ncInterfaces_.Begin(); it != ncInterfaces_.End(); ++it) {
      (*it)->UpdateInterface();
    }
  }

  bool Grid::HasNCI() {
    if (!ncInterfaces_.IsEmpty())
      return true;
    else 
      return false;
  }

  bool Grid::IsSurfacePlanar(const StdVector<SurfElem*>& ifaceElems) const
  {
    std::set<Integer> ifaceNodes;
    std::set<Integer>::iterator it,end;
    Vector<Double> pv1, pv2, pv3;
    Vector<Double> v1, v2, normal, n;
    Double innerProd = 0.0, norm1 = 0.0, norm2 = 0.0;
    Double eps = NORM_EPS;
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
        norm1 = v1.Normalize();
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
        norm2 = v2.Normalize();
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
            v1.CrossProduct(v2, normal);

            // We may not use a zero-normal vector to perform our
            // test for coplanarity.
            if(normal.NormL2() < eps)
            {
              normal[0] = 0.0;
              normal[1] = 0.0;
              normal[2] = 0.0;
              continue;
            }
            normal.Normalize();
            continue;
          }
          v1.CrossProduct(v2, n);

          // If the norm of n is smaller than eps we have multiplied
          // linearly dependant vectors -> a point in the plane
          // has been found.
          if(n.NormL2() < eps)
          {
            continue;
          }

          n.Normalize();

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
    // UInt lastCornerInRegion2;

    for(UInt i=0; i<nElemsRegion1; i++) {
      el = region1Elems[i];
      numCorners = Elem::shapes[el->type].numVertices;
      // lastCornerInRegion2 = 0;
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
        surfEl->type = Elem::ET_LINE2;
        surfEl->extended = new ExtendedElementInfo;

        switch(el->type) {
        case Elem::ET_TRIA6:
        case Elem::ET_QUAD8:
        case Elem::ET_QUAD9:
          surfEl->connect[2] = el->connect[n+numCorners];
          surfEl->type = Elem::ET_LINE3;
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

  StdVector<UInt> Grid::GetRegularDiscretization(RegionIdType region)
  {
    StdVector<UInt> n(3);
    n.Init(0);
    if(!IsRegionRegular(region))
      return n;

    StdVector<double> min(3);
    StdVector<double> max(3);
    UInt dim = this->GetDim();
    min.Init(1e10);
    max.Init(-1e10);

    StdVector <Elem*> elems;
    this->GetElems(elems,region);
    this->SetElementBarycenters(region,false);

    for (UInt i = 0; i < elems.GetSize(); ++i) {
      for (UInt j = 0; j < dim; ++j) {
        if (elems[i]->extended->barycenter[j] > max[j])
          max[j] = elems[i]->extended->barycenter[j];
        if (elems[i]->extended->barycenter[j] < min[j])
          min[j] = elems[i]->extended->barycenter[j];
      }
    }

    // Computes lattice spacing
    StdVector<double> spacing; // the output
    GetElemShapeMap(elems[0], false)->GetEdgeLength(spacing);

    if (dim == 2)
      n[2] = 1;
    for (UInt i = 0; i < dim; ++i)
      n[i] = 1.00001 * (max[i] - min[i]) / spacing[i] + 1;

    LOG_DBG2(grid) << "GB(" << region << ") min=" << min.ToString() << " max=" << max.ToString() << " spacing=" << spacing.ToString() << " -> " << n.ToString();

    return n;
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
      numEdges = Elem::shapes[el->type].numEdges;
      for(UInt n=0; n<numEdges; n++) {
        UInt edgeNum = el->extended->edges[n] < 0 ? -el->extended->edges[n] : el->extended->edges[n];
        edgeCounts[edgeNum]++;
      }
    }

    for(UInt i=0; i<nElemsRegion; i++) {
      el = regionElems[i];
      // get number of edges
      numEdges = Elem::shapes[el->type].numEdges;
      for(UInt n=0; n<numEdges; n++) {
        UInt edgeNum = el->extended->edges[n] < 0 ? -el->extended->edges[n] : el->extended->edges[n];
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
        surfEl->type = Elem::ET_LINE2;
        surfEl->extended = new ExtendedElementInfo;

        switch(el->type) {
          case Elem::ET_TRIA6:
          case Elem::ET_QUAD8:
          case Elem::ET_QUAD9:
            surfEl->connect[2] = el->connect[n+numEdges];
            surfEl->type = Elem::ET_LINE2;
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

  void Grid::GetElemNums( std::set<UInt>& elemNums, 
                          std::set<UInt>& dims,
                          const StdVector<shared_ptr<EntityList> >& entities ) {
    elemNums.clear();
    // Loop over entities
    for( UInt i = 0; i < entities.GetSize(); ++i ) {
      UInt dim = GetEntityDim(entities[i]->GetName());
      dims.insert(dim);
      StdVector<UInt> elemNumVec;
      GetElemNumsByName(elemNumVec, entities[i]->GetName());
      elemNums.insert(elemNumVec.Begin(), elemNumVec.End());
    }
  }
  
  // =======================================================================
  //  ELEMENT / POINT MAPPING
  // =======================================================================

  
  void Grid::MapGlobPointsToLoc( const StdVector<PointElemMatch>& matches,
                                 StdVector<const Elem*>& elems,
                                 StdVector<LocPoint>& lps,
                                 Double tol,
                                 bool printWarnings,
                                 bool updatedGeo) {
    
    
    UInt numMatches = matches.GetSize(); 
    elems.Resize(numMatches);
    lps.Resize(numMatches);
    elems.Init( NULL );

    // loop over matches, perform global->local mapping of coordinates
    // and check, if coordinate is really contained in this element
#pragma omp parallel for num_threads(CFS_NUM_THREADS)
    for( Integer iM = 0; iM < (Integer) numMatches; ++iM ) {
      std::set<const Elem*>::const_iterator it;
      Vector<Double> locCoord;
      const std::set<const Elem*> & mElems = matches[iM].matches;
      StdVector<const Elem*> candidateElem, vagueCandElem;
      StdVector<LocPoint>  candidateLp, vagueCandLp;

      
      // loop over elements
      for( it = mElems.begin(); it != mElems.end(); ++it ) {

        // check, if global point can be mapped to the element
        shared_ptr<ElemShapeMap> esm = GetElemShapeMap(*it, updatedGeo);

        esm->Global2Local(locCoord, matches[iM].globCoord );
        if( esm->CoordIsInsideElem(locCoord, 0.0) ) {
          candidateElem.Push_back( *it );
          candidateLp.Push_back( locCoord );
        } else if ( esm->CoordIsInsideElem(locCoord, tol) ) {
          vagueCandElem.Push_back( *it );
          vagueCandLp.Push_back( locCoord );
        }
      }

      // Check, how many elements have been found
      if ( candidateElem.GetSize() == 0 ) {
        if ( vagueCandElem.GetSize() > 0 ) {
          elems[iM] = vagueCandElem[0];
          lps[iM] = vagueCandLp[0];
        } else if (printWarnings) {
          WARN( "No element found for location "
                << matches[iM].globCoord.ToString() );
        }
      } else {
        elems[iM] = candidateElem[0];
        lps[iM] = candidateLp[0];
      }
    }
  }
  
  void Grid::CreateBBoxFromElement(const Elem* elem,
                                   Double globToler,
                                   Double* bbox,
                                   double updated)
  {
    Vector<Double> p;
    Double& xmin = bbox[0];
    Double& xmax = bbox[3];
    Double& ymin = bbox[1];
    Double& ymax = bbox[4];
    Double& zmin = bbox[2];
    Double& zmax = bbox[5];

    // Create an exact bounding box from all corner nodes.
    GetNodeCoordinate(p, elem->connect[0],updated);
    UInt globalDim = p.GetSize();


    xmin = xmax = p[0];
    ymin = ymax = p[1];
    if(globalDim == 2) {
      zmin = zmax = 0.0;
    }
    else {
      zmin = zmax = p[2];
    }

    for(UInt j = 1, n=elem->connect.GetSize(); j < n; ++j)
    {
      GetNodeCoordinate(p, elem->connect[j],updated);
      xmin = (p[0] < xmin) ? p[0] : xmin;
      xmax = (p[0] > xmax) ? p[0] : xmax;
      ymin = (p[1] < ymin) ? p[1] : ymin;
      ymax = (p[1] > ymax) ? p[1] : ymax;
      if (p.GetSize() == 3) {  // TODO: Linienelemente NC interface probleme, coordinaten probleme // Fr[her statt p.GetSize(), globalDim
        zmin = (p[2] < zmin) ? p[2] : zmin;
        zmax = (p[2] > zmax) ? p[2] : zmax;
      }
    }

    Vector<Double> dia(3);
    dia[0] = xmax - xmin;
    dia[1] = ymax - ymin;
    dia[2] = zmax - zmin;

    // If a two-dimensional element is part of a three dimensional grid we use the maximum diameter of all diameters
    UInt elemDim = Elem::GetShape( Elem::GetShapeType( elem->type) ).dim;
    if (elemDim < globalDim) {
      Double maxDia = dia[0];
      UInt i = dia[1] > dia[2] ? 1 : 2;
      maxDia = dia[i] > maxDia ? dia[i] : maxDia;
      Double thisTol = globToler*maxDia;

      xmin -= thisTol;
      xmax += thisTol;
      ymin -= thisTol;
      ymax += thisTol;
      zmin -= thisTol;
      zmax += thisTol;
    } else {
      xmin -= globToler*dia[0];
      xmax += globToler*dia[0];
      ymin -= globToler*dia[1];
      ymax += globToler*dia[1];
      zmin -= globToler*dia[2];
      zmax += globToler*dia[2];
    }
  }

#ifdef USE_CGAL

  // ========================================================================
  //  C G A L  -  S P E C I F I C   I M P L E M E N T A T I O N
  // ========================================================================
  

  //! Define box handler, which additionally stores an index
  typedef CGAL::Box_intersection_d
      ::Box_with_handle_d<double,3,const UInt*> HandleBox;


  typedef CGAL::Bbox_3 BBox3D;
  
  // Iterator reporter class, returning the two ids of the CGAL-Boxes,
  // the first being the element number, the second being the node index
  template <class OutputIterator>
  struct CGAL_ElemPointIdReporter {
    OutputIterator it;
    CGAL_ElemPointIdReporter(OutputIterator i  ) 
    : it(i) {} // store iterator in object

    // We write the id-number of box a to the output iterator assuming
    // that box b (the query box) is not interesting in the result.
    void operator()( const HandleBox& a, const HandleBox& b) {
      UInt elemNum1 = *a.handle();
      UInt elemNum2 = *b.handle();
      std::pair<UInt, UInt > pair;
      pair.first = elemNum2;
      pair.second = elemNum1;
      *it++ = pair;
    }
  };
  // helper function to create the function object
  template <class Iter> 
  CGAL_ElemPointIdReporter<Iter> elemPointIdReporter(Iter it) 
  { return CGAL_ElemPointIdReporter<Iter>(it); }

  
  

  HandleBox Grid::CreateBoxFromCoord( const Vector<double>& coords, UInt* id,
                                      Double tol )
  {
    if(coords.GetSize()==2){
      return HandleBox(BBox3D(coords[0]-tol/2.0, coords[1]-tol/2.0, 0.0,
                              coords[0]+tol/2.0, coords[1]+tol/2.0, 0.0), id);
    }else{
      return HandleBox(BBox3D(coords[0]-tol/2.0, coords[1]-tol/2.0,
                              coords[2]-tol/2.0, coords[0]+tol/2.0,
                              coords[1]+tol/2.0, coords[2]+tol/2.0), id);
    }
  }

  void Grid::MapPointsToBoundingBoxes( StdVector<PointElemMatch>& matches,
                                       const StdVector<shared_ptr<EntityList> >& srcEntities,
                                       Double tol,
                                       bool updatedGeo ) {
    boost::array<Double,6> bbox;

    // If we haven't initialized the grid bounding boxes yet, do so now!
    if(elemBoxes_.empty() || updatedGeo == 1)
    {
      StdVector<Elem*> elems;
      Vector<Double> p(3);
      GetElems(elems, ALL_REGIONS);
      
      // Loop over dimensions
      for( UInt dim = 1; dim <= GetDim(); ++dim ) {
        std::vector<HandleBox> & boxes = elemBoxes_[dim];
        UInt size = this->GetNumElemOfDim(GetDim());

        boxes.reserve( size );
        for(UInt i = 0; i < elems.GetSize(); i++)       {
          // immediately leave, if the dimension of the element is 
          // lower-dimensional
          if( Elem::shapes[elems[i]->type].dim != dim ) 
            continue;

          CreateBBoxFromElement(elems[i], tol, &bbox[0], updatedGeo);

          HandleBox hbox(BBox3D(bbox[0], bbox[1], bbox[2],
                                bbox[3], bbox[4], bbox[5]),
                         &elems[i]->elemNum);

          boxes.push_back( hbox );
        }
      } //loop: dimension
    }
    
    // Get all element numbers and their dimension
    std::set<UInt> elemNums;
    std::set<UInt> dims;
    if( srcEntities.GetSize() != 0 ) {
      GetElemNums(elemNums, dims, srcEntities);
    } else {
      dims.insert(1);
      dims.insert(GetDim());
      dims.insert(GetDim()-1);
    }
    
    // now set up box list containing the point coordinates
    UInt numPoints = matches.GetSize();
    std::vector<HandleBox> pointBoxes (numPoints);
    
    // create also temporary index array (will be automatically deleted)
    boost::scoped_array<UInt> nodeIndices(new UInt[numPoints]);
    
    for( UInt i = 0; i < numPoints; ++i ) {
      nodeIndices[i] = i;
      pointBoxes[i] = CreateBoxFromCoord( matches[i].globCoord,
                                          &nodeIndices[i], tol );
    }
    
    // Loop over all dimensions
    std::set<UInt>::const_iterator dimIt = dims.begin();
    for( ; dimIt != dims.end(); ++dimIt ) {
      UInt dim = *dimIt;
      std::vector<HandleBox> & boxes = elemBoxes_[dim];

      // run the intersection algorithm and store results in a vector
      std::vector< std::pair<UInt, UInt > > result;
      CGAL::box_intersection_d( boxes.begin(), boxes.end(),
                                pointBoxes.begin(), pointBoxes.end(),
                                elemPointIdReporter( std::back_inserter( result )));

      // now loop over all results and store for each point all candidate elements
      // if they are contained in the desired regions
      std::vector< std::pair<UInt, UInt> >::iterator it = result.begin();
      for(; it != result.end(); ++it ) {
        UInt pointIndex = it->first;
        UInt elemNum = it->second;
        const Elem* ptEl = GetElem(elemNum);
        if( elemNums.size() ) { 
          if( elemNums.find(ptEl->elemNum) != elemNums.end() ) {
            matches[pointIndex].matches.insert(ptEl);
          }
        } else {
          matches[pointIndex].matches.insert(ptEl);
        }
      }
    }
  }

#ifdef USE_EIGEN
  void Grid::SetUpMongeForm(MongeForm& mongeForm, const UInt& degreePolynomFitting, const UInt& degreeMongeCoeffs, 
                        const std::vector<DPoint>& nodeCoordinates) {
    if (dim_ == 2)
      EXCEPTION("Grid::SetUpMongeForm() is only implemented in 3D!");
    // create Monge Form and fit
    MongeViaJetFitting mongeFit;
    mongeForm = mongeFit(nodeCoordinates.begin(), nodeCoordinates.end(), degreePolynomFitting, degreeMongeCoeffs);
  };

  void Grid::ConvertVectorToPoint_3Format(std::vector<DPoint>& pointsAsPoint_3Format, StdVector<Vector<Double>>& pointsAsCfsVectors) {
    UInt numPoints = pointsAsCfsVectors.GetSize();
    pointsAsPoint_3Format.resize(numPoints);
    if (numPoints > 0) {
      UInt dim = pointsAsCfsVectors[0].GetSize();
      if (dim == 1) {
        EXCEPTION("ConvertVectorToPoint_3Format() is untested for 1 dimension");
        for (UInt i = 0; i < numPoints; i++) {
          DPoint p(pointsAsCfsVectors[i][0],0,0);
          pointsAsPoint_3Format[i] = p;
        }
      }
      else if (dim == 2) {
        EXCEPTION("ConvertVectorToPoint_3Format() is untested for 2 dimensions");
        for (UInt i = 0; i < numPoints; i++) {
          DPoint p(pointsAsCfsVectors[i][0],pointsAsCfsVectors[i][1],0);
          pointsAsPoint_3Format[i] = p;
        }
      } 
      else if (dim == 3) {
        for (UInt i = 0; i < numPoints; i++) {
          DPoint p(pointsAsCfsVectors[i][0],pointsAsCfsVectors[i][1],pointsAsCfsVectors[i][2]);
          pointsAsPoint_3Format[i] = p;
        }
      } else
        EXCEPTION("Please provide 1, 2, or 3 dimensions in 'pointsAsCfsVectors'");
    } else
      EXCEPTION("No entries for conversion provided!");
  };

  void Grid::ConvertVectorFromPoint_3Format(StdVector<Vector<Double>>& pointsAsCfsVectors, std::vector<DPoint>& pointsAsPoint_3Format) {
    EXCEPTION("ConvertVectorFromPoint_3Format() is untested yet.");
    size_t numPoints = pointsAsPoint_3Format.size();
    pointsAsCfsVectors.Resize(numPoints);
    if (numPoints > 0) {
      for (UInt i = 0; i < numPoints; i++) {
        pointsAsCfsVectors[i].Resize(3);
        pointsAsCfsVectors[i][0] = pointsAsPoint_3Format[i].x();
        pointsAsCfsVectors[i][1] = pointsAsPoint_3Format[i].y();
        pointsAsCfsVectors[i][2] = pointsAsPoint_3Format[i].z();
      }
    } else
      EXCEPTION("No entries for conversion provided!");
  };
#endif // USE_EIGEN
#elif USE_LIBFBI // USE_CGAL
} // end namespace CoupledField


namespace fbi {
  template<>
  struct Traits<CoupledField::Elem*> : mpl::TraitsGenerator<double, double, double> {};
  
  template<>
  struct Traits< CoupledField::Vector<Double>* > : mpl::TraitsGenerator<double, double, double> {};
} // end namespace fbi

namespace CoupledField {

  struct ElemBoxGenerator
  {
    template <size_t N>
    typename fbi::tuple_element<N,
                                typename fbi::Traits<Elem*>::key_type>::type
    get(const Elem*) const;
    Grid* ptGrid_;
    double globTol_;
    UInt dim_;
    ElemBoxGenerator(Grid* ptGrid, double globTol, UInt dim)
      : ptGrid_(ptGrid), globTol_(globTol), dim_(dim) {}
  };
  
  template <>
  std::pair<double, double>
  ElemBoxGenerator::get<0>(const Elem* elem) const
  {
    // create bounding box array, with the following indices:
    // [xmin, ymin,  zmin, xmax, ymax, zmax]
    //   0     1      2     3     4     5
    boost::array<Double,6> bbox;
    
    ptGrid_->CreateBBoxFromElement(elem, globTol_, &bbox[0], false); // set the updatedGeometry flag to false (similar to the previous implementation)

    return std::make_pair(bbox[0], bbox[3]);
  }
  
  template <>
  std::pair<double, double>
  ElemBoxGenerator::get<1>(const Elem* elem) const
  {
    // create bounding box array, with the following indices:
    // [xmin, ymin,  zmin, xmax, ymax, zmax]
    //   0     1      2     3     4     5
    boost::array<Double,6> bbox;
    
    ptGrid_->CreateBBoxFromElement(elem, globTol_, &bbox[0], false); // set the updatedGeometry flag to false (similar to the previous implementation)

    return std::make_pair(bbox[1], bbox[4]);
  }

  template <>
  std::pair<double, double>
  ElemBoxGenerator::get<2>(const Elem* elem) const
  {
    // create bounding box array, with the following indices:
    // [xmin, ymin,  zmin, xmax, ymax, zmax]
    //   0     1      2     3     4     5
    boost::array<Double,6> bbox;
    
    ptGrid_->CreateBBoxFromElement(elem, globTol_, &bbox[0], false);

    return std::make_pair(bbox[2], bbox[5]);
  }

  struct PointBoxGenerator
  {
    template <size_t N>
    typename fbi::tuple_element<N,
                                typename fbi::Traits< Vector<Double>* >::key_type>::type
    get(const Vector<Double>*) const;
    UInt dim_;
    PointBoxGenerator(UInt dim)
      : dim_(dim) {}
  };
  
  template <>
  std::pair<double, double>
  PointBoxGenerator::get<0>(const Vector<Double>* p) const
  {
    return std::make_pair((*p)[0], (*p)[0]);
  }
  
  template <>
  std::pair<double, double>
  PointBoxGenerator::get<1>(const Vector<Double>* p) const
  {
    return std::make_pair((*p)[1], (*p)[1]);
  }

  template <>
  std::pair<double, double>
  PointBoxGenerator::get<2>(const Vector<Double>* p) const
  {
    if(dim_ == 3) 
    {
      return std::make_pair((*p)[2], (*p)[2]);
    }
    else 
    {
      return std::make_pair(0, 0);
    }    
  }

  void Grid::MapPointsToBoundingBoxes( StdVector<PointElemMatch>& matches,
                                       const StdVector<shared_ptr<EntityList> >& srcEntities,
                                       Double tol,
                                       bool updatedGeo ) {
    WARN("Updated geometry is not used, please implement me!");
    
    std::vector< Elem* > elems;
    std::vector< Vector<Double>* > points;
    
    StdVector<Elem*> allElems;
    GetElems(allElems, ALL_REGIONS);

    if(srcEntities.GetSize()) 
    {
      for( UInt i = 0; i < srcEntities.GetSize(); ++i ) {
        EntityIterator it = srcEntities[i]->GetIterator();
        for( ; !it.IsEnd(); it++ ) {
          const Elem* ptEl = it.GetElem();
          StdVector<Elem*>::iterator elIt;
          elIt = std::find(allElems.Begin(), allElems.End(), ptEl);          
          elems.push_back(*elIt);
        }
      }
    }
    else 
    {
      StdVector<Elem*> volElems;
      GetVolElems(volElems, ALL_REGIONS);
      std::copy(volElems.Begin(), volElems.End(), std::back_inserter(elems));
    } 
    
    UInt numPts = matches.GetSize();
    // Loop over all matches
    for( UInt iPt = 0; iPt < numPts; ++iPt ) {
      Vector<Double>* point = &matches[iPt].globCoord;
      
      points.push_back(point);
    } // loop over points
    
    ElemBoxGenerator ebg(this, tol, GetDim());
    PointBoxGenerator pbg(GetDim());

    // For 2D:
    //  auto adjList = fbi::SetA<Elem*, 0, 1>::SetB<Vector<Double>*, 0, 1>::intersect(
    //    elems, ElemBoxGenerator(this, tol, GetDim()), points, PointBoxGenerator(GetDim()));
    
    fbi::SetA<Elem*, 0, 1, 2>::ResultType adjList;
    adjList = fbi::SetA<Elem*, 0, 1, 2>::SetB<Vector<Double>*, 0, 1, 2>::intersect(
      elems, ebg, points, pbg);
    
    typedef fbi::SetA<Elem*, 0, 1, 2>::IntType LabelType;
    
    for(UInt i=elems.size(), n=adjList.size(), ptIdx = 0;
        i < n;
        i++, ptIdx++ ) {
      
      std::vector<LabelType> queryResultIndexes = adjList[i];
      std::vector<LabelType>::iterator it = queryResultIndexes.begin();
      std::vector<LabelType>::iterator end = queryResultIndexes.end();
      
      // std::cout << "Size of queryResultIndexes for node: " << queryResultIndexes.size() << std::endl;
      for( ; it != end; it++) 
      {
        const Elem* ptEl = elems[(*it)];
        // std::cout << "Elem: " << ptEl->elemNum << std::endl;
        matches[ptIdx].matches.insert(ptEl);
      }    
    }
  }
  
#else // USE_CGAL
  
  // This is a very basic implementation for axis-parallel box intersection. It is just used
  // as internal replacement in case we want to use valgrind and can not use CGAL.
  void Grid::MapPointsToBoundingBoxes( StdVector<PointElemMatch>& matches,
                                       const StdVector<shared_ptr<EntityList> >& srcEntities,
                                       Double tol,
                                       bool updatedGeo) {

    // obtain all volume elements from grid
    StdVector<Elem*> elems;
    GetElems(elems, ALL_REGIONS);
    
    // check, if element boxes are already initialized or if we have to reinitialize them due to the updated geometry
    if( elemBoxes_.size() == 0 || updatedGeo == 1) {

      Vector<Double> p(3);
      
      // Loop over dimensions
      for( UInt dim = 1; dim <= GetDim(); ++dim ) {
        StdVector<BoxType> & boxes = elemBoxes_[dim];
        boxes.Clear(); // for the case of two updated PDEs we need to clear the vector since continuously appending will lead to wrong bounding boxes
        boxes.Reserve(GetNumElemOfDim(dim));

        // loop over all elements
        for(UInt i = 0, m=elems.GetSize(); i < m; i++)       {
          Vector<Double> p(3);

          // only map elements of the given dimension
          if( Elem::shapes[elems[i]->type].dim != dim) 
            continue;

          // create bounding box array, with the following indices:
          // [xmin, ymin,  zmin, xmax, ymax, zmax]
          //   0     1      2     3     4     5
          boost::array<Double,6> bbox;

          CreateBBoxFromElement(elems[i], tol, &bbox[0], updatedGeo);
          
          // assemble tuple of (bounding box, element number)
          boxes.Push_back(BoxType(bbox, elems[i]->elemNum));
          
//          std::cerr << "created box for elem #" << elems[i]->elemNum << " with \n\t"
//              << bbox[0] << ", "
//              << bbox[1] << ", "
//              << bbox[2] << ", "
//              << bbox[3] << ", "
//              << bbox[4] << ", "
//              << bbox[5] << std::endl;
        } //loop: elements
      } //loop: dimensions
    } // check for empty element boxes


    // result vector (pair first: point-index, second: element number)
    std::vector<std::pair<UInt, UInt> > result;
    
    
    // Get all element numbers and their dimension
    std::set<UInt> elemNums;
    std::set<UInt> dims;
    if( srcEntities.GetSize() != 0 ) {
      GetElemNums(elemNums, dims, srcEntities);
    } else {
      dims.insert(1);
      dims.insert(GetDim());
      dims.insert(GetDim()-1);
    }
    
    UInt numPts = matches.GetSize();
    result.reserve(numPts * 2); // assume 2 matches on average
          
    // Loop over all dimensions
    std::set<UInt>::const_iterator dimIt = dims.begin();
    for( ; dimIt != dims.end(); ++dimIt ) {
      UInt dim = *dimIt;
      StdVector<BoxType> & boxes = elemBoxes_[dim];
      
      // Loop over all elements of given dimension
      UInt numElems = boxes.GetSize();

      for( UInt iEl = 0; iEl < numElems; ++iEl ) {
        const BoxType& actBox = boxes[iEl];
        
        const boost::array<Double,6> & eb = actBox.first;
        UInt elemNum = actBox.second;
//        const Elem* ptEl = GetElem(elemNum);
//        std::cerr << "checking elem #" << ptEl->elemNum
//            << " with bbox \n\t "
//            << "x: [" << eb[0] << ", " << eb[3] << "] "
//            << "y: [" << eb[1] << ", " << eb[4] << "] "
//            << "y: [" << eb[2] << ", " << eb[5] << "] "
//            << "\n----------------------------\n";

        // Loop over all points and check bounding box
        for( UInt iPt = 0; iPt < numPts; ++iPt ) {
          const Vector<Double> & point = matches[iPt].globCoord;
          if( point[0] >= eb[0] && point[0] <= eb[3]) {
            if( point[1] >= eb[1] && point[1] <= eb[4]) {
              if( dim == 3 ) {
                if( point[2] >= eb[2] && point[2] <= eb[5]) {
                  result.push_back(std::pair<UInt,UInt>(iPt, elemNum));
                }
              } else {
                result.push_back(std::pair<UInt,UInt>(iPt, elemNum));
              }
            }
          }
        } // loop over points
      } // loop over elements
    } // loop over dimension

    // now loop over all results and store for each point all candidate elements
    // if they are contained in the desired regions
    std::vector< std::pair<UInt, UInt> >::iterator it = result.begin();
    for(; it != result.end(); ++it ) {
      UInt pointIndex = it->first;
      UInt elemNum = it->second;
      const Elem* ptEl = GetElem(elemNum);
      if( elemNums.size() ) { 
        if( elemNums.find(ptEl->elemNum) != elemNums.end() ) {
          matches[pointIndex].matches.insert(ptEl);
        } else {
          //
        }
      } else {
        matches[pointIndex].matches.insert(ptEl);
      }
    }
  }
  

#endif // USE_CGAL
} // end of namespace
