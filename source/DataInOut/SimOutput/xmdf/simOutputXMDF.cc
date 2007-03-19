#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <algorithm>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/exception.hpp>
#include <boost/algorithm/string/trim.hpp>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/exception.hh"
#include "simInputXMDF.hh"

namespace fs = boost::filesystem;

namespace CoupledField {

 #define CHECK_HDF5_ERROR(HID_T, STR) { if(HID_T < 0) {                 \
        std::ostringstream ostr;                                        \
        ostr << STR;                                                    \
        Exception ex(NULL, __FILE__, __LINE__, ostr.str().c_str()); } }


  SimInputXMDF::SimInputXMDF(std::string fileName, ParamNode * inputNode) :
      SimInput(fileName, inputNode)
  {
    ENTER_FCN( "SimInputXMDF::XMDF" );
    mi_ = NULL;
    statsRead_ = false;
    fileName_ = fileName;
    genRegionNodes_ = false;
    
    // Change defaults according to XML file
    if(myParam_->Get("generateRegionNodes", false)->AsString() == "yes")
    {
        genRegionNodes_ = true;
    }

    // Do not print HDF5 exceptions by default
    H5::Exception::dontPrint();

  }


  SimInputXMDF::~SimInputXMDF() {
    ENTER_FCN( "SimInputXMDF::~XMDF" );
    xfCloseFile(fileId_);
  }

  void SimInputXMDF::InitModule(Grid *mi)
  {
    try 
    {
      fs::path fn = fs::system_complete(fileName_);
      fn.normalize();
      baseDir_ = fn.branch_path().native_directory_string();
      baseName_ = (fs::change_extension( fn.leaf(), "" )).native_directory_string();
      if(fs::extension(fn) == "")
      {
        fn = fs::change_extension( fn, ".h5" );
      }
      fileName_ = fn.native_directory_string();
    } catch (fs::filesystem_error& ex)
    {
      MESHIO_ERROR("Received exception: " << ex.what());
      return;
    }
    
    MESHIO_DEBUG("fileName_: " << fileName_);
    MESHIO_DEBUG("baseDir_: " << baseDir_);
    MESHIO_DEBUG("baseName_: " << baseName_);
    
    statsRead_ = false;
    
    mi_ = mi;
    multiStep_ = 0;
    step_ = 0;
    msChange_ = true;

    xid  xGroupId = -1;
    Integer    status;
    std::string pathToMesh = "Mesh";
    std::stringstream strBuf;

    status = xfOpenFile(fileName_.c_str(), &fileId_, true);
    strBuf << "Could not open XMDF file " << fileName_;
    CHECK_HDF5_ERROR(status, strBuf.str());
    
    status = xfOpenGroup(fileId_, pathToMesh.c_str(), &xGroupId);
    strBuf << "Could not open mesh group '" << pathToMesh << "'";
    CHECK_HDF5_ERROR(status, strBuf.str());   

    try
    {
      H5::Group group(xGroupId);
      mainRoot_ = group.openGroup("/");
    }
    catch (H5::Exception& h5ex)
    {
        EXCEPTION(h5ex.getCDetailMsg());
    }
  }

  void SimInputXMDF::ReadMesh()
  {
    Integer status;
    std::stringstream strBuf;
    H5::Group mGroup;
    xid meshGroupId;

    try
    {
      mGroup = mainRoot_.openGroup("Mesh");
    }
    catch (H5::Exception& h5ex)
    {
        EXCEPTION(h5ex.getCDetailMsg());
    }

    meshGroupId = mGroup.getLocId();

    // Get the number of elements, nodes, and Maximum number of nodes per element
    status = xfGetNumberOfElements(meshGroupId, (int*) &fg_nElems);
    if (status >= 0) {
      status = xfGetNumberOfNodes(meshGroupId, (int*) &fg_nNodes);
      if (status >= 0) {
        status = xfGetMaxNodesInElem(meshGroupId, (int*) &fg_nNodesPerElem);
      }
    }
    CHECK_HDF5_ERROR(status, "Unable to read number of elements, nodes or max elemnodes.");

    fg_ElemTypes.resize(fg_nElems);
    MESHIO_DEBUG("fg_nElems: " << fg_nElems);

    status = xfReadElemTypes(meshGroupId, fg_nElems, (int*) &fg_ElemTypes[0]);
    CHECK_HDF5_ERROR(status, "Unable to read element types.");

    // Nodes in each element
    fg_NodesInElem.resize(fg_nElems*fg_nNodesPerElem);
    status = xfReadElemNodeIds(meshGroupId,
                               fg_nElems,
                               fg_nNodesPerElem,
                               (int*) &fg_NodesInElem[0]);
    CHECK_HDF5_ERROR(status, "Unable to read maximum number of nodes per elem.");

    // NodeLocations
    fg_XNodeLocs.resize(fg_nNodes);
    fg_YNodeLocs.resize(fg_nNodes);
    fg_ZNodeLocs.resize(fg_nNodes);

    MESHIO_DEBUG("fg_nNodes: " << fg_nNodes << " " << fg_XNodeLocs.size());

    status = xfReadXNodeLocations(meshGroupId, fg_nNodes, &fg_XNodeLocs[0]);
    if (status >= 0) {
      status = xfReadYNodeLocations(meshGroupId, fg_nNodes, &fg_YNodeLocs[0]);
      if (status >= 0) {
        status = xfReadZNodeLocations(meshGroupId, fg_nNodes, &fg_ZNodeLocs[0]);
      }
    }
    CHECK_HDF5_ERROR(status, "Unable to read nodes.");

    mi_->AddNodes(fg_nNodes);
    Integer node;

    for (node=0; node < fg_nNodes; node++) {
        Point p;

        p[0] = fg_XNodeLocs[node];
        p[1] = fg_YNodeLocs[node];
        p[2] = fg_ZNodeLocs[node];
    
        mi_->SetNodeCoordinate(node+1, p);
    }

    // Release temp memory of nodes
    fg_XNodeLocs.clear();
    fg_YNodeLocs.clear();
    fg_ZNodeLocs.clear();

    FEType eType;
    UInt connectIdx = 0;
    Integer elem;

    mi_->AddElems(fg_nElems);
    
    connectIdx = 0;
    for (elem=0; elem < fg_nElems; elem++) {
      eType = XMDFElemType2ElemType(fg_ElemTypes[elem]);
      if((eType == ET_LINE3) ||
         (eType == ET_TRIA6) ||
         (eType == ET_QUAD8) ||
         (eType == ET_QUAD9))
        this->ReorderConnectivity(eType,
                                  false,
                                  (UInt*)&fg_NodesInElem[connectIdx],
                                  (UInt*)&fg_NodesInElem[connectIdx]);

      mi_->SetElemData(elem+1, eType, 0, (UInt*)&fg_NodesInElem[connectIdx]);
      connectIdx += fg_nNodesPerElem;
    }

    try 
    {
      ReadRegions(mGroup);
      ReadNamedNodes(mGroup);
      ReadNamedElems(mGroup);
    }
    catch (H5::Exception& h5ex)
    {
        EXCEPTION(h5ex.getCDetailMsg());
    }
    catch (Exception& ex)
    {
        EXCEPTION(ex.GetMsg());
    }

    // Release temp memory of elems
    fg_ElemTypes.clear();
    fg_NodesInElem.clear();
  }


  // ======================================================
  // GENERAL MESH INFORMATION
  // ======================================================
  UInt SimInputXMDF::GetDim() {
    MESHIO_WARN("SimInputXMDF::ReadMesh() not implemented");
    return 0;
  }
  
  UInt SimInputXMDF::GetNumNodes(){
    MESHIO_WARN("SimInputXMDF::ReadMesh() not implemented");
    return 0;
  }
    
  UInt SimInputXMDF::GetNumElems(const Integer dim){
    MESHIO_WARN("SimInputXMDF::ReadMesh() not implemented");
    return 0;
  }
  
  UInt SimInputXMDF::GetNumRegions(){
    MESHIO_WARN("SimInputXMDF::ReadMesh() not implemented");
    return 0;
  }

  UInt SimInputXMDF::GetNumNamedNodes(){
    MESHIO_WARN("SimInputXMDF::ReadMesh() not implemented");
    return 0;
  }

  UInt SimInputXMDF::GetNumNamedElems(){
    MESHIO_WARN("SimInputXMDF::ReadMesh() not implemented");
    return 0;
  }
  
  // ======================================================
  // ENTITY NAME ACCESS
  // ======================================================

  void SimInputXMDF::GetAllRegionNames( std::vector<std::string> & regionNames ){
    MESHIO_WARN("SimInputXMDF::ReadMesh() not implemented");
  }

  void SimInputXMDF::GetRegionNamesOfDim( StdVector<std::string> & regionNames,
                                   const UInt dim )
  {
    Error("SimInputXMDF::GetRegionNamesOfDim() not implemented",
          __FILE__,__LINE__);
  }
    
  void SimInputXMDF::GetNodeNames( StdVector<std::string> & nodeNames )
  {
    Error("SimInputXMDF::GetNodeNames() not implemented",
          __FILE__,__LINE__);
  }
  
  void SimInputXMDF::GetElemNames( StdVector<std::string> & elemNames )
  {
    Error("SimInputXMDF::GetElemNames() not implemented",
          __FILE__,__LINE__);
  }

  // =========================================================================
  // MISCELLANEOUS METHODS
  // =========================================================================


  // ***************
  //   Type2ptElem
  // ***************

  FEType SimInputXMDF::XMDFElemType2ElemType( const Integer type ) {

    switch (type) {
    case 100:   return ET_LINE2;   // 1D linear
    case 101:   return ET_LINE3;   // 1D quadratic
    case 200:   return ET_TRIA3;   // 2D linear triangle
    case 201:   return ET_TRIA6;   // 2D quadratic triangle
    case 210:   return ET_QUAD4;   // 2D linear quadrilateral
    case 211:   return ET_QUAD8;   // 2D quadratic quadrilateral
    case 212:   return ET_QUAD9;   // 2D quadratic quadrilateral with center node
    case 300:   return ET_TET4;    // 3D linear tetrahedron
    case 30010: return ET_TET10;   // 3D quadratic tetrahedron
    case 310:   return ET_WEDGE6;  // 3D linear prism
    case 31010: return ET_WEDGE15; // 3D quadratic prism
    case 320:   return ET_HEXA8;   // 3D linear hexahedron
    case 32010: return ET_HEXA20;  // 3D quadratic hexahedron
    case 32011: return ET_HEXA27;  // 3D quadratic hexahedron with center nodes
    case 330:   return ET_PYRA5;   // 3D linear pyramid
    case 33010: return ET_PYRA13;  // 3D quadratic pyramid
    }

    return ET_UNDEF;
    // This place should never be reached!
  }

  Integer SimInputXMDF::ElemType2XMDFElemType( const FEType type )
  {
    switch (type) {
    case ET_LINE2:   return 100;   // 1D linear
    case ET_LINE3:   return 101;   // 1D quadratic
    case ET_TRIA3:   return 200;   // 2D linear triangle
    case ET_TRIA6:   return 201;   // 2D quadratic triangle
    case ET_QUAD4:   return 210;   // 2D linear quadrilateral
    case ET_QUAD8:   return 211;   // 2D quadratic quadrilateral
    case ET_QUAD9:   return 212;   // 2D quadratic quadrilateral with center node
    case ET_TET4:    return 300;   // 3D linear tetrahedron
    case ET_TET10:   return 30010; // 3D linear tetrahedron
    case ET_WEDGE6:  return 310;   // 3D linear prism
    case ET_WEDGE15: return 31010; // 3D linear prism
    case ET_HEXA8:   return 320;   // 3D linear hexahedron
    case ET_HEXA20:  return 32010; // 3D linear hexahedron
    case ET_HEXA27:  return 32011; // 3D linear hexahedron
    case ET_PYRA5:   return 330;   // 3D linear pyramid
    case ET_PYRA13:  return 33010; // 3D linear pyramid
    default: return ET_UNDEF;
    }

    // This place should never be reached!
    return -1;
  }


  void SimInputXMDF::ReorderConnectivity( const Integer eType,
                                  const bool toXMDF,
                                  const UInt* in,
                                  UInt* out) {
    static Integer toXMDFIdxsLine[] = {0, 2, 1};
    static Integer fromXMDFIdxsLine[] = {0, 2, 1};
    static Integer toXMDFIdxsTria[] = {0, 3, 1, 4, 2, 5};
    static Integer fromXMDFIdxsTria[] = {0, 2, 4, 1, 3, 5};
    static Integer toXMDFIdxsQuad[] = {0, 4, 1, 5, 2, 6, 3, 7, 8};
    static Integer fromXMDFIdxsQuad[] = {0, 2, 4, 6, 1, 3, 5, 7, 8};

    std::vector<UInt> tmp;
    Integer* it;

    UInt numNodes = NUM_ELEM_NODES[eType];
    tmp.resize(numNodes);
    memcpy(&tmp[0], in, numNodes*sizeof(UInt));

    switch (eType) {
    case ET_LINE3: // 1D quadratic line
      if(toXMDF)
        it = toXMDFIdxsLine;
      else
        it = fromXMDFIdxsLine;
      break;
    case ET_TRIA6: // 2D quadratic triangle
      if(toXMDF)
        it = toXMDFIdxsTria;
      else
        it = fromXMDFIdxsTria;
      break;
    case ET_QUAD8: // 2D quadratic quadrilateral
    case ET_QUAD9: // 2D quadratic quadrilateral with center node
      if(toXMDF)
        it = toXMDFIdxsQuad;
      else
        it = fromXMDFIdxsQuad;
      break;
    }

    for(UInt i=0; i<numNodes; i++, it++)
    {
      out[i] = tmp[*it];
    }
  }

  void SimInputXMDF::ReadRegions(const H5::Group& meshGroup)
  {
    hid_t status;
    H5::Group regionGroup;
    H5::Group regionElemGroup;
    H5::Group regionNodeGroup;
    std::vector< UInt > regionElems;
    std::vector< UInt > regionNodes;
    std::vector< RegionIdType > regionIds;
    std::set<UInt,
          std::less<UInt>,
          std::allocator<UInt> > regionNodeSet;
    Integer fg_nElems, numInt, numElemNodes, elemIdx;
    StdVector<Elem*> elems;
    std::stringstream nodeNamesStr;

    status = xfGetNumberOfElements(meshGroup.getLocId(), (int*) &fg_nElems);
    CHECK_HDF5_ERROR(status, "Could not read number of elements.");

    if(mi_->GetNumElems() != fg_nElems)
    {
        EXCEPTION ("Number of elements in mesh does not "    \
                   "match numelems in file.");
    }

    if(!statsRead_)
      ReadMeshStats(meshGroup);

    if(numRegions_ == 0)
    {
      mi_->AddRegions(regionNames_, regionIds);
      numRegions_ = 1;
      return;
    }

    try 
    {
      regionGroup = meshGroup.openGroup("Regions");
      regionElemGroup = regionGroup.openGroup("Elements");
      regionNodeGroup = regionGroup.openGroup("Nodes");
    } catch (H5::Exception& h5ex)
    {
      EXCEPTION(h5ex.getCDetailMsg());
    }
    
    // Create the regions
    mi_->AddRegions(regionNames_, regionIds);

    for(Integer i=0; i<numRegions_; i++)
    {
      // Get number of elements in region
      status = xfGetPropertyNumber(regionElemGroup.getLocId(),
                                   regionNames_[i].c_str(),
                                   (int*) &numInt);
      CHECK_HDF5_ERROR(status, "Could not read number of elements in region.");
      
      // Read elements in region
      regionElems.resize(numInt);
      status = xfReadPropertyInt(regionElemGroup.getLocId(),
                                 regionNames_[i].c_str(),
                                 numInt,
                                 (int*) &regionElems[0]);
      CHECK_HDF5_ERROR(status, "Could not read elements in region.");
      

      if(genRegionNodes_)
          regionNodeSet.clear();

      mi_->GetElems(elems, ALL_REGIONS);
      for(int j=0; j<numInt; j++)
      {
          elemIdx = regionElems[j]-1;
          elems[elemIdx]->regionId = regionIds[i];

          if(genRegionNodes_)
          {
              numElemNodes = NUM_ELEM_NODES[XMDFElemType2ElemType(fg_ElemTypes[elemIdx])];
              //          UInt* ptElemNodes = (UInt*)&fg_NodesInElem[elemIdx*fg_nNodesPerElem];
              UInt* ptElemNodes = &elems[elemIdx]->connect[0];
              regionNodeSet.insert( ptElemNodes,
                                    ptElemNodes + numElemNodes);
          }          
      }

      if(genRegionNodes_)
      {
          regionNodes.resize(regionNodeSet.size());
      
          std::copy(regionNodeSet.begin(),
                    regionNodeSet.end(),
                    regionNodes.begin());
          
          nodeNamesStr.str("");
          nodeNamesStr << regionNames_[i] << "_nodes";
          
          mi_->AddNamedNodes(nodeNamesStr.str(), regionNodes);
      }
    }

    regionElemGroup.close();
    regionGroup.close();
  }

  void SimInputXMDF::ReadNamedNodes(const H5::Group& meshGroup)
  {
    H5::Group namedNodesGroup, nNodesGroup;
    hid_t status;
    std::vector< UInt > namedNodes;
    Integer numInt;

    if(!statsRead_)
      ReadMeshStats(meshGroup);

    if(!nodeNames_.size())
      return;
    
    try 
    {
      namedNodesGroup = meshGroup.openGroup("Named Nodes");
      nNodesGroup = namedNodesGroup.openGroup("Nodes");
    } catch (H5::Exception& h5ex)
    {
      EXCEPTION(h5ex.getCDetailMsg ());
    }
    
    for(UInt i=0; i<nodeNames_.size(); i++)
    {
      
      status = xfGetPropertyNumber(nNodesGroup.getLocId(),
                                   nodeNames_[i].c_str(),
                                   (int*) &numInt);
      if(status < 0)
      {
        nNodesGroup.close();
        namedNodesGroup.close();
      }
      CHECK_HDF5_ERROR(status, "Could not read number of named nodes.");

      namedNodes.resize(numInt);

      status = xfReadPropertyInt(nNodesGroup.getLocId(),
                                 nodeNames_[i].c_str(),
                                 numInt,
                                 (int*) &namedNodes[0]);
      if(status < 0)
      {
        nNodesGroup.close();
        namedNodesGroup.close();
      }
      CHECK_HDF5_ERROR(status, "Could not read named nodes.");

      mi_->AddNamedNodes(nodeNames_[i], namedNodes);
    }

    nNodesGroup.close();
    namedNodesGroup.close();
  }

  void SimInputXMDF::ReadNamedElems(const H5::Group& meshGroup)
  {
    H5::Group namedElemsGroup, nElemsGroup;
    hid_t status;
    std::vector< UInt > namedElems;
    Integer numInt;

    if(!statsRead_)
      ReadMeshStats(meshGroup);

    if(!elemNames_.size())
      return;
    
    try 
    {
      namedElemsGroup = meshGroup.openGroup("Named Elements");
      nElemsGroup = namedElemsGroup.openGroup("Elements");
    } catch (H5::Exception& h5ex)
    {
      EXCEPTION(h5ex.getCDetailMsg ());
    }
    
    for(UInt i=0; i<elemNames_.size(); i++)
    {
      status = xfGetPropertyNumber(nElemsGroup.getLocId(),
                                   elemNames_[i].c_str(),
                                   (int*) &numInt);
      if(status < 0)
      {
        nElemsGroup.close();
        namedElemsGroup.close();
      }
      CHECK_HDF5_ERROR(status, "Could not read number of named elems.");

      namedElems.resize(numInt);

      status = xfReadPropertyInt(nElemsGroup.getLocId(),
                                 elemNames_[i].c_str(),
                                 numInt,
                                 (int*) &namedElems[0]);
      if(status < 0)
      {
        nElemsGroup.close();
        namedElemsGroup.close();
      }
      CHECK_HDF5_ERROR(status, "Could not read named elems.");

      mi_->AddNamedElems(elemNames_[i], namedElems);
    }

    nElemsGroup.close();
    namedElemsGroup.close();
  }

  void SimInputXMDF::ReadMeshStats(const H5::Group& meshGroup)
  {
    hsize_t number;
    bool regionGroupExists = true;
    H5::Group regionGroup;
    std::vector<region_desc_type> region_desc;
    
    regionNames_.clear();
    regionDims_.clear();

    try 
    {
      regionGroup = meshGroup.openGroup("Regions");
    } catch (H5::Exception& h5ex)
    {
      MESHIO_WARN("Could not open region group.");
      regionNames_.push_back("XMDF_Region");
      regionDims_.push_back(3);
      numRegions_ = 0;
      regionGroupExists = false;
    }

    if(regionGroupExists)
    {
      try 
      {
        H5::DataSet dataset = regionGroup.openDataSet("Description");
    
        // Create a memory datatype for region_desc_type
        H5::CompType comptype( sizeof(region_desc_type) );

        region_desc.resize(1);
      
        comptype.insertMember( "Name", HOFFSET(region_desc_type, name),
                               H5::StrType(H5::PredType::C_S1,
                                           sizeof(region_desc[0].name)));
        //                             dataset.getCompType().getMemberStrType(0));
        comptype.insertMember( "Dimension", HOFFSET(region_desc_type, dim),
                               H5::PredType::NATIVE_UINT32);

        number = dataset.getStorageSize() / sizeof(region_desc_type);
    
        region_desc.resize(number);
        dataset.read(&region_desc[0], comptype);

      } catch (H5::Exception& h5ex)
      {
        EXCEPTION(h5ex.getCDetailMsg());
      }

      // Push the region names into a vector
      numRegions_ = number;
      regionNames_.resize(numRegions_);
      regionDims_.resize(numRegions_);
      for(int32_t i=0; i<numRegions_; i++)
      {
        region_desc[i].name[31] = 0;
        regionNames_[i] = region_desc[i].name;
        boost::algorithm::trim(regionNames_[i]);
        regionDims_[i] = region_desc[i].dim;
      
        MESHIO_INFO(regionNames_[i] << " " << regionDims_[i]);
      }
    }
    
    // Read Named Nodes Description
    std::vector< named_entity_desc_type > namedNodesDesc;
    bool namedNodesExist = true;
    nodeNames_.clear();
    
    H5::Group namedNodesGroup;
    try 
    {
      namedNodesGroup = meshGroup.openGroup("Named Nodes");
    } catch (H5::Exception& h5ex)
    {
      namedNodesExist = false;
      MESHIO_WARN("Could not open named nodes group.");
    }

    if(namedNodesExist)
    {
      try 
      {
        H5::DataSet dataset = namedNodesGroup.openDataSet("Description");
    
        // Create a memory datatype for region_desc_type
        H5::CompType comptype( sizeof(named_entity_desc_type) );

        comptype.insertMember( "Name", HOFFSET(named_entity_desc_type, name),
                               dataset.getCompType().getMemberStrType(0));

        number = dataset.getStorageSize() / sizeof(named_entity_desc_type);

        namedNodesDesc.resize(number);
        dataset.read(&namedNodesDesc[0], comptype);
      } catch (H5::Exception& h5ex)
      {
        EXCEPTION(h5ex.getCDetailMsg());
      }

      // Push the node names into a vector
      nodeNames_.resize(number);
      for(int32_t i=0; i<number; i++)
      {
        namedNodesDesc[i].name[31] = 0;
        nodeNames_[i] = namedNodesDesc[i].name;
        boost::algorithm::trim(nodeNames_[i]);
      
        MESHIO_INFO("Named Nodes: " << nodeNames_[i]);
      }
    }
    
    // Read Named Elems Description
    std::vector< named_entity_desc_type > namedElemsDesc;
    bool namedElemsExist = true;
    elemNames_.clear();

    H5::Group namedElemsGroup;
    try 
    {
      namedElemsGroup = meshGroup.openGroup("Named Elements");
    } catch (H5::Exception& h5ex)
    {
      namedElemsExist = false;
      MESHIO_WARN("Could not open named elems group.");
    }

    if(namedElemsExist)
    {
      try 
      {
        H5::DataSet dataset = namedElemsGroup.openDataSet("Description");
    
        // Create a memory datatype for region_desc_type
        H5::CompType comptype( sizeof(named_entity_desc_type) );

        comptype.insertMember( "Name", HOFFSET(named_entity_desc_type, name),
                               dataset.getCompType().getMemberStrType(0));

        number = dataset.getStorageSize() / sizeof(named_entity_desc_type);

        namedElemsDesc.resize(number);
        dataset.read(&namedElemsDesc[0], comptype);
      } catch (H5::Exception& h5ex)
      {
        EXCEPTION(h5ex.getCDetailMsg());
      }

      // Push the elem names into a vector
      elemNames_.resize(number);
      for(int32_t i=0; i<number; i++)
      {
        namedElemsDesc[i].name[31] = 0;
        elemNames_[i] = namedElemsDesc[i].name;
        boost::algorithm::trim(elemNames_[i]);
      
        MESHIO_INFO("Named Elems: " << elemNames_[i]);
      }
    }
    
    statsRead_ = true;
  }

}

/// Local Variables:
/// mode: C++
/// c-basic-offset: 2
/// End:
