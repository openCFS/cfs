#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <algorithm>

#include "mod_xmdf.hh"

namespace CoupledField {

#define CHECK_HDF5_ERROR(HID_T, STR) { if(HID_T < 0) { Exception ex(STR, __FILE__, __LINE__); errHandler->Error(ex); } }


  XMDF::XMDF(std::string fileName) : MeshIOModule(fileName) {
    ENTER_FCN( "XMDF::XMDF" );
    mi_ = NULL;
  }


  XMDF::~XMDF() {
    ENTER_FCN( "XMDF::~XMDF" );
  }


  void XMDF::RegisterModule()
  {
  }
    
  
  void XMDF::InitModule(const std::string & baseDir,
                        const std::string & baseName,
                        Grid *mi)
  {
    baseDir_ = baseDir;
    baseName_ = baseName;
    mi_ = mi;
  }

  void XMDF::ReadMesh()
  {
    std::stringstream strBuffer;
    std::string fn;
    std::string pathToMesh = "Mesh";
    strBuffer << this->baseDir_ << baseName_ << ".h5";
    fn = strBuffer.str();

    int32_t fg_nElems;
    int32_t fg_nNodes;
    int32_t fg_nNodesPerElem;
    std::vector<int32_t> fg_ElemTypes;
    std::vector<float64_t> fg_XNodeLocs, fg_YNodeLocs, fg_ZNodeLocs;
    std::vector<int32_t> fg_NodesInElem;

    xid    xFileId, xGroupId;
    int32_t    status;

    // Open the file and the mesh group

    status = xfOpenFile(fn.c_str(), &xFileId, true);
    std::stringstream strBuf;
    strBuf << "Could not read XMDF file " << fn;
    CHECK_HDF5_ERROR(status, strBuf.str());

    if (status >= 0) {
      status = xfOpenGroup(xFileId, pathToMesh.c_str(), &xGroupId);
    }

    // check if group was opened successfully
    if (status < 0)
      xfCloseFile(xFileId);

    strBuf.clear(); strBuf.str("");
    strBuf << "Could not open mesh group '" << pathToMesh << "'";
    CHECK_HDF5_ERROR(status, strBuf.str());

    // Get the number of elements, nodes, and Maximum number of nodes per element
    status = xfGetNumberOfElements(xGroupId, &fg_nElems);
    if (status >= 0) {
      status = xfGetNumberOfNodes(xGroupId, &fg_nNodes);
      if (status >= 0) {
        status = xfGetMaxNodesInElem(xGroupId, &fg_nNodesPerElem);
      }
    }
    if (status < 0)
    {
      xfCloseGroup(xGroupId);
      xfCloseFile(xFileId);
    }
    CHECK_HDF5_ERROR(status, "Unable to read number of elements, nodes or max elemnodes.");

    fg_ElemTypes.resize(fg_nElems);
    MESHIO_DEBUG(fg_nElems);

    status = xfReadElemTypes(xGroupId, fg_nElems, &fg_ElemTypes[0]);
    if (status < 0)
    {
      xfCloseGroup(xGroupId);
      xfCloseFile(xFileId);
    }
    CHECK_HDF5_ERROR(status, "Unable to read element types.");

    // Nodes in each element
    fg_NodesInElem.resize(fg_nElems*fg_nNodesPerElem);
    status = xfReadElemNodeIds(xGroupId, fg_nElems, fg_nNodesPerElem, &fg_NodesInElem[0]);
    if (status < 0)
    {
      xfCloseGroup(xGroupId);
      xfCloseFile(xFileId);
    }
    CHECK_HDF5_ERROR(status, "Unable to read maximum number of nodes per elem.");

    // NodeLocations
    fg_XNodeLocs.resize(fg_nNodes);
    fg_YNodeLocs.resize(fg_nNodes);
    fg_ZNodeLocs.resize(fg_nNodes);

    MESHIO_DEBUG(fg_nNodes << " " << fg_XNodeLocs.size());

    status = xfReadXNodeLocations(xGroupId, fg_nNodes, &fg_XNodeLocs[0]);
    if (status >= 0) {
      status = xfReadYNodeLocations(xGroupId, fg_nNodes, &fg_YNodeLocs[0]);
      if (status >= 0) {
        status = xfReadZNodeLocations(xGroupId, fg_nNodes, &fg_ZNodeLocs[0]);
      }
    }

    if (status < 0)
    {
      xfCloseGroup(xGroupId);
      xfCloseFile(xFileId);
    }
    CHECK_HDF5_ERROR(status, "Unable to read nodes.");

    mi_->AddNodes(fg_nNodes);

    int32_t node, elem;
    uint32_t connectIdx;
    FEType eType;

    MESHIO_INFO(fg_nNodes << " " << fg_XNodeLocs.size());
    for (node=0; node<fg_nNodes; node++) {
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

    mi_->AddElems(fg_nElems);

    connectIdx = 0;
    for (elem=0; elem<fg_nElems; elem++) {
      eType = XMDFElemType2ElemType(fg_ElemTypes[elem]);
      if((eType == ET_LINE3) ||
         (eType == ET_TRIA6) ||
         (eType == ET_QUAD8) ||
         (eType == ET_QUAD9))
        this->ReorderConnectivity(eType,
                                  false,
                                  (uint32_t*)&fg_NodesInElem[connectIdx],
                                  (uint32_t*)&fg_NodesInElem[connectIdx]);

      mi_->SetElemData(elem+1, eType, 0, (uint32_t*)&fg_NodesInElem[connectIdx]);
      connectIdx += NUM_ELEM_NODES[eType];
    }

    // Release temp memory of elems
    fg_ElemTypes.clear();
    fg_NodesInElem.clear();

    ReadRegions(xGroupId);
    ReadNamedNodes(xGroupId);
    ReadNamedElems(xGroupId);

    status = xfCloseGroup(xGroupId);
    if (status < 0)
    {
      xfCloseFile(xFileId);
    }
    CHECK_HDF5_ERROR(status, "Unable to close mesh group.");

    status = xfCloseFile(xFileId);
    CHECK_HDF5_ERROR(status, "Unable to close file.");
  }

  void XMDF::WriteMesh()
  {
    std::stringstream strBuffer;
    std::string fn;
    std::string pathToMesh = "Mesh";

    strBuffer << this->baseDir_ << baseName_ << ".h5";
    fn = strBuffer.str();

    int32_t fg_nElems;
    int32_t fg_nNodes;
    std::vector<int32_t> fg_ElemTypes;
    std::vector<float64_t> fg_XNodeLocs, fg_YNodeLocs, fg_ZNodeLocs;
    std::vector<int32_t> fg_NodesInElem;

    xid    xFileId, xGroupId;
    int32_t    status;
    int32_t    fg_Compression = 1;

    // Create the file and the mesh group
    status = xfCreateFile(fn.c_str(), &xFileId, true);
    if (status < 0) {
      std::stringstream strBuf;
      strBuf << "Could not create XMDF file " << fileName_;
      Exception ex(strBuf.str(), __FILE__, __LINE__);
      errHandler->Error(ex);
    }

    if (status >= 0) {
      status = xfCreateGroupForMesh(xFileId, pathToMesh.c_str(), &xGroupId);
    }

    if (status < 0) {
      // group was not opened successfully
      xfCloseFile(xFileId);

      std::stringstream strBuf;
      strBuf << "Could not create mesh group '" << pathToMesh << "'";
      Exception ex(strBuf.str(), __FILE__, __LINE__);
      errHandler->Error(ex);
    }

    fg_nElems = mi_->GetNumElems();
    fg_nNodes = mi_->GetNumNodes();

    // Set the number of elements and nodes
    status = xfSetNumberOfElements(xGroupId, fg_nElems);
    if (status >= 0) {
      status = xfSetNumberOfNodes(xGroupId, fg_nNodes);
    }

    if (status < 0) {
      std::stringstream strBuf;
      strBuf << "Unable to write number of elements and nodes.";
      Exception ex(strBuf.str(), __FILE__, __LINE__);
      errHandler->Error(ex);
    }

    // NodeLocations
    fg_XNodeLocs.resize(fg_nNodes);
    fg_YNodeLocs.resize(fg_nNodes);
    fg_ZNodeLocs.resize(fg_nNodes);

    for( uint32_t i=0; i < fg_nNodes; i++ )
    {
      Point p;
      
      mi_->GetNodeCoordinate(p, i+1);

      fg_XNodeLocs[i] = p[0];
      fg_YNodeLocs[i] = p[1];
      fg_ZNodeLocs[i] = p[2];
    }

    // Write node locations
    status = xfWriteXNodeLocations(xGroupId, fg_nNodes, &fg_XNodeLocs[0], fg_Compression);
    if (status >= 0) {
      status = xfWriteYNodeLocations(xGroupId, fg_nNodes, &fg_YNodeLocs[0]);
      if (status >= 0) {
        status = xfWriteZNodeLocations(xGroupId, fg_nNodes, &fg_ZNodeLocs[0]);
      }
    }

    if (status < 0) {
      std::stringstream strBuf;
      strBuf << "Unable to write nodes.";
      Exception ex(strBuf.str(), __FILE__, __LINE__);
      errHandler->Error(ex);
    }

    fg_ElemTypes.resize(fg_nElems);
    MESHIO_DEBUG(fg_nElems);

    std::vector<uint32_t> connect;
    std::vector< std::vector<uint32_t> > regionElems;
    std::vector< std::set<uint32_t, std::less<uint32_t>, std::allocator<uint32_t> > > regionNodes;
    uint32_t elemNum, maxNumNodes,numNodes;
    FEType eType;
    RegionIdType region;

    maxNumNodes = mi_->GetMaxNumNodesPerElem();
    connect.resize(maxNumNodes);
    regionElems.resize(mi_->GetNumRegions());
    regionNodes.resize(mi_->GetNumRegions());

    // Build arrays for writing data
    for( uint32_t i=0; i<fg_nElems; i++ )
    {
      elemNum = i+1;
      mi_->GetElemData(elemNum, eType, region, &connect[0]);

      numNodes = NUM_ELEM_NODES[eType];
      fg_ElemTypes[i] = ElemType2XMDFElemType(eType);

      if((eType == ET_LINE3) ||
         (eType == ET_TRIA6) ||
         (eType == ET_QUAD8) ||
         (eType == ET_QUAD9))
        ReorderConnectivity(eType, true, &connect[0], &connect[0]);

      fg_NodesInElem.insert(fg_NodesInElem.end(),
                            connect.begin(),
                            connect.end());
      regionElems[region].push_back(elemNum);
      regionNodes[region].insert(connect.begin(), connect.end());
    }

    // Write element types
    status = xfWriteElemTypes(xGroupId, fg_nElems, &fg_ElemTypes[0], fg_Compression);
    if (status < 0) {
      std::stringstream strBuf;
      strBuf << "Unable to write element types.";
      Exception ex(strBuf.str(), __FILE__, __LINE__);
      errHandler->Error(ex);
    }

    // Write element connectivity
    status = xfWriteElemNodeIds(xGroupId, fg_nElems, maxNumNodes, &fg_NodesInElem[0], fg_Compression);
    if (status < 0) {
      std::stringstream strBuf;
      strBuf << "Unable to write element connectivity.";
      Exception ex(strBuf.str(), __FILE__, __LINE__);
      errHandler->Error(ex);
    }

    WriteRegions(xGroupId, regionElems, regionNodes);
    WriteNamedNodes(xGroupId);
    WriteNamedElems(xGroupId);
  }


  // ======================================================
  // GENERAL MESH INFORMATION
  // ======================================================
  uint32_t XMDF::GetDim() {
    MESHIO_WARN("XMDF::ReadMesh() not implemented");
    return 0;
  }
  
  uint32_t XMDF::GetNumNodes(){
    MESHIO_WARN("XMDF::ReadMesh() not implemented");
    return 0;
  }
    
  uint32_t XMDF::GetNumElems(const uint32_t dim){
    MESHIO_WARN("XMDF::ReadMesh() not implemented");
    return 0;
  }
  
  uint32_t XMDF::GetNumRegions(){
    MESHIO_WARN("XMDF::ReadMesh() not implemented");
    return 0;
  }

  uint32_t XMDF::GetNumNamedNodes(){
    MESHIO_WARN("XMDF::ReadMesh() not implemented");
    return 0;
  }

  uint32_t XMDF::GetNumNamedElems(){
    MESHIO_WARN("XMDF::ReadMesh() not implemented");
    return 0;
  }
  
  // ======================================================
  // ENTITY NAME ACCESS
  // ======================================================

  void XMDF::GetAllRegionNames( std::vector<std::string> & regionNames ){
    MESHIO_WARN("XMDF::ReadMesh() not implemented");
  }
    

  // =========================================================================
  // MISCELLANEOUS METHODS
  // =========================================================================


  // ***************
  //   Type2ptElem
  // ***************

  FEType XMDF::XMDFElemType2ElemType( const int32_t type ) {

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

    // This place should never be reached!
    return ET_UNDEF;
  }

  int32_t XMDF::ElemType2XMDFElemType( const FEType type )
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
    }

    // This place should never be reached!
    return -1;
  }


  void XMDF::ReorderConnectivity( const int32_t eType,
                                  const bool toXMDF,
                                  const uint32_t* in,
                                  uint32_t* out) {
    static int32_t toXMDFIdxsLine[] = {0, 2, 1};
    static int32_t fromXMDFIdxsLine[] = {0, 2, 1};
    static int32_t toXMDFIdxsTria[] = {0, 3, 1, 4, 2, 5};
    static int32_t fromXMDFIdxsTria[] = {0, 2, 4, 1, 3, 5};
    static int32_t toXMDFIdxsQuad[] = {0, 4, 1, 5, 2, 6, 3, 7, 8};
    static int32_t fromXMDFIdxsQuad[] = {0, 2, 4, 6, 1, 3, 5, 7, 8};

    std::vector<uint32_t> tmp;
    int32_t* it;

    uint32_t numNodes = NUM_ELEM_NODES[eType];
    tmp.resize(numNodes);
    memcpy(&tmp[0], in, numNodes*sizeof(uint32_t));

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

    for(uint32_t i=0; i<numNodes; i++, it++)
    {
      out[i] = tmp[*it];
    }
  }

  hid_t XMDF::CreateGroup(const hid_t parentGroup, const std::string name)
  {
    // Add a new Group
    hid_t groupId = H5Gcreate(parentGroup, name.c_str(), -1 );

    std::stringstream strBuf;
    strBuf << "Could not create group " << name;
    CHECK_HDF5_ERROR(groupId, strBuf.str());

    return groupId;
  }

  void XMDF::ReadRegions(const hid_t meshGroup)
  {
    hid_t status, regionGroup;
    hid_t regionElemGroup;
    std::vector< uint32_t > regionElems, regionElDummy;
    std::vector< std::string > regionNames;
    std::vector< RegionIdType > regionIds;
    int32_t fg_nElems, propStringLen, number, numInt;
    char* regnames;
    int fg_Compression = 1;

    status = xfGetNumberOfElements(meshGroup, &fg_nElems);
    CHECK_HDF5_ERROR(status, "Could not read number of elements.");

    if(mi_->GetNumElems() != fg_nElems)
    {
      Exception ex( "Number of elements in mesh does not match numelems in file.", __FILE__, __LINE__);
      errHandler->Error(ex);
    }

    regionGroup = H5Gopen(meshGroup, "Regions");
    if(regionGroup < 0)
    {
      MESHIO_WARN("Could not open region group.");
      regionNames.push_back("XMDF_Region");
      mi_->AddRegions(regionNames, regionIds);
      return;
    }

    regionElemGroup = H5Gopen(regionGroup, "Elements");
    CHECK_HDF5_ERROR(regionElemGroup, "Could not open region elements group.");

    number = 0;
    status = xfGetPropertyStringLength(regionGroup, "Names", &number, &propStringLen);
    CHECK_HDF5_ERROR(status, "Could not read max region names length.");

    regnames = new char[number * propStringLen];
      
    status = xfReadPropertyString(regionGroup, "Names",
                                  number,
                                  propStringLen,
                                  regnames);
    if(status < 0) 
    {
      delete[] regnames;
      H5Gclose(regionElemGroup);
      H5Gclose(regionGroup);
    }
    CHECK_HDF5_ERROR(status, "Could not read region names.");

    // Push the region names into a vector
    for(int32_t i=0; i<number; i++)
    {
      regionNames.push_back(regnames+i*propStringLen);
    }
    delete[] regnames;

    // Create the regions
    mi_->AddRegions(regionNames, regionIds);

    for(int32_t i=0; i<number; i++)
    {
      // Get number of elements in region
      status = xfGetPropertyNumber(regionElemGroup, regionNames[i].c_str(), &numInt);
      if(status < 0) 
      {
        H5Gclose(regionElemGroup);
        H5Gclose(regionGroup);
      }
      CHECK_HDF5_ERROR(status, "Could not read number of elements in region.");

      // Read elements in region
      regionElems.resize(numInt);
      status = xfReadPropertyInt(regionElemGroup, regionNames[i].c_str(), numInt, (int32_t*)&regionElems[0]);
      if(status < 0) 
      {
        H5Gclose(regionElemGroup);
        H5Gclose(regionGroup);
      }
      CHECK_HDF5_ERROR(status, "Could not read elements in region.");

      for(int j=0; j<numInt; j++)
      {
        mi_->SetElemRegion(regionElems[j], regionIds[i]);
      }
    }

    status = H5Gclose(regionElemGroup);
    CHECK_HDF5_ERROR(status, "Could not close region elements group.");

    status = H5Gclose(regionGroup);
    CHECK_HDF5_ERROR(status, "Could not close region group.");
  }

  void XMDF::ReadNamedNodes(const hid_t meshGroup)
  {
    hid_t status, namedNodesGroup;
    hid_t nNodesGroup;
    std::vector< uint32_t > namedNodes;
    int32_t propStringLen, number, numInt;
    char* nodenames;

    namedNodesGroup = H5Gopen(meshGroup, "Named Nodes");
    if(namedNodesGroup < 0)
    {
      MESHIO_WARN("Could not open named nodes group.");
      return;
    }

    nNodesGroup = H5Gopen(namedNodesGroup, "Nodes");
    if(nNodesGroup < 0)
      H5Gclose(namedNodesGroup);
    CHECK_HDF5_ERROR(nNodesGroup, "Could not open named nodes subgroup.");

    number = 0;
    status = xfGetPropertyStringLength(namedNodesGroup, "Names", &number, &propStringLen);
    if(status < 0)
    {
      H5Gclose(namedNodesGroup);
      H5Gclose(nNodesGroup);
      MESHIO_WARN("Obviously no named nodes are defined.");
      return;
    }

    nodenames = new char[number * propStringLen];
      
    status = xfReadPropertyString(namedNodesGroup, "Names",
                                  number,
                                  propStringLen,
                                  nodenames);
    if(status < 0) 
    {
      delete[] nodenames;
      H5Gclose(namedNodesGroup);
      H5Gclose(nNodesGroup);
    }
    CHECK_HDF5_ERROR(status, "Could not read node names.");

    std::string nodename;
    for(int32_t i=0; i<number; i++)
    {
      nodename = nodenames+i*propStringLen;

      status = xfGetPropertyNumber(nNodesGroup, nodename.c_str(), &numInt);
      if(status < 0) 
      {
        delete[] nodenames;
        H5Gclose(nNodesGroup);
        H5Gclose(namedNodesGroup);
      }
      CHECK_HDF5_ERROR(status, "Could not read number of named nodes.");

      namedNodes.resize(numInt);

      status = xfReadPropertyInt(nNodesGroup, nodename.c_str(), numInt, (int32_t*)&namedNodes[0]);
      if(status < 0) 
      {
        delete[] nodenames;
        H5Gclose(nNodesGroup);
        H5Gclose(namedNodesGroup);
      }
      CHECK_HDF5_ERROR(status, "Could not read named nodes.");

      mi_->AddNamedNodes(nodename, namedNodes);
    }
    delete[] nodenames;

    status = H5Gclose(nNodesGroup);
    CHECK_HDF5_ERROR(status, "Could not close named nodes subgroup.");

    status = H5Gclose(namedNodesGroup);
    CHECK_HDF5_ERROR(status, "Could not close named nodes group.");
  }

  void XMDF::ReadNamedElems(const hid_t meshGroup)
  {
    hid_t status, namedElemsGroup;
    hid_t nElemsGroup;
    std::vector< uint32_t > namedElems;
    int32_t propStringLen, number, numInt;
    char* elemnames;

    namedElemsGroup = H5Gopen(meshGroup, "Named Elements");
    if(namedElemsGroup < 0)
    {
      MESHIO_WARN("Could not open named elements group.");
      return;
    }

    nElemsGroup = H5Gopen(namedElemsGroup, "Elements");
    if(nElemsGroup < 0)
      H5Gclose(namedElemsGroup);
    CHECK_HDF5_ERROR(nElemsGroup, "Could not open named elements subgroup.");

    number = 0;
    status = xfGetPropertyStringLength(namedElemsGroup, "Names", &number, &propStringLen);
    if(status < 0)
    {
      H5Gclose(namedElemsGroup);
      H5Gclose(nElemsGroup);
      MESHIO_WARN("Obviously no named elements are defined.");
      return;
    }

    elemnames = new char[number * propStringLen];
      
    status = xfReadPropertyString(namedElemsGroup, "Names",
                                  number,
                                  propStringLen,
                                  elemnames);
    if(status < 0) 
    {
      delete[] elemnames;
      H5Gclose(namedElemsGroup);
      H5Gclose(nElemsGroup);
    }
    CHECK_HDF5_ERROR(status, "Could not read Elem names.");

    std::string elemname;
    for(int32_t i=0; i<number; i++)
    {
      elemname = elemnames+i*propStringLen;

      status = xfGetPropertyNumber(nElemsGroup, elemname.c_str(), &numInt);
      if(status < 0) 
      {
        delete[] elemnames;
        H5Gclose(nElemsGroup);
        H5Gclose(namedElemsGroup);
      }
      CHECK_HDF5_ERROR(status, "Could not read number of named elements.");

      namedElems.resize(numInt);

      status = xfReadPropertyInt(nElemsGroup, elemname.c_str(), numInt, (int32_t*)&namedElems[0]);
      if(status < 0) 
      {
        delete[] elemnames;
        H5Gclose(nElemsGroup);
        H5Gclose(namedElemsGroup);
      }
      CHECK_HDF5_ERROR(status, "Could not read named elements.");

      mi_->AddNamedElems(elemname, namedElems);
    }
    delete[] elemnames;

    status = H5Gclose(nElemsGroup);
    CHECK_HDF5_ERROR(status, "Could not close named elements subgroup.");

    status = H5Gclose(namedElemsGroup);
    CHECK_HDF5_ERROR(status, "Could not close named elements group.");
  }

  void XMDF::WriteRegions(const hid_t meshGroup,
                          const regionElemType regionElems,
                          const regionNodeType regionNodes)
  {
    hid_t status, regionElemGroup, regionNodeGroup;
    std::vector< std::string > regionNames;
    std::vector< RegionIdType > regionIds;
    char* regnames;
    int fg_Compression = 1;
    hid_t regionGroup;
    uint32_t numRegions;

    regionGroup = H5Gopen(meshGroup, "Regions");
    if(regionGroup < 0)
      regionGroup = CreateGroup(meshGroup, "Regions");

    regionElemGroup = H5Gopen(regionGroup, "Elements");
    if(regionElemGroup < 0)
      regionElemGroup = CreateGroup(regionGroup, "Elements");

    regionNodeGroup = H5Gopen(regionGroup, "Nodes");
    if(regionNodeGroup < 0)
      regionNodeGroup = CreateGroup(regionGroup, "Nodes");

    mi_->GetRegionNames(regionNames);
    mi_->GetRegionIds(regionIds);
    numRegions = regionNames.size();
    regnames = new char[numRegions*33];
    memset(regnames, 0, numRegions*33);
    for(uint32_t i=0; i<numRegions; i++)
    {
      strncpy(&regnames[i*33], regionNames[i].c_str(), 32);
      status = xfWritePropertyInt(regionElemGroup,
                                  regionNames[i].c_str(),
                                  regionElems[i].size(),
                                  (int32_t*)&regionElems[i][0],
                                  fg_Compression);
      if(status < 0) delete[] regnames;
      CHECK_HDF5_ERROR(status, "Unable to write region elements.");
      std::vector<uint32_t> nodes(regionNodes[i].begin(),
                                regionNodes[i].end());
      status = xfWritePropertyInt(regionNodeGroup,
                                  regionNames[i].c_str(),
                                  nodes.size(),
                                  (int32_t*)&nodes[0],
                                  fg_Compression);
      if(status < 0) delete[] regnames;
      CHECK_HDF5_ERROR(status, "Unable to write region nodes.");
    }


    if(numRegions)
    {
      status = xfWritePropertyString(regionGroup, "Names", regionNames.size(), 32, regnames);
      delete[] regnames;
      regnames = NULL;
      CHECK_HDF5_ERROR(status, "Unable to write region names.");

      status = xfWritePropertyInt(regionGroup, "Ids", regionIds.size(), (int32_t*)&regionIds[0], fg_Compression);
      CHECK_HDF5_ERROR(status, "Unable to write region Ids.");
    }
    delete[] regnames;
    
    status = H5Gclose(regionElemGroup);
    CHECK_HDF5_ERROR(status, "Could not close region elements group.");

    status = H5Gclose(regionNodeGroup);
    CHECK_HDF5_ERROR(status, "Could not close region nodes group.");

    status = H5Gclose(regionGroup);
    CHECK_HDF5_ERROR(status, "Could not close regions group.");
  }

  void XMDF::WriteNamedNodes(const hid_t meshGroup)
  {
    hid_t status, nNodeGroup;
    std::vector< std::string > nodeNames;
    std::vector< uint32_t > nodes;
    char* nodenames;
    int fg_Compression = 1;
    hid_t namedNodeGroup;
    uint32_t numNamedNodes;

    namedNodeGroup = H5Gopen(meshGroup, "Named Nodes");
    if(namedNodeGroup < 0)
      namedNodeGroup = CreateGroup(meshGroup, "Named Nodes");

    nNodeGroup = H5Gopen(namedNodeGroup, "Nodes");
    if(nNodeGroup < 0)
      nNodeGroup = CreateGroup(namedNodeGroup, "Nodes");

    mi_->GetListNodeNames(nodeNames);
    numNamedNodes = nodeNames.size();
    nodenames = new char[numNamedNodes*33];
    memset(nodenames, 0, numNamedNodes*33);
    for(uint32_t i=0; i<numNamedNodes; i++)
    {
      strncpy(&nodenames[i*33], nodeNames[i].c_str(), 32);
      mi_->GetNodesByName(nodes, nodeNames[i]);
      status = xfWritePropertyInt(nNodeGroup,
                                  nodeNames[i].c_str(),
                                  nodes.size(),
                                  (int32_t*)&nodes[0],
                                  fg_Compression);
      if(status < 0) delete[] nodenames;
      CHECK_HDF5_ERROR(status, "Unable to write named nodes.");
    }

    if(numNamedNodes)
    {
      status = xfWritePropertyString(namedNodeGroup, "Names", nodeNames.size(), 32, nodenames);
      delete[] nodenames;
      nodenames = NULL;
      CHECK_HDF5_ERROR(status, "Unable to write node names.");

      status = xfWritePropertyInt(namedNodeGroup, "Ids", nodeNames.size(), (int32_t*)&nodes[0], fg_Compression);
      CHECK_HDF5_ERROR(status, "Unable to write node Ids.");
    }
    delete[] nodenames;
      
    status = H5Gclose(nNodeGroup);
    CHECK_HDF5_ERROR(status, "Could not close named nodes group.");

    status = H5Gclose(namedNodeGroup);
    CHECK_HDF5_ERROR(status, "Could not close named nodes group.");
  }

  void XMDF::WriteNamedElems(const hid_t meshGroup)
  {
    hid_t status, nElemGroup;
    std::vector< std::string > elemNames;
    std::vector< uint32_t > elems;
    char* elemnames;
    int fg_Compression = 1;
    hid_t namedElemGroup;
    uint32_t numNamedElems;


    namedElemGroup = H5Gopen(meshGroup, "Named Elements");
    if(namedElemGroup < 0)
      namedElemGroup = CreateGroup(meshGroup, "Named Elements");

    nElemGroup = H5Gopen(namedElemGroup, "Elements");
    if(nElemGroup < 0)
      nElemGroup = CreateGroup(namedElemGroup, "Elements");

    mi_->GetListElemNames(elemNames);
    numNamedElems = elemNames.size();
    elemnames = new char[numNamedElems*33];
    memset(elemnames, 0, numNamedElems*33);
    for(uint32_t i=0; i<numNamedElems; i++)
    {
      strncpy(&elemnames[i*33], elemNames[i].c_str(), 32);
      mi_->GetElemNumsByName(elems, elemNames[i]);
      status = xfWritePropertyInt(nElemGroup,
                                  elemNames[i].c_str(),
                                  elems.size(),
                                  (int32_t*)&elems[0],
                                  fg_Compression);
      if(status < 0) delete[] elemnames;
      CHECK_HDF5_ERROR(status, "Unable to write named elements.");
    }

    if(numNamedElems) 
    {
      status = xfWritePropertyString(namedElemGroup, "Names", elemNames.size(), 32, elemnames);
      delete[] elemnames;
      elemnames = NULL;
      CHECK_HDF5_ERROR(status, "Unable to write element names.");

      status = xfWritePropertyInt(namedElemGroup, "Ids", elemNames.size(), (int32_t*)&elems[0], fg_Compression);
      CHECK_HDF5_ERROR(status, "Unable to write Elem Ids.");
    }
    delete[] elemnames;
      
    status = H5Gclose(nElemGroup);
    CHECK_HDF5_ERROR(status, "Could not close named elements group.");

    status = H5Gclose(namedElemGroup);
    CHECK_HDF5_ERROR(status, "Could not close named elements group.");
  }

}

/// Local Variables:
/// mode: C++
/// c-basic-offset: 2
/// End:
