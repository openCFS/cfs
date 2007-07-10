// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

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
#include <boost/tokenizer.hpp>

#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/exception.hh"
#include "simInputXMDF.hh"

namespace fs = boost::filesystem;

namespace CoupledField {

  // declare logging stream
  DECLARE_LOG(simInputXMDF)
  DEFINE_LOG(simInputXMDF, "SimInputXMDF")

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
    ParamNode *readRegionNode = NULL;
    
    // Change defaults according to XML file
    if(myParam_->Get("generateRegionNodes", false)->AsBool())
    {
        genRegionNodes_ = true;
    }

    readRegionNode = myParam_->Get("readRegions", false);
    if(readRegionNode)
    {
       std::string readRegions = myParam_->Get("readRegions", false)->AsString();
       std::cout << std::endl << "readRegions " << readRegions << std::endl;
       typedef boost::tokenizer<char_separator<char> > Tok;
       boost::char_separator<char> sep(";| ");
       Tok t(readRegions, sep);
       Tok::iterator it, end;
       it = t.begin();
       end = t.end();
    
       for( ; it != end; it++)
         readRegions_.push_back(*it);
    }
    else
    {
         readRegions_.push_back("all");
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
      EXCEPTION("Received exception: " << ex.what());
      return;
    }
    
    LOG_TRACE(simInputXMDF) << "fileName_: " << fileName_;
    LOG_TRACE(simInputXMDF) << "baseDir_: " << baseDir_;
    LOG_TRACE(simInputXMDF) << "baseName_: " << baseName_;

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

    // Open mesh group
    try
    {
      mGroup = mainRoot_.openGroup("Mesh");
    }
    catch (H5::Exception& h5ex)
    {
        EXCEPTION(h5ex.getCDetailMsg());
    }

    meshGroupId = mGroup.getLocId();

    // Read infos about mesh
    if(!statsRead_)
      ReadMeshStats(mGroup);

    // If all regions are to be read set list of readRegions accordingly.
    if(readRegions_[0] == "all")
      readRegions_ = regionNames_;
    
    // Check if all readRegions_ can be found in file.
    for(UInt i=0; i<readRegions_.size(); i++)
    {
      if(std::find(regionNames_.begin(),
                   regionNames_.end(),
                   readRegions_[i]) == regionNames_.end())
          EXCEPTION("Region " << readRegions_[i] << " specified for"
                    " reading does not exist." );
    }
    
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
    LOG_TRACE(simInputXMDF) << "fg_nElems: " << fg_nElems;

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

    LOG_TRACE(simInputXMDF) << "fg_nNodes: " << fg_nNodes << " " << fg_XNodeLocs.size();

    status = xfReadXNodeLocations(meshGroupId, fg_nNodes, &fg_XNodeLocs[0]);
    if (status >= 0) {
      status = xfReadYNodeLocations(meshGroupId, fg_nNodes, &fg_YNodeLocs[0]);
      if (status >= 0) {
        status = xfReadZNodeLocations(meshGroupId, fg_nNodes, &fg_ZNodeLocs[0]);
      }
    }
    CHECK_HDF5_ERROR(status, "Unable to read nodes.");

    /*
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
    */

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
    LOG_TRACE(simInputXMDF) << "SimInputXMDF::ReadMesh() not implemented";
    return 0;
  }
  
  UInt SimInputXMDF::GetNumNodes(){
    LOG_TRACE(simInputXMDF) << "SimInputXMDF::ReadMesh() not implemented";
    return 0;
  }
    
  UInt SimInputXMDF::GetNumElems(const Integer dim){
    LOG_TRACE(simInputXMDF) << "SimInputXMDF::ReadMesh() not implemented";
    return 0;
  }
  
  UInt SimInputXMDF::GetNumRegions(){
    LOG_TRACE(simInputXMDF) << "SimInputXMDF::ReadMesh() not implemented";
    return 0;
  }

  UInt SimInputXMDF::GetNumNamedNodes(){
    LOG_TRACE(simInputXMDF) << "SimInputXMDF::ReadMesh() not implemented";
    return 0;
  }

  UInt SimInputXMDF::GetNumNamedElems(){
    LOG_TRACE(simInputXMDF) << "SimInputXMDF::ReadMesh() not implemented";
    return 0;
  }
  
  // ======================================================
  // ENTITY NAME ACCESS
  // ======================================================

  void SimInputXMDF::GetAllRegionNames( std::vector<std::string> & regionNames ){
    LOG_TRACE(simInputXMDF) << "SimInputXMDF::ReadMesh() not implemented";
  }

  void SimInputXMDF::GetRegionNamesOfDim( StdVector<std::string> & regionNames,
                                   const UInt dim )
  {
    EXCEPTION("SimInputXMDF::GetRegionNamesOfDim() not implemented");
  }
    
  void SimInputXMDF::GetNodeNames( StdVector<std::string> & nodeNames )
  {
    EXCEPTION("SimInputXMDF::GetNodeNames() not implemented");
  }
  
  void SimInputXMDF::GetElemNames( StdVector<std::string> & elemNames )
  {
    EXCEPTION("SimInputXMDF::GetElemNames() not implemented");
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

    /*
    status = xfGetNumberOfElements(meshGroup.getLocId(), (int*) &fg_nElems);
    CHECK_HDF5_ERROR(status, "Could not read number of elements.");

    if(mi_->GetNumElems() != fg_nElems)
    {
        EXCEPTION ("Number of elements in mesh does not "    \
                   "match numelems in file.");
    }
    */

    if(!statsRead_)
      ReadMeshStats(meshGroup);

    /*
    if(numRegions_ == 0)
    {
      mi_->AddRegions(regionNames_, regionIds);
      numRegions_ = 1;
      return;
    }
    */

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
    mi_->AddRegions(readRegions_, regionIds);
    UInt numReadRegions = readRegions_.size();
    FEType eType;
    UInt connectIdx = 0;

    for(Integer i=0; i<numReadRegions; i++)
    {
      // Get number of elements in region
      status = xfGetPropertyNumber(regionElemGroup.getLocId(),
                                   readRegions_[i].c_str(),
                                   (int*) &numInt);
      CHECK_HDF5_ERROR(status, "Could not read number of elements in region.");
      
      // Read elements in region
      regionElems.resize(numInt);
      status = xfReadPropertyInt(regionElemGroup.getLocId(),
                                 readRegions_[i].c_str(),
                                 numInt,
                                 (int*) &regionElems[0]);
      CHECK_HDF5_ERROR(status, "Could not read elements in region.");
      
      readElemSet_.insert(regionElems.begin(), regionElems.end());
      
      for(int j=0; j<numInt; j++)
      {
        elemIdx = regionElems[j]-1;
        connectIdx = elemIdx*fg_nNodesPerElem;

        numElemNodes = NUM_ELEM_NODES[XMDFElemType2ElemType(fg_ElemTypes[elemIdx])];
        UInt* ptElemNodes = (UInt*)&fg_NodesInElem[connectIdx];
        readNodeSet_.insert( ptElemNodes,
                             ptElemNodes + numElemNodes);
      }
    }

    UInt nodeOffset = mi_->GetNumNodes();
    UInt numReadNodes = readNodeSet_.size();
    std::set<UInt>::const_iterator nodeIt, nodeEnd;
    UInt node = nodeOffset + 1;
    
    mi_->AddNodes(numReadNodes);

    for(nodeIt = readNodeSet_.begin(),
       nodeEnd = readNodeSet_.end();
        nodeIt != nodeEnd; nodeIt++)
    {
        Point p;
        readNodeMap_[*nodeIt] = node;

        p[0] = fg_XNodeLocs[(*nodeIt) - 1];
        p[1] = fg_YNodeLocs[(*nodeIt) - 1];
        p[2] = fg_ZNodeLocs[(*nodeIt) - 1];
    
        mi_->SetNodeCoordinate(node, p);
        node++;
    }

    // Release temp memory of nodes
    fg_XNodeLocs.clear();
    fg_YNodeLocs.clear();
    fg_ZNodeLocs.clear();

    for(Integer i=0; i<numReadRegions; i++)
    {
      UInt elem = mi_->GetNumElems() + 1;

      // Get number of elements in region
      status = xfGetPropertyNumber(regionElemGroup.getLocId(),
                                   readRegions_[i].c_str(),
                                   (int*) &numInt);
      CHECK_HDF5_ERROR(status, "Could not read number of elements in region.");
      
      // Read elements in region
      regionElems.resize(numInt);
      status = xfReadPropertyInt(regionElemGroup.getLocId(),
                                 readRegions_[i].c_str(),
                                 numInt,
                                 (int*) &regionElems[0]);
      CHECK_HDF5_ERROR(status, "Could not read elements in region.");

      mi_->AddElems(numInt);
    
      if(genRegionNodes_)
          regionNodeSet.clear();

      readElemSet_.insert(regionElems.begin(), regionElems.end());
      
      for(int j=0; j<numInt; j++)
      {
        readElemMap_[regionElems[j]] = elem;
        elemIdx = regionElems[j]-1;
        connectIdx = elemIdx*fg_nNodesPerElem;

        eType = XMDFElemType2ElemType(fg_ElemTypes[elemIdx]);
        if((eType == ET_LINE3) ||
           (eType == ET_TRIA6) ||
           (eType == ET_QUAD8) ||
           (eType == ET_QUAD9))
          this->ReorderConnectivity(eType,
                                    false,
                                    (UInt*)&fg_NodesInElem[connectIdx],
                                    (UInt*)&fg_NodesInElem[connectIdx]);

        numElemNodes = NUM_ELEM_NODES[XMDFElemType2ElemType(fg_ElemTypes[elemIdx])];

        for(int j=0; j<numElemNodes; j++)
          fg_NodesInElem[connectIdx+j] = readNodeMap_[fg_NodesInElem[connectIdx+j]];

        mi_->SetElemData(elem, eType, regionIds[i],
                         (UInt*)&fg_NodesInElem[connectIdx]);
        elem++;
          
        if(genRegionNodes_)
        {
          UInt* ptElemNodes = (UInt*)&fg_NodesInElem[connectIdx];
              
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

    fg_ElemTypes.clear();
    fg_NodesInElem.clear();
    
    regionElemGroup.close();
    regionGroup.close();
  }

  void SimInputXMDF::ReadNamedNodes(const H5::Group& meshGroup)
  {
    H5::Group namedNodesGroup, nNodesGroup;
    hid_t status;
    std::vector< UInt > namedNodes;
    EntitySet namedNodeSet;
    EntitySet result;
    Integer numInt;

    if(!statsRead_)
      ReadMeshStats(meshGroup);

    // Check if named nodes exist.
    if(!nodeNames_.size())
      return;
    
    // Try to open named nodes group.
    try 
    {
      namedNodesGroup = meshGroup.openGroup("Named Nodes");
      nNodesGroup = namedNodesGroup.openGroup("Nodes");
    } catch (H5::Exception& h5ex)
    {
      EXCEPTION(h5ex.getCDetailMsg ());
    }
    
    // Iterate over all named node datasets.
    for(UInt i=0; i<nodeNames_.size(); i++)
    {
      // Get number of named nodes.
      status = xfGetPropertyNumber(nNodesGroup.getLocId(),
                                   nodeNames_[i].c_str(),
                                   (int*) &numInt);
      if(status < 0)
      {
        nNodesGroup.close();
        namedNodesGroup.close();
      }
      CHECK_HDF5_ERROR(status, "Could not read number of named nodes.");

      // Resize array.
      namedNodes.resize(numInt);

      // Read named nodes into array.
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

      // We need to check if the named nodes make sense in the grid.
      // We do this by calculating the intersection of the named node
      // set and the read nodes set.
      namedNodeSet.clear();
      result.clear();

      // Fill named node set.
      namedNodeSet.insert(namedNodes.begin(), namedNodes.end());

      std::insert_iterator<EntitySet> res_ins (result, result.begin());

      // Calculate intersection.
      std::set_intersection (readNodeSet_.begin(), readNodeSet_.end(), 
                             namedNodeSet.begin(), namedNodeSet.end(), res_ins);

      // If there is an intersection add its members as named nodes.
      if(result.size())
      {
        std::copy(result.begin(), result.end(), namedNodes.begin());
        namedNodes.resize(result.size());
        
        for(UInt j=0,
                 n=namedNodes.size();
            j<n; j++)
        {
          namedNodes[j] = readNodeMap_[namedNodes[j]];
        }
        
        mi_->AddNamedNodes(nodeNames_[i], namedNodes);
      }
    }

    // Clear temporary storage.
    readNodeSet_.clear();
    readNodeMap_.clear();
    
    nNodesGroup.close();
    namedNodesGroup.close();
  }

  void SimInputXMDF::ReadNamedElems(const H5::Group& meshGroup)
  {
    H5::Group namedElemsGroup, nElemsGroup;
    hid_t status;
    std::vector< UInt > namedElems;
    EntitySet namedElemSet;
    EntitySet result;
    Integer numInt;

    if(!statsRead_)
      ReadMeshStats(meshGroup);

    // Check if named elems exist.
    if(!elemNames_.size())
      return;
    
    // Try to open named elements group.
    try 
    {
      namedElemsGroup = meshGroup.openGroup("Named Elements");
      nElemsGroup = namedElemsGroup.openGroup("Elements");
    } catch (H5::Exception& h5ex)
    {
      EXCEPTION(h5ex.getCDetailMsg ());
    }
    
    // Iterate over all named elem datasets.
    for(UInt i=0; i<elemNames_.size(); i++)
    {
      // Get number of named elems.
      status = xfGetPropertyNumber(nElemsGroup.getLocId(),
                                   elemNames_[i].c_str(),
                                   (int*) &numInt);
      if(status < 0)
      {
        nElemsGroup.close();
        namedElemsGroup.close();
      }
      CHECK_HDF5_ERROR(status, "Could not read number of named elems.");

      // Resize array.
      namedElems.resize(numInt);

      // Read named elems into array.
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

      // We need to check if the named elems make sense in the grid.
      // We do this by calculating the intersection of the named elem
      // set and the read elems set.
      namedElemSet.clear();
      result.clear();

      // Fill named elem set.
      namedElemSet.insert(namedElems.begin(), namedElems.end());

      std::insert_iterator<EntitySet> res_ins (result, result.begin());

      // Calculate intersection.
      std::set_intersection (readElemSet_.begin(), readElemSet_.end(), 
                             namedElemSet.begin(), namedElemSet.end(), res_ins);

      // If there is an intersection add its members as named elems.
      if(result.size())
      {
        std::copy(result.begin(), result.end(), namedElems.begin());
        namedElems.resize(result.size());
        
        for(UInt j=0,
                 n=namedElems.size();
            j<n; j++)
        {
          namedElems[j] = readElemMap_[namedElems[j]];
        }
        
        mi_->AddNamedElems(elemNames_[i], namedElems);
      }
    }

    // Clear temporary storage.
    readElemSet_.clear();
    readElemMap_.clear();
    
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
      LOG_TRACE(simInputXMDF) << "Could not open region group.";
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
      
        LOG_TRACE(simInputXMDF) << regionNames_[i] << " " << regionDims_[i];
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
      LOG_TRACE(simInputXMDF) << "Could not open named nodes group.";
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
      
        LOG_TRACE(simInputXMDF) << "Named Nodes: " << nodeNames_[i];
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
      LOG_TRACE(simInputXMDF) << "Could not open named elems group.";
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
      
        LOG_TRACE(simInputXMDF) << "Named Elems: " << elemNames_[i];
      }
    }
    
    statsRead_ = true;
  }

}
