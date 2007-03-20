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
#include "simOutputXMDF.hh"

namespace fs = boost::filesystem;

namespace CoupledField {

 #define CHECK_HDF5_ERROR(HID_T, STR) { if(HID_T < 0) {                 \
        std::ostringstream ostr;                                        \
        ostr << STR;                                                    \
        Exception ex(NULL, __FILE__, __LINE__, ostr.str().c_str()); } }


  SimOutputXMDF::SimOutputXMDF(std::string fileName, ParamNode * inputNode) :
      SimOutput(fileName, inputNode)
  {
    ENTER_FCN( "SimOutputXMDF::XMDF" );
    fileName_ = fileName;
    fileId_ = -1;
    formatName_ = "hdf5";
    dirName_ = "simoutput_hdf5";

    capabilities_.insert( MESH );
    //    capabilities_.insert( MESH_RESULTS );
    gridWritten_ = false;
    
    
    // Do not print HDF5 exceptions by default
    H5::Exception::dontPrint();
  }


  SimOutputXMDF::~SimOutputXMDF() {
    ENTER_FCN( "SimOutputXMDF::~XMDF" );
    if(fileId_ > 0)
        xfCloseFile(fileId_);
  }

    void SimOutputXMDF::BeginMultiSequenceStep( UInt step, AnalysisType type ) 
    {
    }
    
    void SimOutputXMDF::RegisterResult( shared_ptr<BaseResult> sol,
                         UInt saveBegin, UInt saveInc,
                         UInt saveEnd ) 
    {
    }

    void SimOutputXMDF::BeginStep( UInt stepNum, Double stepVal ) 
    {
    }

    void SimOutputXMDF::AddResult( shared_ptr<BaseResult> sol ) 
    {
    }

    void SimOutputXMDF::FinishStep( ) 
    {
    }

    void SimOutputXMDF::FinishMultiSequenceStep( ) 
    {
    }

    void SimOutputXMDF::Finalize() 
    {
    }


  void SimOutputXMDF::InitModule()
  {
    std::string pathsep;
    std::string filename;
    std::ostringstream strBuffer;
      
    try 
    {
      fs::create_directory( dirName_ );
      pathsep = fs::path("/").native_directory_string();

    } catch (std::exception &ex)
    {
      EXCEPTION(ex.what());
    }

    strBuffer << dirName_ << pathsep << fileName_ << ".h5";
    filename = strBuffer.str();
    
    xid  xGroupId = -1;
    Integer    status;
    std::string pathToMesh = "Mesh";
    std::stringstream strBuf;

    // Create the file and the mesh group
    status = xfCreateFile(filename.c_str(), &fileId_, true);
    strBuf << "Could not create XMDF file " << fileName_;
    CHECK_HDF5_ERROR(status, strBuf.str());
    
    status = xfCreateGroupForMesh(fileId_, pathToMesh.c_str(), &xGroupId);
    strBuf << "Could not create mesh group '" << pathToMesh << "'";
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

  void SimOutputXMDF::WriteGrid()
  {
    Integer fg_nElems;
    Integer fg_nNodes;
    std::vector<Integer> fg_ElemTypes;
    std::vector<Double> fg_XNodeLocs, fg_YNodeLocs, fg_ZNodeLocs;
    std::vector<Integer> fg_NodesInElem;
    Integer    status;
    Integer    fg_Compression = 1;
    std::stringstream strBuf;
    H5::Group mGroup;
    xid meshGroupId;
    
    if(!gridWritten_)
        InitModule();
    else
        return;
    
    try
    {
      mGroup = mainRoot_.openGroup("Mesh");
    }
    catch (H5::Exception& h5ex)
    {
      EXCEPTION(h5ex.getCDetailMsg());
    }

    meshGroupId = mGroup.getLocId();
    
    fg_nElems = ptGrid_->GetNumElems();
    fg_nNodes = ptGrid_->GetNumNodes();

    // Set the number of elements and nodes
    status = xfSetNumberOfElements(meshGroupId, fg_nElems);
    if (status >= 0) {
      status = xfSetNumberOfNodes(meshGroupId, fg_nNodes);
    }

    if (status < 0) {
       EXCEPTION("Unable to write number of elements and nodes.");
    }

    // NodeLocations
    fg_XNodeLocs.resize(fg_nNodes);
    fg_YNodeLocs.resize(fg_nNodes);
    fg_ZNodeLocs.resize(fg_nNodes);

    for (UInt i = 0; i < fg_nNodes; i++) {
      Point p;
      
      ptGrid_->GetNodeCoordinate(p, i+1);
        
      fg_XNodeLocs[i] = p[0];
      fg_YNodeLocs[i] = p[1];
      fg_ZNodeLocs[i] = p[2];
    }

    // Write node locations
    status = xfWriteXNodeLocations(meshGroupId, fg_nNodes, &fg_XNodeLocs[0], fg_Compression);
    if (status >= 0) {
      status = xfWriteYNodeLocations(meshGroupId, fg_nNodes, &fg_YNodeLocs[0]);
      if (status >= 0) {
        status = xfWriteZNodeLocations(meshGroupId, fg_nNodes, &fg_ZNodeLocs[0]);
      }
    }

    if (status < 0) {
      EXCEPTION("Unable to write nodes.");
    }

    fg_ElemTypes.resize(fg_nElems);

    std::vector<UInt> connect;
    std::vector< std::vector<UInt> > regionElems;
    std::vector< std::set<UInt, std::less<UInt>,
          std::allocator<UInt> > > regionNodes;
    UInt elemNum, maxNumNodes,numNodes;
    UInt numRegions;
    FEType eType;
    RegionIdType region;

    maxNumNodes = ptGrid_->GetMaxNumNodesPerElem();
    connect.resize(maxNumNodes);
    numRegions = ptGrid_->GetNumRegions();
    regionElems.resize(numRegions);
    regionNodes.resize(numRegions);

    // Build arrays for writing data
    for( UInt i=0; i < fg_nElems; i++ )
    {
        elemNum = i+1;
        std::fill(connect.begin(), connect.end(), 0);
      ptGrid_->GetElemData(elemNum, eType, region, &connect[0]);

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
      regionNodes[region].insert(connect.begin(), connect.begin()+numNodes);
    }

    // Write element types
    status = xfWriteElemTypes(meshGroupId,
                              fg_nElems,
                              (int*) &fg_ElemTypes[0],
                              fg_Compression);
    if (status < 0) {
      EXCEPTION("Unable to write element types.");
    }

    // Write element connectivity
    status = xfWriteElemNodeIds(meshGroupId,
                                fg_nElems,
                                maxNumNodes,
                                (int*) &fg_NodesInElem[0],
                                fg_Compression);
    if (status < 0) {
      EXCEPTION("Unable to write element connectivity.");
    }

    try 
    {
      WriteRegions(mGroup, regionElems, regionNodes);
      WriteNamedNodes(mGroup);
      WriteNamedElems(mGroup);
    }
    catch (H5::Exception& h5ex)
    {
      EXCEPTION(h5ex.getCDetailMsg());
    }

    gridWritten_ = true;
    
  }



  // =========================================================================
  // MISCELLANEOUS METHODS
  // =========================================================================


  // ***************
  //   Type2ptElem
  // ***************

  FEType SimOutputXMDF::XMDFElemType2ElemType( const Integer type ) {

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

  Integer SimOutputXMDF::ElemType2XMDFElemType( const FEType type )
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


  void SimOutputXMDF::ReorderConnectivity( const Integer eType,
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

    void SimOutputXMDF::WriteRegions(const H5::Group& meshGroup,
                                   const RegionElemType& regionElems,
                                   const RegionNodeType& regionNodes)

  {
#if 0
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
#endif
  }

  void SimOutputXMDF::WriteNamedNodes(const H5::Group& meshGroup)
  {
#if 0
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
#endif
  }

  void SimOutputXMDF::WriteNamedElems(const H5::Group& meshGroup)
  {
#if 0
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
#endif
  }


}

/// Local Variables:
/// mode: C++
/// c-basic-offset: 2
/// End:
