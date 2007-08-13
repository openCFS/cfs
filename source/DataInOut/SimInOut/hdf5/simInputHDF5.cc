// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/exception.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/tokenizer.hpp>

#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/exception.hh"
#include "simInputHDF5.hh"
#include "hdf5io.hh"

namespace fs = boost::filesystem;

namespace CoupledField {

  // declare logging stream
  DECLARE_LOG(simInputHdf5)
  DEFINE_LOG(simInputHdf5, "simInput.hdf5")

#define H5_CATCH( STR )                                                 \
  catch (H5::Exception& h5Ex ) {                                        \
    EXCEPTION( STR << ":\n" << h5Ex.getCDetailMsg() );                  \
  }

  SimInputHDF5::SimInputHDF5(std::string fileName, ParamNode * inputNode) :
      SimInput(fileName, inputNode)
  {
    ENTER_FCN( "SimInputHDF5::XMDF" );
    mi_ = NULL;
    statsRead_ = false;
    fileName_ = fileName;
    genRegionNodes_ = false;
    ParamNode *readRegionNode = NULL;
    
    // Change defaults according to XML file
    if(myParam_->Get("generateRegionNodes", false)->AsBool()) {
      genRegionNodes_ = true;
    }
    
    readRegionNode = myParam_->Get("readRegions", false);
    if(readRegionNode) {
      std::string readRegions = 
        myParam_->Get("readRegions", false)->AsString();

      typedef boost::tokenizer<char_separator<char> > Tok;
      boost::char_separator<char> sep(";| ");
      Tok t(readRegions, sep);
      Tok::iterator it, end;
      it = t.begin();
      end = t.end();
      
      for( ; it != end; it++)
        readRegions_.push_back(*it);
    } else {
      readRegions_.push_back("all");
    }

    // Do not print HDF5 exceptions by default
    H5::Exception::dontPrint();

  }

  SimInputHDF5::~SimInputHDF5() {
    ENTER_FCN( "SimInputHDF5::~XMDF" );
    
    mainRoot_.close();
    mainFile_.close();
  }

  void SimInputHDF5::InitModule(Grid *mi)
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
    
    LOG_TRACE(simInputHdf5) << "fileName_: " << fileName_;
    LOG_TRACE(simInputHdf5) << "baseDir_: " << baseDir_;
    LOG_TRACE(simInputHdf5) << "baseName_: " << baseName_;

    statsRead_ = false;
    
    mi_ = mi;

    try {
      mainFile_ = H5::H5File( fileName_, H5F_ACC_RDONLY );
    } H5_CATCH( "Could not open XMDF file '" << fileName_ << "'" );

    try {
      mainRoot_ = mainFile_.openGroup("/");
    } H5_CATCH( "Could not open main root" );
  }

  void SimInputHDF5::ReadMesh()
  {
    H5::Group mGroup, nodeGroup;

    // Open mesh group
    try{
      mGroup = mainRoot_.openGroup("Mesh");
    } H5_CATCH( "Could not open mesh group" );

    // Read infos about mesh
    if(!statsRead_)
      ReadMeshStats(mGroup);

    // If all regions are to be read set list of readRegions accordingly.
    if(readRegions_[0] == "all") {
      readRegions_ = regionNames_;
    }
    
    // Check if all readRegions_ can be found in file.
     for(UInt i=0; i<readRegions_.size(); i++) {
       if(std::find(regionNames_.begin(),
                    regionNames_.end(),
                    readRegions_[i]) == regionNames_.end())
         EXCEPTION("Region " << readRegions_[i] << " specified for"
                   " reading does not exist." );
     }
    
     // ========================
     //  READ NODAL INFORMATION
     // ========================

     // get the number of nodes
     UInt numNodes;
     try{ 
       nodeGroup = mGroup.openGroup( "Nodes") ;
     } H5_CATCH( "Could not open Elements / Nodes group" );

     H5IO::ReadAttribute( nodeGroup, "NumNodes", numNodes );

     // read node coordinates
     std::vector<Double> nodeCoords;
     H5IO::ReadArray( nodeGroup, "Coordinates", nodeCoords );

     // add nodes to grid
     mi_->AddNodes( numNodes );

     UInt pos = 0;
     for( UInt i = 0; i < numNodes; i++ ) {
       Vector<Double> p(3);
       p[0] = nodeCoords[pos++];
       p[1] = nodeCoords[pos++];
       p[2] = nodeCoords[pos++];
       mi_->SetNodeCoordinate( i+1, p );
     }
     
     // read region, element and named entity informaion
     ReadRegions(mGroup);
     ReadNamedNodes(mGroup);
     ReadNamedElems(mGroup);
  }


  // ======================================================
  // GENERAL MESH INFORMATION
  // ======================================================
  UInt SimInputHDF5::GetDim() {
    LOG_TRACE(simInputHdf5) << "SimInputHDF5::ReadMesh() not implemented";
    return 0;
  }
  
  UInt SimInputHDF5::GetNumNodes(){
    LOG_TRACE(simInputHdf5) << "SimInputHDF5::ReadMesh() not implemented";
    return 0;
  }
    
  UInt SimInputHDF5::GetNumElems(const Integer dim){
    LOG_TRACE(simInputHdf5) << "SimInputHDF5::ReadMesh() not implemented";
    return 0;
  }
  
  UInt SimInputHDF5::GetNumRegions(){
    LOG_TRACE(simInputHdf5) << "SimInputHDF5::ReadMesh() not implemented";
    return 0;
  }

  UInt SimInputHDF5::GetNumNamedNodes(){
    LOG_TRACE(simInputHdf5) << "SimInputHDF5::ReadMesh() not implemented";
    return 0;
  }

  UInt SimInputHDF5::GetNumNamedElems(){
    LOG_TRACE(simInputHdf5) << "SimInputHDF5::ReadMesh() not implemented";
    return 0;
  }
  
  // ======================================================
  // ENTITY NAME ACCESS
  // ======================================================

  void SimInputHDF5::GetAllRegionNames( std::vector<std::string> & regionNames ){
    LOG_TRACE(simInputHdf5) << "SimInputHDF5::ReadMesh() not implemented";
  }

  void SimInputHDF5::GetRegionNamesOfDim( StdVector<std::string> & regionNames,
                                   const UInt dim )
  {
    EXCEPTION("SimInputHDF5::GetRegionNamesOfDim() not implemented");
  }
    
  void SimInputHDF5::GetNodeNames( StdVector<std::string> & nodeNames )
  {
    EXCEPTION("SimInputHDF5::GetNodeNames() not implemented");
  }
  
  void SimInputHDF5::GetElemNames( StdVector<std::string> & elemNames )
  {
    EXCEPTION("SimInputHDF5::GetElemNames() not implemented");
  }

  // =========================================================================
  // MISCELLANEOUS METHODS
  // =========================================================================
  void SimInputHDF5::ReadRegions(const H5::Group& meshGroup)
  {

    // ==================================
    //  Read General Element Information
    // ==================================
    H5::Group elemGroup;
    try {
      elemGroup = meshGroup.openGroup( "Elements" );
    } H5_CATCH( "Could not open 'Elements' group" );

    
    // read number of elements
    UInt numElems = 0;
    H5IO::ReadAttribute( elemGroup, "NumElems", numElems );

    // read maximum number of nodes per elements
    UInt numNodesPerElem = 
      H5IO::GetArrayDims( elemGroup, "Connectivity")[1];
    
    // read element types
    std::vector<Integer> elemTypes;
    H5IO::ReadArray( elemGroup, "Types", elemTypes );
    
    // read nodes per element
    std::vector<UInt> globConnect;
    H5IO::ReadArray( elemGroup, "Connectivity", globConnect );
    
    elemGroup.close();

    // =========================
    //  Add Elements Per Reigon
    // =========================
    
    // ensure, that region names are already read in
    if( !statsRead_ )
      ReadMeshStats( meshGroup );

    H5::Group regionGroup, actRegion;
    try {
      regionGroup = meshGroup.openGroup( "Regions" );
          } H5_CATCH( "Could not open 'Regions' group" );

    // pass region names to grid and obtain RegionIds
    std::vector<RegionIdType> regionIds;
    mi_->AddRegions(readRegions_, regionIds);

    mi_->AddElems( numElems );

    // iterate over all regions
    for( UInt i = 0; i < readRegions_.size(); i++ ) {

      try {
        actRegion = regionGroup.openGroup( readRegions_[i] );
      } H5_CATCH( "Could not open group for region '" <<
                  regionNames_[i] << "'" );

      // read element numbers for this region
      std::vector<UInt> regionElems;
      H5IO::ReadArray( actRegion, "Elements", regionElems );

      // pass for each element the definition to the grid      
      RegionIdType actRegionId = regionIds[i];
      for( UInt iElem = 0; iElem < regionElems.size(); iElem ++ ) {

        UInt elemNum = regionElems[iElem];
        FEType type = (FEType) elemTypes[elemNum-1];
        UInt * connect = &globConnect[numNodesPerElem*(elemNum-1)];
        mi_->SetElemData( elemNum, type, actRegionId, connect );
      }

       // check, if nodes nodes of region should be additionally added
      // to list of named nodes
      if( genRegionNodes_) {
        
        // read nodes of region
        std::vector<UInt> regionNodes;
        H5IO::ReadArray( actRegion, "Nodes", regionNodes );

        // add nodes as named nodes
        mi_->AddNamedNodes( readRegions_[i]+"_Nodes", regionNodes );
      }
    }
    
    regionGroup.close();
  }

  void SimInputHDF5::ReadNamedNodes(const H5::Group& meshGroup)
  {

    H5::Group namedNodesGroup;
    
    try{
      namedNodesGroup = meshGroup.openGroup( "NamedNodes" );
    } H5_CATCH( "Could not open group for 'NamedNodes'" );

    for( UInt i = 0; i < nodeNames_.size(); i++ ) {

      // read nodes from grid
      std::vector<UInt> nodes;
      H5IO::ReadArray( namedNodesGroup, nodeNames_[i], nodes );

      // add nodes to grid
      mi_->AddNamedNodes( nodeNames_[i], nodes );
    }

    namedNodesGroup.close();
  }

  void SimInputHDF5::ReadNamedElems(const H5::Group& meshGroup)
  {

    H5::Group namedElemGroup;
    
    try{
      namedElemGroup = meshGroup.openGroup( "NamedElems" );
    } H5_CATCH( "Could not open group for 'NamedElems'" );

    for( UInt i = 0; i < elemNames_.size(); i++ ) {

      // read elems from grid
      std::vector<UInt> elems;
      H5IO::ReadArray( namedElemGroup, elemNames_[i], elems );

      // add elems to grid
      mi_->AddNamedElems( elemNames_[i], elems );
    }

    namedElemGroup.close();
  }

  void SimInputHDF5::ReadMeshStats(const H5::Group& meshGroup) {

    // ==================================
    //  Read Region Names and Dimensions
    // ==================================
    H5::Group regionGroup;    
    try {
      regionGroup = meshGroup.openGroup("Regions");
    } H5_CATCH( "Could not open 'Regions' subgroup" );

    regionNames_.clear();
    
    // iterate over all region names
    hsize_t numRegions = regionGroup.getNumObjs();
    for( hsize_t i = 0; i < numRegions; i++ ) {

      // get name
      std::string actName = H5IO::GetObjNameByIdx( regionGroup, i );
      regionNames_.push_back(  actName );

      // get dimension
      UInt dim = 0;
      H5::Group actRegion = regionGroup.openGroup( actName );
      H5IO::ReadAttribute( actRegion, "Dimension", dim );
      regionDims_[actName] = dim;
      actRegion.close();
    }

    regionGroup.close();

    // ==============================
    //  Read Named Nodes Description
    // ==============================
    H5::Group namedNodeGroup;    

    nodeNames_.clear();
    bool namedNodesExist = true;
    try {
      namedNodeGroup = meshGroup.openGroup("NamedNodes");
    } catch (H5::Exception& h5ex) {
      namedNodesExist = false;
      LOG_TRACE(simInputHdf5) << "No named nodes present";
    }
    
    if( namedNodesExist ) {
      
      // iterate over all named nodes' names
      hsize_t numNamedNodeNames = namedNodeGroup.getNumObjs();
      for( hsize_t i = 0; i < numNamedNodeNames; i++ ) {
        
        // get name
        std::string actName = H5IO::GetObjNameByIdx( namedNodeGroup, i );
        nodeNames_.push_back(  actName );
      }
      
      namedNodeGroup.close();
    }

    // ==============================
    //  Read Named Elems Description
    // ==============================
    H5::Group namedElemGroup;    

    elemNames_.clear();
    bool namedElemsExist = true;
    try {
      namedElemGroup = meshGroup.openGroup("NamedElems");
    } catch (H5::Exception& h5ex) {
      namedElemsExist = false;
      LOG_TRACE(simInputHdf5) << "No named elems present";
    }

    if( namedElemsExist ) {
      
      // iterate over all named nodes' names
      hsize_t numNamedElemNames = namedElemGroup.getNumObjs();
      for( hsize_t i = 0; i < numNamedElemNames; i++ ) {
        
        // get name
        std::string actName = H5IO::GetObjNameByIdx( namedElemGroup, i );
        elemNames_.push_back(  actName );
      }
      
      namedElemGroup.close();
    }
    
    statsRead_ = true;
  }

}
