// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/exception.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/exception.hh"
#include "Utils/result.hh"
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
    capabilities_.insert( SimInput::MESH);
    capabilities_.insert( SimInput::MESH_RESULTS);
    
    mi_ = NULL;
    statsRead_ = false;
    fileName_ = fileName;
    genRegionNodes_ = false;
    ParamNode *readRegionNode = NULL;
    
    // Change defaults according to XML file
    readRegionNode = myParam_->Get("generateRegionNodes", false);
    if( readRegionNode ) {
      genRegionNodes_ = readRegionNode->AsBool();
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

  void SimInputHDF5::InitModule()
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
    
    

    try {
      mainFile_ = H5::H5File( fileName_, H5F_ACC_RDONLY );
    } H5_CATCH( "Could not open XMDF file '" << fileName_ << "'" );

    try {
      mainRoot_ = mainFile_.openGroup("/");
    } H5_CATCH( "Could not open main root" );
  }

  void SimInputHDF5::ReadMesh( Grid *mi )
  {
    mi_ = mi;
    
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
  //  GENERAL MESH INFORMATION
  // ======================================================
  UInt SimInputHDF5::GetDim() {
    LOG_TRACE(simInputHdf5) << "SimInputHDF5::ReadMesh() not implemented";

    // Open mesh group
    H5::Group meshGroup;
    try{
      meshGroup = mainRoot_.openGroup("Mesh");
    } H5_CATCH( "Could not open mesh group" );
    
    // Read dimension
    UInt dim;
    H5IO::ReadAttribute( meshGroup, "Dimension", dim );
    meshGroup.close();
    return dim;
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

  void SimInputHDF5::GetAllRegionNames( StdVector<std::string> & regionNames ){
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
  //  GENERAL SOLUTION INFORMATION
  // =========================================================================

  void SimInputHDF5::
  GetNumMultiSequenceSteps( StdVector<AnalysisType>& analysis ) {

    H5::Group gridResultGroup, actMsGroup;
    std::string actAnalysisString;
    AnalysisType actAnalysis;
    analysis.Clear();

    // try to open grid results: if no groups is present,
    // simply return, as this element is optional.
    try{
      gridResultGroup = mainRoot_.openGroup("Results").openGroup("Grid");
    } catch (H5::Exception& h5Ex ) {                                        \
     return;
    }

    UInt numMsSteps = gridResultGroup.getNumObjs();
    analysis.Clear();
    analysis.Resize( numMsSteps );

    // try to find all single multisequence steps and related analysis string
    for( UInt i = 0; i < numMsSteps; i++ ) {
      actMsGroup = H5IO::GetMultiStepGroup( mainFile_, i+1 );

      // get analyisstring
      H5IO::ReadAttribute( actMsGroup, "AnalysisType", actAnalysisString );
      String2Enum( actAnalysisString, actAnalysis );
      analysis[i] = actAnalysis;
      
      actMsGroup.close();
    }
    gridResultGroup.close();
  }

  
  void SimInputHDF5::
  GetStepValues( UInt sequenceStep, 
                 StdVector<Double>& stepVals ) {

    // open corresponding multistep group
    H5::Group actMsGroup = H5IO::GetMultiStepGroup( mainFile_, sequenceStep );

    // iterate over all entries
    // NOTE: Since we have a group called "ResultDescription" within the 
    // multisequence group with all the "Step_?" groups, we have to omit this
    // group, if we want to gather all step values
    UInt numSteps = actMsGroup.getNumObjs()-1;
    
    std::set<UInt> indices;
    std::map<UInt, Double> values;
    
    UInt actStep = 0;
    Double actStepVal = 0.0;
    H5::Group actStepGroup;
    for( UInt i = 0; i < numSteps; i++ ) {
      std::string actGroupName = 
        H5IO::GetObjNameByIdx( actMsGroup, i+1 );
      
      // extract index value by cutting away the "Step_" part
      actStep = boost::lexical_cast<UInt>( std::string( actGroupName, 5 ) );
      
      // open step group and query step value
      actStepGroup = actMsGroup.openGroup( actGroupName );
      H5IO::ReadAttribute( actStepGroup, "StepVal", actStepVal );

      // store values
      indices.insert( actStep );
      values[actStep] = actStepVal;
      
    }
    
    // in the end, copy values to vectors
    stepVals.Clear();
    std::set<UInt>::iterator it = indices.begin();
    for( ; it != indices.end(); it++ ) {
      stepVals.Push_back( values[*it] );
    }

    actMsGroup.close();
  }

  
  void SimInputHDF5::
  GetResultTypes( UInt sequenceStep, 
                  StdVector<shared_ptr<ResultInfo> >& infos ) {
    
    // open ms group and 'Result Description' subgroup
    H5::Group actMsGroup = H5IO::GetMultiStepGroup( mainFile_, sequenceStep );
    H5::Group resInfoGroup;
    try {
      resInfoGroup = actMsGroup.openGroup( "ResultDescription" );
    } H5_CATCH( "Could not open group 'ResultDescription'" );

    UInt numResults = resInfoGroup.getNumObjs();
    
    // iterate over all entries and assemble the resultinfo object
    H5::Group actResInfoGroup;
    std::string actResultName;

    infos.Clear();
    for( UInt i = 0; i < numResults; i++ ) {
      actResultName = H5IO::GetObjNameByIdx( resInfoGroup, i );
      try{
        actResInfoGroup = resInfoGroup.openGroup( actResultName );
      } H5_CATCH( "Could not open description group for result '"
                  << actResultName << "'" );

      // Read resultinfo data
      UInt definedOn, numDOFs, entryType;
      std::string unit;
      std::vector<std::string> dofNames;

      H5IO::ReadArray( actResInfoGroup, "DefinedOn", &definedOn );
      H5IO::ReadArray( actResInfoGroup, "NumDOFs", &numDOFs );
      H5IO::ReadArray( actResInfoGroup, "DOFNames", dofNames );
      H5IO::ReadArray( actResInfoGroup, "EntryType", &entryType );
      H5IO::ReadArray( actResInfoGroup, "Unit", &unit );

      // create new ResultInfo objects
      shared_ptr<ResultInfo> ptInfo( new ResultInfo() );
      SolutionType actResultType;
      
      try{ 
        String2Enum( actResultName, actResultType );
      } H5_CATCH( "Could not convert result '" << 
                  actResultName << "' to a SolutionType ");
      ptInfo->resultType = actResultType;
      ptInfo->resultName = actResultName;
      ptInfo->dofNames = StdVector<std::string>(dofNames);
      ptInfo->unit = unit;
      ptInfo->complexFormat = REAL_IMAG;
      ptInfo->entryType = H5IO::MapEntryType( entryType );
      ptInfo->definedOn = H5IO::MapUnknownType( definedOn );

      infos.Push_back( ptInfo );
    }
    
    resInfoGroup.close();
    actMsGroup.close();
  }

  void SimInputHDF5::
  GetResultEntities( UInt sequenceStep,
                     shared_ptr<ResultInfo> info,
                     StdVector<shared_ptr<EntityList> >& list ) {

    // get resultname from resultinfo object
    std::string resultName;
    Enum2String( info->resultType, resultName );
    
    // open ms group and specific entry in 'ResultDescription'
    H5::Group actMsGroup = H5IO::GetMultiStepGroup( mainFile_, sequenceStep );
    H5::Group resInfoGroup;
    try {
      resInfoGroup = actMsGroup.openGroup( "ResultDescription/" + resultName );
    } H5_CATCH( "Could not open group result description for result "
                << resultName );
    
    // get regions
    std::vector<std::string> regions;
    H5IO::ReadArray( resInfoGroup, "Regions", regions );

    // determine type of list for this result
    EntityList::ListType listType;
    switch( info->definedOn ) {
    case ResultInfo::NODE:
      listType = EntityList::NODE_LIST;
      break;
    case ResultInfo::ELEMENT:
      listType = EntityList::ELEM_LIST;
      break;
    default:
      EXCEPTION( "Only results defined on nodes and elements "
                 << "can be read in from HDF5 file up to now" );
    }

    // iterate over all regions
    list.Clear();
    for( UInt i = 0; i < regions.size(); i++ ) {
      list.Push_back( mi_->GetEntityList( listType, regions[i], 
                                          EntityList::REGION ) );
    }
    resInfoGroup.close();
    actMsGroup.close();
  }
  
  void SimInputHDF5::GetResult( UInt sequenceStep,
                                UInt stepValue,
                                shared_ptr<BaseResult> result ) {

    // open stepgroup, open specific result subgroup
    H5::Group stepGroup = H5IO::GetStepGroup( mainFile_, sequenceStep, 
                                              stepValue );


    // determine region for this results
    std::string regionName =  result->GetEntityList()->GetName();

    // determine entity type string
    std::string entString;
    switch( result->GetResultInfo()->definedOn ) {
    case ResultInfo::NODE:
      entString = "Nodes";
      break;
    case ResultInfo::ELEMENT:
      entString = "Elements";
      break;
    default:
      EXCEPTION( "Currently only results on nodes and elements "
                 << "can be read in from a hdf5 file ");
    }

    std::string groupName = result->GetResultInfo()->resultName;
    groupName += "/" + regionName + "/" + entString;
    
    H5::Group resGroup;
    try {
      resGroup = stepGroup.openGroup( groupName );
    } H5_CATCH( "Unable to open group for result '" 
                << result->GetResultInfo()->resultName
                << "' on '" << regionName << "' in step " << stepValue );

    // read data array
    std::vector<Double> realVals;
    H5IO::ReadArray( resGroup, "Real", realVals );

    // copy data array to result object
    if( result->GetEntryType() == EntryType::DOUBLE ) {
      Vector<Double> & resVec = dynamic_cast<Result<Double>& >(*result).GetVector();
      resVec.Resize( realVals.size() );
      for( UInt i = 0; i < realVals.size(); i++ ) {
        resVec[i] = realVals[i];
      } 
    } else {
      Vector<Complex> & resVec = dynamic_cast<Result<Complex>& >(*result).GetVector();
      std::vector<Double> imagVals;
      H5IO::ReadArray( resGroup, "Imag", imagVals );
      
      resVec.Resize( realVals.size() );
      for( UInt i = 0; i < realVals.size(); i++ ) {
        resVec[i] = Complex( realVals[i], imagVals[i] );
      } 
    }
    resGroup.close();
    stepGroup.close();
  }

  // =========================================================================
  //  MISCELLANEOUS METHODS
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
