// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/exception.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string.hpp>
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
    capabilities_.insert( SimInput::MESH);
    capabilities_.insert( SimInput::MESH_RESULTS);
    
    mi_ = NULL;
    statsRead_ = false;
    hasExternalFiles_ = false;
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
        readRegions_.Push_back(*it);
    } else {
      readRegions_.Push_back("all");
    }

    // Do not print HDF5 exceptions by default
    H5::Exception::dontPrint();
  }

  SimInputHDF5::~SimInputHDF5() {
    
    mainRoot_.close();
    mainFile_.close();
  }

  void SimInputHDF5::InitModule()
  {
    
    std::string baseName;
    try 
    {
      fs::path fn = fs::system_complete(fileName_);
      fn.normalize();
      baseDir_ = fn.branch_path().native_directory_string();
      baseName = (fs::change_extension( fn.leaf(), "" )).native_directory_string();
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
    LOG_TRACE(simInputHdf5) << "baseDir: " << baseDir_;
    LOG_TRACE(simInputHdf5) << "baseName: " << baseName;

    statsRead_ = false;

    try {
      mainFile_ = H5::H5File( fileName_, H5F_ACC_RDONLY );
    } H5_CATCH( "Could not open HDF5 file '" << fileName_ << "'" );

    try {
      mainRoot_ = mainFile_.openGroup("/");
    } H5_CATCH( "Could not open main root" );
    
    // check for use of external files
    try {
      H5::Group meshResGroup = mainRoot_.openGroup("Results/Mesh");
      H5IO::ReadAttribute( meshResGroup, "ExternalFiles",
                           hasExternalFiles_ );
      meshResGroup.close();
    } catch( H5::Exception& ex ) {
      hasExternalFiles_ = false;
    }
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
     for(UInt i=0; i<readRegions_.GetSize(); i++) {
       if(std::find(regionNames_.Begin(),
                    regionNames_.End(),
                    readRegions_[i]) == regionNames_.End())
         EXCEPTION("Region " << readRegions_[i] << " specified for"
                   " reading does not exist." );
     }
    
     // ========================
     //  READ NODAL INFORMATION
     // ========================

     // get the number of nodes
     try{ 
       nodeGroup = mGroup.openGroup( "Nodes") ;
     } H5_CATCH( "Could not open Elements / Nodes group" );

     //     H5IO::ReadAttribute( nodeGroup, "NumNodes", numNodes_ );

     // read node coordinates
     H5IO::ReadArray( nodeGroup, "Coordinates", nodeCoords_ );

     // read region, element and named entity informaion
     ReadNodeElemData(mGroup);
     ReadNodeGroups(mGroup);
     ReadElemGroups(mGroup);
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
  GetNumMultiSequenceSteps( std::map<UInt,AnalysisType>& analysis,
                            std::map<UInt, UInt>& numSteps,
                            bool isHistory ) {

    H5::Group resultGroup, actMsGroup;
    std::string actAnalysisString;
    AnalysisType actAnalysis;
    analysis.clear();
    numSteps.clear();
    UInt actMsNumSteps = 0;
     
    // try to open grid results: if no groups is present,
    // simply return, as this element is optional.
    try{
      if( !isHistory ) {
        resultGroup = mainRoot_.openGroup("Results").openGroup("Mesh");
      } else {
        resultGroup = mainRoot_.openGroup("Results").openGroup("History");
      }
    } catch (H5::Exception& h5Ex ) {
     return;
    }
    
    // Iterate over all children in the specific group and collect the stepvalues
    hsize_t numChildren = resultGroup.getNumObjs();
    std::set<UInt> msStepNums;
    for( UInt i = 0; i < numChildren; i++ ) {
      std::string actName = H5IO::GetObjNameByIdx( resultGroup, i );
      
      // cut away "MultiStep_"-substring and convert  into integer
      boost::erase_all(actName, "MultiStep_");
      msStepNums.insert( boost::lexical_cast<UInt>(actName) );
    }
    
    // try to find all single multisequence steps and related analysis string
    std::set<UInt>::iterator it;
    for( it = msStepNums.begin(); it != msStepNums.end(); it++ ) {
      actMsGroup = H5IO::GetMultiStepGroup( mainFile_, *it, isHistory );

      // get analyisstring
      H5IO::ReadAttribute( actMsGroup, "AnalysisType", actAnalysisString );
      H5IO::ReadAttribute( actMsGroup, "LastStepNum", actMsNumSteps );
      String2Enum( actAnalysisString, actAnalysis );
      analysis[*it] = actAnalysis;
      numSteps[*it] = actMsNumSteps;
      
      actMsGroup.close();
    }
    resultGroup.close();
  }

  
  void SimInputHDF5::
  GetStepValues( UInt sequenceStep,
                 shared_ptr<ResultInfo> info,
                 std::map<UInt, Double>& steps,
                 bool isHistory ) { 

    
    // open corresponding multistep group
    H5::Group actMsGroup = H5IO::GetMultiStepGroup( mainFile_, 
                                                    sequenceStep, isHistory  );
    
    // open result description
    H5::Group resGroup;
    try {
    resGroup = actMsGroup.openGroup("ResultDescription").
               openGroup( info->resultName );
    } H5_CATCH( "Could not open resultdescription for result '"
                 << info->resultName << "'" );
    
    // read stepValues and stepNumbers
    StdVector<Double> values;
    StdVector<UInt> numbers;
    H5IO::ReadArray( resGroup, "StepNumbers", numbers );
    H5IO::ReadArray( resGroup, "StepValues", values );
    
    // sanity check: both vectors need to have the same dimension
    if( values.GetSize() != numbers.GetSize() ) {
      EXCEPTION( "There are not as many stepnumbers as stepvalues" );
    }
    
    // copy to steps-array
    steps.clear();
    for( UInt i = 0, n=numbers.GetSize(); i < n; i++ ) {
      steps[numbers[i]] = values[i];
    }

    resGroup.close();
    actMsGroup.close();
  }

  
  void SimInputHDF5::
  GetResultTypes( UInt sequenceStep, 
                  StdVector<shared_ptr<ResultInfo> >& infos,
                  bool isHistory ) {
    
    // open ms group and 'Result Description' subgroup
    H5::Group actMsGroup = H5IO::GetMultiStepGroup( mainFile_, sequenceStep,
                                                    isHistory );
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
      StdVector<std::string> dofNames;

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
      }  catch (Exception& ex ) {
        actResultType = NO_SOLUTION_TYPE;
      }
      
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
                     StdVector<shared_ptr<EntityList> >& list,
                     bool isHistory ) {

    // get resultname from resultinfo object
    std::string resultName = info->resultName;
    
    // open ms group and specific entry in 'ResultDescription'
    H5::Group actMsGroup = H5IO::GetMultiStepGroup( mainFile_, sequenceStep,
                           isHistory );
    H5::Group resInfoGroup;
    try {
      resInfoGroup = actMsGroup.openGroup( "ResultDescription/" + resultName );
    } H5_CATCH( "Could not open group result description for result "
                << resultName );
    
    // get regions
    StdVector<std::string> regions;
    H5IO::ReadArray( resInfoGroup, "EntityNames", regions );

    // determine type of list for this result
    EntityList::ListType listType;
    EntityList::DefineType defineType;
    switch( info->definedOn ) {
    case ResultInfo::NODE:
      listType = EntityList::NODE_LIST;
      if( isHistory ) 
        defineType = EntityList::NAMED_NODES;
      else
        defineType = EntityList::REGION;
      break;
    case ResultInfo::ELEMENT:
      listType = EntityList::ELEM_LIST;
      if( isHistory )
        defineType = EntityList::NAMED_ELEMS;
      else
        defineType = EntityList::REGION;
      break;
    case ResultInfo::SURF_ELEM:
      listType = EntityList::SURF_ELEM_LIST;
      defineType = EntityList::REGION;
      break;
    case ResultInfo::REGION:
    case ResultInfo::SURF_REGION:
      listType = EntityList::REGION_LIST;
      defineType = EntityList::REGION;
      break;
    default:
      EXCEPTION( "Only results defined on nodes and elements "
                 << "can be read in from HDF5 file up to now" );
    }

    // iterate over all regions
    list.Clear();
    for( UInt i = 0; i < regions.GetSize(); i++ ) {
      list.Push_back( mi_->GetEntityList( listType, regions[i], 
                                          defineType ) );
    }
    resInfoGroup.close();
    actMsGroup.close();
  }
  
  void SimInputHDF5::GetResult( UInt sequenceStep,
                                UInt stepNum,
                                shared_ptr<BaseResult> result,
                                bool isHistory ) {
    
    if( !isHistory ) {
      GetMeshResult( sequenceStep, stepNum, result);
    } else {
      GetHistResult( sequenceStep, stepNum, result );      
    }
    
  }
  
  void SimInputHDF5::GetMeshResult( UInt sequenceStep, UInt stepNum,
                                     shared_ptr<BaseResult> result ) {
                                     
    // open stepgroup, open specific result subgroup
    H5::Group stepGroup = H5IO::GetStepGroup( mainFile_, sequenceStep, 
                                              stepNum );
    
    // check, if results are stored at external file location
    H5::H5File extFile;
    if( hasExternalFiles_ ) {
      std::string extFileString;
      H5IO::ReadAttribute( stepGroup, "ExtHDF5FileName", extFileString);
      std::string pathsep = fs::path("/").native_directory_string();
      std::string extFileNameComplete = baseDir_ + pathsep + extFileString;
      try {
        extFile = H5::H5File( extFileNameComplete, H5F_ACC_RDONLY );
      } H5_CATCH( "Could not open external file '" 
                  << extFileString << "' for result '" 
                  << result->GetResultInfo()->resultName
                  << "' in multisequence step " << sequenceStep
                  << ", analysis step " << stepNum  );

      // replace old step group by new one
      stepGroup.close();
      stepGroup = extFile.openGroup( "/" );
    }
    
    // determine region for this results
    std::string regionName =  result->GetEntityList()->GetName();

    // determine entity type string
    std::string entString;
    switch( result->GetResultInfo()->definedOn ) {
    case ResultInfo::NODE:
      entString = "Nodes";
      break;
    case ResultInfo::ELEMENT:
    case ResultInfo::SURF_ELEM:
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
                << "' on '" << regionName << "' in step " << stepNum );

    // read data array
    StdVector<Double> realVals;
    H5IO::ReadArray( resGroup, "Real", realVals );

    // copy data array to result object
    if( result->GetEntryType() == EntryType::DOUBLE ) {
      Vector<Double> & resVec = dynamic_cast<Result<Double>& >(*result).GetVector();
      resVec.Resize( realVals.GetSize() );
      for( UInt i = 0; i < realVals.GetSize(); i++ ) {
        resVec[i] = realVals[i];
      } 
    } else {
      Vector<Complex> & resVec = dynamic_cast<Result<Complex>& >(*result).GetVector();
      StdVector<Double> imagVals;
      H5IO::ReadArray( resGroup, "Imag", imagVals );
      
      resVec.Resize( realVals.GetSize() );
      for( UInt i = 0; i < realVals.GetSize(); i++ ) {
        resVec[i] = Complex( realVals[i], imagVals[i] );
      } 
    }
    resGroup.close();
    stepGroup.close();
    // close external file for current step
    if( hasExternalFiles_ )
      extFile.close();
  }
  
  void SimInputHDF5::GetHistResult( UInt sequenceStep, UInt stepNum,
                                         shared_ptr<BaseResult> result ) {

    // open multisequence group
    H5::Group actMsGroup = H5IO::GetMultiStepGroup( mainFile_, sequenceStep,
                                                    true);

    // open group for specific result
    std::string resultName = result->GetResultInfo()->resultName;
    H5::Group actResGroup = actMsGroup.openGroup( resultName );

    // determine from definedOn type the correct string representation
    // of the subgroup
    ResultInfo::EntityUnknownType definedOn = result->GetResultInfo()->definedOn;
    std::string entityTypeString = H5IO::MapUnknownTypeAsString( definedOn );
    H5::Group entityGroup = actResGroup.openGroup( entityTypeString );

    // iterate over all entities in the the list
    shared_ptr<EntityList> list = result->GetEntityList();
    EntityIterator it = list->GetIterator();
    UInt numDofs = result->GetResultInfo()->dofNames.GetSize();
    
    // ----------------------------
    // Case 1: Real valued entries
    // ----------------------------
    if( result->GetEntryType() == EntryType::DOUBLE ) {
      Vector<Double> & resVec = 
        dynamic_cast<Result<Double>&>(*result).GetVector();
      resVec.Resize( list->GetSize() * numDofs );
      for( it.Begin(); !it.IsEnd(); it++ ) {
        // open for each entity the specific subgroup
        H5::Group actEntGroup = 
          entityGroup.openGroup(it.GetIdString() );

        // read single part of array and set it in the result vector
        StdVector<Double> vals;
        H5IO::ReadArray( actEntGroup, "Real", vals );

        for( UInt i = 0; i < numDofs; i++ ) {
          resVec[it.GetPos()*numDofs+i] = vals[(stepNum-1)*numDofs+i];
        }
        actEntGroup.close();
      }
    } 
    // ----------------------------
    // Case 2: Complex valued entries
    // ----------------------------
    else {
      Vector<Complex> & resVec = 
        dynamic_cast<Result<Complex>&>(*result).GetVector();
      resVec.Resize( list->GetSize() * numDofs );
      for( it.Begin(); !it.IsEnd(); it++ ) {
        // open for each entity the specific subgroup
        H5::Group actEntGroup = 
          entityGroup.openGroup(it.GetIdString() );

        // read single part of array and set it in the result vector
        StdVector<Double> realVals, imagVals;
        H5IO::ReadArray( actEntGroup, "Real", realVals );
        H5IO::ReadArray( actEntGroup, "Imag", imagVals );

        for( UInt i = 0; i < numDofs; i++ ) {
          resVec[it.GetPos()*numDofs+i] = 
            Complex( realVals[(stepNum-1)*numDofs+i],
                     imagVals[(stepNum-1)*numDofs+i] );
        }
        actEntGroup.close();
      }
    }
    entityGroup.close();
  }
                                       
  // =========================================================================
  //  MISCELLANEOUS METHODS
  // =========================================================================
  void SimInputHDF5::ReadNodeElemData(const H5::Group& meshGroup)
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
    StdVector<Integer> elemTypes;
    H5IO::ReadArray( elemGroup, "Types", elemTypes );
    
    // read nodes per element
    StdVector<UInt> globConnect;
    H5IO::ReadArray( elemGroup, "Connectivity", globConnect );
    
    elemGroup.close();
    
    typedef std::set<UInt> RegionSetType;
    RegionSetType readNodeSet;
    RegionSetType readElemSet;
    
    // map for each element number the related region
    std::map<UInt, RegionIdType> elemRegionMap;
    
    // ensure, that region names are already read in
    if( !statsRead_ )
      ReadMeshStats( meshGroup );
    
    // ================================
    //  Add Elements Per Element Group
    // ================================
    for( UInt i = 0; i < elemNames_.GetSize(); i++ ) {
      H5::Group entityGroup, actEntityGroup;

      try{
        entityGroup = meshGroup.openGroup( "Groups" );
      } catch (H5::Exception& h5ex) {
        LOG_TRACE(simInputHdf5) << "Could not open group for entity groups";
        return;
      }

      for( UInt i = 0; i < elemNames_.GetSize(); i++ ) {
        // open entitygroup with given name
        try {
          actEntityGroup = entityGroup.openGroup( elemNames_[i] );
        } H5_CATCH( "Could not open definition of element group '"
                    << nodeNames_[i] << "'" );

        // read elems from grid
        StdVector<UInt> readElems;
        H5IO::ReadArray( actEntityGroup, "Elements", readElems );
        for( UInt j = 0; j < readElems.GetSize(); j++ ) {
          readElemSet.insert( readElems[j] );
          elemRegionMap[readElems[j]] = NO_REGION_ID;
        }

        // read nodes from grid
        StdVector<UInt> readNodes;
        H5IO::ReadArray( actEntityGroup, "Nodes", readNodes );
        readNodeSet.insert( readNodes.Begin(), readNodes.End( ));

        actEntityGroup.close();
      }
      entityGroup.close();
    }
    
    // =========================
    //  Add Elements Per Region
    // =========================

    H5::Group regionGroup, actRegion;
    try {
      regionGroup = meshGroup.openGroup( "Regions" );
          } H5_CATCH( "Could not open 'Regions' group" );

    // pass region names to grid and obtain RegionIds
    StdVector<RegionIdType> regionIds;
    mi_->AddRegions(readRegions_, regionIds);

    // Read all nodes from regions and initialize mapping from mesh node
    // numbers to grid node numbers accordingly.
    UInt baseNodeNum = mi_->GetNumNodes() + 1;
    UInt baseElemNum = mi_->GetNumElems() + 1;

    for( UInt i = 0; i < readRegions_.GetSize(); i++ ) {
      try {
        actRegion = regionGroup.openGroup( readRegions_[i] );
      } H5_CATCH( "Could not open group for region '" <<
                  regionNames_[i] << "'" );
      
      RegionIdType actRegionId = regionIds[i];
      
      // read node numbers for this region
      StdVector<UInt> regionNodes;
      H5IO::ReadArray( actRegion, "Nodes", regionNodes );
      
      readNodeSet.insert( regionNodes.Begin(), regionNodes.End());

      // read elem numbers for this region
      StdVector<UInt> regionElems;
      H5IO::ReadArray( actRegion, "Elements", regionElems );
      for( UInt j = 0; j < regionElems.GetSize(); j++ ) {
        readElemSet.insert( regionElems[j] );
        elemRegionMap[regionElems[j]] = actRegionId;
      }
    }
    
    // ===================
    //  Add Nodes to Grid
    // ===================
    mi_->AddNodes( readNodeSet.size() );
    UInt idx;
    Vector<Double> p(3);

    RegionSetType::iterator it, end;
    for( it=readNodeSet.begin(), end=readNodeSet.end();
         it != end;
         it++ ) {

      idx = ((*it)-1)*3;
      p[0] = nodeCoords_[idx + 0];
      p[1] = nodeCoords_[idx + 1];
      p[2] = nodeCoords_[idx + 2];
      mi_->SetNodeCoordinate( baseNodeNum, p );

      nodeNumMap_[*it] = baseNodeNum++;
    }

    // ======================
    //  Add Elements To Grid
    // ======================
    // Add only required elements to grid.
    mi_->AddElems( readElemSet.size() );
    idx = baseElemNum;
    for( it=readElemSet.begin(), end=readElemSet.end();
         it != end;
         it++ ) {
      elemNumMap_[*it] = idx++;
    }

    // Remap global connectivity from mesh nodes to grid node numbers
    for( UInt i = 0, n=globConnect.GetSize();
         i < n;
         i++ ) {
      globConnect[i] = nodeNumMap_[globConnect[i]];
    }

    // iterate over all elements
    for( it = readElemSet.begin(), end = readElemSet.end();
         it != end;
         it++ ) {
      FEType type = (FEType) elemTypes[(*it)-1];
      UInt * connect = &globConnect[numNodesPerElem*((*it)-1)];
      RegionIdType actRegionId = elemRegionMap[*it];
      mi_->SetElemData( elemNumMap_[*it], type, actRegionId, connect );
    }
    
    // ======================================
    //  Generate Node Groups For Each Region
    // ======================================
    
    // check, if nodes nodes of region should be additionally added
    // to list of named nodes
    if( genRegionNodes_) {    

      // iterate over all regions
      for( UInt i = 0; i < readRegions_.GetSize(); i++ ) {

        try {
          actRegion = regionGroup.openGroup( readRegions_[i] );
        } H5_CATCH( "Could not open group for region '" <<
                    regionNames_[i] << "'" );

        // read nodes of region
        StdVector<UInt> regionNodes;
        H5IO::ReadArray( actRegion, "Nodes", regionNodes );

        for( UInt j = 0, n=regionNodes.GetSize();
        j < n;
        j++ ) {
          regionNodes[j] = nodeNumMap_[regionNodes[j]];
        }
        // add nodes as named nodes
        mi_->AddNamedNodes( readRegions_[i]+"_Nodes", regionNodes );
        
        regionGroup.close();
      }
    }

  }

  void SimInputHDF5::ReadNodeGroups(const H5::Group& meshGroup)
  {

    H5::Group entityGroup, actEntityGroup;
    
    try{
      entityGroup = meshGroup.openGroup( "Groups" );
    } catch (H5::Exception& h5ex) {
      LOG_TRACE(simInputHdf5) << "Could not open group for entity groups";
      return;
    }

    for( UInt i = 0; i < nodeNames_.GetSize(); i++ ) {

      // open entitygroup with given name
      try {
        actEntityGroup = entityGroup.openGroup( nodeNames_[i] );
      } H5_CATCH( "Could not open definition of node group '"
                  << nodeNames_[i] << "'" );
      
      // read nodes from file and add to grid
      StdVector<UInt> readNodes;
      StdVector<UInt> nodes;
      H5IO::ReadArray( actEntityGroup, "Nodes", readNodes );

      for( UInt j = 0, n=readNodes.GetSize();
           j < n;
           j++ ) {
        if(nodeNumMap_[readNodes[j]] != 0)
          nodes.Push_back(nodeNumMap_[readNodes[j]]);
      }

      // add nodes to grid
      if(nodes.GetSize())
        mi_->AddNamedNodes( nodeNames_[i], nodes );
      
      actEntityGroup.close();
    }

    entityGroup.close();
  }

  void SimInputHDF5::ReadElemGroups(const H5::Group& meshGroup)
  {
    H5::Group entityGroup, actEntityGroup;
     
     try{
       entityGroup = meshGroup.openGroup( "Groups" );
     } catch (H5::Exception& h5ex) {
       LOG_TRACE(simInputHdf5) << "Could not open group for entity groups";
       return;
     }

     for( UInt i = 0; i < elemNames_.GetSize(); i++ ) {
       // open entitygroup with given name
       try {
         actEntityGroup = entityGroup.openGroup( elemNames_[i] );
       } H5_CATCH( "Could not open definition of element group '"
                   << nodeNames_[i] << "'" );
       
       // read elems from grid
       StdVector<UInt> readElems;
       StdVector<UInt> elems;
       H5IO::ReadArray( actEntityGroup, "Elements", readElems );

       for( UInt j = 0, n=readElems.GetSize();
            j < n;
            j++ ) {
         if(elemNumMap_[readElems[j]] != 0)
           elems.Push_back(elemNumMap_[readElems[j]]);
       }

       actEntityGroup.close();

       // add elems to grid
       if(elems.GetSize()) {
         mi_->AddNamedElems( elemNames_[i], elems );
       }
     }

     entityGroup.close();
  }

  void SimInputHDF5::ReadMeshStats(const H5::Group& meshGroup) {

    // ==================================
    //  Read Region Names and Dimensions
    // ==================================
    H5::Group regionGroup;    
    try {
      regionGroup = meshGroup.openGroup("Regions");
    } H5_CATCH( "Could not open 'Regions' subgroup" );

    regionNames_.Clear();
    
    // iterate over all region names
    hsize_t numRegions = regionGroup.getNumObjs();
    for( hsize_t i = 0; i < numRegions; i++ ) {

      // get name
      std::string actName = H5IO::GetObjNameByIdx( regionGroup, i );
      regionNames_.Push_back(  actName );

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
    H5::Group entityGroup;

    nodeNames_.Clear();

    try {
      entityGroup = meshGroup.openGroup("Groups");
    } catch (H5::Exception& h5ex) {
      LOG_TRACE(simInputHdf5) << "No node / elem groups present";
      statsRead_ = true;
      return;
    }
    
    // =================================
    //  Read node / element group names
    // =================================
    nodeNames_.Clear();
    elemNames_.Clear();
    
    hsize_t numGroups = entityGroup.getNumObjs();
    for( hsize_t i = 0; i < numGroups; i++ ) {
      
      // get name of group
      std::string actName = H5IO::GetObjNameByIdx( entityGroup, i );
      
      // open entitygroup and get number of different entity types
      // (nodes, elements) it is defined on
      H5::Group actEntityGroup;
      try {
        actEntityGroup = entityGroup.openGroup( actName );
      } H5_CATCH( "Could not open entity group '" << actName << "'");
      
      hsize_t numTypes = actEntityGroup.getNumObjs();
      bool hasElems = false;
      for( hsize_t iType = 0; iType < numTypes; iType++ ) {
        std::string actType = H5IO::GetObjNameByIdx( actEntityGroup, iType );
        if( actType == "Elements" ) {
          hasElems = true;  
          break;
        }
      } // loop over types
      
      if( hasElems )  {
        elemNames_.Push_back( actName );
      } else {
        nodeNames_.Push_back( actName );
      }

      actEntityGroup.close();
    } // loop over entity groups
    
    entityGroup.close();
    statsRead_ = true;
  }

}
