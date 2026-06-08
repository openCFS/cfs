#include <filesystem>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <hdf5.h>

#include "DataInOut/SimInOut/hdf5/SimInputHDF5.hh"
#include "DataInOut/SimInOut/hdf5/CommonHDF5.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/Exception.hh"
#include "Domain/Results/BaseResults.hh"
#include "Domain/CoordinateSystems/CoordSystem.hh"
#include "Domain/Domain.hh"


namespace fs = std::filesystem;

namespace CoupledField {

  // declare logging stream
  DEFINE_LOG(simInputHdf5, "simInput.hdf5")

  SimInputHDF5::SimInputHDF5(std::string fileName, PtrParamNode inputNode, PtrParamNode infoNode ) :
      SimInput(fileName, inputNode, infoNode )
  {
    capabilities_.insert( SimInput::MESH);
    capabilities_.insert( SimInput::MESH_RESULTS);

    fileName_ = fileName;

    // Change defaults according to XML file
    if(myParam_->Has("generateRegionNodes"))
      genRegionNodes_ = myParam_->Get("generateRegionNodes")->As<bool>();


    boost::char_separator<char> sep(";| ");

    // read, if all entities should be read in
    if(myParam_->Has("readEntities")) 
    {
      std::string readRegions = myParam_->Get("readEntities")->As<std::string>();
      for (const auto& token : boost::tokenizer<boost::char_separator<char>>{readRegions, sep})
        readEntities_.insert(token);
    } 
    if (readEntities_.empty()) {
      readEntities_.insert("__all__");      
      readAllEntities_ = true;
    }
    
    if(myParam_->Has("linearizeEntities")) 
    {
      std::string readRegions =  myParam_->Get("linearizeEntities")->As<std::string>();
      for (const auto& token : boost::tokenizer<boost::char_separator<char>>{readRegions, sep})
        linearizeEntities_.insert(token);
    } 
    
    if(linearizeEntities_.empty()) 
      linearizeEntities_.insert("__none__");
    
    // fetch reference coordinate system
    myParam_->GetValue("coordSysId", coordSysId_, ParamNode::INSERT);
    
    // fetch scale factor
    myParam_->GetValue("scaleFac", scaleFac_, ParamNode::INSERT);

    // this would disable error messages from the library
    // H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
  }

  SimInputHDF5::~SimInputHDF5() 
  {
    // check, if any group is open at all
    if( dbRoot_ > 0 )
      H5Gclose(dbRoot_);

    if( mainRoot_ > 0 )
      H5Gclose(mainRoot_);

    if( mainFile_ > 0 )
      H5Fclose(mainFile_);
  }

  void SimInputHDF5::InitModule()
  {

    // return, if file is already open
    if (mainRoot_ > 0)
      return;
    
    fs::path fn = fs::absolute(fileName_);
    fn = fn.lexically_normal();
    baseDir_ = fn.parent_path().string();
    std::string baseName = fn.filename().replace_extension("").string();
    if(fn.extension() == "")
      fn.replace_extension(".cfs");
    fileName_ = fn.string();

    LOG_TRACE(simInputHdf5) << "fileName_: " << fileName_;
    LOG_TRACE(simInputHdf5) << "baseDir: " << baseDir_;
    LOG_TRACE(simInputHdf5) << "baseName: " << baseName;

    statsRead_ = false;

    if (!fs::exists(fileName_))
      throw Exception("Input file does not exist: " + fileName_);

    hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_file_locking(fapl, false, true); // use_file_locking=false, ignore_when_disabled=true
    mainFile_ = H5Fopen(fileName_.c_str(), H5F_ACC_RDONLY, fapl);
    H5Pclose(fapl);
    if(mainFile_ < 0) 
      throw Exception("Could not open HDF5 file '" + fileName_ +  "' for reading");

    mainRoot_ = H5Gopen2(mainFile_, "/", H5P_DEFAULT);
    if(mainRoot_ < 0) 
      throw Exception("Could not open main root for reading");

    // Read dimension only if present (as it is only written
    // in case a mesh was written before)
    hid_t meshGroup = OpenGroup(mainRoot_, "Mesh");
    if(H5Aexists(meshGroup, "Dimension") > 0) 
      dim_ = ReadAttribute<unsigned int>(meshGroup, ".", "Dimension");
    H5Gclose(meshGroup);

    // check for use of external files
    if(H5Lexists(mainRoot_, "Results/Mesh", H5P_DEFAULT) > 0) 
    {
      hid_t meshResGroup = OpenGroup(mainRoot_, "Results/Mesh");
      hasExternalFiles_ = (ReadAttribute<int>(meshResGroup, ".", "ExternalFiles") != 0);
      H5Gclose(meshResGroup);
    }    
  }

  void SimInputHDF5::ReadMesh( Grid *mi )
  {
    mi_ = mi;

    // Open mesh group
    hid_t mGroup = OpenGroup(mainRoot_, "Mesh");

    // Read infos about mesh
    if(!statsRead_)
      ReadMeshStats(mGroup);

    // If all regions are to be read set list of readRegions accordingly.
    std::set< std::string >::iterator it, end, erase;
    if(*readEntities_.begin() == "__all__") 
    {
      readEntities_.clear();
      readEntities_.insert(regionNames_.Begin(), regionNames_.End());
      readEntities_.insert(nodeNames_.Begin(), nodeNames_.End());
      readEntities_.insert(elemNames_.Begin(), elemNames_.End());
      LOG_DBG(simInputHdf5) << "The following entities will be read:" << ToStringCont(readEntities_);
    }

    // Check if all readEntities_ can be found in file.
    for (const std::string& entity : readEntities_) 
      if (!regionNames_.Contains(entity) && !nodeNames_.Contains(entity) && !elemNames_.Contains(entity))
        throw Exception("Entity " + entity + " specified for reading not found in file.");
    
    // Make sure we have only entities in linearizeEntities_ which are also part of readEntities_
    const std::string& first = *linearizeEntities_.begin();
    if (first == "__none__")
      linearizeEntities_.clear();
    else if (first == "__all__")
      linearizeEntities_ = readEntities_;
    else {
      std::set<std::string> intersection;
      std::set_intersection(linearizeEntities_.begin(), linearizeEntities_.end(),
                            readEntities_.begin(), readEntities_.end(),
                            std::inserter(intersection, intersection.begin()));
      linearizeEntities_ = std::move(intersection);
    }    

    // Remove nodal entities from linearizeEntities_
    for (const auto& name : nodeNames_)
      linearizeEntities_.erase(name); // save for non-existing names


    // ========================
    //  READ NODAL INFORMATION
    // ========================

    hid_t nodeGroup = OpenGroup(mGroup, "Nodes");
    ReadArray(nodeGroup, "Coordinates", nodeCoords_);
    StdVector<unsigned int> dims = GetArrayDims(nodeGroup, "Coordinates");
    // get the number of nodes
    numNodes_ = dims[0];
    H5Gclose(nodeGroup);

    // If a different coordinate system than the default one was specified
    // we map the nodal coordinates into this coordinate system.
    if(coordSysId_ != "default" || scaleFac_ != 1.0) 
    {  
      CoordSystem* coordSys = domain->GetCoordSystem(coordSysId_);
      TransformNodes(*coordSys, scaleFac_);
    }


    // read region, element and named entity information
    ReadNodeElemData(mGroup);
    ReadNodeGroups(mGroup);
    ReadElemGroups(mGroup);

    // Clear / delete any unneeded vector / array
    nodeCoords_.Clear();
    // The mappings file->grid numbers can be deleted in any case, 
    // as it is only needed within the Read...() methods
    f2GElemNumMap_.clear();
    f2GNodeNumMap_.clear();
    H5Gclose(mGroup);
  }



  void SimInputHDF5::GetNumMultiSequenceSteps( std::map<unsigned int,BasePDE::AnalysisType>& analysis,
                            std::map<unsigned int, unsigned int>& numSteps, bool isHistory ) 
  {
    analysis.clear();
    numSteps.clear();

    const char* path = isHistory ? "Results/History" : "Results/Mesh";
    if (H5Lexists(mainRoot_, path, H5P_DEFAULT) <= 0)
        return; // simply return, as this element is optional.
    hid_t resultGroup = OpenGroup(mainRoot_, path);

    // Iterate over all children in the specific group and collect the stepvalues
    hsize_t numChildren = GetInfo(resultGroup).nlinks;
    std::set<unsigned int> msStepNums;
    for( unsigned int i = 0; i < numChildren; i++ ) {
      std::string actName = GetObjNameByIdx( resultGroup, i );

      // cut away "MultiStep_"-substring and convert  into integer
      boost::erase_all(actName, "MultiStep_");
      msStepNums.insert( std::stoul(actName) );
    }

    // try to find all single multisequence steps and related analysis string
    for (unsigned int msStep : msStepNums) {
      hid_t actMsGroup = GetMultiStepGroup(mainFile_, msStep, isHistory);
      analysis[msStep] = BasePDE::analysisType.Parse(ReadAttribute<std::string>(actMsGroup, ".", "AnalysisType"));
      numSteps[msStep] = ReadAttribute<unsigned int>(actMsGroup, ".", "LastStepNum");
      H5Gclose(actMsGroup);
    }    

    H5Gclose(resultGroup);
  }


  void SimInputHDF5::GetStepValues(unsigned int sequenceStep, shared_ptr<ResultInfo> info, std::map<unsigned int, double>& steps, bool isHistory ) 
  {
    // open corresponding multistep group
    hid_t actMsGroup = GetMultiStepGroup( mainFile_, sequenceStep, isHistory );

    // open result description
    hid_t resGroup = OpenGroup(actMsGroup, "ResultDescription/" + info->resultName);

    // read stepValues and stepNumbers
    StdVector<double> values = ReadArray<double>(resGroup, "StepValues");
    StdVector<unsigned int> numbers = ReadArray<unsigned int>(resGroup, "StepNumbers");

    // sanity check: both vectors need to have the same dimension
    if( values.GetSize() != numbers.GetSize() ) 
      throw Exception( "There are not as many step numbers as step values" );
    
    // copy to steps-array
    steps.clear();
    for(unsigned int i = 0; i < numbers.GetSize(); i++)
      steps[numbers[i]] = values[i];

    H5Gclose(resGroup);
    H5Gclose(actMsGroup);
  }

  void SimInputHDF5::GetResultTypes( unsigned int sequenceStep, StdVector<shared_ptr<ResultInfo>>& infos, bool isHistory ) 
  {
    // open ms group and 'Result Description' subgroup
    hid_t actMsGroup = GetMultiStepGroup( mainFile_, sequenceStep, isHistory );
    hid_t resInfoGroup = OpenGroup(actMsGroup, "ResultDescription");

    unsigned int numResults = static_cast<unsigned int>( GetInfo(resInfoGroup).nlinks );

    // iterate over all entries and assemble the resultinfo object
    std::string actResultName;

    infos.Clear();
    for( unsigned int i = 0; i < numResults; i++ ) {
      actResultName = GetObjNameByIdx( resInfoGroup, i );
      hid_t actResInfoGroup = OpenGroup(resInfoGroup, actResultName);

      // Read resultinfo data
      unsigned int definedOn, numDOFs, entryType;
      std::string unit;
      StdVector<std::string> dofNames;

      ReadDataSet( actResInfoGroup, "DefinedOn", &definedOn );
      ReadDataSet( actResInfoGroup, "NumDOFs", &numDOFs );
      ReadArray( actResInfoGroup, "DOFNames", dofNames );
      ReadDataSet( actResInfoGroup, "EntryType", &entryType );
      ReadDataSet( actResInfoGroup, "Unit", &unit );

      // create new ResultInfo objects
      shared_ptr<ResultInfo> ptInfo( new ResultInfo() );
      SolutionType actResultType = SolutionTypeEnum.Parse(actResultName, NO_SOLUTION_TYPE);

      ptInfo->resultType = actResultType;
      ptInfo->resultName = actResultName;
      ptInfo->dofNames = StdVector<std::string>(dofNames);
      ptInfo->unit = unit;
      ptInfo->complexFormat = REAL_IMAG;
      ptInfo->entryType = MapEntryTypeFromInt( entryType );
      ptInfo->definedOn = MapEntityTypeFromInt( definedOn );

      // perform consistency check
      if( ptInfo->entryType == ResultInfo::UNKNOWN ) 
        throw Exception("Result '" + actResultName + "' has no proper EntryType!" );

      if( ptInfo->dofNames.GetSize() == 0 ) 
        throw Exception( "Result '" + actResultName + "' has no degrees of freedoms!");
      
      if(ptInfo->resultName == "" && ptInfo->resultType == NO_SOLUTION_TYPE) 
        throw Exception( "Result has neither a name nor a proper result type!" );
       
      infos.Push_back( ptInfo );
      H5Gclose(actResInfoGroup);
    }

    H5Gclose(resInfoGroup);
    H5Gclose(actMsGroup);
  }

  void SimInputHDF5::GetResultEntities( unsigned int sequenceStep, shared_ptr<ResultInfo> info, StdVector<shared_ptr<EntityList>>& list, bool isHistory ) 
  {
    // get resultname from resultinfo object
    std::string resultName = info->resultName;

    // open ms group and specific entry in 'ResultDescription'
    hid_t actMsGroup = GetMultiStepGroup( mainFile_, sequenceStep, isHistory );
    hid_t resInfoGroup = OpenGroup(actMsGroup, "ResultDescription/" + resultName);

    // get regions
    StdVector<std::string> regions;
    ReadArray( resInfoGroup, "EntityNames", regions );

    // determine type of list for this result
    EntityList::ListType listType;
    switch( info->definedOn ) {
    case ResultInfo::NODE:
    //case ResultInfo::PFEM:
      listType = EntityList::NODE_LIST;
      break;
    case ResultInfo::ELEMENT:
      listType = EntityList::ELEM_LIST;
      break;
    case ResultInfo::SURF_ELEM:
      listType = EntityList::SURF_ELEM_LIST;
      break;
    case ResultInfo::REGION:
      listType = EntityList::REGION_LIST;
      break;
    case ResultInfo::SURF_REGION:
      listType = EntityList::NAME_LIST;
      break;
    default:
      throw Exception( "Only results defined on nodes and elements can be read in from HDF5 file up to now" );
    }

    // iterate over all regions
    list.Clear();
    for( unsigned int i = 0; i < regions.GetSize(); i++ ) 
      list.Push_back( mi_->GetEntityList( listType, regions[i] ) );
    
    H5Gclose(resInfoGroup);
    H5Gclose(actMsGroup);
  }

  void SimInputHDF5::GetResult( unsigned int sequenceStep, unsigned int stepNum, shared_ptr<BaseResult> result, bool isHist) 
  {
    if( !isHist ) 
      GetMeshResult( sequenceStep, stepNum, result);
    else 
      GetHistResult( sequenceStep, stepNum, result );
  }

  void SimInputHDF5::GetMeshResult( unsigned int sequenceStep, unsigned int stepNum, shared_ptr<BaseResult> result ) 
  {
    // open stepgroup, open specific result subgroup
    hid_t stepGroup = GetStepGroup( mainFile_, sequenceStep, stepNum );

    // check, if results are stored at external file location
    hid_t extFile = -1;
    if( hasExternalFiles_ ) {
      std::string extFileString = ReadAttribute<std::string>(stepGroup, ".", "ExtHDF5FileName");
      std::string pathsep = fs::path("/").string();
      std::string extFileNameComplete = baseDir_ + pathsep + extFileString;
      hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
      H5Pset_file_locking(fapl, false, true); // use_file_locking=false, ignore_when_disabled=true
      extFile = H5Fopen(extFileNameComplete.c_str(), H5F_ACC_RDONLY, fapl);
      H5Pclose(fapl);
      if(extFile < 0) EXCEPTION("Could not open external file '" << extFileString
        << "' for result '" << result->GetResultInfo()->resultName
        << "' in multisequence step " << sequenceStep
        << ", analysis step " << stepNum);

      // replace old step group by new one
      H5Gclose(stepGroup);
      stepGroup = H5Gopen(extFile, "/", H5P_DEFAULT);
    }

    // determine region for this results
    std::string regionName =  result->GetEntityList()->GetName();

    // determine entity type string
    std::string entString;
    switch( result->GetResultInfo()->definedOn ) 
    {
      case ResultInfo::NODE:
        entString = "Nodes";
        break;
      case ResultInfo::ELEMENT:
      case ResultInfo::SURF_ELEM:
        entString = "Elements";
        break;
      default:
        throw Exception( "Currently only results on nodes and elements can be read in from a hdf5 file ");
    }

    std::string groupName = result->GetResultInfo()->resultName;
    groupName += "/" + regionName + "/" + entString;

    hid_t resGroup = OpenGroup(stepGroup, groupName);

    // read data array
    StdVector<double> realVals;
    ReadArray( resGroup, "Real", realVals );

    StdVector<unsigned int> idx;
    unsigned int numDofs = result->GetResultInfo()->dofNames.GetSize();
    unsigned int numEntities = result->GetEntityList()->GetSize();
    unsigned int resVecSize =  numEntities * numDofs;

    //std::cout << "resVecSize " << resVecSize << std::endl;

    if(entString == "Nodes") 
      idx = entityNodeMap_[regionName];
    else 
    {
      idx.Resize( numEntities );
      for( unsigned int i = 0; i < numEntities ; i++ )
        idx[i] = i;
    }

    // copy data array to result object
    if(result->GetEntryType() == BaseMatrix::DOUBLE) 
    {
      Vector<double>& resVec = dynamic_cast<Result<double>& >(*result).GetVector();
      resVec.Resize(resVecSize);
      for(unsigned int i = 0; i < numEntities; i++) 
        for(unsigned int iDof = 0; iDof < numDofs; iDof++)  
          resVec[i*numDofs+iDof] = realVals[idx[i]*numDofs+iDof];
    } 
    else 
    {
      Vector<Complex>& resVec = dynamic_cast<Result<Complex>& >(*result).GetVector();
      StdVector<double> imagVals = ReadArray<double>(resGroup, "Imag");

      resVec.Resize(resVecSize);
      for( unsigned int i = 0; i < numEntities; i++ ) 
        for( unsigned int iDof = 0; iDof < numDofs; iDof++ ) 
          resVec[i*numDofs+iDof] = Complex(realVals[idx[i]*numDofs+iDof], imagVals[idx[i]*numDofs+iDof]);
    }
    H5Gclose(resGroup);
    H5Gclose(stepGroup);
    // close external file for current step
    if( hasExternalFiles_ )
      H5Fclose(extFile);
  }

  hid_t SimInputHDF5::OpenHistEntityGroup(hid_t entityGroup,
                                          const std::string& entityTypeString,
                                          const EntityIterator& it) {
    if (entityTypeString == "Nodes") {
      unsigned int fileNodeNum = readAllEntities_ ? it.GetNode() : g2FNodeNumMap_[it.GetNode()];
      return OpenGroup(entityGroup, std::to_string(fileNodeNum));
    } else if (entityTypeString == "Elements") {
      unsigned int fileElemNum = readAllEntities_ ? it.GetElem()->elemNum : g2FElemNumMap_[it.GetElem()->elemNum];
      return OpenGroup(entityGroup, std::to_string(fileElemNum));
    } else {
      return OpenGroup(entityGroup, tempRegionName_ != "" ? tempRegionName_ : it.GetIdString());
    }
  }

  void SimInputHDF5::GetHistResult( unsigned int sequenceStep, unsigned int stepNum, shared_ptr<BaseResult> result) 
  {
    // open multisequence group
    hid_t actMsGroup = GetMultiStepGroup( mainFile_, sequenceStep, true);

    // open group for specific result
    std::string resultName = result->GetResultInfo()->resultName;
    hid_t actResGroup = OpenGroup(actMsGroup, resultName);

    // determine from definedOn type the correct string representation
    // of the subgroup
    ResultInfo::EntityUnknownType definedOn = result->GetResultInfo()->definedOn;
    std::string entityTypeString = entityGroupNameEnum.ToString( definedOn );
    hid_t entityGroup = OpenGroup(actResGroup, entityTypeString);

    // get stepvalues map of result file to extract step indices
    std::map<unsigned int, double> stepVals;
    std::map<shared_ptr<ResultInfo>, std::map<unsigned int, double> > resultSteps;
    shared_ptr<ResultInfo> actRes = result->GetResultInfo();
    this->GetStepValues( sequenceStep, actRes, resultSteps[actRes], true);
    stepVals.insert( resultSteps[actRes].begin(), resultSteps[actRes].end() );
    // Get index of step (because in case of missing steps, stepIdx != stepNum-1)
    unsigned int stepIdx = std::distance(stepVals.begin(),stepVals.find(stepNum));

    // iterate over all entities in the the list
    shared_ptr<EntityList> list = result->GetEntityList();
    EntityIterator it = list->GetIterator();
    unsigned int numDofs = result->GetResultInfo()->dofNames.GetSize();

    // ----------------------------
    // Case 1: Real valued entries
    // ----------------------------
    if( result->GetEntryType() == BaseMatrix::DOUBLE ) 
    {
      Vector<double>& resVec = dynamic_cast<Result<double>&>(*result).GetVector();
      resVec.Resize( list->GetSize() * numDofs );
      for( it.Begin(); !it.IsEnd(); it++ ) {

        hid_t actEntGroup = OpenHistEntityGroup(entityGroup, entityTypeString, it);

        // read single part of array and set it in the result vector
        StdVector<double> vals;
        ReadArray( actEntGroup, "Real", vals );

        for( unsigned int i = 0; i < numDofs; i++ ) {
          resVec[it.GetPos()*numDofs+i] = vals[stepIdx*numDofs+i];
        }
        H5Gclose(actEntGroup);
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
        hid_t actEntGroup = OpenHistEntityGroup(entityGroup, entityTypeString, it);

        // read single part of array and set it in the result vector
        StdVector<double> realVals, imagVals;
        ReadArray( actEntGroup, "Real", realVals );
        ReadArray( actEntGroup, "Imag", imagVals );

        for( unsigned int i = 0; i < numDofs; i++ ) {
          resVec[it.GetPos()*numDofs+i] =
            Complex( realVals[stepIdx*numDofs+i],
                     imagVals[stepIdx*numDofs+i] );
        }
        H5Gclose(actEntGroup);
      }
    }
    H5Gclose(entityGroup);
    H5Gclose(actResGroup);
    H5Gclose(actMsGroup);
  }

  // =========================================================================
  //  MISCELLANEOUS METHODS
  // =========================================================================
  void SimInputHDF5::ReadNodeElemData(hid_t meshGroup)
  {
    // ==================================
    //  Read General Element Information
    // ==================================
    hid_t elemGroup = OpenGroup(meshGroup, "Elements");

    // read number of elements
    unsigned int numElems = ReadAttribute<unsigned int>( elemGroup, ".", "NumElems");

    // read maximum number of nodes per elements
    unsigned int numNodesPerElem = GetArrayDims(elemGroup, "Connectivity")[1];

    // read element types
    StdVector<Integer> elemTypes;
    ReadArray( elemGroup, "Types", elemTypes );

    // read nodes per element
    StdVector<unsigned int> globConnect;
    ReadArray( elemGroup, "Connectivity", globConnect );

    H5Gclose(elemGroup);

    // This sets needs to have defined order - not necessarily sorted. 
    // See comment for SimInputCDB::referencedElems_
    std::set<unsigned int> readNodeSet;
    std::set<unsigned int> readElemSet;

    // map for each element number the related region
    StdVector<RegionIdType> elemRegionMap(numElems);
    elemRegionMap.Init(NO_REGION_ID);

    // ensure, that region names are already read in
    if( !statsRead_ )
      ReadMeshStats( meshGroup );

    // Read all nodes from regions and initialize mapping from mesh node
    // numbers to grid node numbers accordingly.
    unsigned int baseNodeNum = mi_->GetNumNodes() + 1;
    unsigned int baseElemNum = mi_->GetNumElems() + 1;
    if (baseNodeNum > 1 || baseElemNum > 1) {
      readAllEntities_ = false; // We cannot take shortcuts if this is not the only grid
    }

    std::set<std::string>::iterator findIt;

    // ================================
    //  Add Elements Per Element Group
    // ================================
    for( unsigned int i = 0; i < elemNames_.GetSize(); i++ ) 
    {
      hid_t entityGroup = OpenGroup(meshGroup, "Groups");

      for( unsigned int i = 0, n=elemNames_.GetSize(); i < n; i++ ) 
      {
        findIt = readEntities_.find(elemNames_[i]);
        if(findIt == readEntities_.end())
          continue;

        // open entitygroup with given name
        hid_t actEntityGroup = OpenGroup(entityGroup, elemNames_[i]);

        // Check if entity needs to be linearized
        bool linearizeEntity;
        linearizeEntity = linearizeEntities_.find(elemNames_[i]) != linearizeEntities_.end();

        // read elems from grid
        StdVector<unsigned int> readElems;
        ReadArray( actEntityGroup, "Elements", readElems );
        for (unsigned int elemNum : readElems) 
        {
          if(!readAllEntities_)
            readElemSet.insert( elemNum );
          elemRegionMap[elemNum-1] = NO_REGION_ID;
        }

        StdVector<unsigned int> readNodes;
        StdVector<unsigned int> nodeIndices;
        // read nodes from grid
        ReadArray( actEntityGroup, "Nodes", readNodes );
        StdVector<unsigned int> actualNodes(readNodes); // copy

        // Extract nodes needed for linear elements
        if(linearizeEntity) 
        {
          LinearizeElems(readElems, elemTypes, globConnect, actualNodes);
          // To be able to read results for linearized grids we need an index
          // map for each entity (nodes of named elems and regions)
          nodeIndices.Resize(actualNodes.GetSize());
          for( unsigned int j = 0; j < actualNodes.GetSize(); j++ ) 
            nodeIndices[j] = readNodes.Find(actualNodes[j]);
        }
        else 
        {
          // For element lists which need not be linearized we simply use
          // ascending node indices.
          nodeIndices.Resize(actualNodes.GetSize());
          for( unsigned int j = 0; j < actualNodes.GetSize(); j++ ) 
            nodeIndices[j] = j;
        }

        entityNodeMap_[elemNames_[i]] = nodeIndices;
        if( !readAllEntities_)
          readNodeSet.insert( actualNodes.Begin(), actualNodes.End( ));

        H5Gclose(actEntityGroup);
      }
      H5Gclose(entityGroup);
    }

    // =========================
    //  Add Elements Per Region
    // =========================

    hid_t regionGroup = OpenGroup(meshGroup, "Regions");
    hid_t actRegion = -1;

    for( unsigned int i = 0, n=regionNames_.GetSize(); i < n; i++ ) {
      std::string regionName = regionNames_[i];
      findIt = readEntities_.find(regionName);
      if(findIt == readEntities_.end())
        continue;

      actRegion = OpenGroup(regionGroup, regionName);

      // pass region names to grid and obtain RegionIds
      RegionIdType actRegionId = mi_->AddRegion(regionName);

      // Check if entity needs to be linearized
      bool linearizeEntity;
      linearizeEntity = linearizeEntities_.find(regionName) != linearizeEntities_.end();

      // read elem numbers for this region
      StdVector<unsigned int> regionElems;
      ReadArray( actRegion, "Elements", regionElems );
      for( unsigned int elemNum : regionElems ) {
        if( !readAllEntities_ ) 
          readElemSet.insert( elemNum );
        elemRegionMap[elemNum-1] = actRegionId;
      }

      StdVector<unsigned int> regionNodes;
      StdVector<unsigned int> actualNodes;
      StdVector<unsigned int> nodeIndices;

      ReadArray( actRegion, "Nodes", regionNodes );
      actualNodes=regionNodes;

      if(linearizeEntity) 
      {
        LinearizeElems(regionElems, elemTypes, globConnect, actualNodes);
        // To be able to read results for linearized grids we need an index
        // map for each entity (nodes of named elems and regions)
        nodeIndices.Resize(actualNodes.GetSize());
        for( unsigned int j = 0; j < actualNodes.GetSize(); j++ ) 
          nodeIndices[j] = regionNodes.Find(actualNodes[j]);
      }
      else
      {
        // For element list which need not be linearized we simply use
        // ascending node indices.
        nodeIndices.Resize(actualNodes.GetSize());
        for( unsigned int j = 0; j < actualNodes.GetSize(); j++ ) 
          nodeIndices[j] = j;
      }

      entityNodeMap_[regionName] = nodeIndices;
      if( !readAllEntities_ )
        readNodeSet.insert( regionNodes.Begin(), regionNodes.End());
      H5Gclose(actRegion);
    }

    // ===================
    //  Add Nodes to Grid
    // ===================
    
    Vector<double> p(3);

    if( !readAllEntities_ ) {
      numNodes_ = readNodeSet.size();
      mi_->AddNodes( numNodes_ );
      for( unsigned int nodeNum : readNodeSet ) {
        unsigned int idx = (nodeNum-1)*3;
        p[0] = nodeCoords_[idx + 0];
        p[1] = nodeCoords_[idx + 1];
        p[2] = nodeCoords_[idx + 2];
        mi_->SetNodeCoordinate( baseNodeNum, p );

        // perform mapping, if not all entities are read in
        g2FNodeNumMap_[baseNodeNum] = nodeNum;
        f2GNodeNumMap_[nodeNum] = baseNodeNum++;
      }
    } else {
      mi_->AddNodes( numNodes_ );
      for( unsigned int iNode = 0; iNode < numNodes_; ++iNode ) {
        unsigned int idx = iNode * 3;
        p[0] = nodeCoords_[idx + 0];
        p[1] = nodeCoords_[idx + 1];
        p[2] = nodeCoords_[idx + 2];
        mi_->SetNodeCoordinate( iNode+1, p );
      }
    }

    // ======================
    //  Add Elements To Grid
    // ======================
    if( !readAllEntities_ ) {
      // === Add only required elements to grid ===
      mi_->AddElems( readElemSet.size() );
      unsigned int idx = baseElemNum;
      for( unsigned int elemNum : readElemSet ) {
        g2FElemNumMap_[idx] = elemNum;
        f2GElemNumMap_[elemNum] = idx++;
      }

      // Remap global connectivity from mesh nodes to grid node numbers
      for( unsigned int& node : globConnect )
        node = f2GNodeNumMap_[node];

      for( unsigned int elemNum : readElemSet ) {
        Elem::FEType type = (Elem::FEType) elemTypes[elemNum-1];
        unsigned int* connect = &globConnect[numNodesPerElem*(elemNum-1)];
        RegionIdType actRegionId = elemRegionMap[elemNum-1];
        mi_->SetElemData( f2GElemNumMap_[elemNum], type, actRegionId, connect );
      }

    } else {
      // === Add all elements to grid ===
      mi_->AddElems( numElems );
      for( unsigned int iElem = 0; iElem < numElems; ++iElem ) {
        unsigned int elemNum = iElem+1;
        Elem::FEType type = (Elem::FEType) elemTypes[elemNum-1];
        unsigned int* connect = &globConnect[numNodesPerElem*(elemNum-1)];
        RegionIdType actRegionId = elemRegionMap[elemNum-1];
        mi_->SetElemData( elemNum, type, actRegionId, connect );
      }
    }

    // ======================================
    //  Generate Node Groups For Each Region
    // ======================================

    H5Gclose(regionGroup);

    // check, if nodes nodes of region should be additionally added
    // to list of named nodes
    if( genRegionNodes_) {
#if 0
      // iterate over all regions
      for( unsigned int i = 0; i < readEntities_.GetSize(); i++ ) {

        actRegion = OpenGroup(regionGroup, readEntities_[i]);

        // read nodes of region
        StdVector<unsigned int> regionNodes;
        ReadArray( actRegion, "Nodes", regionNodes );

        for( unsigned int j = 0, n=regionNodes.GetSize();
        j < n;
        j++ ) {
          regionNodes[j] = f2GNodeNumMap_[regionNodes[j]];
        }
        // add nodes as named nodes
        mi_->AddNamedNodes( readEntities_[i]+"_Nodes", regionNodes );

        H5Gclose(actRegion);
      }
#else
      EXCEPTION("Generation of nodes per region temporarily not supported!");
#endif
    }

  }

  void SimInputHDF5::ReadNodeGroups(hid_t meshGroup)
  {

    hid_t entityGroup = OpenGroup(meshGroup, "Groups");

    std::set<std::string>::iterator findIt;

    for( unsigned int i = 0; i < nodeNames_.GetSize(); i++ ) {
      findIt = readEntities_.find(nodeNames_[i]);
      if(findIt == readEntities_.end())
        continue;

      // open entitygroup with given name
      hid_t actEntityGroup = OpenGroup(entityGroup, nodeNames_[i]);

      // read nodes from file and add to grid
      StdVector<unsigned int> readNodes;
      StdVector<unsigned int> nodes;
      ReadArray( actEntityGroup, "Nodes", readNodes );

      for( unsigned int nodeNum : readNodes ) {
        if( !readAllEntities_ ) {
          if( f2GNodeNumMap_[nodeNum] != 0 )
            nodes.Push_back( f2GNodeNumMap_[nodeNum] );
        } else {
          nodes.Push_back( nodeNum );
        }
      }

      // add nodes to grid
      if(nodes.GetSize())
        mi_->AddNamedNodes( nodeNames_[i], nodes );

      H5Gclose(actEntityGroup);
    }

    H5Gclose(entityGroup);
  }

  void SimInputHDF5::ReadElemGroups(hid_t meshGroup)
  {
    hid_t entityGroup = OpenGroup(meshGroup, "Groups");

     std::set<std::string>::iterator findIt;

     for( unsigned int i = 0; i < elemNames_.GetSize(); i++ ) {
       findIt = readEntities_.find(elemNames_[i]);
       if(findIt == readEntities_.end())
         continue;

       // open entitygroup with given name
       hid_t actEntityGroup = OpenGroup(entityGroup, elemNames_[i]);

       // read elems from grid
       StdVector<unsigned int> readElems;
       StdVector<unsigned int> elems;
       ReadArray( actEntityGroup, "Elements", readElems );

       for( unsigned int elemNum : readElems ) 
       {
         if( !readAllEntities_ ) {
           if( f2GElemNumMap_[elemNum] != 0 )
             elems.Push_back( f2GElemNumMap_[elemNum] );
         } 
         else 
           elems.Push_back( elemNum );
       }

       H5Gclose(actEntityGroup);

       // add elems to grid
       if(elems.GetSize()) 
         mi_->AddNamedElems( elemNames_[i], elems );
     }

     H5Gclose(entityGroup);
  }

  void SimInputHDF5::ReadMeshStats(hid_t meshGroup) {

    // ==================================
    //  Read Region Names and Dimensions
    // ==================================
    hid_t regionGroup = OpenGroup(meshGroup, "Regions");

    regionNames_.Clear();

    // iterate over all region names
    hsize_t numRegions = GetInfo(regionGroup).nlinks;
    for( hsize_t i = 0; i < numRegions; i++ ) {
      // get name
      std::string actName = GetObjNameByIdx( regionGroup, i );
      regionNames_.Push_back(  actName );
    }

    H5Gclose(regionGroup);

    // ==============================
    //  Read Named Nodes Description
    // ==============================
    if( H5Lexists(meshGroup, "Groups", H5P_DEFAULT) <= 0 ) {
      LOG_TRACE(simInputHdf5) << "Could not open group for entity groups";
      statsRead_ = true;
      return;
    }
    hid_t entityGroup = OpenGroup(meshGroup, "Groups");

    // =================================
    //  Read node / element group names
    // =================================
    nodeNames_.Clear();
    elemNames_.Clear();

    hsize_t numGroups = GetInfo(entityGroup).nlinks;
    for( hsize_t i = 0; i < numGroups; i++ ) 
    {
      // get name of group
      std::string actName = GetObjNameByIdx( entityGroup, i );

      // open entitygroup and get number of different entity types
      // (nodes, elements) it is defined on
      hid_t actEntityGroup = OpenGroup(entityGroup, actName);

      hsize_t numTypes = GetInfo(actEntityGroup).nlinks;
      bool hasElems = false;
      for( hsize_t iType = 0; iType < numTypes; iType++ ) {
        std::string actType = GetObjNameByIdx( actEntityGroup, iType );
        if( actType == "Elements" ) {
          hasElems = true;
          break;
        }
      } // loop over types

      if( hasElems )  
        elemNames_.Push_back( actName );
      else 
        nodeNames_.Push_back( actName );

      H5Gclose(actEntityGroup);
    } // loop over entity groups

    H5Gclose(entityGroup);
    statsRead_ = true;
  }

  void SimInputHDF5::LinearizeElems(const StdVector<unsigned int>& readElems,
                                    StdVector<Integer>& elemTypes,
                                    StdVector<unsigned int>& globConnect,
                                    StdVector<unsigned int>& readNodes) {
    static std::map<Elem::FEType, Elem::FEType> elemTypeMap;
    unsigned int elemIncr = globConnect.GetSize() / elemTypes.GetSize();
    unsigned int numElems = readElems.GetSize();
    std::set<unsigned int> readNodeSet;

    if(!elemTypeMap.size()) {
      elemTypeMap[Elem::ET_LINE2]   = Elem::ET_LINE2;
      elemTypeMap[Elem::ET_LINE3]   = Elem::ET_LINE2;
      elemTypeMap[Elem::ET_TRIA3]   = Elem::ET_TRIA3;
      elemTypeMap[Elem::ET_TRIA6]   = Elem::ET_TRIA3;
      elemTypeMap[Elem::ET_QUAD4]   = Elem::ET_QUAD4;
      elemTypeMap[Elem::ET_QUAD8]   = Elem::ET_QUAD4;
      elemTypeMap[Elem::ET_QUAD9]   = Elem::ET_QUAD4;
      elemTypeMap[Elem::ET_TET4]    = Elem::ET_TET4;
      elemTypeMap[Elem::ET_TET10]   = Elem::ET_TET4;
      elemTypeMap[Elem::ET_HEXA8]   = Elem::ET_HEXA8;
      elemTypeMap[Elem::ET_HEXA20]  = Elem::ET_HEXA8;
      elemTypeMap[Elem::ET_HEXA27]  = Elem::ET_HEXA8;
      elemTypeMap[Elem::ET_PYRA5]   = Elem::ET_PYRA5;
      elemTypeMap[Elem::ET_PYRA13]  = Elem::ET_PYRA5;
      elemTypeMap[Elem::ET_WEDGE6]  = Elem::ET_WEDGE6;
      elemTypeMap[Elem::ET_WEDGE15] = Elem::ET_WEDGE6;
    }

    for(unsigned int i=0; i<numElems; i++) {
      unsigned int elemIdx=readElems[i]-1;
      Elem::FEType newType = elemTypeMap[(Elem::FEType)elemTypes[elemIdx]];
      unsigned int numNodes = Elem::shapes[newType].numNodes;

      elemTypes[elemIdx] = newType;
      readNodeSet.insert(&globConnect[elemIdx*elemIncr],
                         &globConnect[elemIdx*elemIncr+numNodes]);
    }

    readNodes.Resize(readNodeSet.size());
    std::copy(readNodeSet.begin(), readNodeSet.end(), readNodes.Begin());
  }

  void SimInputHDF5::TransformNodes(CoordSystem& coordSys, double scaleFac)
  {
//    if (dim_ != coordSys.GetDim()) {
//      EXCEPTION("Cannot use a " << coordSys.GetDim() << "D coordinate system ("
//                << coordSys.GetName() << ") to transform a "
//                << dim_ << "D mesh (" << fileName_ << ").");
//    }
    
    Vector<double> p, globPoint;
    p.Resize(dim_);
    globPoint.Resize(dim_);

    for(unsigned int i=0; i<numNodes_; ++i) 
    {
      unsigned int idx = i*3;
      p[0] = nodeCoords_[idx + 0];
      p[1] = nodeCoords_[idx + 1];
      if (dim_ == 3) 
        p[2] = nodeCoords_[idx + 2];
      coordSys.Global2LocalCoord(globPoint, p);
      
      nodeCoords_[idx + 0] = globPoint[0] * scaleFac;
      nodeCoords_[idx + 1] = globPoint[1] * scaleFac;
      if (dim_ == 3) 
        nodeCoords_[idx + 2] = globPoint[2] * scaleFac;
    }
  }
  
  void SimInputHDF5::ReadStringFromUserData(const std::string& dSetName, std::string& str) 
  {
    hid_t userDataGroup = OpenGroup(mainRoot_, "UserData");

    StdVector<std::string> userData;
    ReadArray( userDataGroup, dSetName, userData );
    str = userData[0];

    H5Gclose(userDataGroup);
  }


  // ------------------------------------------------------------------------
  //  DATABASE SECTION
  // ------------------------------------------------------------------------
  
  void SimInputHDF5::DB_Init() {
    if(dbRoot_ > 0) 
      H5Gclose(dbRoot_); // don't orphan a previously-opened handle on re-init
    dbRoot_ = OpenGroup(mainRoot_, "DataBase");
  }

  void SimInputHDF5::DB_GetParamFileContent( std::string& params ) {
    hid_t extFiles = OpenGroup(dbRoot_, "InputFiles");

    StdVector<std::string> vec;
    ReadArray( extFiles, "ParameterFile", vec);
    params = vec[0];
    H5Gclose(extFiles);
  }

  void SimInputHDF5::DB_GetMatFileContent( std::string& params ) {
    hid_t extFiles = OpenGroup(dbRoot_, "InputFiles");

    StdVector<std::string> vec;
    ReadArray( extFiles, "MaterialFile", vec);
    params = vec[0];
    H5Gclose(extFiles);
  }

  void SimInputHDF5::DB_GetFeFctCoefs( unsigned int sequenceStep, unsigned int stepNum,
                                       const std::string& pdeName,
                                       const std::string& fctName,
                                       SingleVector * coefs ) {
    
    
    hid_t msGroup = OpenGroup(dbRoot_, "MultiSteps");
    hid_t mGroup  = OpenGroup(msGroup,  std::to_string(sequenceStep));
    hid_t physGroup = OpenGroup(mGroup, pdeName);
    hid_t fctGroup  = OpenGroup(physGroup, fctName);
    std::string stepStr = std::to_string(stepNum);
    hid_t stepGroup = OpenGroup(fctGroup, stepStr);

    // Now aquire coefficients for fefunction
    if( coefs->GetEntryType() == BaseMatrix::DOUBLE ) {
      Vector<double> & rVec = dynamic_cast<Vector<double>& >(*coefs);
      
      // read data array
      StdVector<double> realVals;
      ReadArray( stepGroup, "Real", realVals );
      rVec.Resize( realVals.GetSize() );
      for( unsigned int i = 0; i < realVals.GetSize(); ++i ) {
        rVec[i] = realVals[i];
      }
      
    } else {
      Vector<Complex> & cVec = dynamic_cast<Vector<Complex>& > (*coefs);
      StdVector<double> realVals, imagVals;
      ReadArray( stepGroup, "Real", realVals );
      ReadArray( stepGroup, "Imag", imagVals );
      cVec.Resize( realVals.GetSize() );
      for( unsigned int i = 0; i < realVals.GetSize(); ++i ) {
        cVec[i] = Complex(realVals[i], imagVals[i] );
      }
    }

    H5Gclose(stepGroup);
    H5Gclose(fctGroup);
    H5Gclose(physGroup);
    H5Gclose(mGroup);
    H5Gclose(msGroup);
  }

  void SimInputHDF5::
  DB_GetNumMultiSequenceSteps( std::map<unsigned int, BasePDE::AnalysisType>& analysis,
                               std::map<unsigned int, double>& accTime,
                               std::map<unsigned int, bool>& isFinished ){

    // Initialize DB
    DB_Init();
    
    // Open "MultiStep" section 
    hid_t msGroup = OpenGroup(dbRoot_, "MultiSteps");
    
    // Iterate over all children
    hsize_t numChildren = GetInfo(msGroup).nlinks;
    for( unsigned int i = 0; i < numChildren; ++i ) {
      unsigned int actMsStep = std::stoul(GetObjNameByIdx( msGroup, i ));
      
      // open group and read analysyis type
      hid_t actMsGroup = OpenGroup(msGroup, std::to_string(actMsStep));

      // obtain type of analysis
      std::string analysisString;
      ReadAttribute( actMsGroup, ".", "AnalysisType", analysisString );
      analysis[actMsStep] = BasePDE::analysisType.Parse(analysisString);
      int completedInt = ReadAttribute<int>( actMsGroup, ".", "Completed");
      isFinished[actMsStep] = (completedInt != 0);
      ReadAttribute( actMsGroup, ".", "AccTime", accTime[actMsStep] );
      
      H5Gclose(actMsGroup);
    } // loop: multisequence steps
    
    H5Gclose(msGroup);
  }
  
  void SimInputHDF5::
  DB_GetAvailPdeCoefFcts( unsigned int msStep, 
                          std::map<std::string, 
                          std::set<std::string> >& coefFcts ) {
    hid_t msGroup   = OpenGroup(dbRoot_, "MultiSteps");
    hid_t mGroup    = OpenGroup(msGroup, std::to_string(msStep));

    // Loop over all physic types
    hsize_t numPhysic = GetInfo(mGroup).nlinks;
    for( unsigned int iPhysic = 0; iPhysic < numPhysic; ++iPhysic ) {
      std::string physic = GetObjNameByIdx( mGroup, iPhysic );
      hid_t pGroup = OpenGroup(mGroup, physic);

      // Loop over all coefficient types
      hsize_t numCoefs = GetInfo(pGroup).nlinks;
      for( unsigned int iCoef = 0; iCoef < numCoefs; ++iCoef ) {
        std::string coef = GetObjNameByIdx( pGroup, iCoef );
        coefFcts[physic].insert(coef);
      } // loop: coefs
      
      H5Gclose(pGroup);
    } // loop: physic types
    H5Gclose(mGroup);
    H5Gclose(msGroup);
  }

  void SimInputHDF5::
  DB_GetStepValues( unsigned int sequenceStep,
                    const std::string& pdeName,
                    const std::string& fctName,
                    std::map<unsigned int, double>& stepValues ) {

    hid_t msGroup   = OpenGroup(dbRoot_, "MultiSteps");
    hid_t mGroup    = OpenGroup(msGroup, std::to_string(sequenceStep));
    hid_t physGroup = OpenGroup(mGroup, pdeName);
    hid_t fctGroup  = OpenGroup(physGroup, fctName);

    // Loop over all time steps
    hsize_t numSteps = GetInfo(fctGroup).nlinks;
    for( unsigned int iStep = 0; iStep < numSteps; ++iStep ) {
      unsigned int stepNum = std::stoul(GetObjNameByIdx( fctGroup, iStep ));
      hid_t stepGroup = OpenGroup(fctGroup, std::to_string(stepNum));

      // get time step value
      double stepVal = ReadAttribute<double>( stepGroup, ".", "StepValue");
      stepValues[stepNum] = stepVal;
      H5Gclose(stepGroup);
    } //loop: steps

    H5Gclose(fctGroup);
    H5Gclose(physGroup);
    H5Gclose(mGroup);
    H5Gclose(msGroup);

  }

  void SimInputHDF5::DB_Close() 
  {
    if(dbRoot_ > 0) {
      H5Gclose(dbRoot_);
      dbRoot_ = -1;
    }
  }

  void SimInputHDF5::GetNamedNodeResult(const std::string& nodeName,
                                        const std::string resultName,
                                        StdVector<Complex>& result)
  {
    // Read node numbers for the named group
    if( H5Lexists(mainRoot_, ("Mesh/Groups/" + nodeName).c_str(), H5P_DEFAULT) <= 0 )
      EXCEPTION("Named node '" << nodeName << "' not found in HDF5 file!");

    hid_t groupsGroup = OpenGroup(mainRoot_, "Mesh/Groups");
    hid_t namedNodeGroup = OpenGroup(groupsGroup, nodeName);
    StdVector<unsigned int> nodeNumbers;
    ReadArray(namedNodeGroup, "Nodes", nodeNumbers);
    H5Gclose(namedNodeGroup);
    H5Gclose(groupsGroup);

    for( unsigned int i = 0, n = nodeNumbers.GetSize(); i < n; i++ )
      std::cout << "Node numbers for " << nodeName << ": " << nodeNumbers[i] << std::endl;

    // Try history result first, fall back to mesh result
    std::string histPath = "Results/History/MultiStep_1/" + resultName
                           + "/Nodes/" + std::to_string(nodeNumbers[0]);
    if( H5Lexists(mainRoot_, histPath.c_str(), H5P_DEFAULT) > 0 ) {
      hid_t historyNodeGroup = OpenGroup(mainRoot_, histPath);
      StdVector<double> real, imag;
      ReadArray(historyNodeGroup, "Real", real);
      ReadArray(historyNodeGroup, "Imag", imag);
      H5Gclose(historyNodeGroup);

      result.Resize(1);
      result[0] = Complex(real[0], imag[0]);
      std::cout << "Complex result " << result[0] << std::endl;
    } else {
      StdVector<shared_ptr<ResultInfo>> infos;
      shared_ptr<ResultInfo> myInfo;
      GetResultTypes(1, infos, false);
      for( unsigned int i = 0, n = infos.GetSize(); i < n; i++ ) {
        if( infos[i]->resultName == resultName ) { myInfo = infos[i]; break; }
      }

      StdVector<shared_ptr<EntityList>> list;
      shared_ptr<EntityList> myList;
      GetResultEntities(1, myInfo, list, false);

      unsigned int dist = 0;
      for( unsigned int i = 0, n = list.GetSize(); i < n; i++ ) {
        const StdVector<unsigned int>& nodes = dynamic_cast<NodeList*>(list[i].get())->GetNodes();
        std::string name = dynamic_cast<NodeList*>(list[i].get())->GetName();
        const unsigned int* pt = std::find(&nodes[0], &nodes[nodes.GetSize()-1]+1, nodeNumbers[0]);
        dist = std::distance(&nodes[0], pt);
        std::cout << "Searching for node " << nodeNumbers[0] << " in list " << name << std::endl;
        if( dist < nodes.GetSize() ) { myList = list[i]; std::cout << "Found node!" << std::endl; break; }
      }

      Result<Complex>* resPt = new Result<Complex>;
      shared_ptr<BaseResult> res(resPt);
      res->SetResultInfo(myInfo);
      res->SetEntityList(myList);
      GetMeshResult(1, 1, res);

      Vector<Complex> resVec = resPt->GetVector();
      result.Resize(1);
      result[0] = resVec[dist * myInfo->dofNames.GetSize()];
      std::cout << "Complex result " << result[0] << std::endl;
    }
  }

}
