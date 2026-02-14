#include <filesystem>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/Exception.hh"
#include "Domain/Results/BaseResults.hh"
#include "Domain/CoordinateSystems/CoordSystem.hh"
#include "Domain/Domain.hh"
#include "SimInputHDF5.hh"
#include "hdf5io.hh"

namespace fs = std::filesystem;
using boost::char_separator;

namespace CoupledField {

  // declare logging stream
  DEFINE_LOG(simInputHdf5, "simInput.hdf5")

#define H5_CATCH( STR )                                                 \
  catch (H5::Exception& h5Ex ) {                                        \
    EXCEPTION( STR << ":\n" << h5Ex.getCDetailMsg() );                  \
  }

  SimInputHDF5::SimInputHDF5(std::string fileName, PtrParamNode inputNode,
                             PtrParamNode infoNode ) :
      SimInput(fileName, inputNode, infoNode )
  {
    capabilities_.insert( SimInput::MESH);
    capabilities_.insert( SimInput::MESH_RESULTS);

    mi_ = NULL;
    statsRead_ = false;
    hasExternalFiles_ = false;
    fileName_ = fileName;
    genRegionNodes_ = false;
    readAllEntities_ = true;
    PtrParamNode pNode;
    coordSysId_ = "default";
    scaleFac_ = 1.0;
    numNodes_ = 0;

    // Change defaults according to XML file
    pNode = myParam_->Get("generateRegionNodes", ParamNode::PASS);
    if( pNode ) {
      genRegionNodes_ = pNode->As<bool>();
    }

    // read, if all entities should be read in
    pNode = myParam_->Get("readEntities", ParamNode::PASS);
    if(pNode) {
      std::string readRegions =
        myParam_->Get("readEntities")->As<std::string>();

      typedef boost::tokenizer<boost::char_separator<char> > Tok;
      boost::char_separator<char> sep(";| ");
      Tok t(readRegions, sep);
      Tok::iterator it, end;
      it = t.begin();
      end = t.end();

      for( ; it != end; it++)
        readEntities_.insert(*it);

      // final check: maybe there is no real region
      if( readEntities_.size() == 0 ) {
        readEntities_.insert("__all__");
      }
    } else {
      readEntities_.insert("__all__");
    }
    
    if( *(readEntities_.begin()) == "__all__")
      readAllEntities_ = true;
    else 
      readAllEntities_ = false;

    pNode = myParam_->Get("linearizeEntities", ParamNode::PASS);
    if(pNode) {
      std::string readRegions =
        myParam_->Get("linearizeEntities")->As<std::string>();

      typedef boost::tokenizer<char_separator<char> > Tok;
      boost::char_separator<char> sep(";| ");
      Tok t(readRegions, sep);
      Tok::iterator it, end;
      it = t.begin();
      end = t.end();

      for( ; it != end; it++)
        linearizeEntities_.insert(*it);
      if( linearizeEntities_.size() == 0 ) {
        linearizeEntities_.insert("__none__");
      }
    } else {
      linearizeEntities_.insert("__none__");
    }
    // fetch reference coordinate system
    myParam_->GetValue("coordSysId", coordSysId_, ParamNode::INSERT);
    
    // fetch scale factor
    myParam_->GetValue("scaleFac", scaleFac_, ParamNode::INSERT);

    // Do not print HDF5 exceptions by default
    H5::Exception::dontPrint();
  }

  SimInputHDF5::~SimInputHDF5() {
    
    // check, if any group is open at all
    if( dbRoot_.getLocId() > 0 )
      dbRoot_.close();

    if( mainRoot_.getLocId() > 0 )
      mainRoot_.close();

    if( mainFile_.getLocId() <= 0 )
      return;


    mainRoot_.close();
    mainFile_.close();
  }

  void SimInputHDF5::InitModule()
  {

    // return, if file is already open
    if (mainRoot_.getLocId() > 0)
      return;
    
    std::string baseName;
    try
    {
      fs::path fn = fileName_.empty() ? fs::current_path() : fs::absolute(fileName_);
      fn = fn.lexically_normal();
      baseDir_ = fn.parent_path().string();
      baseName = fn.filename().replace_extension("").string();
      if(fn.extension() == "")
      {
        fn = fn.replace_extension(".h5");
      }
      fileName_ = fn.string();
    } catch (std::filesystem::filesystem_error& ex)
    {
      EXCEPTION("Received exception: " << ex.what());
      return;
    }

    LOG_TRACE(simInputHdf5) << "fileName_: " << fileName_;
    LOG_TRACE(simInputHdf5) << "baseDir: " << baseDir_;
    LOG_TRACE(simInputHdf5) << "baseName: " << baseName;

    statsRead_ = false;

    if ( !fs::exists( fileName_ ))
    {
      EXCEPTION("Input file does not exist:\n" << fileName_ << std::endl);
    }
    try {
      mainFile_ = H5::H5File( fileName_, H5F_ACC_RDONLY );
    } H5_CATCH( "Could not open HDF5 file '" << fileName_ << "'" );

    try {
      mainRoot_ = mainFile_.openGroup("/");
    } H5_CATCH( "Could not open main root" );

    // Open mesh group
    H5::Group meshGroup;
    try{
      meshGroup = mainRoot_.openGroup("Mesh");
    } H5_CATCH( "Could not open mesh group" );

    // Read dimension only if present (as it is only written
    // in case a mesh was written before)
    if(H5IO::AttrExists(meshGroup, "Dimension")) {
      H5IO::ReadAttribute( meshGroup, "Dimension", dim_ );
    }  
    meshGroup.close();

    // check for use of external files
    try {
      H5::Group meshResGroup = mainRoot_.openGroup("Results/Mesh");
      H5IO::ReadAttribute( meshResGroup, "ExternalFiles",
                           hasExternalFiles_ );
      meshResGroup.close();
    } catch( H5::Exception&) {
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
    std::set< std::string >::iterator it, end, erase;
    if(*readEntities_.begin() == "__all__") {
      readEntities_.clear();
      readEntities_.insert(regionNames_.Begin(), regionNames_.End());
      readEntities_.insert(nodeNames_.Begin(), nodeNames_.End());
      readEntities_.insert(elemNames_.Begin(), elemNames_.End());
      LOG_DBG(simInputHdf5) << "The following entities will be read:";
      std::stringstream sstr;
      it=readEntities_.begin();
      end=readEntities_.end();
      for( ; it != end; ++it)
        sstr << (*it) << " ";
      LOG_DBG(simInputHdf5) << sstr.str();
    }

    // Check if all readEntities_ can be found in file.
    it=readEntities_.begin();
    end=readEntities_.end();
    for( ; it != end; ++it) {
      StdVector<std::string>::iterator findIt;

      findIt = std::find(regionNames_.Begin(), regionNames_.End(), *it);
      if( findIt != regionNames_.End()) continue;

      findIt = std::find(nodeNames_.Begin(), nodeNames_.End(), *it);
      if( findIt != nodeNames_.End()) continue;

      findIt = std::find(elemNames_.Begin(), elemNames_.End(), *it);
      if( findIt != elemNames_.End()) continue;

      EXCEPTION("Entity " << (*it) << " specified for"
                " reading does not exist." );
   }

    // Make sure we have only entities in linearizeEntities_ which are
    // also part of readEntities_
    if(!linearizeEntities_.empty()) 
    {      
      if ( *linearizeEntities_.begin() == "__none__" ) {
        linearizeEntities_.clear();
      } else if ( *linearizeEntities_.begin() == "__all__" ) {
        linearizeEntities_.insert(readEntities_.begin(), readEntities_.end());
      } else {
        it = linearizeEntities_.begin();
        end = linearizeEntities_.end();
        for ( ; it != end; ) {
          if ( readEntities_.find(*it) == readEntities_.end() ) {
            linearizeEntities_.erase(it++);
          } else {
            ++it;
          }
        }
      }
    }    

    // Remove nodal entities from linearizeEntities_
    StdVector<std::string>::iterator findIt, endIt = nodeNames_.End();
    it = linearizeEntities_.begin();
    end = linearizeEntities_.end();
    for ( ; it != end; ) {
      findIt = std::find(nodeNames_.Begin(), nodeNames_.End(), *it);
      if ( findIt != endIt ) {
        linearizeEntities_.erase(it++);
      } else {
        ++it;
      }
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
    StdVector<UInt> dims;
    dims = H5IO::GetArrayDims(nodeGroup, "Coordinates");
    numNodes_ = dims[0]; 

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
  }


  // ======================================================
  //  GENERAL MESH INFORMATION
  // ======================================================
  UInt SimInputHDF5::GetDim() {
    return dim_;
  }

  // =========================================================================
  //  GENERAL SOLUTION INFORMATION
  // =========================================================================

  void SimInputHDF5::
  GetNumMultiSequenceSteps( std::map<UInt,BasePDE::AnalysisType>& analysis,
                            std::map<UInt, UInt>& numSteps,
                            bool isHistory ) {

    H5::Group resultGroup, actMsGroup;
    std::string actAnalysisString;
    BasePDE::AnalysisType actAnalysis;
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
    } catch (H5::Exception&) {
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
      actAnalysis = BasePDE::analysisType.Parse(actAnalysisString);
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
    } H5_CATCH( "Could not open result description for result '"
                 << info->resultName << "'" );

    // read stepValues and stepNumbers
    StdVector<Double> values;
    StdVector<UInt> numbers;
    try {
      H5IO::ReadArray( resGroup, "StepNumbers", numbers );
    } H5_CATCH( "Could not read step numbers of result '"
                << info->resultName << "'" );
    try {
      H5IO::ReadArray( resGroup, "StepValues", values );
    } H5_CATCH( "Could not read step values of result '"
                << info->resultName << "'" );

    // sanity check: both vectors need to have the same dimension
    if( values.GetSize() != numbers.GetSize() ) {
      EXCEPTION( "There are not as many step numbers as step values" );
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

    UInt numResults = static_cast<UInt>( resInfoGroup.getNumObjs() );

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
      SolutionType actResultType = SolutionTypeEnum.Parse(actResultName, NO_SOLUTION_TYPE);

      ptInfo->resultType = actResultType;
      ptInfo->resultName = actResultName;
      ptInfo->dofNames = StdVector<std::string>(dofNames);
      ptInfo->unit = unit;
      ptInfo->complexFormat = REAL_IMAG;
      ptInfo->entryType = H5IO::MapEntryType( entryType );
      ptInfo->definedOn = H5IO::MapUnknownType( definedOn );

      // perform consistency check
      if( ptInfo->entryType == ResultInfo::UNKNOWN ) {
        EXCEPTION( "Result '" << actResultName
                   << "' has no proper EntryType!" );
      }

      if( ptInfo->dofNames.GetSize() == 0 ) {
        EXCEPTION( "Result '" << actResultName
                   << "' has no degrees of freedoms!");
      }

      if( ptInfo->resultName == "" &&
          ptInfo->resultType == NO_SOLUTION_TYPE ) {
        EXCEPTION( "Result has neither a name nor a "
                   << "proper result type!" );
      }

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
      EXCEPTION( "Only results defined on nodes and elements "
                 << "can be read in from HDF5 file up to now" );
    }

    // iterate over all regions
    list.Clear();
    for( UInt i = 0; i < regions.GetSize(); i++ ) {
      list.Push_back( mi_->GetEntityList( listType, regions[i] ) );
    }
    resInfoGroup.close();
    actMsGroup.close();
  }

  void SimInputHDF5::GetResult( UInt sequenceStep,
                                UInt stepNum,
                                shared_ptr<BaseResult> result,
                                bool isHist) {

    if( !isHist ) {
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
      std::string pathsep = fs::path("/").string();
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
    //case ResultInfo::PFEM:
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

    StdVector<UInt> idx;
    UInt numDofs = result->GetResultInfo()->dofNames.GetSize();
    UInt numEntities = result->GetEntityList()->GetSize();
    UInt resVecSize =  numEntities * numDofs;

    //std::cout << "resVecSize " << resVecSize << std::endl;

    if(entString == "Nodes") {
      idx = entityNodeMap_[regionName];
      //std::cout << "entityNodeMap_[" << regionName << "] Size " << entityNodeMap_[regionName].GetSize() << std::endl;
    } else {
      idx.Resize( numEntities );

      for( UInt i = 0; i < numEntities ; i++ )
        idx[i] = i;
    }


    // copy data array to result object
    if( result->GetEntryType() == BaseMatrix::DOUBLE ) {
      Vector<Double> & resVec = dynamic_cast<Result<Double>& >(*result).GetVector();
      resVec.Resize( resVecSize );
      for( UInt i = 0; i < numEntities; i++ ) {
        for( UInt iDof = 0; iDof < numDofs; iDof++ ) {
          resVec[i*numDofs+iDof] = realVals[idx[i]*numDofs+iDof];
        }
      }
    } else {
      Vector<Complex> & resVec = dynamic_cast<Result<Complex>& >(*result).GetVector();
      StdVector<Double> imagVals;
      H5IO::ReadArray( resGroup, "Imag", imagVals );

      resVec.Resize( resVecSize );

      for( UInt i = 0; i < numEntities; i++ ) {
        for( UInt iDof = 0; iDof < numDofs; iDof++ ) {
          resVec[i*numDofs+iDof] = Complex( realVals[idx[i]*numDofs+iDof],
                                            imagVals[idx[i]*numDofs+iDof] );
        }
      }
    }
    resGroup.close();
    stepGroup.close();
    // close external file for current step
    if( hasExternalFiles_ )
      extFile.close();
  }

  void SimInputHDF5::GetHistResult( UInt sequenceStep, UInt stepNum,
                                         shared_ptr<BaseResult> result) {

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

    // get stepvalues map of result file to extract step indices
    std::map<UInt, Double> stepVals;
    std::map<shared_ptr<ResultInfo>, std::map<UInt, Double> > resultSteps;
    shared_ptr<ResultInfo> actRes = result->GetResultInfo();
    this->GetStepValues( sequenceStep, actRes, resultSteps[actRes], true);
    stepVals.insert( resultSteps[actRes].begin(), resultSteps[actRes].end() );
    // Get index of step (because in case of missing steps, stepIdx != stepNum-1)
    UInt stepIdx = std::distance(stepVals.begin(),stepVals.find(stepNum));

    // iterate over all entities in the the list
    shared_ptr<EntityList> list = result->GetEntityList();
    EntityIterator it = list->GetIterator();
    UInt numDofs = result->GetResultInfo()->dofNames.GetSize();

    // ----------------------------
    // Case 1: Real valued entries
    // ----------------------------
    if( result->GetEntryType() == BaseMatrix::DOUBLE ) {
      Vector<Double> & resVec =
        dynamic_cast<Result<Double>&>(*result).GetVector();
      resVec.Resize( list->GetSize() * numDofs );
      for( it.Begin(); !it.IsEnd(); it++ ) {

        H5::Group actEntGroup;

        if(entityTypeString == "Nodes") {
          // Find node number in HDF5 file which corresponds to mesh node number
          UInt meshNodeNum = it.GetNode();
          UInt fileNodeNum = meshNodeNum;
          if( !readAllEntities_) {
            fileNodeNum = g2FNodeNumMap_[meshNodeNum];
          }
          std::stringstream sstr;
          sstr << fileNodeNum;
          actEntGroup = entityGroup.openGroup(sstr.str());
        } else if(entityTypeString == "Elements") {
          // Find element number in HDF5 file which corresponds to mesh elem number
          UInt meshElemNum = it.GetElem()->elemNum;
          UInt fileElemNum = meshElemNum;
          if( !readAllEntities_) {
            fileElemNum = g2FElemNumMap_[meshElemNum];
          }
          std::stringstream sstr;
          sstr << fileElemNum;
          actEntGroup = entityGroup.openGroup(sstr.str());
        } else {
          // open for each entity the specific subgroup
          if(tempRegionName_ != ""){
            actEntGroup = entityGroup.openGroup( tempRegionName_ );
          }else{
            actEntGroup = entityGroup.openGroup(it.GetIdString() );
          }

        }

        // read single part of array and set it in the result vector
        StdVector<Double> vals;
        H5IO::ReadArray( actEntGroup, "Real", vals );

        for( UInt i = 0; i < numDofs; i++ ) {
          resVec[it.GetPos()*numDofs+i] = vals[stepIdx*numDofs+i];
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
            Complex( realVals[stepIdx*numDofs+i],
                     imagVals[stepIdx*numDofs+i] );
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

    // This sets needs to have defined order - not necessarily sorted. 
    // See comment for SimInputCDB::referencedElems_
    std::set<UInt> readNodeSet;
    std::set<UInt> readElemSet;

    // map for each element number the related region
    StdVector<RegionIdType> elemRegionMap(numElems);
    elemRegionMap.Init(NO_REGION_ID);

    // ensure, that region names are already read in
    if( !statsRead_ )
      ReadMeshStats( meshGroup );

    // Read all nodes from regions and initialize mapping from mesh node
    // numbers to grid node numbers accordingly.
    UInt baseNodeNum = mi_->GetNumNodes() + 1;
    UInt baseElemNum = mi_->GetNumElems() + 1;
    if (baseNodeNum > 1 || baseElemNum > 1) {
      readAllEntities_ = false; // We cannot take shortcuts if this is not the only grid
    }

    std::set<std::string>::iterator findIt;

    // ================================
    //  Add Elements Per Element Group
    // ================================
    for( UInt i = 0; i < elemNames_.GetSize(); i++ ) {
      H5::Group entityGroup, actEntityGroup;

      try{
        entityGroup = meshGroup.openGroup( "Groups" );
      } catch (H5::Exception&) {
        LOG_TRACE(simInputHdf5) << "Could not open group for entity groups";
        return;
      }

      for( UInt i = 0, n=elemNames_.GetSize(); i < n; i++ ) {
        findIt = readEntities_.find(elemNames_[i]);
        if(findIt == readEntities_.end())
          continue;

        // open entitygroup with given name
        try {
          actEntityGroup = entityGroup.openGroup( elemNames_[i] );
        } H5_CATCH( "Could not open definition of element group '"
                    << nodeNames_[i] << "'" );

        // Check if entity needs to be linearized
        bool linearizeEntity;
        linearizeEntity = linearizeEntities_.find(elemNames_[i]) != linearizeEntities_.end();

        // read elems from grid
        StdVector<UInt> readElems;
        H5IO::ReadArray( actEntityGroup, "Elements", readElems );
        for( UInt j = 0, n2=readElems.GetSize(); j < n2; j++ ) {
          const UInt elemNum = readElems[j];
          if(!readAllEntities_)
            readElemSet.insert( readElems[j] );
          elemRegionMap[elemNum-1] = NO_REGION_ID;
        }

        StdVector<UInt> readNodes;
        StdVector<UInt> actualNodes;
        StdVector<UInt> nodeIndices;
        // read nodes from grid
        H5IO::ReadArray( actEntityGroup, "Nodes", readNodes );
        actualNodes=readNodes;

        // Extract nodes needed for linear elements
        if(linearizeEntity) {
          LinearizeElems(readElems, elemTypes, globConnect, actualNodes);

          // To be able to read results for linearized grids we need an index
          // map for each entity (nodes of named elems and regions)
          nodeIndices.Resize(actualNodes.GetSize());
          for( UInt j = 0, n2=actualNodes.GetSize(); j < n2; j++ ) {
            nodeIndices[j] = readNodes.Find(actualNodes[j]);
          }
        }

        // For element lists which need not be linearized we simply use
        // ascending node indices.
        nodeIndices.Resize(actualNodes.GetSize());
        for( UInt j = 0, n2=actualNodes.GetSize(); j < n2; j++ ) {
          nodeIndices[j] = j;
        }

        entityNodeMap_[elemNames_[i]] = nodeIndices;
        if( !readAllEntities_)
          readNodeSet.insert( actualNodes.Begin(), actualNodes.End( ));

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

    for( UInt i = 0, n=regionNames_.GetSize(); i < n; i++ ) {
      std::string regionName = regionNames_[i];
      findIt = readEntities_.find(regionName);
      if(findIt == readEntities_.end())
        continue;

      try {
        actRegion = regionGroup.openGroup( regionName );
      } H5_CATCH( "Could not open group for region '" <<
                  regionNames_[i] << "'" );


      // pass region names to grid and obtain RegionIds
      RegionIdType actRegionId = mi_->AddRegion(regionName);

      // Check if entity needs to be linearized
      bool linearizeEntity;
      linearizeEntity = linearizeEntities_.find(regionName) != linearizeEntities_.end();

      // read elem numbers for this region
      StdVector<UInt> regionElems;
      H5IO::ReadArray( actRegion, "Elements", regionElems );
      for( UInt j = 0; j < regionElems.GetSize(); j++ ) {
        const UInt elemNum =regionElems[j];
        if( !readAllEntities_ ) 
          readElemSet.insert( elemNum );
        elemRegionMap[elemNum-1] = actRegionId;
      }

      StdVector<UInt> regionNodes;
      StdVector<UInt> actualNodes;
      StdVector<UInt> nodeIndices;

      H5IO::ReadArray( actRegion, "Nodes", regionNodes );
      actualNodes=regionNodes;

      if(linearizeEntity) {
        LinearizeElems(regionElems, elemTypes, globConnect, actualNodes);

        // To be able to read results for linearized grids we need an index
        // map for each entity (nodes of named elems and regions)
        nodeIndices.Resize(actualNodes.GetSize());
        for( UInt j = 0, n2=actualNodes.GetSize(); j < n2; j++ ) {
          nodeIndices[j] = regionNodes.Find(actualNodes[j]);
        }
      }

      // For element list which need not be linearized we simply use
      // ascending node indices.
      nodeIndices.Resize(actualNodes.GetSize());
      for( UInt j = 0, n2=actualNodes.GetSize(); j < n2; j++ ) {
        nodeIndices[j] = j;
      }

      entityNodeMap_[regionName] = nodeIndices;
      if( !readAllEntities_ )
        readNodeSet.insert( regionNodes.Begin(), regionNodes.End());
    }

    // ===================
    //  Add Nodes to Grid
    // ===================
    
    UInt idx;
    Vector<Double> p(3);
    std::set<UInt>::iterator it, end;
    
    if( !readAllEntities_ ) {
      numNodes_ = readNodeSet.size();
      mi_->AddNodes( readNodeSet.size() );
      for( it=readNodeSet.begin(), end=readNodeSet.end();
          it != end;
          it++ ) {

        idx = ((*it)-1)*3;
        p[0] = nodeCoords_[idx + 0];
        p[1] = nodeCoords_[idx + 1];
        p[2] = nodeCoords_[idx + 2];
        mi_->SetNodeCoordinate( baseNodeNum, p );

        // perform mapping, if not all entities are read in
        g2FNodeNumMap_[baseNodeNum] = *it;
        f2GNodeNumMap_[*it] = baseNodeNum++;
      } 
    } else {
      mi_->AddNodes( numNodes_);
      for( UInt iNode = 0; iNode < numNodes_; ++iNode ) {
        idx = iNode * 3;
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
      idx = baseElemNum;
      for( it=readElemSet.begin(), end=readElemSet.end();
          it != end;
          it++ ) {
        g2FElemNumMap_[idx] = *it;
        f2GElemNumMap_[*it] = idx++;
      }

      // Remap global connectivity from mesh nodes to grid node numbers
      for( UInt i = 0, n=globConnect.GetSize();
          i < n;
          i++ ) {
        globConnect[i] = f2GNodeNumMap_[globConnect[i]];
      }

      for( it = readElemSet.begin(), end = readElemSet.end();
          it != end;
          it++ ) {
        UInt elemNum=(*it);
        Elem::FEType type = (Elem::FEType) elemTypes[elemNum-1];
        UInt * connect = &globConnect[numNodesPerElem*(elemNum-1)];
        RegionIdType actRegionId = elemRegionMap[elemNum-1];
        mi_->SetElemData( f2GElemNumMap_[elemNum], type, actRegionId, connect );
      }

    } else {
      // === Add all elements to grid ===
      mi_->AddElems( numElems );
      for( UInt iElem = 0; iElem < numElems; ++iElem ) {
        UInt elemNum= iElem+1;
        Elem::FEType type = (Elem::FEType) elemTypes[elemNum-1];
        UInt * connect = &globConnect[numNodesPerElem*(elemNum-1)];
        RegionIdType actRegionId = elemRegionMap[elemNum-1];
        mi_->SetElemData( elemNum, type, actRegionId, connect );
      }
    }

    // ======================================
    //  Generate Node Groups For Each Region
    // ======================================

    // check, if nodes nodes of region should be additionally added
    // to list of named nodes
    if( genRegionNodes_) {
#if 0
      // iterate over all regions
      for( UInt i = 0; i < readEntities_.GetSize(); i++ ) {

        try {
          actRegion = regionGroup.openGroup( readEntities_[i] );
        } H5_CATCH( "Could not open group for region '" <<
                    regionNames_[i] << "'" );

        // read nodes of region
        StdVector<UInt> regionNodes;
        H5IO::ReadArray( actRegion, "Nodes", regionNodes );

        for( UInt j = 0, n=regionNodes.GetSize();
        j < n;
        j++ ) {
          regionNodes[j] = f2GNodeNumMap_[regionNodes[j]];
        }
        // add nodes as named nodes
        mi_->AddNamedNodes( readEntities_[i]+"_Nodes", regionNodes );

        regionGroup.close();
      }
#else
      EXCEPTION("Generation of nodes per region temporarily not supported!");
#endif
    }

  }

  void SimInputHDF5::ReadNodeGroups(const H5::Group& meshGroup)
  {

    H5::Group entityGroup, actEntityGroup;

    try{
      entityGroup = meshGroup.openGroup( "Groups" );
    } catch (H5::Exception&) {
      LOG_TRACE(simInputHdf5) << "Could not open group for entity groups";
      return;
    }

    std::set<std::string>::iterator findIt;

    for( UInt i = 0; i < nodeNames_.GetSize(); i++ ) {
      findIt = readEntities_.find(nodeNames_[i]);
      if(findIt == readEntities_.end())
        continue;

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
        if( !readAllEntities_ ) {
          if(f2GNodeNumMap_[readNodes[j]] != 0)
            nodes.Push_back(f2GNodeNumMap_[readNodes[j]]);
        } else {
          nodes.Push_back(readNodes[j]);
        }
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
     } catch (H5::Exception&) {
       LOG_TRACE(simInputHdf5) << "Could not open group for entity groups";
       return;
     }

     std::set<std::string>::iterator findIt;

     for( UInt i = 0; i < elemNames_.GetSize(); i++ ) {
       findIt = readEntities_.find(elemNames_[i]);
       if(findIt == readEntities_.end())
         continue;

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
         if( !readAllEntities_ ) {
           if(f2GElemNumMap_[readElems[j]] != 0)
             elems.Push_back(f2GElemNumMap_[readElems[j]]);
         } else {
           elems.Push_back(readElems[j]);
         }
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
    try {
      entityGroup = meshGroup.openGroup("Groups");
    } catch (H5::Exception&) {
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

  void SimInputHDF5::LinearizeElems(const StdVector<UInt>& readElems,
                                    StdVector<Integer>& elemTypes,
                                    StdVector<UInt>& globConnect,
                                    StdVector<UInt>& readNodes) {
    static std::map<Elem::FEType, Elem::FEType> elemTypeMap;
    UInt elemIncr = globConnect.GetSize() / elemTypes.GetSize();
    UInt numElems = readElems.GetSize();
    std::set<UInt> readNodeSet;

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

    for(UInt i=0; i<numElems; i++) {
      UInt elemIdx=readElems[i]-1;
      Elem::FEType newType = elemTypeMap[(Elem::FEType)elemTypes[elemIdx]];
      UInt numNodes = Elem::shapes[newType].numNodes;

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
    
    Vector<Double> p, globPoint;
    p.Resize(dim_);
    globPoint.Resize(dim_);

    for(UInt i=0; i<numNodes_; ++i) 
    {
      UInt idx = i*3;
      p[0] = nodeCoords_[idx + 0];
      p[1] = nodeCoords_[idx + 1];
      if (dim_ == 3) p[2] = nodeCoords_[idx + 2];
      coordSys.Global2LocalCoord(globPoint, p);
      
      nodeCoords_[idx + 0] = globPoint[0] * scaleFac;
      nodeCoords_[idx + 1] = globPoint[1] * scaleFac;
      if (dim_ == 3) nodeCoords_[idx + 2] = globPoint[2] * scaleFac;
    }
  }
  
  void SimInputHDF5::ReadStringFromUserData(const std::string& dSetName,
                                            std::string& str) {
    H5::Group userDataGroup;

    // If it does not exist, create Group for Data.
    try {
      userDataGroup = mainRoot_.openGroup("UserData");
    } H5_CATCH( "Can open 'UserData' group." );

    StdVector<std::string> userData;
    H5IO::ReadArray( userDataGroup, dSetName, userData );
    str = userData[0];

    userDataGroup.close();
  }


  // ------------------------------------------------------------------------
  //  DATABASE SECTION
  // ------------------------------------------------------------------------
  
  void SimInputHDF5::DB_Init() {
    try { 
      dbRoot_= mainRoot_.openGroup("DataBase");
    }  H5_CATCH( "HDF5 file has not database section" );

  }

  void SimInputHDF5::DB_GetParamFileContent( std::string& params ) {
    H5::Group extFiles;
    try {
      extFiles = dbRoot_.openGroup("InputFiles");
    } H5_CATCH( "Could not read group for external files");

    StdVector<std::string> vec;
    H5IO::ReadArray( extFiles, "ParameterFile", vec);
    params = vec[0];
  }

  void SimInputHDF5::DB_GetMatFileContent( std::string& params ) {
    H5::Group extFiles;
    try {
      extFiles = dbRoot_.openGroup("InputFiles");
    } H5_CATCH( "Could not read group for external files");

    StdVector<std::string> vec;
    H5IO::ReadArray( extFiles, "MaterialFile", vec);
    params = vec[0];
  }

  void SimInputHDF5::DB_GetFeFctCoefs( UInt sequenceStep, UInt stepNum,
                                       const std::string& pdeName,
                                       const std::string& fctName,
                                       SingleVector * coefs ) {
    
    
    H5::Group msGroup, mGroup, physGroup, fctGroup, stepGroup;
    // open multistep group
    try {
      msGroup = dbRoot_.openGroup("MultiSteps");
    } H5_CATCH( "Could not open multistep group");
    
    // open current ms step
    std::string msStepStr = lexical_cast<std::string>(sequenceStep);
    try{
      mGroup= msGroup.openGroup(msStepStr);
    } H5_CATCH( "Could not open database group for sequenceStep " 
        << sequenceStep );

    // open group for physics
    try{
      physGroup = mGroup.openGroup(pdeName);
    } H5_CATCH( "Could not open database group for physic " << pdeName );

    //  open group for FeFunction
    try{
      fctGroup = physGroup.openGroup(fctName);
    } H5_CATCH( "Could not open database group for result " 
                << fctName );

    // open group for current step
    std::string stepStr = lexical_cast<std::string>(stepNum);
    try {
        stepGroup = fctGroup.openGroup(stepStr);
    } H5_CATCH("Could not find data of result '" << fctName << "' in step " << stepStr << "\nPerhaps you need to activate <allowPostProc>");

    // Now aquire coefficients for fefunction
    if( coefs->GetEntryType() == BaseMatrix::DOUBLE ) {
      Vector<Double> & rVec = dynamic_cast<Vector<Double>& >(*coefs);
      
      // read data array
      StdVector<Double> realVals;
      H5IO::ReadArray( stepGroup, "Real", realVals );
      rVec.Resize( realVals.GetSize() );
      for( UInt i = 0; i < realVals.GetSize(); ++i ) {
        rVec[i] = realVals[i];
      }
      
    } else {
      Vector<Complex> & cVec = dynamic_cast<Vector<Complex>& > (*coefs);
      StdVector<Double> realVals, imagVals;
      H5IO::ReadArray( stepGroup, "Real", realVals );
      H5IO::ReadArray( stepGroup, "Imag", imagVals );
      cVec.Resize( realVals.GetSize() );
      for( UInt i = 0; i < realVals.GetSize(); ++i ) {
        cVec[i] = Complex(realVals[i], imagVals[i] );
      }
    }

    stepGroup.close();
    fctGroup.close();
    physGroup.close();
    mGroup.close();
    msGroup.close();
  }

  void SimInputHDF5::
  DB_GetNumMultiSequenceSteps( std::map<UInt, BasePDE::AnalysisType>& analysis,
                               std::map<UInt, Double>& accTime,
                               std::map<UInt, bool>& isFinished ){

    // Initialize DB
    DB_Init();
    
    // Open "MultiStep" section 
    H5::Group msGroup;
    // open multistep group
    try {
      msGroup = dbRoot_.openGroup("MultiSteps");
    } H5_CATCH( "Could not open multistep group")
    
    // Iterate over all children
    hsize_t numChildren = msGroup.getNumObjs();
    for( UInt i = 0; i < numChildren; ++i ) {
      UInt actMsStep = lexical_cast<UInt>(H5IO::GetObjNameByIdx( msGroup, i ));
      
      // open group and read analysyis type
      H5::Group actMsGroup;
      try {
        actMsGroup = msGroup.openGroup(lexical_cast<std::string>(actMsStep));
      } H5_CATCH( "Could not open multistep group with number " << actMsStep);

      // obtain type of analysis
      std::string analysisString;
      H5IO::ReadAttribute( actMsGroup, "AnalysisType", analysisString );
      analysis[actMsStep] = BasePDE::analysisType.Parse(analysisString); 
      H5IO::ReadAttribute( actMsGroup, "Completed", isFinished[actMsStep] );
      H5IO::ReadAttribute( actMsGroup, "AccTime", accTime[actMsStep] );
      
      actMsGroup.close();
    } // loop: multisequence steps
    
    msGroup.close();
  }
  
  void SimInputHDF5::
  DB_GetAvailPdeCoefFcts( UInt msStep, 
                          std::map<std::string, 
                          std::set<std::string> >& coefFcts ) {
    H5::Group msGroup, mGroup, physGroup, fctGroup;
    // open multistep group
    try {
      msGroup = dbRoot_.openGroup("MultiSteps");
    } H5_CATCH( "Could not open multistep group");

    // open current ms step
    std::string msStepStr = lexical_cast<std::string>(msStep);
    try{
      mGroup= msGroup.openGroup(msStepStr);
    } H5_CATCH( "Could not open database group for sequenceStep " 
        << msStep )

    // Loop over all physic types
    hsize_t numPhysic = mGroup.getNumObjs();
    for( UInt iPhysic = 0; iPhysic < numPhysic; ++iPhysic ) {
      H5::Group pGroup;
      std::string physic;
      physic = H5IO::GetObjNameByIdx( mGroup, iPhysic );
      try{
        pGroup= mGroup.openGroup(physic);
      } H5_CATCH( "Could not open phsysic database group " << physic);

      // Loop over all coefficient types
      hsize_t numCoefs = pGroup.getNumObjs();
      for( UInt iCoef = 0; iCoef < numCoefs; ++iCoef ) {
        std::string coef;
        coef= H5IO::GetObjNameByIdx( pGroup, iCoef );
        coefFcts[physic].insert(coef);
      } // loop: coefs
      
      pGroup.close();
    } // loop: physic types
    mGroup.close();
    msGroup.close();
  }

  void SimInputHDF5::
  DB_GetStepValues( UInt sequenceStep,
                    const std::string& pdeName,
                    const std::string& fctName,
                    std::map<UInt, Double>& stepValues ) {

    H5::Group msGroup, mGroup, physGroup, fctGroup, stepGroup;
    // open multistep group
    try {
      msGroup = dbRoot_.openGroup("MultiSteps");
    } H5_CATCH( "Could not open multistep group");

    // open current ms step
    std::string msStepStr = lexical_cast<std::string>(sequenceStep);
    try{
      mGroup= msGroup.openGroup(msStepStr);
    } H5_CATCH( "Could not open database group for sequenceStep " 
        << sequenceStep );

    // open group for physics
    try{
      physGroup = mGroup.openGroup(pdeName);
    } H5_CATCH( "Could not open database group for physic " << pdeName );

    //  open group for FeFunction
    try{
      fctGroup = physGroup.openGroup(fctName);
    } H5_CATCH( "Could not open database group for result " 
        << fctName );

    // Loop over all time steps
    hsize_t numSteps = fctGroup.getNumObjs();
    for( UInt iStep = 0; iStep < numSteps; ++iStep ) {
      H5::Group stepGroup;
      UInt stepNum = 
          lexical_cast<UInt>(H5IO::GetObjNameByIdx( fctGroup, iStep ));
      try{
        stepGroup = fctGroup.openGroup(lexical_cast<std::string>(stepNum));
      } H5_CATCH( "Could not read step " << stepNum << " for result " 
                  << fctName << " in database"); 

      // get time step value
      Double stepVal;
      H5IO::ReadAttribute( stepGroup, "StepValue", stepVal );
      stepValues[stepNum] = stepVal;
      stepGroup.close();
    } //loop: steps

    fctGroup.close();
    physGroup.close();
    mGroup.close();
    msGroup.close();

  }

  void SimInputHDF5::DB_Close() {
    try { 
      dbRoot_.close();
    }  H5_CATCH( "Could not close database section" );

  }

  void SimInputHDF5::GetNamedNodeResult(const std::string& nodeName,
                                        const std::string resultName,
                                        StdVector<Complex>& result)
  {    
    StdVector<UInt> nodeNumbers;
    StdVector<Double> real, imag;
    
    try {
      H5::Group groupsGroup = mainRoot_.openGroup("Mesh/Groups");

      H5::Group namedNodeGroup;
      try {      
        namedNodeGroup = groupsGroup.openGroup(nodeName);
      } catch (H5::Exception& h5Ex ) {
        EXCEPTION("Named node '" << nodeName << "' not found in HDF5 file!");
      }      

      H5IO::ReadArray( namedNodeGroup, "Nodes", nodeNumbers);
      
      namedNodeGroup.close();
      groupsGroup.close();
      
      for(UInt i=0, n=nodeNumbers.GetSize(); i<n; i++) 
      {
        std::cout << "Node numbers for " << nodeName << ": " << nodeNumbers[i] << std::endl;
      }
      
      std::stringstream sstr;
      
      try {
        // Try to read node results from history group.
        sstr << "Results/History/MultiStep_1/" << resultName << "/Nodes/" << nodeNumbers[0];

        H5::Group historyNodeGroup = mainRoot_.openGroup(sstr.str());
        H5IO::ReadArray( historyNodeGroup, "Real", real);
        H5IO::ReadArray( historyNodeGroup, "Imag", imag);
        historyNodeGroup.close();

        std::cout << "real.size " << real.GetSize() << std::endl;
        std::cout << "imag.size " << imag.GetSize() << std::endl;

        result.Resize(1);
        result[0] = Complex(real[0], imag[0]);

        std::cout << "Complex result " << result[0] << std::endl;

      }  catch (H5::Exception& h5Ex ) {
        // Try to read node results from history group.

        sstr.clear(); sstr.str("");
        sstr << "Results/Mesh/MultiStep_1/ResultDescription/" << resultName;
        H5::Group resDescGroup = mainRoot_.openGroup(sstr.str());

        StdVector< shared_ptr<ResultInfo> > infos;
        shared_ptr<ResultInfo> myInfo;
        GetResultTypes( 1, infos, false );
        
        for(UInt i=0, n= infos.GetSize(); i<n; i++)
        {
          if(infos[i]->resultName == resultName) 
          {
            myInfo = infos[i];
            break;
          }          
        }
        

        StdVector<shared_ptr<EntityList> > list;
        shared_ptr<EntityList> myList;

        GetResultEntities( 1, myInfo, list, false );

        UInt dist = 0;
        for(UInt i=0, n= list.GetSize(); i<n; i++)
        {
          const StdVector<UInt>& nodes = dynamic_cast<NodeList*>(list[i].get())->GetNodes();
          std::string name = dynamic_cast<NodeList*>(list[i].get())->GetName();

          UInt numNodes = nodes.GetSize();
          const UInt* pt = std::find(&nodes[0], &nodes[nodes.GetSize()-1]+1, nodeNumbers[0]);
          dist = std::distance(&nodes[0], pt);

          std::cout << "Searching for node " << nodeNumbers[0] << " in list " << name << std::endl;
          
          if(dist < numNodes) 
          {
            myList = list[i];
            std::cout << "Found node!" << std::endl;
            break;
          }
        }
        

        Result<Complex>* resPt = new Result<Complex>;
        shared_ptr<BaseResult> res( resPt );
        res->SetResultInfo(myInfo);
        res->SetEntityList(myList);
        
        GetMeshResult( 1, 1, res );

        Vector<Complex> resVec =  resPt->GetVector();
        result.Resize(1);
        result[0] = resVec[dist * myInfo->dofNames.GetSize()];

        std::cout << "Complex result " << result[0] << std::endl;
        // EXCEPTION( STR << ":\n" << h5Ex.getCDetailMsg() );
      }

    }  H5_CATCH( "Could not close database section" );
  }

}
