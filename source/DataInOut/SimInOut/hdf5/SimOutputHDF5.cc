// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <algorithm>

#include <boost/filesystem/fstream.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#include <def_cfs_stats.hh>

#include "DataInOut/ProgramOptions.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "General/Exception.hh"
#include "SimOutputHDF5.hh"
#include "Domain/ElemMapping/Elem.hh"
#include "FeBasis/BaseFE.hh"
#include "hdf5io.hh"

namespace CoupledField {

#define H5_EXCEPTION(STR, EX)                                           \
  EXCEPTION( STR, EX.getCDetailMsg() );

#define H5_CATCH( STR )                                                 \
  catch (H5::Exception& h5Ex ) {                                        \
    EXCEPTION( STR << ":\n" << h5Ex.getCDetailMsg() );                  \
  }

  DEFINE_LOG(h5Out, "hdf5Out")



  SimOutputHDF5::SimOutputHDF5(std::string fileName, PtrParamNode inputNode,
                               PtrParamNode infoNode, bool isRestart ) :
    SimOutput(fileName, inputNode, infoNode, isRestart) {

    LOG_DBG(h5Out) << "SO fn=" << fileName << " rs=" << isRestart;

    fileName_ = fileName;
    formatName_ = "hdf5";
    isInitialized_ = false;
    
    std::string dirString = "results_" + formatName_; 
    inputNode->GetValue("directory", dirString, ParamNode::PASS );
    dirName_ = dirString; 
    
    useDataBase_ = false;
    
    capabilities_.insert( MESH );
    capabilities_.insert( MESH_RESULTS );
    capabilities_.insert( HISTORY );
    capabilities_.insert( DATABASE );

    gridWritten_ = false;
    externalFiles_ = false;
    printGridOnly_ = false;

    currMS_ = 0;
    currStep_ = 0;

    // Initialize data layout (compression, maxChunkSize)
    // using values specified by user
    UInt compressionLevel = 6;
    UInt maxChunkSize = 100;
    myParam_->GetValue("compressionLevel", compressionLevel, ParamNode::PASS );
    if( compressionLevel > 9) {
      EXCEPTION( "Value for compressionLevel must be between 0 and 9" );
    }
    myParam_->GetValue("maxChunkSize", maxChunkSize, ParamNode::PASS );
    dPropList_ = H5::DSetCreatPropList::DEFAULT;
    if (maxChunkSize > 0 || compressionLevel > 0) {
      dPropList_.setLayout( H5D_CHUNKED );
    } else {
      dPropList_.setLayout( H5D_CONTIGUOUS );
    }
    if (compressionLevel > 0) {
      dPropList_.setDeflate( compressionLevel );
      if (maxChunkSize == 0) {
        EXCEPTION("HDF5 compression level > 0 requires a maxChunckSize > 0");
      }
    } else if (maxChunkSize > 0) {
      H5IO::SetMaxChunkSize( maxChunkSize );
    }

    // Change defaults according to XML file
    myParam_->GetValue("externalFiles", externalFiles_, ParamNode::PASS);

    H5::Exception::dontPrint();
    
    std::string extString = "cfs";
    inputNode->GetValue("extension", extString, ParamNode::PASS );

    std::string fName = fileName_ + "." + extString;
    currFileName_ = fs::path(dirName_ / fName).string();
  }


  SimOutputHDF5::~SimOutputHDF5() {
    CloseFile();
  }

  void SimOutputHDF5::Init( Grid* ptGrid, bool printGridOnly ) {
    LOG_DBG(h5Out) << "Init";
    ptGrid_ = ptGrid;
    printGridOnly_ = printGridOnly;
    WriteGrid();
  }

  void SimOutputHDF5::WriteFileInfoHeader() {
    LOG_DBG(h5Out) << "WFIH";
    H5::Group infoGroup;
    try {
      infoGroup = mainGroup_.openGroup( "FileInfo" );
    } H5_CATCH( "Could not open group for FileInfo" );

    // write file version
    std::stringstream version;
    version << CFS_HDF5_FORMAT_MAJOR << "." << CFS_HDF5_FORMAT_MINOR;
    std::string versionString = version.str();
    H5IO::Write1DArray( infoGroup, "Version", 1, &versionString, dPropList_ );

    // write date / time information
    using namespace boost::posix_time;
    using namespace boost::gregorian;
    std::string now = to_simple_string( second_clock::local_time() );
    H5IO::Write1DArray( infoGroup, "Date", 1, &now, dPropList_ );

    // write creator
    std::stringstream creator;
    // we omit CFS_GIT_COMMIT to prevent recompiles, the info is also in the info.xml and the branch tells probably more
    creator << "openCFS " << CFS_VERSION << ", " << CFS_NAME << " (" << CFS_GIT_BRANCH << ")";
    std::string creatorString = creator.str();
    H5IO::Write1DArray( infoGroup, "Creator", 1, &creatorString, dPropList_ );

    // create group for content
    StdVector<Integer> content;
    std::set<Capability>::iterator it;
    for( it = usedCapabilities_.begin();
         it != usedCapabilities_.end(); it++ ) {
      content.Push_back( H5IO::MapCapabilityType( *it ) );
    }
    H5IO::Write1DArray( infoGroup, "Content", content.GetSize(),
                        &content[0], dPropList_ );

  }


  void SimOutputHDF5::BeginMultiSequenceStep( UInt step,
                                              BasePDE::AnalysisType type,
                                              UInt numSteps  ) {
    LOG_DBG(h5Out) << "BMSS step=" << step << " at=" << BasePDE::analysisType.ToString(type) << " ns=" << numSteps;
    std::stringstream msName;
    H5::Group resultDescGroup;
    
    // acquire lock
    LockFile();

    currMSNumSteps_ = numSteps;

    // If it does not exist, create Group for Grid / Volume data
    try {
      resultsGroup_ = mainGroup_.openGroup("Results");
    } H5_CATCH( "Could not open group for results" );

    // Map analysistype
    std::string analysisType = BasePDE::analysisType.ToString(type);
    currAnalysisType_ = type;

    // Assemble name of multistep
    msName << "MultiStep_" << step;

    // 1) Write group for mesh results
    if( registeredMeshResults_.size() != 0 ) {
      meshResultsGroup_ = H5IO::OpenCreateGroup(resultsGroup_, "Mesh");
      // write attribute indicating use of external files for simlation
      // steps
      H5IO::WriteAttribute( meshResultsGroup_, "ExternalFiles", externalFiles_ );

      try {
        // create new multistep group.
        currMSMeshGroup_ = H5IO::OpenCreateGroup(meshResultsGroup_, msName.str());

        // add analysistype and number of steps to group
        H5IO::WriteAttribute( currMSMeshGroup_, "AnalysisType", analysisType );
        H5IO::WriteAttribute( currMSMeshGroup_, "LastStepNum", (UInt) 0 );
        H5IO::WriteAttribute( currMSMeshGroup_, "LastStepValue", (Double) 0.0);

        // add a group for the result description datasets.
        resultDescGroup = H5IO::OpenCreateGroup(currMSMeshGroup_, "ResultDescription");

        // write result meta information to file
        WriteResultDescriptions( resultDescGroup, numSteps, false);
        resultDescGroup.close();
      } H5_CATCH( "Could not create group for multi sequence step " << step );
    }

    // 2) Write group for history results
    if( registeredHistResults_.size() != 0 ) {
      histResultsGroup_ = H5IO::OpenCreateGroup(resultsGroup_, "History" );

      try {
        currMSHistGroup_ = H5IO::OpenCreateGroup(histResultsGroup_, msName.str() );

        // add analysistype and number of steps to group
        H5IO::WriteAttribute( currMSHistGroup_, "AnalysisType", analysisType );
        H5IO::WriteAttribute( currMSHistGroup_, "LastStepNum", (UInt) 0  );
        H5IO::WriteAttribute( currMSHistGroup_, "LastStepValue", (Double) 0.0 );
        
        // add a group for the result description datasets.
        resultDescGroup = H5IO::OpenCreateGroup(currMSHistGroup_,"ResultDescription");

        // write result meta information to file
        WriteResultDescriptions( resultDescGroup, numSteps, true );
        resultDescGroup.close();

        // iterate over all results
        ResDescType::iterator it;
        for( it = registeredHistResults_.begin();
            it != registeredHistResults_.end();
            it++ ) {

          // create for each result a group within the ms group
          H5::Group resultGroup = H5IO::OpenCreateGroup(currMSHistGroup_, it->first);

          // create subgroup for entitytype
          ResultInfo::EntityUnknownType definedOn
          = it->second[0]->GetResultInfo()->definedOn;
          std::string entityString = H5IO::MapUnknownTypeAsString(definedOn );
          H5::Group entityTypeGroup = H5IO::OpenCreateGroup( resultGroup,  entityString );

          // iterate over all entitylists of result and create sub-subgroup
          std::vector<shared_ptr<BaseResult> > const & lists = it->second;
          for( UInt iList = 0; iList < lists.size(); iList++ )  {

            // iterate over all entities in this list
            EntityIterator entIt = lists[iList]->GetEntityList()->GetIterator();
            for( entIt.Begin(); !entIt.IsEnd(); entIt++ ) {
              H5::Group entityGroup; 
              try {
                entityGroup = entityTypeGroup.openGroup( entIt.GetIdString() );
                // In the restart case it is okay to have the group already present
                if( !isRestart_)
                WARN("You are trying to add history entity '" << entIt.GetIdString()
                     << "' under group '"
                     << "History/" << msName.str() << "/" << it->first << "/" << entityString 
                     << "'\nwhich already exists under a different name! Please check your mesh and XML files.");
                continue;
              } catch( H5::Exception& h5Ex ) {
                entityGroup = entityTypeGroup.createGroup( entIt.GetIdString() );
              }
            }
          }
          entityTypeGroup.close();
        }
      } H5_CATCH( "Could not create group for multi sequence step " << step );
    }
    currMS_ = step;

    // release lock
    UnlockFile();
  }

  void SimOutputHDF5::RegisterResult( shared_ptr<BaseResult> sol,
                                      UInt saveBegin, UInt saveInc,
                                      UInt saveEnd,
                                      bool isHistory ) {
    std::string resultName = sol->GetResultInfo()->resultName;

    LOG_DBG(h5Out) << "RS sol=" << resultName << " sb=" << saveBegin << " se=" << saveEnd << " inc=" << saveInc << " hist=" << isHistory;

    if( !isHistory ) {
      registeredMeshResults_[resultName].push_back(sol);
      meshResultSaveBegin_[resultName] = saveBegin;
      meshResultSaveEnd_[resultName] =  saveEnd;
      meshResultSaveInc_[resultName] = saveInc;
      usedCapabilities_.insert(MESH_RESULTS);
    } else {
      registeredHistResults_[resultName].push_back(sol);
      histResultSaveBegin_[resultName] = saveBegin;
      histResultSaveEnd_[resultName] = saveEnd;
      histResultSaveInc_[resultName] = saveInc;
      usedCapabilities_.insert(HISTORY);
    }
  }

  void SimOutputHDF5::BeginStep( UInt stepNum, Double stepVal ) {
    LOG_DBG(h5Out) << "BS num=" << stepNum << " v=" << stepVal;
    currStep_ = stepNum;
    currStepValue_ = stepVal;
  }


  void SimOutputHDF5::AddResult( shared_ptr<BaseResult> sol )
  {
    // acquire lock
    LockFile();
    
    std::string resultName = sol->GetResultInfo()->resultName;

    LOG_DBG(h5Out) << "AR sol=" << resultName;

    // try to determine, if current result is a history or mesh result
    bool isHistory = false;
    if( registeredHistResults_.find(resultName) !=
      registeredHistResults_.end() ) {
      if( std::find( registeredHistResults_[resultName].begin(),
                     registeredHistResults_[resultName].end(),
                     sol) !=
                       registeredHistResults_[resultName].end() ) {
        isHistory = true;
      }
    }

    if( !isHistory) {
      AddMeshResult( sol );
    } else {
      AddHistResult( sol );
    }

    // continuously update the "lastStep" value and number
    // to ensure a consistent hdf5 file throughout the 
    // simulation, even in case of a Ctrl-C action
    if( registeredMeshResults_.size() > 0 ) {
      H5IO::WriteAttribute( currMSMeshGroup_, "LastStepNum", currStep_ );
      H5IO::WriteAttribute( currMSMeshGroup_, "LastStepValue", currStepValue_ );
    }

    if( registeredHistResults_.size() > 0 ) {
      H5IO::WriteAttribute( currMSHistGroup_, "LastStepNum", currStep_ );
      H5IO::WriteAttribute( currMSHistGroup_, "LastStepValue", currStepValue_ );

    }

    // release lock
    UnlockFile();
  }

  void SimOutputHDF5::AddMeshResult( shared_ptr<BaseResult> sol) {

    LOG_DBG(h5Out) << "AMR sol=" << sol->GetResultInfo()->resultName;
    // No need to aquire lock, as this method just gets called from
    // SimOutputHDF5::AddResult()
    
    H5::Group resultGroup, subGroup, regionGroup;
    UInt numDOFs;
    Vector<Double> realVec, imagVec;
    std::vector<std::string> resultNames;

    std::string regionName = sol->GetEntityList()->GetName();
    shared_ptr<ResultInfo> resInfo = sol->GetResultInfo();
    std::string resultName = resInfo->resultName;
    numDOFs = resInfo->dofNames.GetSize();

    // check, if step group is already open
    if( currMeshStepGroup_.getId() <= 0 ) {
      std::stringstream stepName;
      stepName << "Step_" << currStep_;

      // Create new step group.
      try {
        currMeshStepGroup_= currMSMeshGroup_.openGroup(stepName.str());
      } catch (H5::Exception& ) {
        try {
          // TODO: strieben - Creation of step groups fails when converting from h5 to h5 (see cube2d). Do something about it!
          currMeshStepGroup_= currMSMeshGroup_.createGroup(stepName.str());  
          H5IO::WriteAttribute( currMeshStepGroup_, "StepValue", currStepValue_ );

          if(externalFiles_ )
            CreateExternalFile();
        } H5_CATCH( "Can not create dataset for step " << currStep_ );
      }
    }

    // Add current stepvalue to related group in result description,
    // if not yet present
    bool writeStep = false;
    StdVector<UInt> & myStepNums = meshResultStepNums_[resultName];
    StdVector<Double> & myStepVals = meshResultStepVal_[resultName];
    if( myStepNums.GetSize() == 0 ) {
      writeStep = true;
    } else {
      if( myStepNums.Last() != currStep_ )
        writeStep = true;
    }

    if( writeStep ) {
      myStepNums.Push_back( currStep_ );
      myStepVals.Push_back( currStepValue_ );
      try {
        H5::Group resDescGroup =
          currMSMeshGroup_.openGroup( "ResultDescription").openGroup(resultName);

        // "extend" the StepValues and StepNumbers array
        StdVector<Double> tmp_double(1);
        tmp_double[0] =  currStepValue_;
        H5IO::Extend1DArray( resDescGroup, "StepValues", myStepVals.GetSize(), 
                             tmp_double.GetPointer(), dPropList_ );

        StdVector<Double> tmp_uint(1);
        tmp_uint[0]  = currStep_;
        H5IO::Extend1DArray( resDescGroup, "StepNumbers", myStepNums.GetSize(), 
                                     tmp_uint.GetPointer(), dPropList_ );
        resDescGroup.close();
      } H5_CATCH( "Could not write current step value for result '"
                  << resultName << "'" );
    }

    // check, if result was already written
    try {
      resultGroup = currMeshStepGroup_.openGroup( resultName );
    } catch( H5::Exception& ) {
      resultGroup = currMeshStepGroup_.createGroup( resultName );
    }

    // determine, on which type of entity the result is defined
    std::string entityString;
    switch( resInfo->definedOn ) {
    case ResultInfo::NODE:
      entityString = "Nodes";
      break;
    case ResultInfo::ELEMENT:
    case ResultInfo::SURF_ELEM:
      entityString = "Elements";
      break;
    default:
      EXCEPTION( "Result of type '" << resInfo->resultName
                  << "'can not be written as mesh result" );
    }

    // try to create regionGroup
    try {
      regionGroup = resultGroup.openGroup( regionName );
    } catch( H5::Exception& ) {
      regionGroup = resultGroup.createGroup( regionName );
    }

    // try to create subgroup for entity
    try {
      LOG_DBG2(h5Out) << "Create subgroup " << entityString
          << " for result " << resultName << " on region " << regionName
          << " in step " << currStep_;
      subGroup = regionGroup.createGroup( entityString );
    } H5_CATCH( "Could not create subgroup " << entityString
                << " for result " << resultName << " on region "
                << regionName << ". Maybe the group already exists.");

    if( sol->GetEntryType() == BaseMatrix::DOUBLE ) {

      Vector<Double> & resultVec = dynamic_cast<Result<Double>&>
      (*sol).GetVector();

      WriteResults(subGroup, resultVec, numDOFs, false);
    } else {
      Vector<Complex> & resultVec = dynamic_cast<Result<Complex>&>
      (*sol).GetVector();

      realVec.Resize( resultVec.GetSize() );
      imagVec.Resize( resultVec.GetSize() );

      for(UInt i = 0; i < resultVec.GetSize(); i++) {
        realVec[i] = resultVec[i].real();
        imagVec[i] = resultVec[i].imag();
      }
      WriteResults(subGroup, realVec, numDOFs, false);
      WriteResults(subGroup, imagVec, numDOFs, true);

    }

    // close groups
    subGroup.close();
    regionGroup.close();
    resultGroup.close();
  }

  void SimOutputHDF5::AddHistResult( shared_ptr<BaseResult> sol) {

    // No need to aquire lock, as this method just gets called from
    // SimOutputHDF5::AddResult()
    
    shared_ptr<ResultInfo> resInfo = sol->GetResultInfo();
    std::string resultName = resInfo->resultName;
    UInt numDofs = resInfo->dofNames.GetSize();
    std::string entityString = H5IO::MapUnknownTypeAsString(resInfo->definedOn );

    LOG_DBG(h5Out) << "AHR sol=" << resultName;


    // Add current stepvalue to related group in result description,
    // if not yet present
    bool writeStep = false;
    StdVector<UInt> & myStepNums = histResultStepNums_[resultName];
    StdVector<Double> & myStepVals = histResultStepVal_[resultName];
    if( myStepNums.GetSize() == 0 ) {
      writeStep = true;
    } else {
      if( myStepNums.Last() != currStep_ )
        writeStep = true;
    }
    if( writeStep ) {
      myStepNums.Push_back( currStep_ );
      myStepVals.Push_back( currStepValue_ );

      try {
        H5::Group resDescGroup =
            currMSHistGroup_.openGroup( "ResultDescription").openGroup(resultName);
        StdVector<Double> tmp_double(1);
        tmp_double[0] =  currStepValue_;
        H5IO::Extend1DArray( resDescGroup, "StepValues", myStepVals.GetSize(), 
                             tmp_double.GetPointer(), dPropList_ );

        StdVector<Double> tmp_uint(1);
        tmp_uint[0]  = currStep_;
        H5IO::Extend1DArray( resDescGroup, "StepNumbers", myStepNums.GetSize(), 
                             tmp_uint.GetPointer(), dPropList_ );
        
        resDescGroup.close();
      } H5_CATCH( "Could not write current step value for result '"
                  << resultName << "'" );
    }
    // ---------------------
    //  Write result itself
    // ---------------------
    try {
      H5::Group resultGroup = currMSHistGroup_.openGroup( resultName);
      H5::Group entityTypeGroup = resultGroup.openGroup( entityString );


      // iterate over all entities in this list
      EntityIterator entIt = sol->GetEntityList()->GetIterator();
      UInt pos = 0;
      for( entIt.Begin(); !entIt.IsEnd(); entIt++ ) {
        H5::Group entityGroup =
          entityTypeGroup.openGroup( entIt.GetIdString() );

        if( sol->GetEntryType() == BaseMatrix::DOUBLE ){
          Vector<Double> & resultVec = dynamic_cast<Result<Double>&>
          (*sol).GetVector();
          H5IO::Extend2DArray( entityGroup, "Real", myStepNums.GetSize(),
                               numDofs,  &resultVec[pos] );
          pos += numDofs;
        } else {
          Vector<Complex> & resultVec = dynamic_cast<Result<Complex>&>
          (*sol).GetVector();
          Vector<Double> realVec(numDofs), imagVec(numDofs);

          for( UInt i = 0; i < numDofs; i++ ) {
            realVec[i] = resultVec[pos+i].real();
            imagVec[i] = resultVec[pos+i].imag();
          }

          H5IO::Extend2DArray( entityGroup, "Real", myStepNums.GetSize(),
                               numDofs, &realVec[0] );
          H5IO::Extend2DArray( entityGroup, "Imag", myStepNums.GetSize(),
                               numDofs, &imagVec[0] );
          pos += numDofs;
        }
        entityGroup.close();
      }
      entityTypeGroup.close();
      resultGroup.close();
    } H5_CATCH( "Could not open history result group for result '"
                << resultName << "'" );
  }


  void SimOutputHDF5::FinishStep( )
  {
    LOG_DBG(h5Out) << "FS";

    LockFile();
    
    if(externalFiles_ && myInfo_)
    {
      PtrParamNode in = myInfo_->Get("analysis/output/externalFile");
      try {
        in->Get("name")->SetValue(currStepFile_.getFileName());
        in->Get("size")->SetValue((int) currStepFile_.getFileSize());
      } catch (H5::FileIException &h5ex) {}
    }

   // we close everything here, so that the file is in a consistent state
    currStepFile_.close();

    if( currMeshStepGroup_.getId() > 0 )
      currMeshStepGroup_.close();

    if( currHistStepGroup_.getId() > 0 )
      currHistStepGroup_.close();

    // Release lock
    UnlockFile();
  }

  void SimOutputHDF5::FinishMultiSequenceStep( ) {
    LOG_DBG(h5Out) << "FMSS";
    // Acquire lock
    LockFile();

    // close groups, which were opened in BeginMultiSequenceStep()
    if( registeredMeshResults_.size() > 0 ) {
      currMSMeshGroup_.close();
      meshResultsGroup_.close();
      registeredMeshResults_.clear();
    }

    if( registeredHistResults_.size() > 0 ) {
      currMSHistGroup_.close();
      histResultsGroup_.close();
      registeredHistResults_.clear();
    }
    resultsGroup_.close();

    // reset all data per sequence step
    meshResultSaveBegin_.clear();
    meshResultSaveEnd_.clear();
    meshResultSaveInc_.clear();
    histResultSaveBegin_.clear();
    histResultSaveEnd_.clear();
    histResultSaveInc_.clear();
    meshResultStepNums_.clear();
    histResultStepNums_.clear();
    meshResultStepVal_.clear();
    histResultStepVal_.clear();

    // Release lock
    UnlockFile();
  }

  void SimOutputHDF5::Finalize() {
    LOG_DBG(h5Out) << "Finalize";
    // return in case of an already restarted file
    if (isRestart_)
      return;
    
    // Acquire lock
    LockFile();
        
    // Write file header
    WriteFileInfoHeader();

    // return, if only the grid is to be printed
    if( printGridOnly_ )
      return;

    // return, if no commandLine handler or
    // global root ParaemNode are present
    if( !progOpts || !myParam_ )
      return;

    std::vector<std::string> fileNames;
    std::vector<std::string> dataSetNames;
    std::ifstream fin;
    std::ostringstream dumpStr;

    fileNames.push_back(progOpts->GetParamFileStr());
    fileNames.push_back(myParam_->GetRoot()->Get("fileFormats")
                        ->Get("materialData")
                        ->Get("file")->As<std::string>());

    dataSetNames.push_back("ParameterFile");
    dataSetNames.push_back("MaterialFile");

    for(UInt i=0; i<fileNames.size(); i++)
    {
      fin.open( fileNames[i].c_str(), std::ios::binary );

      if(fin.fail())
        EXCEPTION("Cannot open file '" << fileNames[i]
                  <<"' to dump into HDF5!");

      // seek to the end of the file
      fin.seekg (0, std::ios::end);
      UInt numBytes = fin.tellg();
      fin.seekg (0, std::ios::beg);

      std::string str;
      str.resize(numBytes);
      fin.read(&str[0], numBytes);
      WriteStringToUserData(dataSetNames[i], str);
      fin.close();
    }

    progOpts->PrintVersion( dumpStr, false );
    WriteStringToUserData( "ProgramStats", dumpStr.str() );

    // Release lock
    UnlockFile();
  }

  void SimOutputHDF5::InitModule() {
    LOG_DBG(h5Out) << "IM";
    if( isInitialized_)
      return;

    // concatenate output file name
    try {
      fs::create_directory( dirName_ );
    } catch (std::exception &ex) {
      EXCEPTION(ex.what());
    }
    
    // In case of re-start, we simply append information
    bool truncate = !isRestart_;
    OpenFile(truncate);
    isInitialized_ = true;
  }
  
  void SimOutputHDF5::OpenFile(bool truncate){
    LOG_DBG(h5Out) << "OF truncate=" << truncate;
    // create main file and obtain main group
    try {
      mainFile_ = H5::H5File (currFileName_, truncate ? H5F_ACC_TRUNC : H5F_ACC_RDWR );
    } H5_CATCH( "Could not create hdf5 file '" << currFileName_ << "' : " );

    mainGroup_ = mainFile_.openGroup( "/" );
    if(truncate){
      meshGroup_ = mainGroup_.createGroup( "Mesh" );
      mainGroup_.createGroup( "FileInfo" ).close();

      mainGroup_.createGroup( "UserData" ).close();
      mainGroup_.createGroup( "Results" ).close();
    }else{
      meshGroup_ = mainGroup_.openGroup("Mesh");
    }    
  }
  
  void SimOutputHDF5::CloseFile(){
    LOG_DBG(h5Out) << "CF";
    if(currStepFile_.getLocId() > 0) {
      if (currStepFile_.getObjCount( H5F_OBJ_DATASET |
                                     H5F_OBJ_GROUP |
                                     H5F_OBJ_DATATYPE | H5F_OBJ_ATTR) > 0 ) {
        std::cerr << "There are still objects open in the hdf5 file "
                  << currStepFile_.getFileName() << "\n\n";
        H5IO::CheckOpenObjects(currStepFile_, true);
      }
      
      currStepFile_.close();
    }

//todo    if(resultsGroup_.getLocId() > 0){
//      resultsGroup_.close();
//    }
    if(meshGroup_.getLocId() > 0){
      meshGroup_.close();
    }

    // check, if any group is open at all
    if( dbGroup_.getLocId() > 0 )
      dbGroup_.close();
    
    // check, if any group is open at all
    if( currMsDbGroup_.getLocId() > 0 )
      currMsDbGroup_.close();
    
    // check, if any group is open at all
    if( mainGroup_.getLocId() > 0 )
      mainGroup_.close();

    // check for open groups, datasets etc.
    if (mainFile_.getLocId() > 0 )
    {
      if (mainFile_.getObjCount( H5F_OBJ_DATASET |
                                 H5F_OBJ_GROUP |
                                 H5F_OBJ_DATATYPE | H5F_OBJ_ATTR) > 0 ) {
        std::cerr << "There are still objects open in the hdf5 file "
                  << mainFile_.getFileName() << "\n\n";
        H5IO::CheckOpenObjects(mainFile_, true);
      }

      mainFile_.close();
    }
  }

  void SimOutputHDF5::WriteGrid() {
    LOG_DBG(h5Out) << "WG";
    // ensure that grid gets only written once
    if(!gridWritten_)
      InitModule();
    else
      return;

    // in case of restart, we do not re-write the grid!
    if( isRestart_) {
      gridWritten_ = true;
      return;
    }
    
    // Lock the file
    LockFile();
    
    // write the dimension of the grid.
    H5IO::WriteAttribute( meshGroup_, "Dimension", ptGrid_->GetDim() );

    // ================
    //  Node Locations
    // ================
    UInt nNodes = ptGrid_->GetNumNodes();
    H5::Group nodeGroup;
    try {
      nodeGroup = meshGroup_.createGroup( "Nodes" );
    } H5_CATCH( "Could not create node group" );

    H5IO::WriteAttribute( nodeGroup, "NumNodes", nNodes );

    // collect all nodal coordinates
    std::vector<Double> locs( nNodes * 3 );
    for (UInt i = 0; i < nNodes; i++) {
      Vector<Double> p;
      ptGrid_->GetNodeCoordinate3D(p, i+1);
      locs[i*3+0] = p[0];
      locs[i*3+1] = p[1];
      locs[i*3+2] = p[2];
    }

    H5IO::Write2DArray( nodeGroup, "Coordinates", nNodes,
                        3, &locs[0], dPropList_ );

    nodeGroup.close();


    // =====================
    //  Element definitions
    // =====================
    UInt nElems = ptGrid_->GetNumElems();
    H5::Group elemGroup;
    try{
      elemGroup = meshGroup_.createGroup("Elements");
    } H5_CATCH( "Could not create element group" );

    UInt maxNumNodes = ptGrid_->GetMaxNumNodesPerElem();
    std::vector<UInt> connect (nElems * maxNumNodes);
    std::vector<UInt> elConnect;
    std::vector<Integer> feTypes (nElems);
    std::vector< UInt > numElemsOfDim ( 4 );

    // Fill connectivity array
    std::fill( connect.begin(), connect.end(), 0 );
    UInt offset;
    Elem::FEType eType;
    RegionIdType region;

    // iterate over all elements
    for( UInt i = 0; i < nElems; i++ ) {
      elConnect.resize( maxNumNodes );
      std::fill(elConnect.begin(), elConnect.end(),
                0 );
      ptGrid_->GetElemData( i+1, eType, region, &elConnect[0] );
      numElemsOfDim[(Elem::shapes[eType].dim)]++;
      feTypes[i] = eType;

      // insert connectivity into global array
      offset = i * maxNumNodes;
      for( UInt j = 0; j < elConnect.size(); j++ ) {
        connect[offset + j] = elConnect[j];
      }
    }

    // write connectivity
    H5IO::Write2DArray( elemGroup, "Connectivity", nElems,
                        maxNumNodes, &connect[0], dPropList_ );

    // write element types
    H5IO::Write1DArray( elemGroup, "Types", nElems,
                        &feTypes[0], dPropList_ );


    // ==========================
    //  Grid Meta Information
    // ==========================

    H5IO::WriteAttribute( elemGroup, "NumElems", nElems );

    H5IO::WriteAttribute( elemGroup, "QuadraticElems",
                          ptGrid_->IsQuadratic() );

    // number of elements per dimension
    for(UInt i=0; i<4; i++) {
      std::stringstream attrName;
      attrName << "Num" << (i) << "DElems";
      H5IO::WriteAttribute( elemGroup, attrName.str(), numElemsOfDim[i] );
    }

    // number of elements per type
    UInt numElemTypes = Elem::feType.map.size();
    for(UInt i=0; i<numElemTypes; i++) {
      std::stringstream attrName;
      attrName << "Num_" << Elem::feType.ToString((Elem::FEType)i);
      H5IO::WriteAttribute( elemGroup, attrName.str(),
                            ptGrid_->GetNumElemOfType((Elem::FEType)i) );
    }

    // close element group
    elemGroup.close();

    // ============================================
    //  Write Regions, Node Groups, Element Groups
    // ============================================

    WriteRegions( meshGroup_ );

    // create new group for entity groups
    H5::Group groupsGroup;
    try{
      groupsGroup = meshGroup_.createGroup("Groups");
    } H5_CATCH( "Could not create mesh regiongroup" );

    WriteNodeGroups( groupsGroup );
    WriteElemGroups( groupsGroup );
    groupsGroup.close();

    gridWritten_ = true;
    usedCapabilities_.insert(MESH);

    // close meshGroup
    meshGroup_.close();
    
    // release lock
    UnlockFile();
  }


  void SimOutputHDF5::WriteRegions(const H5::Group& meshGroup) {
    LOG_DBG(h5Out) << "WR";
    H5::Group regionListGroup;
    StdVector< std::string > regionNames;
    StdVector< UInt > regionDims;
    StdVector< std::vector<UInt> > regionElems;
    StdVector< StdVector<UInt> > regionNodes;
    StdVector<Elem*> elems;
    StdVector<RegionIdType> surfRegionIds, volRegionIds;
    UInt dim, numRegions;
    Integer idx;

    // create region group
    try{
      regionListGroup = meshGroup.createGroup("Regions");
    } H5_CATCH( "Could not create region group" );

    numRegions = ptGrid_->GetNumRegions();
    if(!numRegions)
      return;

    dim = ptGrid_->GetDim();
    ptGrid_->GetVolRegionIds(volRegionIds);
    ptGrid_->GetSurfRegionIds(surfRegionIds);
    ptGrid_->GetRegionNames(regionNames);
    regionDims.Resize(numRegions);
    regionElems.Resize(numRegions);
    regionNodes.Resize(numRegions);
    regionDims.Resize(numRegions);

    // obtain nodes and elements per surface region
    for(UInt i=0, n=surfRegionIds.GetSize(); i<n; i++) {
      idx = surfRegionIds[i];
      regionDims[idx] = dim-1;

      ptGrid_->GetElems(elems, idx);
      UInt nElems = elems.GetSize();
      regionElems[idx].resize(nElems);
      for(UInt j=0; j<nElems; j++) {
        regionElems[idx][j] = elems[j]->elemNum;
      }

      ptGrid_->GetNodesByRegion(regionNodes[idx], idx);

      // determine dimensionality
      if( ptGrid_->GetDim() == 3 ) {
        regionDims[idx] = 2;
      } else {
        regionDims[idx] = 1;
      }

    }

    // obtain nodes and elements per volume region
    for(UInt i=0, n=volRegionIds.GetSize(); i<n; i++) {
      idx = volRegionIds[i];
      regionDims[idx] = dim;

      ptGrid_->GetElems(elems, idx);
      UInt nElems = elems.GetSize();
      regionElems[idx].resize(nElems);
      for(UInt j=0; j<nElems; j++) {
        regionElems[idx][j] = elems[j]->elemNum;
      }

      ptGrid_->GetNodesByRegion(regionNodes[idx], idx);

      // determine dimensionality
      regionDims[idx] = ptGrid_->GetDim();
    }


    // loop over regions and write out nodes and elements
    for(UInt i = 0; i < numRegions; i++)
    {
      // create new region group
      H5::Group actRegionGroup;
      try {
        actRegionGroup = regionListGroup.createGroup(regionNames[i] );
      } H5_CATCH( "Could not create region group for region '"
                  << regionNames[i] << "'" );
      LOG_DBG2(h5Out) << "WR: write region " << regionNames[i];
      H5IO::WriteAttribute( actRegionGroup, "Dimension",
                            regionDims[i] );

      // create new node group
      H5IO::Write1DArray<UInt>( actRegionGroup, "Nodes",
                          regionNodes[i].GetSize(),
                          (const UInt*)&regionNodes[i][0], dPropList_ );

      // create new element group
      H5IO::Write1DArray( actRegionGroup,
                          "Elements",
                          regionElems[i].size(),
                          (const Integer*)&regionElems[i][0],
                          dPropList_);

      // create new face group
      // .. to be implemented ..

      // create new edge group
      // .. to be implemented ..

      // close current region group
      actRegionGroup.close();
    }

    // close regionlist group
    regionListGroup.close();


  }

  void SimOutputHDF5::WriteNodeGroups(const H5::Group& meshGroup) {
    LOG_DBG(h5Out) << "WNG";
    H5::Group myGroup;
    StdVector< UInt > nodes;
    StdVector<std::string> nodeNames;
    UInt numNodeGroups = 0;

    // obtain list with names of nodes
    ptGrid_->GetListNodeNames(nodeNames);
    numNodeGroups = nodeNames.GetSize();

    for(UInt i = 0; i < numNodeGroups; i++ ) {
      ptGrid_->GetNodesByName(nodes, nodeNames[i]);

      // try to open group with given name
      try {
        myGroup = meshGroup.openGroup( nodeNames[i] );
      } catch (H5::Exception& ) {
        myGroup = meshGroup.createGroup( nodeNames[i] );
      }
      H5IO::WriteAttribute( myGroup, "Dimension", (Integer) 0 );
      H5IO::Write1DArray( myGroup, "Nodes",
                          nodes.GetSize(), &nodes[0], dPropList_ );

      // close nodes array of current group
      myGroup.close();
    }
  }

  void SimOutputHDF5::WriteElemGroups(const H5::Group& meshGroup) {
    LOG_DBG(h5Out) << "WEG";
    H5::Group myGroup;
    StdVector< UInt > elemNums, elemNodes;
    StdVector<Elem*> elems;
    StdVector<std::string> elemNames;
    std::set<UInt> nodeSet;

    // obtain list with names of elements
    ptGrid_->GetListElemNames(elemNames);
    UInt numElemGroups = elemNames.GetSize();

    for(UInt i = 0; i < numElemGroups; i++ ) {
      ptGrid_->GetElemsByName(elems, elemNames[i]);
      elemNums.Resize( elems.GetSize() );
      nodeSet.clear();
      
      std::set<UInt> dims;

      for( UInt j = 0; j < elems.GetSize(); j++ ) {
        elemNums[j] = elems[j]->elemNum;
        dims.insert( Elem::shapes[elems[j]->type].dim);
        nodeSet.insert( elems[j]->connect.Begin(),
                        elems[j]->connect.End() );
      }
      if( dims.size() > 1 ) {
        EXCEPTION( "Element group '" << elemNames[i]
                    << "' contains elements of different dimensions" );
      }
      try {
        myGroup = meshGroup.openGroup( elemNames[i] );
      } catch (H5::Exception& ) {
        myGroup = meshGroup.createGroup( elemNames[i] );
      }
      H5IO::WriteAttribute( myGroup, "Dimension", *dims.begin() );
      H5IO::Write1DArray( myGroup, "Elements",
                          elemNums.GetSize(), &elemNums[0],
                          dPropList_);

      // Write nodes of element group
      elemNodes.Resize( nodeSet.size() );
      std::copy( nodeSet.begin(), nodeSet.end(), elemNodes.Begin() );
      H5IO::Write1DArray( myGroup, "Nodes",
                          elemNodes.GetSize(), &elemNodes[0],
                          dPropList_);

      // close nodes array of current group
      myGroup.close();
    }
  }

  void SimOutputHDF5::WriteResultDescriptions( const H5::Group& resGroup,
                                               UInt numSteps,
                                               bool isHistory ) {
    LOG_DBG(h5Out) << "WRD num=" << numSteps << " h=" << isHistory;
    std::string resultName, unit;
    UInt definedOn, numDofs, entryType; 
    //UInt saveBegin, saveEnd, saveInc;
    std::vector<std::string> entityNames;
    ResDescType::const_iterator it, end;
    std::vector< boost::shared_ptr<BaseResult> >::const_iterator solIt, solEnd;
    boost::shared_ptr<ResultInfo> resInfo, actResInfo;

    if( !isHistory ) {
      it = registeredMeshResults_.begin();
      end = registeredMeshResults_.end();
    } else {
      it = registeredHistResults_.begin();
      end = registeredHistResults_.end();
    }

    for( ; it != end; it++ ) {
      resultName = it->first;

      solIt = it->second.begin();
      solEnd = it->second.end();

      resInfo = (*solIt)->GetResultInfo();
      numDofs = resInfo->dofNames.GetSize();
      unit = resInfo->unit;
      definedOn = H5IO::MapUnknownType( resInfo->definedOn );
      entryType = H5IO::MapEntryType( resInfo->entryType );

      // Generate list of entityNames for the current result.
      entityNames.clear();
      for( ; solIt != solEnd; solIt++ ) {
        actResInfo = (*solIt)->GetResultInfo();

        entityNames.push_back((*solIt)->GetEntityList()->GetName());
      }

      // Reset solIt to beginning of result vector
      solIt = it->second.begin();

      try {

        // Generate compound datatype

        /* // First version: Compound data type
           H5IO::CompoundType resInfo;
           typedef std::pair<std::string, boost::any> CEntryType;
           resInfo.push_back( CEntryType( "DefinedOn", definedOn ) );
           resInfo.push_back( CEntryType( "Regions", regions ) );
           resInfo.push_back( CEntryType( "NumDOFs", numDOFs ) );
           resInfo.push_back( CEntryType( "DOFNames", dofNames ) );
           resInfo.push_back( CEntryType( "EntryType", entryType ) );
           resInfo.push_back( CEntryType( "Unit", unit ) );

           H5IO::WriteCompound( currAttrDescGroup_, resNames[0], resInfo );
        */

        // Check, if the group for the group for the result has to be created.
        // which is either the case if the simulation is not restarted or
        // if the simulation is restarted, but the restarted sequenceStep
        // is not the current one.
        if( !isRestart_ || (isRestart_ && !H5IO::GroupExists(resGroup, resultName) ) )
        {
          H5::Group actGroup = resGroup.createGroup(resultName );

          H5IO::Write1DArray( actGroup, "DefinedOn", 1, &definedOn, dPropList_ );
          H5IO::Write1DArray( actGroup, "EntityNames", entityNames.size(), &entityNames[0], dPropList_ );
          H5IO::Write1DArray( actGroup, "NumDOFs", 1, &numDofs, dPropList_ );
          H5IO::Write1DArray( actGroup, "DOFNames", resInfo->dofNames.GetSize(), &(resInfo->dofNames[0]), dPropList_ );
          H5IO::Write1DArray( actGroup, "EntryType", 1, &entryType, dPropList_ );
          H5IO::Write1DArray( actGroup, "Unit", 1, &unit, dPropList_ );
          // In order to write a valid entry, we also set the initial stepNumber/
          // stepValues array
          UInt dummyStepNum = 0;
          H5IO::Extend1DArray(actGroup, "StepNumbers", 0, &dummyStepNum, dPropList_ );
          Double dummyStepVal = 0.0;
          H5IO::Extend1DArray(actGroup, "StepValues", 0, &dummyStepVal, dPropList_ );
          
          actGroup.close();
        }
        else
        { // general restart case
          H5::Group actGroup = resGroup.openGroup(resultName );
          
          // Obtain already written stepValue and stepNumbers and
          // store them in the meshResultStepNums_, meshResultStepVal_.
          // These array may not be created yet.
          if( H5IO::DatasetExists(actGroup, "StepNumbers") && H5IO::DatasetExists(actGroup, "StepValues") )
          {
            if( isHistory ) {

              StdVector<UInt> & oldStepNums = histResultStepNums_[resultName];
              H5IO::ReadArray( actGroup, "StepNumbers", oldStepNums);

              StdVector<Double> & oldStepVals= histResultStepVal_[resultName];
              H5IO::ReadArray( actGroup, "StepValues", oldStepVals);
            } else {
              // Obtain already written stepValue and stepNumbers and
              // store them in the meshResultStepNums_, meshResultStepVal_
              StdVector<UInt> & oldStepNums = meshResultStepNums_[resultName];
              H5IO::ReadArray( actGroup, "StepNumbers", oldStepNums);


              StdVector<Double> & oldStepVals= meshResultStepVal_[resultName];
              H5IO::ReadArray( actGroup, "StepValues", oldStepVals);

            }
          }
          actGroup.close();
        }
      } H5_CATCH( "Could not write result description for result '" << resultName << "'" );
    } //loop: registered mesh / history results
  }

  void SimOutputHDF5::WriteResults( H5::Group& resultGroup,
                                    Vector<Double>& resultVals,
                                    const UInt numDOFs,
                                    const bool isImag ) {
    LOG_DBG(h5Out) << "WR n=" << resultVals.GetSize() << " nd=" << numDOFs << " imag=" << isImag;
    // create dataset with related name
    std::string name;
    if( !isImag )
      name = "Real";
    else
      name = "Imag";

    UInt numEntities = (UInt) resultVals.GetSize() / numDOFs;
    
    H5IO::Write2DArray( resultGroup, name,
                        numEntities, numDOFs, &resultVals[0],
                        dPropList_ );
  }

  void SimOutputHDF5::CreateExternalFile() {
    LOG_DBG(h5Out) << "CEF";

    std::stringstream fName, masterGroup;
    std::string pathsep, fn;

    try {

      // open external file
      pathsep = fs::path("/").string();

      fName << fileName_ << "_ms" << currMS_ << "_step"
            << currStep_ << ".h5";
      fn = fName.str();
      fName.str("");
      fName << dirName_.string() << pathsep << fn;
         
      currStepFile_ = H5::H5File(fName.str(), H5F_ACC_TRUNC);

      // Write reference to external file to main file
      H5IO::WriteAttribute( currMeshStepGroup_, "ExtHDF5FileName", fn );

      // set current step group to external file
      currMeshStepGroup_.close();
      currMeshStepGroup_ = currStepFile_.openGroup("/");

      // Store reference to master file in external file
      fName.str("");
      fName << fileName_ << ".h5";
      fn = fName.str();
      H5IO::WriteAttribute( currMeshStepGroup_, "MasterHDF5FileName", fn );

      // Store reference to main group in external file
      masterGroup << "/Results/Mesh/MultiStep_" << currMS_
                  << "/Step_" << currStep_;
      H5IO::WriteAttribute( currMeshStepGroup_, "MasterGroup", masterGroup.str() );
    } H5_CATCH( "Could not create external file" );
  }

  void SimOutputHDF5::WriteStringToUserData(const std::string& dSetName,
                                            const std::string& str) {
    H5::Group userDataGroup;

    // If it does not exist, create Group for Data.
    userDataGroup = H5IO::OpenCreateGroup(mainGroup_, "UserData" );

    H5IO::Write1DArray( userDataGroup, dSetName, 1, &str, dPropList_ );
    userDataGroup.close();
  }
  
  
  void SimOutputHDF5::LockFile() {

    // disable temporarily the Ctr+c key
    //signal(SIGINT, SIG_IGN);
  }

  void SimOutputHDF5::UnlockFile() {
    // flush file to get a consistent state
    mainFile_.flush(H5F_SCOPE_GLOBAL);

    // re-enable default behavior for Ctr+C
    //signal(SIGINT, SIG_DFL); 
  }

  
  // ------------------------------------------------------------------------
  //  DATABASE SECTION
  // ------------------------------------------------------------------------
  void SimOutputHDF5::DB_Init() {
    InitModule();
    useDataBase_ = true;

    // insert capability
    usedCapabilities_.insert(DATABASE);

    if( isRestart_) {
      dbGroup_ = mainGroup_.openGroup("DataBase");
    } else {
      // add "DataBase" group
      try {
        dbGroup_ = mainGroup_.createGroup("DataBase");
        dbGroup_.createGroup("MultiSteps").close();
      } H5_CATCH( "Could not create section for internal database")
    }

  }
  
  void SimOutputHDF5::DB_WriteXmlFiles( fs::path simFile, fs::path matFile ) {
    if( isRestart_)
      return;
    
    fs::ifstream fin;
    std::ostringstream dumpStr;
    H5::Group extFiles;
    try
    {
      extFiles = dbGroup_.createGroup("InputFiles");
    } H5_CATCH( "Could not create group for external files");
    
    // open external Files
    StdVector<fs::path> filePaths(2);
    StdVector<std::string> setNames(2);
    filePaths[0] = simFile;
    setNames[0] = "ParameterFile";
    filePaths[1] = matFile;
    setNames[1] = "MaterialFile";

    for(UInt i=0; i<filePaths.GetSize(); i++)
    {
      fin.open( filePaths[i], std::ios::binary );

      if(fin.fail())
        EXCEPTION("Cannot open file '" << filePaths[i]
                  <<"' to dump into HDF5!");

      // seek to the end of the file
      fin.seekg (0, std::ios::end);
      UInt numBytes = fin.tellg();
      fin.seekg (0, std::ios::beg);

      std::string str;
      str.resize(numBytes);
      fin.read(&str[0], numBytes);
      fin.close();
      
      // now write out string
      H5IO::Write1DArray( extFiles, setNames[i], 1, &str, dPropList_ );
    }
    extFiles.close();
  }

  void SimOutputHDF5::DB_WritePythonFile( fs::path pythonFile) {
    if(isRestart_)
      return;

    fs::ifstream fin;
    std::ostringstream dumpStr;
    H5::Group extFiles;
    meshResultsGroup_ = H5IO::OpenCreateGroup(dbGroup_, "InputFiles");

    // open external file
    fs::path filePath;
    std::string setName;
    filePath = pythonFile;
    setName = "PythonFile";

    fin.open( filePath, std::ios::binary );

    if(fin.fail())
      EXCEPTION("Cannot open file '" << filePath << "' to dump into HDF5!");

    // seek to the end of the file
    fin.seekg (0, std::ios::end);
    UInt numBytes = fin.tellg();
    fin.seekg (0, std::ios::beg);

    std::string str;
    str.resize(numBytes);
    fin.read(&str[0], numBytes);
    fin.close();

    // now write out string
    H5IO::Write1DArray( extFiles, setName, 1, &str, dPropList_ );
    extFiles.close();
  }

  void SimOutputHDF5::DB_BeginMultiSequenceStep( UInt step,
                                                 BasePDE::AnalysisType type ) {
    currMS_ = step;
    currAnalysisType_ = type;
    std::string stepStr = lexical_cast<std::string>(step);
    H5::Group msGroup = dbGroup_.openGroup("MultiSteps");
    try {  
      currMsDbGroup_ = msGroup.openGroup( stepStr );
    } catch (H5::Exception&) {
      try {
        currMsDbGroup_ = msGroup.createGroup( stepStr );
      } H5_CATCH( "Could not open database group for step" << step );
    } 
    
    // set attributes for number of steps and analysis types
    std::string analysisType = BasePDE::analysisType.ToString(type);
    H5IO::WriteAttribute( currMsDbGroup_, "AnalysisType", analysisType );
    H5IO::WriteAttribute( currMsDbGroup_, "AccTime", 0.0 );
    //! Only write attribute, if it not exists yet
    if( !H5IO::AttrExists( currMsDbGroup_, "Completed" ) ) {
      H5IO::WriteAttribute( currMsDbGroup_, "Completed", false );
    }
  }

  void SimOutputHDF5::DB_BeginStep( UInt stepNum, Double stepVal ) {
    currStepDb_ = stepNum;
    currStepValueDb_ = stepVal;
    DB_BeginMultiSequenceStep(currMS_, currAnalysisType_ );
  }

  void SimOutputHDF5::DB_WriteFeFunction( const std::string& pdeCplName,
                                          const std::string& fctName,
                                          SingleVector* coefs ) {


    H5::Group physGroup, fctGroup, stepGroup;
    // Create / open group for physics
    try{
      physGroup = currMsDbGroup_.openGroup(pdeCplName);
    } catch (H5::Exception&) {
      try {
        physGroup = currMsDbGroup_.createGroup(pdeCplName);
      } H5_CATCH( "Could not create database group for physic " << pdeCplName );
    }

    // Create / open group for FeFunction
    try{
      fctGroup = physGroup.openGroup(fctName);
    } catch (H5::Exception&) {
      try {
        fctGroup = physGroup.createGroup(fctName);
      } H5_CATCH( "Could not create database group for physic " << pdeCplName );
    }

    // Create group for current step
    std::string stepStr = lexical_cast<std::string>(currStepDb_);
    stepGroup = fctGroup.createGroup(stepStr);
    H5IO::WriteAttribute( stepGroup, "StepValue", currStepValueDb_ );
    
    // Now write coefficients for fefunction
    if( coefs->GetEntryType() == BaseMatrix::DOUBLE ) {
      Vector<Double> & rVec = dynamic_cast<Vector<Double>& >(*coefs);
      WriteResults( stepGroup, rVec, 1, false );
      
    } else {
      Vector<Complex> & cVec = dynamic_cast<Vector<Complex>& > (*coefs);
      Vector<Double> realVec, imagVec;
      realVec.Resize( cVec.GetSize() );
      imagVec.Resize( cVec.GetSize() );

      for(UInt i = 0; i < cVec.GetSize(); i++) {
        realVec[i] = cVec[i].real();
        imagVec[i] = cVec[i].imag();
      }
      WriteResults(stepGroup, realVec, 1, false);
      WriteResults(stepGroup, imagVec, 1, true);

    }
    
    // now close groups
    stepGroup.close();
    fctGroup.close();
    physGroup.close();
  }
  
  void SimOutputHDF5::DB_FinishMultiSequenceStep( bool completed, 
                                                  Double accTime ) {
    H5IO::WriteAttribute( currMsDbGroup_, "Completed", completed );
    H5IO::WriteAttribute( currMsDbGroup_, "AccTime", accTime );
    
  }
   
  

} // end of namespace
