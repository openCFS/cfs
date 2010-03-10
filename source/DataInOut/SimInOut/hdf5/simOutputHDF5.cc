// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <algorithm>

#include "boost/date_time/posix_time/posix_time.hpp"
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/exception.hpp>

#include <def_cfs_stats.hh>

#include "DataInOut/programOptions.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/exception.hh"
#include "simOutputHDF5.hh"
#include "Domain/elem.hh"
#include "Elements/basefe.hh"
#include "hdf5io.hh"

namespace fs = boost::filesystem;

namespace CoupledField {

#define H5_EXCEPTION(STR, EX)                                           \
  EXCEPTION( STR, EX.getCDetailMsg() );

#define H5_CATCH( STR )                                                 \
  catch (H5::Exception& h5Ex ) {                                        \
    EXCEPTION( STR << ":\n" << h5Ex.getCDetailMsg() );                  \
  }


  SimOutputHDF5::SimOutputHDF5(std::string fileName, ParamNode * inputNode) :
    SimOutput(fileName, inputNode) {

    fileName_ = fileName;
    formatName_ = "hdf5";
    dirName_ = "results_" + formatName_;
    inputNode->Get("directory", dirName_, false );
    
    capabilities_.insert( MESH );
    capabilities_.insert( MESH_RESULTS );
    capabilities_.insert( HISTORY );

    gridWritten_ = false;
    externalFiles_ = false;
    printGridOnly_ = false;

    currMS_ = 0;
    currStep_ = 0;

    // Initialize data layout (compression, maxChunkSize)
    // using values specified by user
    UInt compressionLevel = 6;
    UInt maxChunkSize = 100;
    myParam_->Get("compressionLevel", compressionLevel, false );
    if( compressionLevel > 9) {
      EXCEPTION( "Value for compressionLevel must be betwen 1 and 9" );
    }
    myParam_->Get("maxChunkSize", maxChunkSize, false );
    dPropList_ = H5::DSetCreatPropList::DEFAULT;
    dPropList_.setLayout( H5D_CHUNKED );
    dPropList_.setDeflate( compressionLevel );
    H5IO::SetMaxChunkSize( maxChunkSize );

    // Change defaults according to XML file
    myParam_->Get("externalFiles", externalFiles_, false);

    // Do not print HDF5 exceptions by default
    H5::Exception::dontPrint();
  }


  SimOutputHDF5::~SimOutputHDF5() {
    // check for open groups, datasets etc. in current step file
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

    // check, if any group is open at all
    if( mainGroup_.getLocId() > 0 )
      mainGroup_.close();

    // check for open groups, datasets etc.
    if (mainFile_.getObjCount( H5F_OBJ_DATASET |
                               H5F_OBJ_GROUP |
                               H5F_OBJ_DATATYPE | H5F_OBJ_ATTR) > 0 ) {
      std::cerr << "There are still objects open in the hdf5 file "
                << mainFile_.getFileName() << "\n\n";
      H5IO::CheckOpenObjects(mainFile_, true);
    }

    mainFile_.close();
  }

  void SimOutputHDF5::Init( Grid* ptGrid, bool printGridOnly ) {
    ptGrid_ = ptGrid;
    printGridOnly_ = printGridOnly;
    WriteGrid();
  }

  void SimOutputHDF5::WriteFileInfoHeader() {
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
    creator << "CFS++ " << CFS_VERSION << ", " << CFS_NAME << " ( r" << CFS_SUBVERSION_REV << " )";
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
    std::stringstream msName;
    H5::Group resultDescGroup;


    currMSNumSteps_ = numSteps;

    // If it does not exist, create Group for Grid / Volume data
    try {
      resultsGroup_ = mainGroup_.openGroup("Results");
    } H5_CATCH( "Could open group for results" );

    // Map analysistype
    std::string analysisType = BasePDE::analysisType.ToString(type);
    currAnalysisType_ = type;

    // Assemble name of multistep
    msName << "MultiStep_" << step;

    // 1) Write group for mesh results
    if( registeredMeshResults_.size() != 0 ) {
      try{
        meshResultsGroup_ = resultsGroup_.openGroup("Mesh");
      } catch (H5::Exception&) {
        try {
          meshResultsGroup_ = resultsGroup_.createGroup("Mesh");
          // write attribute indicating use of external files for simlation
          // steps
          H5IO::WriteAttribute( meshResultsGroup_, "ExternalFiles", externalFiles_ );
        } H5_CATCH( "Could not create group for mesh results" );
      }

      try {
        // create new multistep group.
        currMSMeshGroup_ = meshResultsGroup_.createGroup(msName.str());

        // add analysistype and number of steps to group
        H5IO::WriteAttribute( currMSMeshGroup_, "AnalysisType", analysisType );

        // add a group for the result description datasets.
        resultDescGroup = currMSMeshGroup_.createGroup("ResultDescription");

        // write result meta information to file
        WriteResultDescriptions( resultDescGroup, numSteps, false);
        resultDescGroup.close();
      } H5_CATCH( "Could not create group for multi sequence step " << step );
    }

    // 2) Write group for history results
    if( registeredHistResults_.size() != 0 ) {
      try{
        histResultsGroup_ = resultsGroup_.openGroup("History");
      } catch (H5::Exception&) {
        try {
          histResultsGroup_ = resultsGroup_.createGroup("History");
        } H5_CATCH( "Could not create group for history results" );
      }

      try {
        // create new multistep group.
        currMSHistGroup_ = histResultsGroup_.createGroup(msName.str());

        // add analysistype and number of steps to group
        H5IO::WriteAttribute( currMSHistGroup_, "AnalysisType", analysisType );

        // add a group for the result description datasets.
        resultDescGroup = currMSHistGroup_.createGroup("ResultDescription");

        // write result meta information to file
        WriteResultDescriptions( resultDescGroup, numSteps, true );
        resultDescGroup.close();

        // iterate over all results
        ResDescType::iterator it;
        for( it = registeredHistResults_.begin();
             it != registeredHistResults_.end();
             it++ ) {

          // create for each result a group within the ms group
          H5::Group resultGroup;
          try {
            resultGroup = currMSHistGroup_.openGroup(it->first);
          } catch( H5::Exception& ) {
            resultGroup = currMSHistGroup_.createGroup(it->first);
          }

          // calculate number of "real steps"
          UInt saveBegin, saveEnd, saveInc, numRealSteps;
          saveBegin = std::max( (UInt) 1, histResultSaveBegin_[it->first] );
          saveEnd = std::min( currMSNumSteps_, histResultSaveEnd_[it->first] );
          saveInc = histResultSaveInc_[it->first];
          numRealSteps = (UInt) (saveEnd-saveBegin) / saveInc + 1;

          // create subgroup for entitytype
          ResultInfo::EntityUnknownType definedOn
              = it->second[0]->GetResultInfo()->definedOn;
          std::string entityString = H5IO::MapUnknownTypeAsString(definedOn );
          H5::Group entityTypeGroup = resultGroup.createGroup( entityString );

          // iterate over all entitylists of result and create sub-subgroup
          std::vector<shared_ptr<BaseResult> > const & lists = it->second;
          for( UInt iList = 0; iList < lists.size(); iList++ )  {

            // iterate over all entities in this list
            EntityIterator entIt = lists[iList]->GetEntityList()->GetIterator();
            UInt numDofs = lists[iList]->GetResultInfo()->dofNames.GetSize();
            for( entIt.Begin(); !entIt.IsEnd(); entIt++ ) {
              H5::Group entityGroup; 
              try {
                entityGroup = entityTypeGroup.openGroup( entIt.GetIdString() );
                WARN("You are trying to add history entity '" << entIt.GetIdString()
                     << "' under group '"
                     << "History/" << msName.str() << "/" << it->first << "/" << entityString 
                     << "'\nwhich already exists under a different name! Please check your mesh and XML files.");
                entityGroup.close();

                continue;
              } catch( H5::Exception& h5Ex ) {
                entityGroup = entityTypeGroup.createGroup( entIt.GetIdString() );
              }

              H5IO::Reserve2DArray<Double>(entityGroup, "Real", numRealSteps,
                                          numDofs, dPropList_ );
              if( lists[iList]->GetEntryType() == BaseMatrix::COMPLEX){
                H5IO::Reserve2DArray<Double>(entityGroup, "Imag", numRealSteps,
                                             numDofs, dPropList_ );
              }
              entityGroup.close();
            }
          }
          entityTypeGroup.close();
        }


      } H5_CATCH( "Could not create group for multi sequence step " << step );
    }
    currMS_ = step;
  }

  void SimOutputHDF5::RegisterResult( shared_ptr<BaseResult> sol,
                                      UInt saveBegin, UInt saveInc,
                                      UInt saveEnd,
                                      bool isHistory ) {

    std::string resultName = sol->GetResultInfo()->resultName;

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

    currStep_ = stepNum;
    currStepValue_ = stepVal;

  }


  void SimOutputHDF5::AddResult( shared_ptr<BaseResult> sol )
  {
    std::string resultName = sol->GetResultInfo()->resultName;

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
  }

  void SimOutputHDF5::AddMeshResult( shared_ptr<BaseResult> sol) {

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
        // TODO: strieben - Creation of step groups fails when converting from h5 to h5 (see cube2d). Do something about it!
        currMeshStepGroup_= currMSMeshGroup_.createGroup(stepName.str());
        H5IO::WriteAttribute( currMeshStepGroup_, "StepValue", currStepValue_ );

        if(externalFiles_ )
          CreateExternalFile();
      } H5_CATCH( "Can not create dataset for step " << currStep_ );
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
        H5IO::SetEntries1DArray(resDescGroup, "StepValues", myStepNums.GetSize()-1,
                                myStepNums.GetSize()-1,
                                &myStepVals[myStepNums.GetSize()-1] );
        H5IO::SetEntries1DArray(resDescGroup, "StepNumbers", myStepNums.GetSize()-1,
                                        myStepNums.GetSize()-1,
                                        &myStepNums[myStepNums.GetSize()-1] );
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

    // determine, on whicht type of entity the result is defined
    std::string entityString;
    switch( resInfo->definedOn ) {
    case ResultInfo::NODE:
    case ResultInfo::PFEM:
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
      subGroup = regionGroup.createGroup( entityString );
    } H5_CATCH( "Could not create subgroup " << entityString
                << " for result " << resultName << " on region "
                << regionName );

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

    shared_ptr<ResultInfo> resInfo = sol->GetResultInfo();
    std::string resultName = resInfo->resultName;
    UInt numDofs = resInfo->dofNames.GetSize();
    std::string entityString = H5IO::MapUnknownTypeAsString(resInfo->definedOn );


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
        H5IO::SetEntries1DArray(resDescGroup, "StepValues", myStepNums.GetSize()-1,
                                myStepNums.GetSize()-1,
                                &myStepVals[myStepNums.GetSize()-1] );
        H5IO::SetEntries1DArray(resDescGroup, "StepNumbers", myStepNums.GetSize()-1,
                                        myStepNums.GetSize()-1,
                                        &myStepNums[myStepNums.GetSize()-1] );
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
          H5IO::SetEntries2DArray( entityGroup, "Real", myStepNums.GetSize()-1,
                                   myStepNums.GetSize()-1, 0, numDofs-1,
                                   &resultVec[pos] );
          pos += numDofs;
        } else {
          Vector<Complex> & resultVec = dynamic_cast<Result<Complex>&>
          (*sol).GetVector();
          Vector<Double> realVec(numDofs), imagVec(numDofs);

          for( UInt i = 0; i < numDofs; i++ ) {
            realVec[i] = resultVec[pos+i].real();
            imagVec[i] = resultVec[pos+i].imag();
          }

          H5IO::SetEntries2DArray( entityGroup, "Real", myStepNums.GetSize()-1,
                                   myStepNums.GetSize()-1, 0, numDofs-1,
                                   &realVec[0] );
          H5IO::SetEntries2DArray( entityGroup, "Imag", myStepNums.GetSize()-1,
                                   myStepNums.GetSize()-1, 0, numDofs-1,
                                   &imagVec[0] );
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
    if(externalFiles_)
    {
      InfoNode* in = info->Get("analysis/output/externalFile");
      try {
        in->Get("name")->SetValue(currStepFile_.getFileName());
        in->Get("size")->SetValue((int) currStepFile_.getFileSize());
        info->ToFile();
      } catch (H5::FileIException &h5ex) {}
      currStepFile_.close();
    }

    if( currMeshStepGroup_.getId() > 0 )
      currMeshStepGroup_.close();

    if( currHistStepGroup_.getId() > 0 )
      currHistStepGroup_.close();

  }

  void SimOutputHDF5::FinishMultiSequenceStep( ) {

    // close groups, which were opened in BeginMultiSequenceStep()
    if( registeredMeshResults_.size() > 0 ) {

      // write attributes containing stepNumber and stepValue
      // of last simulation step
      H5IO::WriteAttribute( currMSMeshGroup_, "LastStepNum", currStep_ );
      H5IO::WriteAttribute( currMSMeshGroup_, "LastStepValue", currStepValue_ );

      currMSMeshGroup_.close();
      meshResultsGroup_.close();
      registeredMeshResults_.clear();
    }

    if( registeredHistResults_.size() > 0 ) {
      // write attributes containing stepNumber and stepValue
      // of last simulation step
      H5IO::WriteAttribute( currMSHistGroup_, "LastStepNum", currStep_ );
      H5IO::WriteAttribute( currMSHistGroup_, "LastStepValue", currStepValue_ );

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

  }

  void SimOutputHDF5::Finalize() {

    // Write file header
    WriteFileInfoHeader();

    // return, if only the grid is to be printed
    if( printGridOnly_ )
      return;

    // return, if no commandLine handler or
    // global root ParaemNode are present
    if( !progOpts || !param )
      return;

    std::vector<std::string> fileNames;
    std::vector<std::string> dataSetNames;
    std::ifstream fin;
    std::ostringstream dumpStr;

    fileNames.push_back(progOpts->GetParamFileStr());
    fileNames.push_back(param->Get("fileFormats")
                        ->Get("materialData")
                        ->Get("file")->AsString());

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

    progOpts->GetVersionString( dumpStr, false );
    WriteStringToUserData( "ProgramStats", dumpStr.str() );
  }

  void SimOutputHDF5::InitModule() {
    std::string pathsep;
    std::string fileName;
    std::ostringstream strBuffer;

    // concatenate output file name
    try {
      fs::create_directory( dirName_ );
      pathsep = fs::path("/").native_directory_string();
    } catch (std::exception &ex) {
      EXCEPTION(ex.what());
    }

    strBuffer << dirName_ << pathsep << fileName_ << ".h5";
    fileName = strBuffer.str();

    // create main file and obtain main group
    try {
      mainFile_ = H5::H5File (fileName, H5F_ACC_TRUNC );
    } H5_CATCH( "Could not create hdf5 file '" << fileName << "' : " );

    mainGroup_ = mainFile_.openGroup( "/" );
    meshGroup_ = mainGroup_.createGroup( "Mesh" );
    mainGroup_.createGroup( "FileInfo" ).close();

    mainGroup_.createGroup( "UserData" ).close();
    mainGroup_.createGroup( "Results" ).close();

  }

  void SimOutputHDF5::WriteGrid() {

    // ensure that grid gets only written once
    if(!gridWritten_)
      InitModule();
    else
      return;

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
      Point p;
      ptGrid_->GetNodeCoordinate(p, i+1);
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
    std::vector< UInt > numElemsOfDim ( 3 );

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
      numElemsOfDim[Elem::GetElemDim(eType)-1]++;
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
    for(UInt i=0; i<3; i++) {
      std::stringstream attrName;
      attrName << "Num" << (i+1) << "DElems";
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
  }


  void SimOutputHDF5::WriteRegions(const H5::Group& meshGroup) {
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
        dims.insert( elems[j]->ptElem->GetDim() );
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
    std::string resultName, unit;
    UInt definedOn, numDofs, entryType, saveBegin, saveEnd, saveInc;
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
      if( !isHistory ) {
        saveBegin = std::max( (UInt)1, meshResultSaveBegin_[resultName]);
        saveEnd = std::min( currMSNumSteps_, meshResultSaveEnd_[resultName] );
        saveInc = meshResultSaveInc_[resultName];
      } else {
        saveBegin = std::max( (UInt) 1, histResultSaveBegin_[resultName] );
        saveEnd = std::min( currMSNumSteps_, histResultSaveEnd_[resultName] );
        saveInc = histResultSaveInc_[resultName];
      }

      // Generate list of entityNames for the current result.
      entityNames.clear();
      for( ; solIt != solEnd; solIt++ ) {
        actResInfo = (*solIt)->GetResultInfo();

        entityNames.
          push_back((*solIt)->GetEntityList()->GetName());
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

        // === Second version: Separate datasets for each entry
        H5::Group actGroup = resGroup.createGroup(resultName );

        H5IO::Write1DArray( actGroup, "DefinedOn", 1, &definedOn, dPropList_ );
        H5IO::Write1DArray( actGroup, "EntityNames", entityNames.size(),
                            &entityNames[0], dPropList_ );
        H5IO::Write1DArray( actGroup, "NumDOFs", 1, &numDofs, dPropList_ );
        H5IO::Write1DArray( actGroup, "DOFNames", resInfo->dofNames.GetSize(),
                            &(resInfo->dofNames[0]), dPropList_ );
        H5IO::Write1DArray( actGroup, "EntryType", 1, &entryType, dPropList_ );
        H5IO::Write1DArray( actGroup, "Unit", 1, &unit, dPropList_ );

        UInt numRealSteps = (UInt) (saveEnd-saveBegin) / saveInc + 1;

        // do not reserve 1D array for StepValues and StepNumbers but initialize with 0
        // otherwise strange memory bugs occur! :(
        StdVector<Double> tmp_double;
        tmp_double.Resize(numRealSteps, 0.0);
        H5IO::Write1DArray( actGroup, "StepValues", numRealSteps, tmp_double.GetPointer(), dPropList_ );

        StdVector<UInt> tmp_uint;
        tmp_uint.Resize(numRealSteps, 0);
        H5IO::Write1DArray( actGroup, "StepNumbers", numRealSteps, tmp_uint.GetPointer(), dPropList_ );

        actGroup.close();


      } H5_CATCH( "Could not write result description for result '"
                  << resultName << "'" );
    }
  }

  void SimOutputHDF5::WriteResults( H5::Group& resultGroup,
                                    Vector<Double>& resultVals,
                                    const UInt numDOFs,
                                    const bool isImag ) {

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
    std::stringstream fName, masterGroup;
    std::string pathsep, fn;

    try {

      // open external file
      pathsep = fs::path("/").native_directory_string();

      fName << fileName_ << "_ms" << currMS_ << "_step"
            << currStep_ << ".h5";
      fn = fName.str();
      fName.str("");
      fName << dirName_ << pathsep << fn;
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
    try {
      userDataGroup = mainGroup_.openGroup("UserData");
    } H5_CATCH( "Can not write meta information to hdf5 file" );

    H5IO::Write1DArray( userDataGroup, dSetName, 1, &str, dPropList_ );
    userDataGroup.close();
  }

} // end of namespace
