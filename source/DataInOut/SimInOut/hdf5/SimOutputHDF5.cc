// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <fstream>
#include <filesystem>
#include <boost/container/flat_set.hpp>
#include <def_cfs_stats.hh>
#include <hdf5.h>
#include <hdf5_hl.h>

#include "DataInOut/SimInOut/hdf5/SimOutputHDF5.hh"
#include "DataInOut/SimInOut/hdf5/CommonHDF5.hh"
#include "DataInOut/ProgramOptions.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "General/Exception.hh"
#include "SimOutputHDF5.hh"
#include "Domain/ElemMapping/Elem.hh"
#include "Utils/Timer.hh"

namespace CoupledField 
{
  DEFINE_LOG(h5Out, "hdf5Out")

  SimOutputHDF5::SimOutputHDF5(std::string fileName, PtrParamNode inputNode, PtrParamNode infoNode, bool isRestart ) :
    SimOutput(fileName, inputNode, infoNode, isRestart) {

    LOG_DBG(h5Out) << "SO fn=" << fileName << " rs=" << isRestart;

    fileName_ = fileName;
    formatName_ = "hdf5";
   
    std::string dirString = "results_" + formatName_; 
    inputNode->GetValue("directory", dirString, ParamNode::PASS );
    dirName_ = dirString; 
    
    capabilities_.insert( MESH );
    capabilities_.insert( MESH_RESULTS );
    capabilities_.insert( HISTORY );
    capabilities_.insert( DATABASE );
    
    // note that for some tests SimState::SetOutputHdf5Writer() calls us with a ParamNode not from a schema, 
    // hence we have no default values set in that case (check e.g. Coil3DExtJ)
    myParam_->GetValue("externalFiles", externalFiles_, ParamNode::PASS);
    myParam_->GetValue("flushSeconds", autoFlushSeconds_, ParamNode::PASS);
    myParam_->GetValue("compressionLevel", compressionLevel_, ParamNode::PASS);
    std::string extString = "cfs"; // sensible default
    inputNode->GetValue("extension", extString, ParamNode::PASS);    
    
    std::string fName = fileName_ + "." + extString;
    currFileName_ = fs::path(dirName_ / fName).string();

    initTimer_ = make_shared<Timer>(timer);
    if(progOpts && progOpts->DoDetailedInfo()) // not for cfsdat
      myInfo_->Get("init/timer")->SetValue(initTimer_);

    // this would disable error messages from the library
    // H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
  }


  SimOutputHDF5::~SimOutputHDF5() {
    CloseFile();
  }

  void SimOutputHDF5::Init( Grid* ptGrid, bool printGridOnly ) {
    initTimer_->Start();
    LOG_DBG(h5Out) << "Init";
    ptGrid_ = ptGrid;
    printGridOnly_ = printGridOnly;
    WriteGrid();
    initTimer_->Stop();
  }

  void SimOutputHDF5::WriteFileInfoHeader() 
  {
    LOG_DBG(h5Out) << "WFIH";
    hid_t infoGroup = -1;
    infoGroup = OpenGroup( mainGroup_, "FileInfo" );

    // write file version
    std::stringstream version;
    version << 0 << "." << 9; // not of any real use
    std::string versionString = version.str();
    WriteSingleDataSet( infoGroup, "Version", versionString);

    // write date / time information
    std::string now = Timer::TimeStampYYYYmmDD();
    WriteSingleDataSet( infoGroup, "Date", now );

    // write creator
    std::string creator = "openCFS " + std::string(CFS_VERSION) + ", " + std::string(CFS_NAME);
    WriteSingleDataSet( infoGroup, "Creator", creator );

    // create group for content
    StdVector<int> content;
    for(const Capability& cap : usedCapabilities_ ) 
      content.Push_back((int) cap);
    
    WriteDataSet1D( infoGroup, "Content", content.GetPointer(), content.GetSize(), compressionLevel_ );
    H5Gclose(infoGroup);
  }


  void SimOutputHDF5::BeginMultiSequenceStep(unsigned int step, BasePDE::AnalysisType type, unsigned int numSteps) 
  {
    LOG_DBG(h5Out) << "BMSS step=" << step << " at=" << BasePDE::analysisType.ToString(type) << " ns=" << numSteps;
    
    currMSNumSteps_ = numSteps;

    // If it does not exist, create Group for Grid / Volume data
    resultsGroup_ = OpenGroup(mainGroup_, "Results");

    // Map analysistype
    std::string analysisType = BasePDE::analysisType.ToString(type);
    currAnalysisType_ = type;

    // Assemble name of multistep
    std::string msName = "MultiStep_" + std::to_string(step);

    // 1) Write group for mesh results
    if( registeredMeshResults_.size() != 0) 
    {
      meshResultsGroup_ = CreateGroup(resultsGroup_, "Mesh", true);
      // write attribute indicating use of external files for simlation
      // steps
      WriteAttribute( meshResultsGroup_, "ExternalFiles", externalFiles_ );

      // create new multistep group.
      currMSMeshGroup_ = CreateGroup(meshResultsGroup_, msName, true);

      // add analysistype and number of steps to group
      WriteAttribute( currMSMeshGroup_, "AnalysisType", analysisType );
      WriteAttribute( currMSMeshGroup_, "LastStepNum", 0u );
      WriteAttribute( currMSMeshGroup_, "LastStepValue", 0.0);

      // add a group for the result description datasets.
      hid_t resultDescGroup = CreateGroup(currMSMeshGroup_, "ResultDescription", true); // use existing

      // write result meta information to file
      WriteResultDescriptions( resultDescGroup, numSteps, false);
      H5Gclose(resultDescGroup);
    }

    // 2) Write group for history results
    if( registeredHistResults_.size() != 0 ) 
    {
      histResultsGroup_ = CreateGroup(resultsGroup_, "History", true);
      currMSHistGroup_ = CreateGroup(histResultsGroup_, msName, true);

      // add analysistype and number of steps to group
      WriteAttribute( currMSHistGroup_, "AnalysisType", analysisType );
      WriteAttribute( currMSHistGroup_, "LastStepNum", 0u  );
      WriteAttribute( currMSHistGroup_, "LastStepValue", 0.0 );
      
      // add a group for the result description datasets.
      hid_t resultDescGroup = CreateGroup(currMSHistGroup_, "ResultDescription", true);

      // write result meta information to file
      WriteResultDescriptions( resultDescGroup, numSteps, true );
      H5Gclose(resultDescGroup);

      // iterate over all results
      for (auto& [name, results] : registeredHistResults_) 
      {
        // create for each result a group within the ms group
        hid_t resultGroup = CreateGroup(currMSHistGroup_, name, true);

        // create subgroup for entitytype
        ResultInfo::EntityUnknownType definedOn = results[0]->GetResultInfo()->definedOn;
        std::string entityString = entityGroupNameEnum.ToString(definedOn);
        hid_t entityTypeGroup = CreateGroup(resultGroup, entityString, true);

        // iterate over all entitylists of result and create sub-subgroup
        for (const shared_ptr<BaseResult>& sol : results) 
        {
          EntityIterator entIt = sol->GetEntityList()->GetIterator();          
          
          for( entIt.Begin(); !entIt.IsEnd(); entIt++ ) 
          {
            // we make sure, that the group exists at the end
            if(H5Lexists(entityTypeGroup, entIt.GetIdString().c_str(), H5P_DEFAULT) > 0)
            {
              if(!isRestart_) // don't warn in the restart case - do nothing otherwise
                WARN("You are trying to add history entity '" << entIt.GetIdString() << "' under group '"
                      << "History/" << msName << "/" << name << "/" << entityString 
                      << "'\nwhich already exists under a different name! Please check your mesh and XML files.");
            }
            else 
            {
              hid_t entityGroup = CreateGroup(entityTypeGroup, entIt.GetIdString());
              H5Gclose(entityGroup);
            }
          }
        }
        H5Gclose(entityTypeGroup);
        H5Gclose(resultGroup);
      }
    }
    currMS_ = step;
  }

  void SimOutputHDF5::RegisterResult( shared_ptr<BaseResult> sol,
                                      unsigned int saveBegin, unsigned int saveInc,
                                      unsigned int saveEnd,
                                      bool isHistory ) {

    std::string resultName = sol->GetResultInfo()->resultName;

    LOG_DBG(h5Out) << "RS sol=" << resultName << " sb=" << saveBegin << " se=" << saveEnd << " inc=" << saveInc << " hist=" << isHistory;

    if( !isHistory ) {
      registeredMeshResults_[resultName].Push_back(sol);
      usedCapabilities_.insert(MESH_RESULTS);
    } else {
      registeredHistResults_[resultName].Push_back(sol);
      usedCapabilities_.insert(HISTORY);
    }
  }

  void SimOutputHDF5::BeginStep( unsigned int stepNum, Double stepVal ) {
    LOG_DBG(h5Out) << "BS num=" << stepNum << " v=" << stepVal;
    currStep_ = stepNum;
    currStepValue_ = stepVal;
  }


  void SimOutputHDF5::AddResult( shared_ptr<BaseResult> sol )
  {
    std::string resultName = sol->GetResultInfo()->resultName;

    // try to determine, if current result is a history or mesh result
    bool isHistory = registeredHistResults_.count(resultName) > 0 
                  && registeredHistResults_[resultName].Contains(sol);
    
    LOG_DBG(h5Out) << "AR sol=" << resultName << " history=" << isHistory; 
                  
    if(isHistory)
      AddHistResult(sol);
    else
      AddMeshResult(sol);
    
    // continuously update the "lastStep" value and number
    // to ensure a consistent hdf5 file throughout the 
    // simulation, even in case of a Ctrl-C action
    if( registeredMeshResults_.size() > 0 ) {
      WriteAttribute( currMSMeshGroup_, "LastStepNum", currStep_ );
      WriteAttribute( currMSMeshGroup_, "LastStepValue", currStepValue_ );
    }

    if( registeredHistResults_.size() > 0 ) {
      WriteAttribute( currMSHistGroup_, "LastStepNum", currStep_ );
      WriteAttribute( currMSHistGroup_, "LastStepValue", currStepValue_ );
    }
  }

  void SimOutputHDF5::AddMeshResult(shared_ptr<BaseResult> sol) 
  {
    LOG_DBG(h5Out) << "AMR: sol=" << sol->GetResultInfo()->resultName;
    
    hid_t resultGroup = -1, subGroup = -1, regionGroup = -1;

    std::string regionName = sol->GetEntityList()->GetName();
    shared_ptr<ResultInfo> resInfo = sol->GetResultInfo();
    std::string resultName = resInfo->resultName;
    unsigned int numDOFs = resInfo->dofNames.GetSize();

    // check, if step group is already open
    if( currMeshStepGroup_ < 0 ) 
    {
      std::string stepName = "Step_" + std::to_string(currStep_);
      
      if(H5Lexists(currMSMeshGroup_, stepName.c_str(), H5P_DEFAULT) > 0)
        currMeshStepGroup_ = OpenGroup(currMSMeshGroup_, stepName);
      else
      {
        currMeshStepGroup_ = CreateGroup(currMSMeshGroup_, stepName);
        WriteAttribute( currMeshStepGroup_, "StepValue", currStepValue_ );

        if(externalFiles_ )
          CreateExternalFile();
      }
    }

    // Add current stepvalue to related group in result description,
    // if not yet present
    StdVector<unsigned int> & myStepNums = meshResultStepNums_[resultName];
    StdVector<Double> & myStepVals = meshResultStepVal_[resultName];
    bool writeStep = myStepNums.GetSize() == 0 || myStepNums.Last() != currStep_;

    if( writeStep ) 
    {
      myStepNums.Push_back( currStep_ );
      myStepVals.Push_back( currStepValue_ );
      
      hid_t msDescGroup = OpenGroup(currMSMeshGroup_, "ResultDescription");
      hid_t resDescGroup = OpenGroup(msDescGroup, resultName);
      H5Gclose(msDescGroup);

      WriteGrowingDataSet1D<Double>(resDescGroup, "StepValues", &myStepVals.Last(), 1);
      unsigned int stepNum = myStepNums.Last();
      WriteGrowingDataSet1D<unsigned int>(resDescGroup, "StepNumbers", &stepNum, 1);
      H5Gclose(resDescGroup);
    }

    // check, if result was already written
    resultGroup = CreateGroup(currMeshStepGroup_, resultName, true); // use existing
    
    // determine, on which type of entity the result is defined
    assert(resInfo->definedOn  == ResultInfo::NODE || resInfo->definedOn == ResultInfo::ELEMENT || resInfo->definedOn == ResultInfo::SURF_ELEM);
    std::string entityString = resInfo->definedOn == ResultInfo::NODE ? "Nodes" : "Elements";

    // try to create regionGroup
    regionGroup = CreateGroup(resultGroup, regionName, true);

    // try to create subgroup for entity
    LOG_DBG2(h5Out) << "AMR: Create subgroup " << entityString << " for result " << resultName << " on region " << regionName << " in step " << currStep_;
    subGroup = CreateGroup(regionGroup, entityString);

    if( sol->GetEntryType() == BaseMatrix::DOUBLE ) 
    {
      Vector<double>& resultVec = dynamic_cast<Result<double>*>(sol.get())->GetVector();
      WriteResults(subGroup, resultVec, numDOFs, false);
    } 
    else 
    {
      Vector<Complex> & resultVec = dynamic_cast<Result<Complex>*>(sol.get())->GetVector();

      Vector<double> realVec(resultVec.GetSize());
      Vector<double> imagVec(resultVec.GetSize());
      for(unsigned int i = 0; i < resultVec.GetSize(); i++) 
      {
        realVec[i] = resultVec[i].real();
        imagVec[i] = resultVec[i].imag();
      }
      WriteResults(subGroup, realVec, numDOFs, false);
      WriteResults(subGroup, imagVec, numDOFs, true);
    }

    // close groups
    H5Gclose(subGroup);
    H5Gclose(regionGroup);
    H5Gclose(resultGroup);
  }

  void SimOutputHDF5::AddHistResult( shared_ptr<BaseResult> sol) 
  {
    shared_ptr<ResultInfo> resInfo = sol->GetResultInfo();
    std::string resultName = resInfo->resultName;
    unsigned int numDofs = resInfo->dofNames.GetSize();
    std::string entityString = entityGroupNameEnum.ToString(resInfo->definedOn);

    LOG_DBG(h5Out) << "AHR sol=" << resultName;

    // Add current stepvalue to related group in result description,
    // if not yet present
    StdVector<unsigned int>& myStepNums = histResultStepNums_[resultName];
    StdVector<Double>& myStepVals = histResultStepVal_[resultName];

    bool writeStep = myStepNums.GetSize() == 0 || myStepNums.Last() != currStep_;
    
    if( writeStep ) 
    {
      myStepNums.Push_back( currStep_ );
      myStepVals.Push_back( currStepValue_ );

      hid_t msDescGroup = OpenGroup(currMSHistGroup_, "ResultDescription");
      hid_t resDescGroup = OpenGroup(msDescGroup, resultName);
      H5Gclose(msDescGroup);

      WriteGrowingDataSet1D<Double>(resDescGroup, "StepValues", &myStepVals.Last(), 1);
      unsigned int stepNum = myStepNums.Last();
      WriteGrowingDataSet1D<unsigned int>(resDescGroup, "StepNumbers", &stepNum, 1);
      H5Gclose(resDescGroup);
    }
    // ---------------------
    //  Write result itself
    // ---------------------
    hid_t resultGroup = OpenGroup(currMSHistGroup_, resultName);
    hid_t entityTypeGroup = OpenGroup(resultGroup, entityString);


    // iterate over all entities in this list
    EntityIterator entIt = sol->GetEntityList()->GetIterator();
    unsigned int pos = 0;
    for( entIt.Begin(); !entIt.IsEnd(); entIt++ ) 
    {
      hid_t entityGroup = OpenGroup(entityTypeGroup, entIt.GetIdString());

      if( sol->GetEntryType() == BaseMatrix::DOUBLE ){
        Vector<Double> & resultVec = dynamic_cast<Result<Double>&>(*sol).GetVector();
        WriteGrowingDataSet2D<Double>(entityGroup, "Real", 1, numDofs, &resultVec[pos]);
        pos += numDofs;
      } else {
        Vector<Complex> & resultVec = dynamic_cast<Result<Complex>&>(*sol).GetVector();
        Vector<Double> realVec(numDofs), imagVec(numDofs);

        for( unsigned int i = 0; i < numDofs; i++ ) {
          realVec[i] = resultVec[pos+i].real();
          imagVec[i] = resultVec[pos+i].imag();
        }

        WriteGrowingDataSet2D<Double>(entityGroup, "Real", 1, numDofs, realVec.GetPointer());
        WriteGrowingDataSet2D<Double>(entityGroup, "Imag", 1, numDofs, imagVec.GetPointer());
        pos += numDofs;
      }
      H5Gclose(entityGroup);
    }
    H5Gclose(entityTypeGroup);
    H5Gclose(resultGroup);
  }


  void SimOutputHDF5::FinishStep( )
  {
    LOG_DBG(h5Out) << "FS";

    if(externalFiles_ && myInfo_)
    {
      PtrParamNode in = myInfo_->Get("analysis/output/externalFile");

      std::string fn(H5Fget_name(currStepFile_, nullptr, 0), '\0');
      H5Fget_name(currStepFile_, fn.data(), fn.size() + 1);
      in->Get("name")->SetValue(fn);

      hsize_t _sz = 0; 
      H5Fget_filesize(currStepFile_, &_sz);
      in->Get("size")->SetValue((int)_sz);
    }

    // we close everything here, so that the file is in a consistent state
    if (currStepFile_ >= 0) 
      H5Fclose(currStepFile_); 

    if(currMeshStepGroup_ >= 0)
      H5Gclose(currMeshStepGroup_); 

    if(currHistStepGroup_ >= 0) 
      H5Gclose(currHistStepGroup_); 

    currStepFile_ = -1;   
    currMeshStepGroup_ = -1; 
    currHistStepGroup_ = -1;       

    // flush file to allow to read by ParaView -> don't do too often, it slows down a lot
    // TODO: Change to Single Writer Multiple Reader with hdf5 1.10
    H5Fflush(mainFile_, H5F_SCOPE_GLOBAL);
  }

  void SimOutputHDF5::FinishMultiSequenceStep() 
  {
    LOG_DBG(h5Out) << "FMSS";
    // close groups, which were opened in BeginMultiSequenceStep()
    if( registeredMeshResults_.size() > 0 ) {
      H5Gclose(currMSMeshGroup_);  currMSMeshGroup_ = -1;
      H5Gclose(meshResultsGroup_); meshResultsGroup_ = -1;
      registeredMeshResults_.clear();
    }

    if( registeredHistResults_.size() > 0 ) {
      H5Gclose(currMSHistGroup_);  currMSHistGroup_ = -1;
      H5Gclose(histResultsGroup_); histResultsGroup_ = -1;
      registeredHistResults_.clear();
    }
    H5Gclose(resultsGroup_); resultsGroup_ = -1;

    // reset all data per sequence step
    meshResultStepNums_.clear();
    histResultStepNums_.clear();
    meshResultStepVal_.clear();
    histResultStepVal_.clear();

    AutoFlush();
  }

  void SimOutputHDF5::Finalize() {
    LOG_DBG(h5Out) << "Finalize";
    // return in case of an already restarted file
    if (isRestart_)
      return;
    
    // Write file header
    WriteFileInfoHeader();

    // return, if only the grid is to be printed
    if( printGridOnly_ )
      return;

    // return, if no commandLine handler or
    // global root ParaemNode are present
    if( !progOpts || !myParam_ )
      return;

    StdVector<std::string> fileNames;
    fileNames.push_back(progOpts->GetParamFileStr());   
    fileNames.push_back(myParam_->GetRoot()->Get("fileFormats")->Get("materialData")->GetAsFilePath("file"));
    StdVector<std::string> dataSetNames = {"ParameterFile", "MaterialFile"};

    for(unsigned int i=0; i<fileNames.GetSize(); i++)
    {
      std::ifstream fin( fileNames[i].c_str(), std::ios::binary );

      if(fin.fail())
        throw Exception("Cannot open file '" + fileNames[i] + "' to dump into HDF5!");

      std::string str;        
      str.resize(fs::file_size(fileNames[i]));
      fin.read(str.data(), str.size());        

      WriteStringToUserData(dataSetNames[i], str);
    }

    std::ostringstream dumpStr;
    progOpts->PrintVersion( dumpStr, false );
    WriteStringToUserData( "ProgramStats", dumpStr.str() );
  }

  void SimOutputHDF5::InitModule() 
  {
    LOG_DBG(h5Out) << "IM";
    if( isInitialized_)
      return;

    // concatenate output file name
    fs::create_directory( dirName_ ); // throws std::exception which is caught in CFS.cc
    
    // In case of re-start, we simply append information
    bool truncate = !isRestart_;
    OpenFile(truncate);
    isInitialized_ = true;
  }
  
  void SimOutputHDF5::OpenFile(bool truncate)
  {
    LOG_DBG(h5Out) << "OF truncate=" << truncate;
    // create main file and obtain main group
    if (truncate)
      mainFile_ = H5Fcreate(currFileName_.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    else
      mainFile_ = H5Fopen(currFileName_.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    
    if (mainFile_ < 0)
      throw Exception("Could not open/create hdf5 file '" + currFileName_ + "'");

    mainGroup_ = H5Gopen2(mainFile_, "/", H5P_DEFAULT);
    if(truncate)
    {
      meshGroup_ = CreateGroup(mainGroup_, "Mesh");
      H5Gclose(CreateGroup(mainGroup_, "FileInfo"));
      H5Gclose(CreateGroup(mainGroup_, "UserData"));
      H5Gclose(CreateGroup(mainGroup_, "Results"));
    }
    else
      meshGroup_ = OpenGroup(mainGroup_, "Mesh"); 
  }
  
  void ReportOpenHDF5Objects(hid_t fid)
  {
    assert(fid >= 0);
    const unsigned int types = H5F_OBJ_DATASET | H5F_OBJ_GROUP | H5F_OBJ_DATATYPE | H5F_OBJ_ATTR;
    ssize_t nOpen = H5Fget_obj_count(fid, types);
    if (nOpen <= 0)
      return;

    std::string fn(H5Fget_name(fid, nullptr, 0), '\0');
    H5Fget_name(fid, fn.data(), fn.size() + 1);
    std::string msg = std::to_string(nOpen) + " object(s) were open on closing " + fn + ": ";

    std::vector<hid_t> ids(nOpen);
    H5Fget_obj_ids(fid, types, nOpen, ids.data());
    for (hid_t obj : ids)
    {
      std::string typeName;
      switch (H5Iget_type(obj))
      {
        case H5I_GROUP:    typeName = "group";    break;
        case H5I_DATASET:  typeName = "dataset";  break;
        case H5I_DATATYPE: typeName = "datatype"; break;
        case H5I_ATTR:     typeName = "attr";     break;
        default:           typeName = "unknown";  break;
      }
      ssize_t nameLen = H5Iget_name(obj, nullptr, 0);
      std::string name(nameLen, '\0');
      H5Iget_name(obj, name.data(), nameLen + 1);
      msg += " " + typeName + ": " + name;
    }
    WARN(msg);
  }

  void SimOutputHDF5::CloseFile()
  {
    LOG_DBG(h5Out) << "CF";
    if(currStepFile_ >= 0) {
      ReportOpenHDF5Objects(currStepFile_);
      H5Fclose(currStepFile_);
      currStepFile_ = -1;
    }

    // check, if any group is open at all
    if(resultsGroup_ >= 0) 
      H5Gclose(resultsGroup_); 

    if(meshGroup_ >= 0)
      H5Gclose(meshGroup_); 

    if(dbGroup_ >= 0) 
      H5Gclose(dbGroup_); 
    
    if(currMsDbGroup_ >= 0) 
      H5Gclose(currMsDbGroup_); 
    
    if(mainGroup_ >= 0) 
      H5Gclose(mainGroup_); 

    resultsGroup_ = -1; 
    meshGroup_ = -1;
    dbGroup_ = -1;
    currMsDbGroup_ = -1;
    mainGroup_ = -1;

    // check for open groups, datasets etc.
    if (mainFile_ >= 0) {
      ReportOpenHDF5Objects(mainFile_);
      H5Fclose(mainFile_);
      mainFile_ = -1;
    }
  }

  void SimOutputHDF5::WriteGrid() 
  {
    LOG_DBG(h5Out) << "WG";
    
    // ensure that grid gets only written once
    if(gridWritten_)  
      return; 

    lastFlush_ = std::chrono::steady_clock::now();

    InitModule();      
 
    // in case of restart, we do not re-write the grid!
    if( isRestart_) {
      gridWritten_ = true;
      return;
    }
    
    // write the dimension of the grid.
    WriteAttribute( meshGroup_, "Dimension", ptGrid_->GetDim());

    // ================
    //  Node Locations
    // ================
    hid_t nodeGroup = CreateGroup( meshGroup_, "Nodes" );
     
    unsigned int nNodes = ptGrid_->GetNumNodes();
    WriteAttribute( nodeGroup, "NumNodes", nNodes );

    // collect all nodal coordinates
    StdVector<double> locs(nNodes * 3);
    Vector<double> p(3);
    for(unsigned int i = 0; i < nNodes; i++) {
      ptGrid_->GetNodeCoordinate3D(p, i+1);
      std::copy(p.GetPointer(), p.GetPointer() + 3, locs.GetPointer() + i*3);
    }

    WriteDataSet2D( nodeGroup, "Coordinates", nNodes, 3, locs.GetPointer(), compressionLevel_ );

    H5Gclose(nodeGroup); 
    nodeGroup = -1;


    // =====================
    //  Element definitions
    // =====================
    hid_t elemGroup = CreateGroup(meshGroup_, "Elements");
    unsigned int nElems = ptGrid_->GetNumElems();
    unsigned int maxNumNodes = ptGrid_->GetMaxNumNodesPerElem();
    StdVector<unsigned int> connect (nElems * maxNumNodes);
    StdVector<unsigned int> elConnect (maxNumNodes);
    StdVector<int> feTypes (nElems);
    StdVector<unsigned int> numElemsOfDim (3);

    // Fill connectivity array
    connect.Init(0);

    unsigned int offset;
    Elem::FEType eType;
    RegionIdType region;

    // iterate over all elements and fill hdf5 contiguous data arrays
    for( unsigned int i = 0; i < nElems; i++ ) 
    {
      elConnect.Resize(maxNumNodes,0);
      ptGrid_->GetElemData( i+1, eType, region, elConnect.GetPointer() );

      numElemsOfDim[(Elem::shapes[eType].dim)-1]++;
      feTypes[i] = eType;

      // insert connectivity into global array
      offset = i * maxNumNodes;
      for( unsigned int j = 0; j < elConnect.GetSize(); j++ ) 
        connect[offset + j] = elConnect[j];
    }

    // write connectivity
    WriteDataSet2D(elemGroup, "Connectivity", nElems, maxNumNodes, connect.GetPointer(), compressionLevel_);

    // write element types
    WriteDataSet1D(elemGroup, "Types", feTypes.GetPointer(), nElems, compressionLevel_);

    // ==========================
    //  Grid Meta Information
    // ==========================

    WriteAttribute( elemGroup, "NumElems", nElems );

    WriteAttribute( elemGroup, "QuadraticElems", (int) ptGrid_->IsQuadratic());

    // number of elements per dimension
    for(unsigned int i=0; i<3; i++) {
      std::string attrName = "Num" + std::to_string(i+1) + "DElems";
      WriteAttribute( elemGroup, attrName, numElemsOfDim[i] );
    }

    // number of elements per type
    unsigned int numElemTypes = Elem::feType.map.size();
    for(unsigned int i=0; i<numElemTypes; i++) {
      std::string attrName = "Num_" + Elem::feType.ToString((Elem::FEType)i);
      WriteAttribute( elemGroup, attrName, ptGrid_->GetNumElemOfType((Elem::FEType)i) );
    }

    // close element group
    H5Gclose(elemGroup); 
    elemGroup = -1;

    // ============================================
    //  Write Regions, Node Groups, Element Groups
    // ============================================

    WriteRegions( meshGroup_ );

    // create new group for entity groups
    hid_t groupsGroup = CreateGroup(meshGroup_, "Groups");

    WriteNodeGroups( groupsGroup );
    WriteElemGroups( groupsGroup );
    H5Gclose(groupsGroup);

    gridWritten_ = true;
    usedCapabilities_.insert(MESH);

    // close meshGroup
    H5Gclose(meshGroup_); 
    meshGroup_ = -1;

    AutoFlush(); // flush only of writing this grid took really long
  }


  void SimOutputHDF5::WriteRegions(hid_t meshGroup) 
  {
    LOG_DBG(h5Out) << "WR";
    StdVector< std::string > regionNames;
    StdVector< unsigned int > regionDims;
    StdVector< StdVector<unsigned int> > regionElems;
    StdVector< StdVector<unsigned int> > regionNodes;
    StdVector<Elem*> elems;
    StdVector<RegionIdType> surfRegionIds, volRegionIds;

    // create region group
    hid_t regionListGroup = CreateGroup(meshGroup, "Regions");

    unsigned int numRegions = ptGrid_->GetNumRegions();
    if(numRegions == 0)
      return;

    unsigned int dim = ptGrid_->GetDim();
    ptGrid_->GetVolRegionIds(volRegionIds);
    ptGrid_->GetSurfRegionIds(surfRegionIds);
    ptGrid_->GetRegionNames(regionNames);
    regionDims.Resize(numRegions);
    regionElems.Resize(numRegions);
    regionNodes.Resize(numRegions);
    regionDims.Resize(numRegions);

    // obtain nodes and elements per surface region
    int idx = -1;
    for(unsigned int i=0, n=surfRegionIds.GetSize(); i<n; i++) {
      idx = surfRegionIds[i];
      regionDims[idx] = dim-1;

      ptGrid_->GetElems(elems, idx);
      unsigned int nElems = elems.GetSize();
      regionElems[idx].Resize(nElems);
      for(unsigned int j=0; j<nElems; j++) 
        regionElems[idx][j] = elems[j]->elemNum;
      
      ptGrid_->GetNodesByRegion(regionNodes[idx], idx);

      // determine dimensionality
      regionDims[idx] = ptGrid_->GetDim() == 3 ? 2 : 1;
    }

    // obtain nodes and elements per volume region
    for(unsigned int i=0, n=volRegionIds.GetSize(); i<n; i++) {
      idx = volRegionIds[i];
      regionDims[idx] = dim;

      ptGrid_->GetElems(elems, idx);
      unsigned int nElems = elems.GetSize();
      regionElems[idx].Resize(nElems);
      for(unsigned int j=0; j<nElems; j++) {
        regionElems[idx][j] = elems[j]->elemNum;
      }

      ptGrid_->GetNodesByRegion(regionNodes[idx], idx);

      // determine dimensionality
      regionDims[idx] = ptGrid_->GetDim();
    }


    // loop over regions and write out nodes and elements
    for(unsigned int i = 0; i < numRegions; i++)
    {
      // create new region group
      hid_t actRegionGroup = CreateGroup(regionListGroup, regionNames[i] );
      LOG_DBG2(h5Out) << "WR: write region " << regionNames[i];
      WriteAttribute( actRegionGroup, "Dimension", regionDims[i] );

      // create new node group
      WriteDataSet1D<unsigned int>( actRegionGroup, "Nodes", regionNodes[i].GetPointer(), regionNodes[i].GetSize(), compressionLevel_);

      // create new element group
      WriteDataSet1D( actRegionGroup, "Elements", (const int*) regionElems[i].GetPointer(), regionElems[i].GetSize(), compressionLevel_);

      // create new face group
      // .. to be implemented ..

      // create new edge group
      // .. to be implemented ..

      // close current region group
      H5Gclose(actRegionGroup);
    }

    // close regionlist group
    H5Gclose(regionListGroup);

  }

  void SimOutputHDF5::WriteNodeGroups(hid_t meshGroup) 
  {
    LOG_DBG(h5Out) << "WNG";
    StdVector<unsigned int> nodes;
    StdVector<std::string> nodeNames;

    // obtain list with names of nodes
    ptGrid_->GetListNodeNames(nodeNames);
    unsigned int numNodeGroups = nodeNames.GetSize();

    for(unsigned int i = 0; i < numNodeGroups; i++ ) {
      ptGrid_->GetNodesByName(nodes, nodeNames[i]);

      // try to open group with given name
      hid_t myGroup = CreateGroup( meshGroup, nodeNames[i], true ); // true -> use old if exits     
      WriteAttribute( myGroup, "Dimension", (int) 0 );
      WriteDataSet1D( myGroup, "Nodes", nodes.GetPointer(), nodes.GetSize(), compressionLevel_);

      // close nodes array of current group
      H5Gclose(myGroup); myGroup = -1;
    }
  }

  void SimOutputHDF5::WriteElemGroups(hid_t meshGroup) 
  {
    LOG_DBG(h5Out) << "WEG";
    StdVector<unsigned int> elemNums, elemNodes;
    StdVector<Elem*> elems;
    StdVector<std::string> elemNames;
    
    // find unique nodes by element group -> keep order for test cases
    boost::container::flat_set<unsigned int> nodeSet;

    // obtain list with names of elements
    ptGrid_->GetListElemNames(elemNames);
    unsigned int numElemGroups = elemNames.GetSize();

    for(unsigned int i = 0; i < numElemGroups; i++ ) 
    {
      // get all elements by name
      ptGrid_->GetElemsByName(elems, elemNames[i]);
      elemNums.Resize( elems.GetSize() );
      nodeSet.clear();
      nodeSet.reserve(elems.GetSize() * 10); // rough upper estimate
      
      unsigned int dim = Elem::shapes[elems[0]->type].dim;
      for( unsigned int j = 0; j < elems.GetSize(); j++ ) 
      {
        if (Elem::shapes[elems[j]->type].dim != dim)
           EXCEPTION( "Element group '" << elemNames[i] << "' contains elements of different dimensions" );
        elemNums[j] = elems[j]->elemNum;
        nodeSet.insert( elems[j]->connect.Begin(), elems[j]->connect.End() );
      }
      
      hid_t myGroup = CreateGroup( meshGroup, elemNames[i], true ); // use existing
      WriteAttribute( myGroup, "Dimension", dim);
      WriteDataSet1D( myGroup, "Elements", elemNums.GetPointer(), elemNums.GetSize(), compressionLevel_);

      // Write nodes of element group
      elemNodes.Resize( nodeSet.size()); // we need a contiguous array
      std::copy( nodeSet.begin(), nodeSet.end(), elemNodes.Begin() );
      WriteDataSet1D( myGroup, "Nodes", elemNodes.GetPointer(), elemNodes.GetSize(), compressionLevel_);

      // close nodes array of current group
      H5Gclose(myGroup); 
    }
  }

  void SimOutputHDF5::WriteResultDescriptions( hid_t resGroup, unsigned int numSteps, bool isHistory )  
  {
    LOG_DBG(h5Out) << "WRD num=" << numSteps << " h=" << isHistory;
    std::string unit;
    unsigned int definedOn, numDofs, entryType; 
    //unsigned int saveBegin, saveEnd, saveInc;
    std::vector<std::string> entityNames;
    shared_ptr<ResultInfo> resInfo;

    const auto& reg = isHistory ? registeredHistResults_ : registeredMeshResults_;
    for (const auto& [resultName, solVec] : reg) 
    {
      resInfo = solVec[0]->GetResultInfo();
      numDofs = resInfo->dofNames.GetSize();
      unit = resInfo->unit;
      definedOn = MapEntityTypeToInt(resInfo->definedOn);
      entryType = MapEntryTypeToInt(resInfo->entryType);

      // Generate list of entityNames for the current result.
      entityNames.clear();
      for(const auto& sol : solVec) 
        entityNames.push_back(sol->GetEntityList()->GetName());

      // Generate compound datatype

      // Check, if the group for the group for the result has to be created.
      // which is either the case if the simulation is not restarted or
      // if the simulation is restarted, but the restarted sequenceStep
      // is not the current one.
      bool exists = H5Lexists(resGroup, resultName.c_str(), H5P_DEFAULT) > 0;
      if( !isRestart_ || (isRestart_ && !exists) )
      {
        hid_t actGroup = CreateGroup(resGroup, resultName );

        WriteSingleDataSet( actGroup, "DefinedOn", definedOn );
        WriteStringArray( actGroup, "EntityNames", entityNames );
        WriteSingleDataSet( actGroup, "NumDOFs", numDofs );
        { // DOFNames: convert StdVector<std::string> to std::vector<std::string>
          std::vector<std::string> dofNamesVec;
          for (const std::string& name : resInfo->dofNames)
            dofNamesVec.push_back(name);
          WriteStringArray( actGroup, "DOFNames", dofNamesVec );
        }
        WriteSingleDataSet( actGroup, "EntryType", entryType );
        WriteSingleDataSet( actGroup, "Unit", unit );
        // In order to write a valid entry, we also set the initial stepNumber/
        // stepValues array
        WriteGrowingDataSet1D<unsigned int>(actGroup, "StepNumbers"); // no data yet!
        WriteGrowingDataSet1D<double>(actGroup, "StepValues");
        
        H5Gclose(actGroup);
      }
      else
      { // general restart case
        hid_t actGroup = OpenGroup(resGroup, resultName );
        
        // Obtain already written stepValue and stepNumbers and
        // store them in the meshResultStepNums_, meshResultStepVal_.
        // These array may not be created yet.
        if( H5Lexists(actGroup, "StepNumbers", H5P_DEFAULT) > 0 &&
            H5Lexists(actGroup, "StepValues",  H5P_DEFAULT) > 0 )
        {
          if( isHistory ) {
            ReadArray(actGroup, "StepNumbers", histResultStepNums_[resultName]);
            ReadArray(actGroup, "StepValues",  histResultStepVal_[resultName]);
          } else {
            // Obtain already written stepValue and stepNumbers and
            // store them in the meshResultStepNums_, meshResultStepVal_
            ReadArray(actGroup, "StepNumbers", meshResultStepNums_[resultName]);
            ReadArray(actGroup, "StepValues",  meshResultStepVal_[resultName]);
          }
        }
        H5Gclose(actGroup);
      }
    } //loop: registered mesh / history results
  }

  void SimOutputHDF5::WriteResults( hid_t resultGroup, Vector<double>& resultVals, const unsigned int numDOFs, const bool isImag ) 
  {
    LOG_DBG(h5Out) << "WR n=" << resultVals.GetSize() << " nd=" << numDOFs << " imag=" << isImag;
    // create dataset with related name
    std::string name = isImag ? "Imag" : "Real";

    unsigned int numEntities = resultVals.GetSize() / numDOFs;
    
    WriteDataSet2D(resultGroup, name, numEntities, numDOFs, resultVals.GetPointer(), compressionLevel_);
  }

  void SimOutputHDF5::CreateExternalFile() 
  {
    LOG_DBG(h5Out) << "CEF";

    // open external file
    string fn = fileName_ + "_ms" + std::to_string(currMS_) + "_step" + std::to_string(currStep_) + ".cfs";
    std::string fullPath = (dirName_ / fn).string(); // std::filesystem operator magic

    currStepFile_ = H5Fcreate(fullPath.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // Write reference to external file to main file
    WriteAttribute( currMeshStepGroup_, "ExtHDF5FileName", fn );

    // set current step group to external file
    H5Gclose(currMeshStepGroup_);
    currMeshStepGroup_ = H5Gopen2(currStepFile_, "/", H5P_DEFAULT);

    // Store reference to master file in external file
    WriteAttribute( currMeshStepGroup_, "MasterHDF5FileName", fileName_ + ".cfs" );

    // Store reference to main group in external file
    std::string masterGroup = "/Results/Mesh/MultiStep_" + std::to_string(currMS_) + "/Step_" + std::to_string(currStep_);
    WriteAttribute( currMeshStepGroup_, "MasterGroup", masterGroup );    
  }

  void SimOutputHDF5::WriteStringToUserData(const std::string& dSetName, const std::string& str) 
  {
    // If it does not exist, create Group for Data.
    hid_t userDataGroup = CreateGroup(mainGroup_, "UserData", true);
    WriteSingleDataSet( userDataGroup, dSetName, str );
    H5Gclose(userDataGroup);
  }
  
  
  // ------------------------------------------------------------------------
  //  DATABASE SECTION
  // ------------------------------------------------------------------------
  void SimOutputHDF5::DB_Init() 
  {
    InitModule();
    useDataBase_ = true;

    // insert capability
    usedCapabilities_.insert(DATABASE);

    if( isRestart_) 
      dbGroup_ = OpenGroup(mainGroup_, "DataBase");
    else 
    {
      // add "DataBase" group
      dbGroup_ = CreateGroup(mainGroup_, "DataBase");
      H5Gclose(CreateGroup(dbGroup_, "MultiSteps"));
    }
  }
  
  void SimOutputHDF5::DB_WriteXmlFiles( fs::path simFile, fs::path matFile ) 
  {
    if( isRestart_)
      return;
    
    std::ifstream fin;
    hid_t extFiles = CreateGroup(dbGroup_, "InputFiles");
    
    // open external Files
    StdVector<fs::path> filePaths = {simFile, matFile};
    StdVector<std::string> setNames = {"ParameterFile", "MaterialFile"};

    for(unsigned int i = 0; i < filePaths.GetSize(); i++)
    {
      fin.open( filePaths[i], std::ios::binary );

      if(fin.fail())
        throw Exception("Cannot open external file '" + filePaths[i].string() + "' to dump into HDF5!");

      // seek to the end of the file
      fin.seekg (0, std::ios::end);
      unsigned int numBytes = fin.tellg();
      fin.seekg (0, std::ios::beg);

      std::string str;
      str.resize(numBytes);
      fin.read(&str[0], numBytes);
      fin.close();
      
      // now write out string
      WriteSingleDataSet( extFiles, setNames[i], str );
    }
    H5Gclose(extFiles);
  }

  void SimOutputHDF5::DB_BeginMultiSequenceStep( unsigned int step, BasePDE::AnalysisType type ) 
  {
    currMS_ = step;
    currAnalysisType_ = type;
    std::string stepStr = std::to_string(step);
    hid_t msGroup = OpenGroup(dbGroup_, "MultiSteps");

    if(currMsDbGroup_ >= 0)  
      H5Gclose(currMsDbGroup_);
    currMsDbGroup_ = CreateGroup(msGroup, stepStr, true); // use existing
    H5Gclose(msGroup);
    
    // set attributes for number of steps and analysis types
    std::string analysisType = BasePDE::analysisType.ToString(type);
    WriteAttribute( currMsDbGroup_, "AnalysisType", analysisType );
    WriteAttribute( currMsDbGroup_, "AccTime", 0.0 );
    //! Only write attribute, if it not exists yet
    if( H5Aexists(currMsDbGroup_, "Completed") == 0 ) {
      WriteAttribute( currMsDbGroup_, "Completed", (int)false );
    }
  }

  void SimOutputHDF5::DB_BeginStep( unsigned int stepNum, Double stepVal ) 
  {
    currStepDb_ = stepNum;
    currStepValueDb_ = stepVal;
    DB_BeginMultiSequenceStep(currMS_, currAnalysisType_ );
  }

  void SimOutputHDF5::DB_WriteFeFunction( const std::string& pdeCplName,
                                          const std::string& fctName,
                                          SingleVector* coefs ) {


    hid_t physGroup = -1, fctGroup = -1, stepGroup = -1;
    // Create / open group for physics
    physGroup = CreateGroup(currMsDbGroup_, pdeCplName, true); // use existing

    // Create / open group for FeFunction
    fctGroup = CreateGroup(physGroup, fctName, true); // use existing
    
    // Create group for current step
    std::string stepStr = std::to_string(currStepDb_);
    stepGroup = CreateGroup(fctGroup, stepStr);
    WriteAttribute( stepGroup, "StepValue", currStepValueDb_ );
    
    // Now write coefficients for fefunction
    if( coefs->GetEntryType() == BaseMatrix::DOUBLE ) {
      Vector<Double> & rVec = dynamic_cast<Vector<Double>& >(*coefs);
      WriteResults( stepGroup, rVec, 1, false );
      
    } else {
      Vector<Complex> & cVec = dynamic_cast<Vector<Complex>& > (*coefs);
      Vector<Double> realVec, imagVec;
      realVec.Resize( cVec.GetSize() );
      imagVec.Resize( cVec.GetSize() );

      for(unsigned int i = 0; i < cVec.GetSize(); i++) {
        realVec[i] = cVec[i].real();
        imagVec[i] = cVec[i].imag();
      }
      WriteResults(stepGroup, realVec, 1, false);
      WriteResults(stepGroup, imagVec, 1, true);

    }
    
    // now close groups
    H5Gclose(stepGroup);
    H5Gclose(fctGroup);
    H5Gclose(physGroup);
  }
  
  void SimOutputHDF5::DB_FinishMultiSequenceStep( bool completed, double accTime ) 
  {
    WriteAttribute( currMsDbGroup_, "Completed", (int)completed );
    WriteAttribute( currMsDbGroup_, "AccTime", accTime );
  }
     
void SimOutputHDF5::AutoFlush() 
{
  auto now = std::chrono::steady_clock::now();
  LOG_DBG(h5Out) << "AF: time since last event: " << std::chrono::duration_cast<std::chrono::seconds>(now - lastFlush_).count() << "s";
  if(now - lastFlush_ > std::chrono::duration<double>(autoFlushSeconds_))
  {
    // flush file to allow to read by ParaView -> don't do too often, it slows down a lot
    // TODO: Change to Single Writer Multiple Reader with hdf5 1.10
    H5Fflush(mainFile_, H5F_SCOPE_GLOBAL);
    lastFlush_ = now;
    LOG_DBG(h5Out) << "AF: flushed file at " << Timer::TimeStamp();
  }
}

} // end of namespace
