// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/exception.hpp>
namespace fs=boost::filesystem;

#include <def_cfs_stats.hh>
#include <def_use_hdf5.hh>

#include "General/exception.hh"
#include "Settings.hh"
#include "General/environment.hh"
#include "DataInOut/SimInOut/hdf5/hdf5io.hh"
#include "cplreader/CouplingHandler.hh"


#include "OutputWriter_HDF5.hh"

#define H5_EXCEPTION(STR, EX)                   \
  EXCEPTION( STR, EX.getCDetailMsg() );

#define H5_CATCH( STR )                                 \
  catch (H5::Exception& h5Ex ) {                        \
    EXCEPTION( STR << ":\n" << h5Ex.getCDetailMsg() );  \
  }

namespace CoupledField
{

  OutputWriter_HDF5::OutputWriter_HDF5() {
    actStepNum_ = 0;
    actStepValue_ = 0;
  }

  OutputWriter_HDF5::~OutputWriter_HDF5() {
    Settings& settings = Settings::Instance();    

    // check for open groups, datasets etc. in current step file
    if(currStepFile_.getLocId() > 0) {
      if (currStepFile_.getObjCount( H5F_OBJ_DATASET |
                                     H5F_OBJ_GROUP |
                                     H5F_OBJ_DATATYPE | H5F_OBJ_ATTR) > 0 ) {
        if(settings.GetInt("verbose"))
          std::cerr << "There are still objects open in the hdf5 file "
                    << currStepFile_.getFileName() << "\n\n";
        H5IO::CheckOpenObjects(currStepFile_, settings.GetInt("verbose") != 0);
      }
      
      currStepFile_.close();
    }

    // check, if any group is open at all
    if( mainGroup_.getLocId() > 0 )
      mainGroup_.close();

    // check for open groups, datasets etc.
    if (mainGroup_.getLocId() > 0 &&
        mainFile_.getObjCount( H5F_OBJ_DATASET |
                               H5F_OBJ_GROUP |
                               H5F_OBJ_DATATYPE | H5F_OBJ_ATTR) > 0 ) {
      if(settings.GetInt("verbose"))
        std::cerr << "There are still objects open in the hdf5 file "
                  << mainFile_.getFileName() << "\n\n";
      H5IO::CheckOpenObjects(mainFile_, settings.GetInt("verbose") != 0);
    }

    mainFile_.close();

    if(!doneWriting_) {
      try {
        fs::remove( mainFileName_ );
      } catch (std::exception &ex) {
        EXCEPTION("An error occured while trying to delete incomplete HDF5 file '"
                  << mainFileName_ << "'\n"
                  << ex.what());
      }
    }
  }

  void OutputWriter_HDF5::Init(int argc, char** argv,
                               const std::string& outputPath)
  {
    hdf5DirName_ = outputPath;

    InitHDF5();
    WriteFileInfoHeader();
  }

  void OutputWriter_HDF5::WriteNodalCoords(const std::vector<Double>& coords) {
    // write the dimension of the grid.
    H5IO::WriteAttribute( meshGroup_, "Dimension", dim_ );

    // ================
    //  Node Locations
    // ================
    UInt nNodes = coords.size()/3;
    H5::Group nodeGroup;
    try {
      nodeGroup = meshGroup_.createGroup( "Nodes" );
    } H5_CATCH( "Could not create node group" );

    H5IO::WriteAttribute( nodeGroup, "NumNodes", nNodes );

    H5IO::Write2DArray( nodeGroup, "Coordinates", nNodes,
                        3, &coords[0], dPropList_ );

    nodeGroup.close();
  }

  void OutputWriter_HDF5::WriteTopology(UInt maxNumElemNodes,
                                        const std::vector<UInt>& elemTypes,
                                        const std::vector<UInt>& elems) {
    // =====================
    //  Element definitions
    // =====================
    UInt nElems = elemTypes.size();

    H5::Group elemGroup;
    try{
      elemGroup = meshGroup_.createGroup("Elements");
    } H5_CATCH( "Could not create element group" );

    // write connectivity
    H5IO::Write2DArray( elemGroup, "Connectivity", nElems,
                        maxNumElemNodes, &elems[0], dPropList_ );

    // write element types
    H5IO::Write1DArray( elemGroup, "Types", nElems,
                        &elemTypes[0], dPropList_ );

    std::cout << "done.\nWriting grid meta info... ";

    // ==========================
    //  Grid Meta Information
    // ==========================

    H5IO::WriteAttribute( elemGroup, "NumElems", nElems );
    std::vector< UInt > numElemsOfDim ( 3 );
    UInt numElemTypes = Elem::feType.map.size();
    std::map<Elem::FEType, UInt> numElemsOfType;
    std::vector<UInt>::const_iterator it, end;

    it = elemTypes.begin();
    end = elemTypes.end();
    UInt quadrElems = 0;

    for( ; it != end; it++ )
    {
      numElemsOfDim[ Elem::GetElemDim((Elem::FEType)*it)-1 ]++;
      quadrElems &= static_cast<UInt>(Elem::GetElemQuadratic((Elem::FEType)*it));
      numElemsOfType[(Elem::FEType)*it]++;
    }

    // This has still to to be checked
    H5IO::WriteAttribute( elemGroup, "QuadraticElems",
                          quadrElems );

    // number of elements per dimension
    for(UInt i=0; i<3; i++) {
      std::stringstream attrName;
      attrName << "Num" << (i+1) << "DElems";
      H5IO::WriteAttribute( elemGroup, attrName.str(), numElemsOfDim[i] );
    }

    // number of elements per type
    for(UInt i=0; i<numElemTypes; i++) {
      std::stringstream attrName;
      attrName << "Num_" << Elem::feType.ToString((Elem::FEType)i);
      H5IO::WriteAttribute( elemGroup, attrName.str(),
                            numElemsOfType[(Elem::FEType)i] );
    }

    std::cout << "done.\n";

    // close element group
    elemGroup.close();

  }

  void OutputWriter_HDF5::WriteRegion(std::string regionName,
                                      const std::vector<UInt>& regionElems,
                                      const std::vector<UInt>& regionNodes,
                                      UInt regionDim)
  {
    regionNames_.push_back(regionName);

    H5::Group regionListGroup;

    // Try to open Regions group. If it does not exist try to create it.
    try{
      regionListGroup = meshGroup_.openGroup("Regions");
    } catch(H5::GroupIException&) {;
      // create region group
      try{
        regionListGroup = meshGroup_.createGroup("Regions");
      } H5_CATCH( "Could not create region group" );
    }

    // create new region group
    H5::Group actRegionGroup;
    try {
      actRegionGroup = regionListGroup.createGroup(regionName);
    } H5_CATCH( "Could not create region group for region '"
                << regionName << "'" );
    H5IO::WriteAttribute( actRegionGroup, "Dimension",
                          regionDim );

    // create new node group
    H5IO::Write1DArray<UInt>( actRegionGroup, "Nodes",
                              regionNodes.size(),
                              (const UInt*)&regionNodes[0], dPropList_ );

    // create new element group
    H5IO::Write1DArray( actRegionGroup,
                        "Elements",
                        regionElems.size(),
                        (const Integer*)&regionElems[0],
                        dPropList_);

    try{
      actRegionGroup.close();
      regionListGroup.close();
    } H5_CATCH( "Could not close region group" );
  }

  void OutputWriter_HDF5::WriteNodeGroups(const std::map< std::string, std::vector<UInt> >&
                                          nodeGroups) {
    H5::Group myGroup;
    std::string nodesName;
    std::map< std::string, std::vector<UInt> >::const_iterator it, end;
    it = nodeGroups.begin();
    end = nodeGroups.end();

    for(; it != end; it++) {
      nodesName = it->first;

      std::cout << "Writing node group " << nodesName << "... ";

      // try to open group with given name
      try {
        myGroup = groupsGroup_.openGroup( nodesName );
      } catch (H5::Exception&) {
        try {
          myGroup = groupsGroup_.createGroup( nodesName );
        }
        H5_CATCH("Could not create node group '" << nodesName << "'");
      }
      H5IO::WriteAttribute( myGroup, "Dimension", (Integer)0 );
      H5IO::Write1DArray( myGroup, "Nodes",
                          it->second.size(), &it->second[0], dPropList_ );

      // close nodes array of current group
      myGroup.close();

      std::cout << "done." << std::endl;
    }
  }

  void OutputWriter_HDF5::WriteElemGroups(const std::map< std::string, std::vector<UInt> >&
                                          elemGroups,
                                          const std::map< std::string, std::vector<UInt> >&
                                          elemGroupNodes,
                                          const std::map< std::string, UInt >&
                                          groupDims)
  {
    Settings& settings = Settings::Instance();

    H5::Group myGroup;
    std::map <std::string, std::vector<UInt> >::const_iterator it, end, nodesIt;
    std::map <std::string, UInt >::const_iterator dimsIt;
    it = elemGroups.begin();
    end = elemGroups.end();
    nodesIt = elemGroupNodes.begin();
    dimsIt = groupDims.begin();

    for( ; it != end; it++, nodesIt++, dimsIt++ ) {
      const std::string& elemsName = it->first;
      std::cout << "Writing element group " << elemsName << "... ";

      try {
        myGroup = groupsGroup_.openGroup( elemsName );
      } catch (H5::Exception&) {
        try {
          myGroup = groupsGroup_.createGroup( elemsName );
        }
        H5_CATCH("Could not create element group '" << elemsName << "'");
      }
      H5IO::WriteAttribute( myGroup, "Dimension", dimsIt->second );
      H5IO::Write1DArray( myGroup, "Elements",
                          it->second.size(), &it->second[0],
                          dPropList_);

      // Write nodes of element group
      H5IO::Write1DArray( myGroup, "Nodes",
                          nodesIt->second.size(),
                          &nodesIt->second[0],
                          dPropList_);

      // close nodes array of current group
      myGroup.close();
      
      std::cout << "done." << std::endl;
    }

    groupsGroup_.close();

    if(settings.GetInt("justmesh"))
      doneWriting_ = true;
  }

  void OutputWriter_HDF5::WriteFlowData(CouplingHandler* cplHandler,
                                        UInt actRegion,
                                        const std::vector<std::string>& outputFields)
  {
    // Iterate over all flow datasets for this partition
    // and write them to HDF5
    UInt numDOFs;
    std::string groupName;
    H5::Group currResultGroup;

    FlowDataType::const_iterator fdIt, fdEnd;
    fdIt = (*flowData_)[actRegion].begin();
    fdEnd = (*flowData_)[actRegion].end();

    for( ; fdIt != fdEnd; fdIt++)
    {
      const FlowDataPartStruct& fdps = fdIt->second;
      // SolutionType st = fdIt->first;
      
      // Only write required datasets
      if(!fdps.isActive)
        continue;
      
      if(*outputFields.begin() != "all")
      {
        if(std::find(outputFields.begin(),
                     outputFields.end(),
                     fdps.resultName) == outputFields.end())
        {
          continue;
        }
      }
      
      // Create result groups in HDF5 file and write result.
      try {
        // Open or create subgroup for result
        try {
          currResultGroup = currMeshStepGroup_.openGroup( fdps.resultName );
        } catch (H5::GroupIException&)
        {
          try {
            currResultGroup = currMeshStepGroup_.createGroup( fdps.resultName );
          }
          H5_CATCH("Could not create results group '"
                   << fdps.resultName << "'");
        }
        
        // Create subgroup for region
        groupName = regionNames_[actRegion];
        currResultGroup = currResultGroup.createGroup(groupName);
        
        // Create subgroup for Nodes
        groupName = H5IO::MapUnknownTypeAsString(fdps.definedOn);
        currResultGroup = currResultGroup.createGroup(groupName);
        
        // Write result dataset
        numDOFs = fdps.dofNames.size();
        std::cout << "Writing result: " << fdps.resultName
                  << " on region " << regionNames_[actRegion] << "... ";
        
        std::vector<Double> output;
        cplHandler->ShrinkNodalVector(actRegion, numDOFs, fdps.data, output);
        
        WriteResults(currResultGroup, output, numDOFs, false);
        std::cout << "done." << std::endl;

        // Close nodes subgroup
        currResultGroup.close();

        if(stepNums_[fdps.resultName].empty() ||
           (*stepNums_[fdps.resultName].rbegin()) != actStepNum_)
        {
          stepNums_[fdps.resultName].push_back(actStepNum_);
          stepValues_[fdps.resultName].push_back(actStepValue_);
        }
      } H5_CATCH("Failed to write result '" << fdps.resultName
                 << "' in step " << actStepNum_);
    }
  }
  
  void OutputWriter_HDF5::WriteUserData(const std::map<std::string, std::string>& userData) {
    std::map<std::string, std::string>::const_iterator udIt, udEnd;

    udIt = userData.begin();
    udEnd = userData.end();
    for( ; udIt != udEnd; udIt++ )
      WriteStringToUserData(udIt->first, udIt->second);
  }

  void OutputWriter_HDF5::BeginResults(const std::vector<FlowDataType>* flowData)
  {
    flowData_ = flowData;
    
    // Initialize results tree in HDF5 file
    InitResultsGroup();    
  }
  
  void OutputWriter_HDF5::BeginStep(UInt stepNum, Double stepValue) {
    OutputWriter::BeginStep(stepNum, stepValue);
    Settings& settings = Settings::Instance();
    UInt externalFiles = settings.GetInt("extfiles");

    // Check, if step group is already open and create it if not
    if( currMeshStepGroup_.getId() <= 0 ) {
      std::stringstream stepName;
      stepName << "/Results/Mesh/MultiStep_1/Step_" << stepNum;
      
      // Create new step group.
      try {
        currMeshStepGroup_= mainGroup_.createGroup(stepName.str());
        H5IO::WriteAttribute( currMeshStepGroup_, "StepValue", actStepValue_ );
        
        if(externalFiles)
          CreateExternalFile(stepNum);
      } H5_CATCH( "Can not create dataset for step " << actStepNum_ );
    }
  }

  void OutputWriter_HDF5::EndStep() {
    OutputWriter::EndStep();

    // Close current step group in HDF5 file
    try {
      currMeshStepGroup_.close();
    } H5_CATCH( "Could close current step group" );
  }

  void OutputWriter_HDF5::EndResults()
  {
    // Finally write out result descriptions. This must be done in
    // the end since we do not know how many steps there are in advance.
    WriteResultDescriptions();
    
    doneWriting_ = true;
  }

  void OutputWriter_HDF5::WriteFileInfoHeader() {
    Settings& settings = Settings::Instance();
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
    creator << "cplreader-cfs " << CFS_VERSION << " ( "
            << CFS_SUBVERSION_REPOS
            << " rev. " << CFS_SUBVERSION_REV
            << ", date " << CFS_CONF_DATE
            << " )\n"
            << "Built on: " << CFS_BUILD_HOST << "\n"
            << "Platform: " << CFS_DISTRO << "_" << CFS_DISTRO_VER << "_" << CFS_ARCH << "\n"
            << "Built by: " << CFS_BUILD_USER << "\n"
            << "HDF5 library version: " << CFS_HDF5_VERSION;
    std::string creatorString = creator.str();
    H5IO::Write1DArray( infoGroup, "Creator", 1, &creatorString, dPropList_ );

    // create group for content
    StdVector<Integer> content;
    content.Push_back( H5IO::MapCapabilityType( SimOutput::USERDATA ) );
    if(!settings.GetInt("justinit"))
      content.Push_back( H5IO::MapCapabilityType( SimOutput::MESH ) );
    if(!settings.GetInt("justmesh") && !settings.GetInt("justinit"))
      content.Push_back( H5IO::MapCapabilityType( SimOutput::MESH_RESULTS ) );
    H5IO::Write1DArray( infoGroup, "Content", content.GetSize(),
                        &content[0], dPropList_ );

  }

  void OutputWriter_HDF5::InitHDF5() {
    Settings& settings = Settings::Instance();
    std::string pathsep;
    std::string fileName  = settings.GetString("name");
    std::ostringstream strBuffer;

    // Set compression and chunk size parameters for HDF5
    UInt compressionLevel = 9;
    UInt maxChunkSize = 4096;
    dPropList_ = H5::DSetCreatPropList::DEFAULT;
    dPropList_.setLayout( H5D_CHUNKED );
    dPropList_.setDeflate( compressionLevel );
    H5IO::SetMaxChunkSize( maxChunkSize );

    // Do not print HDF5 exceptions by default
    H5::Exception::dontPrint();

    // concatenate output file name
    try {
      fs::create_directory( hdf5DirName_ );
      pathsep = fs::path("/").native_directory_string();
    } catch (std::exception &ex) {
      EXCEPTION(ex.what());
    }

    strBuffer.clear();
    strBuffer.str("");
    strBuffer << hdf5DirName_ << pathsep << fileName << ".h5";
    mainFileName_ = strBuffer.str();

    // create main file and obtain main group
    try {
      mainFile_ = H5::H5File (mainFileName_, H5F_ACC_TRUNC );
    } H5_CATCH( "Could not create hdf5 file '" << fileName << "' : " );

    mainGroup_ = mainFile_.openGroup( "/" );
    meshGroup_ = mainGroup_.createGroup( "Mesh" );
    groupsGroup_ = meshGroup_.createGroup( "Groups" );
    mainGroup_.createGroup( "FileInfo" ).close();

    mainGroup_.createGroup( "UserData" );
    mainGroup_.createGroup( "Results" );
  }

  void OutputWriter_HDF5::WriteResultDescriptions() {
    std::string resultName, unit;
    UInt definedOn, numDofs, entryType;
    std::vector<std::string> dofNames;
    FlowDataType::const_iterator it, end, it_tmp, end_tmp;
    H5::Group resultDescGroup;
    H5::Group msGroup;
    std::vector<std::string> resultRegions;

    try {
      // open the group for the result description datasets.
      msGroup = mainGroup_.openGroup("/Results/Mesh/MultiStep_1");
      resultDescGroup = msGroup.openGroup("ResultDescription");

      H5IO::WriteAttribute( msGroup, "LastStepNum", actStepNum_ );
      H5IO::WriteAttribute( msGroup, "LastStepValue", actStepValue_ );
      msGroup.close();

    } H5_CATCH( "Could not open result description group" );

    // Extract active regions
    UInt numRegions = (*flowData_).size();

    // Write result descriptions
    it = (*flowData_)[0].begin();
    end = (*flowData_)[0].end();

    for( ; it != end; it++ ) {
      if(!it->second.isActive)
        continue;

      resultName = it->second.resultName;

      resultRegions.clear();
      for(UInt iRegion = 0; iRegion < numRegions; ++iRegion)
      {
        it_tmp = (*flowData_)[iRegion].begin();
        end_tmp = (*flowData_)[iRegion].end();

        for( ; it_tmp != end_tmp; it_tmp++ ) {
          if (resultName == it_tmp->second.resultName)
          {
            resultRegions.push_back(regionNames_[iRegion]);
          }
        }
      }
      
      // if stepNums_[resultName] is empty,
      // then this result is no output field
      if (stepNums_[resultName].empty())
        continue;
      
      definedOn = H5IO::MapUnknownType(it->second.definedOn);
      dofNames = it->second.dofNames;
      unit = it->second.unit;
      numDofs = dofNames.size();
      entryType = H5IO::MapEntryType(it->second.entryType);

      try {
        // === Second version: Separate datasets for each entry
        H5::Group actGroup = resultDescGroup.createGroup(resultName);

        H5IO::Write1DArray( actGroup, "DefinedOn", 1, &definedOn, dPropList_ );
        H5IO::Write1DArray( actGroup, "EntityNames", resultRegions.size(),
            &resultRegions[0], dPropList_ );
        H5IO::Write1DArray( actGroup, "NumDOFs", 1, &numDofs, dPropList_ );
        H5IO::Write1DArray( actGroup, "DOFNames", dofNames.size(),
            &dofNames[0], dPropList_ );
        H5IO::Write1DArray( actGroup, "EntryType", 1, &entryType, dPropList_ );
        H5IO::Write1DArray( actGroup, "Unit", 1, &unit, dPropList_ );

        H5IO::Write1DArray<Double>( actGroup, "StepValues",
            stepValues_[resultName].size(), &stepValues_[resultName][0], dPropList_ );
        H5IO::Write1DArray<UInt>( actGroup, "StepNumbers",
            stepNums_[resultName].size(), &stepNums_[resultName][0], dPropList_ );

        actGroup.close();

      } H5_CATCH( "Could not write result description for result '"
          << resultName << "'" );
    }

    try {
      resultDescGroup.close();
    } H5_CATCH( "Could not close result description group" );
  }

  void OutputWriter_HDF5::WriteResults( H5::Group& resultGroup,
                                        const std::vector<Double>& resultVals,
                                        const UInt numDOFs,
                                        const bool isImag ) {
    Settings& settings = Settings::Instance();
    std::string outputPrecision = settings.GetString("outprec");

    // create dataset with related name
    std::string name;
    if( !isImag )
      name = "Real";
    else
      name = "Imag";

    UInt numEntities = (UInt) resultVals.size() / numDOFs;

    if(outputPrecision == "double")
    {
      H5IO::Write2DArray( resultGroup, name,
                          numEntities, numDOFs, &resultVals[0],
                          dPropList_ );
    }
    else
    {
      std::vector<Float> floatResultVals(resultVals.size());

      std::vector<Double>::const_iterator it = resultVals.begin();
      std::vector<Double>::const_iterator end = resultVals.end();
      std::vector<Float>::iterator flIt = floatResultVals.begin();

      for( ; it != end; it++ )
      {
        *flIt = static_cast<Float>(*it);
        flIt++;
      }

      H5IO::Write2DArray( resultGroup, name,
                          numEntities, numDOFs, &floatResultVals[0],
                          dPropList_ );
    }
  }

  void OutputWriter_HDF5::CreateExternalFile(UInt timeStep) {
    Settings& settings = Settings::Instance();
    std::string fileName = settings.GetString("name");
    std::stringstream fName, masterGroup;
    std::string pathsep, fn;

    try {

      // open external file
      pathsep = fs::path("/").native_directory_string();

      fName << fileName << "_ms" << 1 << "_step"
            << timeStep << ".h5";
      fn = fName.str();
      fName.str("");
      fName << hdf5DirName_ << pathsep << fn;
      currStepFileName_ = fName.str();
      
      currStepFile_ = H5::H5File(currStepFileName_, H5F_ACC_TRUNC);

      // Write reference to external file to main file
      H5IO::WriteAttribute( currMeshStepGroup_, "ExtHDF5FileName", fn );

      // set current step group to external file
      currMeshStepGroup_.close();
      currMeshStepGroup_ = currStepFile_.openGroup("/");

      // Store reference to master file in external file
      fName.str("");
      fName << fileName << ".h5";
      fn = fName.str();
      H5IO::WriteAttribute( currMeshStepGroup_, "MasterHDF5FileName", fn );

      // Store reference to main group in external file
      masterGroup << "/Results/Mesh/MultiStep_" << 1
                  << "/Step_" << timeStep;
      H5IO::WriteAttribute( currMeshStepGroup_, "MasterGroup", masterGroup.str() );
    } H5_CATCH( "Could not create external file" );
  }

  void OutputWriter_HDF5::WriteStringToUserData(const std::string& dSetName,
                                                const std::string& str) {
    H5::Group userDataGroup;

    // If it does not exist, create Group for Data.
    try {
      userDataGroup = mainGroup_.openGroup("UserData");
    } H5_CATCH( "Can not write meta information to hdf5 file" );

    H5IO::Write1DArray( userDataGroup, dSetName, 1, &str, dPropList_ );
    userDataGroup.close();

  }

  void OutputWriter_HDF5::InitResultsGroup()
  {
    Settings& settings = Settings::Instance();
    UInt externalFiles = settings.GetInt("extfiles");
    std::string analysisType = "transient";

    try {
      resultsGroup_ = mainGroup_.openGroup("Results");
    } H5_CATCH( "Could open group for results" );

    try {
      meshResultsGroup_ = resultsGroup_.createGroup("Mesh");
      H5IO::WriteAttribute( meshResultsGroup_, "ExternalFiles", externalFiles );

      H5::Group msGroup = meshResultsGroup_.createGroup("MultiStep_1");

      H5IO::WriteAttribute( msGroup, "AnalysisType", analysisType );

      msGroup.createGroup("ResultDescription").close();
      msGroup.close();
      meshResultsGroup_.close();
    } H5_CATCH( "Could init group for results" );
  }
}
