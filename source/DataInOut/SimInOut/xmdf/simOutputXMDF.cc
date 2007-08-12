// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

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

#include "DataInOut/CommandLine/BaseCommandLineHandler.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/exception.hh"
#include "simOutputXMDF.hh"
#include "hdf5io.hh"

namespace fs = boost::filesystem;

namespace CoupledField {
  
#define H5_EXCEPTION(STR, EX)                                           \
  EXCEPTION( STR, EX.getCDetailMsg() );

#define H5_CATCH( STR )                                                 \
  catch (H5::Exception& h5Ex ) {                                        \
    EXCEPTION( STR << ":\n" << h5Ex.getCDetailMsg() );                  \
  }

  
  SimOutputXMDF::SimOutputXMDF(std::string fileName, ParamNode * inputNode) :
    SimOutput(fileName, inputNode) 
  {
    ENTER_FCN( "SimOutputXMDF::XMDF" );
    fileName_ = fileName;
    formatName_ = "hdf5";
    dirName_ = "simoutput_hdf5";
    
    capabilities_.insert( MESH );
    capabilities_.insert( MESH_RESULTS );
    
    gridWritten_ = false;
    externalFiles_ = false;
    printGridOnly_ = false;

    currMS_ = 0;
    currStep_ = 0;
    
    // Change defaults according to XML file
    myParam_->Get("externalFiles", externalFiles_, false);
    
    // Do not print HDF5 exceptions by default
    H5::Exception::dontPrint();
  }


  SimOutputXMDF::~SimOutputXMDF() 
  {
    ENTER_FCN( "SimOutputXMDF::~XMDF" );
    
    // close groups
    dataGroup_.close();
    volDataGroup_.close();
    mainGroup_.close();

    // check for open groups, datasets etc.
    if (mainFile_.getObjCount( H5F_OBJ_DATASET | 
                               H5F_OBJ_GROUP | 
                               H5F_OBJ_DATATYPE | H5F_OBJ_ATTR) > 0 ) {
      std::cerr << "There are still objects open in the hdf5 file\n\n";
      CheckOpenObjects();
    }

    mainFile_.close();
  }


  void SimOutputXMDF::CheckOpenObjects() {
    // check for open groups, datasets etc.
    std::cerr << "Number of open objects:\n"
              << "--------------------------";
    std::cerr << "Datasets: "<<  mainFile_.getObjCount( H5F_OBJ_DATASET) << std::endl;
    std::cerr << "Groups: "<<  mainFile_.getObjCount( H5F_OBJ_GROUP ) << std::endl;
    std::cerr << "DataTypes: "<<  mainFile_.getObjCount( H5F_OBJ_DATATYPE) << std::endl;
    std::cerr << "Attributes: "<<  mainFile_.getObjCount( H5F_OBJ_ATTR) << std::endl;


  }
  
  void SimOutputXMDF::Init( Grid* ptGrid, bool printGridOnly ) 
  {
    ptGrid_ = ptGrid;
    printGridOnly_ = printGridOnly;
    WriteGrid();
  }

  void SimOutputXMDF::BeginMultiSequenceStep( UInt step,
                                              AnalysisType type,
                                              UInt numSteps  ) 
  {
    std::stringstream msName;
    std::string analysisType;
    H5::Group resultDescGroup;    
    // If it does not exist, create Group for Grid / Volume data
    try {
      dataGroup_ = mainGroup_.openGroup("Results");
      volDataGroup_ = dataGroup_.openGroup("Grid");
    } catch (H5::Exception& h5ex) {
      try {
        volDataGroup_ = dataGroup_.createGroup("Grid");
      } H5_CATCH( "Could not create group for grid results" );
    }
    
    // Assemble name of multistep
    msName << "MultiStep_" << step;
    
    try {
      // create new multistep group.
      currMSGroup_ = volDataGroup_.createGroup(msName.str());

      // add the attribute AnalysisType to multisequence step.
      Enum2String(type, analysisType);
      currAnalysisType_ = type;
      H5IO::WriteAttribute( currMSGroup_, "AnalysisType", analysisType );
      
      // add a group for the result description datasets.
      resultDescGroup= currMSGroup_.createGroup("ResultDescription");
      
      // write result meta information to file
      WriteResultDescriptions( resultDescGroup );
      resultDescGroup.close();

    } H5_CATCH( "Could not create group for multi sequence step " << step );
    currMS_ = step;
  }

  void SimOutputXMDF::RegisterResult( shared_ptr<BaseResult> sol,
                                      UInt saveBegin, UInt saveInc,
                                      UInt saveEnd ) {
    registeredResults_[sol->GetResultInfo()->resultName].push_back(sol);
  }
  
  void SimOutputXMDF::BeginStep( UInt stepNum, Double stepVal ) 
  {
    std::stringstream stepName;
    
    // Assemble name of multistep
    stepName << "Step_" << stepNum;
    currStep_ = stepNum;

    try {
      // Create new step group.
      currStepGroup_ = currMSGroup_.createGroup(stepName.str());
      H5IO::WriteAttribute( currStepGroup_, "StepVal", stepVal );
      
      if(externalFiles_)
        CreateExternalFile();
    } H5_CATCH( "Can not create dataset for step " << currStep_ );
    
  }

  void SimOutputXMDF::AddResult( shared_ptr<BaseResult> sol ) 
  {
    H5::Group resultGroup, subGroup, regionGroup;
    UInt numDOFs;
    Vector<Double> realVec, imagVec;
    std::vector<std::string> resultNames;
    
    std::string regionName = sol->GetEntityList()->GetName();
    shared_ptr<ResultInfo> resInfo = sol->GetResultInfo();
    std::string resultName = resInfo->resultName;
    numDOFs = resInfo->dofNames.GetSize();

    // check if group for result exists already
    try {
      resultGroup = currStepGroup_.openGroup( resultName );
    } catch( H5::Exception& h5Ex ) {
      resultGroup = currStepGroup_.createGroup( resultName );
    }

    // determine, on whicht type of entity the result is defined
    std::string entityString;
    switch( resInfo->definedOn ) {
    case ResultInfo::NODE:
      entityString = "Nodes";
      break;
    case ResultInfo::ELEMENT:
      entityString = "Elements";
      break;
    default:
      EXCEPTION( "Only nodal and element results can be stored to hdf5 "
                 << "files up to now ");
    }

    // try to create regionGroup
    try {
      regionGroup = resultGroup.openGroup( regionName );
    } catch( H5::Exception& h5Ex ) {
      regionGroup = resultGroup.createGroup( regionName );
    }

    // try to create subgroup for entity
    try {
      subGroup = regionGroup.createGroup( entityString );
    } H5_CATCH( "Could not create subgroup " << entityString
                << " for result " << resultName << " on region "
                << regionName );
    
    if( sol->GetEntryType() == EntryType::DOUBLE ) {

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

  void SimOutputXMDF::FinishStep( ) 
  {
    if(externalFiles_)
      currStepFile_.close();

    currStepGroup_.close();
  }

  void SimOutputXMDF::FinishMultiSequenceStep( ) 
  {
    registeredResults_.clear();
    currMSGroup_.close();
  }

  void SimOutputXMDF::Finalize() 
  {

    // return, if only the grid is to be printed
    if( printGridOnly_ )
      return;

    std::vector<std::string> fileNames;
    std::vector<std::string> dataSetNames;
    std::ifstream fin;
    std::string dumpStr;

    fileNames.push_back(commandLine->GetParamFile());
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

    commandLine->GetDumpString(dumpStr);
    WriteStringToUserData("ProgramStats", dumpStr);
  }

  void SimOutputXMDF::InitModule()
  {
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
    
    std::string newNames [3];
    H5IO::ReadArray( mainGroup_, "Names", newNames );
    for( UInt i = 0; i < 3; i++ ) {
      std::cerr << newNames[i] << std::endl;
    }

  }

  void SimOutputXMDF::WriteGrid() {
  
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
    FEType eType;
    RegionIdType region;

    // iterate over all elements
    for( UInt i = 0; i < nElems; i++ ) {
      elConnect.resize( maxNumNodes );
      std::fill(elConnect.begin(), elConnect.end(),
                0 );
      ptGrid_->GetElemData( i+1, eType, region, &elConnect[0] );
      numElemsOfDim[ELEM_DIM[eType]-1]++;
      feTypes[i] = eType;

      // insert connectivity into global array
      offset = i * maxNumNodes;
      for( UInt j = 0; j < elConnect.size(); j++ ) {
        connect[offset + j] = elConnect[j];
      }
    }

    // write connectivity
    H5IO::Write2DArray( elemGroup, "Connectivity", nElems, 
                        maxNumNodes, &connect[0] );

    // write element types
    H5IO::Write1DArray( elemGroup, "Types", nElems,
                        &feTypes[0] );


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
    UInt numElemTypes = sizeof(ELEM_DIM) / sizeof(UInt);
    for(UInt i=0; i<numElemTypes; i++) {
      std::stringstream attrName;
      attrName << "Num_" << ELEM_TYPE_NAMES[i].substr(3);
      H5IO::WriteAttribute( elemGroup, attrName.str(),
                            ptGrid_->GetNumElemOfType((FEType)i) );
    }

    // close element group
    elemGroup.close();
    
    // ============================================
    //  Write Regions, Named Elements, Named Nodes
    // ============================================
    WriteRegions( meshGroup_ );
    WriteNamedElems( meshGroup_ );
    gridWritten_ = true;

    // close meshGroup
    meshGroup_.close();
  }


  void SimOutputXMDF::WriteRegions(const H5::Group& meshGroup)
  {
    H5::Group regionListGroup;
    std::vector< std::string > regionNames;
    std::vector< UInt > regionDims;
    std::vector< std::vector<UInt> > regionElems;
    std::vector< StdVector<UInt> > regionNodes;
    StdVector<Elem*> elems;
    std::vector<RegionIdType> surfRegionIds, volRegionIds;
    UInt dim, idx, numRegions;

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
    regionDims.resize(numRegions);
    regionElems.resize(numRegions);
    regionNodes.resize(numRegions);
    regionDims.resize(numRegions);

    // obtain nodes and elements per surface region
    for(UInt i=0; i<surfRegionIds.size(); i++) {
      idx = surfRegionIds[i];
      regionDims[idx] = dim-1;
      
      ptGrid_->GetElems(elems, idx);
      regionElems[idx].resize(elems.GetSize());
      for(UInt j=0; j<elems.GetSize(); j++) {
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
    for(UInt i=0; i<volRegionIds.size(); i++) {
      idx = volRegionIds[i];
      regionDims[idx] = dim;

      ptGrid_->GetElems(elems, idx);
      regionElems[idx].resize(elems.GetSize());
      for(UInt j=0; j<elems.GetSize(); j++) {
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
                          (const UInt*)&regionNodes[i][0] );
      
      // create new element group
      H5IO::Write1DArray( actRegionGroup,
                          "Elements",
                          regionElems[i].size(),
                          (const Integer*)&regionElems[i][0] );

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

  void SimOutputXMDF::WriteNamedNodes(const H5::Group& meshGroup)
  {
    H5::Group namedNodeGroup;
    std::vector< UInt > nodes;
    std::vector<std::string> nodeNames;
    UInt numNamedNodes = 0;

    // create group for named nodes
    try {
      namedNodeGroup = meshGroup.createGroup( "NamedNodes" );
    } H5_CATCH( "Could not create group for named nodes" );
    
    // obtain list with names of nodes
    ptGrid_->GetListNodeNames(nodeNames);
    numNamedNodes = nodeNames.size();
    
    for(UInt i = 0; i < numNamedNodes; i++ ) {
      ptGrid_->GetNodesByName(nodes, nodeNames[i]);      
      H5IO::Write1DArray( namedNodeGroup, nodeNames[i],
                          nodes.size(), &nodes[0] );
    }

    // close named node group
    namedNodeGroup.close();
  }

  void SimOutputXMDF::WriteNamedElems(const H5::Group& meshGroup)
  {
    H5::Group namedElemGroup;
    std::vector< UInt > elemNums;
    StdVector<Elem*> elems;
    std::vector<std::string> elemNames;

    // create group for named elements
    try {
      namedElemGroup = meshGroup.createGroup( "NamedElems" );
    } H5_CATCH( "Could not create group for named elems" );
    
    // obtain list with names of elements
    ptGrid_->GetListElemNames(elemNames);
    UInt numNamedElems = elemNames.size();

    for(UInt i = 0; i < numNamedElems; i++ ) {
      ptGrid_->GetElemsByName(elems, elemNames[i]);
      elemNums.resize( elems.GetSize() );
      for( UInt j = 0; j < elems.GetSize(); j++ ) {
        elemNums[j] = elems[j]->elemNum;
      }
      
      H5IO::Write1DArray( namedElemGroup, elemNames[i],
                          elemNums.size(), &elemNums[0] );
    }

    // write namedElemsGroup
    namedElemGroup.close();
  }

  void SimOutputXMDF::WriteResultDescriptions( const H5::Group& resGroup )
  {
    std::string resultName;
    UInt definedOn;
    UInt numDOFs;
    UInt entryType;
    std::string unit;
    std::vector<std::string> regions;
    ResDescType::const_iterator it, end;
    std::vector< boost::shared_ptr<BaseResult> >::const_iterator solIt, solEnd;
    boost::shared_ptr<ResultInfo> resInfo, actResInfo;

    it = registeredResults_.begin();
    end = registeredResults_.end();

    for( ; it != end; it++ ) {
      resultName = it->first;
      
      solIt = it->second.begin();
      solEnd = it->second.end();
      
      resInfo = (*solIt)->GetResultInfo();
      numDOFs = resInfo->dofNames.GetSize();
      unit = resInfo->unit;
      definedOn = H5IO::MapUnknownType( resInfo->definedOn );
      entryType = H5IO::MapEntryType( resInfo->entryType );
      
      // Generate list of regions for the current result.
      regions.clear();
      for( ; solIt != solEnd; solIt++ ) {
        actResInfo = (*solIt)->GetResultInfo();
        
        regions.
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
        
        H5IO::Write1DArray( actGroup, "DefinedOn", 1, &definedOn );
        H5IO::Write1DArray( actGroup, "Regions", regions.size(), &regions[0] );
        H5IO::Write1DArray( actGroup, "NumDOFs", 1, &numDOFs );
        H5IO::Write1DArray( actGroup, "DOFNames", resInfo->dofNames.GetSize(), 
                            &(resInfo->dofNames[0]) );
        H5IO::Write1DArray( actGroup, "EntryType", 1, &entryType );
        H5IO::Write1DArray( actGroup, "Unit", 1, &unit );

        actGroup.close();


      } H5_CATCH( "Could not write result description for result '"
                  << resultName << "'" );
    }
  }

  void SimOutputXMDF::WriteResults( H5::Group& resultGroup,
                                    Vector<Double>& resultVals,
                                    const UInt numDOFs,
                                    const bool isImag ) 
  {

    // create dataset with related name
    std::string name;
    if( !isImag )
      name = "Real";
    else
      name = "Imag";

    UInt numEntities = (UInt) resultVals.GetSize() / numDOFs;

    H5IO::Write2DArray( resultGroup, name, 
                        numEntities, numDOFs, &resultVals[0] );
  }

  void SimOutputXMDF::CreateExternalFile()
  {
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
      H5IO::WriteAttribute( currStepGroup_, "ExtHDF5FileName", fn );

      // set current step group to external file
      currStepGroup_.close();
      currStepGroup_ = currStepFile_.openGroup("/");
      
      // Store reference to master file in external file
      fName.str("");
      fName << fileName_ << ".h5";
      fn = fName.str();
      H5IO::WriteAttribute( currStepGroup_, "MasterHDF5FileName", fn );
      
      // Store reference to main group in external file
      masterGroup << "/Results/Grid/MultiStep_" << currMS_
                  << "/Step_" << currStep_;
      H5IO::WriteAttribute( currStepGroup_, "MasterGroup", masterGroup.str() );
    } H5_CATCH( "Could not create external file" );
  }
  
  void SimOutputXMDF::WriteStringToUserData(const std::string& dSetName, 
                                            const std::string& str) {
    H5::Group userDataGroup;
    
    // If it does not exist, create Group for Data.
    try {
      userDataGroup = mainGroup_.openGroup("UserData");
    } H5_CATCH( "Can not write meta information to hdf5 file" );
    
    H5IO::Write1DArray( userDataGroup, dSetName, 1, &str );
    userDataGroup.close();
  }

} // end of namespace

