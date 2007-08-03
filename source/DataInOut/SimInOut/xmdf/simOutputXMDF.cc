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
  
#define CHECK_HDF5_ERROR(HID_T, STR) { if(HID_T < 0) {                  \
      std::ostringstream ostr;                                          \
      ostr << STR;                                                      \
      Exception ex(NULL, __FILE__, __LINE__, ostr.str().c_str()); } }

#define H5_EXCEPTION(STR, EX)                                           \
  EXCEPTION( STR, EX.getCDetailMsg() );

#define H5_CATCH( STR )                                                 \
  catch (H5::Exception& h5Ex ) {                                        \
    EXCEPTION( STR << h5Ex.getCDetailMsg() );                           \
  }

  
  SimOutputXMDF::SimOutputXMDF(std::string fileName, ParamNode * inputNode) :
    SimOutput(fileName, inputNode) 
  {
    ENTER_FCN( "SimOutputXMDF::XMDF" );
    fileName_ = fileName;
    fileId_ = -1;
    formatName_ = "hdf5";
    dirName_ = "simoutput_hdf5";

    capabilities_.insert( MESH );
    capabilities_.insert( MESH_RESULTS );

    gridWritten_ = false;
    externalFiles_ = false;
    printGridOnly_ = false;
    
    // Change defaults according to XML file
    myParam_->Get("externalFiles", externalFiles_, false);

    // Do not print HDF5 exceptions by default
    H5::Exception::dontPrint();
  }


  SimOutputXMDF::~SimOutputXMDF() {
    ENTER_FCN( "SimOutputXMDF::~XMDF" );

    mainFile_.close();
  }

  
  void SimOutputXMDF::Init( Grid* ptGrid, bool printGridOnly ) {
    ptGrid_ = ptGrid;
    printGridOnly_ = printGridOnly;
    WriteGrid();
  }

  void SimOutputXMDF::BeginMultiSequenceStep( UInt step,
                                              AnalysisType type,
                                              UInt numSteps  ) 
  {
    std::stringstream msName;
    H5::Attribute attrib;
    H5::DataSpace space;
    std::string analysisType;

    // If it does not exist, create Group for Data.
    try 
    {
      dataGroup_ = mainGroup_.openGroup("Data");
      volDataGroup_ = dataGroup_.openGroup("Volume");
    }
    catch (H5::Exception& h5ex)
    {
      try 
      {
        dataGroup_ = mainGroup_.createGroup("Data");
        volDataGroup_ = dataGroup_.createGroup("Volume");
      }
      catch (H5::Exception& h5ex)
      {
        EXCEPTION(h5ex.getCDetailMsg());
      }
    }

    // Assemble name of multistep
    msName << "Multistep " << step;

    try
    {
      // Create new multistep group.
      currMSGroup_ = volDataGroup_.createGroup(msName.str());

      // Add the attribute AnalysisType to multisequence step.
      Enum2String(type, analysisType);
      currAnalysisType_ = type;
      H5IO::WriteAttribute( currMSGroup_, "AnalysisType", analysisType );

      // Add a group for the attribute description datasets.
      currAttrDescGroup_ = currMSGroup_.createGroup("ResultDescription");

      std::cerr << "before wriring AttributeDescriptions\n";
      WriteAttribDescriptions();

    }
    catch (H5::Exception& h5ex)
    {
      EXCEPTION(h5ex.getCDetailMsg());
    }

    currMS_ = step;
  }

  void SimOutputXMDF::RegisterResult( shared_ptr<BaseResult> sol,
                                      UInt saveBegin, UInt saveInc,
                                      UInt saveEnd ) 
  {
    registeredResults_[sol->GetResultInfo()->resultName].push_back(sol);
  }

  void SimOutputXMDF::BeginStep( UInt stepNum, Double stepVal ) 
  {
    std::stringstream stepName;
    H5::Attribute attrib;
    H5::DataSpace space;
    Double timeValue;
    Double freqValue;

    switch(currAnalysisType_)
    {
    case STATIC:

    case TRANSIENT:
      timeValue = stepVal;
      freqValue = 0;
      break;

    case HARMONIC:

    case EIGENFREQUENCY:
      timeValue = 0;
      freqValue = stepVal;
      break;

    case TRANSIENTHARMONIC:
      EXCEPTION("Don't know how to treat TRANSIENTHARMONIC in SimOutputXMDF.");
      break;

    case MULTI_SEQUENCE:
      EXCEPTION("Don't know how to treat MULTI_SEQUENCE in SimOutputXMDF.");
      break;
    }


    // Assemble name of multistep
    stepName << "Step " << stepNum;
    currStep_ = stepNum;

    try
    {
      // Create new step group.
      currStepGroup_ = currMSGroup_.createGroup(stepName.str());
      H5IO::WriteAttribute( currStepGroup_, "StepFrequency",
                            freqValue );

      H5IO::WriteAttribute( currStepGroup_, "StepTime",
                            timeValue );

      if(externalFiles_)
        CreateExternalFile();
    }
    catch (H5::Exception& h5ex)
    {
      EXCEPTION(h5ex.getCDetailMsg());
    }
  }

  void SimOutputXMDF::AddResult( shared_ptr<BaseResult> sol ) 
  {
    std::string regionName = sol->GetEntityList()->GetName();
    std::string resultName = (*sol->GetResultInfo()).resultName;
    std::stringstream regionStr;
    UInt regionNum;
    UInt numDOFs;
    std::vector< Vector<Double>* > results;
    Vector<Double> complexVec1;
    Vector<Double> complexVec2;
    std::vector<std::string> resultNames;


    if( sol->GetEntryType() == EntryType::DOUBLE ) {

      Vector<Double> & resultVec = dynamic_cast<Result<Double>&>
                                   (*sol).GetVector();

      numDOFs = (*sol->GetResultInfo()).dofNames.GetSize();

      results.push_back(&resultVec);
      resultNames.push_back(resultName);
    } else {
      Vector<Complex> & resultVec = dynamic_cast<Result<Complex>&>
                                    (*sol).GetVector();

      numDOFs = (*sol->GetResultInfo()).dofNames.GetSize();
      UInt numEntities = resultVec.GetSize() / numDOFs;

      complexVec1.Resize(resultVec.GetSize());
      complexVec2.Resize(resultVec.GetSize());

      for(UInt i=0; i<numEntities; i++)
      {
        complexVec1[i] = resultVec[i].real();
        complexVec2[i] = resultVec[i].imag();        
      }

      results.push_back(&complexVec1);
      results.push_back(&complexVec2);

      std::stringstream resName;
      if((*sol->GetResultInfo()).complexFormat == REAL_IMAG)
      {
        resName << resultName << "_real";
        resultNames.push_back(resName.str());
        resName.str("");
        resName << resultName << "_imag";
        resultNames.push_back(resName.str());
      }
      else 
      {
        resName << resultName << "_ampl";
        resultNames.push_back(resName.str());
        resName.str("");
        resName << resultName << "_phase";
        resultNames.push_back(resName.str());
      }
    }


    regionNum = ptGrid_->RegionNameToId(regionName)+1;
    regionStr << "Region " << regionNum;

    WriteAttributes(resultNames, results, regionStr.str(), numDOFs);
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
    currAttrDescGroup_.close();
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

    dataSetNames.push_back("Parameter File");
    dataSetNames.push_back("Material File");
    
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
      //WriteStringToUserData(dataSetNames[i], str);
      fin.close();
    }

    commandLine->GetDumpString(dumpStr);
    //WriteStringToUserData("Program Stats", dumpStr);
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
    for (UInt i = 0; i < nNodes; i+= 3) {
      Point p;
      ptGrid_->GetNodeCoordinate(p, i+1);
      locs[i+0] = p[0];
      locs[i+1] = p[1];
      locs[i+2] = p[2];
    }

    // write nodal coordinates to file
    H5IO::Write2DArray( nodeGroup, "Coordinates", 
                        nNodes, 3, &locs[0] );
    
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
    std::vector< UInt > numElemsOfDim (ptGrid_->GetDim());

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
      attrName << "NUM_" << ELEM_TYPE_NAMES[i].substr(3);
      H5IO::WriteAttribute( elemGroup, attrName.str(),
                            ptGrid_->GetNumElemOfType((FEType)i) );
    }
    
    // ============================================
    //  Write Regions, Named Elements, Named Nodes
    // ============================================
    WriteRegions( meshGroup_ );
    WriteNamedNodes( meshGroup_ );
    WriteNamedElems( meshGroup_ );
    
    gridWritten_ = true;
  }


  void SimOutputXMDF::WriteRegions(const H5::Group& meshGroup)
  {
    H5::Group regionListGroup, actRegionGroup;
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
    }


    // loop over regions and write out nodes and elements
    for(UInt i=0; i<numRegions; i++)
    {
      // create new region group 
      try {
        actRegionGroup = regionListGroup.createGroup(regionNames[i] );
      } H5_CATCH( "Could not create region group for region '" 
                  << regionNames[i] << "'" );
      H5IO::WriteAttribute( actRegionGroup, "Dimension", 
                            ptGrid_->GetDim() );
      
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

    }
      

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
  }

  void SimOutputXMDF::WriteAttribDescriptions()
  {
    std::string resultName;
    UInt definedOn;
    UInt numDOFs;
    UInt entryType;
    hsize_t numRegions;
    std::string unit;
    std::vector<std::string> dofNames;
    std::vector<UInt> regions;
    RegionIdType regionId;
    ResDescType::const_iterator it, end;
    std::vector< boost::shared_ptr<BaseResult> >::const_iterator solIt, solEnd;
    boost::shared_ptr<ResultInfo> resInfo, actResInfo;

    it = registeredResults_.begin();
    end = registeredResults_.end();

    for( ; it != end; it++ )
    {
      resultName = it->first;

      solIt = it->second.begin();
      solEnd = it->second.end();

      resInfo = (*solIt)->GetResultInfo();
      numDOFs = resInfo->dofNames.GetSize();
      unit = resInfo->unit;

      switch(resInfo->definedOn)
      {
        case ResultInfo::NODE:
          definedOn = 1;
          break;
        case ResultInfo::EDGE:
          definedOn = 2;
          break;
        case ResultInfo::FACE:
          definedOn = 3;
          break;
        case ResultInfo::ELEMENT:
          definedOn = 4;
          break;
        case ResultInfo::SURF_ELEM:
          definedOn = 5;
          break;
        case ResultInfo::PFEM:
          definedOn = 6;
          break;
        case ResultInfo::REGION:
          definedOn = 7;
          break;
        case ResultInfo::SURF_REGION:
          definedOn = 8;
          break;
        case ResultInfo::NODELIST:
          definedOn = 9;
          break;
        case ResultInfo::COIL:
          definedOn = 10;
          break;
        case ResultInfo::FREE:
          definedOn = 11;
          break;
      }

      switch(resInfo->entryType)
      {
        case ResultInfo::UNKNOWN:
          entryType = 0;
          break;
        case ResultInfo::SCALAR:
          entryType = 1;
          break;
        case ResultInfo::VECTOR:
          entryType = 3;
          break;
        case ResultInfo::TENSOR:
          entryType = 6;
          break;
        case ResultInfo::STRING:
          entryType = 32;
          break;
      }

      dofNames.resize( resInfo->dofNames.GetSize() );
      for(UInt i = 0; i<numDOFs; i++) {
        dofNames[i] = resInfo->dofNames[i];
      }

      // Generate list of regions for the current result.
      regions.clear();
      for( ; solIt != solEnd; solIt++ )
      {
        actResInfo = (*solIt)->GetResultInfo();

        regionId = ptGrid_->RegionNameToId(
            (*solIt)->GetEntityList()->GetName());
        regions.push_back(++regionId);
      }

      // Reset solIt to beginning of result vector
      solIt = it->second.begin();

      try
      {
        numRegions = regions.size();

     

        std::vector< std::string > resNames;

        if( (*solIt)->GetEntryType() == EntryType::DOUBLE ) {
          resNames.push_back(resultName);
        }
        else
        {
          std::string resName;

          if(resInfo->complexFormat == REAL_IMAG)
          {
            resName = resultName; resName.append("_real");
            resNames.push_back(resName);
            resName = resultName; resName.append("_imag");
            resNames.push_back(resName);
          }
          else
          {
            resName = resultName; resName.append("_ampl");
            resNames.push_back(resName);
            resName = resultName; resName.append("_phase");
            resNames.push_back(resName);
          }
        }

        // Generate compound datatype
        H5IO::CompoundType resInfo;
        typedef std::pair<std::string, boost::any> CEntryType;
        resInfo.push_back( CEntryType( "DefinedOn", definedOn ) );
        resInfo.push_back( CEntryType( "Regions", regions ) );
        resInfo.push_back( CEntryType( "NumDOFs", numDOFs ) );
        resInfo.push_back( CEntryType( "DOFNames", dofNames ) );
        resInfo.push_back( CEntryType( "EntryType", entryType ) );
        resInfo.push_back( CEntryType( "Unit", unit ) );

        H5IO::WriteCompound( currAttrDescGroup_, resNames[0], resInfo );
        


      } catch (H5::Exception& h5ex2)
      {
        EXCEPTION(h5ex2.getCDetailMsg());
      }
    }
  }

  void SimOutputXMDF::WriteAttributes(const std::vector<std::string>& resultNames,
                                      std::vector< Vector<Double>* >& results,
                                      const std::string& regionStr,
                                      const UInt numDOFs)
  {
    H5::Group resultGroup;
    UInt numResults = resultNames.size();

    for(UInt i=0; i<numResults; i++)
    {
      UInt numEntities = results[i]->GetSize() / numDOFs;

      try
      {
        // Try to open result group.
        resultGroup = currStepGroup_.openGroup(resultNames[i]);
      }
      catch (H5::Exception& h5ex)
      {
        try
        {
          resultGroup = currStepGroup_.createGroup(resultNames[i]);
        }
        catch (H5::Exception& h5ex)
        {
          EXCEPTION(h5ex.getCDetailMsg());
        }
      }

      try
      {
        // dim sizes of ds (on disk)
        hsize_t fdim[] = {numEntities, numDOFs}; 
        H5::DataSpace fspace( 2, fdim );

        // Create dataset and write it into the file.
        H5::DataSet dataset;
        dataset = resultGroup.createDataSet( regionStr,
                                             H5::PredType::NATIVE_DOUBLE,
                                             fspace );
        dataset.write( &(*results[i])[0], H5::PredType::IEEE_F64LE);

        dataset.close();
        fspace.close();
      }
      catch (H5::Exception& h5ex)
      {
        EXCEPTION(h5ex.getCDetailMsg());
      }
    }

    resultGroup.close();
  }

  void SimOutputXMDF::CreateExternalFile()
  {
    H5::Attribute attrib;
    H5::DataSpace space;
    std::stringstream fName;
    std::string pathsep;
    std::string fn;
    std::stringstream masterGroup;
    std::string mg;

    try 
    {

      pathsep = fs::path("/").native_directory_string();

      fName << fileName_ << "_ms" << currMS_ << "_step"
          << currStep_ << ".h5";
      fn = fName.str();
      fName.str("");
      fName << dirName_ << pathsep << fn;

      currStepFile_ = H5::H5File(fName.str(), H5F_ACC_TRUNC);

      // Add the attribute StepFrequency to step.
      attrib = currStepGroup_.createAttribute("ExtHDF5FileName",
                                              H5::StrType(H5::PredType::C_S1,
                                                  fn.size()),
                                              space);
      attrib.write(H5::StrType(H5::PredType::C_S1,
                   fn.size()),
      fn.c_str());
      attrib.close();
      currStepGroup_.close();

      currStepGroup_ = currStepFile_.openGroup("/");

      fName.str("");
      fName << fileName_ << ".h5";
      fn = fName.str();

      // Add the attribute StepFrequency to step.
      attrib = currStepGroup_.createAttribute("MasterHDF5FileName",
                                              H5::StrType(H5::PredType::C_S1,
                                                  fn.size()),
                                              space);
      attrib.write(H5::StrType(H5::PredType::C_S1,
                   fn.size()),
      fn.c_str());
      attrib.close();

      masterGroup << "/Data/Volume/Multistep " << currMS_
          << "/Step " << currStep_;
      mg = masterGroup.str();

      // Add the attribute StepFrequency to step.
      attrib = currStepGroup_.createAttribute("MasterGroup",
                                              H5::StrType(H5::PredType::C_S1,
                                                  mg.size()),
                                              space);
      attrib.write(H5::StrType(H5::PredType::C_S1,
                   mg.size()),
      mg.c_str());
      attrib.close();
    }
    catch (H5::Exception& h5ex)
    {
      EXCEPTION(h5ex.getCDetailMsg());
    }
  }

  void SimOutputXMDF::WriteStringToUserData(const std::string& dSetName, 
                                            const std::string& str)
  {
    H5::Group userDataGroup;
    
    // If it does not exist, create Group for Data.
    try 
    {
      userDataGroup = dataGroup_.openGroup("User Data");
    }
    catch (H5::Exception& h5ex)
    {
      try 
      {
        userDataGroup = dataGroup_.createGroup("User Data");
      }
      catch (H5::Exception& h5ex)
      {
        EXCEPTION(h5ex.getCDetailMsg());
      }
    }

    H5IO::WriteAttribute( userDataGroup, dSetName, str );
    H5::DataSpace space;

    H5::DataSet dataset;
    dataset = userDataGroup.createDataSet(dSetName,
                                          H5::StrType(H5::PredType::C_S1,
                                                      str.length()),
                                          space);
    
    dataset.write(str.c_str(), H5::StrType(H5::PredType::C_S1, str.length()));
    
    userDataGroup.close();
  }

} // end of namespace

