// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>

#define BOOST_FILESYSTEM_VERSION 3

// #include "boost/filesystem.hpp"
// #include "boost/filesystem/operations.hpp"
// #include "boost/filesystem/path.hpp"
// #include "boost/filesystem/convenience.hpp"
// #include "boost/filesystem/exception.hpp"

#include "boost/algorithm/string/trim.hpp"
#include "boost/algorithm/string.hpp"
#include "boost/tokenizer.hpp"
#include "boost/lexical_cast.hpp"

#include "boost/version.hpp"

//
#include "hdf5Reader.h"

// namespace fs = boost::filesystem;

#ifdef WIN32
#define PATH_SEPARATOR '\\'
#else
#define PATH_SEPARATOR '/'
#endif

namespace H5CFS {


//#define H5_CATCH( STR )                                                 \
//  catch (H5::Exception& h5Ex ) {                                        \
//    H5CFS_EXCEPTION( STR << ":\n" << h5Ex.getCDetailMsg() );                  \
//  }

  Hdf5Reader::Hdf5Reader( ) 
  {

      statsRead_ = false;
      hasExternalFiles_ = false;

    // Do not print HDF5 exceptions by default
    H5::Exception::dontPrint();
  }

  Hdf5Reader::~Hdf5Reader() {
    if( mainFile_.getLocId() <= 0 )
      return;
    

    mainRoot_.close();
    mainFile_.close();
  }

  void Hdf5Reader::LoadFile( const std::string& fileName )
  {
    fileName_ = fileName;
    
    std::string::size_type len = fileName.rfind(PATH_SEPARATOR);

    baseDir_ = fileName.substr(0, len);
    

#if 0
    try 
    {
      //      fs::path fn = fs::system_complete(fileName_);
      //      fn.normalize();
      //      baseDir_ = fn.branch_path().native();
      //      baseName = (fs::change_extension( fn.leaf(), "" )).native();
      //      fileName_ = fn.native();
      //      if(fs::extension(fn) == "")
      //      {
      //   fn = fs::change_extension( fn, ".h5" );
      // }
    } catch (fs::filesystem_error& ex)
    {
      H5CFS_EXCEPTION("Received exception: " << ex.what());
      return;
    }
    ;
#endif

    //! open file and store main group and mesh group
    try {
      mainFile_ = H5::H5File( fileName_, H5F_ACC_RDONLY );
    } H5_CATCH( "Could not open HDF5 file '" << fileName_ << "'" );

    try {
      mainRoot_ = mainFile_.openGroup("/");
    } H5_CATCH( "Could not open main root" );
    
    try{
      meshRoot_ = mainRoot_.openGroup("Mesh");
    } H5_CATCH( "Could not open mesh group" );
    
    // check for use of external files
    try {
      H5::Group meshResGroup = mainRoot_.openGroup("Results/Mesh");
      Hdf5Common::ReadAttribute( meshResGroup, "ExternalFiles",
                           hasExternalFiles_ );
      meshResGroup.close();
    } catch( H5::Exception&) {
      hasExternalFiles_ = false;
    }

    // read general mesh information
    ReadMeshStats( meshRoot_ );
  }
/**
  void Hdf5Reader::ReadMesh( Grid *mi )
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
    if(*readEntities_.begin() == "all") {
      readEntities_.clear();
      readEntities_.insert(regionNames_.Begin(), regionNames_.End());
      readEntities_.insert(nodeNames_.Begin(), nodeNames_.End());
      readEntities_.insert(elemNames_.Begin(), elemNames_.End());
      LOG_DBG(Hdf5Reader) << "The following entities will be read:";
      std::stringstream sstr;
      it=readEntities_.begin();
      end=readEntities_.end();
      for( ; it != end; it++)
        sstr << (*it) << " ";
      LOG_DBG(Hdf5Reader) << sstr.str();
    }
    
    // Check if all readEntities_ can be found in file.
    it=readEntities_.begin();
    end=readEntities_.end();
    for( ; it != end; it++) {
      std::vector<std::string>::iterator findIt;
      
      findIt = std::find(regionNames_.Begin(), regionNames_.End(), *it); 
      if( findIt != regionNames_.End()) continue;
       
      findIt = std::find(nodeNames_.Begin(), nodeNames_.End(), *it); 
      if( findIt != nodeNames_.End()) continue;
       
      findIt = std::find(elemNames_.Begin(), elemNames_.End(), *it); 
      if( findIt != elemNames_.End()) continue;

      H5CFS_EXCEPTION("Entity " << (*it) << " specified for"
                " reading does not exist." );
   }

    // Make sure we have only entities in linearizeEntities_ which are
    // also part of readEntities_
    if(*linearizeEntities_.begin() == "none") {
      linearizeEntities_.clear();
    } else if(*linearizeEntities_.begin() == "all") {
      linearizeEntities_.insert(readEntities_.begin(), readEntities_.end());
    } else {
      it=linearizeEntities_.begin();
      end=linearizeEntities_.end();
      for( ; it != end; ) {
        if(readEntities_.find(*it) == readEntities_.end()) {
          erase = it; it++;
//          std::cout << "Erasing nonexistant entity " << (*erase) << std::endl;
          linearizeEntities_.erase(erase);
        }
        it++;
      }      
    }

    // Remove nodal entities from linearizeEntities_
    it=linearizeEntities_.begin();
    end=linearizeEntities_.end();
    for( ; it != end; ) {
      std::vector<std::string>::iterator findIt;
      
      findIt = std::find(nodeNames_.Begin(), nodeNames_.End(), *it); 
      if( findIt != nodeNames_.End()) {
        erase = it; it++;
//        std::cout << "Erasing nodal entity " << (*erase) << std::endl;
        linearizeEntities_.erase(erase);
      }
      it++;
    }      

//  TODO: strieben - Remove these lines!   
//    std::cout << "linearizeEntities" << std::endl;
//    typedef std::ostream_iterator<std::string> string_os_iter;
//    std::copy (linearizeEntities_.begin(),
//               linearizeEntities_.end(),
//               string_os_iter (std::cout, " "));
//    std::cout << std::endl;
   
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

**/
  // ======================================================
  //  GENERAL MESH INFORMATION
  // ======================================================
  unsigned int Hdf5Reader::GetDim() {

    // Open mesh group
    H5::Group meshGroup;
    try{
      meshGroup = mainRoot_.openGroup("Mesh");
    } H5_CATCH( "Could not open mesh group" );
    
    // Read dimension
    unsigned int dim;
    Hdf5Common::ReadAttribute( meshGroup, "Dimension", dim );
    meshGroup.close();
    return dim;
  }
  
  unsigned int Hdf5Reader::GetNumNodes(){
    return numNodes_;
  }
    
  unsigned int Hdf5Reader::GetNumElems(const int dim){
    return numElems_;
  }
  
  unsigned int Hdf5Reader::GetNumRegions(){
    return regionNames_.size();
  }

  unsigned int Hdf5Reader::GetNumNamedNodes(){
    return nodeNames_.size();
  }

  unsigned int Hdf5Reader::GetNumNamedElems(){
    return elemNames_.size();
  }
  
  // ======================================================
  // ENTITY NAME ACCESS
  // ======================================================

  void Hdf5Reader::GetAllRegionNames( std::vector<std::string> & regionNames ){
    regionNames = regionNames_;
  }

  void Hdf5Reader::GetRegionNamesOfDim( std::vector<std::string> & regionNames,
                                   const unsigned int dim ) const
  {
    for (unsigned int i = 0; i < regionNames_.size(); ++i)
    {
      if (regionDims_.find(regionNames_[i])->second == dim)
      {
        regionNames.push_back(regionNames_[i]);
      }
    }
  }
    
  void Hdf5Reader::GetNodeNames( std::vector<std::string> & nodeNames )
  {
    nodeNames = nodeNames_;
  }
  
  void Hdf5Reader::GetElemNames( std::vector<std::string> & elemNames )
  {
    elemNames = elemNames_;
  }


  // =======================================================================
  //  MESH ENTITY ACCESS
  // =======================================================================
  
  void Hdf5Reader::GetNodeCoords( std::vector<std::vector<double> >& coords ) {
    
    // open node group
    H5::Group nodeGroup;
    try{ 
      nodeGroup = meshRoot_.openGroup( "Nodes") ;
    } H5_CATCH( "Could not open Elements / Nodes group" );

    // read nodal coordinates
    std::vector<double> coordVec ;
    unsigned int numNodes = Hdf5Common::GetArrayDims( nodeGroup, "Coordinates" )[0];
    Hdf5Common::ReadArray( nodeGroup, "Coordinates", coordVec );

    coords.resize( numNodes );
    unsigned int index = 0;
    for( unsigned int i = 0; i < numNodes; i++ ) {
      coords[i].resize(3);
      coords[i][0] = coordVec[index];
      coords[i][1] = coordVec[index+1];
      coords[i][2] = coordVec[index+2];
      index += 3;
    }
    
    nodeGroup.close();
  }


  void Hdf5Reader::GetNodeCoord( unsigned int nodeNum, std::vector<double>& coord ) {
    H5CFS_EXCEPTION( "Not yet implemented" );
  }

  void Hdf5Reader::GetElems( std::vector<ElemType>& types, 
                             std::vector<std::vector<unsigned int> >& connect ) {
    H5::Group elemGroup;
    try {
      elemGroup = meshRoot_.openGroup( "Elements" );
    } H5_CATCH( "Could not open 'Elements' group" );
    
    // read number of elements
    //unsigned int numElems = 0;
    //Hdf5Common::ReadAttribute( elemGroup, "NumElems", numElems );
    unsigned int numElems = 
      Hdf5Common::GetArrayDims( elemGroup, "Connectivity")[0];
      
    // read maximum number of nodes per elements
    unsigned int numNodesPerElem = 
      Hdf5Common::GetArrayDims( elemGroup, "Connectivity")[1];

    // read element types
    std::vector<int> elemTypes;
    Hdf5Common::ReadArray( elemGroup, "Types", elemTypes );

    // read nodes per element
    std::vector<unsigned int> globConnect;
    Hdf5Common::ReadArray( elemGroup, "Connectivity", globConnect );
    
    // add element definition
    types.resize( numElems );
    connect.resize( numElems );
    std::vector<unsigned int>::const_iterator it1, it2;
    it1 = it2 = globConnect.begin();
    for( unsigned int i = 0; i < numElems; i++ ) {
      unsigned numElemNodes = NUM_ELEM_NODES[elemTypes[i]];
      it2 = it1 + numElemNodes;
      connect[i] = std::vector<unsigned int>(it1, it2);
      types[i] = (ElemType) elemTypes[i];
      it1 += numNodesPerElem;
    }

    elemGroup.close();
  }
    
  void Hdf5Reader::GetElem( unsigned int elemNum, ElemType& type, 
                            std::vector<unsigned int>& connect ) {
    H5CFS_EXCEPTION( "Not yet implemented" );                        
  }
  
  void Hdf5Reader::GetElemsOfRegion( const std::string& regionName,
                                     std::vector<unsigned int>& elemNums ) {
    if(std::find(regionNames_.begin(), regionNames_.end(), regionName) != regionNames_.end()) 
    {
      elemNums = regionElems_[regionName];
    }
    else
    {
      H5CFS_EXCEPTION("No elements present for region '" << regionName << "'.");
    }
  }
  
  void Hdf5Reader::GetNodesOfRegion( const std::string& regionName,
                                     std::vector<unsigned int>& nodeNums ) {
    if(std::find(regionNames_.begin(), regionNames_.end(), regionName) != regionNames_.end()) 
    {
      nodeNums = regionNodes_[regionName];
    }
    else
    {
      H5CFS_EXCEPTION("No nodes present for region '" << regionName << "'.");
    }
  }


  // =========================================================================
  //  GENERAL SOLUTION INFORMATION
  // =========================================================================

  void Hdf5Reader::
  GetNumMultiSequenceSteps( std::map<unsigned int,AnalysisType>& analysis,
                            std::map<unsigned int, unsigned int>& numSteps,
                            bool isHistory ) {

    H5::Group resultGroup, actMsGroup;
    std::string actAnalysisString;
    analysis.clear();
    numSteps.clear();
    unsigned int actMsNumSteps = 0;
     
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
    //std::cerr << "trying to get number of multisequence steps ... \n";
    
    // Iterate over all children in the specific group and collect the stepvalues
    hsize_t numChildren = resultGroup.getNumObjs();
    std::set<unsigned int> msStepNums;
    for( unsigned int i = 0; i < numChildren; i++ ) {
      std::string actName = Hdf5Common::GetObjNameByIdx( resultGroup, i );
      
      // cut away "MultiStep_"-substring and convert  into int
      boost::erase_all(actName, "MultiStep_");
      msStepNums.insert( boost::lexical_cast<unsigned int>(actName) );
    }
    
    // try to find all single multisequence steps and related analysis string
    std::set<unsigned int>::iterator it;
    for( it = msStepNums.begin(); it != msStepNums.end(); it++ ) {
      actMsGroup = Hdf5Common::GetMultiStepGroup( mainFile_, *it, isHistory );

      // get analyisstring
      Hdf5Common::ReadAttribute( actMsGroup, "AnalysisType", actAnalysisString );
      Hdf5Common::ReadAttribute( actMsGroup, "LastStepNum", actMsNumSteps );
      AnalysisType actAnalysis;
      if( actAnalysisString == "static") {
        actAnalysis = STATIC;
      } else if( actAnalysisString == "transient" ) {
        actAnalysis = TRANSIENT;
      } else if( actAnalysisString == "harmonic" ) {
        actAnalysis = HARMONIC;
      } else if( actAnalysisString == "eigenFrequency" ) {
        actAnalysis = EIGENFREQUENCY;
      } else {
        H5CFS_EXCEPTION( "Unknown analysistype found in hdf file");
      }
      analysis[*it] = actAnalysis;
      numSteps[*it] = actMsNumSteps;
      
      actMsGroup.close();
    }
    resultGroup.close();
  }

  
  void Hdf5Reader::
  GetStepValues( unsigned int sequenceStep,
                 shared_ptr<ResultInfo> info,
                 std::map<unsigned int, double>& steps,
                 bool isHistory ) { 

    
    // open corresponding multistep group
    H5::Group actMsGroup = Hdf5Common::GetMultiStepGroup( mainFile_, 
                                                    sequenceStep, isHistory  );
    
    // open result description
    H5::Group resGroup;
    try {
    resGroup = actMsGroup.openGroup("ResultDescription").
               openGroup( info->name );
    } H5_CATCH( "Could not open resultdescription for result '"
                 << info->name << "'" );
    
    // read stepValues and stepNumbers
    std::vector<double> values;
    std::vector<unsigned int> numbers;
    Hdf5Common::ReadArray( resGroup, "StepNumbers", numbers );
    Hdf5Common::ReadArray( resGroup, "StepValues", values );
    
    // sanity check: both vectors need to have the same dimension
    if( values.size() != numbers.size() ) {
      H5CFS_EXCEPTION( "There are not as many stepnumbers as stepvalues" );
    }
    
    // copy to steps-array
    // make it robust to handle old optimization files which are corrupt
    // stop reading if the value is not larger than the last
    steps.clear();
    bool corrupt = false;
    for( unsigned int i = 0, n=numbers.size(); !corrupt && i < n; i++ )
    {
      if(!corrupt && (i == 0 || values[i] > values[i-1]))
        steps[numbers[i]] = values[i];
      else
      {
        std::cout << "======== Assume corrupt data and cancel reading ==== " << std::endl
                  << "last good was step " << i << " with stepnumer " << numbers[i-1]
                  << " and value " << values[i-1] << std::endl
                  << "corrupt is stepnumber " << numbers[i] << " and value " << values[i] << std::endl
                  << "ignore a total of " << (n-i) << " steps!" << std::endl;
        corrupt = true;
      }
    }

    resGroup.close();
    actMsGroup.close();
  }


  void Hdf5Reader::
  GetResultTypes( unsigned int sequenceStep, 
                  std::vector<shared_ptr<ResultInfo> >& infos,
                  bool isHistory ) {
    
    // open ms group and 'Result Description' subgroup
    H5::Group actMsGroup = Hdf5Common::GetMultiStepGroup( mainFile_, sequenceStep,
                                                    isHistory );
    H5::Group resInfoGroup;
    try {
      resInfoGroup = actMsGroup.openGroup( "ResultDescription" );
    } H5_CATCH( "Could not open group 'ResultDescription'" );
    unsigned int numResults = static_cast<unsigned int>( resInfoGroup.getNumObjs() );
    
    // iterate over all entries and assemble the resultinfo object
    H5::Group actResInfoGroup;
    std::string actResultName;

    infos.clear();
    for( unsigned int i = 0; i < numResults; i++ ) {
      actResultName = Hdf5Common::GetObjNameByIdx( resInfoGroup, i );
      try{
        actResInfoGroup = resInfoGroup.openGroup( actResultName );
      } H5_CATCH( "Could not open description group for result '"
                  << actResultName << "'" );

      // Read resultinfo data
      unsigned int definedOn, numDOFs, entryType;
      std::string unit;
      std::vector<std::string> dofNames, entities;
      Hdf5Common::ReadArray( actResInfoGroup, "DefinedOn", &definedOn );
      Hdf5Common::ReadArray( actResInfoGroup, "NumDOFs", &numDOFs );
      Hdf5Common::ReadArray( actResInfoGroup, "DOFNames", dofNames );
      Hdf5Common::ReadArray( actResInfoGroup, "EntityNames", entities );
      Hdf5Common::ReadArray( actResInfoGroup, "EntryType", &entryType );
      Hdf5Common::ReadArray( actResInfoGroup, "Unit", &unit );

      // create new ResultInfo objects
      shared_ptr<ResultInfo> ptInfo( new ResultInfo() );
      ptInfo->name = actResultName;
      ptInfo->dofNames = dofNames;
      ptInfo->unit = unit;
      ptInfo->entryType = (EntryType) entryType;
      ptInfo->listType = (EntityType)  definedOn;
      
      // perform consistency check
      if( ptInfo->entryType == UNKNOWN ) {
        H5CFS_EXCEPTION( "Result '" << actResultName 
                   << "' has no proper EntryType!" );
      }
      
      if( ptInfo->dofNames.size() == 0 ) {
        H5CFS_EXCEPTION( "Result '" << actResultName 
                   << "' has no degrees of freedoms!");
      }
      
      if( ptInfo->name == "" ) {
        H5CFS_EXCEPTION( "Result has neither a name nor a "
                   << "proper result type!" );
      }
      
      // now, run over all regions and create for each region
      // a new result info
      for( unsigned int iRegion = 0; iRegion < entities.size(); iRegion++) {
        shared_ptr<ResultInfo> actInfo(new ResultInfo(*ptInfo));
        actInfo->listName = entities[iRegion];
        infos.push_back( actInfo );
      }
      
    }
    
    resInfoGroup.close();
    actMsGroup.close();
  }



  void Hdf5Reader::GetEntities( EntityType type, const std::string& name,
                              std::vector<unsigned int>& entities ) {
    if( type == NODE ) {
      GetNodesOfRegion( name, entities );
    }
    if( type == ELEMENT || type == SURF_ELEM ) {
      GetElemsOfRegion( name, entities );
    }
  }
  

/**
  void Hdf5Reader::
  GetResultEntities( unsigned int sequenceStep,
                     shared_ptr<ResultInfo> info,
                     std::vector<shared_ptr<EntityList> >& list,
                     bool isHistory ) {

    // get resultname from resultinfo object
    std::string resultName = info->resultName;
    
    // open ms group and specific entry in 'ResultDescription'
    H5::Group actMsGroup = Hdf5Common::GetMultiStepGroup( mainFile_, sequenceStep,
                           isHistory );
    H5::Group resInfoGroup;
    try {
      resInfoGroup = actMsGroup.openGroup( "ResultDescription/" + resultName );
    } H5_CATCH( "Could not open group result description for result "
                << resultName );
    
    // get regions
    std::vector<std::string> regions;
    Hdf5Common::ReadArray( resInfoGroup, "EntityNames", regions );

    // determine type of list for this result
    EntityList::ListType listType;
    EntityList::DefineType defineType;
    switch( info->definedOn ) {
    case ResultInfo::NODE:
    case ResultInfo::PFEM:
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
      H5CFS_EXCEPTION( "Only results defined on nodes and elements "
                 << "can be read in from HDF5 file up to now" );
    }

    // iterate over all regions
    list.Clear();
    for( unsigned int i = 0; i < regions.GetSize(); i++ ) {
      list.Push_back( mi_->GetEntityList( listType, regions[i], 
                                          defineType ) );
    }
    resInfoGroup.close();
    actMsGroup.close();
  }
  **/
  void Hdf5Reader::GetResult( unsigned int sequenceStep,
                              unsigned int stepNum,
                              shared_ptr<Result> result,
                              bool isHistory ) {
    
    if( !isHistory ) {
      GetMeshResult( sequenceStep, stepNum, result);
    } else {
      GetHistResult( sequenceStep, stepNum, result );      
    }
    
  }
  void Hdf5Reader::GetMeshResult( unsigned int sequenceStep, unsigned int stepNum,
                                  shared_ptr<Result> result ) {
                                     
    // open stepgroup, open specific result subgroup
    H5::Group stepGroup = Hdf5Common::GetStepGroup( mainFile_, sequenceStep, 
                                              stepNum );
    
    // check, if results are stored at external file location
    H5::H5File extFile;
    if( hasExternalFiles_ ) {
      std::string extFileString;
      Hdf5Common::ReadAttribute( stepGroup, "ExtHDF5FileName", extFileString);
      std::string extFileNameComplete = baseDir_ + PATH_SEPARATOR + extFileString;
      try {
        extFile = H5::H5File( extFileNameComplete, H5F_ACC_RDONLY );
      } H5_CATCH( "Could not open external file '" 
                  << extFileString << "' for result '" 
                  << result->resultInfo->name
                  << "' in multisequence step " << sequenceStep
                  << ", analysis step " << stepNum  );

      // replace old step group by new one
      stepGroup.close();
      stepGroup = extFile.openGroup( "/" );
    }
    
    // determine region for this results
    std::string regionName =  result->resultInfo->listName;

    // determine entity type string
    std::string entString;
    switch( result->resultInfo->listType ) {
    case NODE:
      entString = "Nodes";
      break;
    case ELEMENT:
    case SURF_ELEM:
      entString = "Elements";
      break;
    default:
      H5CFS_EXCEPTION( "Currently only results on nodes and elements "
                 << "can be read in from a hdf5 file ");
    }

    std::string groupName = result->resultInfo->name;
    groupName += "/" + regionName + "/" + entString;
    
    H5::Group resGroup;
    try {
      resGroup = stepGroup.openGroup( groupName );
    } H5_CATCH( "Unable to open group for result '" 
                << result->resultInfo->name
                << "' on '" << regionName << "' in step " << stepNum );

    // read data array
    std::vector<double> realVals;
    Hdf5Common::ReadArray( resGroup, "Real", realVals );

    std::vector<unsigned int> idx, entities;
    unsigned int numDofs = result->resultInfo->dofNames.size();
    this->GetEntities( result->resultInfo->listType, result->resultInfo->listName,
                       entities );
    unsigned int numEntities = entities.size(); 
    unsigned int resVecSize =  numEntities * numDofs;    
    
    // copy data array to result object
    // REAL part
    std::vector<double> & resRealVec= result->realVals;
    resRealVec.resize( resVecSize );
    for( unsigned int i = 0; i < numEntities; i++ ) {
      for( unsigned int iDof = 0; iDof < numDofs; iDof++ ) {
        resRealVec[i*numDofs+iDof] = realVals[i*numDofs+iDof];
      }
    }
    
    // check if also imaginary values are present
    if( resGroup.getNumObjs() > 1 ) {
      result->isComplex = true;
      std::vector<double> imagVals;
      Hdf5Common::ReadArray( resGroup, "Imag", imagVals );
      std::vector<double> & resImagVec= result->imagVals;
      resImagVec.resize( resVecSize );
      for( unsigned int i = 0; i < numEntities; i++ ) {
        for( unsigned int iDof = 0; iDof < numDofs; iDof++ ) {
          resImagVec[i*numDofs+iDof] = imagVals[i*numDofs+iDof];
        }
      }
    } else {
      result->isComplex = false;
    }
      
    resGroup.close();
    stepGroup.close();
    // close external file for current step
    if( hasExternalFiles_ )
      extFile.close();
  }

  void Hdf5Reader::GetHistResult( unsigned int sequenceStep, 
                                  unsigned int stepNum,
                                  shared_ptr<Result> result )
  {
  }

  /**
  void Hdf5Reader::GetHistResult( unsigned int sequenceStep, unsigned int stepNum,
                                         shared_ptr<BaseResult> result ) {

    // open multisequence group
    H5::Group actMsGroup = Hdf5Common::GetMultiStepGroup( mainFile_, sequenceStep,
                                                    true);

    // open group for specific result
    std::string resultName = result->GetResultInfo()->resultName;
    H5::Group actResGroup = actMsGroup.openGroup( resultName );

    // determine from definedOn type the correct string representation
    // of the subgroup
    ResultInfo::EntityUnknownType definedOn = result->GetResultInfo()->definedOn;
    std::string entityTypeString = Hdf5Common::MapUnknownTypeAsString( definedOn );
    H5::Group entityGroup = actResGroup.openGroup( entityTypeString );

    // iterate over all entities in the the list
    shared_ptr<EntityList> list = result->GetEntityList();
    EntityIterator it = list->GetIterator();
    unsigned int numDofs = result->GetResultInfo()->dofNames.GetSize();
    
    // ----------------------------
    // Case 1: Real valued entries
    // ----------------------------
    if( result->GetEntryType() == EntryType::double ) {
      std::vector<double> & resVec = 
        dynamic_cast<Result<double>&>(*result).Getstd::vector();
      resVec.Resize( list->GetSize() * numDofs );
      for( it.Begin(); !it.IsEnd(); it++ ) {
        
        H5::Group actEntGroup;
        
        if(entityTypeString == "Nodes") {
          // Find node number in HDF5 file which corresponds to mesh node number
          unsigned int meshNodeNum = it.GetNode();
          unsigned int fileNodeNum = g2FNodeNumMap_[meshNodeNum];
          std::stringstream sstr;
          sstr << fileNodeNum;
          actEntGroup = entityGroup.openGroup(sstr.str());
        } else if(entityTypeString == "Elements") {
          // Find element number in HDF5 file which corresponds to mesh elem number
          unsigned int meshElemNum = it.GetElem()->elemNum;
          unsigned int fileElemNum = g2FElemNumMap_[meshElemNum];
          std::stringstream sstr;
          sstr << fileElemNum;
          actEntGroup = entityGroup.openGroup(sstr.str());
        } else {
          // open for each entity the specific subgroup
          actEntGroup = entityGroup.openGroup(it.GetIdString() );
        }

        // read single part of array and set it in the result vector
        std::vector<double> vals;
        Hdf5Common::ReadArray( actEntGroup, "Real", vals );

        for( unsigned int i = 0; i < numDofs; i++ ) {
          resVec[it.GetPos()*numDofs+i] = vals[(stepNum-1)*numDofs+i];
        }
        actEntGroup.close();
      }
    } 
    // ----------------------------
    // Case 2: Complex valued entries
    // ----------------------------
    else {
      std::vector<Complex> & resVec = 
        dynamic_cast<Result<Complex>&>(*result).Getstd::vector();
      resVec.Resize( list->GetSize() * numDofs );
      for( it.Begin(); !it.IsEnd(); it++ ) {
        // open for each entity the specific subgroup
        H5::Group actEntGroup = 
          entityGroup.openGroup(it.GetIdString() );

        // read single part of array and set it in the result vector
        std::vector<double> realVals, imagVals;
        Hdf5Common::ReadArray( actEntGroup, "Real", realVals );
        Hdf5Common::ReadArray( actEntGroup, "Imag", imagVals );

        for( unsigned int i = 0; i < numDofs; i++ ) {
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
  void Hdf5Reader::ReadNodeElemData(const H5::Group& meshGroup)
  {

    // ==================================
    //  Read General Element Information
    // ==================================
    H5::Group elemGroup;
    try {
      elemGroup = meshGroup.openGroup( "Elements" );
    } H5_CATCH( "Could not open 'Elements' group" );

    
    // read number of elements
    unsigned int numElems = 0;
    Hdf5Common::ReadAttribute( elemGroup, "NumElems", numElems );

    // read maximum number of nodes per elements
    unsigned int numNodesPerElem = 
      Hdf5Common::GetArrayDims( elemGroup, "Connectivity")[1];
    
    // read element types
    std::vector<int> elemTypes;
    Hdf5Common::ReadArray( elemGroup, "Types", elemTypes );
    
    // read nodes per element
    std::vector<unsigned int> globConnect;
    Hdf5Common::ReadArray( elemGroup, "Connectivity", globConnect );
    
    elemGroup.close();
    
    typedef std::set<unsigned int> RegionSetType;
    RegionSetType readNodeSet;
    RegionSetType readElemSet;
    
    // map for each element number the related region
    std::map<unsigned int, RegionIdType> elemRegionMap;
    
    // ensure, that region names are already read in
    if( !statsRead_ )
      ReadMeshStats( meshGroup );

    std::set<std::string>::iterator findIt;
    
    // ================================
    //  Add Elements Per Element Group
    // ================================
    for( unsigned int i = 0; i < elemNames_.GetSize(); i++ ) {
      H5::Group entityGroup, actEntityGroup;

      try{
        entityGroup = meshGroup.openGroup( "Groups" );
      } catch (H5::Exception&) {
        std::cout << "Could not open group for entity groups";
        return;
      }

      for( unsigned int i = 0, n=elemNames_.GetSize(); i < n; i++ ) {
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
        std::vector<unsigned int> readElems;
        Hdf5Common::ReadArray( actEntityGroup, "Elements", readElems );
        for( unsigned int j = 0, n2=readElems.GetSize(); j < n2; j++ ) {
          readElemSet.insert( readElems[j] );
          elemRegionMap[readElems[j]] = NO_REGION_ID;
        }

        std::vector<unsigned int> readNodes;
        std::vector<unsigned int> actualNodes;
        std::vector<unsigned int> nodeIndices;
        // read nodes from grid
        Hdf5Common::ReadArray( actEntityGroup, "Nodes", readNodes );
        actualNodes=readNodes;
        
        // Extract nodes needed for linear elements
        if(linearizeEntity) {
          LinearizeElems(readElems, elemTypes, globConnect, actualNodes);
          
          // To be able to read results for linearized grids we need an index
          // map for each entity (nodes of named elems and regions)
          nodeIndices.Resize(actualNodes.GetSize());
          for( unsigned int j = 0, n2=actualNodes.GetSize(); j < n2; j++ ) {
            nodeIndices[j] = readNodes.Find(actualNodes[j]);
          }          
        }
        
        // For element lists which need not be linearized we simply use
        // ascending node indices.
        nodeIndices.Resize(actualNodes.GetSize());
        for( unsigned int j = 0, n2=actualNodes.GetSize(); j < n2; j++ ) {
          nodeIndices[j] = j;
        }
        
        entityNodeMap_[elemNames_[i]] = nodeIndices;
        
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


    // Read all nodes from regions and initialize mapping from mesh node
    // numbers to grid node numbers accordingly.
    unsigned int baseNodeNum = mi_->GetNumNodes() + 1;
    unsigned int baseElemNum = mi_->GetNumElems() + 1;

    for( unsigned int i = 0, n=regionNames_.GetSize(); i < n; i++ ) {
      std::string regionName = regionNames_[i];
      findIt = readEntities_.find(regionName);
      if(findIt == readEntities_.end())
        continue;
      
      try {
        actRegion = regionGroup.openGroup( regionName );
      } H5_CATCH( "Could not open group for region '" <<
                  regionNames_[i] << "'" );

      
      // pass region names to grid and obtain RegionIds
      RegionIdType actRegionId;
      mi_->AddRegion(regionName, actRegionId);

      // Check if entity needs to be linearized
      bool linearizeEntity;
      linearizeEntity = linearizeEntities_.find(regionName) != linearizeEntities_.end();

      // read elem numbers for this region
      std::vector<unsigned int> regionElems;
      Hdf5Common::ReadArray( actRegion, "Elements", regionElems );
      for( unsigned int j = 0; j < regionElems.GetSize(); j++ ) {
        readElemSet.insert( regionElems[j] );
        elemRegionMap[regionElems[j]] = actRegionId;
      }
      
      std::vector<unsigned int> regionNodes;
      std::vector<unsigned int> actualNodes;
      std::vector<unsigned int> nodeIndices;

      Hdf5Common::ReadArray( actRegion, "Nodes", regionNodes );
      actualNodes=regionNodes;

      if(linearizeEntity) {
        LinearizeElems(regionElems, elemTypes, globConnect, actualNodes);
        
        // To be able to read results for linearized grids we need an index
        // map for each entity (nodes of named elems and regions)
        nodeIndices.Resize(actualNodes.GetSize());
        for( unsigned int j = 0, n2=actualNodes.GetSize(); j < n2; j++ ) {
          nodeIndices[j] = regionNodes.Find(actualNodes[j]);
        }
      }

      // For element list which need not be linearized we simply use
      // ascending node indices.
      nodeIndices.Resize(actualNodes.GetSize());
      for( unsigned int j = 0, n2=actualNodes.GetSize(); j < n2; j++ ) {
        nodeIndices[j] = j;
      }
      
      entityNodeMap_[regionName] = nodeIndices;
      
      readNodeSet.insert( regionNodes.Begin(), regionNodes.End());
    }
    
    // ===================
    //  Add Nodes to Grid
    // ===================
    mi_->AddNodes( readNodeSet.size() );
    unsigned int idx;
    std::vector<double> p(3);

    RegionSetType::iterator it, end;
    for( it=readNodeSet.begin(), end=readNodeSet.end();
         it != end;
         it++ ) {

      idx = ((*it)-1)*3;
      p[0] = nodeCoords_[idx + 0];
      p[1] = nodeCoords_[idx + 1];
      p[2] = nodeCoords_[idx + 2];
      mi_->SetNodeCoordinate( baseNodeNum, p );

      g2FNodeNumMap_[baseNodeNum] = *it;
      f2GNodeNumMap_[*it] = baseNodeNum++;
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
      g2FElemNumMap_[idx] = *it;
      f2GElemNumMap_[*it] = idx++;
    }

    // Remap global connectivity from mesh nodes to grid node numbers
    for( unsigned int i = 0, n=globConnect.GetSize();
         i < n;
         i++ ) {
      globConnect[i] = f2GNodeNumMap_[globConnect[i]];
    }

    // iterate over all elements
    for( it = readElemSet.begin(), end = readElemSet.end();
         it != end;
         it++ ) {
      unsigned int elemNum=(*it);
      FEType type = (FEType) elemTypes[elemNum-1];
      unsigned int * connect = &globConnect[numNodesPerElem*(elemNum-1)];
      RegionIdType actRegionId = elemRegionMap[elemNum];
      mi_->SetElemData( f2GElemNumMap_[elemNum], type, actRegionId, connect );
    }
    
    // ======================================
    //  Generate Node Groups For Each Region
    // ======================================
    
    // check, if nodes nodes of region should be additionally added
    // to list of named nodes
    if( genRegionNodes_) {    
#if 0      
      // iterate over all regions
      for( unsigned int i = 0; i < readEntities_.GetSize(); i++ ) {

        try {
          actRegion = regionGroup.openGroup( readEntities_[i] );
        } H5_CATCH( "Could not open group for region '" <<
                    regionNames_[i] << "'" );

        // read nodes of region
        std::vector<unsigned int> regionNodes;
        Hdf5Common::ReadArray( actRegion, "Nodes", regionNodes );

        for( unsigned int j = 0, n=regionNodes.GetSize();
        j < n;
        j++ ) {
          regionNodes[j] = f2GNodeNumMap_[regionNodes[j]];
        }
        // add nodes as named nodes
        mi_->AddNamedNodes( readEntities_[i]+"_Nodes", regionNodes );
        
        regionGroup.close();
      }
#else
      H5CFS_EXCEPTION("Generation of nodes per region temporarily not supported!");
#endif
    }

  }
*/
  void Hdf5Reader::GetNamedNodes(const std::string& name,
                                 std::vector<unsigned int>& nodeNums)
  {
    if(std::find(nodeNames_.begin(), nodeNames_.end(), name) != nodeNames_.end())
    {
      nodeNums = entityNodes_[name];
    }
    else
    {
      H5CFS_EXCEPTION("No nodes present for named entity '" << name << "'.");
    }
  }
  
                           
  void Hdf5Reader::GetNamedElems(const std::string& name,
                                 std::vector<unsigned int>& elemNums)
  {
    if(std::find(elemNames_.begin(), elemNames_.end(), name) != elemNames_.end())
    {
      elemNums = entityElems_[name];
    }
    else
    {
      H5CFS_EXCEPTION("No elements present for named entity '" << name << "'.");
    }
  }
  
  void Hdf5Reader::ReadMeshStats(const H5::Group& meshGroup) {


    // ========================
    //  READ NODAL INFORMATION
    // ========================
    H5::Group nodeGroup;
    // get the number of nodes
    try{ 
      nodeGroup = meshRoot_.openGroup( "Nodes") ;
    } H5_CATCH( "Could not open Elements / Nodes group" );
    
    numNodes_ = Hdf5Common::GetArrayDims( nodeGroup, "Coordinates" )[0];

    // ==========================
    //  READ ELEMENT INFORMATION
    // ==========================
    H5::Group elemGroup;
    // get the number of nodes
    try{ 
      elemGroup = meshRoot_.openGroup( "Elements") ;
    } H5_CATCH( "Could not open Elements / Nodes group" );
    
    numElems_ = Hdf5Common::GetArrayDims( elemGroup, "Connectivity" )[0];

    // ==================================
    //  Read Region Names and Dimensions
    // ==================================
    H5::Group regionGroup;    
    try {
      regionGroup = meshGroup.openGroup("Regions");
    } H5_CATCH( "Could not open 'Regions' subgroup" );

    regionNames_.clear();
    
    // Collect all nodes elements which belong to regions
    std::set<unsigned> nodesInRegions;
    std::set<unsigned> elemsInRegions;

    // iterate over all region names
    hsize_t numRegions = regionGroup.getNumObjs();
    for( hsize_t i = 0; i < numRegions; i++ ) {

      // get name
      std::string actName = Hdf5Common::GetObjNameByIdx( regionGroup, i );
      regionNames_.push_back(  actName );

      // get dimension
      unsigned int dim = 0;
      H5::Group actRegion = regionGroup.openGroup( actName );
      Hdf5Common::ReadAttribute( actRegion, "Dimension", dim );
      regionDims_[actName] = dim;
      // read elem numbers for this region
      Hdf5Common::ReadArray( actRegion, "Elements", regionElems_[actName] );
      // read node numbers for this region
      Hdf5Common::ReadArray( actRegion, "Nodes", regionNodes_[actName] );
      elemsInRegions.insert(regionElems_[actName].begin(),
                            regionElems_[actName].end());

      nodesInRegions.insert(regionNodes_[actName].begin(),
                            regionNodes_[actName].end());

      actRegion.close();
    }

    regionGroup.close();

    // Get all elements which are not part of a region (for NACS)
    unsigned int numElemsInRegions = elemsInRegions.size();
    std::vector<unsigned> regionLessElems;
    if( numElemsInRegions < numElems_ ) 
    {
      std::set<unsigned> allElems;
      for(unsigned i=0; i<numElems_; i++) allElems.insert(i+1);
      regionLessElems.resize(numElems_ - numElemsInRegions);
      std::vector<unsigned>::iterator it;

      it = std::set_difference (allElems.begin(),
                                allElems.end(), 
                                elemsInRegions.begin(),
                                elemsInRegions.end(),
                                regionLessElems.begin());

      regionNames_.push_back("OtherEntities");
      regionElems_["OtherEntities"] = regionLessElems;
      regionDims_["OtherEntities"] = 0;
    }

    // ==============================
    //  Read Named Nodes Description
    // ==============================
    H5::Group entityGroup;
    try {
      entityGroup = meshGroup.openGroup("Groups");
    } catch (H5::Exception&) {
      std::cout << "No node / elem groups present";
      statsRead_ = true;
      return;
    }
    
    // =================================
    //  Read node / element group names
    // =================================
    nodeNames_.clear();
    elemNames_.clear();
    std::set<unsigned> nodesInEntities;
    
    hsize_t numGroups = entityGroup.getNumObjs();
    for( hsize_t i = 0; i < numGroups; i++ ) {
      
      // get name of group
      std::string actName = Hdf5Common::GetObjNameByIdx( entityGroup, i );
      
      // open entitygroup and get number of different entity types
      // (nodes, elements) it is defined on
      H5::Group actEntityGroup;
      try {
        actEntityGroup = entityGroup.openGroup( actName );
      } H5_CATCH( "Could not open entity group '" << actName << "'");
      
      hsize_t numTypes = actEntityGroup.getNumObjs();
      bool hasElems = false;
      for( hsize_t iType = 0; iType < numTypes; iType++ ) {
        std::string actType = Hdf5Common::GetObjNameByIdx( actEntityGroup, iType );
        if( actType == "Elements" ) {
          hasElems = true;  
          break;
        }
      } // loop over types
      
      if( hasElems )  {
        elemNames_.push_back( actName );

        // read nodes and elements from file and add to grid
        Hdf5Common::ReadArray( actEntityGroup, "Nodes", entityNodes_[actName] );
        Hdf5Common::ReadArray( actEntityGroup, "Elements", entityElems_[actName] );

        if( numElemsInRegions < numElems_ ) 
        {
          std::vector<unsigned> v(entityElems_[actName].size());
          std::vector<unsigned>::iterator it;

          it = std::set_intersection (regionLessElems.begin(),
                                      regionLessElems.end(), 
                                      entityElems_[actName].begin(),
                                      entityElems_[actName].end(),
                                      v.begin());

          if(std::distance(v.begin(), it) > 0) 
          {
            unsigned dim = 0;
            unsigned regDim = regionDims_["OtherEntities"];

            Hdf5Common::ReadAttribute( actEntityGroup, "Dimension", dim );
            regionDims_["OtherEntities"] = regDim > dim ? regDim : dim;

            nodesInEntities.insert(entityNodes_[actName].begin(), entityNodes_[actName].end());
          }
        }
      } else {
        nodeNames_.push_back( actName );

        // read nodes from file and add to grid
        Hdf5Common::ReadArray( actEntityGroup, "Nodes", entityNodes_[actName] );
      }

      actEntityGroup.close();
    } // loop over entity groups

    std::copy(nodesInEntities.begin(),
              nodesInEntities.end(),
              std::back_inserter(regionNodes_["OtherEntities"]));
    
    entityGroup.close();
    statsRead_ = true;
  }
/**
  void Hdf5Reader::LinearizeElems(const std::vector<unsigned int>& readElems,
                                    std::vector<int>& elemTypes, 
                                    std::vector<unsigned int>& globConnect, 
                                    std::vector<unsigned int>& readNodes) {
    static std::map<FEType, FEType> elemTypeMap;
    unsigned int elemIncr = globConnect.GetSize() / elemTypes.GetSize();
    unsigned int numElems = readElems.GetSize();
    std::set<unsigned int> readNodeSet;
    
    if(!elemTypeMap.size()) {
      elemTypeMap[ET_LINE2] = ET_LINE2;
      elemTypeMap[ET_LINE3] = ET_LINE2;
      elemTypeMap[ET_TRIA3] = ET_TRIA3;
      elemTypeMap[ET_TRIA6] = ET_TRIA3;
      elemTypeMap[ET_QUAD4] = ET_QUAD4;
      elemTypeMap[ET_QUAD8] = ET_QUAD4;
      elemTypeMap[ET_TET4] = ET_TET4;
      elemTypeMap[ET_TET10] = ET_TET4;
      elemTypeMap[ET_HEXA8] = ET_HEXA8;
      elemTypeMap[ET_HEXA20] = ET_HEXA8;
      elemTypeMap[ET_PYRA5] = ET_PYRA5;
      elemTypeMap[ET_PYRA13] = ET_PYRA5;
      elemTypeMap[ET_WEDGE6] = ET_WEDGE6;
      elemTypeMap[ET_WEDGE15] = ET_WEDGE6;
    }

    for(unsigned int i=0; i<numElems; i++) {
      unsigned int elemIdx=readElems[i]-1;
      FEType newType = elemTypeMap[(FEType)elemTypes[elemIdx]];
      unsigned int numNodes = NUM_ELEM_NODES[newType];
        
      elemTypes[elemIdx] = newType;
      readNodeSet.insert(&globConnect[elemIdx*elemIncr],
                         &globConnect[elemIdx*elemIncr+numNodes]);
    }
    
    readNodes.Resize(readNodeSet.size());
    std::copy(readNodeSet.begin(), readNodeSet.end(), readNodes.Begin());
  }
**/  
}
